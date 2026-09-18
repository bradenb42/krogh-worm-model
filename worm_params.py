"""Krogh-cylinder O2 model of C. elegans body size.

Parameter record (WormParams), driving concentration, radius cap (r_max), length caps
(length_caps), metabolic headroom (scope), minimum tolerable ambient O2 for a
given geometry (pO2_min), the combined result (max_worm_size), radial O2
profile (o2_profile), and a pO2 x q sweep (size_envelope). SI units except pressure, which is in atm.
"""

from dataclasses import dataclass, field, replace
from math import inf, isfinite, sqrt
from typing import Dict, Sequence, Tuple

__all__ = [
    "WormParams", "NoAerobicSize", "Driving", "SizeCap", "LengthCaps", "SizeResult", "Profile",
    "driving_concentration", "r_max", "length_caps", "max_worm_size", "o2_profile",
    "resistance", "scope", "pO2_min", "invert_for_tolerance", "size_envelope",
]


@dataclass(frozen=True)
class WormParams:
    """All inputs to max_worm_size, o2_profile and invert_for_tolerance.

    Units are SI except pressure, which is in atm. P = inf means no cuticle resistance; every
    formula uses R / (2 * P), which evaluates to 0.0 for P = inf.
    """

    D: float = field(default=2e-9, metadata={"units": "m^2/s", "doc": "O2 diffusivity in tissue"})
    alpha: float = field(default=1.4, metadata={"units": "mol/m^3/atm", "doc": "O2 solubility"})
    pO2: float = field(default=0.21, metadata={"units": "atm", "doc": "ambient O2 partial pressure, 0..1"})
    q: float = field(default=0.02, metadata={"units": "mol/m^3/s", "doc": "volumetric O2 consumption"})
    C_crit: float = field(default=0.0, metadata={"units": "mol/m^3", "doc": "minimum tolerable core O2, >= 0"})
    P: float = field(default=inf, metadata={"units": "m/s", "doc": "cuticle permeability; inf = no cuticle"})
    AR: float = field(default=15.0, metadata={"units": "1", "doc": "length / diameter aspect ratio"})
    lam: float = field(default=1e-3, metadata={"units": "m", "doc": "neural passive length constant (lambda)"})
    k: float = field(default=2.0, metadata={"units": "1", "doc": "neural length multiplier, L_max_neural = k * lam"})
    R_obs: float = field(default=35e-6, metadata={"units": "m", "doc": "observed body radius"})

    def __post_init__(self) -> None:
        strictly_positive = ("D", "alpha", "q", "AR", "lam", "k", "R_obs")
        for name in strictly_positive:
            v = getattr(self, name)
            if not isfinite(v) or v <= 0:
                raise ValueError(f"{name} must be finite and > 0, got {v!r}")
        if not (isfinite(self.C_crit) and self.C_crit >= 0):
            raise ValueError(f"C_crit must be finite and >= 0, got {self.C_crit!r}")
        if not (self.P > 0):  # inf passes; nan, 0, negative fail
            raise ValueError(f"P must be > 0 (inf allowed), got {self.P!r}")
        if not (isfinite(self.pO2) and 0.0 <= self.pO2 <= 1.0):
            raise ValueError(f"pO2 must be in [0, 1] atm, got {self.pO2!r}")

    def cuticle_resistance(self, R: float) -> float:
        """R / (2 P) [s]. Returns 0.0 when P is inf."""
        return R / (2.0 * self.P)

    def diffusion_resistance(self, R: float) -> float:
        """R^2 / (4 D) [s]."""
        return R * R / (4.0 * self.D)

    def resistance(self, R: float) -> float:
        """Total series resistance R^2/(4D) + R/(2P) [s]."""
        return self.diffusion_resistance(R) + self.cuticle_resistance(R)


class NoAerobicSize(ValueError):
    """Raised when C_amb - C_crit <= 0: no radius satisfies the core condition."""


@dataclass(frozen=True)
class Driving:
    C_amb: float  # alpha * pO2 [mol/m^3]
    dC: float     # C_amb - C_crit [mol/m^3], > 0


def driving_concentration(p: WormParams) -> Driving:
    """Return C_amb = alpha * pO2 and dC = C_amb - C_crit.

    Raises NoAerobicSize when dC <= 0. Radius, scope and profile
    calculations obtain C_amb and dC from here.
    """
    C_amb = p.alpha * p.pO2
    dC = C_amb - p.C_crit
    if dC <= 0.0:
        raise NoAerobicSize(
            f"C_amb - C_crit = {C_amb:.4g} - {p.C_crit:.4g} = {dC:.4g} mol/m^3 <= 0"
        )
    return Driving(C_amb=C_amb, dC=dC)


@dataclass(frozen=True)
class SizeCap:
    R_max: float  # [m]
    d_max: float  # 2 * R_max [m]


def r_max(p: WormParams) -> SizeCap:
    """Largest radius whose core stays at or above C_crit.

    P infinite: closed form R = sqrt(4 D dC / q).
    P finite:   positive root of q R^2/(4D) + q R/(2P) - dC = 0, written as
                R = 2 dC / (b + sqrt(b^2 + 4 a dC)) with a = q/(4D), b = q/(2P),
                which is finite for every P > 0 and tends to 0 as P -> 0.
    """
    d = driving_concentration(p)
    if p.P == inf:
        R = sqrt(4.0 * p.D * d.dC / p.q)
    else:
        a = p.q / (4.0 * p.D)
        b = p.q / (2.0 * p.P)
        R = 2.0 * d.dC / (b + sqrt(b * b + 4.0 * a * d.dC))
    return SizeCap(R_max=R, d_max=2.0 * R)


@dataclass(frozen=True)
class LengthCaps:
    L_max_iso: float     # AR * d_max [m]
    L_max_neural: float  # k * lam [m]
    L_max: float         # min of the two [m]


def length_caps(p: WormParams, d_max: float) -> LengthCaps:
    """Length caps from a diameter cap and the neural length scale.

    L_max_iso    = AR * d_max   (isometric growth at fixed aspect ratio)
    L_max_neural = k * lam      (graded signalling over a passive neurite)
    L_max        = min(L_max_iso, L_max_neural)

    d_max comes from any diameter model; this function does not call r_max.
    """
    if not (isfinite(d_max) and d_max > 0):
        raise ValueError(f"d_max must be finite and > 0, got {d_max!r}")
    L_iso = p.AR * d_max
    L_neu = p.k * p.lam
    return LengthCaps(L_max_iso=L_iso, L_max_neural=L_neu, L_max=min(L_iso, L_neu))


@dataclass(frozen=True)
class SizeResult:
    R_max: float          # [m]
    d_max: float          # [m]
    L_max_iso: float      # AR * d_max [m]
    L_max_neural: float   # k * lam [m]
    L_max: float          # min of the two [m]
    scope: float          # q_max / q at R_obs, dimensionless
    pO2_min: float        # lowest ambient pO2 sustaining R_obs at q [atm]; may exceed 1


def max_worm_size(p: WormParams) -> SizeResult:
    """Compose r_max, length_caps, scope and pO2_min into one result record.

    Raises NoAerobicSize when dC <= 0. scope < 1 means R_obs exceeds R_max.
    """
    cap = r_max(p)
    L = length_caps(p, cap.d_max)
    return SizeResult(
        R_max=cap.R_max,
        d_max=cap.d_max,
        L_max_iso=L.L_max_iso,
        L_max_neural=L.L_max_neural,
        L_max=L.L_max,
        scope=scope(p),
        pO2_min=pO2_min(p),
    )


@dataclass(frozen=True)
class Profile:
    r: Tuple[float, ...]           # radial positions 0..R [m], N+1 points
    C: Tuple[float, ...]           # O2 concentration at r [mol/m^3]
    below_crit: Tuple[bool, ...]   # C[i] < C_crit
    feasible: bool                 # no point below C_crit
    C_surf: float                  # C_amb - q R/(2P) [mol/m^3]
    C_core: float                  # C[0] [mol/m^3]


def o2_profile(R: float, p: WormParams, N: int = 100) -> Profile:
    """Radial O2 concentration inside a cylinder of radius R.

    C_surf = C_amb - q R/(2P);  C(r) = C_surf - q (R^2 - r^2)/(4D) on N+1
    points r_i = R i/N. below_crit[i] is True where C(r_i) < C_crit; because
    C rises monotonically with r, the flagged points form a prefix starting at
    the axis, and feasible is False exactly when R > R_max(p).
    Raises NoAerobicSize when C_amb - C_crit <= 0 (no R is feasible).
    """
    if not (isfinite(R) and R > 0):
        raise ValueError(f"R must be finite and > 0, got {R!r}")
    if N < 1:
        raise ValueError(f"N must be >= 1, got {N!r}")
    d = driving_concentration(p)
    C_surf = d.C_amb - p.q * p.cuticle_resistance(R)
    r = tuple(R * i / N for i in range(N + 1))
    C = tuple(C_surf - p.q * (R * R - ri * ri) / (4.0 * p.D) for ri in r)
    below = tuple(c < p.C_crit for c in C)
    return Profile(r=r, C=C, below_crit=below, feasible=not any(below),
                   C_surf=C_surf, C_core=C[0])


def resistance(p: WormParams, R: float) -> float:
    """R^2/(4D) + R/(2P) [s]. Shared by scope and pO2_min."""
    return p.resistance(R)


def scope(p: WormParams) -> float:
    """Metabolic headroom at R_obs: dC / (q * resistance(R_obs)).

    Factor by which q can rise before the core reaches C_crit. Equals
    (R_max / R_obs)^2 when P = inf. < 1 means R_obs already exceeds R_max.
    Raises NoAerobicSize when dC <= 0.
    """
    d = driving_concentration(p)
    return d.dC / (p.q * resistance(p, p.R_obs))


def pO2_min(p: WormParams) -> float:
    """Lowest ambient pO2 [atm] keeping R_obs aerobic at consumption q.

    (q * resistance(R_obs) + C_crit) / alpha. Exact inverse of r_max: replacing
    pO2 with this value makes R_max == R_obs. Independent of the current pO2;
    may exceed 1 atm.
    """
    return (p.q * resistance(p, p.R_obs) + p.C_crit) / p.alpha


invert_for_tolerance = pO2_min  # backward-compatible alias


def size_envelope(
    p: WormParams,
    pO2_values: Sequence[float] = (0.005, 0.01, 0.02, 0.05, 0.08, 0.21),
    q_values: Sequence[float] = (0.01, 0.02, 0.05),
) -> Dict[float, Dict[float, SizeResult | None]]:
    """max_worm_size over a pO2 x q grid, other parameters from p.

    A cell is None where dC <= 0 (no aerobic size).
    """
    out: Dict[float, Dict[float, SizeResult | None]] = {}
    for pO2 in pO2_values:
        row: Dict[float, SizeResult | None] = {}
        for q in q_values:
            try:
                row[q] = max_worm_size(replace(p, pO2=pO2, q=q))
            except NoAerobicSize:
                row[q] = None
        out[pO2] = row
    return out
