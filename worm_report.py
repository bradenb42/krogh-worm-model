"""Display helpers for the Krogh worm model: unit conversion, tables, plots.

Nothing here is imported by worm_params. Inputs use model units (SI, with pressure in atm); conversion
happens only at the point of formatting.
"""

from dataclasses import replace
from typing import Iterable, List, Optional, Sequence, Tuple

from worm_params import WormParams, max_worm_size, o2_profile, r_max, scope, NoAerobicSize

# Unit conversion (SI -> display)

def m_to_um(x: float) -> float:
    return x * 1e6


def m_to_mm(x: float) -> float:
    return x * 1e3


def molm3_to_mM(x: float) -> float:
    return x  # 1 mol/m^3 == 1 mmol/L


def atm_to_pct(x: float) -> float:
    return x * 100.0


def pct_to_atm(x: float) -> float:
    return x / 100.0


# Tables (rows as tuples; markdown rendering separate)

def diameter_cap_rows(
    p: WormParams,
    pO2_values: Sequence[float] = (0.21, 0.08, 0.02, 0.01, 0.005),
) -> List[Tuple[float, float, float, float]]:
    """(pO2 [%], C_amb [mM], R_max [um], d_max [um]) per ambient level.

    Levels with no aerobic size are omitted.
    """
    rows = []
    for pO2 in pO2_values:
        params = replace(p, pO2=pO2)
        try:
            cap = r_max(params)
        except NoAerobicSize:
            continue
        rows.append((atm_to_pct(pO2), molm3_to_mM(params.alpha * pO2),
                     m_to_um(cap.R_max), m_to_um(cap.d_max)))
    return rows


def scope_rows(
    p: WormParams,
    pO2_values: Sequence[float] = (0.21, 0.08, 0.01, 0.005),
) -> List[Tuple[float, float]]:
    """(pO2 [%], scope) at R_obs per ambient level. Levels with no aerobic size are omitted."""
    rows = []
    for pO2 in pO2_values:
        try:
            rows.append((atm_to_pct(pO2), scope(replace(p, pO2=pO2))))
        except NoAerobicSize:
            continue
    return rows


def _md_table(header: Sequence[str], rows: Iterable[Sequence[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def diameter_cap_table(p: WormParams, pO2_values: Sequence[float] = (0.21, 0.08, 0.02, 0.01, 0.005)) -> str:
    rows = diameter_cap_rows(p, pO2_values)
    return _md_table(
        ["ambient O2 (%)", "C_amb (mM)", "R_max (µm)", "d_max (µm)"],
        [(f"{pc:g}", f"{c:.3g}", f"{r:.0f}", f"{d:.0f}") for pc, c, r, d in rows],
    )


def scope_table(p: WormParams, pO2_values: Sequence[float] = (0.21, 0.08, 0.01, 0.005)) -> str:
    rows = scope_rows(p, pO2_values)
    return _md_table(
        ["ambient O2 (%)", "scope at R_obs"],
        [(f"{pc:g}", f"{s:.3g}") for pc, s in rows],
    )


def result_summary(p: WormParams) -> str:
    """One-line human-readable summary of max_worm_size(p)."""
    r = max_worm_size(p)
    return (f"R_max {m_to_um(r.R_max):.0f} µm, d_max {m_to_um(r.d_max):.0f} µm, "
            f"L_max_iso {m_to_mm(r.L_max_iso):.1f} mm, L_max_neural {m_to_mm(r.L_max_neural):.1f} mm, "
            f"L_max {m_to_mm(r.L_max):.1f} mm, scope {r.scope:.3g}, "
            f"pO2_min {atm_to_pct(r.pO2_min):.2g}%")


# Plots (matplotlib imported lazily so the module loads without it)

def plot_dmax_vs_pO2(
    p: WormParams,
    q_values: Sequence[float] = (0.01, 0.02, 0.05),
    pO2_values: Optional[Sequence[float]] = None,
    ax=None,
):
    """d_max [um] against ambient O2 [%] on a log-log axis, one line per q.

    Draws the observed diameter 2*R_obs as a horizontal reference.
    Returns the matplotlib Axes.
    """
    import matplotlib.pyplot as plt

    if pO2_values is None:
        pO2_values = [0.002 * 1.15 ** i for i in range(35)]  # 0.2% .. ~23%
        pO2_values = [v for v in pO2_values if v <= 1.0]
    if ax is None:
        _, ax = plt.subplots()
    for q in q_values:
        xs, ys = [], []
        for pO2 in pO2_values:
            try:
                cap = r_max(replace(p, pO2=pO2, q=q))
            except NoAerobicSize:
                continue
            xs.append(atm_to_pct(pO2))
            ys.append(m_to_um(cap.d_max))
        ax.plot(xs, ys, label=f"q = {q:g} mol m⁻³ s⁻¹")
    ax.axhline(m_to_um(2 * p.R_obs), linestyle="--", color="k", label="observed d")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("ambient O2 (%)")
    ax.set_ylabel("d_max (µm)")
    ax.legend()
    return ax


def plot_o2_profiles(
    p: WormParams,
    R_values: Optional[Sequence[float]] = None,
    N: int = 200,
    ax=None,
):
    """C(r) [mM] against r [um] for several radii. Dashed where C < C_crit.

    Default radii: 0.5, 1, 1.5 times R_max. Returns the matplotlib Axes.
    """
    import matplotlib.pyplot as plt

    if R_values is None:
        Rm = r_max(p).R_max
        R_values = (0.5 * Rm, Rm, 1.5 * Rm)
    if ax is None:
        _, ax = plt.subplots()
    for R in R_values:
        prof = o2_profile(R, p, N)
        r_um = [m_to_um(x) for x in prof.r]
        C_mM = [molm3_to_mM(c) for c in prof.C]
        line, = ax.plot(r_um, C_mM, label=f"R = {m_to_um(R):.0f} µm")
        if not prof.feasible:
            k = prof.below_crit.index(False)
            ax.plot(r_um[:k + 1], C_mM[:k + 1], linestyle="--", color=line.get_color())
    ax.axhline(molm3_to_mM(p.C_crit), linestyle=":", color="k", label="C_crit")
    ax.set_xlabel("r (µm)")
    ax.set_ylabel("O2 (mM)")
    ax.legend()
    return ax
