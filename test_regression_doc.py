"""A change to the physics in worm_params must break at least one assertion here.
Agreement criterion: |computed - doc| <= half a unit in the doc value's second
significant figure (i.e. the doc value is the computed value rounded to 2 s.f.).
"""

from dataclasses import replace
from math import floor, log10

import pytest

from worm_params import WormParams, r_max, scope, pO2_min, length_caps, max_worm_size

NOMINAL = WormParams()  # D=2e-9, alpha=1.4, q=0.02, C_crit=0, P=inf, AR=15, lam=1e-3, k=2, R_obs=35e-6


def sf2(computed: float, doc: float) -> None:
    """Assert computed agrees with doc to two significant figures."""
    unit = 10 ** (floor(log10(abs(doc))) - 1)
    assert abs(computed - doc) <= 0.5 * unit + 1e-15, f"{computed!r} vs doc {doc!r} (2 s.f.)"


# ---- §7 table: ambient -> C_amb, R_max, d_max ---------------------------------
# pO2 [atm], C_amb [mol/m^3], R_max [um], d_max [um]
SEC7_TABLE = [
    (0.21,  0.29,  340, 690),
    (0.08,  0.11,  210, 420),
    (0.02,  0.028, 106, 210),
    (0.01,  0.014,  75, 150),
    (0.005, 0.007,  53, 106),
]


@pytest.mark.parametrize("pO2,C_amb,R_um,d_um", SEC7_TABLE)
def test_sec7_diameter_table(pO2, C_amb, R_um, d_um):
    p = replace(NOMINAL, pO2=pO2)
    sf2(p.alpha * p.pO2, C_amb)
    cap = r_max(p)
    sf2(cap.R_max * 1e6, R_um)
    sf2(cap.d_max * 1e6, d_um)


# ---- §7 scope table ------------------------------------------------------------
SEC7_SCOPE = [(0.21, 96), (0.08, 37), (0.01, 4.6), (0.005, 2.3)]


@pytest.mark.parametrize("pO2,s", SEC7_SCOPE)
def test_sec7_scope_table(pO2, s):
    sf2(scope(replace(NOMINAL, pO2=pO2)), s)


def test_sec7_observed_diameter_about_10x_below_air_limit():
    ratio = r_max(NOMINAL).d_max / (2 * NOMINAL.R_obs)
    assert 9.0 <= ratio <= 11.0  # doc says "about 10x"


def test_sec7_pO2_min():
    sf2(pO2_min(NOMINAL), 0.0022)          # atm
    sf2(pO2_min(NOMINAL) * 100, 0.22)      # percent


def test_sec7_q_scaling_100x_gives_10x_diameter():
    d1 = r_max(NOMINAL).d_max
    d2 = r_max(replace(NOMINAL, q=NOMINAL.q / 100)).d_max
    assert d2 / d1 == pytest.approx(10.0, rel=1e-12)


# ---- §8 length caps ------------------------------------------------------------
# pO2 [atm], L_max_iso [mm]
SEC8_ISO = [(0.21, 10), (0.08, 6.3), (0.01, 2.2), (0.005, 1.6)]


@pytest.mark.parametrize("pO2,L_mm", SEC8_ISO)
def test_sec8_isometric_length(pO2, L_mm):
    p = replace(NOMINAL, pO2=pO2)
    L = length_caps(p, r_max(p).d_max)
    sf2(L.L_max_iso * 1e3, L_mm)
    assert L.L_max_iso == pytest.approx(p.AR * r_max(p).d_max, rel=1e-12)


def test_sec8_neural_cap():
    # L_max_neural = k * lambda; k = 1..3, lambda = 1 mm -> 1..3 mm; nominal k=2 -> 2.0 mm
    sf2(length_caps(NOMINAL, 1e-3).L_max_neural * 1e3, 2.0)
    sf2(length_caps(replace(NOMINAL, k=1), 1e-3).L_max_neural * 1e3, 1.0)
    sf2(length_caps(replace(NOMINAL, k=3), 1e-3).L_max_neural * 1e3, 3.0)


@pytest.mark.parametrize("pO2", [0.21, 0.08, 0.02])
def test_sec8_neural_cap_binds_above_1pct(pO2):
    r = max_worm_size(replace(NOMINAL, pO2=pO2))
    assert r.L_max == r.L_max_neural
    assert r.L_max_iso > r.L_max_neural


def test_sec8_radial_model_independent_of_length():
    # §8: "A 1 mm and a 10 mm worm of equal diameter have identical O2 profiles."
    a = r_max(replace(NOMINAL, AR=15))
    b = r_max(replace(NOMINAL, AR=150))
    assert a == b


# ---- §9 scaling laws stated in the doc and used in §7/§8 -----------------------
def test_sec9_scaling_exponents():
    base = r_max(NOMINAL).d_max
    assert r_max(replace(NOMINAL, q=2 * NOMINAL.q)).d_max / base == pytest.approx(2 ** -0.5, rel=1e-12)
    assert r_max(replace(NOMINAL, pO2=NOMINAL.pO2 / 2)).d_max / base == pytest.approx(2 ** -0.5, rel=1e-12)
