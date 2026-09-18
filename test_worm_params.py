import math
from dataclasses import fields, replace
from math import inf
from pathlib import Path

import pytest

import worm_params
from worm_params import (
    NoAerobicSize,
    SizeResult,
    WormParams,
    driving_concentration,
    invert_for_tolerance,
    length_caps,
    max_worm_size,
    o2_profile,
    pO2_min,
    r_max,
    resistance,
    scope,
    size_envelope,
)

NOMINAL = WormParams()  # D=2e-9, alpha=1.4, pO2=0.21, q=0.02, C_crit=0, P=inf


def test_closed_form_nominal_is_340um():
    cap = r_max(NOMINAL)
    assert cap.R_max == pytest.approx(340e-6, abs=5e-6)
    assert cap.d_max == 2 * cap.R_max


def test_closed_form_matches_formula_exactly():
    cap = r_max(NOMINAL)
    d = driving_concentration(NOMINAL)
    assert cap.R_max == pytest.approx(math.sqrt(4 * NOMINAL.D * d.dC / NOMINAL.q), rel=1e-12)


def test_large_finite_P_converges_to_closed_form():
    ref = r_max(NOMINAL).R_max
    errs = [abs(r_max(WormParams(P=P)).R_max - ref) / ref for P in (1e-2, 1e0, 1e2, 1e4)]
    assert errs[-1] < 1e-6
    assert all(e2 < e1 for e1, e2 in zip(errs, errs[1:]))  # monotone convergence


def test_finite_root_satisfies_quadratic():
    p = WormParams(P=1e-5)
    R = r_max(p).R_max
    d = driving_concentration(p)
    resid = p.q * R * R / (4 * p.D) + p.q * R / (2 * p.P) - d.dC
    assert abs(resid) < 1e-12 * d.dC


def test_P_to_zero_sends_R_to_zero():
    Rs = [r_max(WormParams(P=P)).R_max for P in (1e-4, 1e-6, 1e-8, 1e-10)]
    assert all(r2 < r1 for r1, r2 in zip(Rs, Rs[1:]))
    assert Rs[-1] < 1e-8  # P=1e-10 gives R ~ 2.9e-9 m
    # asymptote R ~ 2 P dC / q
    p = WormParams(P=1e-12)
    assert r_max(p).R_max == pytest.approx(2 * p.P * driving_concentration(p).dC / p.q, rel=1e-6)


def test_no_aerobic_size_raises():
    with pytest.raises(NoAerobicSize):
        r_max(WormParams(pO2=0.0))


def test_scope_equals_ratio_squared_when_P_inf():
    res = max_worm_size(NOMINAL)
    assert res.scope == pytest.approx((res.R_max / NOMINAL.R_obs) ** 2, rel=1e-12)
    assert res.scope >= 1


def test_scope_below_one_when_R_obs_exceeds_R_max():
    cap = r_max(NOMINAL).R_max
    res = max_worm_size(WormParams(R_obs=2 * cap))
    assert res.scope == pytest.approx(0.25, rel=1e-12)


def test_length_caps():
    d_max = 4e-4
    L = length_caps(NOMINAL, d_max)
    assert L.L_max_iso == pytest.approx(NOMINAL.AR * d_max)
    assert L.L_max_neural == pytest.approx(NOMINAL.k * NOMINAL.lam)
    assert L.L_max == min(L.L_max_iso, L.L_max_neural)
    # neural cap binds at nominal; isometric cap binds when lam is large
    assert L.L_max == L.L_max_neural
    L2 = length_caps(WormParams(lam=1.0), d_max)
    assert L2.L_max == L2.L_max_iso
    with pytest.raises(ValueError):
        length_caps(NOMINAL, 0.0)


def test_max_worm_size_uses_length_caps():
    res = max_worm_size(NOMINAL)
    L = length_caps(NOMINAL, res.d_max)
    assert (res.L_max_iso, res.L_max_neural, res.L_max) == (L.L_max_iso, L.L_max_neural, L.L_max)


def test_invert_is_exact_inverse_of_r_max():
    for P in (inf, 1e-4, 1e-6):
        p = WormParams(P=P)
        cap = r_max(p).R_max
        q = replace(p, R_obs=cap)
        assert invert_for_tolerance(q) == pytest.approx(p.pO2, rel=1e-10)


def test_profile_core_matches_closed_form():
    p = WormParams(P=1e-5)
    prof = o2_profile(p.R_obs, p, N=10)
    d = driving_concentration(p)
    R = p.R_obs
    assert prof.C_surf == pytest.approx(d.C_amb - p.q * R / (2 * p.P))
    assert prof.C_core == pytest.approx(prof.C_surf - p.q * R * R / (4 * p.D))
    assert prof.C[-1] == pytest.approx(prof.C_surf)
    assert len(prof.r) == 11 and prof.r[-1] == R
    assert all(c2 >= c1 for c1, c2 in zip(prof.C, prof.C[1:]))
    assert prof.feasible and not any(prof.below_crit)
    assert len(prof.below_crit) == len(prof.C)


def test_profile_at_R_max_has_core_at_C_crit():
    p = WormParams(P=1e-5, C_crit=0.01)
    prof = o2_profile(r_max(p).R_max, p)
    assert prof.C_core == pytest.approx(p.C_crit, abs=1e-12)


def test_profile_flags_infeasible_radius():
    p = WormParams(C_crit=0.01)
    Rm = r_max(p).R_max
    ok = o2_profile(0.9 * Rm, p, N=20)
    bad = o2_profile(1.5 * Rm, p, N=20)
    assert ok.feasible
    assert not bad.feasible
    assert bad.below_crit[0]            # axis is the worst point
    assert not bad.below_crit[-1]       # surface is above C_crit
    # flagged points form a prefix
    flags = bad.below_crit
    first_ok = flags.index(False)
    assert all(flags[:first_ok]) and not any(flags[first_ok:])


def test_profile_rejects_bad_args():
    with pytest.raises(ValueError):
        o2_profile(0.0, NOMINAL)
    with pytest.raises(ValueError):
        o2_profile(1e-5, NOMINAL, N=0)
    with pytest.raises(NoAerobicSize):
        o2_profile(1e-5, WormParams(pO2=0.0))


def test_envelope_shape_and_monotonicity():
    env = size_envelope(NOMINAL)
    assert list(env) == [0.005, 0.01, 0.02, 0.05, 0.08, 0.21]
    for pO2, row in env.items():
        assert list(row) == [0.01, 0.02, 0.05]
        Rs = [r.R_max for r in row.values()]
        assert Rs[0] > Rs[1] > Rs[2]  # R_max falls with q
    assert env[0.21][0.02].R_max > env[0.01][0.02].R_max


def test_envelope_marks_no_aerobic_cells():
    env = size_envelope(WormParams(C_crit=0.05), pO2_values=(0.01, 0.21), q_values=(0.02,))
    assert env[0.01][0.02] is None
    assert env[0.21][0.02] is not None


def test_validation_rejects_bad_inputs():
    for kw in (dict(D=0), dict(q=-1), dict(P=0), dict(P=float("nan")), dict(pO2=1.5),
               dict(pO2=-0.1), dict(R_obs=0), dict(C_crit=-1), dict(alpha=inf)):
        with pytest.raises(ValueError):
            WormParams(**kw)
    WormParams(P=inf)  # allowed


def test_resistance_shared_helper():
    p = WormParams(P=1e-5)
    R = p.R_obs
    assert resistance(p, R) == pytest.approx(R * R / (4 * p.D) + R / (2 * p.P), rel=1e-12)
    assert resistance(WormParams(), R) == pytest.approx(R * R / (4 * p.D), rel=1e-12)


def test_scope_is_ratio_squared_when_P_inf():
    p = WormParams()
    assert scope(p) == pytest.approx((r_max(p).R_max / p.R_obs) ** 2, rel=1e-12)


def test_scope_not_ratio_squared_when_P_finite():
    p = WormParams(P=1e-6)
    assert scope(p) != pytest.approx((r_max(p).R_max / p.R_obs) ** 2, rel=1e-3)
    # but scope is exactly dC / (q * resistance)
    assert scope(p) == pytest.approx(driving_concentration(p).dC / (p.q * resistance(p, p.R_obs)), rel=1e-12)


def test_pO2_min_roundtrips_through_r_max():
    for P in (inf, 1e-4, 1e-6):
        for C_crit in (0.0, 0.01):
            p = WormParams(P=P, C_crit=C_crit)
            p2 = replace(p, pO2=pO2_min(p))
            assert r_max(p2).R_max == pytest.approx(p.R_obs, rel=1e-10)
            assert scope(p2) == pytest.approx(1.0, rel=1e-10)


def test_pO2_min_independent_of_current_pO2():
    a, b = WormParams(pO2=0.21), WormParams(pO2=0.05)
    assert pO2_min(a) == pO2_min(b)


def test_result_record_fields_and_composition():
    names = [f.name for f in fields(SizeResult)]
    assert names == ["R_max", "d_max", "L_max_iso", "L_max_neural", "L_max", "scope", "pO2_min"]
    res = max_worm_size(NOMINAL)
    cap = r_max(NOMINAL)
    L = length_caps(NOMINAL, cap.d_max)
    assert (res.R_max, res.d_max) == (cap.R_max, cap.d_max)
    assert (res.L_max_iso, res.L_max_neural, res.L_max) == (L.L_max_iso, L.L_max_neural, L.L_max)
    assert res.scope == scope(NOMINAL)
    assert res.pO2_min == pO2_min(NOMINAL)


def test_max_worm_size_prints_nothing(capsys):
    max_worm_size(NOMINAL)
    out, err = capsys.readouterr()
    assert out == "" and err == ""


def test_report_module_not_imported_by_model():
    src = Path(worm_params.__file__).read_text(encoding="utf-8")
    assert "worm_report" not in src and "matplotlib" not in src
