# Krogh-cylinder O2 model of C. elegans body size

Computes the largest body diameter a nematode can sustain by radial O2 diffusion alone, the length caps that follow from aspect ratio and neural signalling range, and the metabolic headroom at the observed size. Pure Python, no dependencies beyond the standard library; matplotlib only for the optional plots.

## Files

| File | Contents |
| --- | --- |
| [krogh_worm_model.md](krogh_worm_model.md) | The model: derivation, parameter values, predicted numbers |
| [worm_params.py](worm_params.py) | Parameter record and all numeric functions |
| [worm_report.py](worm_report.py) | Unit conversion, the tables from §7 of the doc, plots |
| [test_worm_params.py](test_worm_params.py) | Unit tests for each function |
| [test_regression_doc.py](test_regression_doc.py) | Nominal results in §7–§8 and scaling checks from §9 |

## Quick start

Requires Python 3.10 or newer. Run the examples from the repository root.

```python
from worm_params import WormParams, max_worm_size

p = WormParams(pO2=0.08)          # bacterial lawn; other fields at nominal
r = max_worm_size(p)
r.d_max        # 4.2e-4 m
r.L_max        # 2.0e-3 m  (neural cap binds)
r.scope        # 36.6      (q can rise 37x before the core goes anoxic)
r.pO2_min      # 0.0022 atm (lowest ambient O2 the observed radius tolerates)
```

Lengths, times, and concentrations use SI units (m, s, mol m⁻³); oxygen partial pressure uses atm. Convert only for display, via `worm_report`.

## Parameters

`WormParams` is a frozen dataclass. Defaults are the nominal set used throughout the doc.

| Field | Units | Default | Meaning |
| --- | --- | --- | --- |
| `D` | m² s⁻¹ | 2e-9 | O2 diffusivity in tissue |
| `alpha` | mol m⁻³ atm⁻¹ | 1.4 | O2 solubility at 20 °C |
| `pO2` | atm | 0.21 | ambient partial pressure, must lie in [0, 1] |
| `q` | mol m⁻³ s⁻¹ | 0.02 | volumetric O2 consumption |
| `C_crit` | mol m⁻³ | 0 | minimum tolerable core O2, ≥ 0 |
| `P` | m s⁻¹ | inf | cuticle permeability; `inf` means no cuticle resistance |
| `AR` | 1 | 15 | length / diameter |
| `lam` | m | 1e-3 | neural passive length constant |
| `k` | 1 | 2 | `L_max_neural = k * lam` |
| `R_obs` | m | 35e-6 | observed body radius |

Validation runs in `__post_init__`. `P = inf` is the only non-finite value accepted; every formula uses `R / (2P)`, which is `0.0` for `P = inf`, so the resistance helpers need no special case for an infinitely permeable cuticle.

## Functions

All take a `WormParams`. `driving_concentration`, `r_max`, `scope`, `max_worm_size`, and `o2_profile` raise `NoAerobicSize` when `C_amb − C_crit ≤ 0`. `size_envelope` records those cases as `None`; `resistance`, `length_caps`, and `pO2_min` do not require positive driving concentration.

- `driving_concentration(p)` → `Driving(C_amb, dC)`. The one place `alpha * pO2` and `C_amb − C_crit` are computed.
- `resistance(p, R)` = `R²/(4D) + R/(2P)`, in seconds. Shared by everything below.
- `r_max(p)` → `SizeCap(R_max, d_max)`. Closed form `sqrt(4 D dC / q)` when `P = inf`; otherwise the positive root of `q R²/(4D) + q R/(2P) = dC`, written as `2 dC / (b + sqrt(b² + 4 a dC))` so it is stable as `P → 0`.
- `length_caps(p, d_max)` → `LengthCaps(L_max_iso, L_max_neural, L_max)`. Takes `d_max` as an argument and does not call `r_max`, so the diameter model can be swapped.
- `scope(p)` = `dC / (q · resistance(R_obs))`. Equals `(R_max / R_obs)²` when `P = inf`. Below 1 means `R_obs` already exceeds `R_max`.
- `pO2_min(p)` = `(q · resistance(R_obs) + C_crit) / alpha`. Exact inverse of `r_max`: setting `pO2` to this value makes `R_max == R_obs`. Independent of the current `pO2`; may exceed 1 atm.
- `max_worm_size(p)` → `SizeResult(R_max, d_max, L_max_iso, L_max_neural, L_max, scope, pO2_min)`. Composes the above; no printing.
- `o2_profile(R, p, N=100)` → `Profile(r, C, below_crit, feasible, C_surf, C_core)`. `C(r)` on N+1 points with the cuticle drop folded into `C_surf`. `below_crit[i]` flags points under `C_crit`; because the profile is monotone in r these form a prefix from the axis, and `feasible` is False exactly when `R > R_max`.
- `size_envelope(p, pO2_values, q_values)` → nested dict of `SizeResult` over a pO2 × q grid, `None` where no aerobic size exists.

## Display

Tables and summaries use only the standard library. For plots, install matplotlib:

```sh
python -m pip install matplotlib
```

```python
from worm_report import diameter_cap_table, scope_table, result_summary, plot_dmax_vs_pO2, plot_o2_profiles

print(diameter_cap_table(p))   # markdown, µm and mM
print(scope_table(p))
print(result_summary(p))       # one line

import matplotlib.pyplot as plt
plot_dmax_vs_pO2(p)            # d_max vs ambient O2 (%), one line per q, log-log
plot_o2_profiles(p)            # C(r) at 0.5, 1, 1.5 × R_max; dashed where C < C_crit
plt.show()
```

`worm_params.py` never imports `worm_report.py` or matplotlib; a test enforces this.

## Tests

```sh
python -m pip install pytest
python -m pytest -q
```

48 tests. `test_worm_params.py` checks each function against its defining formula, the two branches of `r_max` and their limits, the `pO2_min` round trip, and the profile flags. `test_regression_doc.py` checks the tabulated nominal results and selected claims in §7–§9 of `krogh_worm_model.md` at nominal parameters and asserts agreement to two significant figures (half a unit in the second digit). A factor-of-2 change to the closed form fails 10 of its 22 tests. If you change the physics on purpose, update the doc and that file together.

## Nominal results

At `WormParams()` defaults (air, `q = 0.02`):

- `R_max` 343 µm, `d_max` 686 µm
- `L_max_iso` 10.3 mm, `L_max_neural` 2.0 mm, `L_max` 2.0 mm
- `scope` 96
- `pO2_min` 0.22%

The observed diameter (60 to 80 µm) sits at the diffusion limit for 0.5 to 1% ambient O2, not for air. Details and the parameter provenance are in the doc.
