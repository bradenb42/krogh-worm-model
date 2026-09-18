# Oxygen-diffusion limit on C. elegans body size

A Krogh-cylinder model. Sets a hard cap on body diameter and, through aspect ratio and neural signalling range, a cap on length.

## 1. Notation

| **Symbol** | **Meaning** | **Units** | **Value used** |
| --- | --- | --- | --- |
| R | body radius | m | 35e-6 observed |
| d = 2R | body diameter | m | 60 to 80e-6 observed |
| L | body length | m | 1e-3 observed |
| AR = L/d | aspect ratio | 1 | 15 |
| D | O2 diffusivity in tissue | m² s⁻¹ | 2e-9 |
| α | O2 solubility | mol m⁻³ atm⁻¹ | 1.4 at 20 °C |
| pO2 | ambient partial pressure | atm | 0.21 air, 0.05 to 0.10 lawn |
| C_amb = α·pO2 | ambient dissolved O2 | mol m⁻³ | 0.29 in air |
| C_surf | O2 just inside the cuticle | mol m⁻³ | ≤ C_amb |
| C_crit | minimum tolerable core O2 | mol m⁻³ | 0 (mitochondrial K_m ≈ 1e-3) |
| ΔC = C_amb − C_crit | driving concentration | mol m⁻³ | |
| q | volumetric O2 consumption | mol m⁻³ s⁻¹ | 0.02 nominal, 0.01 to 0.05 range |
| P | cuticle permeability | m s⁻¹ | ∞ if ignored |
| λ | neural passive length constant | m | 1e-3 |
| k | neural length multiplier | 1 | 1 to 3 |

`D` is diffusivity throughout. Diameter is `d`.

## 2. Geometry and transport

The worm is a solid cylinder of radius R and length L with L ≫ R (AR ≈ 15), so axial diffusion is negligible and O2 enters only radially. There is no circulation. Consumption q is uniform. The gut lumen is not an O2 source: pumped bacteria respire and lumen pO2 is low.

## 3. Radial profile

Steady state:

```text
D (1/r) d/dr ( r dC/dr ) = q
```

with dC/dr = 0 at r = 0 and C = C_surf at r = R. Solution:

```text
C(r) = C_surf − q (R² − r²) / (4D)
```

The minimum is on the axis:

```text
C(0) = C_surf − q R² / (4D)
```

## 4. Cuticle

Consumption per unit length is πR²q. Surface area per unit length is 2πR. Flux through the cuticle is therefore

```text
J = q R / 2          [mol m⁻² s⁻¹]
```

and the drop across a cuticle of permeability P is

```text
C_amb − C_surf = J / P = q R / (2P)
```

Combining with §3, the core condition C(0) ≥ C_crit is

```text
q R² / (4D) + q R / (2P) ≤ ΔC
```

The two terms are diffusion resistance R²/(4D) and cuticle resistance R/(2P). Their ratio is R P / (2D). The cuticle dominates when

```text
P < 2D / R
```

For R = 35e-6 and D = 2e-9 that threshold is P ≈ 1e-4 m s⁻¹. Measured C. elegans cuticle permeabilities to small solutes are not well constrained; P = ∞ is used below unless stated.

## 5. Size cap

With P = ∞:

```text
R_max = sqrt( 4 D ΔC / q ),    d_max = 2 R_max
```

With finite P, solve the quadratic in R:

```text
a = q / (4D),  b = q / (2P),  c = −ΔC
R_max = ( −b + sqrt(b² − 4ac) ) / (2a)
```

Metabolic headroom at the observed radius, defined as the factor by which q can rise before C(0) hits C_crit:

```text
q_max = ΔC / ( R_obs² / (4D) + R_obs / (2P) )
scope = q_max / q
```

With P = ∞ this equals (R_max / R_obs)².
Lowest ambient O2 the observed geometry can sustain at consumption q:

```text
pO2_min = ( q R_obs² / (4D) + q R_obs / (2P) + C_crit ) / α
```

## 6. Parameter values

- D: 1.5 to 2e-9 m² s⁻¹. O2 in soft tissue is about two thirds of the value in water.
- α: 1.4 mol m⁻³ atm⁻¹ at 20 °C. Air-saturated C_amb = 0.29 mol m⁻³. A bacterial lawn sits near 5 to 10% O2, C_amb = 0.07 to 0.14.
- q: adult O2 consumption is 2 to 8 pmol min⁻¹ per worm depending on method and feeding state. Body volume for L = 1 mm, d = 60 µm is πR²L = 2.8 nL. That gives q = 0.012 to 0.048 mol m⁻³ s⁻¹ (0.7 to 2.9 mM min⁻¹). Nominal q = 0.02.
- C_crit: mitochondrial O2 K_m is below 1 µM, so C_crit = 0 gives the outer bound.

## 7. Predicted diameter cap (D = 2e-9, q = 0.02, P = ∞)

| **ambient** | **C_amb (mol m⁻³)** | **R_max (µm)** | **d_max (µm)** |
| --- | --- | --- | --- |
| air, 21% | 0.29 | 340 | 690 |
| lawn, 8% | 0.11 | 210 | 420 |
| 2% | 0.028 | 106 | 210 |
| 1% | 0.014 | 75 | 150 |
| 0.5% | 0.007 | 53 | 106 |

Observed d = 60 to 80 µm sits about 10× below the air limit and at the limit for 0.5 to 1% O2. Two observations line up with this. C. elegans tolerates hypoxia down to about 0.5% O2 before entering suspended animation, and it prefers 5 to 12% O2 over air.
Headroom at R_obs = 35 µm:

| **ambient** | **scope = (R_max/R_obs)²** |
| --- | --- |
| air | ≈ 96 |
| lawn, 8% | ≈ 37 |
| 1% | ≈ 4.6 |
| 0.5% | ≈ 2.3 |

Inverting: at q = 0.02, R_obs = 35 µm, C_crit = 0, pO2_min = 0.0022 atm = 0.22%.
R_max ∝ q^(−1/2), so lowering aerobic q by 100× raises d_max by 10×. That is the regime of the large nematodes: Ascaris (d ≈ 5 mm) and Placentonema (L ≈ 8 m) are gut parasites with largely anaerobic metabolism.

## 8. Length

The radial model does not depend on L. A 1 mm and a 10 mm worm of equal diameter have identical O2 profiles. Length is capped only through a coupled constraint.
Isometric growth at fixed AR:

```text
L_max_iso = AR · d_max
```

With AR = 15: 10 mm in air, 6.3 mm on a lawn, 2.2 mm at 1% O2, 1.6 mm at 0.5%. With k = 2 and λ = 1 mm the neural cap is 2.0 mm, so it binds at every ambient level above 1% O2. Observed L ≈ 1 mm matches the hypoxic bound, not the normoxic one.
Neural signalling range. C. elegans neurons do not fire sodium action potentials; they use graded and plateau potentials, so signal range is the passive length constant λ = sqrt(r_m / r_i). For neurites of this diameter λ is a few hundred µm to about 1 mm. A neuron that must span the body passively is limited to

```text
L_max_neural = k · λ,   k ≈ 1 to 3
```

Combined:

```text
L_max = min( L_max_iso, L_max_neural )
```

Undulatory mechanics add a weaker constraint. Bending stiffness of the body scales as R⁴ while muscle torque scales as R³ times muscle fraction. Growth in L at fixed R lowers the achievable wave frequency for a given wavelength.

## 9. Scaling summary

```text
d_max        = 2 sqrt( 4 D ΔC / q )               ∝ (pO2 / q)^(1/2)
L_max_iso    = AR · d_max                         ∝ (pO2 / q)^(1/2)
L_max_neural = k λ
scope        = ΔC / ( q ( R_obs²/(4D) + R_obs/(2P) ) )
```

Doubling q shrinks d_max by 1/√2. Halving pO2 does the same.

## Appendix: pseudocode

```text
MaxWormSize(D, alpha, pO2, q, C_crit, P, AR, lambda, k, R_obs):
  C_amb <- alpha * pO2
  dC    <- C_amb - C_crit
  if dC <= 0: return "no aerobic size possible"
  if P == INF:
    R_max <- sqrt(4*D*dC / q)
  else:
    a <- q/(4*D); b <- q/(2*P); c <- -dC
    R_max <- (-b + sqrt(b*b - 4*a*c)) / (2*a)
  d_max        <- 2*R_max
  L_max_iso    <- AR * d_max
  L_max_neural <- k * lambda
  L_max        <- min(L_max_iso, L_max_neural)
  resist <- R_obs^2/(4*D) + R_obs/(2*P)      // second term 0 if P == INF
  scope  <- dC / (q * resist)
  return R_max, d_max, L_max_iso, L_max_neural, L_max, scope

O2Profile(R, D, C_amb, q, P, N):
  C_surf <- C_amb - q*R/(2*P)
  for i in 0..N:
    r[i] <- R*i/N
    C[i] <- C_surf - q*(R^2 - r[i]^2)/(4*D)
  return r[], C[]

InvertForTolerance(R_obs, D, alpha, q, C_crit, P):
  return (q*R_obs^2/(4*D) + q*R_obs/(2*P) + C_crit) / alpha
```
