# Manuscript Supplement Experiments
Date: 2026-03-27 21:48
Families: 25 (17 standard + 8 adversarial)
Reference REPS: 80



========================================================================
# Exp 1: Counter-Factual S_MD — Self-Selection Test
========================================================================

For each Lorentzian center family X, compute S_MD centered on X
and check whether X uniquely ranks #1 across all 25 families.

| N | Center | Center-Rank | #1 Family | Margin | Self-Select? |
|---|--------|-------------|-----------|--------|--------------|
| 16 | Lor2D | 1 | Lor2D | +5.0520 | YES |
| 16 | Lor3D | 1 | Lor3D | +2.6051 | YES |
| 16 | Lor4D | 1 | Lor4D | +3.6522 | YES |
| 16 | Lor5D | 1 | Lor5D | +12.4538 | YES |
| 28 | Lor2D | 1 | Lor2D | +6.3994 | YES |
| 28 | Lor3D | 1 | Lor3D | +14.3601 | YES |
| 28 | Lor4D | 1 | Lor4D | +4.4656 | YES |
| 28 | Lor5D | 1 | Lor5D | +11.5811 | YES |
| 48 | Lor2D | 1 | Lor2D | +29.6296 | YES |
| 48 | Lor3D | 1 | Lor3D | +21.3970 | YES |
| 48 | Lor4D | 1 | Lor4D | +4.4141 | YES |
| 48 | Lor5D | 1 | Lor5D | +19.7793 | YES |
| 64 | Lor2D | 1 | Lor2D | +37.4490 | YES |
| 64 | Lor3D | 1 | Lor3D | +31.3780 | YES |
| 64 | Lor4D | 1 | Lor4D | +20.0927 | YES |
| 64 | Lor5D | 1 | Lor5D | +22.7525 | YES |
| 96 | Lor2D | 1 | Lor2D | +31.8156 | YES |
| 96 | Lor3D | 1 | Lor3D | +72.0002 | YES |
| 96 | Lor4D | 1 | Lor4D | +26.3260 | YES |
| 96 | Lor5D | 1 | Lor5D | +70.7637 | YES |
| 128 | Lor2D | 1 | Lor2D | +54.7496 | YES |
| 128 | Lor3D | 1 | Lor3D | +88.6509 | YES |
| 128 | Lor4D | 1 | Lor4D | +40.4729 | YES |
| 128 | Lor5D | 1 | Lor5D | +60.8419 | YES |


========================================================================
# Exp 2: Scaling Laws with Bootstrap Uncertainties
========================================================================

## Power-law fits: quantity ~ N^exponent

| Quantity | Exponent | Bootstrap SE | R² | Data Points |
|----------|----------|--------------|----|-------------|
| det(Sigma) | -3.314 ± 0.144 | 0.1440 | 0.9832 | 10 |
| V_eff = sqrt(det) | -1.657 ± 0.072 | 0.0720 | 0.9832 | 10 |
| Fisher I_F = tr(Sigma^-1) | 1.123 ± 0.073 | 0.0732 | 0.9649 | 10 |
| Delta_hist (gap) | 1.218 ± 0.100 | 0.1000 | 0.8924 | 10 |

N values used: [16, 20, 28, 36, 48, 64, 96, 128, 192, 256]

Consistency: det exponent = -3.314, V_eff exponent = -1.657, ½·det_exp = -1.657 (should match V_eff)


========================================================================
# Exp 3: S_BD vs S_MD Functional Separation — Extended
========================================================================

Pearson r and Spearman rho between S_BD and S_MD across families,
with bootstrap 95% CI.

| N | Pearson r | 95% CI | Spearman rho | 95% CI | n_families |
|---|-----------|--------|--------------|--------|------------|
| 16 | -0.586 | [-0.830, +0.014] | -0.207 | [-0.593, +0.305] | 25 |
| 18 | -0.332 | [-0.627, +0.162] | -0.033 | [-0.433, +0.387] | 25 |
| 20 | -0.570 | [-0.802, -0.127] | -0.294 | [-0.649, +0.225] | 25 |
| 24 | -0.065 | [-0.419, +0.525] | +0.175 | [-0.253, +0.565] | 25 |
| 28 | -0.140 | [-0.585, +0.564] | +0.100 | [-0.383, +0.563] | 25 |
| 32 | -0.145 | [-0.681, +0.492] | +0.082 | [-0.383, +0.557] | 25 |
| 36 | +0.235 | [-0.362, +0.557] | +0.244 | [-0.192, +0.597] | 25 |
| 48 | +0.059 | [-0.500, +0.639] | +0.064 | [-0.379, +0.527] | 25 |
| 64 | -0.069 | [-0.442, +0.451] | +0.091 | [-0.327, +0.588] | 25 |
| 96 | -0.386 | [-0.663, +0.044] | -0.132 | [-0.533, +0.305] | 25 |
| 128 | +0.114 | [-0.262, +0.512] | +0.194 | [-0.238, +0.581] | 25 |


========================================================================
# Exp 4: Physical vs Statistical Gap Decomposition
========================================================================

Separating the Mahalanobis gap into:
  - Euclidean distance (absolute feature deviation)
  - Precision amplification (Sigma^-1 scaling)
If Euclidean distance also grows → genuine physical separation.
If only Mahalanobis grows → statistical precision artifact.

| N | Runner-Up | Euclid_dist | Mahal_gap | Sigma_amplif | Euclid_trend |
|---|-----------|-------------|-----------|--------------|-------------|
| 16 | Lor5D | 0.510652 | 3.6522 | 14.01 | — |
| 20 | Lor5D | 0.482315 | 4.6515 | 20.00 | — |
| 28 | Lor5D | 0.369265 | 4.4656 | 32.75 | — |
| 36 | Lor5D | 0.443352 | 14.2345 | 72.42 | — |
| 48 | Lor5D | 0.375122 | 4.4141 | 31.37 | — |
| 64 | Lor5D | 0.488254 | 20.0927 | 84.28 | — |
| 96 | Lor5D | 0.508074 | 26.3260 | 101.98 | — |
| 128 | Lor5D | 0.459235 | 40.4729 | 191.91 | — |
| 192 | Lor5D | 0.456442 | 71.4218 | 342.82 | — |
| 256 | Lor5D | 0.451549 | 92.5763 | 454.04 | — |

Euclidean distance scaling: exponent = 0.007, R² = 0.0032
→ Euclidean distance does NOT vanish — physical separation is genuine
Mahalanobis gap scaling: exponent = 1.218, R² = 0.8924


========================================================================
# Exp 5: Feature Space Ablation
========================================================================

Test S_MD with all subsets of the 3 features:
  Full: (d_eff, C1/C0, w/N)
  Drop-1: each pair of 2 features
  Single: each feature alone

| N | Feature Set | Lor4D Rank | Runner-Up | Margin |
|---|-------------|------------|-----------|--------|
| 28 | d_eff                | 2 | KR_2layer    | +0.5879 |
| 28 | C1/C0                | 1 | Lor4D        | +0.1823 |
| 28 | w/N                  | 2 | KR_like      | +0.9094 |
| 28 | d_eff+C1/C0          | 1 | Lor4D        | +1.6702 |
| 28 | d_eff+w/N            | 1 | Lor4D        | +4.3149 |
| 28 | C1/C0+w/N            | 1 | Lor4D        | +0.1925 |
| 28 | ALL 3                | 1 | Lor4D        | +4.4656 |
| 64 | d_eff                | 2 | KR_2layer    | +0.7917 |
| 64 | C1/C0                | 1 | Lor4D        | +12.0921 |
| 64 | w/N                  | 2 | KR_4layer    | +0.8638 |
| 64 | d_eff+C1/C0          | 1 | Lor4D        | +14.0967 |
| 64 | d_eff+w/N            | 1 | Lor4D        | +14.9914 |
| 64 | C1/C0+w/N            | 1 | Lor4D        | +19.2660 |
| 64 | ALL 3                | 1 | Lor4D        | +20.0927 |
| 128 | d_eff                | 2 | KR_2layer    | +0.1741 |
| 128 | C1/C0                | 1 | Lor4D        | +7.1713 |
| 128 | w/N                  | 1 | Lor4D        | +0.0417 |
| 128 | d_eff+C1/C0          | 1 | Lor4D        | +36.1915 |
| 128 | d_eff+w/N            | 1 | Lor4D        | +28.7820 |
| 128 | C1/C0+w/N            | 1 | Lor4D        | +8.3109 |
| 128 | ALL 3                | 1 | Lor4D        | +40.4729 |


========================================================================
# Exp 6: Turn-On CI and Gap Error Bars (Multi-Seed)
========================================================================

10 independent seeds × REPS=80, full 25-family library

| N | #1_rate | Mean_margin | Margin_SE | Min_margin | 95% CI_margin | Runner_census |
|---|---------|-------------|-----------|------------|---------------|---------------|
| 12 | 10/10 | +0.4757 | 0.1491 | +0.0010 | [+0.183, +0.768] | — |
| 14 | 9/10 | +0.8632 | 0.2488 | -0.0943 | [+0.376, +1.351] | Lor5D |
| 16 | 10/10 | +1.5290 | 0.1723 | +0.6480 | [+1.191, +1.867] | — |
| 18 | 10/10 | +1.9378 | 0.3209 | +0.9746 | [+1.309, +2.567] | — |
| 20 | 10/10 | +1.8878 | 0.3365 | +0.2871 | [+1.228, +2.547] | — |
| 24 | 10/10 | +2.5100 | 0.2219 | +1.7199 | [+2.075, +2.945] | — |
| 28 | 10/10 | +3.2498 | 0.3044 | +1.2072 | [+2.653, +3.846] | — |
| 32 | 10/10 | +4.1237 | 0.3494 | +3.1198 | [+3.439, +4.809] | — |


---
Total runtime: 1044s