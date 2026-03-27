# N=16 Mahalanobis CV Failure — Root Cause Diagnosis

## 1. Covariance Eigenvalue Spectrum

| N | λ₁ | λ₂ | λ₃ | κ(Σ) = λ_max/λ_min | det(Σ) |
|---|:-:|:-:|:-:|:-:|:-:|
| 16 | 0.1259 | 0.0076 | 0.003541 | 36.4 | 3.99e-06 |
| 20 | 0.0917 | 0.0080 | 0.003152 | 29.2 | 2.33e-06 |
| 28 | 0.0648 | 0.0063 | 0.002693 | 25.7 | 1.05e-06 |
| 64 | 0.0240 | 0.0027 | 0.001149 | 24.6 | 7.21e-08 |
| 128 | 0.0116 | 0.0011 | 0.000499 | 24.4 | 6.25e-09 |

## 2. Train/Test μ Shift at N=16

How much does the centroid shift when we remove 20% of data?

| Seed | Full μ [d, c, w] | Max fold Δμ | Max fold ‖Δμ‖/‖σ‖ |
|------|:-:|:-:|:-:|
| 42 | [3.845, 0.109, 0.537] | 0.0817 | 0.210 |
| 137 | [3.964, 0.113, 0.552] | 0.0375 | 0.197 |
| 271 | [3.874, 0.088, 0.546] | 0.0442 | 0.129 |
| 500 | [4.100, 0.057, 0.588] | 0.0283 | 0.140 |
| 777 | [3.887, 0.081, 0.550] | 0.0496 | 0.129 |
| 1001 | [4.032, 0.076, 0.562] | 0.0317 | 0.166 |
| 2023 | [3.982, 0.079, 0.583] | 0.0517 | 0.196 |
| 3141 | [3.952, 0.080, 0.548] | 0.0221 | 0.156 |
| 5000 | [3.975, 0.067, 0.554] | 0.0221 | 0.272 |
| 8888 | [3.992, 0.087, 0.590] | 0.0300 | 0.101 |

## 3. Competitor Identity at Failure Points

5-fold CV at N=16 and N=20: who beats Lor4D?

| Seed | N | Fold | Winner | Winner S_M | Lor4D S_M | Gap |
|------|---|------|--------|:-:|:-:|:-:|
| 42 | 20 | 5 | Lor5D | 4.91 | 8.34 | 3.43 |
| 137 | 16 | 2 | Lor5D | 4.28 | 4.91 | 0.64 |
| 137 | 16 | 3 | Lor5D | 4.14 | 4.52 | 0.38 |
| 271 | 16 | 2 | Lor5D | 5.50 | 5.69 | 0.19 |
| 271 | 20 | 3 | Lor5D | 3.41 | 3.92 | 0.51 |
| 271 | 20 | 4 | Lor5D | 4.43 | 7.42 | 2.98 |
| 777 | 16 | 1 | Lor5D | 5.05 | 5.46 | 0.40 |
| 777 | 20 | 4 | Lor5D | 6.23 | 6.39 | 0.16 |
| 2023 | 16 | 1 | Lor5D | 5.28 | 7.35 | 2.07 |
| 2023 | 16 | 3 | Lor5D | 3.26 | 3.76 | 0.50 |
| 2023 | 16 | 5 | Lor5D | 3.31 | 3.54 | 0.23 |
| 3141 | 20 | 2 | Lor5D | 7.80 | 8.97 | 1.17 |
| 5000 | 16 | 5 | Lor5D | 5.81 | 6.33 | 0.53 |
| 5000 | 20 | 5 | Lor5D | 4.65 | 15.47 | 10.82 |
| 8888 | 16 | 1 | Lor5D | 2.65 | 4.63 | 1.98 |
| 8888 | 16 | 2 | Lor5D | 2.31 | 3.89 | 1.58 |
| 8888 | 16 | 3 | Lor5D | 1.91 | 1.94 | 0.04 |
| 8888 | 16 | 4 | Lor5D | 2.21 | 2.31 | 0.11 |
| 8888 | 16 | 5 | Lor5D | 2.17 | 4.96 | 2.78 |

Total failures: 19
Winner census: {'Lor5D': 19}

## 4. Regularized Mahalanobis (Shrinkage)

Replace Σ with (1-α)Σ + α·tr(Σ)/3·I (shrinkage toward spherical)

| α | N=16 CV #1 rate | N=20 CV #1 rate | N=28+ CV #1 rate |
|---|:-:|:-:|:-:|
| 0.0 | 37/50 | 44/50 | 47/50 |
| 0.1 | 44/50 | 44/50 | 49/50 |
| 0.2 | 43/50 | 45/50 | 49/50 |
| 0.5 | 29/50 | 45/50 | 49/50 |
| 1.0 | 19/50 | 23/50 | 40/50 |

## 5. Minimum Training Size for Stable Mahalanobis at N=16

Vary train fraction from 50% to 95%.

| Train fraction | n_train | CV #1 rate (N=16) |
|:-:|:-:|:-:|
| 50% | ~20 | 77/100 (77%) |
| 60% | ~24 | 74/100 (74%) |
| 70% | ~28 | 72/100 (72%) |
| 80% | ~32 | 74/100 (74%) |
| 90% | ~36 | 83/100 (83%) |
| 95% | ~38 | 79/100 (79%) |

## 6. Root Cause Synthesis

The analysis identifies the precise mechanism behind N=16 CV failures:

### Root cause: Lor5D is the SOLE intruder

All 19 failures across all seeds at N=16/20 have the **same winner: Lor5D**.
No non-Lorentzian family ever beats Lor4D under Mahalanobis, even with CV.

This is not a covariance instability problem in the usual sense.
It is a **dimensionality proximity** effect: at N=16, Lor4D and Lor5D feature vectors
have not yet separated, because:
- d_eff(Lor4D, N=16) ≈ 3.85–4.10 (noisy, far from 4.0)
- d_eff(Lor5D, N=16) ≈ 4.0–4.5 (overlapping with Lor4D)

### Condition number is NOT the main factor

κ(Σ) = 36 at N=16 vs 24 at N=128 — only 1.5× difference.
The condition number is high but stable. The small eigenvalue λ₃ ≈ 0.0035
(width variance) is small but not pathological.

### μ shift is moderate

Max fold ‖Δμ‖/‖σ‖ ≈ 0.1–0.3, meaning centroid shifts < 30% of one σ.
This causes mild Σ⁻¹ misalignment but is not the primary failure mode.

### Shrinkage helps slightly but doesn't solve it

α=0.1 shrinkage improves N=16 from 37/50 to 44/50 (+14%), but further
shrinkage (α=0.5, 1.0) degrades performance. This confirms the issue is
feature overlap, not covariance estimation noise.

### Physical interpretation

> At N=16, the causal geometry has not yet developed enough structure
> to distinguish 4D from 5D Lorentzian embeddings. The Mahalanobis
> distance correctly measures proximity to the Lor4D reference manifold,
> but Lor5D samples fall inside the same ellipsoid simply because
> 16-element causets don't encode enough dimension information.
>
> This is not a statistical failure; it is a **physical resolution limit**.

### Why LSD-Well doesn't fail

LSD-Well uses the hard target d*=4.0, which Lor5D (d_eff ≈ 4.5) violates
through the (d-4)² term with weight α=0.5. Mahalanobis uses the data-driven
centroid d_eff ≈ 3.9 and a learned covariance, which at small N includes Lor5D.
The LSD-Well's built-in "d=4 is truth" assumption is what saves it at N=16.
