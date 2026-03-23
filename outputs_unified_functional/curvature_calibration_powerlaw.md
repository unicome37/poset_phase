# Curvature Calibration: Power-Law Fitting

## Motivation

The naive ratio c = R_hat/R_dS grows with H, indicating R_hat is NOT
simply proportional to H². We fit R_hat = c·H^α and check if α→2.


## Strategy A: Fit mean(R_hat) vs H per (d, N)

Using bias-subtracted R_hat: R_hat_corr = R_hat - R_hat(H=0)

| d | N | α (exponent) | log(c) | R² | α_raw (no bias sub) |
|---|---|-------------|--------|-----|---------------------|
| 2 | 256 | **3.359** ± 0.345 | 4.93 | 0.9896 | nan |
| 2 | 512 | **3.809** ± 0.095 | 4.76 | 0.9994 | 3.502 |
| 2 | 1024 | **3.842** ± 0.056 | 4.73 | 0.9998 | 3.589 |
| 2 | 2048 | **4.013** ± 0.035 | 4.68 | 0.9999 | nan |
| 3 | 256 | **3.069** ± 0.924 | 4.49 | 0.9168 | 2.960 |
| 3 | 512 | **3.517** ± 0.000 | 4.28 | 1.0000 | nan |
| 3 | 1024 | **3.592** ± 0.019 | 4.46 | 1.0000 | 3.259 |
| 3 | 2048 | **4.310** ± 0.280 | 4.10 | 0.9958 | 3.407 |
| 4 | 256 | **1.423** ± 0.020 | 3.84 | 0.9998 | nan |
| 4 | 512 | **4.588** ± 1.296 | 3.18 | 0.9261 | 3.823 |
| 4 | 1024 | **2.933** ± 0.264 | 3.37 | 0.9919 | 2.500 |
| 4 | 2048 | **2.930** ± 0.000 | 3.62 | 1.0000 | nan |

## α(d, N) Convergence

Does α → 2.0 (linear in R_dS) as N → ∞?


### d=2

α(N): 3.359 → 3.809 → 3.842 → 4.013
- Mean α: 3.756 ± 0.242
- Target: α = 2.0 (proportional to H²)
- Deviation from 2: +1.756
- Spearman(N, α): ρ=+1.000 (p=0.000e+00)

### d=3

α(N): 3.069 → 3.517 → 3.592 → 4.310
- Mean α: 3.622 ± 0.445
- Target: α = 2.0 (proportional to H²)
- Deviation from 2: +1.622
- Spearman(N, α): ρ=+1.000 (p=0.000e+00)

### d=4

α(N): 1.423 → 4.588 → 2.933 → 2.930
- Mean α: 2.968 ± 1.120
- Target: α = 2.0 (proportional to H²)
- Deviation from 2: +0.968
- Spearman(N, α): ρ=+0.200 (p=8.000e-01)
- → α stable (no N trend)

## Strategy B: Linear Slope in H² (recap from Layer 2)

If R_hat = c·H² + bias, the slope is the calibration constant.

- d=2, N=256: slope(R_hat vs H²) = 428.7, R² = 0.913, c_eff = slope/2 = **214.4**
- d=2, N=512: slope(R_hat vs H²) = 417.1, R² = 0.930, c_eff = slope/2 = **208.6**
- d=2, N=1024: slope(R_hat vs H²) = 438.2, R² = 0.944, c_eff = slope/2 = **219.1**
- d=2, N=2048: slope(R_hat vs H²) = 455.9, R² = 0.971, c_eff = slope/2 = **227.9**
- d=3, N=256: slope(R_hat vs H²) = 286.2, R² = 0.377, c_eff = slope/6 = **47.7**
- d=3, N=512: slope(R_hat vs H²) = 219.9, R² = 0.818, c_eff = slope/6 = **36.6**
- d=3, N=1024: slope(R_hat vs H²) = 272.1, R² = 0.898, c_eff = slope/6 = **45.4**
- d=3, N=2048: slope(R_hat vs H²) = 283.2, R² = 0.959, c_eff = slope/6 = **47.2**
- d=4, N=256: slope(R_hat vs H²) = -36.9, R² = 0.128, c_eff = slope/12 = **-3.1**
- d=4, N=512: slope(R_hat vs H²) = 89.6, R² = 0.154, c_eff = slope/12 = **7.5**
- d=4, N=1024: slope(R_hat vs H²) = 52.3, R² = 0.133, c_eff = slope/12 = **4.4**
- d=4, N=2048: slope(R_hat vs H²) = 76.1, R² = 0.816, c_eff = slope/12 = **6.3**

## Strategy C: Per-H Calibration at N=2048

Best estimate using largest N.

| d | H | R_dS | mean(R_hat) | bias-corrected | c = corr/R_dS |
|---|---|------|-------------|----------------|---------------|
| 2 | 0.25 | 0.12 | +0.3 | -3.6 | -28.9 |
| 2 | 0.5 | 0.50 | +10.6 | +6.6 | +13.2 |
| 2 | 1.0 | 2.00 | +115.2 | +111.2 | +55.6 |
| 2 | 2.0 | 8.00 | +1726.5 | +1722.6 | +215.3 |
| 3 | 0.25 | 0.38 | -1.8 | -8.7 | -23.2 |
| 3 | 0.5 | 1.50 | +9.6 | +2.7 | +1.8 |
| 3 | 1.0 | 6.00 | +82.1 | +75.3 | +12.5 |
| 3 | 2.0 | 24.00 | +1074.2 | +1067.3 | +44.5 |
| 4 | 0.25 | 0.75 | -0.3 | -10.2 | -13.6 |
| 4 | 0.5 | 3.00 | +7.0 | -2.9 | -1.0 |
| 4 | 1.0 | 12.00 | +47.1 | +37.2 | +3.1 |
| 4 | 2.0 | 48.00 | +293.7 | +283.9 | +5.9 |

## Summary

- **d=2**: α(N=2048) = **4.013** (target=2.0, deviation=+2.013, R²=1.000)
- **d=3**: α(N=2048) = **4.310** (target=2.0, deviation=+2.310, R²=0.996)
- **d=4**: α(N=2048) = **2.930** (target=2.0, deviation=+0.930, R²=1.000)

**Interpretation**: If α > 2, R_hat responds super-linearly to curvature —
the discrete proxy amplifies strong curvature. If α < 2, it under-responds.
Convergence α → 2 as N → ∞ would confirm correct continuum scaling.
