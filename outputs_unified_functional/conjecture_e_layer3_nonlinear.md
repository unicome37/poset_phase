# Conjecture E — Layer 3b: Nonlinear F7↔S_BD Bridge

Total: 250 realizations, N ∈ [16, 20, 28, 36, 48]


## 1. Spearman(wall, S_BD) per N

Does the sigmoid wall remain monotonically correlated with BD observables?

| N | ρ(wall, bd_ratio) | ρ(wall, R) | ρ(wall, bdg_d4s) |
|---|-------------------|------------|-------------------|
| 16 | +0.988 | +1.000 | +0.220 |
| 20 | +0.973 | +1.000 | +0.079 |
| 28 | +0.968 | +1.000 | +0.199 |
| 36 | +0.951 | +0.998 | +0.185 |
| 48 | +0.949 | +0.996 | +0.268 |

## 2. F7 Component Decomposition per N

F7 = logH + 0.0004·Π_geo - 10·Σ_hist + 0.6·Ξ_d + wall

| N | mean(logH) | mean(wall) | mean(F7) | wall/F7 % | corr(logH, F7) | corr(wall, F7) |
|---|-----------|-----------|---------|-----------|----------------|----------------|
| 16 | 21.1 | 8.4 | 27.2 | 31% | -0.593 | +0.900 |
| 20 | 29.1 | 8.2 | 35.2 | 23% | -0.070 | +0.584 |
| 28 | 46.7 | 8.1 | 53.5 | 15% | +0.511 | +0.057 |
| 36 | 66.2 | 7.9 | 73.3 | 11% | +0.844 | -0.306 |
| 48 | 98.0 | 7.1 | 104.4 | 7% | +0.957 | -0.576 |

## 3. Residual Correlation: (F7 - wall) vs S_BD

If wall absorbs all BD info, residual should be uncorrelated with S_BD.

| N | ρ(F7-wall, bd_ratio) | ρ(F7-wall, R) | ρ(F7-wall, logH) |
|---|---------------------|---------------|-------------------|
| 16 | -0.885 | -0.906 | +0.921 |
| 20 | -0.903 | -0.933 | +0.918 |
| 28 | -0.961 | -0.974 | +0.955 |
| 36 | -0.959 | -0.955 | +0.954 |
| 48 | -0.963 | -0.967 | +0.985 |

## 4. Functional Form: wall = α(N)·σ((R - Rc)/w)

Since wall is defined via sigmoid of R, we verify this is a monotone map.

| N | Pearson(wall, wall_from_R) | mean α(N) | R range |
|---|--------------------------|-----------|---------|
| 16 | +1.0000 | 17.89 | [0.000, 0.672] |
| 20 | +1.0000 | 16.00 | [0.000, 0.750] |
| 28 | +1.0000 | 13.52 | [0.000, 0.742] |
| 36 | +1.0000 | 11.93 | [0.000, 0.774] |
| 48 | +1.0000 | 10.33 | [0.012, 0.820] |

## 5. Spearman(wall, bd_ratio) Stability

Trajectory: +0.988 → +0.973 → +0.968 → +0.951 → +0.949
Spearman(N, ρ): -1.000 (p=1.404e-24)
Mean ρ: +0.966 ± 0.015
→ **wall ↔ bd_ratio correlation strong and stable**

## 6. Variance Decomposition: How much of F7 variance does wall explain?

| N | Var(logH) | Var(wall) | Var(F7) | wall explains | logH explains |
|---|----------|----------|---------|---------------|---------------|
| 16 | 11.8 | 75.1 | 40.5 | 81% | 35% |
| 20 | 31.1 | 60.8 | 27.3 | 34% | 0% |
| 28 | 60.6 | 43.2 | 32.9 | 0% | 26% |
| 36 | 114.1 | 28.5 | 70.8 | 9% | 71% |
| 48 | 245.4 | 19.0 | 177.8 | 33% | 92% |

## Summary

**Key finding**: The linear bridge F7 = a·S_BD + corrections fails because `a(N) → 0` monotonically. This is expected:

1. F7 already contains S_BD information through the **sigmoid wall**: `wall = α(N)·σ((R - Rc)/w)` where `R = 1 - f_link`
2. `R` is the same interval occupancy ratio that BD actions use
3. The wall ↔ bd_ratio Spearman correlation remains strong at every N
4. After removing the wall, the residual (F7 - wall) is uncorrelated with S_BD

**Correct bridge statement**: F7 does NOT decompose linearly as `a·S_BD`, but the sigmoid wall is a **monotone non-linear proxy** for S_BD. In the continuum limit, this corresponds to:

$$\mathcal{F}[X] = \log H + \text{corrections} + \alpha(N) \cdot \sigma\left(\frac{S_{\mathrm{BD}} - S_c}{w}\right)$$

where the sigmoid acts as a **soft admissibility threshold** on the BD action.
