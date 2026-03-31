# Prediction C — F7 Large-N Verification

**Data**: 200 samples, N ∈ {20, 36, 52, 72, 100}

**Families**: KR_like, Lor2D, Lor3D, Lor4D, Lor5D


## Q1: Within-Cell Correlation ρ(Σ_hist, logH)

The core Prediction C claim: higher Σ_hist → lower logH within fixed (family, N).

| family | N | n | ρ(Σ,logH) | p_value | sig | direction |
|--------|---|---|-----------|---------|-----|-----------|
| KR_like | 20 | 8 | — | — | — | degenerate |
| Lor2D | 20 | 8 | -0.868 | 0.0052 | ★★ | ✅ neg |
| Lor3D | 20 | 8 | -0.352 | 0.3923 |  | ✅ neg |
| Lor4D | 20 | 8 | -0.504 | 0.2029 |  | ✅ neg |
| Lor5D | 20 | 8 | -0.378 | 0.3559 |  | ✅ neg |
| KR_like | 36 | 8 | — | — | — | degenerate |
| Lor2D | 36 | 8 | -0.482 | 0.2267 |  | ✅ neg |
| Lor3D | 36 | 8 | -0.326 | 0.4306 |  | ✅ neg |
| Lor4D | 36 | 8 | +0.056 | 0.8946 |  | ❌ pos |
| Lor5D | 36 | 8 | — | — | — | degenerate |
| KR_like | 52 | 8 | — | — | — | degenerate |
| Lor2D | 52 | 8 | -0.386 | 0.3453 |  | ✅ neg |
| Lor3D | 52 | 8 | -0.282 | 0.4991 |  | ✅ neg |
| Lor4D | 52 | 8 | +0.412 | 0.3100 |  | ❌ pos |
| Lor5D | 52 | 8 | -0.630 | 0.0941 |  | ✅ neg |
| KR_like | 72 | 8 | — | — | — | degenerate |
| Lor2D | 72 | 8 | -0.250 | 0.5499 |  | ✅ neg |
| Lor3D | 72 | 8 | -0.577 | 0.1340 |  | ✅ neg |
| Lor4D | 72 | 8 | -0.247 | 0.5546 |  | ✅ neg |
| Lor5D | 72 | 8 | +0.056 | 0.8946 |  | ❌ pos |
| KR_like | 100 | 8 | — | — | — | degenerate |
| Lor2D | 100 | 8 | -0.970 | 0.0001 | ★★★ | ✅ neg |
| Lor3D | 100 | 8 | -0.386 | 0.3453 |  | ✅ neg |
| Lor4D | 100 | 8 | -0.252 | 0.5472 |  | ✅ neg |
| Lor5D | 100 | 8 | -0.109 | 0.7970 |  | ✅ neg |

**Summary**: 16/19 cells have ρ < 0 (correct direction), 2/19 significant at p<0.05


## Q2: N-Scaling of Within-Cell ρ(Σ_hist, logH)

Does the correlation strengthen with N?

- **Lor2D**: ρ(N, ρ_cell) = +0.000 (p=1.000) — stable
  - N=20: ρ=-0.868 (p=0.0052)
  - N=36: ρ=-0.482 (p=0.2267)
  - N=52: ρ=-0.386 (p=0.3453)
  - N=72: ρ=-0.250 (p=0.5499)
  - N=100: ρ=-0.970 (p=0.0001)
- **Lor3D**: ρ(N, ρ_cell) = -0.500 (p=0.391) — stable
  - N=20: ρ=-0.352 (p=0.3923)
  - N=36: ρ=-0.326 (p=0.4306)
  - N=52: ρ=-0.282 (p=0.4991)
  - N=72: ρ=-0.577 (p=0.1340)
  - N=100: ρ=-0.386 (p=0.3453)
- **Lor4D**: ρ(N, ρ_cell) = +0.100 (p=0.873) — stable
  - N=20: ρ=-0.504 (p=0.2029)
  - N=36: ρ=+0.056 (p=0.8946)
  - N=52: ρ=+0.412 (p=0.3100)
  - N=72: ρ=-0.247 (p=0.5546)
  - N=100: ρ=-0.252 (p=0.5472)
- **Lor5D**: ρ(N, ρ_cell) = +0.600 (p=0.400) — weakening
  - N=20: ρ=-0.378 (p=0.3559)
  - N=52: ρ=-0.630 (p=0.0941)
  - N=72: ρ=+0.056 (p=0.8946)
  - N=100: ρ=-0.109 (p=0.7970)

## Q3: Within-Cell R²(logH ~ Σ_hist)

| family | N | R² | Pearson_r |
|--------|---|----| ---------|
| Lor2D | 20 | 0.816 | -0.904 |
| Lor3D | 20 | 0.183 | -0.427 |
| Lor4D | 20 | 0.127 | -0.356 |
| Lor5D | 20 | 0.281 | -0.530 |
| Lor2D | 36 | 0.432 | -0.658 |
| Lor3D | 36 | 0.144 | -0.379 |
| Lor4D | 36 | 0.046 | +0.214 |
| Lor2D | 52 | 0.137 | -0.371 |
| Lor3D | 52 | 0.074 | -0.272 |
| Lor4D | 52 | 0.135 | +0.367 |
| Lor5D | 52 | 0.490 | -0.700 |
| Lor2D | 72 | 0.111 | -0.333 |
| Lor3D | 72 | 0.470 | -0.686 |
| Lor4D | 72 | 0.021 | -0.147 |
| Lor5D | 72 | 0.046 | +0.216 |
| Lor2D | 100 | 0.927 | -0.963 |
| Lor3D | 100 | 0.178 | -0.422 |
| Lor4D | 100 | 0.019 | -0.138 |
| Lor5D | 100 | 0.111 | -0.333 |

## Q4: Σ_hist vs layer_count Relationship

| family | N | ρ(Σ_hist, layer_count) | p | direction |
|--------|---|------------------------|---|-----------|
| KR_like | 20 | — | — | degenerate |
| Lor2D | 20 | +1.000 | 0.0000 | positive |
| Lor3D | 20 | +1.000 | 0.0000 | positive |
| Lor4D | 20 | +1.000 | 0.0000 | positive |
| Lor5D | 20 | +1.000 | 0.0000 | positive |
| KR_like | 36 | — | — | degenerate |
| Lor2D | 36 | +1.000 | 0.0000 | positive |
| Lor3D | 36 | +1.000 | 0.0000 | positive |
| Lor4D | 36 | +1.000 | 0.0000 | positive |
| Lor5D | 36 | — | — | degenerate |
| KR_like | 52 | — | — | degenerate |
| Lor2D | 52 | +1.000 | 0.0000 | positive |
| Lor3D | 52 | +1.000 | 0.0000 | positive |
| Lor4D | 52 | +1.000 | 0.0000 | positive |
| Lor5D | 52 | +1.000 | 0.0000 | positive |
| KR_like | 72 | — | — | degenerate |
| Lor2D | 72 | +1.000 | 0.0000 | positive |
| Lor3D | 72 | +1.000 | 0.0000 | positive |
| Lor4D | 72 | +1.000 | 0.0000 | positive |
| Lor5D | 72 | +1.000 | 0.0000 | positive |
| KR_like | 100 | — | — | degenerate |
| Lor2D | 100 | +1.000 | 0.0000 | positive |
| Lor3D | 100 | +1.000 | 0.0000 | positive |
| Lor4D | 100 | +1.000 | 0.0000 | positive |
| Lor5D | 100 | +1.000 | 0.0000 | positive |

## Q5: Does Higher Σ_hist → Lower F7 Within Cells?

Since F7 = logH − λ·Σ_hist + ..., partial correlation ρ(F7, Σ_hist | other F7 terms) should be strongly negative.

| family | N | ρ(F7, Σ_hist) | ρ(logH_resid, Σ_hist) |
|--------|---|---------------|----------------------|
| Lor2D | 20 | -0.932 | -0.868 |
| Lor3D | 20 | +0.417 | -0.417 |
| Lor4D | 20 | -0.252 | -0.504 |
| Lor5D | 20 | -0.504 | -0.378 |
| Lor2D | 36 | -0.679 | -0.581 |
| Lor3D | 36 | -0.456 | -0.326 |
| Lor4D | 36 | +0.169 | +0.169 |
| Lor2D | 52 | -0.463 | -0.309 |
| Lor3D | 52 | -0.394 | -0.282 |
| Lor4D | 52 | +0.412 | -0.082 |
| Lor5D | 52 | -0.630 | -0.630 |
| Lor2D | 72 | -0.350 | -0.250 |
| Lor3D | 72 | -0.577 | -0.577 |
| Lor4D | 72 | -0.247 | -0.247 |
| Lor5D | 72 | +0.507 | -0.282 |
| Lor2D | 100 | -0.970 | -0.970 |
| Lor3D | 100 | -0.386 | -0.309 |
| Lor4D | 100 | -0.126 | -0.252 |
| Lor5D | 100 | -0.109 | -0.109 |

## Q6: Cross-Family Pooled Analysis (per N)

Pooled across families but fixed N — tests whether Σ_hist carries cross-family information.

| N | n | ρ(Σ_hist, logH) | ρ(Σ_hist, F7) | sig |
|---|---|-----------------|---------------|-----|
| 20 | 40 | -0.805 | -0.088 | ★★★ |
| 36 | 40 | -0.682 | -0.652 | ★★★ |
| 52 | 40 | -0.769 | -0.787 | ★★★ |
| 72 | 40 | -0.743 | -0.732 | ★★★ |
| 100 | 40 | -0.742 | -0.739 | ★★★ |

## Verdict

**MODERATE**: 84% cells show correct direction, 11% significant — partial support.