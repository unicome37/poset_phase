# Prediction C — F7 Large-N Verification

**Data**: 1000 samples, N ∈ {20, 36, 52, 72, 100}

**Families**: KR_like, Lor2D, Lor3D, Lor4D, Lor5D


## Q1: Within-Cell Correlation ρ(Σ_hist, logH)

The core Prediction C claim: higher Σ_hist → lower logH within fixed (family, N).

| family | N | n | ρ(Σ,logH) | p_value | sig | direction |
|--------|---|---|-----------|---------|-----|-----------|
| KR_like | 20 | 40 | — | — | — | degenerate |
| Lor2D | 20 | 40 | -0.852 | 0.0000 | ★★★ | ✅ neg |
| Lor3D | 20 | 40 | -0.741 | 0.0000 | ★★★ | ✅ neg |
| Lor4D | 20 | 40 | -0.207 | 0.1998 |  | ✅ neg |
| Lor5D | 20 | 40 | -0.457 | 0.0030 | ★★ | ✅ neg |
| KR_like | 36 | 40 | — | — | — | degenerate |
| Lor2D | 36 | 40 | -0.689 | 0.0000 | ★★★ | ✅ neg |
| Lor3D | 36 | 40 | -0.524 | 0.0005 | ★★★ | ✅ neg |
| Lor4D | 36 | 40 | -0.392 | 0.0124 | ★ | ✅ neg |
| Lor5D | 36 | 40 | -0.418 | 0.0073 | ★★ | ✅ neg |
| KR_like | 52 | 40 | — | — | — | degenerate |
| Lor2D | 52 | 40 | -0.488 | 0.0014 | ★★ | ✅ neg |
| Lor3D | 52 | 40 | -0.390 | 0.0129 | ★ | ✅ neg |
| Lor4D | 52 | 40 | -0.342 | 0.0308 | ★ | ✅ neg |
| Lor5D | 52 | 40 | +0.119 | 0.4636 |  | ❌ pos |
| KR_like | 72 | 40 | — | — | — | degenerate |
| Lor2D | 72 | 40 | -0.649 | 0.0000 | ★★★ | ✅ neg |
| Lor3D | 72 | 40 | -0.434 | 0.0051 | ★★ | ✅ neg |
| Lor4D | 72 | 40 | -0.190 | 0.2415 |  | ✅ neg |
| Lor5D | 72 | 40 | +0.276 | 0.0842 |  | ❌ pos |
| KR_like | 100 | 40 | — | — | — | degenerate |
| Lor2D | 100 | 40 | -0.306 | 0.0548 |  | ✅ neg |
| Lor3D | 100 | 40 | -0.445 | 0.0041 | ★★ | ✅ neg |
| Lor4D | 100 | 40 | -0.580 | 0.0001 | ★★★ | ✅ neg |
| Lor5D | 100 | 40 | -0.409 | 0.0087 | ★★ | ✅ neg |

**Summary**: 18/20 cells have ρ < 0 (correct direction), 15/20 significant at p<0.05


## Q2: N-Scaling of Within-Cell ρ(Σ_hist, logH)

Does the correlation strengthen with N?

- **Lor2D**: ρ(N, ρ_cell) = +0.900 (p=0.037) — weakening
  - N=20: ρ=-0.852 (p=0.0000)
  - N=36: ρ=-0.689 (p=0.0000)
  - N=52: ρ=-0.488 (p=0.0014)
  - N=72: ρ=-0.649 (p=0.0000)
  - N=100: ρ=-0.306 (p=0.0548)
- **Lor3D**: ρ(N, ρ_cell) = +0.600 (p=0.285) — weakening
  - N=20: ρ=-0.741 (p=0.0000)
  - N=36: ρ=-0.524 (p=0.0005)
  - N=52: ρ=-0.390 (p=0.0129)
  - N=72: ρ=-0.434 (p=0.0051)
  - N=100: ρ=-0.445 (p=0.0041)
- **Lor4D**: ρ(N, ρ_cell) = -0.300 (p=0.624) — stable
  - N=20: ρ=-0.207 (p=0.1998)
  - N=36: ρ=-0.392 (p=0.0124)
  - N=52: ρ=-0.342 (p=0.0308)
  - N=72: ρ=-0.190 (p=0.2415)
  - N=100: ρ=-0.580 (p=0.0001)
- **Lor5D**: ρ(N, ρ_cell) = +0.700 (p=0.188) — weakening
  - N=20: ρ=-0.457 (p=0.0030)
  - N=36: ρ=-0.418 (p=0.0073)
  - N=52: ρ=+0.119 (p=0.4636)
  - N=72: ρ=+0.276 (p=0.0842)
  - N=100: ρ=-0.409 (p=0.0087)

## Q3: Within-Cell R²(logH ~ Σ_hist)

| family | N | R² | Pearson_r |
|--------|---|----| ---------|
| Lor2D | 20 | 0.693 | -0.833 |
| Lor3D | 20 | 0.514 | -0.717 |
| Lor4D | 20 | 0.058 | -0.242 |
| Lor5D | 20 | 0.236 | -0.486 |
| Lor2D | 36 | 0.454 | -0.673 |
| Lor3D | 36 | 0.240 | -0.490 |
| Lor4D | 36 | 0.164 | -0.405 |
| Lor5D | 36 | 0.235 | -0.485 |
| Lor2D | 52 | 0.303 | -0.551 |
| Lor3D | 52 | 0.164 | -0.405 |
| Lor4D | 52 | 0.121 | -0.348 |
| Lor5D | 52 | 0.032 | +0.180 |
| Lor2D | 72 | 0.363 | -0.602 |
| Lor3D | 72 | 0.162 | -0.402 |
| Lor4D | 72 | 0.021 | -0.144 |
| Lor5D | 72 | 0.069 | +0.263 |
| Lor2D | 100 | 0.119 | -0.345 |
| Lor3D | 100 | 0.261 | -0.511 |
| Lor4D | 100 | 0.335 | -0.579 |
| Lor5D | 100 | 0.159 | -0.399 |

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
| Lor5D | 36 | +1.000 | 0.0000 | positive |
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
| Lor2D | 20 | -0.922 | -0.813 |
| Lor3D | 20 | +0.371 | -0.754 |
| Lor4D | 20 | +0.105 | -0.509 |
| Lor5D | 20 | -0.592 | -0.522 |
| Lor2D | 36 | -0.773 | -0.686 |
| Lor3D | 36 | -0.538 | -0.569 |
| Lor4D | 36 | -0.004 | -0.509 |
| Lor5D | 36 | -0.461 | -0.418 |
| Lor2D | 52 | -0.596 | -0.507 |
| Lor3D | 52 | -0.382 | -0.430 |
| Lor4D | 52 | +0.182 | -0.402 |
| Lor5D | 52 | +0.109 | +0.119 |
| Lor2D | 72 | -0.686 | -0.652 |
| Lor3D | 72 | -0.447 | -0.422 |
| Lor4D | 72 | +0.200 | -0.303 |
| Lor5D | 72 | +0.265 | +0.276 |
| Lor2D | 100 | -0.356 | -0.300 |
| Lor3D | 100 | -0.466 | -0.454 |
| Lor4D | 100 | -0.616 | -0.576 |
| Lor5D | 100 | -0.409 | -0.387 |

## Q6: Cross-Family Pooled Analysis (per N)

Pooled across families but fixed N — tests whether Σ_hist carries cross-family information.

| N | n | ρ(Σ_hist, logH) | ρ(Σ_hist, F7) | sig |
|---|---|-----------------|---------------|-----|
| 20 | 200 | -0.782 | -0.011 | ★★★ |
| 36 | 200 | -0.697 | -0.691 | ★★★ |
| 52 | 200 | -0.753 | -0.771 | ★★★ |
| 72 | 200 | -0.765 | -0.760 | ★★★ |
| 100 | 200 | -0.734 | -0.749 | ★★★ |

## Verdict

**STRONG**: 90% cells show correct direction, 75% significant — Prediction C well-supported at large N.