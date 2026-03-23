# Conjecture E: ρ_BD vs R(p) Curvature Cross-Check

## Question

Does the BD action density (ρ_BD, Layer 3 of Conjecture E) co-vary with
the discrete curvature proxy (bdg_d4, Layer 2) — both globally across
poset instances and locally within individual posets?

## Method

**Strategy A (cross-poset):** For each poset instance, compute global
bd_ratio and global bdg_d4. Spearman-correlate across instances.

**Strategy B (within-poset):** For deep-layer families (Lor2D, TransPerc),
use a sliding window (center ± 1 layer) to compute local observables.
Correlate across windows within each poset.

Note: High-d families (Lor4D/5D, KR_like) have only ~3 layers, with
boundary layers having bd_ratio=0. Within-poset correlation is impossible.

## Strategy A: Cross-poset results

### Within-family Spearman ρ(bd_ratio, bdg_d4)

| family | ρ | p-value | n | mean bd | mean d4 |
|---|---:|---:|---:|---:|---:|
| Lor2D | +0.133 | 0.2381 | 80 | 1.2147 | +0.939 |
| Lor3D | +0.105 | 0.3544 | 80 | 0.6526 | +1.339 |
| Lor4D | +0.245 | 0.0286 | 80 | 0.2100 | +0.517 |
| Lor5D | +0.290 | 0.0092 | 80 | 0.0633 | +0.029 |
| KR_like | +0.370 | 0.0007 | 80 | 0.4405 | -3.031 |
| TransPerc | -0.060 | 0.5987 | 80 | 0.9035 | +1.797 |

### Cross-family pooled: ρ(bd_ratio, bdg_d4) = +0.242 (p = 7.56e-08, n = 480)
### Cross-family pooled: ρ(bd_ratio, link_frac) = -0.990 (p = 0.00e+00)

## Strategy B: Within-poset results

| family | N | mean ρ(bd,d4) | mean ρ(bd,lf) | n_posets |
|---|---:|---:|---:|---:|
| Lor2D | 36 | +0.320 | -0.977 | 20 |
| Lor2D | 48 | +0.354 | -0.969 | 20 |
| TransPerc | 36 | +0.690 | -0.999 | 18 |
| TransPerc | 48 | +0.781 | -0.997 | 20 |

## Interpretation

- Positive ρ(bd_ratio, bdg_d4): regions with richer interval structure
  also have higher BDG d=4 action — the BD-to-EH bridge is locally consistent.
- Negative ρ(bd_ratio, link_frac): more non-trivial intervals (fewer links)
  → richer BD structure — bd_ratio measures 'causal thickness'.
- Cross-family correlations reflect that different families have
  systematically different bd_ratio and curvature profiles.
