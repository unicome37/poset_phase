# E2 v3: Direction-Agnostic Structural Complexity

**Date:** 2026-03-23

**Design:** 5 families × 3 N × 3 T × 30 reps = 1350


## Reformulated Prediction

E2 (original): "Forward augmentation deepens layer structure faster" — WEAK (layer count is symmetric)

E2 (v3): "Forward augmentation preserves/enhances structural complexity better than backward"

- Height ratio h = n_layers/N: forward should maintain depth-per-element

- Width entropy S_w: forward should produce richer, less concentrated layer profiles

- Order fraction R: forward should preserve causal density

- Σ_hist: forward should increase historical sedimentation


## Height Ratio (n_layers / N)

**Prediction:** fwd deepens more

### Per-family (all T pooled)

| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | verdict |
|--------|-----------|-----------|-----------|-----------|------------|---------|
| Lor2D | +0.1013 | +0.1018 | -0.0005 | 36.3% | 7.666e-01 | ❌ |
| Lor3D | +0.1331 | +0.1316 | +0.0016 | 37.8% | 6.966e-01 | ❌ |
| Lor4D | +0.1466 | +0.1469 | -0.0004 | 40.0% | 7.381e-01 | ❌ |
| KR_like | +0.1510 | +0.1476 | +0.0034 | 41.1% | 1.438e-01 | ❌ |
| TransPerc | +0.1264 | +0.1257 | +0.0007 | 38.5% | 8.147e-01 | ❌ |

### Per-T (all families pooled)

| T | mean diff | fwd>bwd % | Wilcoxon p |
|---|-----------|-----------|------------|
| 4 | +0.0016 | 35.1% | 4.748e-01 |
| 8 | +0.0004 | 38.7% | 8.357e-01 |
| 16 | +0.0010 | 42.4% | 5.170e-01 |

### Family × T (fwd>bwd rate)

| Family | T=4 | T=8 | T=16 |
|--------|---:|---:|---:|
| Lor2D | 33% | 39% | 37% |
| Lor3D | 31% | 37% | 46% |
| Lor4D | 37% | 41% | 42% |
| KR_like | 34% | 38% | 51% |
| TransPerc | 40% | 39% | 37% |

## Width Entropy S_w

**Prediction:** fwd produces richer layer structure

### Per-family (all T pooled)

| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | verdict |
|--------|-----------|-----------|-----------|-----------|------------|---------|
| Lor2D | +0.5544 | +0.5780 | -0.0235 | 35.6% | 1.270e-05 | ❌ |
| Lor3D | +0.6958 | +0.7701 | -0.0743 | 24.4% | 5.871e-22 | ❌ |
| Lor4D | +0.7788 | +0.9128 | -0.1341 | 15.9% | 8.237e-35 | ❌ |
| KR_like | +0.7636 | +0.8393 | -0.0758 | 21.5% | 1.684e-26 | ❌ |
| TransPerc | +0.6966 | +0.7960 | -0.0994 | 19.6% | 5.349e-29 | ❌ |

### Per-T (all families pooled)

| T | mean diff | fwd>bwd % | Wilcoxon p |
|---|-----------|-----------|------------|
| 4 | -0.0925 | 21.3% | 1.161e-44 |
| 8 | -0.0852 | 22.7% | 1.862e-37 |
| 16 | -0.0666 | 26.2% | 6.147e-31 |

### Family × T (fwd>bwd rate)

| Family | T=4 | T=8 | T=16 |
|--------|---:|---:|---:|
| Lor2D | 33% | 36% | 38% |
| Lor3D | 22% | 22% | 29% |
| Lor4D | 14% | 13% | 20% |
| KR_like | 19% | 18% | 28% |
| TransPerc | 18% | 24% | 17% |

## Order Fraction R

**Prediction:** fwd preserves causal density better

### Per-family (all T pooled)

| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | verdict |
|--------|-----------|-----------|-----------|-----------|------------|---------|
| Lor2D | +0.1115 | +0.1158 | -0.0044 | 34.8% | 3.754e-07 | ❌ |
| Lor3D | +0.1397 | +0.1595 | -0.0198 | 14.1% | 8.588e-35 | ❌ |
| Lor4D | +0.1465 | +0.1784 | -0.0319 | 10.7% | 4.748e-39 | ❌ |
| KR_like | +0.1385 | +0.1382 | +0.0003 | 47.0% | 7.382e-01 | ❌ |
| TransPerc | +0.1297 | +0.1866 | -0.0570 | 1.1% | 9.986e-46 | ❌ |

### Per-T (all families pooled)

| T | mean diff | fwd>bwd % | Wilcoxon p |
|---|-----------|-----------|------------|
| 4 | -0.0132 | 25.1% | 2.508e-33 |
| 8 | -0.0222 | 20.7% | 1.086e-43 |
| 16 | -0.0322 | 18.9% | 3.154e-52 |

### Family × T (fwd>bwd rate)

| Family | T=4 | T=8 | T=16 |
|--------|---:|---:|---:|
| Lor2D | 34% | 39% | 31% |
| Lor3D | 22% | 14% | 6% |
| Lor4D | 17% | 10% | 6% |
| KR_like | 49% | 40% | 52% |
| TransPerc | 3% | 0% | 0% |

## Σ_hist (Historical Sedimentation)

**Prediction:** fwd produces deeper sedimentation

### Per-family (all T pooled)

| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | verdict |
|--------|-----------|-----------|-----------|-----------|------------|---------|
| Lor2D | +0.1520 | +0.1527 | -0.0007 | 36.3% | 7.726e-01 | ❌ |
| Lor3D | +0.1997 | +0.1973 | +0.0024 | 37.8% | 6.879e-01 | ❌ |
| Lor4D | +0.2199 | +0.2204 | -0.0005 | 40.0% | 7.395e-01 | ❌ |
| KR_like | +0.2265 | +0.2214 | +0.0051 | 41.1% | 1.420e-01 | ❌ |
| TransPerc | +0.1896 | +0.1885 | +0.0011 | 38.5% | 8.028e-01 | ❌ |

### Per-T (all families pooled)

| T | mean diff | fwd>bwd % | Wilcoxon p |
|---|-----------|-----------|------------|
| 4 | +0.0023 | 35.1% | 4.617e-01 |
| 8 | +0.0005 | 38.7% | 8.566e-01 |
| 16 | +0.0015 | 42.4% | 5.053e-01 |

### Family × T (fwd>bwd rate)

| Family | T=4 | T=8 | T=16 |
|--------|---:|---:|---:|
| Lor2D | 33% | 39% | 37% |
| Lor3D | 31% | 37% | 46% |
| Lor4D | 37% | 41% | 42% |
| KR_like | 34% | 38% | 51% |
| TransPerc | 40% | 39% | 37% |

## Max Width Fraction

**Prediction:** lower = deeper structure

### Per-family (all T pooled)

| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | verdict |
|--------|-----------|-----------|-----------|-----------|------------|---------|
| Lor2D | -0.0655 | -0.0732 | +0.0077 | 28.1% | 3.335e-08 | ❌ |
| Lor3D | -0.1053 | -0.1322 | +0.0269 | 51.9% | 4.292e-20 | ✅ |
| Lor4D | -0.1336 | -0.1881 | +0.0545 | 73.3% | 1.938e-31 | ✅ |
| KR_like | -0.1301 | -0.1570 | +0.0269 | 47.0% | 1.316e-22 | ❌ |
| TransPerc | -0.1260 | -0.1850 | +0.0589 | 80.4% | 7.153e-37 | ✅ |

### Per-T (all families pooled)

| T | mean diff | fwd>bwd % | Wilcoxon p |
|---|-----------|-----------|------------|
| 4 | +0.0381 | 53.6% | 1.398e-38 |
| 8 | +0.0365 | 56.9% | 3.707e-39 |
| 16 | +0.0305 | 58.0% | 2.755e-40 |

### Family × T (fwd>bwd rate)

| Family | T=4 | T=8 | T=16 |
|--------|---:|---:|---:|
| Lor2D | 26% | 27% | 32% |
| Lor3D | 52% | 48% | 56% |
| Lor4D | 69% | 77% | 74% |
| KR_like | 39% | 57% | 46% |
| TransPerc | 82% | 77% | 82% |

## Grand Summary

| Metric | Lorentzian fwd>bwd % | KR fwd>bwd % | TransPerc fwd>bwd % | Signal? |
|--------|---------------------|-------------|--------------------|---------|

| Height Ratio (n_layers / N) | 38.0% | 41.1% | 38.5% | ❌ |
| Width Entropy S_w | 25.3% | 21.5% | 19.6% | ❌ |
| Order Fraction R | 19.9% | 47.0% | 1.1% | ❌ |
| Σ_hist (Historical Sedimentation) | 38.0% | 41.1% | 38.5% | ❌ |
| Max Width Fraction | 51.1% | 47.0% | 80.4% | ⚠️ |

## Interpretation

If forward augmentation consistently produces HIGHER structural complexity

(height ratio, width entropy, Σ_hist) than backward, this confirms E2:

the causal future direction is structurally richer than the past direction.

This is complementary to E1 (entropy asymmetry) and strengthens the

"time arrow from structural asymmetry" thesis.
