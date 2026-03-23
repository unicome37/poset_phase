# Enhanced E2: Causal Depth Directional Growth

**Date:** 2026-03-23

**Design:** 5 families × 4 N × 4 T × 30 reps = 2400 experiments

**Families:** Lor2D, Lor3D, Lor4D, KR_like, TransPerc

**N values:** [16, 24, 36, 48]

**T values:** [4, 8, 16, 32]


## Metric 1: Layer Count (original E2 metric)

### Per-family (all T pooled)

| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | n |
|--------|-----------|-----------|-----------|-----------|------------|---|
| Lor2D | 9.96 | 9.97 | -0.019 | 38.3% | 9.187e-01 | 480 |
| Lor3D | 10.11 | 10.16 | -0.054 | 36.0% | 5.182e-01 | 480 |
| Lor4D | 10.02 | 10.13 | -0.115 | 34.4% | 2.498e-01 | 480 |
| KR_like | 10.31 | 10.40 | -0.085 | 33.5% | 1.541e-01 | 480 |
| TransPerc | 10.08 | 10.07 | +0.013 | 36.5% | 9.731e-01 | 480 |

## Metric 2: Mean Element Depth (continuous-valued)

### Per-family (all T pooled)

| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | n |
|--------|-----------|-----------|-----------|-----------|------------|---|
| Lor2D | 3.263 | 6.879 | -3.6167 | 1.0% | 3.777e-80 | 480 |
| Lor3D | 2.808 | 7.368 | -4.5598 | 0.0% | 2.333e-80 | 480 |
| Lor4D | 2.617 | 7.466 | -4.8481 | 0.0% | 2.333e-80 | 480 |
| KR_like | 2.589 | 7.628 | -5.0390 | 0.0% | 2.333e-80 | 480 |
| TransPerc | 3.049 | 7.332 | -4.2833 | 0.0% | 2.333e-80 | 480 |

### Per-T breakdown (mean element depth, all families pooled)

| T | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | n |
|---|-----------|-----------|-----------|-----------|------------|---|
| 4 | 0.502 | 2.320 | -1.8188 | 0.2% | 6.246e-100 | 600 |
| 8 | 1.193 | 4.296 | -3.1033 | 0.5% | 6.522e-100 | 600 |
| 16 | 2.860 | 8.046 | -5.1855 | 0.2% | 6.035e-100 | 600 |
| 32 | 6.906 | 14.676 | -7.7699 | 0.0% | 5.978e-100 | 600 |

### Per-N breakdown (mean element depth, all families pooled)

| N | mean diff | fwd>bwd % | Wilcoxon p | n |
|---|-----------|-----------|------------|---|
| 16 | -3.0343 | 0.8% | 9.002e-100 | 600 |
| 24 | -4.0646 | 0.0% | 5.975e-100 | 600 |
| 36 | -5.0722 | 0.0% | 5.977e-100 | 600 |
| 48 | -5.7064 | 0.0% | 5.977e-100 | 600 |

### Family × T (mean diff in mean_depth)

| Family | T=4 | T=8 | T=16 | T=32 |
|--------|---:|---:|---:|---:|
| Lor2D | -1.454 | -2.478 | -4.234 | -6.301 |
| Lor3D | -1.884 | -3.176 | -5.246 | -7.933 |
| Lor4D | -1.970 | -3.349 | -5.601 | -8.473 |
| KR_like | -2.112 | -3.618 | -5.782 | -8.643 |
| TransPerc | -1.674 | -2.896 | -5.063 | -7.500 |

### Family × T (fwd > bwd rate for mean_depth)

| Family | T=4 | T=8 | T=16 | T=32 |
|--------|---:|---:|---:|---:|
| Lor2D | 1% | 2% | 1% | 0% |
| Lor3D | 0% | 0% | 0% | 0% |
| Lor4D | 0% | 0% | 0% | 0% |
| KR_like | 0% | 0% | 0% | 0% |
| TransPerc | 0% | 0% | 0% | 0% |

## Global Summary

| Group | mean diff | fwd>bwd % | Wilcoxon p | n |
|-------|-----------|-----------|------------|---|
| Lorentzian (2D/3D/4D) | -4.3415 | 0.3% | 7.620e-237 | 1440 |
| KR_like (control) | -5.0390 | 0.0% | 2.333e-80 | 480 |
| TransPerc | -4.2833 | 0.0% | 2.333e-80 | 480 |

## Signal Scaling with T

Does the directional asymmetry grow with augmentation steps?

| Family | Spearman(T, mean_diff) | p-value |
|--------|----------------------|---------|
| Lor2D | -1.000 | 0.000e+00 |
| Lor3D | -1.000 | 0.000e+00 |
| Lor4D | -1.000 | 0.000e+00 |
| KR_like | -1.000 | 0.000e+00 |
| TransPerc | -1.000 | 0.000e+00 |

## Conclusion

(Filled based on results above)
