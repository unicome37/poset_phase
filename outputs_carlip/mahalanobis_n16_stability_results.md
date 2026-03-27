# N=16 Stability Test — Mahalanobis Variants

Date: 2026-03-27 11:04
N: 16
Families: 25
Seeds tested: 10
Reps per family per seed: 80

## 1. Winner Census (by seed)

### Baseline (full)

| Winner | Count |
|---|---:|
| Lor4D | 10 |

### Anchored (full)

| Winner | Count |
|---|---:|
| Lor4D | 9 |
| Lor5D | 1 |

### Anchored (diag)

| Winner | Count |
|---|---:|
| Lor4D | 10 |

## 2. Lor4D Rank and Margin Summary

| Method | Mean rank | #1 rate | Mean margin |
|---|---:|---:|---:|
| baseline_full | 1.00 | 100% | 1.3753 |
| anchored_full | 1.10 | 90% | 1.2475 |
| anchored_diag | 1.00 | 100% | 1.7275 |

## 3. Per-Seed Detail

| seed | mu_d | mu_c | mu_w | winner baseline | winner anchor_full | winner anchor_diag | Lor4D rank baseline | rank anchor_full | rank anchor_diag |
|---:|---:|---:|---:|---|---|---|---:|---:|---:|
| 42 | 3.9873 | 0.0738 | 0.5680 | Lor4D | Lor4D | Lor4D | 1 | 1 | 1 |
| 137 | 3.9748 | 0.1021 | 0.5469 | Lor4D | Lor4D | Lor4D | 1 | 1 | 1 |
| 271 | 3.9841 | 0.0915 | 0.5500 | Lor4D | Lor4D | Lor4D | 1 | 1 | 1 |
| 500 | 3.9991 | 0.0974 | 0.5477 | Lor4D | Lor4D | Lor4D | 1 | 1 | 1 |
| 777 | 3.9579 | 0.0836 | 0.5648 | Lor4D | Lor4D | Lor4D | 1 | 1 | 1 |
| 1001 | 3.9335 | 0.0868 | 0.5523 | Lor4D | Lor4D | Lor4D | 1 | 1 | 1 |
| 2023 | 3.9479 | 0.0801 | 0.5586 | Lor4D | Lor4D | Lor4D | 1 | 1 | 1 |
| 3141 | 3.9529 | 0.0863 | 0.5664 | Lor4D | Lor5D | Lor4D | 1 | 2 | 1 |
| 5000 | 4.0254 | 0.1046 | 0.5563 | Lor4D | Lor4D | Lor4D | 1 | 1 | 1 |
| 8888 | 3.9660 | 0.0804 | 0.5594 | Lor4D | Lor4D | Lor4D | 1 | 1 | 1 |

Runtime: 62.2s
