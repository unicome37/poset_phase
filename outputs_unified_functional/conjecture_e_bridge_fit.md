# Conjecture E bridge fit

- Input: `outputs_unified_functional/bd_actions.csv`
- Output CSV: `outputs_unified_functional/conjecture_e_bridge_fit.csv`

## Observable correlations

| metric | Spearman(F5) | Spearman(log_H) | mean focus | mean control |
|---|---:|---:|---:|---:|
| bdg_d2_link_norm | -0.2912 | -0.4551 | -2.6511 | -5.2407 |
| bdg_d2_corrected_norm | -0.3746 | -0.4320 | -2.2134 | -4.8401 |
| bdg_d4_standard_norm | -0.0348 | -0.0446 | 0.4865 | -1.8676 |
| bd_d4_trunc_norm | -0.1886 | -0.3425 | -0.3303 | -1.7989 |
| bd_ratio | -0.3395 | -0.1429 | 0.1452 | 0.6681 |

## Best lambda per observable

| metric | best lambda | best win rate |
|---|---:|---:|
| bdg_d2_link_norm | 0 | 1.000 |
| bdg_d2_corrected_norm | 0 | 1.000 |
| bdg_d4_standard_norm | 0 | 1.000 |
| bd_d4_trunc_norm | 0 | 1.000 |
| bd_ratio | 0 | 1.000 |

## Residual bridge

Baseline model: `F5 ~ N + family` with R² = 0.9692

| metric | R² ext | ΔR² | Spearman(metric, residual F5) |
|---|---:|---:|---:|
| bdg_d2_link_norm | 0.9693 | +0.0001 | +0.0397 |
| bdg_d2_corrected_norm | 0.9705 | +0.0013 | -0.0971 |
| bdg_d4_standard_norm | 0.9694 | +0.0002 | -0.0865 |
| bd_d4_trunc_norm | 0.9693 | +0.0001 | +0.0435 |
| bd_ratio | 0.9869 | +0.0178 | -0.0694 |

## Interpretation

The bridge test asks whether a BD-like scalar can improve Lorentzian-vs-KR ordering when added to the current F5 baseline. A high best win rate with O(1) lambda is the first finite-size signal that the bridge is not just a tiny perturbation.
