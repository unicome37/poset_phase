# Conjecture E bridge candidate compare

- Input: `outputs_unified_functional/bd_actions.csv`
- Output CSV: `outputs_unified_functional/conjecture_e_bridge_candidate_compare.csv`

## Residual bridge comparison

Baseline model: `F5 ~ N + family` with R² = 0.9692

| metric | R² ext | ΔR² | coef | corr(metric, residual F5) |
|---|---:|---:|---:|---:|
| bd_ratio | 0.9869 | +0.0178 | -18.5279 | -0.2786 |
| bdg_d2_corrected_norm | 0.9705 | +0.0013 | -1.3998 | -0.1079 |

## Alpha scan

Scoring convention: `score = -1.0·F5 + alpha·metric`

| metric | scope | best alpha | best win rate | mean gap (control - focus) |
|---|---|---:|---:|---:|
| bd_ratio | ALL | +0.00 | 1.000 | +29.0512 |
| bd_ratio | 16 | +0.00 | 1.000 | +29.7487 |
| bd_ratio | 20 | +0.00 | 1.000 | +26.7292 |
| bd_ratio | 28 | +0.00 | 1.000 | +29.7380 |
| bd_ratio | 36 | +0.00 | 1.000 | +29.9891 |
| bdg_d2_corrected_norm | ALL | +0.00 | 1.000 | +29.0512 |
| bdg_d2_corrected_norm | 16 | +0.00 | 1.000 | +29.7487 |
| bdg_d2_corrected_norm | 20 | +0.00 | 1.000 | +26.7292 |
| bdg_d2_corrected_norm | 28 | +0.00 | 1.000 | +29.7380 |
| bdg_d2_corrected_norm | 36 | +0.00 | 1.000 | +29.9891 |

## Takeaway

The main question is whether any candidate both preserves a nonzero residual bridge and can improve the Lor4D-vs-KR_like ordering at finite alpha. If one metric keeps the larger ΔR² while both best alphas stay pinned near zero, it is the better bridge proxy for the current stage.
