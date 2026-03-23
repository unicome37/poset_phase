# Conjecture E bd_ratio fine scan

- Input: `outputs_unified_functional/bd_actions.csv`
- Output CSV: `outputs_unified_functional/conjecture_e_bd_ratio_scan.csv`

## Residual bridge

Baseline model: `F5 ~ N + family` with R² = 0.9692
Extended model: `F5 ~ N + family + bd_ratio` with R² = 0.9869
ΔR² = +0.0178
Fitted bd_ratio coefficient (raw units) = -18.5279
Corr(bd_ratio, residual F5) = -0.2786
Corr(z(bd_ratio), residual F5) = -0.2786

## Alpha scan

Scoring convention: `score = -1.0·F5 + alpha·bd_ratio`

| scope | best alpha | best win rate | mean gap (control - focus) |
|---|---:|---:|---:|
| ALL | +0.00 | 1.000 | +29.0512 |
| 16 | +0.00 | 1.000 | +29.7487 |
| 20 | +0.00 | 1.000 | +26.7292 |
| 28 | +0.00 | 1.000 | +29.7380 |
| 36 | +0.00 | 1.000 | +29.9891 |

## Interpretation

This scan treats `bd_ratio` as the strongest bridge proxy and checks whether a finite `alpha` can actually improve the current Lor4D-vs-KR_like ordering. If the best alpha stays near zero while the residual R² gain stays positive, that is a good sign that `bd_ratio` is a residual-corrective bridge rather than a rank-flipping replacement for F5.
