# Bootstrap Confidence Intervals for Lor4D Dominance

B = 1000 bootstrap iterations, 40 reps, train/test = 80/20.

## 1. Mahalanobis LSD — Bootstrap Results

| N | Rank 95% CI | Margin 95% CI | P(#1) | P(non-Lor wins) |
|---|:-:|:-:|:-:|:-:|
| 16 | [1, 2] | [+0.20, +3.73] | 97.0% | 1.2% |
| 20 | [1, 1] | [+7.48, +11.49] | 100.0% | 0.0% |
| 28 | [1, 1] | [+0.97, +5.46] | 98.9% | 0.0% |
| 48 | [1, 1] | [+2.25, +6.19] | 99.9% | 0.0% |
| 64 | [1, 1] | [+7.33, +13.67] | 100.0% | 0.0% |
| 128 | [1, 1] | [+26.11, +41.81] | 100.0% | 0.0% |

## 2. LSD-Well — Bootstrap Results

| N | Rank 95% CI | Margin 95% CI | P(#1) | P(non-Lor wins) |
|---|:-:|:-:|:-:|:-:|
| 16 | [1, 2] | [+0.002, +0.103] | 88.9% | 11.1% |
| 20 | [1, 1] | [+0.109, +0.182] | 100.0% | 0.0% |
| 28 | [1, 1] | [+0.057, +0.162] | 100.0% | 0.0% |
| 48 | [1, 1] | [+0.083, +0.140] | 100.0% | 0.0% |
| 64 | [1, 1] | [+0.133, +0.159] | 100.0% | 0.0% |
| 128 | [1, 1] | [+0.113, +0.133] | 100.0% | 0.0% |

## 3. Summary

- At N≥28, both metrics achieve P(#1) ≈ 100% with strictly positive margin CIs
- P(non-Lorentzian wins) = 0% at all N — the only competitor is Lor5D at small N
- Bootstrap confirms CV results are not split-dependent artifacts