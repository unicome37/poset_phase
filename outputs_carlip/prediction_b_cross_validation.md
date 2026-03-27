# Prediction B — Cross-Validation Test

**Goal**: Verify Lor4D's #1 ranking is NOT due to overfitting.

**Design**: 5-fold CV, 40 reps, 3 seeds.

- **Oracle mode**: μ, Σ from ALL Lor4D samples; score on same set
- **CV mode**: μ, Σ from 80% train; score on 20% test (Lor4D test only)

## 1. Oracle vs CV: Lor4D Rank

| Seed | N | Oracle LSD | Oracle Mahal | CV LSD (5 folds) | CV Mahal (5 folds) |
|------|---|:----------:|:-----------:|:----------------:|:-----------------:|
| 42 | 16 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,2] |
| 42 | 20 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 28 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 36 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 48 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 64 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 96 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 128 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 16 | #1 | #1 | [1,1,1,2,1] | [1,2,3,2,1] |
| 777 | 20 | #1 | #1 | [1,1,1,1,1] | [2,1,1,1,1] |
| 777 | 28 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 36 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 48 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 64 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 96 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 128 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 3141 | 16 | #1 | #1 | [1,1,1,2,2] | [1,1,1,3,1] |
| 3141 | 20 | #1 | #1 | [1,1,1,1,1] | [1,2,1,1,1] |
| 3141 | 28 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 3141 | 36 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 3141 | 48 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 3141 | 64 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 3141 | 96 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 3141 | 128 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |

## 2. Summary

| Mode | LSD-Well #1 rate | Mahalanobis #1 rate |
|------|:----------------:|:-------------------:|
| Oracle | 24/24 (100%) | 24/24 (100%) |
| CV (5-fold) | 117/120 (98%) | 113/120 (94%) |

## 3. Margin Degradation: Oracle vs CV

| N | Oracle LSD margin | CV LSD margin (mean±std) | Oracle Mahal margin | CV Mahal margin (mean±std) |
|---|:-:|:-:|:-:|:-:|
| 16 | +0.056 | +0.065±0.040 | +1.8 | +1.7±1.0 |
| 20 | +0.118 | +0.116±0.039 | +2.6 | +2.2±1.0 |
| 28 | +0.133 | +0.132±0.025 | +5.5 | +5.0±1.4 |
| 36 | +0.131 | +0.131±0.018 | +8.2 | +8.0±1.7 |
| 48 | +0.123 | +0.122±0.015 | +10.6 | +10.1±4.0 |
| 64 | +0.125 | +0.125±0.026 | +14.7 | +14.5±5.8 |
| 96 | +0.140 | +0.140±0.006 | +23.4 | +23.3±2.7 |
| 128 | +0.138 | +0.138±0.007 | +30.1 | +30.3±7.3 |

## 4. Overfitting Diagnosis

**LSD-Well**: CV #1 rate = 98%. Partial degradation from oracle.

**Mahalanobis**: CV #1 rate = 94%. Minor degradation at small N expected.

Mean margin retention (CV/Oracle):
- LSD-Well: 105.7%
- Mahalanobis: 95.1%

## 5. Conclusion

Train/test separation confirms that Lor4D's dominance is a genuine
structural signal, not an artifact of using the same data to estimate
centroids and to score. The margin retention shows the scoring
generalizes from training to unseen test data.