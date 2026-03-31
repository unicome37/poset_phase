# Prediction B — Cross-Validation Test

**Goal**: Verify Lor4D's #1 ranking is NOT due to overfitting.

**Design**: 5-fold CV, 40 reps, 3 seeds.

- **Oracle mode**: μ, Σ from ALL Lor4D samples; score on same set
- **CV mode**: μ, Σ from 80% train; score on 20% test (Lor4D test only)

## 1. Oracle vs CV: Lor4D Rank

| Seed | N | Oracle LSD | Oracle Mahal | CV LSD (5 folds) | CV Mahal (5 folds) |
|------|---|:----------:|:-----------:|:----------------:|:-----------------:|
| 42 | 16 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 20 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,2] |
| 42 | 28 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 36 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 48 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 64 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 96 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 42 | 128 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 16 | #1 | #1 | [1,1,1,1,1] | [3,1,1,1,1] |
| 777 | 20 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 28 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 36 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 48 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 64 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 96 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 777 | 128 | #1 | #1 | [1,1,1,1,1] | [1,1,1,1,1] |
| 3141 | 16 | #1 | #1 | [1,1,1,1,1] | [2,1,1,1,1] |
| 3141 | 20 | #1 | #1 | [1,1,1,1,1] | [3,1,1,1,1] |
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
| CV (5-fold) | 120/120 (100%) | 116/120 (97%) |

## 3. Margin Degradation: Oracle vs CV

| N | Oracle LSD margin | CV LSD margin (mean±std) | Oracle Mahal margin | CV Mahal margin (mean±std) |
|---|:-:|:-:|:-:|:-:|
| 16 | +0.095 | +0.092±0.033 | +3.3 | +3.0±1.4 |
| 20 | +0.096 | +0.093±0.032 | +1.9 | +2.0±1.0 |
| 28 | +0.115 | +0.114±0.019 | +3.4 | +3.0±1.2 |
| 36 | +0.129 | +0.128±0.022 | +5.7 | +5.0±1.9 |
| 48 | +0.114 | +0.113±0.015 | +8.3 | +8.2±2.9 |
| 64 | +0.123 | +0.123±0.008 | +11.5 | +11.3±2.3 |
| 96 | +0.133 | +0.133±0.010 | +22.4 | +22.6±5.3 |
| 128 | +0.133 | +0.133±0.009 | +32.0 | +32.3±6.1 |

## 4. Overfitting Diagnosis

**LSD-Well**: Zero overfitting. CV #1 rate = 100%. ✅

**Mahalanobis**: CV #1 rate = 97%. Minor degradation at small N expected.

Mean margin retention (CV/Oracle):
- LSD-Well: 98.8%
- Mahalanobis: 97.5%

## 5. Conclusion

Train/test separation confirms that Lor4D's dominance is a genuine
structural signal, not an artifact of using the same data to estimate
centroids and to score. The margin retention shows the scoring
generalizes from training to unseen test data.