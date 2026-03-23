# Prediction E: Time Arrow from Structural Asymmetry

## E1: Augmentation Entropy Asymmetry

A_entropy = Δlog H(forward) - Δlog H(backward)

Prediction: A_entropy > 0 (forward augmentation increases entropy more)

### Per-family results

| family | mean A | std | A>0 % | n |
|---|---:|---:|---:|---:|
| Lor2D | +0.1597 | 0.8123 | 58.9% | 180 |
| Lor3D | +0.2721 | 0.8972 | 61.1% | 180 |
| Lor4D | +0.3703 | 0.8576 | 67.8% | 180 |
| KR_like | +0.0734 | 0.7720 | 48.9% | 180 |
| TransPerc | +0.9502 | 1.0567 | 87.8% | 180 |

### Family × k table (mean A_entropy)

| family | k=1 | k=2 | k=4 |
|---|---:|---:|---:|
| Lor2D | +0.036 | +0.074 | +0.369 |
| Lor3D | +0.105 | +0.203 | +0.508 |
| Lor4D | +0.129 | +0.327 | +0.655 |
| KR_like | +0.097 | +0.043 | +0.080 |
| TransPerc | +0.406 | +0.839 | +1.606 |

## E2: Causal Depth Directional Growth

R_depth = (layers_forward - layers_base) / (layers_backward - layers_base)

Prediction: R_depth > 1

| family | mean R | std | R>1 % | n |
|---|---:|---:|---:|---:|
| Lor2D | 1.001 | 0.401 | 21.5% | 177 |
| Lor3D | 1.097 | 0.431 | 31.3% | 179 |
| Lor4D | 1.071 | 0.415 | 32.6% | 178 |
| KR_like | 1.037 | 0.349 | 27.2% | 180 |
| TransPerc | 1.130 | 0.501 | 33.5% | 176 |

## E3: Coarse-graining Robustness

Sign preservation rate: 57/75 (76.0%)

| family | preserved | total |
|---|---:|---:|
| Lor2D | 10 | 15 |
| Lor3D | 11 | 15 |
| Lor4D | 12 | 15 |
| KR_like | 11 | 15 |
| TransPerc | 13 | 15 |

## Interpretation

If A_entropy > 0 consistently across families and N values, this supports
the claim that time direction emerges from structural growth asymmetry:
augmenting toward the future (adding successors of maxima) increases
combinatorial entropy more than augmenting toward the past.

R_depth > 1 would confirm that forward growth deepens causal structure
faster than backward growth — a complementary structural arrow of time.
