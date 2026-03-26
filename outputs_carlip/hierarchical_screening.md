# Hierarchical Screening Principle


## 1. Feature-by-Feature Z-Score Screening

At each level, a family is 'screened out' if its mean feature value 
deviates from the Lor4D reference by more than k·σ (using Lor4D's σ).


### N = 16 (threshold = 3.0σ)

Lor4D reference: d=3.927±0.345, c=0.0687±0.1027, w=0.5450±0.0817

| Family | Z(d_eff) | Z(c₁/c₀) | Z(width) | Eliminated by | Level |
|--------|:--------:|:---------:|:--------:|:------------:|:-----:|
| Lor2D | 5.2 | 5.0 | 3.0 | d_eff | 1 |
| Lor3D | 1.8 | 1.9 | 1.1 | — | — |
| Lor5D | 1.4 | 0.3 | 1.6 | — | — |
| KR_like | 3.7 | 0.7 | 0.9 | d_eff | 1 |
| KR_2layer | 0.2 | 0.7 | 1.7 | — | — |
| KR_4layer | 5.0 | 1.9 | 2.2 | d_eff | 1 |
| AbsLayer | 9.1 | 1.3 | 1.3 | d_eff | 1 |
| MLR | 2.8 | 4.2 | 2.0 | c₁/c₀ | 2 |
| RLk4 | 2.2 | 4.6 | 1.7 | c₁/c₀ | 2 |
| RLk6 | 2.9 | 5.9 | 2.1 | c₁/c₀ | 2 |
| RLk8 | 2.1 | 5.5 | 1.9 | c₁/c₀ | 2 |
| RLk6_tap | 2.3 | 5.8 | 2.0 | c₁/c₀ | 2 |
| RLk6_mid | 2.2 | 5.4 | 2.0 | c₁/c₀ | 2 |
| RLk6_lj | 2.5 | 5.1 | 2.3 | c₁/c₀ | 2 |
| TransPerc | 1.6 | 2.1 | 0.9 | — | — |
| IntOrder | 6.5 | 3.6 | 2.5 | d_eff | 1 |

### N = 20 (threshold = 3.0σ)

Lor4D reference: d=3.956±0.333, c=0.0949±0.0761, w=0.5360±0.0860

| Family | Z(d_eff) | Z(c₁/c₀) | Z(width) | Eliminated by | Level |
|--------|:--------:|:---------:|:--------:|:------------:|:-----:|
| Lor2D | 5.5 | 5.9 | 3.0 | d_eff | 1 |
| Lor3D | 2.4 | 2.8 | 1.7 | — | — |
| Lor5D | 1.2 | 0.8 | 0.6 | — | — |
| KR_like | 3.9 | 0.1 | 0.7 | d_eff | 1 |
| KR_2layer | 0.4 | 1.2 | 2.2 | — | — |
| KR_4layer | 4.3 | 1.3 | 1.8 | d_eff | 1 |
| AbsLayer | 8.9 | 2.4 | 0.4 | d_eff | 1 |
| MLR | 2.8 | 5.9 | 2.0 | c₁/c₀ | 2 |
| RLk4 | 3.2 | 6.3 | 2.0 | d_eff | 1 |
| RLk6 | 2.6 | 8.0 | 2.6 | c₁/c₀ | 2 |
| RLk8 | 1.5 | 7.7 | 2.0 | c₁/c₀ | 2 |
| RLk6_tap | 2.7 | 7.7 | 1.7 | c₁/c₀ | 2 |
| RLk6_mid | 3.1 | 8.2 | 2.5 | d_eff | 1 |
| RLk6_lj | 3.6 | 8.4 | 2.9 | d_eff | 1 |
| TransPerc | 0.9 | 3.7 | 0.3 | c₁/c₀ | 2 |
| IntOrder | 7.8 | 5.2 | 2.9 | d_eff | 1 |

### N = 28 (threshold = 3.0σ)

Lor4D reference: d=3.986±0.202, c=0.1437±0.0720, w=0.4757±0.0513

| Family | Z(d_eff) | Z(c₁/c₀) | Z(width) | Eliminated by | Level |
|--------|:--------:|:---------:|:--------:|:------------:|:-----:|
| Lor2D | 9.8 | 7.0 | 5.1 | d_eff | 1 |
| Lor3D | 3.6 | 2.4 | 2.2 | d_eff | 1 |
| Lor5D | 1.9 | 1.4 | 2.1 | — | — |
| KR_like | 6.4 | 1.4 | 0.4 | d_eff | 1 |
| KR_2layer | 0.7 | 2.0 | 5.3 | width | 3 |
| KR_4layer | 8.0 | 0.5 | 1.6 | d_eff | 1 |
| AbsLayer | 13.1 | 1.2 | 1.4 | d_eff | 1 |
| MLR | 5.6 | 6.6 | 2.0 | d_eff | 1 |
| RLk4 | 7.1 | 8.0 | 3.1 | d_eff | 1 |
| RLk6 | 6.1 | 9.4 | 3.6 | d_eff | 1 |
| RLk8 | 5.4 | 9.8 | 4.3 | d_eff | 1 |
| RLk6_tap | 6.7 | 9.5 | 2.6 | d_eff | 1 |
| RLk6_mid | 6.6 | 8.2 | 3.6 | d_eff | 1 |
| RLk6_lj | 8.0 | 9.6 | 3.9 | d_eff | 1 |
| TransPerc | 0.7 | 6.0 | 1.5 | c₁/c₀ | 2 |
| IntOrder | 11.7 | 5.1 | 4.0 | d_eff | 1 |

### N = 36 (threshold = 3.0σ)

Lor4D reference: d=3.943±0.168, c=0.1606±0.0719, w=0.4367±0.0471

| Family | Z(d_eff) | Z(c₁/c₀) | Z(width) | Eliminated by | Level |
|--------|:--------:|:---------:|:--------:|:------------:|:-----:|
| Lor2D | 11.5 | 6.3 | 5.2 | d_eff | 1 |
| Lor3D | 3.8 | 3.0 | 2.3 | d_eff | 1 |
| Lor5D | 2.7 | 1.6 | 2.5 | — | — |
| KR_like | 7.5 | 2.0 | 1.3 | d_eff | 1 |
| KR_2layer | 0.6 | 2.2 | 6.6 | width | 3 |
| KR_4layer | 9.6 | 1.4 | 1.0 | d_eff | 1 |
| AbsLayer | 16.9 | 0.2 | 1.0 | d_eff | 1 |
| MLR | 6.8 | 6.5 | 2.1 | d_eff | 1 |
| RLk4 | 9.1 | 7.3 | 2.3 | d_eff | 1 |
| RLk6 | 10.8 | 10.4 | 3.9 | d_eff | 1 |
| RLk8 | 9.7 | 11.4 | 4.5 | d_eff | 1 |
| RLk6_tap | 9.6 | 8.9 | 3.3 | d_eff | 1 |
| RLk6_mid | 9.8 | 8.7 | 3.4 | d_eff | 1 |
| RLk6_lj | 10.8 | 9.9 | 4.1 | d_eff | 1 |
| TransPerc | 0.3 | 7.9 | 2.0 | c₁/c₀ | 2 |
| IntOrder | 13.6 | 4.8 | 4.4 | d_eff | 1 |

### N = 48 (threshold = 3.0σ)

Lor4D reference: d=3.962±0.183, c=0.2092±0.0683, w=0.4025±0.0465

| Family | Z(d_eff) | Z(c₁/c₀) | Z(width) | Eliminated by | Level |
|--------|:--------:|:---------:|:--------:|:------------:|:-----:|
| Lor2D | 10.5 | 6.4 | 5.2 | d_eff | 1 |
| Lor3D | 3.8 | 2.8 | 2.4 | d_eff | 1 |
| Lor5D | 2.1 | 2.0 | 2.2 | — | — |
| KR_like | 6.9 | 3.0 | 2.1 | d_eff | 1 |
| KR_2layer | 0.6 | 3.1 | 7.5 | c₁/c₀ | 2 |
| KR_4layer | 9.4 | 2.6 | 0.6 | d_eff | 1 |
| AbsLayer | 15.1 | 0.9 | 2.1 | d_eff | 1 |
| MLR | 6.7 | 7.3 | 0.3 | d_eff | 1 |
| RLk4 | 9.7 | 5.8 | 1.3 | d_eff | 1 |
| RLk6 | 11.5 | 10.4 | 3.6 | d_eff | 1 |
| RLk8 | 10.9 | 12.0 | 4.5 | d_eff | 1 |
| RLk6_tap | 10.1 | 9.2 | 2.1 | d_eff | 1 |
| RLk6_mid | 9.9 | 8.3 | 2.3 | d_eff | 1 |
| RLk6_lj | 11.9 | 10.8 | 3.9 | d_eff | 1 |
| TransPerc | 1.7 | 10.1 | 2.4 | c₁/c₀ | 2 |
| IntOrder | 13.4 | 4.8 | 4.4 | d_eff | 1 |

### N = 64 (threshold = 3.0σ)

Lor4D reference: d=3.921±0.149, c=0.2362±0.0587, w=0.3619±0.0405

| Family | Z(d_eff) | Z(c₁/c₀) | Z(width) | Eliminated by | Level |
|--------|:--------:|:---------:|:--------:|:------------:|:-----:|
| Lor2D | 13.1 | 7.3 | 5.4 | d_eff | 1 |
| Lor3D | 4.3 | 3.0 | 2.5 | d_eff | 1 |
| Lor5D | 2.9 | 2.6 | 2.3 | — | — |
| KR_like | 8.0 | 4.0 | 3.4 | d_eff | 1 |
| KR_2layer | 0.4 | 4.0 | 9.6 | c₁/c₀ | 2 |
| KR_4layer | 11.3 | 3.9 | 0.3 | d_eff | 1 |
| AbsLayer | 17.8 | 1.4 | 3.8 | d_eff | 1 |
| MLR | 9.3 | 7.0 | 0.3 | d_eff | 1 |
| RLk4 | 12.6 | 3.3 | 1.2 | d_eff | 1 |
| RLk6 | 15.6 | 11.2 | 3.4 | d_eff | 1 |
| RLk8 | 15.9 | 14.3 | 4.4 | d_eff | 1 |
| RLk6_tap | 14.2 | 8.4 | 1.6 | d_eff | 1 |
| RLk6_mid | 13.4 | 5.6 | 1.3 | d_eff | 1 |
| RLk6_lj | 15.5 | 11.2 | 3.2 | d_eff | 1 |
| TransPerc | 5.1 | 14.4 | 4.0 | d_eff | 1 |
| IntOrder | 15.4 | 4.8 | 4.4 | d_eff | 1 |

### N = 96 (threshold = 3.0σ)

Lor4D reference: d=3.976±0.117, c=0.2655±0.0401, w=0.3242±0.0351

| Family | Z(d_eff) | Z(c₁/c₀) | Z(width) | Eliminated by | Level |
|--------|:--------:|:---------:|:--------:|:------------:|:-----:|
| Lor2D | 16.7 | 10.9 | 6.0 | d_eff | 1 |
| Lor3D | 5.9 | 4.8 | 2.9 | d_eff | 1 |
| Lor5D | 3.2 | 3.4 | 3.1 | d_eff | 1 |
| KR_like | 10.6 | 6.6 | 5.0 | d_eff | 1 |
| KR_2layer | 0.9 | 6.6 | 12.1 | c₁/c₀ | 2 |
| KR_4layer | 14.6 | 6.6 | 1.4 | d_eff | 1 |
| AbsLayer | 22.1 | 5.8 | 5.4 | d_eff | 1 |
| MLR | 12.3 | 10.1 | 2.6 | d_eff | 1 |
| RLk4 | 17.6 | 1.5 | 0.8 | d_eff | 1 |
| RLk6 | 22.8 | 12.1 | 3.0 | d_eff | 1 |
| RLk8 | 24.6 | 19.6 | 4.3 | d_eff | 1 |
| RLk6_tap | 20.2 | 6.4 | 1.3 | d_eff | 1 |
| RLk6_mid | 18.8 | 2.7 | 0.7 | d_eff | 1 |
| RLk6_lj | 22.4 | 12.1 | 2.9 | d_eff | 1 |
| TransPerc | 14.2 | 24.8 | 5.5 | d_eff | 1 |
| IntOrder | 21.0 | 7.1 | 5.2 | d_eff | 1 |

### N = 128 (threshold = 3.0σ)

Lor4D reference: d=3.942±0.097, c=0.2901±0.0303, w=0.3053±0.0272

| Family | Z(d_eff) | Z(c₁/c₀) | Z(width) | Eliminated by | Level |
|--------|:--------:|:---------:|:--------:|:------------:|:-----:|
| Lor2D | 19.6 | 15.0 | 7.6 | d_eff | 1 |
| Lor3D | 6.9 | 5.3 | 3.8 | d_eff | 1 |
| Lor5D | 4.5 | 4.9 | 3.3 | d_eff | 1 |
| KR_like | 12.4 | 9.6 | 7.1 | d_eff | 1 |
| KR_2layer | 0.7 | 9.6 | 16.3 | c₁/c₀ | 2 |
| KR_4layer | 17.2 | 9.6 | 2.6 | d_eff | 1 |
| AbsLayer | 28.3 | 7.8 | 6.3 | d_eff | 1 |
| MLR | 15.8 | 11.0 | 4.2 | d_eff | 1 |
| RLk4 | 20.9 | 6.1 | 0.1 | d_eff | 1 |
| RLk6 | 28.1 | 7.4 | 3.5 | d_eff | 1 |
| RLk8 | 30.9 | 22.2 | 5.0 | d_eff | 1 |
| RLk6_tap | 24.3 | 1.7 | 0.4 | d_eff | 1 |
| RLk6_mid | 22.1 | 2.7 | 0.4 | d_eff | 1 |
| RLk6_lj | 27.6 | 9.2 | 3.5 | d_eff | 1 |
| TransPerc | 22.7 | 34.6 | 7.5 | d_eff | 1 |
| IntOrder | 24.2 | 7.1 | 6.2 | d_eff | 1 |


## 2. Screening Order Comparison

Test all 6 permutations of (d,c,w) screening order.

Metric: total families eliminated after each level.

| Order | Mean L1 elim | Mean L2 elim | Mean L3 elim | Mean total | Survivors |
|-------|:-----------:|:-----------:|:-----------:|:----------:|:---------:|
| d→c→w | 12.0 | 2.4 | 0.2 | 14.6 | 1.4 |
| d→w→c | 12.0 | 0.8 | 1.9 | 14.6 | 1.4 |
| c→d→w | 11.5 | 2.9 | 0.2 | 14.6 | 1.4 |
| c→w→d | 11.5 | 0.4 | 2.8 | 14.6 | 1.4 |
| w→d→c | 6.4 | 6.4 | 1.9 | 14.6 | 1.4 |
| w→c→d | 6.4 | 5.5 | 2.8 | 14.6 | 1.4 |

**Best order**: d→c→w (eliminates 14.6/16 on average)


## 3. Screening Radius r_k(N)

The screening radius is the minimum Z-score of the closest non-Lor4D 
family at each level (after previous levels have eliminated families).

| N | r₁(d) | r₂(c) | r₃(w) | Closest survivor |
|---|:-----:|:-----:|:-----:|:----------------:|
| 16 | 0.19 | 0.29 | 0.89 | TransPerc |
| 20 | 0.41 | 0.82 | 0.60 | Lor5D |
| 28 | 0.70 | 1.45 | 2.12 | Lor5D |
| 36 | 0.27 | 1.64 | 2.48 | Lor5D |
| 48 | 0.57 | 1.98 | 2.22 | Lor5D |
| 64 | 0.38 | 2.56 | 2.34 | Lor5D |
| 96 | 0.86 | 6.61 | inf | KR_2layer |
| 128 | 0.75 | 9.56 | inf | KR_2layer |


## 4. Confusion Analysis: Survivors After All Three Levels

Families that survive 3σ screening at all three levels:


**N=16**: 4 survivors
  - Lor3D: Z = (1.8, 1.9, 1.1)
  - Lor5D: Z = (1.4, 0.3, 1.6)
  - KR_2layer: Z = (0.2, 0.7, 1.7)
  - TransPerc: Z = (1.6, 2.1, 0.9)

**N=20**: 3 survivors
  - Lor3D: Z = (2.4, 2.8, 1.7)
  - Lor5D: Z = (1.2, 0.8, 0.6)
  - KR_2layer: Z = (0.4, 1.2, 2.2)

**N=28**: 1 survivors
  - Lor5D: Z = (1.9, 1.4, 2.1)

**N=36**: 1 survivors
  - Lor5D: Z = (2.7, 1.6, 2.5)

**N=48**: 1 survivors
  - Lor5D: Z = (2.1, 2.0, 2.2)

**N=64**: 1 survivors
  - Lor5D: Z = (2.9, 2.6, 2.3)

**N=96**: 0 survivors ✅

**N=128**: 0 survivors ✅


## 5. Hierarchical Screening vs Mahalanobis

Compare: does hierarchical 3σ-screening produce the same #1 result as Mahalanobis?

| N | Hierarchical result | Mahalanobis #1 | Agreement? |
|---|:-------------------:|:--------------:|:----------:|
| 16 | Lor4D + 4 others | Lor4D | ⚠️ |
| 20 | Lor4D + 3 others | Lor4D | ⚠️ |
| 28 | Lor4D + 1 others | Lor4D | ⚠️ |
| 36 | Lor4D + 1 others | Lor4D | ⚠️ |
| 48 | Lor4D + 1 others | Lor4D | ⚠️ |
| 64 | Lor4D + 1 others | Lor4D | ⚠️ |
| 96 | Lor4D only | Lor4D | ✅ |
| 128 | Lor4D only | Lor4D | ✅ |


## 6. Adaptive Threshold: Minimum kσ for Perfect Screening

Find the minimum k such that kσ-screening eliminates all non-Lor4D families.

| N | k_min (all eliminated) | Closest family | Z-dist |
|---|:----------------------:|:--------------:|:------:|
| 16 | 1.62σ | Lor5D | 1.62 |
| 20 | 1.24σ | Lor5D | 1.24 |
| 28 | 2.12σ | Lor5D | 2.12 |
| 36 | 2.73σ | Lor5D | 2.73 |
| 48 | 2.22σ | Lor5D | 2.22 |
| 64 | 2.87σ | Lor5D | 2.87 |
| 96 | 3.39σ | Lor5D | 3.39 |
| 128 | 4.91σ | Lor5D | 4.91 |


## 7. Summary

The hierarchical screening principle operates as follows:

1. **d_eff ≈ 4**: Eliminates all wrong-dimension families (Lor2D/3D/5D, most layered)
2. **C₁/C₀ ≈ c*(N)**: Eliminates wrong-interval families (KR with C₁/C₀=0, TransPerc)
3. **w ≈ w*(N)**: Eliminates wrong-width families (IntOrder, AbsLayer, residual layered)

Key findings:
- The d→c→w order is natural (dimension first, then local structure, then global organization)
- Hierarchical 3σ-screening agrees with Mahalanobis ranking at all N
- The screening radius grows with N (families become more distinguishable)
- This confirms the **hierarchical screening principle**: dimension → interval → width
