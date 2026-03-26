# LSD-Well Large-N Scalability Test

N = [16, 20, 28, 36, 48, 64, 96, 128, 192, 256]
reps = 15 (N≤64), 10 (N>64)
Total samples: 2210
Generation time: 25.8s


## 1. Lor4D Centroid Evolution (Extended Range)

| N | c*(N) | σ(c) | w*(N) | σ(w) | d_eff | σ(d) |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| 16 | 0.1009 | 0.0827 | 0.5458 | 0.0868 | 4.005 | 0.182 |
| 20 | 0.1042 | 0.0745 | 0.5000 | 0.0516 | 3.871 | 0.208 |
| 28 | 0.1438 | 0.0725 | 0.4810 | 0.0701 | 3.930 | 0.273 |
| 36 | 0.1638 | 0.0861 | 0.4481 | 0.0737 | 3.993 | 0.223 |
| 48 | 0.2109 | 0.0565 | 0.3972 | 0.0399 | 3.925 | 0.215 |
| 64 | 0.2570 | 0.0519 | 0.3823 | 0.0445 | 3.860 | 0.161 |
| 96 | 0.2696 | 0.0233 | 0.3250 | 0.0322 | 3.986 | 0.098 |
| 128 | 0.2854 | 0.0336 | 0.3172 | 0.0145 | 3.956 | 0.103 |
| 192 | 0.3064 | 0.0217 | 0.2786 | 0.0240 | 3.987 | 0.091 |
| 256 | 0.3318 | 0.0231 | 0.2543 | 0.0230 | 3.982 | 0.069 |

**Full fit (all N)**: c*(N) = 0.3147 + -4.06/N
**Full fit (all N)**: w*(N) = 0.2791 + 4.75/N
**Small fit (N≤64)**: c*(N) = 0.2759 + -3.18/N
**Small fit (N≤64)**: w*(N) = 0.3390 + 3.39/N

### Extrapolation Accuracy (N≤64 fit → large N)

| N | c pred | c actual | |Δc| | w pred | w actual | |Δw| |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| 96 | 0.2428 | 0.2696 | 0.0268 | 0.3743 | 0.3250 | 0.0493 |
| 128 | 0.2511 | 0.2854 | 0.0344 | 0.3655 | 0.3172 | 0.0483 |
| 192 | 0.2593 | 0.3064 | 0.0470 | 0.3566 | 0.2786 | 0.0780 |
| 256 | 0.2635 | 0.3318 | 0.0684 | 0.3522 | 0.2543 | 0.0979 |


## 2. Fixed Weights (α=0.5, β=1.0, γ=5.0)

### MODE A: Oracle (same-N centroids)

**Beat rate**: 100.0%

| N | Lor4D Rank | Margin to runner-up |
|---|:----------:|:-------------------:|
| 16 | #1/17 | 0.1274 |
| 20 | #1/17 | 0.2360 |
| 28 | #1/17 | 0.2847 |
| 36 | #1/17 | 0.3591 |
| 48 | #1/17 | 0.6389 |
| 64 | #1/17 | 0.7158 |
| 96 | #1/17 | 0.9746 |
| 128 | #1/17 | 1.0164 |
| 192 | #1/17 | 1.1256 |
| 256 | #1/17 | 1.1941 |

### MODE B: Extrapolated (N≤64 fit applied to all N)

**Beat rate**: 100.0%

| N | Lor4D Rank | Margin to runner-up |
|---|:----------:|:-------------------:|
| 16 | #1/17 | 0.1145 |
| 20 | #1/17 | 0.2195 |
| 28 | #1/17 | 0.2500 |
| 36 | #1/17 | 0.3153 |
| 48 | #1/17 | 0.5947 |
| 64 | #1/17 | 0.6644 |
| 96 | #1/17 | 0.7507 |
| 128 | #1/17 | 0.7879 |
| 192 | #1/17 | 0.8085 |
| 256 | #1/17 | 0.8099 |


## 3. Grid Search: Optimal Weights at Large N

Do optimal weights shift when large-N data is included?

**Best weights (all N)**: α=0.5, β=0.5, γ=1.0
**Beat rate**: 100.0%

| N | Lor4D Rank |
|---|:----------:|
| 16 | #1/17 |
| 20 | #1/17 |
| 28 | #1/17 |
| 36 | #1/17 |
| 48 | #1/17 |
| 64 | #1/17 |
| 96 | #1/17 |
| 128 | #1/17 |
| 192 | #1/17 |
| 256 | #1/17 |


## 4. Full Rankings at N=128

### MODE A (Oracle, α=0.5 β=1.0 γ=5.0)

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.008459 |
| 2 | Lor5D ◆ | Lorentzian | 0.121734 |
| 3 | Lor3D ◆ | Lorentzian | 0.407408 |
| 4 | KR_2layer ● | KR-family | 1.024860 |
| 5 | KR_like ● | KR-family | 1.043036 |
| 6 | KR_4layer ● | KR-family | 1.595145 |
| 7 | MLR | Layered | 1.934933 |
| 8 | RLk4 | Layered | 2.153205 |
| 9 | RLk6_mid | Layered | 2.427423 |
| 10 | Lor2D ◆ | Lorentzian | 2.545134 |
| 11 | RLk6_tap | Layered | 2.743182 |
| 12 | IntOrder | Other | 3.319107 |
| 13 | RLk6_lj | Layered | 3.801393 |
| 14 | TransPerc | Other | 3.814501 |
| 15 | RLk6 | Layered | 4.058620 |
| 16 | AbsLayer | Layered | 4.107243 |
| 17 | RLk8 | Layered | 5.254997 |

### MODE B (Extrapolated, α=0.5 β=1.0 γ=5.0)

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.021291 |
| 2 | Lor5D ◆ | Lorentzian | 0.087015 |
| 3 | Lor3D ◆ | Lorentzian | 0.486562 |
| 4 | KR_2layer ● | KR-family | 0.809143 |
| 5 | KR_like ● | KR-family | 0.948004 |
| 6 | KR_4layer ● | KR-family | 1.560460 |
| 7 | MLR | Layered | 1.914706 |
| 8 | RLk4 | Layered | 2.160773 |
| 9 | RLk6_mid | Layered | 2.444455 |
| 10 | Lor2D ◆ | Lorentzian | 2.698610 |
| 11 | RLk6_tap | Layered | 2.764742 |
| 12 | IntOrder | Other | 3.441295 |
| 13 | RLk6_lj | Layered | 3.881865 |
| 14 | TransPerc | Other | 4.000177 |
| 15 | AbsLayer | Layered | 4.003791 |
| 16 | RLk6 | Layered | 4.137207 |
| 17 | RLk8 | Layered | 5.386513 |


## 4. Full Rankings at N=256

### MODE A (Oracle, α=0.5 β=1.0 γ=5.0)

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.005746 |
| 2 | Lor5D ◆ | Lorentzian | 0.140370 |
| 3 | Lor3D ◆ | Lorentzian | 0.366285 |
| 4 | KR_like ● | KR-family | 1.199893 |
| 5 | KR_2layer ● | KR-family | 1.346377 |
| 6 | KR_4layer ● | KR-family | 1.671575 |
| 7 | MLR | Layered | 2.354208 |
| 8 | RLk4 | Layered | 2.375010 |
| 9 | Lor2D ◆ | Lorentzian | 2.469961 |
| 10 | RLk6_mid | Layered | 2.613987 |
| 11 | IntOrder | Other | 3.161782 |
| 12 | RLk6_tap | Layered | 3.234474 |
| 13 | RLk6_lj | Layered | 4.137121 |
| 14 | AbsLayer | Layered | 4.205630 |
| 15 | RLk6 | Layered | 4.236519 |
| 16 | RLk8 | Layered | 5.454561 |
| 17 | TransPerc | Other | 7.100325 |

### MODE B (Extrapolated, α=0.5 β=1.0 γ=5.0)

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.058353 |
| 2 | Lor5D ◆ | Lorentzian | 0.087472 |
| 3 | Lor3D ◆ | Lorentzian | 0.546709 |
| 4 | KR_2layer ● | KR-family | 0.868262 |
| 5 | KR_like ● | KR-family | 0.966561 |
| 6 | KR_4layer ● | KR-family | 1.560636 |
| 7 | MLR | Layered | 2.281082 |
| 8 | RLk4 | Layered | 2.353554 |
| 9 | RLk6_mid | Layered | 2.582844 |
| 10 | Lor2D ◆ | Lorentzian | 2.765109 |
| 11 | RLk6_tap | Layered | 3.220831 |
| 12 | IntOrder | Other | 3.393907 |
| 13 | AbsLayer | Layered | 3.980331 |
| 14 | RLk6_lj | Layered | 4.221664 |
| 15 | RLk6 | Layered | 4.311876 |
| 16 | RLk8 | Layered | 5.621818 |
| 17 | TransPerc | Other | 7.499131 |


## 5. Margin Scaling: How Does Discrimination Improve with N?

| N | Margin A (Oracle) | Margin B (Extrap) | Margin ratio B/A |
|---|:-:|:-:|:-:|
| 16 | 0.127412 | 0.114519 | 0.899 |
| 20 | 0.235964 | 0.219506 | 0.930 |
| 28 | 0.284741 | 0.250020 | 0.878 |
| 36 | 0.359094 | 0.315278 | 0.878 |
| 48 | 0.638924 | 0.594655 | 0.931 |
| 64 | 0.715809 | 0.664438 | 0.928 |
| 96 | 0.974611 | 0.750670 | 0.770 |
| 128 | 1.016401 | 0.787852 | 0.775 |
| 192 | 1.125631 | 0.808536 | 0.718 |
| 256 | 1.194147 | 0.809909 | 0.678 |


## 6. Feature Space at N=256

| Family | d_eff | c₁/c₀ | width | F (Oracle) |
|--------|:-:|:-:|:-:|:-:|
| AbsLayer | 1.284 | 0.0000 | 0.4918 | 4.205630 |
| IntOrder | 1.556 | 0.5119 | 0.0961 | 3.161782 |
| KR_2layer | 3.877 | 0.0000 | 0.7500 | 1.346377 |
| KR_4layer | 2.275 | 0.0000 | 0.3750 | 1.671575 |
| KR_like | 2.745 | 0.0000 | 0.5000 | 1.199893 |
| Lor2D | 1.949 | 0.7631 | 0.0668 | 2.469961 |
| Lor3D | 3.267 | 0.5310 | 0.1516 | 0.366285 |
| Lor4D | 3.982 | 0.3318 | 0.2543 | 0.005746 |
| Lor5D | 4.397 | 0.1868 | 0.3418 | 0.140370 |
| MLR | 2.160 | 0.6375 | 0.4254 | 2.354208 |
| RLk4 | 1.874 | 0.0027 | 0.2840 | 2.375010 |
| RLk6 | 1.119 | 0.0758 | 0.1953 | 4.236519 |
| RLk6_lj | 1.145 | 0.1318 | 0.1938 | 4.137121 |
| RLk6_mid | 1.764 | 0.0213 | 0.2965 | 2.613987 |
| RLk6_tap | 1.489 | 0.0654 | 0.2848 | 3.234474 |
| RLk8 | 0.719 | 0.4682 | 0.1562 | 5.454561 |
| TransPerc | 0.629 | 1.4346 | 0.0547 | 7.100325 |


## 7. Summary

- MODE A (Oracle): Lor4D #1 at ALL N? **YES**
- MODE B (Extrap): Lor4D #1 at ALL N? **YES**
- Grid search best: α=0.5, β=0.5, γ=1.0 (was α=0.5, β=1.0, γ=5.0)
- Finite-size scaling prediction error at N=256: Δc=0.0684, Δw=0.0979