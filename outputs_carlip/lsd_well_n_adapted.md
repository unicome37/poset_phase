# LSD-Well N-Adapted: N-dependent Well Centers

**Upgrade**: well centers c*(N), w*(N) no longer fixed constants.

N = [16, 20, 28, 36, 48, 64], reps = 15


## 1. Lor4D Centroid Evolution

| N | c*(N) = C₁/C₀ | σ(c) | w*(N) = width | σ(w) | d_eff | σ(d) |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| 16 | 0.0940 | 0.0851 | 0.5750 | 0.1075 | 4.005 | 0.354 |
| 20 | 0.1647 | 0.1010 | 0.4833 | 0.0394 | 3.858 | 0.240 |
| 28 | 0.1441 | 0.0938 | 0.4857 | 0.0484 | 3.914 | 0.268 |
| 36 | 0.1556 | 0.0613 | 0.4315 | 0.0485 | 3.931 | 0.209 |
| 48 | 0.2078 | 0.0499 | 0.4167 | 0.0437 | 3.948 | 0.176 |
| 64 | 0.2290 | 0.0437 | 0.3688 | 0.0411 | 3.911 | 0.139 |

**Finite-size fit**: c*(N) = 0.2485 + -2.33/N
**Finite-size fit**: w*(N) = 0.3255 + 3.80/N
**d* = 4 (fixed, physics input)**


## 2. MODE A: Oracle N-Adapted (same-N centroids)

F(N) = α·(d_eff − 4)² + β·(C₁/C₀ − c*(N))² + γ·(w − w*(N))²

where c*(N), w*(N) = Lor4D centroid at each N.

**Best**: α=0.5, β=1.0, γ=5.0
**Lor4D beats 100.0% of non-Lor across all N**

| N | Lor4D Rank | c*(N) | w*(N) |
|---|:----------:|:-----:|:-----:|
| 16 | #1/17 | 0.0940 | 0.5750 |
| 20 | #1/17 | 0.1647 | 0.4833 |
| 28 | #1/17 | 0.1441 | 0.4857 |
| 36 | #1/17 | 0.1556 | 0.4315 |
| 48 | #1/17 | 0.2078 | 0.4167 |
| 64 | #1/17 | 0.2290 | 0.3688 |

### MODE A: N=16

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.1277 |
| 2 | KR_2layer ● | KR-family | 0.1278 |
| 3 | Lor5D ◆ | Lorentzian | 0.1788 |
| 4 | TransPerc | Other | 0.2250 |
| 5 | RLk6 | Layered | 0.5494 |
| 6 | Lor3D ◆ | Lorentzian | 0.5823 |
| 7 | RLk6_tap | Layered | 0.8715 |
| 8 | KR_like ● | KR-family | 0.9736 |
| 9 | RLk6_mid | Layered | 0.9896 |
| 10 | MLR | Layered | 1.0068 |
| 11 | RLk4 | Layered | 1.1005 |
| 12 | RLk8 | Layered | 1.1338 |
| 13 | RLk6_lj | Layered | 1.3178 |
| 14 | KR_4layer ● | KR-family | 1.8242 |
| 15 | Lor2D ◆ | Lorentzian | 2.4043 |
| 16 | IntOrder | Other | 4.3708 |
| 17 | AbsLayer | Layered | 5.0751 |

### MODE A: N=20

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0569 |
| 2 | TransPerc | Other | 0.2072 |
| 3 | Lor5D ◆ | Lorentzian | 0.2217 |
| 4 | Lor3D ◆ | Lorentzian | 0.3729 |
| 5 | KR_2layer ● | KR-family | 0.3852 |
| 6 | RLk8 | Layered | 0.6963 |
| 7 | MLR | Layered | 0.7725 |
| 8 | KR_like ● | KR-family | 0.8796 |
| 9 | RLk6 | Layered | 0.9497 |
| 10 | RLk6_lj | Layered | 1.1002 |
| 11 | RLk4 | Layered | 1.1368 |
| 12 | RLk6_mid | Layered | 1.1646 |
| 13 | RLk6_tap | Layered | 1.2272 |
| 14 | KR_4layer ● | KR-family | 1.3351 |
| 15 | Lor2D ◆ | Lorentzian | 2.2962 |
| 16 | IntOrder | Other | 3.3102 |
| 17 | AbsLayer | Layered | 4.9718 |

### MODE A: N=28

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0600 |
| 2 | Lor5D ◆ | Lorentzian | 0.1371 |
| 3 | TransPerc | Other | 0.3741 |
| 4 | KR_2layer ● | KR-family | 0.3777 |
| 5 | Lor3D ◆ | Lorentzian | 0.4640 |
| 6 | KR_like ● | KR-family | 0.9167 |
| 7 | MLR | Layered | 0.9777 |
| 8 | KR_4layer ● | KR-family | 1.4100 |
| 9 | RLk6_tap | Layered | 1.4383 |
| 10 | RLk4 | Layered | 1.5500 |
| 11 | RLk8 | Layered | 1.7494 |
| 12 | RLk6_mid | Layered | 1.8722 |
| 13 | RLk6 | Layered | 1.9271 |
| 14 | RLk6_lj | Layered | 2.2537 |
| 15 | Lor2D ◆ | Lorentzian | 3.0200 |
| 16 | IntOrder | Other | 3.5039 |
| 17 | AbsLayer | Layered | 5.1099 |

### MODE A: N=36

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0398 |
| 2 | Lor5D ◆ | Lorentzian | 0.1675 |
| 3 | Lor3D ◆ | Lorentzian | 0.3848 |
| 4 | KR_2layer ● | KR-family | 0.5489 |
| 5 | TransPerc | Other | 0.5802 |
| 6 | KR_like ● | KR-family | 0.8872 |
| 7 | MLR | Layered | 1.3668 |
| 8 | KR_4layer ● | KR-family | 1.3688 |
| 9 | RLk4 | Layered | 1.6690 |
| 10 | RLk6_tap | Layered | 1.9306 |
| 11 | RLk6_mid | Layered | 1.9633 |
| 12 | RLk8 | Layered | 2.0123 |
| 13 | RLk6 | Layered | 2.2669 |
| 14 | Lor2D ◆ | Lorentzian | 2.3793 |
| 15 | RLk6_lj | Layered | 2.5927 |
| 16 | IntOrder | Other | 3.0258 |
| 17 | AbsLayer | Layered | 5.0034 |

### MODE A: N=48

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0289 |
| 2 | Lor5D ◆ | Lorentzian | 0.1349 |
| 3 | Lor3D ◆ | Lorentzian | 0.4362 |
| 4 | KR_2layer ● | KR-family | 0.6106 |
| 5 | KR_like ● | KR-family | 0.8976 |
| 6 | TransPerc | Other | 0.9707 |
| 7 | MLR | Layered | 1.3306 |
| 8 | KR_4layer ● | KR-family | 1.5951 |
| 9 | RLk4 | Layered | 1.8472 |
| 10 | RLk6_mid | Layered | 2.1422 |
| 11 | Lor2D ◆ | Lorentzian | 2.4267 |
| 12 | RLk6_tap | Layered | 2.4594 |
| 13 | RLk6 | Layered | 2.9290 |
| 14 | RLk6_lj | Layered | 3.0143 |
| 15 | IntOrder | Other | 3.2618 |
| 16 | RLk8 | Layered | 3.3511 |
| 17 | AbsLayer | Layered | 4.0237 |

### MODE A: N=64

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0239 |
| 2 | Lor5D ◆ | Lorentzian | 0.1330 |
| 3 | Lor3D ◆ | Lorentzian | 0.3907 |
| 4 | KR_2layer ● | KR-family | 0.7901 |
| 5 | KR_like ● | KR-family | 0.9493 |
| 6 | MLR | Layered | 1.1626 |
| 7 | TransPerc | Other | 1.4283 |
| 8 | KR_4layer ● | KR-family | 1.5737 |
| 9 | RLk4 | Layered | 1.9468 |
| 10 | RLk6_mid | Layered | 2.2996 |
| 11 | Lor2D ◆ | Lorentzian | 2.3906 |
| 12 | RLk6_tap | Layered | 2.6245 |
| 13 | IntOrder | Other | 3.1903 |
| 14 | RLk6 | Layered | 3.4414 |
| 15 | RLk6_lj | Layered | 3.4792 |
| 16 | RLk8 | Layered | 4.0024 |
| 17 | AbsLayer | Layered | 4.5171 |


## 3. MODE B: Extrapolated N-Adapted

Fit c*(N) = c_∞ + b_c/N and w*(N) = w_∞ + b_w/N from N ≤ 36,

then predict well centers for N = 48, 64.


**Train fit (N ≤ 36)**: c*(N) = 0.2010 + -1.40/N
**Train fit (N ≤ 36)**: w*(N) = 0.3379 + 3.55/N

| N | c*(N) pred | c*(N) actual | Δc | w*(N) pred | w*(N) actual | Δw |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| 16 | 0.1138 | 0.0940 | 0.0198 | 0.5595 | 0.5750 | 0.0155 |
| 20 | 0.1312 | 0.1647 | 0.0335 | 0.5152 | 0.4833 | 0.0318 |
| 28 | 0.1512 | 0.1441 | 0.0070 | 0.4645 | 0.4857 | 0.0212 |
| 36 | 0.1623 | 0.1556 | 0.0066 | 0.4364 | 0.4315 | 0.0049 |
| 48 (test) | 0.1720 | 0.2078 | 0.0359 | 0.4118 | 0.4167 | 0.0049 |
| 64 (test) | 0.1792 | 0.2290 | 0.0498 | 0.3933 | 0.3688 | 0.0245 |

**Best**: α=0.5, β=1.0, γ=5.0
**Lor4D beats 100.0% of non-Lor across all N**

| N | Lor4D Rank |
|---|:----------:|
| 16 | #1/17 |
| 20 | #1/17 |
| 28 | #1/17 |
| 36 | #1/17 |
| 48 | #1/17 |
| 64 | #1/17 |

### MODE B: N=16

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.1293 |
| 2 | KR_2layer ● | KR-family | 0.1519 |
| 3 | Lor5D ◆ | Lorentzian | 0.1947 |
| 4 | TransPerc | Other | 0.2188 |
| 5 | RLk6 | Layered | 0.5093 |
| 6 | Lor3D ◆ | Lorentzian | 0.5590 |
| 7 | RLk6_tap | Layered | 0.8211 |
| 8 | RLk6_mid | Layered | 0.9375 |
| 9 | KR_like ● | KR-family | 0.9571 |
| 10 | MLR | Layered | 0.9650 |
| 11 | RLk4 | Layered | 1.0544 |
| 12 | RLk8 | Layered | 1.0814 |
| 13 | RLk6_lj | Layered | 1.2624 |
| 14 | KR_4layer ● | KR-family | 1.7840 |
| 15 | Lor2D ◆ | Lorentzian | 2.3463 |
| 16 | IntOrder | Other | 4.3184 |
| 17 | AbsLayer | Layered | 5.0583 |

### MODE B: N=20

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0631 |
| 2 | Lor5D ◆ | Lorentzian | 0.1821 |
| 3 | TransPerc | Other | 0.2188 |
| 4 | KR_2layer ● | KR-family | 0.2997 |
| 5 | Lor3D ◆ | Lorentzian | 0.4064 |
| 6 | RLk8 | Layered | 0.7727 |
| 7 | MLR | Layered | 0.8387 |
| 8 | KR_like ● | KR-family | 0.8796 |
| 9 | RLk6 | Layered | 1.0276 |
| 10 | RLk6_lj | Layered | 1.1876 |
| 11 | RLk4 | Layered | 1.2163 |
| 12 | RLk6_mid | Layered | 1.2625 |
| 13 | RLk6_tap | Layered | 1.3107 |
| 14 | KR_4layer ● | KR-family | 1.3694 |
| 15 | Lor2D ◆ | Lorentzian | 2.3957 |
| 16 | IntOrder | Other | 3.3945 |
| 17 | AbsLayer | Layered | 4.9788 |

### MODE B: N=28

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0623 |
| 2 | Lor5D ◆ | Lorentzian | 0.1551 |
| 3 | TransPerc | Other | 0.3552 |
| 4 | KR_2layer ● | KR-family | 0.4371 |
| 5 | Lor3D ◆ | Lorentzian | 0.4389 |
| 6 | KR_like ● | KR-family | 0.9230 |
| 7 | MLR | Layered | 0.9444 |
| 8 | KR_4layer ● | KR-family | 1.3931 |
| 9 | RLk6_tap | Layered | 1.4006 |
| 10 | RLk4 | Layered | 1.5124 |
| 11 | RLk8 | Layered | 1.6960 |
| 12 | RLk6_mid | Layered | 1.8175 |
| 13 | RLk6 | Layered | 1.8742 |
| 14 | RLk6_lj | Layered | 2.1966 |
| 15 | Lor2D ◆ | Lorentzian | 2.9571 |
| 16 | IntOrder | Other | 3.4548 |
| 17 | AbsLayer | Layered | 5.1078 |

### MODE B: N=36

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0400 |
| 2 | Lor5D ◆ | Lorentzian | 0.1632 |
| 3 | Lor3D ◆ | Lorentzian | 0.3881 |
| 4 | KR_2layer ● | KR-family | 0.5355 |
| 5 | TransPerc | Other | 0.5775 |
| 6 | KR_like ● | KR-family | 0.8858 |
| 7 | MLR | Layered | 1.3620 |
| 8 | KR_4layer ● | KR-family | 1.3722 |
| 9 | RLk4 | Layered | 1.6687 |
| 10 | RLk6_tap | Layered | 1.9274 |
| 11 | RLk6_mid | Layered | 1.9639 |
| 12 | RLk8 | Layered | 2.0123 |
| 13 | RLk6 | Layered | 2.2666 |
| 14 | Lor2D ◆ | Lorentzian | 2.3848 |
| 15 | RLk6_lj | Layered | 2.5925 |
| 16 | IntOrder | Other | 3.0314 |
| 17 | AbsLayer | Layered | 4.9984 |

### MODE B: N=48

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0303 |
| 2 | Lor5D ◆ | Lorentzian | 0.1327 |
| 3 | Lor3D ◆ | Lorentzian | 0.4464 |
| 4 | KR_2layer ● | KR-family | 0.6135 |
| 5 | KR_like ● | KR-family | 0.8886 |
| 6 | TransPerc | Other | 1.0191 |
| 7 | MLR | Layered | 1.3660 |
| 8 | KR_4layer ● | KR-family | 1.5811 |
| 9 | RLk4 | Layered | 1.8688 |
| 10 | RLk6_mid | Layered | 2.1729 |
| 11 | Lor2D ◆ | Lorentzian | 2.4499 |
| 12 | RLk6_tap | Layered | 2.4957 |
| 13 | RLk6 | Layered | 2.9721 |
| 14 | RLk6_lj | Layered | 3.0580 |
| 15 | IntOrder | Other | 3.2753 |
| 16 | RLk8 | Layered | 3.4007 |
| 17 | AbsLayer | Layered | 4.0157 |

### MODE B: N=64

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0294 |
| 2 | Lor5D ◆ | Lorentzian | 0.1002 |
| 3 | Lor3D ◆ | Lorentzian | 0.4430 |
| 4 | KR_2layer ● | KR-family | 0.6792 |
| 5 | KR_like ● | KR-family | 0.8999 |
| 6 | MLR | Layered | 1.1927 |
| 7 | KR_4layer ● | KR-family | 1.5556 |
| 8 | TransPerc | Other | 1.5645 |
| 9 | RLk4 | Layered | 1.9834 |
| 10 | RLk6_mid | Layered | 2.3692 |
| 11 | Lor2D ◆ | Lorentzian | 2.5001 |
| 12 | RLk6_tap | Layered | 2.6969 |
| 13 | IntOrder | Other | 3.2633 |
| 14 | RLk6 | Layered | 3.5447 |
| 15 | RLk6_lj | Layered | 3.5848 |
| 16 | RLk8 | Layered | 4.1366 |
| 17 | AbsLayer | Layered | 4.4840 |


## 4. Baseline: Constant-Center LSD-W2

For comparison: use N=48 centroid as fixed well center.

Fixed centers: c* = 0.2078, w* = 0.4167, d* = 4

**Best**: α=0.5, β=1.0, γ=1.0
**Lor4D beats 100.0% of non-Lor across all N**

| N | Lor4D Rank |
|---|:----------:|
| 16 | #1/17 |
| 20 | #1/17 |
| 28 | #1/17 |
| 36 | #1/17 |
| 48 | #1/17 |
| 64 | #1/17 |


## 5. Head-to-Head Comparison

| Mode | Description | Beat% | Worst Rank |
|------|------------|:-----:|:----------:|
| A (Oracle) | c*(N), w*(N) from same-N | 100.0% | #1/17 |
| B (Extrap) | c*(N), w*(N) power-law fit N≤36 | 100.0% | #1/17 |
| C (Const)  | Fixed c*, w* from N=48 | 100.0% | #1/17 |

### Per-N Rank Comparison

| N | A (Oracle) | B (Extrap) | C (Const) |
|---|:----------:|:----------:|:---------:|
| 16 | #1 | #1 | #1 |
| 20 | #1 | #1 | #1 |
| 28 | #1 | #1 | #1 |
| 36 | #1 | #1 | #1 |
| 48 | #1 | #1 | #1 |
| 64 | #1 | #1 | #1 |

### Lor4D vs KR_2layer (Oracle mode)

| N | Lor4D F | KR_2layer F | Gap | Win? |
|---|:------:|:----------:|:---:|:----:|
| 16 | 0.1277 | 0.1278 | +0.0001 | ✅ |
| 20 | 0.0569 | 0.3852 | +0.3283 | ✅ |
| 28 | 0.0600 | 0.3777 | +0.3177 | ✅ |
| 36 | 0.0398 | 0.5489 | +0.5090 | ✅ |
| 48 | 0.0289 | 0.6106 | +0.5817 | ✅ |
| 64 | 0.0239 | 0.7901 | +0.7662 | ✅ |


## 6. Conclusion

N-adapted wells match constant-center performance.

Extrapolated mode matches or exceeds constant baseline → power-law scaling is viable for prediction.

**Final form:**

```
F_LSD(N) = α·(d_eff − 4)² + β·(C₁/C₀ − c*(N))² + γ·(w − w*(N))²
  d* = 4 (fixed)
  c*(N) = 0.2485 + -2.33/N
  w*(N) = 0.3255 + 3.80/N
  Best weights (oracle): α=0.5, β=1.0, γ=5.0
```