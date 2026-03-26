# LSD-Well: Lorentzian Structural Discriminator with Well Potentials

**No logH. Wells in interval-shape and transverse-structure space.**

N = [16, 20, 28, 36, 48, 64], reps = 15


## LSD-W1: α·(C₁/C₀ − c*)² + β·(width − w*)²

**Fixed well centers from Lor4D large-N centroid.**

Reference (N=48 centroid): c*=0.2127, w*=0.4083, d*=3.94

**Best**: c*=0.1701, w*=0.4083, α=10.0, β=0.5
**Lor4D beats 97.4% of non-Lor across all N**

| N | Lor4D Rank |
|---|:----------:|
| 16 | #2/17 |
| 20 | #2/17 |
| 28 | #1/17 |
| 36 | #1/17 |
| 48 | #1/17 |
| 64 | #1/17 |

### LSD-W1 N=16

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | KR_like ● | KR-family | 0.0601 |
| 2 | Lor4D ◆ | Lorentzian | 0.1720 |
| 3 | Lor3D ◆ | Lorentzian | 0.2410 |
| 4 | Lor5D ◆ | Lorentzian | 0.2734 |
| 5 | KR_4layer ● | KR-family | 0.3167 |
| 6 | KR_2layer ● | KR-family | 0.3282 |
| 7 | TransPerc | Other | 0.5302 |
| 8 | IntOrder | Other | 1.3253 |
| 9 | RLk6_tap | Layered | 1.7422 |
| 10 | MLR | Layered | 1.8157 |
| 11 | RLk8 | Layered | 1.8768 |
| 12 | RLk6_lj | Layered | 1.9712 |
| 13 | Lor2D ◆ | Lorentzian | 2.0334 |
| 14 | RLk4 | Layered | 2.2616 |
| 15 | RLk6 | Layered | 2.3266 |
| 16 | RLk6_mid | Layered | 2.3474 |
| 17 | AbsLayer | Layered | 6.5574 |

### LSD-W1 N=20

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | KR_4layer ● | KR-family | 0.0506 |
| 2 | Lor4D ◆ | Lorentzian | 0.0913 |
| 3 | KR_like ● | KR-family | 0.0940 |
| 4 | Lor3D ◆ | Lorentzian | 0.2077 |
| 5 | Lor5D ◆ | Lorentzian | 0.2481 |
| 6 | KR_2layer ● | KR-family | 0.3425 |
| 7 | TransPerc | Other | 0.8005 |
| 8 | MLR | Layered | 1.5852 |
| 9 | IntOrder | Other | 2.1379 |
| 10 | Lor2D ◆ | Lorentzian | 2.3211 |
| 11 | RLk4 | Layered | 2.4102 |
| 12 | RLk6 | Layered | 2.7497 |
| 13 | RLk6_lj | Layered | 2.9168 |
| 14 | RLk6_mid | Layered | 2.9172 |
| 15 | RLk6_tap | Layered | 2.9560 |
| 16 | RLk8 | Layered | 3.3788 |
| 17 | AbsLayer | Layered | 5.5302 |

### LSD-W1 N=28

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0629 |
| 2 | KR_4layer ● | KR-family | 0.0807 |
| 3 | Lor5D ◆ | Lorentzian | 0.1668 |
| 4 | KR_like ● | KR-family | 0.1782 |
| 5 | KR_2layer ● | KR-family | 0.3470 |
| 6 | Lor3D ◆ | Lorentzian | 0.4406 |
| 7 | IntOrder | Other | 1.3860 |
| 8 | Lor2D ◆ | Lorentzian | 2.0327 |
| 9 | MLR | Layered | 2.4172 |
| 10 | TransPerc | Other | 2.4609 |
| 11 | RLk6_mid | Layered | 3.5777 |
| 12 | RLk4 | Layered | 4.0102 |
| 13 | RLk6 | Layered | 4.0783 |
| 14 | RLk6_tap | Layered | 4.2376 |
| 15 | RLk6_lj | Layered | 4.4372 |
| 16 | AbsLayer | Layered | 4.6478 |
| 17 | RLk8 | Layered | 4.8023 |

### LSD-W1 N=36

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0532 |
| 2 | Lor5D ◆ | Lorentzian | 0.1155 |
| 3 | KR_4layer ● | KR-family | 0.1330 |
| 4 | KR_like ● | KR-family | 0.2318 |
| 5 | AbsLayer | Layered | 0.2936 |
| 6 | KR_2layer ● | KR-family | 0.3478 |
| 7 | Lor3D ◆ | Lorentzian | 0.5150 |
| 8 | IntOrder | Other | 1.2061 |
| 9 | Lor2D ◆ | Lorentzian | 1.9409 |
| 10 | RLk4 | Layered | 2.9433 |
| 11 | TransPerc | Other | 3.5936 |
| 12 | MLR | Layered | 3.6208 |
| 13 | RLk6_mid | Layered | 3.6279 |
| 14 | RLk6 | Layered | 4.9215 |
| 15 | RLk6_lj | Layered | 5.1093 |
| 16 | RLk6_tap | Layered | 5.3232 |
| 17 | RLk8 | Layered | 6.3721 |

### LSD-W1 N=48

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0877 |
| 2 | Lor5D ◆ | Lorentzian | 0.0923 |
| 3 | KR_4layer ● | KR-family | 0.2061 |
| 4 | KR_like ● | KR-family | 0.2785 |
| 5 | KR_2layer ● | KR-family | 0.3478 |
| 6 | Lor3D ◆ | Lorentzian | 0.5118 |
| 7 | RLk4 | Layered | 1.3879 |
| 8 | IntOrder | Other | 1.5502 |
| 9 | Lor2D ◆ | Lorentzian | 2.4701 |
| 10 | RLk6_mid | Layered | 3.1648 |
| 11 | MLR | Layered | 3.2582 |
| 12 | RLk6_tap | Layered | 4.4534 |
| 13 | RLk6 | Layered | 5.3296 |
| 14 | TransPerc | Other | 6.1653 |
| 15 | RLk6_lj | Layered | 6.3107 |
| 16 | RLk8 | Layered | 6.6612 |
| 17 | AbsLayer | Layered | 8.7928 |

### LSD-W1 N=64

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0903 |
| 2 | Lor5D ◆ | Lorentzian | 0.1018 |
| 3 | KR_4layer ● | KR-family | 0.2781 |
| 4 | KR_like ● | KR-family | 0.2928 |
| 5 | AbsLayer | Layered | 0.3012 |
| 6 | KR_2layer ● | KR-family | 0.3478 |
| 7 | RLk4 | Layered | 0.5978 |
| 8 | Lor3D ◆ | Lorentzian | 0.7856 |
| 9 | IntOrder | Other | 1.9239 |
| 10 | RLk6_mid | Layered | 2.4183 |
| 11 | Lor2D ◆ | Lorentzian | 2.8108 |
| 12 | RLk6_tap | Layered | 3.5503 |
| 13 | MLR | Layered | 4.4605 |
| 14 | RLk6 | Layered | 4.9512 |
| 15 | RLk6_lj | Layered | 5.6198 |
| 16 | RLk8 | Layered | 7.9097 |
| 17 | TransPerc | Other | 9.2254 |


## LSD-W2: α·(C₁/C₀−c*)² + β·(w−w*)² + γ·(d_eff−d*)²

**Triple well: interval shape + width + dimension.**

**Best**: c*=0.2127, w*=0.4083, d*=3.75, α=1.0, β=1.0, γ=0.5
**Lor4D beats 100.0% of non-Lor across all N**

| N | Lor4D Rank |
|---|:----------:|
| 16 | #1/17 |
| 20 | #1/17 |
| 28 | #1/17 |
| 36 | #1/17 |
| 48 | #1/17 |
| 64 | #1/17 |

### LSD-W2 N=16

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.1174 |
| 2 | KR_2layer ● | KR-family | 0.1434 |
| 3 | Lor3D ◆ | Lorentzian | 0.2749 |
| 4 | Lor5D ◆ | Lorentzian | 0.3431 |
| 5 | RLk6_tap | Layered | 0.3459 |
| 6 | RLk8 | Layered | 0.4034 |
| 7 | TransPerc | Other | 0.4261 |
| 8 | RLk6_mid | Layered | 0.4485 |
| 9 | RLk6_lj | Layered | 0.4527 |
| 10 | MLR | Layered | 0.4970 |
| 11 | KR_like ● | KR-family | 0.5754 |
| 12 | RLk6 | Layered | 0.7732 |
| 13 | RLk4 | Layered | 0.8114 |
| 14 | KR_4layer ● | KR-family | 1.1542 |
| 15 | Lor2D ◆ | Lorentzian | 1.8617 |
| 16 | IntOrder | Other | 2.7421 |
| 17 | AbsLayer | Layered | 4.7891 |

### LSD-W2 N=20

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.1159 |
| 2 | Lor3D ◆ | Lorentzian | 0.1435 |
| 3 | KR_2layer ● | KR-family | 0.1612 |
| 4 | TransPerc | Other | 0.2634 |
| 5 | Lor5D ◆ | Lorentzian | 0.3941 |
| 6 | MLR | Layered | 0.4529 |
| 7 | RLk6 | Layered | 0.5442 |
| 8 | RLk6_tap | Layered | 0.5533 |
| 9 | RLk8 | Layered | 0.6264 |
| 10 | RLk6_mid | Layered | 0.6385 |
| 11 | KR_like ● | KR-family | 0.6386 |
| 12 | RLk6_lj | Layered | 0.7987 |
| 13 | RLk4 | Layered | 0.8003 |
| 14 | KR_4layer ● | KR-family | 0.9357 |
| 15 | Lor2D ◆ | Lorentzian | 1.9699 |
| 16 | IntOrder | Other | 2.5335 |
| 17 | AbsLayer | Layered | 4.4365 |

### LSD-W2 N=28

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0745 |
| 2 | KR_2layer ● | KR-family | 0.1667 |
| 3 | Lor3D ◆ | Lorentzian | 0.2050 |
| 4 | Lor5D ◆ | Lorentzian | 0.2470 |
| 5 | TransPerc | Other | 0.3220 |
| 6 | KR_like ● | KR-family | 0.5894 |
| 7 | MLR | Layered | 0.6680 |
| 8 | RLk6 | Layered | 0.9673 |
| 9 | RLk6_mid | Layered | 0.9756 |
| 10 | KR_4layer ● | KR-family | 0.9851 |
| 11 | RLk8 | Layered | 0.9867 |
| 12 | RLk6_tap | Layered | 1.0107 |
| 13 | RLk4 | Layered | 1.0758 |
| 14 | RLk6_lj | Layered | 1.3144 |
| 15 | Lor2D ◆ | Lorentzian | 1.7055 |
| 16 | IntOrder | Other | 2.5613 |
| 17 | AbsLayer | Layered | 4.2529 |

### LSD-W2 N=36

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0630 |
| 2 | KR_2layer ● | KR-family | 0.1693 |
| 3 | Lor3D ◆ | Lorentzian | 0.2072 |
| 4 | Lor5D ◆ | Lorentzian | 0.2537 |
| 5 | TransPerc | Other | 0.3649 |
| 6 | KR_like ● | KR-family | 0.6043 |
| 7 | MLR | Layered | 0.9117 |
| 8 | KR_4layer ● | KR-family | 1.0088 |
| 9 | RLk6_mid | Layered | 1.3212 |
| 10 | RLk4 | Layered | 1.3311 |
| 11 | Lor2D ◆ | Lorentzian | 1.4389 |
| 12 | RLk6_tap | Layered | 1.6491 |
| 13 | RLk6 | Layered | 1.7542 |
| 14 | RLk8 | Layered | 1.7941 |
| 15 | RLk6_lj | Layered | 2.0062 |
| 16 | IntOrder | Other | 2.5550 |
| 17 | AbsLayer | Layered | 3.8325 |

### LSD-W2 N=48

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0423 |
| 2 | KR_2layer ● | KR-family | 0.1707 |
| 3 | Lor3D ◆ | Lorentzian | 0.1711 |
| 4 | Lor5D ◆ | Lorentzian | 0.2531 |
| 5 | KR_like ● | KR-family | 0.5785 |
| 6 | TransPerc | Other | 0.6096 |
| 7 | MLR | Layered | 0.9562 |
| 8 | KR_4layer ● | KR-family | 1.1960 |
| 9 | RLk4 | Layered | 1.3503 |
| 10 | RLk6_mid | Layered | 1.6951 |
| 11 | Lor2D ◆ | Lorentzian | 1.8449 |
| 12 | RLk6_tap | Layered | 1.9112 |
| 13 | RLk8 | Layered | 2.2006 |
| 14 | RLk6_lj | Layered | 2.3953 |
| 15 | RLk6 | Layered | 2.3991 |
| 16 | IntOrder | Other | 2.6357 |
| 17 | AbsLayer | Layered | 3.8074 |

### LSD-W2 N=64

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.0238 |
| 2 | KR_2layer ● | KR-family | 0.1677 |
| 3 | Lor3D ◆ | Lorentzian | 0.2272 |
| 4 | Lor5D ◆ | Lorentzian | 0.2604 |
| 5 | KR_like ● | KR-family | 0.5948 |
| 6 | MLR | Layered | 1.1625 |
| 7 | KR_4layer ● | KR-family | 1.1962 |
| 8 | TransPerc | Other | 1.2407 |
| 9 | RLk4 | Layered | 1.4989 |
| 10 | Lor2D ◆ | Lorentzian | 1.7017 |
| 11 | RLk6_mid | Layered | 1.7692 |
| 12 | RLk6_tap | Layered | 2.2546 |
| 13 | IntOrder | Other | 2.5234 |
| 14 | RLk6_lj | Layered | 2.9147 |
| 15 | RLk6 | Layered | 2.9334 |
| 16 | RLk8 | Layered | 3.1950 |
| 17 | AbsLayer | Layered | 3.3626 |


## LSD-W3: N·[α·(C₁/C₀−c*)² + β·(w−w*)² + γ·(d_eff−d*)²]

**N-scaled well: O(N) penalty ensures wells compete at all scales.**

**Best**: c*=0.2127, w*=0.4083, d*=3.75, α=0.5, β=0.5, γ=0.2
**Lor4D beats 100.0% of non-Lor across all N**

| N | Lor4D Rank |
|---|:----------:|
| 16 | #1/17 |
| 20 | #1/17 |
| 28 | #1/17 |
| 36 | #1/17 |
| 48 | #1/17 |
| 64 | #1/17 |

### LSD-W3 N=16

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.8399 |
| 2 | KR_2layer ● | KR-family | 1.1141 |
| 3 | Lor3D ◆ | Lorentzian | 1.8004 |
| 4 | Lor5D ◆ | Lorentzian | 2.3616 |
| 5 | RLk6_tap | Layered | 2.4493 |
| 6 | RLk8 | Layered | 2.8436 |
| 7 | TransPerc | Other | 2.9200 |
| 8 | RLk6_lj | Layered | 3.1700 |
| 9 | RLk6_mid | Layered | 3.1998 |
| 10 | MLR | Layered | 3.4364 |
| 11 | KR_like ● | KR-family | 3.7067 |
| 12 | RLk6 | Layered | 5.2792 |
| 13 | RLk4 | Layered | 5.5028 |
| 14 | KR_4layer ● | KR-family | 7.4254 |
| 15 | Lor2D ◆ | Lorentzian | 12.2083 |
| 16 | IntOrder | Other | 17.7524 |
| 17 | AbsLayer | Layered | 31.6686 |

### LSD-W3 N=20

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.9959 |
| 2 | Lor3D ◆ | Lorentzian | 1.1854 |
| 3 | KR_2layer ● | KR-family | 1.5921 |
| 4 | TransPerc | Other | 2.2604 |
| 5 | Lor5D ◆ | Lorentzian | 3.3765 |
| 6 | MLR | Layered | 3.8965 |
| 7 | RLk6 | Layered | 4.8433 |
| 8 | RLk6_tap | Layered | 4.9447 |
| 9 | KR_like ● | KR-family | 5.1570 |
| 10 | RLk8 | Layered | 5.6204 |
| 11 | RLk6_mid | Layered | 5.6272 |
| 12 | RLk4 | Layered | 6.8227 |
| 13 | RLk6_lj | Layered | 6.9110 |
| 14 | KR_4layer ● | KR-family | 7.5034 |
| 15 | Lor2D ◆ | Lorentzian | 16.1969 |
| 16 | IntOrder | Other | 20.6633 |
| 17 | AbsLayer | Layered | 36.5973 |

### LSD-W3 N=28

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.8848 |
| 2 | KR_2layer ● | KR-family | 2.3165 |
| 3 | Lor3D ◆ | Lorentzian | 2.3953 |
| 4 | Lor5D ◆ | Lorentzian | 2.9211 |
| 5 | TransPerc | Other | 4.2103 |
| 6 | KR_like ● | KR-family | 6.7083 |
| 7 | MLR | Layered | 8.0734 |
| 8 | KR_4layer ● | KR-family | 11.0774 |
| 9 | RLk6_mid | Layered | 11.8415 |
| 10 | RLk6 | Layered | 11.8832 |
| 11 | RLk8 | Layered | 12.3039 |
| 12 | RLk6_tap | Layered | 12.3856 |
| 13 | RLk4 | Layered | 13.0565 |
| 14 | RLk6_lj | Layered | 15.8697 |
| 15 | Lor2D ◆ | Lorentzian | 19.6749 |
| 16 | IntOrder | Other | 29.0805 |
| 17 | AbsLayer | Layered | 48.9608 |

### LSD-W3 N=36

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.9428 |
| 2 | KR_2layer ● | KR-family | 3.0215 |
| 3 | Lor3D ◆ | Lorentzian | 3.1496 |
| 4 | Lor5D ◆ | Lorentzian | 3.7938 |
| 5 | TransPerc | Other | 6.4014 |
| 6 | KR_like ● | KR-family | 8.8659 |
| 7 | MLR | Layered | 14.2928 |
| 8 | KR_4layer ● | KR-family | 14.6171 |
| 9 | RLk4 | Layered | 20.1040 |
| 10 | RLk6_mid | Layered | 20.2360 |
| 11 | Lor2D ◆ | Lorentzian | 21.4400 |
| 12 | RLk6_tap | Layered | 25.4960 |
| 13 | RLk6 | Layered | 26.9337 |
| 14 | RLk8 | Layered | 28.0208 |
| 15 | RLk6_lj | Layered | 30.6157 |
| 16 | IntOrder | Other | 37.2435 |
| 17 | AbsLayer | Layered | 55.3812 |

### LSD-W3 N=48

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.8561 |
| 2 | Lor3D ◆ | Lorentzian | 3.5210 |
| 3 | KR_2layer ● | KR-family | 4.0544 |
| 4 | Lor5D ◆ | Lorentzian | 4.9858 |
| 5 | KR_like ● | KR-family | 11.3548 |
| 6 | TransPerc | Other | 14.4544 |
| 7 | MLR | Layered | 19.7402 |
| 8 | KR_4layer ● | KR-family | 23.1331 |
| 9 | RLk4 | Layered | 26.5032 |
| 10 | RLk6_mid | Layered | 33.9287 |
| 11 | Lor2D ◆ | Lorentzian | 36.7075 |
| 12 | RLk6_tap | Layered | 38.6281 |
| 13 | RLk8 | Layered | 45.3212 |
| 14 | RLk6 | Layered | 48.4821 |
| 15 | RLk6_lj | Layered | 48.8535 |
| 16 | IntOrder | Other | 51.4054 |
| 17 | AbsLayer | Layered | 77.4319 |

### LSD-W3 N=64

| Rank | Family | Category | F |
|------|--------|----------|:-:|
| 1 | Lor4D ◆ | Lorentzian | 0.6584 |
| 2 | KR_2layer ● | KR-family | 5.3299 |
| 3 | Lor3D ◆ | Lorentzian | 6.3414 |
| 4 | Lor5D ◆ | Lorentzian | 6.8423 |
| 5 | KR_like ● | KR-family | 15.5692 |
| 6 | KR_4layer ● | KR-family | 30.9093 |
| 7 | MLR | Layered | 32.3703 |
| 8 | TransPerc | Other | 37.4274 |
| 9 | RLk4 | Layered | 38.7275 |
| 10 | Lor2D ◆ | Lorentzian | 45.5348 |
| 11 | RLk6_mid | Layered | 46.6739 |
| 12 | RLk6_tap | Layered | 59.7792 |
| 13 | IntOrder | Other | 65.9482 |
| 14 | RLk6_lj | Layered | 78.0002 |
| 15 | RLk6 | Layered | 78.0933 |
| 16 | AbsLayer | Layered | 86.5236 |
| 17 | RLk8 | Layered | 86.6881 |


## Overall Winner: LSD-W2 (100.0%)

### Lor4D vs KR_2layer across N

| N | Lor4D F | KR_2layer F | Lor4D wins? |
|---|:-------:|:-----------:|:-----------:|
| 16 | 0.1174 | 0.1434 | ✅ |
| 20 | 0.1159 | 0.1612 | ✅ |
| 28 | 0.0745 | 0.1667 | ✅ |
| 36 | 0.0630 | 0.1693 | ✅ |
| 48 | 0.0423 | 0.1707 | ✅ |
| 64 | 0.0238 | 0.1677 | ✅ |

### Lor4D vs ALL at N=64

| Competitor | Category | Lor4D F | Comp F | Win? |
|------------|----------|:-------:|:------:|:----:|
| AbsLayer | Layered | 0.0238 | 3.3626 | ✅ |
| IntOrder | Other | 0.0238 | 2.5234 | ✅ |
| KR_2layer | KR-family | 0.0238 | 0.1677 | ✅ |
| KR_4layer | KR-family | 0.0238 | 1.1962 | ✅ |
| KR_like | KR-family | 0.0238 | 0.5948 | ✅ |
| Lor2D | Lorentzian | 0.0238 | 1.7017 | ✅ |
| Lor3D | Lorentzian | 0.0238 | 0.2272 | ✅ |
| Lor5D | Lorentzian | 0.0238 | 0.2604 | ✅ |
| MLR | Layered | 0.0238 | 1.1625 | ✅ |
| RLk4 | Layered | 0.0238 | 1.4989 | ✅ |
| RLk6 | Layered | 0.0238 | 2.9334 | ✅ |
| RLk6_lj | Layered | 0.0238 | 2.9147 | ✅ |
| RLk6_mid | Layered | 0.0238 | 1.7692 | ✅ |
| RLk6_tap | Layered | 0.0238 | 2.2546 | ✅ |
| RLk8 | Layered | 0.0238 | 3.1950 | ✅ |
| TransPerc | Other | 0.0238 | 1.2407 | ✅ |


## Comparison: All Generations

| Functional | Architecture | Lor4D best rank | N=64 rank | KR_2layer? | Score |
|------------|-------------|:---------------:|:---------:|:----------:|:-----:|
| F7 | logH + wall | #1 (N=16) | #8+ | ❌ N≥28 | ~60% |
| F10 | logH + d_eff-well | #1 (all N) | #1 | ✅ | ~95% |
| F_link | S_link + wall | #3 (N=16) | #8 | ❌ always | ~40% |
| LSD2 | C₁/C₀ + wall | #1 (N≤28) | #10 | ❌ N≥48 | 89.7% |
| LSD-W1 | shape+width well | ? | #1 | ? | 97.4% |
| LSD-W2 | triple well | ? | #1 | ? | 100.0% |
| LSD-W3 | N·triple well | ? | #1 | ? | 100.0% |