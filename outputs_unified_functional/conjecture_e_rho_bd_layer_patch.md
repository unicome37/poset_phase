# Conjecture E: ρ_BD Layer Block Patch Locality Test

## Method

Partition poset into consecutive layer blocks of size `k` (=2,3,4 layers),
compute bd_ratio on each block's induced sub-poset, test additivity:

  additivity_ratio = (node-weighted mean of local bd_ratios) / global bd_ratio

A ratio near 1.0 means bd_ratio decomposes into local densities ρ_BD(x).

## Block size = 2

Valid: 139/160, mean=0.0000, median=0.0000, std=0.0000

| family | add_ratio | std | n | blocks | coverage |
|---|---:|---:|---:|---:|---:|
| Lor2D | 0.0000 | 0.0000 | 32 | 3.7 | 0.989 |
| Lor3D | 0.0000 | 0.0000 | 31 | 2.2 | 0.997 |
| Lor4D | 0.0000 | 0.0000 | 27 | 2.0 | 0.998 |
| Lor5D | 0.0000 | 0.0000 | 17 | 2.0 | 0.987 |
| KR_like | 0.0000 | 0.0000 | 32 | 2.0 | 1.000 |

| N | add_ratio | std | n |
|---:|---:|---:|---:|
| 16 | 0.0000 | 0.0000 | 29 |
| 20 | 0.0000 | 0.0000 | 34 |
| 28 | 0.0000 | 0.0000 | 37 |
| 36 | 0.0000 | 0.0000 | 39 |

## Block size = 3

Valid: 139/160, mean=0.7133, median=1.0000, std=0.3582

| family | add_ratio | std | n | blocks | coverage |
|---|---:|---:|---:|---:|---:|
| Lor2D | 0.1903 | 0.0649 | 32 | 2.6 | 0.991 |
| Lor3D | 0.6049 | 0.2854 | 31 | 1.7 | 0.993 |
| Lor4D | 0.9586 | 0.1221 | 27 | 1.1 | 0.998 |
| Lor5D | 0.9664 | 0.1343 | 17 | 1.1 | 0.998 |
| KR_like | 1.0000 | 0.0000 | 32 | 1.0 | 1.000 |

| N | add_ratio | std | n |
|---:|---:|---:|---:|
| 16 | 0.7534 | 0.3363 | 29 |
| 20 | 0.7276 | 0.3404 | 34 |
| 28 | 0.7279 | 0.3750 | 37 |
| 36 | 0.6573 | 0.3661 | 39 |

## Block size = 4

Valid: 139/160, mean=0.8400, median=1.0000, std=0.2811

| family | add_ratio | std | n | blocks | coverage |
|---|---:|---:|---:|---:|---:|
| Lor2D | 0.3775 | 0.1878 | 32 | 2.1 | 0.995 |
| Lor3D | 0.9254 | 0.1641 | 31 | 1.2 | 0.997 |
| Lor4D | 1.0000 | 0.0000 | 27 | 1.0 | 1.000 |
| Lor5D | 1.0000 | 0.0000 | 17 | 1.0 | 1.000 |
| KR_like | 1.0000 | 0.0000 | 32 | 1.0 | 1.000 |

| N | add_ratio | std | n |
|---:|---:|---:|---:|
| 16 | 0.8815 | 0.2185 | 29 |
| 20 | 0.8527 | 0.2582 | 34 |
| 28 | 0.8262 | 0.3047 | 37 |
| 36 | 0.8113 | 0.3122 | 39 |

## Comparison with Alexandrov interval baseline

| method | mean add_ratio |
|---|---:|
| Alexandrov interval | 0.117 |
| **Layer block (bs=4)** | **0.8400** |

## Interpretation

If layer blocks achieve additivity_ratio >> 0.12 (the Alexandrov baseline),
this validates that bd_ratio can be decomposed into a local density ρ_BD(x)
defined on layer patches — a necessary step for the continuum limit.

Perfect additivity (ratio=1) is not expected because boundary effects
between blocks create extra/missing pairs. The key test is whether the
ratio is stable across N (finite-size scaling) and across families.
