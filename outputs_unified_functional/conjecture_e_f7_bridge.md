# Conjecture E: F7↔S_BD Bridge (sigmoid wall upgrade)

## F7 main model (§5.10.7)

F7 = logH + 0.0004·Π_geo - 10·Σ_hist + 0.6·Ξ_d + α(N)·σ((R-Rc)/w)

α₀=16, q=-0.5, Rc=0.25, w=0.015, N₀=20

## F7 vs F5 correlation

Spearman(F7, F5) = +0.7539

## F7 vs BD observables

| pair | Spearman |
|---|---:|
| F7 vs bd_ratio | +0.2192 |
| F7 vs R | +0.1804 |
| wall vs R | +0.8215 |
| wall vs bd_ratio | +0.7903 |

## Residual bridge: F7 ~ N + family (+BD)

Baseline R² = 0.9373

| metric | R² ext | ΔR² | coef | ρ(resid) |
|---|---:|---:|---:|---:|
| bd_ratio | 0.9589 | +0.0216 | -8.8147 | +0.0968 |
| bdg_d2_corrected_norm | 0.9391 | +0.0018 | -0.7104 | -0.1370 |
| bdg_d4_standard_norm | 0.9373 | +0.0000 | -0.0385 | -0.0120 |
| R | 0.9373 | +0.0000 | +1.2941 | +0.1361 |

## F7 ordering: Lor4D vs KR_like

| N | F7 win | F5 win |
|---:|---:|---:|
| 16 | 1.000 | 0.000 |
| 20 | 1.000 | 0.000 |
| 28 | 1.000 | 0.000 |
| 36 | 0.875 | 0.000 |

## ρ_BD local density additivity

Samples with patches: 139/160

| family | add_ratio | std | n |
|---|---:|---:|---:|
| Lor2D | 0.092 | 0.033 | 32 |
| Lor3D | 0.222 | 0.076 | 31 |
| Lor4D | 0.324 | 0.027 | 27 |
| Lor5D | 0.333 | 0.000 | 17 |
| KR_like | 0.333 | 0.000 | 32 |

## Per-family means

| family | F7 | F5 | R | wall | bd_ratio |
|---|---:|---:|---:|---:|---:|
| Lor2D | 41.53 | 60.14 | 0.6481 | 14.83 | 2.0966 |
| Lor3D | 49.91 | 88.83 | 0.3382 | 11.39 | 0.5806 |
| Lor4D | 44.87 | 111.86 | 0.1199 | 0.43 | 0.1452 |
| Lor5D | 48.71 | 126.91 | 0.0543 | 0.80 | 0.0632 |
| KR_like | 54.49 | 82.81 | 0.3269 | 14.59 | 0.6681 |

## Interpretation

The F7↔S_BD bridge tests whether the sigmoid wall (which uses R = 1 - f_link) 
is already capturing the same information as the BD action's interval statistics. 
A high Spearman(wall, bd_ratio) confirms that the sigmoid wall in F7 is a 
monotone proxy for the same causal interval richness that BD actions encode.

The ρ_BD local density additivity check tests whether bd_ratio can be 
decomposed into local layer-block patches — a prerequisite for writing it as 
ρ_BD(x) in the continuum limit.
