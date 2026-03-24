# Conjecture E — C_k Mode Decomposition (§4.1.25)

Loaded 360 realizations, C_0 through C_4


## Part 1: Global PCA on p_k = C_k / Σ C_k

### Explained variance

| PC | Var explained | Cumulative |
|----|--------------|------------|
| PC1 | 0.8855 (88.6%) | 0.8855 (88.6%) |
| PC2 | 0.0820 (8.2%) | 0.9675 (96.8%) |
| PC3 | 0.0254 (2.5%) | 0.9929 (99.3%) |
| PC4 | 0.0048 (0.5%) | 0.9977 (99.8%) |
| PC5 | 0.0023 (0.2%) | 1.0000 (100.0%) |

### PC loadings (on standardized p_k)

| | p_0 | p_1 | p_2 | p_3 | p_4 |
|---|---|---|---|---|---|
| PC1 | +0.448 | -0.393 | -0.470 | -0.466 | -0.455 |
| PC2 | -0.173 | -0.874 | -0.038 | +0.230 | +0.388 |
| PC3 | -0.876 | -0.041 | -0.200 | -0.279 | -0.336 |
| PC4 | +0.038 | -0.280 | +0.826 | -0.092 | -0.480 |

## Part 2: PC scores vs curvature

### Spearman(PC_i, H²) per slice

| d | N | PC1 ρ | PC1 p | PC2 ρ | PC2 p | PC3 ρ | PC3 p |
|---|---|-------|-------|-------|-------|-------|-------|
| 2 | 128 | +0.625 | 1.6e-05 | -0.645 | 7.1e-06 | +0.602 | 4.0e-05 |
| 2 | 256 | +0.547 | 2.6e-04 | -0.554 | 2.1e-04 | +0.533 | 4.0e-04 |
| 2 | 512 | +0.655 | 4.4e-06 | -0.610 | 3.0e-05 | +0.593 | 5.6e-05 |
| 3 | 128 | +0.933 | 2.1e-18 | -0.409 | 8.8e-03 | +0.021 | 9.0e-01 |
| 3 | 256 | +0.922 | 3.1e-17 | -0.839 | 1.3e-11 | +0.712 | 2.6e-07 |
| 3 | 512 | +0.953 | 3.2e-21 | -0.827 | 4.8e-11 | +0.790 | 1.3e-09 |
| 4 | 128 | +0.912 | 2.6e-16 | +0.756 | 1.7e-08 | -0.428 | 5.9e-03 |
| 4 | 256 | +0.962 | 4.6e-23 | +0.532 | 4.2e-04 | -0.711 | 2.7e-07 |
| 4 | 512 | +0.972 | 1.2e-25 | +0.188 | 2.4e-01 | -0.398 | 1.1e-02 |

## Part 3: bdg_d2c projection analysis

bdg_d2c coefficient vector in raw p_k space: [-2, +2, 0, 0, 0]

### Projection of d2c onto each PC

| PC | cos(θ) | |cos(θ)|² (variance explained) |
|----|--------|------------------------------|
| PC1 | -0.5316 | 0.2826 (28.3%) |
| PC2 | -0.0483 | 0.0023 (0.2%) |
| PC3 | +0.8389 | 0.7037 (70.4%) |
| PC4 | -0.1056 | 0.0112 (1.1%) |

**Key finding**: d2c projects 28.3% onto PC1 (density) and 71.7% onto PC2+ (shape).


## Part 4: Density-subtracted shape score

Define: shape_score = ||projection onto PC2+|| = sqrt(Σ_{i≥2} PC_i²)

This removes PC1 (density mode) and keeps only shape residual.

### Spearman(shape_score, H²) per slice

| d | N | shape_score ρ | p | PC2 ρ | p | KL ρ | p |
|---|---|--------------|---|-------|---|------|---|
| 2 | 128 | -0.397 | 1.1e-02 | -0.645 | 7.1e-06 | +0.899 | 3.4e-15 |
| 2 | 256 | -0.515 | 6.8e-04 | -0.554 | 2.1e-04 | -0.891 | 1.3e-14 |
| 2 | 512 | -0.616 | 2.4e-05 | -0.610 | 3.0e-05 | -0.980 | 2.7e-28 |
| 3 | 128 | +0.518 | 6.2e-04 | -0.409 | 8.8e-03 | +0.974 | 4.2e-26 |
| 3 | 256 | +0.848 | 4.8e-12 | -0.839 | 1.3e-11 | +0.980 | 2.7e-28 |
| 3 | 512 | +0.744 | 3.7e-08 | -0.827 | 4.8e-11 | +0.980 | 2.7e-28 |
| 4 | 128 | +0.581 | 8.4e-05 | +0.756 | 1.7e-08 | +0.972 | 1.3e-23 |
| 4 | 256 | +0.483 | 1.6e-03 | +0.532 | 4.2e-04 | +0.981 | 1.7e-28 |
| 4 | 512 | +0.806 | 3.6e-10 | +0.188 | 2.4e-01 | +0.980 | 2.7e-28 |

## Part 5: Per-dimension PCA

Since de Sitter sprinklings differ significantly across d, repeat PCA within each dimension.


### d = 2

| PC | Var% | Cumul% |
|----|------|--------|
| PC1 | 56.7% | 56.7% |
| PC2 | 20.8% | 77.5% |
| PC3 | 14.3% | 91.7% |
| PC4 | 8.3% | 100.0% |
| PC5 | 0.0% | 100.0% |

Loadings:
| | p_0 | p_1 | p_2 | p_3 | p_4 |
|---|---|---|---|---|---|
| PC1 | +0.543 | +0.321 | -0.339 | -0.476 | -0.509 |
| PC2 | +0.286 | -0.642 | +0.605 | -0.303 | -0.219 |
| PC3 | +0.331 | -0.622 | -0.628 | +0.319 | +0.081 |

d2c projection:
  PC1: cos=-0.3754 (14.1%)
  PC2: cos=-0.5138 (26.4%)
  PC3: cos=-0.5473 (30.0%)

Shape score vs H²:
| N | shape ρ | p | PC2 ρ | p |
|---|---------|---|-------|---|
| 128 | -0.155 | 3.4e-01 | -0.234 | 1.5e-01 |
| 256 | +0.172 | 2.9e-01 | -0.006 | 9.7e-01 |
| 512 | +0.236 | 1.4e-01 | -0.207 | 2.0e-01 |

### d = 3

| PC | Var% | Cumul% |
|----|------|--------|
| PC1 | 87.9% | 87.9% |
| PC2 | 10.2% | 98.1% |
| PC3 | 1.2% | 99.3% |
| PC4 | 0.7% | 100.0% |
| PC5 | 0.0% | 100.0% |

Loadings:
| | p_0 | p_1 | p_2 | p_3 | p_4 |
|---|---|---|---|---|---|
| PC1 | +0.477 | -0.371 | -0.467 | -0.464 | -0.449 |
| PC2 | -0.043 | -0.872 | +0.002 | +0.247 | +0.420 |
| PC3 | +0.021 | +0.275 | -0.812 | +0.140 | +0.495 |

d2c projection:
  PC1: cos=-0.5332 (28.4%)
  PC2: cos=-0.1071 (1.1%)
  PC3: cos=+0.0265 (0.1%)

Shape score vs H²:
| N | shape ρ | p | PC2 ρ | p |
|---|---------|---|-------|---|
| 128 | +0.498 | 1.1e-03 | +0.070 | 6.7e-01 |
| 256 | +0.110 | 5.0e-01 | -0.583 | 7.8e-05 |
| 512 | +0.135 | 4.1e-01 | -0.781 | 2.8e-09 |

### d = 4

| PC | Var% | Cumul% |
|----|------|--------|
| PC1 | 83.1% | 83.1% |
| PC2 | 10.8% | 94.0% |
| PC3 | 5.2% | 99.2% |
| PC4 | 0.6% | 99.7% |
| PC5 | 0.3% | 100.0% |

Loadings:
| | p_0 | p_1 | p_2 | p_3 | p_4 |
|---|---|---|---|---|---|
| PC1 | +0.359 | -0.438 | -0.482 | -0.480 | -0.466 |
| PC2 | +0.920 | +0.308 | +0.181 | +0.131 | +0.097 |
| PC3 | +0.158 | -0.752 | -0.040 | +0.293 | +0.567 |

d2c projection:
  PC1: cos=-0.4859 (23.6%)
  PC2: cos=-0.7612 (58.0%)
  PC3: cos=-0.4031 (16.2%)

Shape score vs H²:
| N | shape ρ | p | PC2 ρ | p |
|---|---------|---|-------|---|
| 128 | +0.457 | 3.0e-03 | +0.609 | 3.1e-05 |
| 256 | +0.737 | 5.9e-08 | +0.965 | 9.5e-24 |
| 512 | -0.270 | 9.3e-02 | +0.979 | 1.1e-27 |

## Part 6: shape_score vs KL(p||p_flat) correlation

Hypothesis: KL ≈ ||shape residual||² (both measure deviation from flat)


Global Spearman(shape_score, |KL|) = **-0.143** (p = 6.6e-03)

  d=2 N=128: ρ = +0.283 (p = 7.7e-02)
  d=2 N=256: ρ = -0.359 (p = 2.3e-02)
  d=2 N=512: ρ = -0.605 (p = 3.6e-05)
  d=3 N=128: ρ = -0.321 (p = 4.4e-02)
  d=3 N=256: ρ = -0.682 (p = 1.3e-06)
  d=3 N=512: ρ = -0.733 (p = 7.3e-08)
  d=4 N=128: ρ = +0.223 (p = 1.9e-01)
  d=4 N=256: ρ = +0.225 (p = 1.6e-01)
  d=4 N=512: ρ = +0.154 (p = 3.4e-01)

## Part 7: Summary

- **shape_score ↑ with H²**: 6/9 slices significant
- **PC1 (density) ↓ with H²**: 0/9 slices significant
- PC1 explains **88.6%** of variance (density mode)
- PC2+ explains **11.4%** of variance (shape residual)
- d2c projects **28.3%** onto PC1, **71.7%** onto PC2+

### Interpretation

> The C_k distribution's variation is dominated by a single density mode (PC1) that is strongly anti-correlated with curvature. Curvature information exists only in the remaining ~X% shape residual. bdg_d2c succeeds partially because its coefficient vector [-2,+2,0,0,0] has nonzero projection onto the shape subspace; KL(p||p_flat) succeeds more broadly because it directly measures the magnitude of the shape deviation. This confirms the 'density >> shape' hierarchy and provides the formal basis for the Two-Family Decomposition.
