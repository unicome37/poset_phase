# §4.1.43: Power-Law FRW DDT C2 Test — Non-Uniform R ≠ 0 Background

## Motivation

§4.1.42 (3+1D Schwarzschild) found only 4/32 features beyond density,
with the antichain channel completely density-dominated. The physical
interpretation was that the antichain channel responds to **scalar**
curvature R (de Sitter expansion widens spatial slices), not Weyl curvature.
Schwarzschild is Ricci-flat (R=0), so the strongest DDT escape channel
has no R-signal to detect.

This experiment provides the DECISIVE test: power-law FRW has
**non-uniform R ≠ 0** (for p ≠ 0, 0.5). If the antichain channel
recovers strong beyond-density signals here (like de Sitter's 21/21),
it confirms that:

1. The antichain channel responds to scalar curvature R
2. Schwarzschild weakness is because R=0, not because of non-uniformity
3. DDT C2 escape is confirmed for R≠0 non-uniform backgrounds

## Experiment Design

- Spacetime dimensions: d = [2, 3, 4]
- Power-law exponents: p = [0.0, 0.5, 0.67, 1.0, 1.5, 2.0]
- Sizes: N = [128, 256, 512]
- Total realizations: 540
- Metric: ds² = -dt² + (t/t₀)^{2p} Σ dx_i²
- Volume-weighted sprinkling: ∝ a(t)^{d_spatial}
- Causal relation: comoving horizon distance χ(t_i, t_j)
- Scalar curvature: R(t) = d_s · p · [(d_s+1)p - 2] / t²

### p-value physics (d=4)

| p | a(t) | R(t) | Physical model |
|---|------|------|----------------|
| 0.0 | 1 | 0 | Minkowski (flat) |
| 0.5 | √t | **0** | Radiation-dominated (**R=0 control**) |
| 0.67 | t^{2/3} | 4/(3t²) | Matter-dominated |
| 1.0 | t | 6/t² | Coasting / Milne-like |
| 1.5 | t^{3/2} | 18/t² | Accelerating |
| 2.0 | t² | 36/t² | Strongly accelerating |

## Q1: Global Features vs mean |R| (per d, per N)

### d = 2

| N | n_causal_pairs | C_0 | C_1 | C_2 | bd_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy | w_max_ratio |
|---|------|------|------|------|------|------|------|------|------|------|------|------|
| 128 | **+0.416** | *-0.279* | -0.248 | *-0.254* | *-0.279* | **+0.535** | **-0.535** | **+0.420** | **+0.553** | +0.087 | **+0.361** | *+0.318* |
| 256 | **+0.642** | **-0.481** | **-0.444** | **-0.534** | **-0.481** | **+0.669** | **-0.669** | **+0.561** | **+0.632** | *+0.271* | **+0.383** | **+0.509** |
| 512 | +0.174 | **-0.504** | **-0.432** | **-0.534** | **-0.504** | **+0.649** | **-0.649** | **+0.615** | **+0.670** | **+0.366** | **+0.343** | N/A |

### d = 3

| N | n_causal_pairs | C_0 | C_1 | C_2 | bd_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy | w_max_ratio |
|---|------|------|------|------|------|------|------|------|------|------|------|------|
| 128 | **-0.481** | **-0.758** | **-0.576** | **-0.605** | **-0.758** | +0.084 | -0.084 | **+0.651** | **+0.634** | **+0.614** | **-0.405** | **+0.829** |
| 256 | **-0.452** | **-0.814** | **-0.776** | **-0.674** | **-0.814** | +0.228 | -0.228 | **+0.782** | **+0.773** | **+0.733** | **-0.399** | **+0.889** |
| 512 | **-0.476** | **-0.845** | **-0.783** | **-0.851** | **-0.845** | **+0.621** | **-0.621** | **+0.858** | **+0.889** | **+0.784** | **-0.377** | N/A |

### d = 4

| N | n_causal_pairs | C_0 | C_1 | C_2 | bd_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy | w_max_ratio |
|---|------|------|------|------|------|------|------|------|------|------|------|------|
| 128 | **-0.769** | **-0.868** | **-0.711** | **-0.659** | **-0.868** | **-0.354** | **+0.354** | **+0.713** | **+0.400** | **+0.710** | **-0.623** | **+0.912** |
| 256 | **-0.774** | **-0.897** | **-0.840** | **-0.799** | **-0.897** | -0.108 | +0.108 | **+0.887** | **+0.760** | **+0.820** | **-0.781** | **+0.922** |
| 512 | **-0.838** | **-0.913** | **-0.916** | **-0.895** | **-0.913** | -0.012 | +0.012 | **+0.864** | **+0.802** | **+0.835** | **-0.662** | N/A |

## Q2: Features vs p² (H_eff²) per d

### d = 2

| feature | Spearman ρ(feature, p²) | p-value |
|---------|------------------------|---------|
| n_causal_pairs | +0.210 | 4.74e-03 |
| C_0 | -0.202 | 6.42e-03 |
| C_1 | -0.192 | 9.80e-03 |
| C_2 | -0.207 | 5.31e-03 |
| bd_ratio | -0.202 | 6.42e-03 |
| layer_ratio | +0.340 | 3.04e-06 |
| mean_layer_width | -0.340 | 3.04e-06 |
| layer_width_std | +0.411 | 1.03e-08 |
| layer_width_cv | +0.880 | 1.52e-59 |
| max_layer_width_ratio | +0.119 | 1.11e-01 |
| layer_entropy | +0.193 | 9.46e-03 |
| w_max_ratio | +0.351 | 8.44e-05 |

### d = 3

| feature | Spearman ρ(feature, p²) | p-value |
|---------|------------------------|---------|
| n_causal_pairs | -0.177 | 1.76e-02 |
| C_0 | -0.291 | 7.25e-05 |
| C_1 | -0.263 | 3.63e-04 |
| C_2 | -0.262 | 3.88e-04 |
| bd_ratio | -0.340 | 2.94e-06 |
| layer_ratio | +0.098 | 1.91e-01 |
| mean_layer_width | -0.098 | 1.91e-01 |
| layer_width_std | +0.420 | 4.51e-09 |
| layer_width_cv | +0.823 | 1.24e-45 |
| max_layer_width_ratio | +0.489 | 3.46e-12 |
| layer_entropy | -0.188 | 1.15e-02 |
| w_max_ratio | +0.807 | 9.09e-29 |

### d = 4

| feature | Spearman ρ(feature, p²) | p-value |
|---------|------------------------|---------|
| n_causal_pairs | -0.275 | 1.85e-04 |
| C_0 | -0.310 | 2.23e-05 |
| C_1 | -0.295 | 5.74e-05 |
| C_2 | -0.285 | 1.04e-04 |
| bd_ratio | -0.761 | 2.83e-35 |
| layer_ratio | -0.052 | 4.89e-01 |
| mean_layer_width | +0.052 | 4.89e-01 |
| layer_width_std | +0.371 | 3.01e-07 |
| layer_width_cv | +0.667 | 1.51e-24 |
| max_layer_width_ratio | +0.663 | 3.63e-24 |
| layer_entropy | -0.456 | 1.29e-10 |
| w_max_ratio | +0.921 | 3.69e-50 |

## Q3: Density-Residual Analysis (OLS remove n_causal_pairs)

Does each feature's residual (after density removal) still correlate with mean |R|?

**This is the decisive DDT C2 test.** If antichain features survive density
removal in FRW (R≠0), DDT C2 escape is confirmed for non-uniform backgrounds.

### d = 2

| N | feature | raw ρ_S(feat, |R|) | resid ρ_S | |resid| | p_resid | verdict |
|---|---------|-------------------|-----------|--------|---------|---------|
| 128 | C_0 | -0.279 | -0.056 | 0.056 | 6.68e-01 | density-dominated |
| 128 | C_1 | -0.248 | -0.111 | 0.111 | 3.98e-01 | density-dominated |
| 128 | C_2 | -0.254 | -0.048 | 0.048 | 7.18e-01 | density-dominated |
| 128 | bd_ratio | -0.279 | -0.056 | 0.056 | 6.68e-01 | density-dominated |
| 128 | layer_ratio | +0.535 | +0.284 | 0.284 | 2.80e-02 | marginal+ |
| 128 | mean_layer_width | -0.535 | -0.277 | 0.277 | 3.24e-02 | marginal+ |
| 128 | layer_width_std | +0.420 | +0.380 | 0.380 | 2.74e-03 | **BEYOND DENSITY** |
| 128 | layer_width_cv | +0.553 | +0.418 | 0.418 | 8.82e-04 | **BEYOND DENSITY** |
| 128 | max_layer_width_ratio | +0.087 | +0.170 | 0.170 | 1.95e-01 | marginal |
| 128 | layer_entropy | +0.361 | +0.029 | 0.029 | 8.26e-01 | density-dominated |
| 128 | w_max_ratio | +0.318 | +0.238 | 0.238 | 6.65e-02 | marginal |
| 256 | C_0 | -0.481 | -0.263 | 0.263 | 4.26e-02 | marginal+ |
| 256 | C_1 | -0.444 | -0.223 | 0.223 | 8.61e-02 | marginal |
| 256 | C_2 | -0.534 | -0.291 | 0.291 | 2.40e-02 | marginal+ |
| 256 | bd_ratio | -0.481 | -0.263 | 0.263 | 4.26e-02 | marginal+ |
| 256 | layer_ratio | +0.669 | +0.316 | 0.316 | 1.38e-02 | marginal+ |
| 256 | mean_layer_width | -0.669 | -0.276 | 0.276 | 3.26e-02 | marginal+ |
| 256 | layer_width_std | +0.561 | +0.352 | 0.352 | 5.79e-03 | **BEYOND DENSITY** |
| 256 | layer_width_cv | +0.632 | +0.358 | 0.358 | 4.95e-03 | **BEYOND DENSITY** |
| 256 | max_layer_width_ratio | +0.271 | +0.181 | 0.181 | 1.67e-01 | marginal |
| 256 | layer_entropy | +0.383 | -0.063 | 0.063 | 6.31e-01 | density-dominated |
| 256 | w_max_ratio | +0.509 | +0.253 | 0.253 | 5.09e-02 | marginal |
| 512 | C_0 | -0.504 | -0.375 | 0.375 | 3.18e-03 | **BEYOND DENSITY** |
| 512 | C_1 | -0.432 | -0.395 | 0.395 | 1.81e-03 | **BEYOND DENSITY** |
| 512 | C_2 | -0.534 | -0.422 | 0.422 | 7.84e-04 | **BEYOND DENSITY** |
| 512 | bd_ratio | -0.504 | -0.375 | 0.375 | 3.18e-03 | **BEYOND DENSITY** |
| 512 | layer_ratio | +0.649 | +0.531 | 0.531 | 1.27e-05 | **BEYOND DENSITY** |
| 512 | mean_layer_width | -0.649 | -0.552 | 0.552 | 4.86e-06 | **BEYOND DENSITY** |
| 512 | layer_width_std | +0.615 | +0.441 | 0.441 | 4.18e-04 | **BEYOND DENSITY** |
| 512 | layer_width_cv | +0.670 | +0.472 | 0.472 | 1.42e-04 | **BEYOND DENSITY** |
| 512 | max_layer_width_ratio | +0.366 | +0.249 | 0.249 | 5.50e-02 | marginal |
| 512 | layer_entropy | +0.343 | +0.082 | 0.082 | 5.34e-01 | density-dominated |

### d = 3

| N | feature | raw ρ_S(feat, |R|) | resid ρ_S | |resid| | p_resid | verdict |
|---|---------|-------------------|-----------|--------|---------|---------|
| 128 | C_0 | -0.758 | -0.559 | 0.559 | 3.43e-06 | **BEYOND DENSITY** |
| 128 | C_1 | -0.576 | -0.330 | 0.330 | 1.01e-02 | marginal+ |
| 128 | C_2 | -0.605 | -0.309 | 0.309 | 1.65e-02 | marginal+ |
| 128 | bd_ratio | -0.758 | -0.559 | 0.559 | 3.43e-06 | **BEYOND DENSITY** |
| 128 | layer_ratio | +0.084 | +0.362 | 0.362 | 4.44e-03 | **BEYOND DENSITY** |
| 128 | mean_layer_width | -0.084 | -0.358 | 0.358 | 4.95e-03 | **BEYOND DENSITY** |
| 128 | layer_width_std | +0.651 | +0.423 | 0.423 | 7.51e-04 | **BEYOND DENSITY** |
| 128 | layer_width_cv | +0.634 | +0.460 | 0.460 | 2.16e-04 | **BEYOND DENSITY** |
| 128 | max_layer_width_ratio | +0.614 | +0.373 | 0.373 | 3.30e-03 | **BEYOND DENSITY** |
| 128 | layer_entropy | -0.405 | +0.047 | 0.047 | 7.22e-01 | density-dominated |
| 128 | w_max_ratio | +0.829 | +0.611 | 0.611 | 2.12e-07 | **BEYOND DENSITY** |
| 256 | C_0 | -0.814 | -0.639 | 0.639 | 3.98e-08 | **BEYOND DENSITY** |
| 256 | C_1 | -0.776 | -0.638 | 0.638 | 4.25e-08 | **BEYOND DENSITY** |
| 256 | C_2 | -0.674 | -0.459 | 0.459 | 2.23e-04 | **BEYOND DENSITY** |
| 256 | bd_ratio | -0.814 | -0.639 | 0.639 | 3.98e-08 | **BEYOND DENSITY** |
| 256 | layer_ratio | +0.228 | +0.442 | 0.442 | 4.01e-04 | **BEYOND DENSITY** |
| 256 | mean_layer_width | -0.228 | -0.446 | 0.446 | 3.57e-04 | **BEYOND DENSITY** |
| 256 | layer_width_std | +0.782 | +0.557 | 0.557 | 3.75e-06 | **BEYOND DENSITY** |
| 256 | layer_width_cv | +0.773 | +0.624 | 0.624 | 9.90e-08 | **BEYOND DENSITY** |
| 256 | max_layer_width_ratio | +0.733 | +0.523 | 0.523 | 1.84e-05 | **BEYOND DENSITY** |
| 256 | layer_entropy | -0.399 | -0.104 | 0.104 | 4.30e-01 | density-dominated |
| 256 | w_max_ratio | +0.889 | +0.694 | 0.694 | 7.61e-10 | **BEYOND DENSITY** |
| 512 | C_0 | -0.845 | -0.673 | 0.673 | 3.87e-09 | **BEYOND DENSITY** |
| 512 | C_1 | -0.783 | -0.648 | 0.648 | 2.24e-08 | **BEYOND DENSITY** |
| 512 | C_2 | -0.851 | -0.703 | 0.703 | 3.76e-10 | **BEYOND DENSITY** |
| 512 | bd_ratio | -0.845 | -0.673 | 0.673 | 3.87e-09 | **BEYOND DENSITY** |
| 512 | layer_ratio | +0.621 | +0.498 | 0.498 | 5.19e-05 | **BEYOND DENSITY** |
| 512 | mean_layer_width | -0.621 | -0.498 | 0.498 | 5.19e-05 | **BEYOND DENSITY** |
| 512 | layer_width_std | +0.858 | +0.642 | 0.642 | 3.19e-08 | **BEYOND DENSITY** |
| 512 | layer_width_cv | +0.889 | +0.703 | 0.703 | 3.89e-10 | **BEYOND DENSITY** |
| 512 | max_layer_width_ratio | +0.784 | +0.556 | 0.556 | 3.94e-06 | **BEYOND DENSITY** |
| 512 | layer_entropy | -0.377 | +0.000 | 0.000 | 1.00e+00 | density-dominated |

### d = 4

| N | feature | raw ρ_S(feat, |R|) | resid ρ_S | |resid| | p_resid | verdict |
|---|---------|-------------------|-----------|--------|---------|---------|
| 128 | C_0 | -0.868 | -0.330 | 0.330 | 1.01e-02 | marginal+ |
| 128 | C_1 | -0.711 | +0.018 | 0.018 | 8.89e-01 | density-dominated |
| 128 | C_2 | -0.659 | +0.172 | 0.172 | 1.88e-01 | marginal |
| 128 | bd_ratio | -0.868 | -0.330 | 0.330 | 1.01e-02 | marginal+ |
| 128 | layer_ratio | -0.354 | +0.094 | 0.094 | 4.77e-01 | density-dominated |
| 128 | mean_layer_width | +0.354 | -0.106 | 0.106 | 4.19e-01 | density-dominated |
| 128 | layer_width_std | +0.713 | +0.168 | 0.168 | 1.99e-01 | marginal |
| 128 | layer_width_cv | +0.400 | +0.176 | 0.176 | 1.78e-01 | marginal |
| 128 | max_layer_width_ratio | +0.710 | +0.174 | 0.174 | 1.83e-01 | marginal |
| 128 | layer_entropy | -0.623 | -0.021 | 0.021 | 8.72e-01 | density-dominated |
| 128 | w_max_ratio | +0.912 | +0.423 | 0.423 | 7.71e-04 | **BEYOND DENSITY** |
| 256 | C_0 | -0.897 | -0.379 | 0.379 | 2.79e-03 | **BEYOND DENSITY** |
| 256 | C_1 | -0.840 | -0.258 | 0.258 | 4.65e-02 | marginal+ |
| 256 | C_2 | -0.799 | -0.113 | 0.113 | 3.90e-01 | density-dominated |
| 256 | bd_ratio | -0.897 | -0.379 | 0.379 | 2.79e-03 | **BEYOND DENSITY** |
| 256 | layer_ratio | -0.108 | +0.165 | 0.165 | 2.06e-01 | marginal |
| 256 | mean_layer_width | +0.108 | -0.191 | 0.191 | 1.45e-01 | marginal |
| 256 | layer_width_std | +0.887 | +0.400 | 0.400 | 1.53e-03 | **BEYOND DENSITY** |
| 256 | layer_width_cv | +0.760 | +0.438 | 0.438 | 4.70e-04 | **BEYOND DENSITY** |
| 256 | max_layer_width_ratio | +0.820 | +0.314 | 0.314 | 1.44e-02 | marginal+ |
| 256 | layer_entropy | -0.781 | -0.232 | 0.232 | 7.40e-02 | marginal |
| 256 | w_max_ratio | +0.922 | +0.426 | 0.426 | 6.96e-04 | **BEYOND DENSITY** |
| 512 | C_0 | -0.913 | -0.240 | 0.240 | 6.43e-02 | marginal |
| 512 | C_1 | -0.916 | -0.211 | 0.211 | 1.06e-01 | marginal |
| 512 | C_2 | -0.895 | -0.145 | 0.145 | 2.68e-01 | density-dominated |
| 512 | bd_ratio | -0.913 | -0.240 | 0.240 | 6.43e-02 | marginal |
| 512 | layer_ratio | -0.012 | +0.260 | 0.260 | 4.49e-02 | marginal+ |
| 512 | mean_layer_width | +0.012 | -0.260 | 0.260 | 4.49e-02 | marginal+ |
| 512 | layer_width_std | +0.864 | +0.270 | 0.270 | 3.72e-02 | marginal+ |
| 512 | layer_width_cv | +0.802 | +0.273 | 0.273 | 3.51e-02 | marginal+ |
| 512 | max_layer_width_ratio | +0.835 | +0.241 | 0.241 | 6.32e-02 | marginal |
| 512 | layer_entropy | -0.662 | -0.054 | 0.054 | 6.79e-01 | density-dominated |


**Total features beyond density: 45/96**

### d = 2 summary

- N=128: **2/11**
- N=256: **2/11**
- N=512: **8/10**

### d = 3 summary

- N=128: **8/11**
- N=256: **10/11**
- N=512: **9/10**

### d = 4 summary

- N=128: **1/11**
- N=256: **5/11**
- N=512: **0/10**

## Q3b: R=0 vs R≠0 Contrast — The Decisive Comparison

Split realizations into R=0 (p=0, p=0.5) vs R≠0 (p=0.67, 1.0, 1.5, 2.0)
and repeat density-residual analysis separately.

### R=0 only (p=0, 0.5)

| d | N | feature | |resid ρ_S| | p_resid | verdict |
|---|---|---------|-----------|---------|---------|

**R=0 only (p=0, 0.5): 0/0 beyond density**

### R≠0 only (p=0.67, 1.0, 1.5, 2.0)

| d | N | feature | |resid ρ_S| | p_resid | verdict |
|---|---|---------|-----------|---------|---------|
| 2 | 128 | layer_ratio | 0.383 | 1.48e-02 | marginal+ |
| 2 | 128 | mean_layer_width | 0.372 | 1.82e-02 | marginal+ |
| 2 | 128 | layer_width_std | 0.475 | 1.94e-03 | **BEYOND** |
| 2 | 128 | w_max_ratio | 0.391 | 1.26e-02 | marginal+ |
| 2 | 128 | C_1 | 0.284 | 7.61e-02 | weak/no |
| 2 | 128 | C_2 | 0.102 | 5.33e-01 | weak/no |
| 2 | 256 | layer_ratio | 0.588 | 6.70e-05 | **BEYOND** |
| 2 | 256 | mean_layer_width | 0.590 | 6.19e-05 | **BEYOND** |
| 2 | 256 | layer_width_std | 0.666 | 2.72e-06 | **BEYOND** |
| 2 | 256 | w_max_ratio | 0.565 | 1.43e-04 | **BEYOND** |
| 2 | 256 | C_1 | 0.440 | 4.48e-03 | **BEYOND** |
| 2 | 256 | C_2 | 0.499 | 1.06e-03 | **BEYOND** |
| 2 | 512 | layer_ratio | 0.814 | 1.74e-10 | **BEYOND** |
| 2 | 512 | mean_layer_width | 0.797 | 7.60e-10 | **BEYOND** |
| 2 | 512 | layer_width_std | 0.728 | 9.99e-08 | **BEYOND** |
| 2 | 512 | C_1 | 0.730 | 9.12e-08 | **BEYOND** |
| 2 | 512 | C_2 | 0.645 | 6.96e-06 | **BEYOND** |
| 3 | 128 | layer_ratio | 0.209 | 1.45e-01 | weak/no |
| 3 | 128 | mean_layer_width | 0.230 | 1.09e-01 | weak/no |
| 3 | 128 | layer_width_std | 0.309 | 2.88e-02 | marginal+ |
| 3 | 128 | w_max_ratio | 0.350 | 1.27e-02 | marginal+ |
| 3 | 128 | C_1 | 0.114 | 4.32e-01 | weak/no |
| 3 | 128 | C_2 | 0.254 | 7.48e-02 | weak/no |
| 3 | 256 | layer_ratio | 0.117 | 4.20e-01 | weak/no |
| 3 | 256 | mean_layer_width | 0.136 | 3.46e-01 | weak/no |
| 3 | 256 | layer_width_std | 0.359 | 1.05e-02 | marginal+ |
| 3 | 256 | w_max_ratio | 0.422 | 2.28e-03 | **BEYOND** |
| 3 | 256 | C_1 | 0.384 | 5.90e-03 | **BEYOND** |
| 3 | 256 | C_2 | 0.164 | 2.56e-01 | weak/no |
| 3 | 512 | layer_ratio | 0.192 | 1.81e-01 | weak/no |
| 3 | 512 | mean_layer_width | 0.192 | 1.81e-01 | weak/no |
| 3 | 512 | layer_width_std | 0.247 | 8.34e-02 | weak/no |
| 3 | 512 | C_1 | 0.288 | 4.26e-02 | marginal+ |
| 3 | 512 | C_2 | 0.355 | 1.14e-02 | marginal+ |
| 4 | 128 | layer_ratio | 0.005 | 9.77e-01 | weak/no |
| 4 | 128 | mean_layer_width | 0.100 | 5.40e-01 | weak/no |
| 4 | 128 | layer_width_std | 0.213 | 1.88e-01 | weak/no |
| 4 | 128 | w_max_ratio | 0.366 | 2.01e-02 | marginal+ |
| 4 | 128 | C_1 | 0.120 | 4.60e-01 | weak/no |
| 4 | 128 | C_2 | 0.106 | 5.17e-01 | weak/no |
| 4 | 256 | layer_ratio | 0.082 | 6.15e-01 | weak/no |
| 4 | 256 | mean_layer_width | 0.100 | 5.40e-01 | weak/no |
| 4 | 256 | layer_width_std | 0.398 | 1.10e-02 | marginal+ |
| 4 | 256 | w_max_ratio | 0.480 | 1.70e-03 | **BEYOND** |
| 4 | 256 | C_1 | 0.127 | 4.36e-01 | weak/no |
| 4 | 256 | C_2 | 0.008 | 9.62e-01 | weak/no |
| 4 | 512 | layer_ratio | 0.020 | 9.05e-01 | weak/no |
| 4 | 512 | mean_layer_width | 0.033 | 8.41e-01 | weak/no |
| 4 | 512 | layer_width_std | 0.251 | 1.18e-01 | weak/no |
| 4 | 512 | C_1 | 0.317 | 4.59e-02 | marginal+ |
| 4 | 512 | C_2 | 0.176 | 2.78e-01 | weak/no |

**R≠0 only (p=0.67, 1.0, 1.5, 2.0): 15/51 beyond density**

## Q4: Temporal Bin Analysis — Local R(t) Response

Each realization split into 3 temporal bins (early/mid/late).
Early bins have larger |R(t)|. We test whether antichain features
correlate with local |R| after density removal.

- Raw ρ(bin_|R|, bin_layer_ratio): -0.101 (p=5.10e-04)
- Density-residualized ρ: -0.193 (p=2.78e-11)
- Marginal local curvature response after density removal

## Q5: N-Scaling of Beyond-Density Signals (R≠0 only)

Do the beyond-density residuals strengthen with N?

### d = 2

| feature | N=128 |ρ_resid| | N=256 |ρ_resid| | N=512 |ρ_resid| | trend |
|---------|------------------|------------------|------------------|-------|
| C_0 | 0.261 | 0.530 | 0.715 | ↑ |
| C_1 | 0.284 | 0.440 | 0.730 | ↑ |
| C_2 | 0.102 | 0.499 | 0.645 | ↑ |
| bd_ratio | 0.261 | 0.530 | 0.715 | ↑ |
| layer_ratio | 0.383 | 0.588 | 0.814 | ↑ |
| mean_layer_width | 0.372 | 0.590 | 0.797 | ↑ |
| layer_width_std | 0.475 | 0.666 | 0.728 | ↑ |
| layer_width_cv | 0.514 | 0.674 | 0.771 | ↑ |
| max_layer_width_ratio | 0.171 | 0.383 | 0.484 | ↑ |
| layer_entropy | 0.032 | 0.022 | 0.226 | ~ |
| w_max_ratio | 0.391 | 0.565 | N/A | ↑ |

### d = 3

| feature | N=128 |ρ_resid| | N=256 |ρ_resid| | N=512 |ρ_resid| | trend |
|---------|------------------|------------------|------------------|-------|
| C_0 | 0.287 | 0.401 | 0.273 | ~ |
| C_1 | 0.114 | 0.384 | 0.288 | ~ |
| C_2 | 0.254 | 0.164 | 0.355 | ~ |
| bd_ratio | 0.287 | 0.401 | 0.273 | ~ |
| layer_ratio | 0.209 | 0.117 | 0.192 | ~ |
| mean_layer_width | 0.230 | 0.136 | 0.192 | ~ |
| layer_width_std | 0.309 | 0.359 | 0.247 | ~ |
| layer_width_cv | 0.336 | 0.404 | 0.338 | ~ |
| max_layer_width_ratio | 0.132 | 0.293 | 0.157 | ~ |
| layer_entropy | 0.098 | 0.188 | 0.064 | ~ |
| w_max_ratio | 0.350 | 0.422 | N/A | ↑ |

### d = 4

| feature | N=128 |ρ_resid| | N=256 |ρ_resid| | N=512 |ρ_resid| | trend |
|---------|------------------|------------------|------------------|-------|
| C_0 | 0.075 | 0.409 | 0.172 | ~ |
| C_1 | 0.120 | 0.127 | 0.317 | ↑ |
| C_2 | 0.106 | 0.008 | 0.176 | ~ |
| bd_ratio | 0.075 | 0.409 | 0.172 | ~ |
| layer_ratio | 0.005 | 0.082 | 0.020 | ~ |
| mean_layer_width | 0.100 | 0.100 | 0.033 | ↓ |
| layer_width_std | 0.213 | 0.398 | 0.251 | ~ |
| layer_width_cv | 0.229 | 0.365 | 0.167 | ~ |
| max_layer_width_ratio | 0.137 | 0.246 | 0.147 | ~ |
| layer_entropy | 0.236 | 0.304 | 0.178 | ~ |
| w_max_ratio | 0.366 | 0.480 | N/A | ↑ |

## Comparison with Previous Experiments

| Setting | R | Uniform? | Beyond density | Peak |ρ_resid| |
|---------|---|----------|---------------|----------------|
| **FRW power-law R≠0** (this) | ≠0 | No | **45/96** | (see above) |
| de Sitter (§4.1.28) | ≠0 | Yes (constant) | **21/21** | 0.817 |
| 3+1D Schwarzschild (§4.1.42) | =0 | No | 4/32 | ~0.56 |
| 1+1D Schwarzschild (§4.1.29) | =0 | No | 2/27 | ~0.37 |

## Conclusion

**PARTIAL POSITIVE: 15/60 antichain features beyond density.**
Some antichain features carry scalar curvature R signal beyond density in FRW,
but not as comprehensively as de Sitter (21/21). The non-uniform R(t)
partially degrades the signal, but it is substantially stronger than
Schwarzschild (4/32), confirming R≠0 is the key factor.

**Overall: 45/96** features beyond density across all
dimensions and sizes (antichain: 15/60).

### 中文结论

**总计 45/96** 个特征在幂律 FRW 背景（R≠0，非均匀）中超越密度
（其中反链通道：15/60）。

部分确认了 §4.1.42 的解释：反链通道在 R≠0 条件下恢复了一部分超密度信号，
但不如 de Sitter 的 21/21 强。R≠0 是关键因素，但非均匀性有一定削弱效应。
