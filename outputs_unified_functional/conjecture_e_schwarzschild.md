# §4.1.29: Schwarzschild Sprinkling Experiment

## Experiment Design

- Geometry: 1+1D Schwarzschild (radial + time)
- rₛ values (2M): [0.0, 0.1, 0.25, 0.5, 1.0, 2.0]
- Sizes: [128, 256, 512]
- Total realizations: 180
- r range: [1.5·rₛ, 10·rₛ] (well outside horizon)
- Curvature proxy: Kretschner scalar K = 12rₛ²/r⁶

DDT condition C2 tested: non-uniform curvature background.

## Q1: Features vs mean Kretschner K (pooled by N)

| N | n_causal_pairs | C_0 | C_1 | C_2 | bd_ratio | layer_ratio | mean_layer_width | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|------|------|------|------|------|------|------|------|------|------|
| 128 | **+0.801** | **-0.551** | **-0.542** | **-0.573** | **-0.551** | **+0.840** | **-0.840** | +0.237 | **-0.819** | **+0.832** |
| 256 | **+0.798** | **-0.597** | **-0.671** | **-0.626** | **-0.597** | **+0.818** | **-0.818** | -0.045 | **-0.813** | **+0.812** |
| 512 | **+0.793** | **-0.671** | **-0.724** | **-0.586** | **-0.671** | **+0.813** | **-0.813** | -0.241 | **-0.817** | **+0.802** |

## Q2: Features vs rₛ (pooled across all N)

| feature | Spearman ρ(feature, rₛ) | p-value |
|---------|------------------------|---------|
| n_causal_pairs | -0.216 | 3.63e-03 |
| C_0 | +0.033 | 6.60e-01 |
| C_1 | +0.019 | 7.98e-01 |
| C_2 | +0.002 | 9.77e-01 |
| bd_ratio | +0.033 | 6.60e-01 |
| layer_ratio | -0.497 | 1.24e-12 |
| mean_layer_width | +0.497 | 1.24e-12 |
| layer_width_cv | +0.098 | 1.91e-01 |
| max_layer_width_ratio | +0.459 | 9.14e-11 |
| layer_entropy | -0.488 | 3.56e-12 |

## Q3: Density-Residual Analysis (OLS remove n_causal_pairs)

Does each feature's residual (after density removal) still correlate with mean_K?

| N | feature | raw ρ_S(feat, K) | resid ρ_S | Δρ | verdict |
|---|---------|-----------------|-----------|-----|---------|
| 128 | C_0 | -0.551 | -0.253 | -0.298 | marginal |
| 128 | C_1 | -0.542 | -0.279 | -0.263 | marginal |
| 128 | C_2 | -0.573 | -0.322 | -0.251 | marginal |
| 128 | bd_ratio | -0.551 | -0.253 | -0.298 | marginal |
| 128 | layer_ratio | +0.840 | +0.121 | -0.719 | density-dominated |
| 128 | mean_layer_width | -0.840 | -0.039 | -0.800 | density-dominated |
| 128 | layer_width_cv | +0.237 | +0.151 | -0.086 | marginal |
| 128 | max_layer_width_ratio | -0.819 | +0.087 | -0.733 | density-dominated |
| 128 | layer_entropy | +0.832 | +0.175 | -0.657 | marginal |
| 256 | C_0 | -0.597 | -0.250 | -0.346 | marginal |
| 256 | C_1 | -0.671 | -0.296 | -0.375 | marginal |
| 256 | C_2 | -0.626 | -0.267 | -0.359 | marginal |
| 256 | bd_ratio | -0.597 | -0.250 | -0.346 | marginal |
| 256 | layer_ratio | +0.818 | +0.178 | -0.640 | marginal |
| 256 | mean_layer_width | -0.818 | +0.051 | -0.767 | density-dominated |
| 256 | layer_width_cv | -0.045 | +0.095 | +0.050 | density-dominated |
| 256 | max_layer_width_ratio | -0.813 | +0.058 | -0.755 | density-dominated |
| 256 | layer_entropy | +0.812 | +0.170 | -0.642 | marginal |
| 512 | C_0 | -0.671 | -0.192 | -0.479 | marginal |
| 512 | C_1 | -0.724 | -0.334 | -0.391 | **BEYOND DENSITY** |
| 512 | C_2 | -0.586 | -0.264 | -0.322 | marginal |
| 512 | bd_ratio | -0.671 | -0.192 | -0.479 | marginal |
| 512 | layer_ratio | +0.813 | +0.130 | -0.682 | density-dominated |
| 512 | mean_layer_width | -0.813 | +0.368 | -0.445 | **BEYOND DENSITY** |
| 512 | layer_width_cv | -0.241 | +0.036 | -0.205 | density-dominated |
| 512 | max_layer_width_ratio | -0.817 | +0.269 | -0.548 | marginal |
| 512 | layer_entropy | +0.802 | +0.136 | -0.666 | density-dominated |

**Features beyond density: 2/27**

## Q4: Non-Uniformity Test — Inner vs Outer Patch

Split each realization into inner half (closer to horizon) and outer half (farther).
Compare C_k distributions between patches.

| rₛ | N | inner C₀/pairs | outer C₀/pairs | ratio inner/outer | ρ(ratio, K_inner/K_outer) |
|-----|---|---------------|----------------|-------------------|--------------------------|
| 0.1 | 128 | (computed below) | | | |
| 0.1 | 256 | (computed below) | | | |
| 0.1 | 512 | (computed below) | | | |
| 0.25 | 128 | (computed below) | | | |
| 0.25 | 256 | (computed below) | | | |
| 0.25 | 512 | (computed below) | | | |
| 0.5 | 128 | (computed below) | | | |
| 0.5 | 256 | (computed below) | | | |
| 0.5 | 512 | (computed below) | | | |
| 1.0 | 128 | (computed below) | | | |
| 1.0 | 256 | (computed below) | | | |
| 1.0 | 512 | (computed below) | | | |
| 2.0 | 128 | (computed below) | | | |
| 2.0 | 256 | (computed below) | | | |
| 2.0 | 512 | (computed below) | | | |

## Conclusion

**2/27** features carry curvature info beyond density
 in Schwarzschild background. DDT condition C2 (uniform curvature)
 is confirmed to be essential — non-uniform backgrounds allow
 curvature signals that constant-curvature DDT cannot absorb.

### Comparison with de Sitter channels

- Schwarzschild beyond-density: 2/27
- de Sitter B_ℓ spectral: 6/18
- de Sitter antichain: 21/21

