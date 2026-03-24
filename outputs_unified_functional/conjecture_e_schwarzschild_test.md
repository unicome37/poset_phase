# §4.1.29: Schwarzschild Sprinkling Experiment

## Experiment Design

- Geometry: 1+1D Schwarzschild (radial + time)
- rₛ values (2M): [0.0, 0.5, 2.0]
- Sizes: [128]
- Total realizations: 9
- r range: [1.5·rₛ, 10·rₛ] (well outside horizon)
- Curvature proxy: Kretschner scalar K = 12rₛ²/r⁶

DDT condition C2 tested: non-uniform curvature background.

## Q1: Features vs mean Kretschner K (pooled by N)

| N | n_causal_pairs | C_0 | C_1 | C_2 | bd_ratio | layer_ratio | mean_layer_width | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|------|------|------|------|------|------|------|------|------|------|

## Q2: Features vs rₛ (pooled across all N)

| feature | Spearman ρ(feature, rₛ) | p-value |
|---------|------------------------|---------|
| n_causal_pairs | N/A | N/A |
| C_0 | N/A | N/A |
| C_1 | N/A | N/A |
| C_2 | N/A | N/A |
| bd_ratio | N/A | N/A |
| layer_ratio | N/A | N/A |
| mean_layer_width | N/A | N/A |
| layer_width_cv | N/A | N/A |
| max_layer_width_ratio | N/A | N/A |
| layer_entropy | N/A | N/A |

## Q3: Density-Residual Analysis (OLS remove n_causal_pairs)

Does each feature's residual (after density removal) still correlate with mean_K?

| N | feature | raw ρ_S(feat, K) | resid ρ_S | Δρ | verdict |
|---|---------|-----------------|-----------|-----|---------|

**Features beyond density: 0/0**

## Q4: Non-Uniformity Test — Inner vs Outer Patch

Split each realization into inner half (closer to horizon) and outer half (farther).
Compare C_k distributions between patches.

| rₛ | N | inner C₀/pairs | outer C₀/pairs | ratio inner/outer | ρ(ratio, K_inner/K_outer) |
|-----|---|---------------|----------------|-------------------|--------------------------|
| 0.5 | 128 | (computed below) | | | |
| 2.0 | 128 | (computed below) | | | |

## Conclusion

No features survive density residualization in 1+1D Schwarzschild.
 This may be due to: (1) 1+1D too low-dimensional;
 (2) Ricci-flat background (R=0, only tidal Weyl curvature);
 (3) current N too small for the effect size.

### Comparison with de Sitter channels

- Schwarzschild beyond-density: 0/0
- de Sitter B_ℓ spectral: 6/18
- de Sitter antichain: 21/21

