# §4.1.28: Antichain Structure Experiment

## Experiment Design

- Dimensions: [2, 3, 4]
- Sizes: [128, 256, 512]
- Hubble values: [0.0, 0.25, 0.5, 1.0, 2.0]
- Total realizations: 360

Antichain (transverse) statistics measure the 'width' of the poset,
orthogonal to the 'longitudinal' interval counts {C_k} closed by DDT.

## Q1: Antichain Features vs H² (Spearman, pooled by d)

| d | w_max_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|------|------|------|------|------|------|------|
| 2 | **+0.720** | **-0.478** | **+0.478** | **+0.689** | **+0.668** | **+0.601** | **-0.528** |
| 3 | **+0.888** | **-0.501** | **+0.501** | **+0.724** | **+0.745** | **+0.836** | **-0.765** |
| 4 | **+0.940** | **-0.490** | **+0.490** | **+0.719** | **+0.760** | **+0.886** | **-0.854** |

## Q2: Per-Slice Spearman (d, N) for Key Features

| d | N | w_max_ratio | layer_ratio | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|---|------|------|------|------|------|
| 2 | 128 | **+0.902** | **-0.844** | **+0.539** | **+0.833** | **-0.870** |
| 2 | 256 | **+0.958** | **-0.886** | **+0.741** | **+0.857** | **-0.936** |
| 2 | 512 | **+0.973** | **-0.916** | **+0.701** | **+0.932** | **-0.940** |
| 3 | 128 | **+0.945** | **-0.879** | **+0.668** | **+0.888** | **-0.946** |
| 3 | 256 | **+0.970** | **-0.893** | **+0.802** | **+0.961** | **-0.954** |
| 3 | 512 | **+0.977** | **-0.949** | **+0.786** | **+0.962** | **-0.980** |
| 4 | 128 | **+0.977** | **-0.872** | **+0.708** | **+0.944** | **-0.949** |
| 4 | 256 | **+0.979** | **-0.886** | **+0.802** | **+0.954** | **-0.953** |
| 4 | 512 | **+0.981** | **-0.922** | **+0.787** | **+0.970** | **-0.969** |

## Q3: Antichain Features vs n_causal_pairs (Density)

If |ρ| ≈ 1.0, the feature is density-dominated.

| d | w_max_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|------|------|------|------|------|------|------|
| 2 | -0.812 | -0.628 | +0.628 | +0.307 | -0.296 | -0.866 | +0.926 |
| 3 | -0.816 | -0.338 | +0.338 | +0.033 | -0.439 | -0.844 | +0.922 |
| 4 | -0.888 | -0.106 | +0.106 | -0.190 | -0.549 | -0.905 | +0.945 |

## Q5: Density-Residual Analysis

After OLS-removing n_causal_pairs (density proxy ≈ ΣC_k),
does each antichain feature's residual still correlate with H²?

| d | feature | raw ρ_S | residual ρ_S | Δρ | verdict |
|---|---------|--------|-------------|-----|---------|
| 2 | w_max_ratio | +0.720 | +0.696 | -0.024 | **BEYOND DENSITY** |
| 2 | layer_ratio | -0.478 | -0.794 | +0.315 | **BEYOND DENSITY** |
| 2 | mean_layer_width | +0.478 | +0.817 | +0.338 | **BEYOND DENSITY** |
| 2 | layer_width_std | +0.689 | +0.765 | +0.075 | **BEYOND DENSITY** |
| 2 | layer_width_cv | +0.668 | +0.606 | -0.061 | **BEYOND DENSITY** |
| 2 | max_layer_width_ratio | +0.601 | +0.570 | -0.031 | **BEYOND DENSITY** |
| 2 | layer_entropy | -0.528 | -0.559 | +0.031 | **BEYOND DENSITY** |
| 3 | w_max_ratio | +0.888 | +0.683 | -0.204 | **BEYOND DENSITY** |
| 3 | layer_ratio | -0.501 | -0.755 | +0.254 | **BEYOND DENSITY** |
| 3 | mean_layer_width | +0.501 | +0.665 | +0.164 | **BEYOND DENSITY** |
| 3 | layer_width_std | +0.724 | +0.632 | -0.093 | **BEYOND DENSITY** |
| 3 | layer_width_cv | +0.745 | +0.620 | -0.125 | **BEYOND DENSITY** |
| 3 | max_layer_width_ratio | +0.836 | +0.637 | -0.199 | **BEYOND DENSITY** |
| 3 | layer_entropy | -0.765 | -0.625 | -0.139 | **BEYOND DENSITY** |
| 4 | w_max_ratio | +0.940 | +0.740 | -0.201 | **BEYOND DENSITY** |
| 4 | layer_ratio | -0.490 | -0.730 | +0.240 | **BEYOND DENSITY** |
| 4 | mean_layer_width | +0.490 | +0.585 | +0.095 | **BEYOND DENSITY** |
| 4 | layer_width_std | +0.719 | +0.604 | -0.115 | **BEYOND DENSITY** |
| 4 | layer_width_cv | +0.760 | +0.556 | -0.205 | **BEYOND DENSITY** |
| 4 | max_layer_width_ratio | +0.886 | +0.639 | -0.247 | **BEYOND DENSITY** |
| 4 | layer_entropy | -0.854 | -0.655 | -0.199 | **BEYOND DENSITY** |

**Features with |residual ρ| > 0.3: 21/21**

### Interpretation

Antichain features that survive density removal carry 
information **beyond** the {C_k} family. 
This would validate DDT escape via condition C1 
(transverse structure ≠ interval counts).


## Conclusion

- **d=2**: 7 features BEYOND DENSITY: w_max_ratio(ρ_resid=+0.696), layer_ratio(ρ_resid=-0.794), mean_layer_width(ρ_resid=+0.817), layer_width_std(ρ_resid=+0.765), layer_width_cv(ρ_resid=+0.606), max_layer_width_ratio(ρ_resid=+0.570), layer_entropy(ρ_resid=-0.559)
- **d=3**: 7 features BEYOND DENSITY: w_max_ratio(ρ_resid=+0.683), layer_ratio(ρ_resid=-0.755), mean_layer_width(ρ_resid=+0.665), layer_width_std(ρ_resid=+0.632), layer_width_cv(ρ_resid=+0.620), max_layer_width_ratio(ρ_resid=+0.637), layer_entropy(ρ_resid=-0.625)
- **d=4**: 7 features BEYOND DENSITY: w_max_ratio(ρ_resid=+0.740), layer_ratio(ρ_resid=-0.730), mean_layer_width(ρ_resid=+0.585), layer_width_std(ρ_resid=+0.604), layer_width_cv(ρ_resid=+0.556), max_layer_width_ratio(ρ_resid=+0.639), layer_entropy(ρ_resid=-0.655)

**VERDICT: Antichain structure encodes curvature BEYOND {C_k}.**
 21/21 features survive density residualization.
 DDT condition C1 escaped via transverse (antichain) observables.

### Comparison with B_ℓ Spectral (§4.1.27)

- Antichain beyond-density: 21/21
- B_ℓ spectral beyond-density: 6/18 (d=4: 4/6)

Both transverse and spectral channels provide independent curvature info.
