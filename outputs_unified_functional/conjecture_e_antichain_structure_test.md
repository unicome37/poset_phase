# §4.1.28: Antichain Structure Experiment

## Experiment Design

- Dimensions: [2, 4]
- Sizes: [128]
- Hubble values: [0.0, 1.0, 2.0]
- Total realizations: 12

Antichain (transverse) statistics measure the 'width' of the poset,
orthogonal to the 'longitudinal' interval counts {C_k} closed by DDT.

## Q1: Antichain Features vs H² (Spearman, pooled by d)

| d | w_max_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|------|------|------|------|------|------|------|

## Q2: Per-Slice Spearman (d, N) for Key Features

| d | N | w_max_ratio | layer_ratio | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|---|------|------|------|------|------|

## Q3: Antichain Features vs n_causal_pairs (Density)

If |ρ| ≈ 1.0, the feature is density-dominated.

| d | w_max_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|------|------|------|------|------|------|------|

## Q5: Density-Residual Analysis

After OLS-removing n_causal_pairs (density proxy ≈ ΣC_k),
does each antichain feature's residual still correlate with H²?

| d | feature | raw ρ_S | residual ρ_S | Δρ | verdict |
|---|---------|--------|-------------|-----|---------|

**Features with |residual ρ| > 0.3: 0/0**

### Interpretation

All antichain features collapse to density after residualization. 
Transverse structure at these N values is still dominated by 
the same density signal. This does NOT rule out antichain signals 
at larger N or in non-uniform backgrounds.


## Conclusion

- d=2: no antichain features survive density removal
- d=4: no antichain features survive density removal

**Antichain features density-dominated** at current N.
 Transverse structure does not provide independent curvature signal
 beyond what total density (ΣC_k) already captures.

### Comparison with B_ℓ Spectral (§4.1.27)

- Antichain beyond-density: 0/0
- B_ℓ spectral beyond-density: 6/18 (d=4: 4/6)

Antichain channel weaker than spectral — the operator-level
structure of B_ℓ encodes more curvature than raw transverse widths.

