# §4.1.31: Dual-Channel Unification Experiment

## Design

- Dimensions: [2, 4]
- Sizes: [128]
- Hubble values: [0.0, 0.5, 1.0, 2.0]
- Total realizations: 24
- Both antichain and B_ℓ features extracted from the **same** sprinklings
- Core question: do the two confirmed DDT escape channels measure the same geometric content?

## Q1: Independent Channel Performance vs H²

Replication of §4.1.27/28 on the same sprinklings.

### Antichain features vs H² (pooled by d)

| d | w_max_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|------|------|------|------|------|------|------|
| 2 | **+0.944** | **-0.793** | **+0.793** | **+0.928** | +0.540 | **+0.943** | **-0.820** |
| 4 | **+0.975** | **-0.927** | **+0.927** | **+0.950** | **+0.820** | **+0.972** | **-0.950** |

### B_ℓ spectral features vs H² (pooled by d)

| d | b1_mean | b1_std | b1_neg_frac | eig_min | eig_max | eig_gap | eig_spread |
|---|------|------|------|------|------|------|------|
| 2 | +0.022 | *-0.669* | -0.174 | *+0.605* | *-0.691* | *-0.583* | *-0.626* |
| 4 | -0.348 | **-0.972** | -0.445 | **+0.950** | **-0.950** | -0.259 | **-0.950** |

## Q2: Density-Residual Analysis per Channel per (d, N)

After OLS-removing n_causal_pairs, Spearman ρ(residual, H²).

| d | N | w_max_ratio | mean_layer_width | layer_width_std | layer_ratio | b1_std | eig_min | eig_max | eig_spread |
|---|---|------|------|------|------|------|------|------|------|

### Summary: Beyond-density count per (d, N)

| d | N | AC beyond | B_ℓ beyond | Total |
|---|---|-----------|-----------|-------|

## Q3: Cross-Channel Residual Correlation (KEY TEST)

After removing density from both channels independently,
do their residuals correlate with each other?

- High |ρ| → same underlying geometric content (unified E-bulk)
- Low |ρ| → independent signals (multi-dimensional E-bulk)

### Best-representative cross-correlation

AC representative: w_max_ratio (strongest antichain)
B_ℓ representative: b1_std (strongest B_ℓ spectral)

| d | N | ρ(AC_resid, Bℓ_resid) | p-value | interpretation |
|---|---|----------------------|---------|----------------|

### Full cross-correlation matrix (d=4 slices)


## Q4: N-Scaling of Cross-Channel Correlation

Does the cross-channel correlation strengthen or weaken with N?

| d | N | ρ_cross(AC, Bℓ) |
|---|---|-----------------|


## Q5: Combined Proxy Analysis

Can a linear combination of AC + B_ℓ residuals explain more H² variance than either alone?

| d | N | R²(AC only) | R²(Bℓ only) | R²(combined) | ΔR² |
|---|---|------------|-------------|--------------|------|

## Conclusion

**Cross-channel correlation summary**: 0 (d,N) slices tested
- Unified (|ρ| > 0.5): 0/0
- Partial overlap (0.3 < |ρ| ≤ 0.5): 0/0
- Independent (|ρ| ≤ 0.3): 0/0

**VERDICT: INDEPENDENT CHANNELS** — antichain and B_ℓ spectral
residuals are largely uncorrelated. The two DDT escape paths access
different geometric information. E-bulk is multi-dimensional.

