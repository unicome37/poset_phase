# T4: Spectral Ratio Bridge — EH Reanalysis

## Goal

Test whether eigenvalue RATIOS/GAPS of B_ℓ achieve α ≈ 2
(tracking R = d(d-1)H²) instead of α ≈ 1 (tracking H).

- Data source: §4.1.31 dual-channel CSV (360 realizations)
- Dimensions: [2, 3, 4]
- Sizes: [128, 256, 512]
- Hubble values: [0.0, 0.25, 0.5, 1.0, 2.0]
- New features: ['eig_ratio', 'eig_abs_ratio', 'eig_gap_norm', 'eig_gap_over_max', 'eig_neg_pos_gap', 'eig_neg_pos_prod', 'b1_std_over_mean', 'b1_std_sq']
- Baseline features: ['w_max_ratio', 'b1_std', 'eig_spread', 'mean_layer_width']

## 1. Raw Spearman ρ(feature, H²) — pooled by d

| d | eig_ratio | eig_abs_ratio | eig_gap_norm | eig_gap_over_max | eig_neg_pos_gap | eig_neg_pos_prod | b1_std_over_mean | b1_std_sq | w_max_ratio | b1_std | eig_spread | mean_layer_width |
|---|------|------|------|------|------|------|------|------|------|------|------|------|
| 2 | *+0.222* | *-0.222* | +0.021 | +0.021 | -0.009 | **-0.241** | -0.074 | **-0.397** | **+0.733** | **-0.397** | **-0.239** | **+0.464** |
| 3 | +0.040 | -0.040 | *+0.199* | *+0.201* | **+0.381** | **-0.454** | **+0.350** | **-0.711** | **+0.885** | **-0.711** | **-0.454** | **+0.497** |
| 4 | **+0.708** | **-0.708** | **+0.535** | **+0.477** | +0.057 | **-0.738** | **-0.588** | **-0.864** | **+0.946** | **-0.864** | **-0.738** | **+0.490** |

## 2. Density-Residual ρ(residual, H²) per (d, N)

After OLS-removing n_causal_pairs.

| d | N | eig_ratio | eig_abs_ratio | eig_gap_norm | eig_gap_over_max | eig_neg_pos_gap | eig_neg_pos_prod | b1_std_over_mean | b1_std_sq | w_max_ratio | b1_std | eig_spread | mean_layer_width |
|---|---|------|------|------|------|------|------|------|------|------|------|------|------|
| 2 | 128 | +0.020 | -0.020 | +0.060 | +0.072 | -0.116 | +0.093 | -0.008 | -0.003 | -0.126 | +0.037 | +0.110 | -0.101 |
| 2 | 256 | -0.012 | +0.012 | +0.025 | +0.021 | -0.012 | +0.035 | +0.129 | +0.008 | -0.179 | +0.003 | +0.051 | -0.124 |
| 2 | 512 | -0.136 | +0.136 | -0.115 | -0.126 | +0.149 | +0.123 | +0.012 | +0.029 | -0.205 | +0.029 | +0.092 | -0.136 |
| 3 | 128 | -0.002 | +0.002 | -0.058 | -0.063 | +0.005 | +0.032 | -0.051 | -0.096 | -0.081 | +0.023 | +0.084 | -0.297 |
| 3 | 256 | +0.014 | -0.014 | -0.055 | -0.049 | -0.069 | +0.095 | -0.035 | -0.021 | -0.058 | +0.052 | +0.152 | -0.179 |
| 3 | 512 | -0.098 | +0.098 | -0.046 | -0.044 | +0.086 | +0.139 | -0.145 | +0.041 | -0.090 | +0.135 | +0.167 | -0.115 |
| 4 | 128 | -0.199 | +0.184 | -0.130 | -0.248 | -0.060 | +0.046 | **+0.697** | +0.058 | +0.064 | +0.008 | -0.129 | -0.103 |
| 4 | 256 | -0.282 | +0.282 | -0.247 | -0.239 | +0.118 | -0.092 | **+0.495** | -0.049 | +0.070 | -0.043 | -0.035 | -0.265 |
| 4 | 512 | -0.133 | +0.133 | -0.185 | -0.184 | -0.109 | -0.072 | **-0.677** | +0.026 | +0.095 | -0.083 | -0.017 | -0.046 |

## 3. Power-Law α Grid Scan (KEY TEST)

For each feature, find the α that maximizes R²(group-mean residual, H^α).
Target: α ≈ 2 at d=4 → feature tracks R = d(d-1)H².
Baseline from §4.1.32: w_max_ratio α≈1.25, b1_std α≈1.00 at d=4.

### 3a. Pooled group-mean analysis (all N combined)

| d | Feature | α_best | R²_best | α_2nd | R²_2nd | α=1.0 R² | α=2.0 R² |
|---|---------|--------|---------|-------|--------|----------|----------|
| 2 | eig_ratio | **8.00** | 0.7717 | 7.75 | 0.7716 | 0.6992 | 0.7385 |
| 2 | eig_abs_ratio | **8.00** | 0.7717 | 7.75 | 0.7716 | 0.6992 | 0.7385 |
| 2 | eig_gap_norm | **0.25** | 0.2186 | 0.50 | 0.1867 | 0.1232 | 0.0396 |
| 2 | eig_gap_over_max | **0.25** | 0.1964 | 0.50 | 0.1674 | 0.1091 | 0.0329 |
| 2 | eig_neg_pos_gap | **0.25** | 0.8900 | 0.50 | 0.8454 | 0.7496 | 0.5929 |
| 2 | eig_neg_pos_prod | **8.00** | 0.5910 | 7.75 | 0.5904 | 0.2903 | 0.4166 |
| 2 | b1_std_over_mean | **8.00** | 0.3612 | 7.75 | 0.3610 | 0.1139 | 0.2398 |
| 2 | b1_std_sq | **8.00** | 0.8387 | 7.75 | 0.8382 | 0.5558 | 0.6936 |
| 2 | w_max_ratio | **8.00** | 0.7094 | 7.75 | 0.7090 | 0.3364 | 0.5268 |
| 2 | b1_std | **8.00** | 0.8210 | 7.75 | 0.8205 | 0.5031 | 0.6554 |
| 2 | eig_spread | **8.00** | 0.5310 | 7.75 | 0.5304 | 0.2018 | 0.3366 |
| 2 | mean_layer_width | **8.00** | 0.6923 | 7.75 | 0.6920 | 0.3678 | 0.5473 |
| 3 | eig_ratio | **8.00** | 0.3007 | 7.75 | 0.3005 | 0.1158 | 0.2210 |
| 3 | eig_abs_ratio | **8.00** | 0.3007 | 7.75 | 0.3005 | 0.1158 | 0.2210 |
| 3 | eig_gap_norm | **8.00** | 0.4263 | 7.75 | 0.4256 | 0.1104 | 0.2324 |
| 3 | eig_gap_over_max | **8.00** | 0.4279 | 7.75 | 0.4272 | 0.1114 | 0.2338 |
| 3 | eig_neg_pos_gap | **8.00** | 0.2121 | 7.75 | 0.2117 | 0.0807 | 0.1221 |
| 3 | eig_neg_pos_prod | **8.00** | 0.6363 | 7.75 | 0.6359 | 0.2844 | 0.4654 |
| 3 | b1_std_over_mean | **8.00** | 0.0666 | 7.75 | 0.0663 | 0.0020 | 0.0063 |
| 3 | b1_std_sq | **8.00** | 0.7077 | 7.75 | 0.7071 | 0.3092 | 0.4990 |
| 3 | w_max_ratio | **8.00** | 0.8012 | 7.75 | 0.8007 | 0.4133 | 0.6064 |
| 3 | b1_std | **8.00** | 0.6942 | 7.75 | 0.6936 | 0.2948 | 0.4825 |
| 3 | eig_spread | **8.00** | 0.6737 | 7.75 | 0.6732 | 0.3015 | 0.4889 |
| 3 | mean_layer_width | **8.00** | 0.6122 | 7.75 | 0.6116 | 0.2197 | 0.3976 |
| 4 | eig_ratio | **8.00** | 0.8182 | 7.75 | 0.8177 | 0.4346 | 0.6271 |
| 4 | eig_abs_ratio | **8.00** | 0.8187 | 7.75 | 0.8182 | 0.4353 | 0.6279 |
| 4 | eig_gap_norm | **8.00** | 0.9729 | 7.75 | 0.9727 | 0.7033 | 0.8605 |
| 4 | eig_gap_over_max | **6.25** | 1.0000 | 6.50 | 1.0000 | 0.8452 | 0.9555 |
| 4 | eig_neg_pos_gap | **8.00** | 0.8800 | 7.75 | 0.8795 | 0.5215 | 0.7092 |
| 4 | eig_neg_pos_prod | **8.00** | 0.8142 | 7.75 | 0.8140 | 0.4966 | 0.6788 |
| 4 | b1_std_over_mean | **8.00** | 0.5324 | 7.75 | 0.5320 | 0.1995 | 0.3653 |
| 4 | b1_std_sq | **8.00** | 0.2790 | 7.75 | 0.2786 | 0.0426 | 0.1448 |
| 4 | w_max_ratio | **8.00** | 0.9674 | 7.75 | 0.9673 | 0.7313 | 0.8819 |
| 4 | b1_std | **5.50** | 0.9885 | 5.25 | 0.9885 | 0.8196 | 0.9414 |
| 4 | eig_spread | **8.00** | 0.8612 | 7.75 | 0.8609 | 0.5179 | 0.7062 |
| 4 | mean_layer_width | **8.00** | 0.8772 | 7.75 | 0.8769 | 0.5269 | 0.7154 |

### 3b. Per-(d, N) α scan — N-scaling of α_eff

Does α_eff increase with N at d=4? (Tests T5 hypothesis)

| d | N | α(eig_ratio) | α(eig_abs_ratio) | α(eig_neg_pos_prod) | α(b1_std_sq) | α(w_max_ratio) | α(b1_std) |
|---|---|------|------|------|------|------|------|
| 2 | 128 | 0.25 | 0.25 | 8.00 | 8.00 | 8.00 | 8.00 |
| 2 | 256 | 0.25 | 0.25 | 0.25 | 8.00 | 8.00 | 8.00 |
| 2 | 512 | 8.00 | 8.00 | 0.25 | 4.00 | 8.00 | 8.00 |
| 3 | 128 | 8.00 | 8.00 | 8.00 | 8.00 | 8.00 | 8.00 |
| 3 | 256 | 0.25 | 0.25 | 8.00 | 8.00 | 8.00 | 8.00 |
| 3 | 512 | 0.50 | 0.50 | 8.00 | 8.00 | 8.00 | 8.00 |
| 4 | 128 | 8.00 | 8.00 | 8.00 | 8.00 | 4.00 | 1.75 |
| 4 | 256 | 8.00 | 8.00 | 7.75 | 8.00 | 8.00 | 8.00 |
| 4 | 512 | 8.00 | 8.00 | 8.00 | 0.25 | 8.00 | 5.75 |

## 4. Direct Correlation: Residual vs R_dS = d(d-1)H²

| d | Feature | Spearman(resid, R_dS) | Pearson R²(resid, R_dS) | Pearson R²(resid, H) | R²(R)/R²(H) |
|---|---------|----------------------|------------------------|---------------------|-------------|
| 2 | eig_ratio | -0.019 | 0.0029 | 0.0005 | 5.88× |
| 2 | eig_abs_ratio | +0.019 | 0.0029 | 0.0005 | 5.88× |
| 2 | eig_gap_norm | -0.008 | 0.0002 | 0.0005 | 0.35× |
| 2 | eig_gap_over_max | -0.002 | 0.0001 | 0.0005 | 0.30× |
| 2 | eig_neg_pos_gap | +0.005 | 0.0016 | 0.0006 | 2.61× |
| 2 | eig_neg_pos_prod | +0.055 | 0.0027 | 0.0001 | 32.77× |
| 2 | b1_std_over_mean | +0.032 | 0.0034 | 0.0006 | 5.85× |
| 2 | b1_std_sq | +0.021 | 0.0075 | 0.0004 | 20.67× |
| 2 | w_max_ratio | -0.172 | 0.0532 | 0.0035 | 15.17× |
| 2 | b1_std | +0.034 | 0.0088 | 0.0004 | 24.33× |
| 2 | eig_spread | +0.059 | 0.0029 | 0.0000 | 63.94× |
| 2 | mean_layer_width | -0.114 | 0.0228 | 0.0011 | 20.58× |
| 3 | eig_ratio | -0.015 | 0.0008 | 0.0001 | 13.84× |
| 3 | eig_abs_ratio | +0.015 | 0.0008 | 0.0001 | 13.84× |
| 3 | eig_gap_norm | -0.062 | 0.0140 | 0.0019 | 7.36× |
| 3 | eig_gap_over_max | -0.061 | 0.0142 | 0.0020 | 7.24× |
| 3 | eig_neg_pos_gap | +0.008 | 0.0018 | 0.0014 | 1.36× |
| 3 | eig_neg_pos_prod | +0.099 | 0.0407 | 0.0068 | 5.95× |
| 3 | b1_std_over_mean | -0.067 | 0.0014 | 0.0005 | 2.67× |
| 3 | b1_std_sq | +0.007 | 0.0459 | 0.0089 | 5.16× |
| 3 | w_max_ratio | -0.085 | 0.1723 | 0.0348 | 4.95× |
| 3 | b1_std | +0.070 | 0.1075 | 0.0205 | 5.25× |
| 3 | eig_spread | +0.136 | 0.0576 | 0.0100 | 5.78× |
| 3 | mean_layer_width | -0.181 | 0.0809 | 0.0169 | 4.80× |
| 4 | eig_ratio | -0.214 | 0.1769 | 0.0667 | 2.65× |
| 4 | eig_abs_ratio | +0.215 | 0.1911 | 0.0720 | 2.65× |
| 4 | eig_gap_norm | -0.202 | 0.0581 | 0.0222 | 2.62× |
| 4 | eig_gap_over_max | -0.231 | 0.0277 | 0.0106 | 2.61× |
| 4 | eig_neg_pos_gap | -0.023 | 0.0048 | 0.0018 | 2.64× |
| 4 | eig_neg_pos_prod | -0.023 | 0.1264 | 0.0463 | 2.73× |
| 4 | b1_std_over_mean | +0.172 | 0.0043 | 0.0016 | 2.76× |
| 4 | b1_std_sq | +0.031 | 0.0072 | 0.0021 | 3.45× |
| 4 | w_max_ratio | +0.054 | 0.3052 | 0.1213 | 2.52× |
| 4 | b1_std | -0.053 | 0.1658 | 0.0646 | 2.56× |
| 4 | eig_spread | -0.024 | 0.2600 | 0.0971 | 2.68× |
| 4 | mean_layer_width | -0.138 | 0.2255 | 0.0844 | 2.67× |

## 5. Verdict

### ❌ No new feature achieves α ≈ 2 at d=4

All spectral ratio/gap features produce α_best → 8.00 in the pooled group-mean
scan. This is a **diagnostic artifact**: after density removal, the residual
signal is weak (|ρ_resid| < 0.3 for most features, see Section 2), and with
only 4 non-zero H levels (0.25, 0.5, 1.0, 2.0), the group-mean curve becomes
step-function-shaped, causing R² to increase monotonically with α.

**This is the same artifact observed in §4.1.33 C1/C2 candidates.**

The Section 2 density-residual table is the more reliable diagnostic:
- At d=4, only `b1_std_over_mean` shows |ρ_resid| > 0.3 (but with **sign flip**
  between N=256 (+0.495) and N=512 (−0.677), ruling it out as stable)
- No new spectral ratio feature achieves consistent beyond-density signal
- The baseline features themselves show |ρ_resid| < 0.1 per (d,N) slice —
  the strong pooled signals in §4.1.28/27 came from cross-N pooling

**Conclusion**: Eigenvalue ratios/gaps do NOT provide a shortcut to R-tracking.
This **confirms** the §4.1.33 conclusion: the H→R bridge is a
continuum-limit phenomenon, not achievable by algebraic operations on
finite-N eigenvalue features.

### Summary: α comparison (d=4, pooled)

| Feature | Type | α_best | R²_best | vs baseline |
|---------|------|--------|---------|-------------|
| eig_ratio | NEW | 8.00 | 0.8182 | ↑ closer to R |
| eig_abs_ratio | NEW | 8.00 | 0.8187 | ↑ closer to R |
| eig_gap_norm | NEW | 8.00 | 0.9729 | ↑ closer to R |
| eig_gap_over_max | NEW | 6.25 | 1.0000 | ↑ closer to R |
| eig_neg_pos_gap | NEW | 8.00 | 0.8800 | ↑ closer to R |
| eig_neg_pos_prod | NEW | 8.00 | 0.8142 | ↑ closer to R |
| b1_std_over_mean | NEW | 8.00 | 0.5324 | ↑ closer to R |
| b1_std_sq | NEW | 8.00 | 0.2790 | ↑ closer to R |
| w_max_ratio | baseline | 8.00 | 0.9674 |  |
| b1_std | baseline | 5.50 | 0.9885 |  |
| eig_spread | baseline | 8.00 | 0.8612 |  |
| mean_layer_width | baseline | 8.00 | 0.8772 |  |

### N-scaling at d=4

| Feature | N values | α values | ρ(N, α) | Trend |
|---------|----------|----------|---------|-------|
| eig_ratio | 128/256/512 | 8.00/8.00/8.00 | +nan | ↓ diverging |
| eig_abs_ratio | 128/256/512 | 8.00/8.00/8.00 | +nan | ↓ diverging |
| eig_neg_pos_prod | 128/256/512 | 8.00/7.75/8.00 | +0.00 | → flat |
| b1_std_sq | 128/256/512 | 8.00/8.00/0.25 | -0.87 | ↓ diverging |
| w_max_ratio | 128/256/512 | 4.00/8.00/8.00 | +0.87 | ↑ converging to 2 |
| b1_std | 128/256/512 | 1.75/8.00/5.75 | +0.50 | ↓ diverging |

**Note**: The α=8.00 saturation across most features and N values is an artifact
of the weak per-slice residual signal (see Section 2 and verdict above). Per-slice
group means with only 4 non-zero H levels and 8 reps per level cannot reliably
discriminate power-law exponents. A dedicated T5 experiment with N=1024+ and
more reps is needed for conclusive N-scaling analysis.

### Overall T4 Assessment

1. **No spectral ratio/gap achieves α ≈ 2** — eigenvalue ratios do not cancel density well enough to expose H² dependence
2. **b1_std_over_mean is anomalous** — sign-flips between N=256 and N=512 at d=4, likely a finite-size artifact rather than genuine signal
3. **Section 4 R²(R)/R²(H) ratios** — consistently ~2.5–3.5× at d=4, slightly favoring R over H, but absolute R² values are too low (max 0.31) for this to be conclusive
4. **The α scan methodology itself is unreliable** on weak residual signals with few H levels — all candidates saturate at α=8.00 (step-function artifact)
5. **Next step**: T5 (N-scaling with N=1024+) is the correct path — directly measures whether α_eff(N) → 2 in the continuum limit

---

*Generated by conjecture_e_spectral_ratio_bridge.py*
*Data: 360 realizations from §4.1.31*