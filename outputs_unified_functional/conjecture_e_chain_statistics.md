# §4.1.30: Causal Chain Statistics Experiment

## Experiment Design

- Dimensions: [2, 3, 4]
- Sizes: [128, 256, 512]
- Hubble values: [0.0, 0.25, 0.5, 1.0, 1.5, 2.0]
- Total realizations: 432
- Chain features: N_2..N_5, longest chain, chain ratios, chain entropy
- DDT path I.4: do k-chain counts carry curvature info beyond {C_k}?

## Q1: Chain Features vs H² (pooled by d)

| d | N_3 | N_4 | N_5 | longest_chain | longest_chain_ratio | N3_over_N2 | N4_over_N3 | N5_over_N4 | N4_over_N2sq | chain_entropy | N3_per_N | N4_per_N |
|---|------|------|------|------|------|------|------|------|------|------|------|------|
| 2 | **-0.395** | **-0.450** | **-0.480** | **-0.515** | **-0.484** | **-0.551** | **-0.560** | **-0.563** | **-0.982** | **+0.561** | **-0.539** | **-0.546** |
| 3 | **-0.707** | **-0.746** | **-0.770** | **-0.760** | **-0.522** | **-0.833** | **-0.829** | **-0.797** | **-0.974** | **-0.367** | **-0.830** | **-0.831** |
| 4 | **-0.841** | **-0.848** | **-0.781** | **-0.831** | **-0.538** | **-0.904** | **-0.844** | **-0.674** | **-0.934** | **-0.894** | **-0.907** | **-0.888** |

## Q2: Chain Features vs Density (n_causal_pairs)

| d | N_3 | N_4 | N_5 | longest_chain | longest_chain_ratio | N3_over_N2 | N4_over_N3 | N5_over_N4 | N4_over_N2sq | chain_entropy | N3_per_N | N4_per_N |
|---|------|------|------|------|------|------|------|------|------|------|------|------|
| 2 | +0.995 | +0.985 | +0.977 | +0.949 | -0.636 | +0.959 | +0.954 | +0.951 | +0.358 | -0.953 | +0.963 | +0.960 |
| 3 | +0.982 | +0.968 | +0.954 | +0.909 | -0.327 | +0.916 | +0.910 | +0.890 | +0.592 | +0.181 | +0.921 | +0.918 |
| 4 | +0.987 | +0.967 | +0.889 | +0.941 | -0.010 | +0.948 | +0.930 | +0.824 | +0.812 | +0.946 | +0.958 | +0.949 |

## Q5: Density-Residual Analysis (OLS remove n_causal_pairs)

After removing density, does the residual still correlate with H²?

| d | N | feature | raw ρ(feat,H²) | ρ(feat,density) | resid ρ | p_resid | verdict |
|---|---|---------|----------------|-----------------|---------|---------|---------|
| 2 | 128 | N_3 | -0.985 | +0.997 | -0.021 | 8.87e-01 | density-dominated |
| 2 | 128 | N_4 | -0.984 | +0.995 | +0.017 | 9.10e-01 | density-dominated |
| 2 | 128 | N_5 | -0.982 | +0.991 | +0.046 | 7.57e-01 | density-dominated |
| 2 | 128 | longest_chain | -0.888 | +0.904 | -0.077 | 6.01e-01 | density-dominated |
| 2 | 128 | longest_chain_ratio | -0.888 | +0.904 | -0.077 | 6.01e-01 | density-dominated |
| 2 | 128 | N3_over_N2 | -0.984 | +0.994 | +0.000 | 1.00e+00 | density-dominated |
| 2 | 128 | N4_over_N3 | -0.980 | +0.985 | -0.009 | 9.53e-01 | density-dominated |
| 2 | 128 | N5_over_N4 | -0.972 | +0.975 | -0.006 | 9.67e-01 | density-dominated |
| 2 | 128 | N4_over_N2sq | -0.975 | +0.976 | +0.019 | 8.96e-01 | density-dominated |
| 2 | 128 | chain_entropy | +0.979 | -0.984 | +0.041 | 7.80e-01 | density-dominated |
| 2 | 128 | N3_per_N | -0.985 | +0.997 | -0.021 | 8.87e-01 | density-dominated |
| 2 | 128 | N4_per_N | -0.984 | +0.995 | +0.017 | 9.10e-01 | density-dominated |
| 2 | 256 | N_3 | -0.986 | +0.997 | -0.070 | 6.34e-01 | density-dominated |
| 2 | 256 | N_4 | -0.986 | +0.996 | -0.048 | 7.44e-01 | density-dominated |
| 2 | 256 | N_5 | -0.986 | +0.993 | -0.011 | 9.43e-01 | density-dominated |
| 2 | 256 | longest_chain | -0.942 | +0.944 | -0.069 | 6.43e-01 | density-dominated |
| 2 | 256 | longest_chain_ratio | -0.942 | +0.944 | -0.069 | 6.43e-01 | density-dominated |
| 2 | 256 | N3_over_N2 | -0.986 | +0.995 | +0.020 | 8.91e-01 | density-dominated |
| 2 | 256 | N4_over_N3 | -0.985 | +0.988 | -0.005 | 9.72e-01 | density-dominated |
| 2 | 256 | N5_over_N4 | -0.984 | +0.985 | -0.060 | 6.86e-01 | density-dominated |
| 2 | 256 | N4_over_N2sq | -0.984 | +0.982 | +0.018 | 9.01e-01 | density-dominated |
| 2 | 256 | chain_entropy | +0.985 | -0.987 | -0.079 | 5.92e-01 | density-dominated |
| 2 | 256 | N3_per_N | -0.986 | +0.997 | -0.070 | 6.34e-01 | density-dominated |
| 2 | 256 | N4_per_N | -0.986 | +0.996 | -0.048 | 7.44e-01 | density-dominated |
| 2 | 512 | N_3 | -0.986 | +0.998 | -0.145 | 3.24e-01 | density-dominated |
| 2 | 512 | N_4 | -0.986 | +0.997 | -0.142 | 3.36e-01 | density-dominated |
| 2 | 512 | N_5 | -0.986 | +0.994 | -0.066 | 6.56e-01 | density-dominated |
| 2 | 512 | longest_chain | -0.968 | +0.945 | -0.089 | 5.48e-01 | density-dominated |
| 2 | 512 | longest_chain_ratio | -0.968 | +0.945 | -0.089 | 5.48e-01 | density-dominated |
| 2 | 512 | N3_over_N2 | -0.986 | +0.996 | -0.135 | 3.61e-01 | density-dominated |
| 2 | 512 | N4_over_N3 | -0.986 | +0.991 | -0.094 | 5.24e-01 | density-dominated |
| 2 | 512 | N5_over_N4 | -0.986 | +0.987 | -0.085 | 5.64e-01 | density-dominated |
| 2 | 512 | N4_over_N2sq | -0.986 | +0.984 | -0.104 | 4.82e-01 | density-dominated |
| 2 | 512 | chain_entropy | +0.986 | -0.989 | -0.103 | 4.86e-01 | density-dominated |
| 2 | 512 | N3_per_N | -0.986 | +0.998 | -0.145 | 3.24e-01 | density-dominated |
| 2 | 512 | N4_per_N | -0.986 | +0.997 | -0.142 | 3.36e-01 | density-dominated |
| 3 | 128 | N_3 | -0.985 | +0.994 | +0.093 | 5.28e-01 | density-dominated |
| 3 | 128 | N_4 | -0.980 | +0.988 | +0.192 | 1.91e-01 | marginal |
| 3 | 128 | N_5 | -0.965 | +0.973 | +0.299 | 3.93e-02 | marginal |
| 3 | 128 | longest_chain | -0.906 | +0.891 | -0.019 | 8.96e-01 | density-dominated |
| 3 | 128 | longest_chain_ratio | -0.906 | +0.891 | -0.019 | 8.96e-01 | density-dominated |
| 3 | 128 | N3_over_N2 | -0.977 | +0.984 | -0.010 | 9.48e-01 | density-dominated |
| 3 | 128 | N4_over_N3 | -0.966 | +0.973 | -0.030 | 8.40e-01 | density-dominated |
| 3 | 128 | N5_over_N4 | -0.889 | +0.897 | +0.010 | 9.50e-01 | density-dominated |
| 3 | 128 | N4_over_N2sq | -0.949 | +0.951 | -0.024 | 8.73e-01 | density-dominated |
| 3 | 128 | chain_entropy | -0.870 | +0.870 | +0.055 | 7.08e-01 | density-dominated |
| 3 | 128 | N3_per_N | -0.985 | +0.994 | +0.093 | 5.28e-01 | density-dominated |
| 3 | 128 | N4_per_N | -0.980 | +0.988 | +0.192 | 1.91e-01 | marginal |
| 3 | 256 | N_3 | -0.986 | +0.997 | +0.072 | 6.26e-01 | density-dominated |
| 3 | 256 | N_4 | -0.986 | +0.990 | +0.136 | 3.55e-01 | density-dominated |
| 3 | 256 | N_5 | -0.984 | +0.982 | +0.163 | 2.69e-01 | marginal |
| 3 | 256 | longest_chain | -0.913 | +0.887 | -0.088 | 5.52e-01 | density-dominated |
| 3 | 256 | longest_chain_ratio | -0.913 | +0.887 | -0.088 | 5.52e-01 | density-dominated |
| 3 | 256 | N3_over_N2 | -0.986 | +0.990 | +0.062 | 6.77e-01 | density-dominated |
| 3 | 256 | N4_over_N3 | -0.984 | +0.976 | -0.016 | 9.15e-01 | density-dominated |
| 3 | 256 | N5_over_N4 | -0.959 | +0.942 | -0.026 | 8.63e-01 | density-dominated |
| 3 | 256 | N4_over_N2sq | -0.982 | +0.971 | +0.052 | 7.26e-01 | density-dominated |
| 3 | 256 | chain_entropy | -0.380 | +0.357 | +0.004 | 9.81e-01 | density-dominated |
| 3 | 256 | N3_per_N | -0.986 | +0.997 | +0.072 | 6.26e-01 | density-dominated |
| 3 | 256 | N4_per_N | -0.986 | +0.990 | +0.136 | 3.55e-01 | density-dominated |
| 3 | 512 | N_3 | -0.986 | +0.997 | +0.007 | 9.62e-01 | density-dominated |
| 3 | 512 | N_4 | -0.986 | +0.995 | +0.022 | 8.82e-01 | density-dominated |
| 3 | 512 | N_5 | -0.986 | +0.992 | +0.053 | 7.21e-01 | density-dominated |
| 3 | 512 | longest_chain | -0.962 | +0.955 | +0.031 | 8.35e-01 | density-dominated |
| 3 | 512 | longest_chain_ratio | -0.962 | +0.955 | +0.031 | 8.35e-01 | density-dominated |
| 3 | 512 | N3_over_N2 | -0.986 | +0.993 | +0.052 | 7.26e-01 | density-dominated |
| 3 | 512 | N4_over_N3 | -0.986 | +0.990 | -0.006 | 9.67e-01 | density-dominated |
| 3 | 512 | N5_over_N4 | -0.986 | +0.987 | +0.002 | 9.91e-01 | density-dominated |
| 3 | 512 | N4_over_N2sq | -0.986 | +0.986 | +0.061 | 6.82e-01 | density-dominated |
| 3 | 512 | chain_entropy | +0.628 | -0.621 | +0.059 | 6.90e-01 | density-dominated |
| 3 | 512 | N3_per_N | -0.986 | +0.997 | +0.007 | 9.62e-01 | density-dominated |
| 3 | 512 | N4_per_N | -0.986 | +0.995 | +0.022 | 8.82e-01 | density-dominated |
| 4 | 128 | N_3 | -0.978 | +0.979 | +0.295 | 4.18e-02 | marginal |
| 4 | 128 | N_4 | -0.950 | +0.945 | +0.479 | 5.72e-04 | **BEYOND DENSITY** |
| 4 | 128 | N_5 | -0.733 | +0.703 | +0.625 | 2.03e-06 | **BEYOND DENSITY** |
| 4 | 128 | longest_chain | -0.908 | +0.896 | -0.114 | 4.42e-01 | density-dominated |
| 4 | 128 | longest_chain_ratio | -0.908 | +0.896 | -0.114 | 4.42e-01 | density-dominated |
| 4 | 128 | N3_over_N2 | -0.969 | +0.961 | -0.005 | 9.72e-01 | density-dominated |
| 4 | 128 | N4_over_N3 | -0.914 | +0.904 | +0.155 | 3.38e-01 | marginal |
| 4 | 128 | N5_over_N4 | -0.611 | +0.525 | -0.016 | 9.36e-01 | density-dominated |
| 4 | 128 | N4_over_N2sq | -0.908 | +0.897 | +0.341 | 1.78e-02 | marginal |
| 4 | 128 | chain_entropy | -0.966 | +0.956 | -0.166 | 2.58e-01 | marginal |
| 4 | 128 | N3_per_N | -0.978 | +0.979 | +0.295 | 4.18e-02 | marginal |
| 4 | 128 | N4_per_N | -0.950 | +0.945 | +0.479 | 5.72e-04 | **BEYOND DENSITY** |
| 4 | 256 | N_3 | -0.986 | +0.994 | +0.250 | 8.65e-02 | marginal |
| 4 | 256 | N_4 | -0.977 | +0.975 | +0.452 | 1.27e-03 | **BEYOND DENSITY** |
| 4 | 256 | N_5 | -0.915 | +0.906 | +0.506 | 2.41e-04 | **BEYOND DENSITY** |
| 4 | 256 | longest_chain | -0.919 | +0.919 | -0.093 | 5.28e-01 | density-dominated |
| 4 | 256 | longest_chain_ratio | -0.919 | +0.919 | -0.093 | 5.28e-01 | density-dominated |
| 4 | 256 | N3_over_N2 | -0.975 | +0.979 | +0.195 | 1.85e-01 | marginal |
| 4 | 256 | N4_over_N3 | -0.950 | +0.944 | +0.095 | 5.37e-01 | density-dominated |
| 4 | 256 | N5_over_N4 | -0.722 | +0.710 | -0.148 | 4.11e-01 | density-dominated |
| 4 | 256 | N4_over_N2sq | -0.951 | +0.945 | +0.178 | 2.26e-01 | marginal |
| 4 | 256 | chain_entropy | -0.952 | +0.950 | -0.179 | 2.24e-01 | marginal |
| 4 | 256 | N3_per_N | -0.986 | +0.994 | +0.250 | 8.65e-02 | marginal |
| 4 | 256 | N4_per_N | -0.977 | +0.975 | +0.452 | 1.27e-03 | **BEYOND DENSITY** |
| 4 | 512 | N_3 | -0.986 | +0.997 | +0.081 | 5.84e-01 | density-dominated |
| 4 | 512 | N_4 | -0.987 | +0.989 | +0.100 | 4.97e-01 | density-dominated |
| 4 | 512 | N_5 | -0.969 | +0.969 | +0.187 | 2.04e-01 | marginal |
| 4 | 512 | longest_chain | -0.929 | +0.937 | -0.151 | 3.04e-01 | marginal |
| 4 | 512 | longest_chain_ratio | -0.929 | +0.937 | -0.151 | 3.04e-01 | marginal |
| 4 | 512 | N3_over_N2 | -0.986 | +0.990 | +0.223 | 1.28e-01 | marginal |
| 4 | 512 | N4_over_N3 | -0.976 | +0.971 | -0.015 | 9.20e-01 | density-dominated |
| 4 | 512 | N5_over_N4 | -0.963 | +0.963 | +0.102 | 5.21e-01 | density-dominated |
| 4 | 512 | N4_over_N2sq | -0.963 | +0.957 | +0.043 | 7.71e-01 | density-dominated |
| 4 | 512 | chain_entropy | -0.980 | +0.979 | -0.176 | 2.31e-01 | marginal |
| 4 | 512 | N3_per_N | -0.986 | +0.997 | +0.081 | 5.84e-01 | density-dominated |
| 4 | 512 | N4_per_N | -0.987 | +0.989 | +0.100 | 4.97e-01 | density-dominated |

**Total: 6/108 features BEYOND DENSITY**

## Conclusion

**6/108** (6%) chain features carry
 curvature information beyond density.

### Comparison with other DDT escape channels

| Channel | Path | Beyond density |
|---------|------|---------------|
| B_ℓ spectral | I.1 | 6/18 |
| Antichain transverse | I.2 | **21/21** |
| Schwarzschild | I.3 | 2/27 |
| **Chain statistics** | **I.4** | **6/108** |

### Physical Interpretation

Some chain features carry curvature signal beyond density.
 The surviving features likely encode path-correlation structure
 that goes beyond pairwise interval counting.

