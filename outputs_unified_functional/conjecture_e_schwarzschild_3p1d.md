# §4.1.42: 3+1D Schwarzschild Sprinkling Experiment

## Motivation

§4.1.29 tested DDT condition C2 in 1+1D Schwarzschild, finding only 2/27
features beyond density. The 1+1D limitation is severe: no transverse
spatial degrees of freedom exist. In 3+1D, the antichain (transverse)
channel — which is the strongest DDT escape path in de Sitter (21/21,
§4.1.28) — has full spatial structure to respond to non-uniform tidal
curvature.

## Experiment Design

- Geometry: **3+1D Schwarzschild** (t, r, θ, φ)
- rₛ values (2M): [0.0, 0.5, 1.0, 2.0, 4.0]
- Sizes: [128, 256, 512]
- Total realizations: 150
- r range: [2·rₛ, 12·rₛ] (well outside horizon)
- Volume-weighted sprinkling: √|g| = r²sinθ / √(1-rₛ/r)
- Causal relation: Schwarzschild null cone (tortoise + angular)
- Curvature proxy: Kretschner scalar K = 12rₛ²/r⁶
- New: **radial binning** (inner/mid/outer) for local curvature analysis

### Key physics

- Schwarzschild is **Ricci-flat** (R=0, Rμν=0)
- Kretschner scalar K = 48M²/r⁶ (tidal/Weyl curvature)
- DDT was proven under constant curvature → condition C2 relaxed here
- §4.1.28 showed antichain channel escapes DDT C1 (21/21 in de Sitter)
- Question: does the transverse channel also respond to **non-uniform** Weyl curvature?

## Q1: Global Features vs mean Kretschner K (pooled by N)

| N | n_causal_pairs | C_0 | C_1 | C_2 | bd_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy | w_max_ratio |
|---|------|------|------|------|------|------|------|------|------|------|------|------|
| 128 | **+0.389** | **+0.454** | **+0.384** | *+0.353* | **+0.454** | **+0.506** | **-0.506** | **-0.434** | *-0.314* | **-0.496** | **+0.413** | **-0.454** |
| 256 | **+0.392** | **+0.394** | **+0.398** | **+0.392** | **+0.394** | **+0.470** | **-0.470** | **-0.493** | *-0.314* | **-0.414** | **+0.403** | **-0.405** |
| 512 | **+0.375** | **+0.413** | **+0.375** | **+0.381** | **+0.413** | **+0.457** | **-0.457** | **-0.456** | **-0.430** | **-0.463** | **+0.383** | N/A |

## Q2: Features vs rₛ (pooled across all N)

| feature | Spearman ρ(feature, rₛ) | p-value |
|---------|------------------------|---------|
| n_causal_pairs | -0.739 | 3.26e-27 |
| C_0 | -0.669 | 8.66e-21 |
| C_1 | -0.742 | 1.63e-27 |
| C_2 | -0.763 | 7.46e-30 |
| bd_ratio | -0.814 | 1.00e-36 |
| layer_ratio | -0.576 | 1.24e-14 |
| mean_layer_width | +0.576 | 1.24e-14 |
| layer_width_std | +0.735 | 9.89e-27 |
| layer_width_cv | +0.680 | 1.16e-21 |
| max_layer_width_ratio | +0.812 | 2.09e-36 |
| layer_entropy | -0.819 | 1.33e-37 |
| w_max_ratio | +0.822 | 1.04e-25 |

## Q3: Density-Residual Analysis (OLS remove n_causal_pairs)

Does each feature's residual (after density removal) still correlate with mean_K?

| N | feature | raw ρ_S(feat, K) | resid ρ_S | |resid| | p_resid | verdict |
|---|---------|-----------------|-----------|--------|---------|---------|
| 128 | C_0 | +0.454 | -0.145 | 0.145 | 3.16e-01 | density-dominated |
| 128 | C_1 | +0.384 | -0.515 | 0.515 | 1.31e-04 | **BEYOND DENSITY** |
| 128 | C_2 | +0.353 | -0.562 | 0.562 | 2.14e-05 | **BEYOND DENSITY** |
| 128 | bd_ratio | +0.454 | -0.145 | 0.145 | 3.16e-01 | density-dominated |
| 128 | layer_ratio | +0.506 | +0.028 | 0.028 | 8.47e-01 | density-dominated |
| 128 | mean_layer_width | -0.506 | +0.149 | 0.149 | 3.02e-01 | density-dominated |
| 128 | layer_width_std | -0.434 | +0.098 | 0.098 | 4.97e-01 | density-dominated |
| 128 | layer_width_cv | -0.314 | +0.197 | 0.197 | 1.70e-01 | marginal |
| 128 | max_layer_width_ratio | -0.496 | +0.012 | 0.012 | 9.32e-01 | density-dominated |
| 128 | layer_entropy | +0.413 | -0.109 | 0.109 | 4.52e-01 | density-dominated |
| 128 | w_max_ratio | -0.454 | +0.062 | 0.062 | 6.67e-01 | density-dominated |
| 256 | C_0 | +0.394 | -0.084 | 0.084 | 5.60e-01 | density-dominated |
| 256 | C_1 | +0.398 | -0.276 | 0.276 | 5.24e-02 | marginal |
| 256 | C_2 | +0.392 | -0.425 | 0.425 | 2.07e-03 | **BEYOND DENSITY** |
| 256 | bd_ratio | +0.394 | -0.084 | 0.084 | 5.60e-01 | density-dominated |
| 256 | layer_ratio | +0.470 | -0.155 | 0.155 | 2.81e-01 | marginal |
| 256 | mean_layer_width | -0.470 | +0.245 | 0.245 | 8.64e-02 | marginal |
| 256 | layer_width_std | -0.493 | +0.142 | 0.142 | 3.26e-01 | density-dominated |
| 256 | layer_width_cv | -0.314 | +0.115 | 0.115 | 4.28e-01 | density-dominated |
| 256 | max_layer_width_ratio | -0.414 | +0.123 | 0.123 | 3.94e-01 | density-dominated |
| 256 | layer_entropy | +0.403 | -0.113 | 0.113 | 4.34e-01 | density-dominated |
| 256 | w_max_ratio | -0.405 | +0.033 | 0.033 | 8.22e-01 | density-dominated |
| 512 | C_0 | +0.413 | -0.025 | 0.025 | 8.65e-01 | density-dominated |
| 512 | C_1 | +0.375 | -0.291 | 0.291 | 4.00e-02 | marginal+ |
| 512 | C_2 | +0.381 | -0.419 | 0.419 | 2.46e-03 | **BEYOND DENSITY** |
| 512 | bd_ratio | +0.413 | -0.025 | 0.025 | 8.65e-01 | density-dominated |
| 512 | layer_ratio | +0.457 | -0.070 | 0.070 | 6.27e-01 | density-dominated |
| 512 | mean_layer_width | -0.457 | +0.095 | 0.095 | 5.13e-01 | density-dominated |
| 512 | layer_width_std | -0.456 | +0.021 | 0.021 | 8.82e-01 | density-dominated |
| 512 | layer_width_cv | -0.430 | +0.102 | 0.102 | 4.81e-01 | density-dominated |
| 512 | max_layer_width_ratio | -0.463 | +0.022 | 0.022 | 8.79e-01 | density-dominated |
| 512 | layer_entropy | +0.383 | -0.117 | 0.117 | 4.18e-01 | density-dominated |

**Features beyond density: 4/32**

- N=128: 2/11
- N=256: 1/11
- N=512: 1/10

## Q4: Radial Bin Analysis — Local Curvature Response

Each realization is split into 3 equal-count radial bins (inner/mid/outer).
Inner bins have stronger tidal curvature K(r). We test whether
antichain features differ between bins after controlling for local density.

| rₛ | N | inner K̄ | mid K̄ | outer K̄ | inner_lr | mid_lr | outer_lr |
|-----|---|---------|-------|---------|----------|--------|----------|
| 0.5 | 128 | 4.31e-02 | 2.92e-04 | 9.68e-05 | 0.143 | 0.128 | 0.121 |
| 0.5 | 256 | 3.76e-02 | 3.07e-04 | 9.53e-05 | 0.089 | 0.080 | 0.077 |
| 0.5 | 512 | 5.03e-02 | 3.00e-04 | 9.57e-05 | 0.054 | 0.049 | 0.047 |
| 1.0 | 128 | 2.95e-03 | 2.08e-05 | 6.29e-06 | 0.100 | 0.091 | 0.081 |
| 1.0 | 256 | 3.18e-03 | 1.85e-05 | 6.18e-06 | 0.056 | 0.055 | 0.051 |
| 1.0 | 512 | 2.96e-03 | 1.87e-05 | 6.14e-06 | 0.035 | 0.036 | 0.032 |
| 2.0 | 128 | 1.35e-04 | 1.10e-06 | 3.70e-07 | 0.071 | 0.067 | 0.060 |
| 2.0 | 256 | 1.61e-04 | 1.24e-06 | 3.87e-07 | 0.039 | 0.036 | 0.037 |
| 2.0 | 512 | 1.96e-04 | 1.14e-06 | 3.81e-07 | 0.024 | 0.025 | 0.022 |
| 4.0 | 128 | 1.26e-05 | 7.44e-08 | 2.36e-08 | 0.048 | 0.047 | 0.049 |
| 4.0 | 256 | 1.23e-05 | 7.66e-08 | 2.42e-08 | 0.025 | 0.028 | 0.027 |
| 4.0 | 512 | 1.13e-05 | 7.23e-08 | 2.33e-08 | 0.018 | 0.018 | 0.018 |

### Q4b: Pooled bin-level correlation

- Raw ρ(bin_K, bin_layer_ratio): +0.589 (p=5.98e-35)
- Density-residualized ρ: -0.045 (p=3.97e-01)
- Local curvature response not significant after density removal

## Q5: N-Scaling of Beyond-Density Signals

Do the beyond-density residuals strengthen with N?

| feature | N=128 |ρ_resid| | N=256 |ρ_resid| | N=512 |ρ_resid| | trend |
|---------|-----------------|-----------------|-----------------|-------|
| C_0 | 0.145 | 0.084 | 0.025 | ↓ |
| C_1 | 0.515 | 0.276 | 0.291 | ~ |
| C_2 | 0.562 | 0.425 | 0.419 | ↓ |
| bd_ratio | 0.145 | 0.084 | 0.025 | ↓ |
| layer_ratio | 0.028 | 0.155 | 0.070 | ~ |
| mean_layer_width | 0.149 | 0.245 | 0.095 | ~ |
| layer_width_std | 0.098 | 0.142 | 0.021 | ~ |
| layer_width_cv | 0.197 | 0.115 | 0.102 | ↓ |
| max_layer_width_ratio | 0.012 | 0.123 | 0.022 | ~ |
| layer_entropy | 0.109 | 0.113 | 0.117 | ↑ |
| w_max_ratio | 0.062 | 0.033 | N/A | ↓ |

## Comparison

| Setting | Beyond density | Peak |ρ_resid| |
|---------|--------------|----------------|
| **3+1D Schwarzschild** (this) | **4/32** | (see above) |
| 1+1D Schwarzschild (§4.1.29) | 2/27 | ~0.33 |
| de Sitter antichain (§4.1.28) | 21/21 | 0.817 |
| de Sitter B_ℓ spectral (§4.1.27) | 6/18 | 0.703 |

## Conclusion

**4/32** features carry curvature information
beyond density in 3+1D Schwarzschild. While modest, this is an improvement
over 1+1D (2/27). The weaker signal compared to de Sitter (21/21) is
expected: Schwarzschild is Ricci-flat (R=0), so only tidal (Weyl)
curvature exists — the antichain channel's primary response may be
to scalar curvature R, not Weyl.

### Physical Interpretation

- Schwarzschild: R=0 (Ricci-flat), K=48M²/r⁶ (Weyl/tidal)
- de Sitter: R=d(d-1)H² (constant scalar curvature), K∝R²
- The antichain channel in de Sitter (§4.1.28) responds to the
  de Sitter expansion that WIDENS spatial slices — a scalar R effect
- In Schwarzschild, there is no uniform expansion; curvature is tidal
- If beyond-density signals are found, they indicate the antichain
  channel also encodes Weyl curvature, strengthening DDT C2 escape

### 中文结论

**4/32** 个特征在 3+1D Schwarzschild 中超越密度。
信号弱于 de Sitter（21/21），因为 Schwarzschild 是 Ricci 平坦的——
只有潮汐（Weyl）曲率，而反链通道可能主要响应标量曲率 R。
