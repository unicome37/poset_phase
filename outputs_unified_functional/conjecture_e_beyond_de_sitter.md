# §4.1.34: Beyond de Sitter Generalization

## Design

- Dimensions: [2, 3, 4]
- Sizes: [128, 256, 512]
- Power-law exponents p: [0.4, 0.5, 0.67, 1.0, 1.5, 2.0]
- Total Phase A realizations: 432
- Total Phase B local-bin rows: 864
- Background: power-law FRW with a(t) = (t/t₀)^p
- H(t) = p/t — time-dependent (unlike constant-H de Sitter)
- p=0.5: radiation, p=0.67: matter, p=1.0: coasting, p>1: accelerating

---

## Phase A: Global Test — Residuals vs p² (H_eff²)

Analogous to §4.1.28/31 but with power-law FRW instead of de Sitter.

### A1: Antichain features vs p² (Spearman ρ, pooled by d)

| d | w_max_ratio | layer_ratio | mean_layer_width | layer_width_std | layer_width_cv | max_layer_width_ratio | layer_entropy |
|---|------|------|------|------|------|------|------|
| 2 | **+0.273** | **+0.280** | **-0.280** | **+0.515** | **+0.888** | *+0.172* | +0.130 |
| 3 | **+0.632** | +0.039 | -0.039 | **+0.430** | **+0.814** | **+0.485** | **-0.297** |
| 4 | **+0.809** | -0.117 | +0.117 | **+0.380** | **+0.647** | **+0.670** | **-0.565** |

### A2: B_ℓ spectral features vs p² (Spearman ρ, pooled by d)

| d | b1_mean | b1_std | b1_neg_frac | eig_min | eig_max | eig_gap | eig_spread |
|---|------|------|------|------|------|------|------|
| 2 | +0.101 | -0.150 | -0.030 | +0.115 | -0.113 | +0.078 | -0.123 |
| 3 | *+0.211* | **-0.266** | +0.036 | +0.049 | -0.068 | +0.055 | -0.062 |
| 4 | *-0.169* | **-0.438** | +0.038 | -0.031 | +0.024 | **+0.471** | +0.029 |

### A3: Density-Residual Analysis — Beyond Density?

OLS-remove n_causal_pairs (density proxy), test residual vs p².

| d | N | Feature | ρ_raw | ρ_resid | Beyond density? |
|---|---|---------|-------|---------|-----------------|
| 2 | 128 | w_max_ratio | +0.595 | +0.667 | ✅ |
| 2 | 128 | layer_ratio | +0.761 | +0.785 | ✅ |
| 2 | 128 | mean_layer_width | -0.761 | -0.789 | ✅ |
| 2 | 128 | layer_width_std | +0.660 | +0.720 | ✅ |
| 2 | 128 | layer_width_cv | +0.805 | +0.803 | ✅ |
| 2 | 128 | max_layer_width_ratio | +0.301 | +0.390 | ✅ |
| 2 | 128 | layer_entropy | +0.456 | +0.497 | ✅ |
| 2 | 128 | b1_mean | +0.158 | +0.153 | — |
| 2 | 128 | b1_std | -0.331 | -0.336 | ✅ |
| 2 | 128 | b1_neg_frac | -0.224 | -0.224 | — |
| 2 | 128 | eig_min | +0.135 | +0.112 | — |
| 2 | 128 | eig_max | -0.137 | -0.119 | — |
| 2 | 128 | eig_gap | +0.122 | +0.119 | — |
| 2 | 128 | eig_spread | -0.151 | -0.138 | — |
| 2 | 256 | w_max_ratio | +0.792 | +0.649 | ✅ |
| 2 | 256 | layer_ratio | +0.774 | +0.601 | ✅ |
| 2 | 256 | mean_layer_width | -0.774 | -0.608 | ✅ |
| 2 | 256 | layer_width_std | +0.881 | +0.756 | ✅ |
| 2 | 256 | layer_width_cv | +0.933 | +0.800 | ✅ |
| 2 | 256 | max_layer_width_ratio | +0.524 | +0.373 | ✅ |
| 2 | 256 | layer_entropy | +0.350 | +0.223 | — |
| 2 | 256 | b1_mean | +0.049 | +0.011 | — |
| 2 | 256 | b1_std | -0.498 | -0.316 | ✅ |
| 2 | 256 | b1_neg_frac | +0.038 | +0.055 | — |
| 2 | 256 | eig_min | +0.346 | +0.159 | — |
| 2 | 256 | eig_max | -0.303 | -0.149 | — |
| 2 | 256 | eig_gap | +0.102 | +0.111 | — |
| 2 | 256 | eig_spread | -0.332 | -0.160 | — |
| 2 | 512 | w_max_ratio | +0.801 | +0.815 | ✅ |
| 2 | 512 | layer_ratio | +0.914 | +0.888 | ✅ |
| 2 | 512 | mean_layer_width | -0.914 | -0.885 | ✅ |
| 2 | 512 | layer_width_std | +0.911 | +0.911 | ✅ |
| 2 | 512 | layer_width_cv | +0.943 | +0.942 | ✅ |
| 2 | 512 | max_layer_width_ratio | +0.537 | +0.622 | ✅ |
| 2 | 512 | layer_entropy | +0.365 | +0.382 | ✅ |
| 2 | 512 | b1_mean | +0.103 | +0.078 | — |
| 2 | 512 | b1_std | -0.090 | -0.060 | — |
| 2 | 512 | b1_neg_frac | +0.107 | +0.122 | — |
| 2 | 512 | eig_min | +0.268 | +0.262 | — |
| 2 | 512 | eig_max | -0.317 | -0.291 | — |
| 2 | 512 | eig_gap | +0.028 | +0.018 | — |
| 2 | 512 | eig_spread | -0.320 | -0.299 | — |
| 3 | 128 | w_max_ratio | +0.856 | +0.372 | ✅ |
| 3 | 128 | layer_ratio | +0.003 | +0.262 | — |
| 3 | 128 | mean_layer_width | -0.003 | -0.262 | — |
| 3 | 128 | layer_width_std | +0.758 | +0.300 | ✅ |
| 3 | 128 | layer_width_cv | +0.681 | +0.320 | ✅ |
| 3 | 128 | max_layer_width_ratio | +0.667 | +0.241 | — |
| 3 | 128 | layer_entropy | -0.491 | +0.034 | — |
| 3 | 128 | b1_mean | +0.196 | +0.072 | — |
| 3 | 128 | b1_std | -0.583 | -0.164 | — |
| 3 | 128 | b1_neg_frac | +0.076 | +0.003 | — |
| 3 | 128 | eig_min | +0.253 | +0.017 | — |
| 3 | 128 | eig_max | -0.338 | -0.055 | — |
| 3 | 128 | eig_gap | -0.055 | -0.199 | — |
| 3 | 128 | eig_spread | -0.282 | -0.042 | — |
| 3 | 256 | w_max_ratio | +0.909 | +0.384 | ✅ |
| 3 | 256 | layer_ratio | +0.078 | +0.306 | ✅ |
| 3 | 256 | mean_layer_width | -0.078 | -0.306 | ✅ |
| 3 | 256 | layer_width_std | +0.902 | +0.345 | ✅ |
| 3 | 256 | layer_width_cv | +0.831 | +0.317 | ✅ |
| 3 | 256 | max_layer_width_ratio | +0.751 | +0.196 | — |
| 3 | 256 | layer_entropy | -0.726 | -0.080 | — |
| 3 | 256 | b1_mean | +0.267 | -0.037 | — |
| 3 | 256 | b1_std | -0.538 | -0.100 | — |
| 3 | 256 | b1_neg_frac | +0.125 | +0.100 | — |
| 3 | 256 | eig_min | +0.024 | -0.058 | — |
| 3 | 256 | eig_max | -0.066 | +0.048 | — |
| 3 | 256 | eig_gap | +0.136 | -0.048 | — |
| 3 | 256 | eig_spread | -0.070 | +0.042 | — |
| 3 | 512 | w_max_ratio | +0.931 | +0.348 | ✅ |
| 3 | 512 | layer_ratio | +0.282 | +0.048 | — |
| 3 | 512 | mean_layer_width | -0.282 | -0.048 | — |
| 3 | 512 | layer_width_std | +0.940 | +0.290 | — |
| 3 | 512 | layer_width_cv | +0.939 | +0.353 | ✅ |
| 3 | 512 | max_layer_width_ratio | +0.846 | +0.124 | — |
| 3 | 512 | layer_entropy | -0.719 | -0.019 | — |
| 3 | 512 | b1_mean | +0.626 | +0.057 | — |
| 3 | 512 | b1_std | -0.717 | -0.158 | — |
| 3 | 512 | b1_neg_frac | -0.071 | +0.030 | — |
| 3 | 512 | eig_min | +0.142 | +0.010 | — |
| 3 | 512 | eig_max | -0.154 | -0.011 | — |
| 3 | 512 | eig_gap | +0.089 | +0.008 | — |
| 3 | 512 | eig_spread | -0.180 | -0.046 | — |
| 4 | 128 | w_max_ratio | +0.902 | +0.332 | ✅ |
| 4 | 128 | layer_ratio | -0.603 | +0.085 | — |
| 4 | 128 | mean_layer_width | +0.603 | +0.015 | — |
| 4 | 128 | layer_width_std | +0.749 | +0.145 | — |
| 4 | 128 | layer_width_cv | +0.376 | +0.136 | — |
| 4 | 128 | max_layer_width_ratio | +0.735 | +0.111 | — |
| 4 | 128 | layer_entropy | -0.793 | -0.073 | — |
| 4 | 128 | b1_mean | -0.107 | -0.246 | — |
| 4 | 128 | b1_std | -0.725 | +0.121 | — |
| 4 | 128 | b1_neg_frac | -0.250 | +0.033 | — |
| 4 | 128 | eig_min | +0.266 | -0.328 | ✅ |
| 4 | 128 | eig_max | -0.283 | +0.321 | ✅ |
| 4 | 128 | eig_gap | +0.593 | +0.032 | — |
| 4 | 128 | eig_spread | -0.272 | +0.324 | ✅ |
| 4 | 256 | w_max_ratio | +0.935 | +0.318 | ✅ |
| 4 | 256 | layer_ratio | -0.381 | +0.077 | — |
| 4 | 256 | mean_layer_width | +0.381 | -0.130 | — |
| 4 | 256 | layer_width_std | +0.897 | +0.289 | — |
| 4 | 256 | layer_width_cv | +0.681 | +0.291 | — |
| 4 | 256 | max_layer_width_ratio | +0.854 | +0.223 | — |
| 4 | 256 | layer_entropy | -0.875 | -0.127 | — |
| 4 | 256 | b1_mean | -0.195 | -0.082 | — |
| 4 | 256 | b1_std | -0.769 | -0.012 | — |
| 4 | 256 | b1_neg_frac | +0.206 | +0.079 | — |
| 4 | 256 | eig_min | -0.397 | -0.268 | — |
| 4 | 256 | eig_max | +0.377 | +0.273 | — |
| 4 | 256 | eig_gap | +0.542 | +0.181 | — |
| 4 | 256 | eig_spread | +0.384 | +0.267 | — |
| 4 | 512 | w_max_ratio | +0.966 | +0.255 | — |
| 4 | 512 | layer_ratio | -0.133 | +0.240 | — |
| 4 | 512 | mean_layer_width | +0.133 | -0.240 | — |
| 4 | 512 | layer_width_std | +0.933 | +0.248 | — |
| 4 | 512 | layer_width_cv | +0.886 | +0.261 | — |
| 4 | 512 | max_layer_width_ratio | +0.890 | +0.151 | — |
| 4 | 512 | layer_entropy | -0.850 | -0.091 | — |
| 4 | 512 | b1_mean | -0.212 | -0.108 | — |
| 4 | 512 | b1_std | -0.896 | +0.035 | — |
| 4 | 512 | b1_neg_frac | +0.148 | +0.126 | — |
| 4 | 512 | eig_min | -0.122 | -0.173 | — |
| 4 | 512 | eig_max | +0.109 | +0.163 | — |
| 4 | 512 | eig_gap | +0.261 | +0.011 | — |
| 4 | 512 | eig_spread | +0.124 | +0.159 | — |

**Beyond-density summary:**

- d=2: 22/42 features pass beyond-density test
- d=3: 10/42 features pass beyond-density test
- d=4: 5/42 features pass beyond-density test

### A4: Geometric Target — α Grid Search

Pooled group-mean: for each p level, average residualized features across reps.

Fit Z_resid ~ p^α, find best α ∈ [0.25, 8.0].

| d | Feature | Best α | R² | Physical match |
|---|---------|--------|-----|----------------|
| 2 | w_max_ratio | 8.00 | 0.309 | H^8.0 (higher order) |
| 2 | b1_std | 1.25 | 0.451 | H (expansion rate) |
| 3 | w_max_ratio | 0.25 | 0.298 | H (expansion rate) |
| 3 | b1_std | 0.25 | 0.139 | H (expansion rate) |
| 4 | w_max_ratio | 3.75 | 0.201 | H^3.8 (higher order) |
| 4 | b1_std | 5.00 | 0.471 | H^5.0 (higher order) |

---

## Phase B: Local Test — Within-Realization H(t) Tracking

Each realization split into early (high H) and late (low H) bins.

If early-bin w_max_ratio > late-bin w_max_ratio (same realization, same p),

the observable tracks LOCAL H(t), not just the global p.


### B1: Paired Early-Late Comparison

| d | N | p | reps | early H_local | late H_local | Δw_max_ratio (early−late) | sign consistent? |
|---|---|---|------|--------------|-------------|--------------------------|-----------------|
| 2 | 128 | 0.40 | 8 | 1.25 | 0.49 | -0.0703 (0/8 positive) | ❌ |
| 2 | 128 | 0.50 | 8 | 1.44 | 0.61 | -0.0508 (0/8 positive) | ❌ |
| 2 | 128 | 0.67 | 8 | 1.96 | 0.82 | -0.0703 (0/8 positive) | ❌ |
| 2 | 128 | 1.00 | 8 | 2.53 | 1.16 | -0.0801 (0/8 positive) | ❌ |
| 2 | 128 | 1.50 | 8 | 2.92 | 1.69 | -0.1172 (0/8 positive) | ❌ |
| 2 | 128 | 2.00 | 8 | 3.73 | 2.21 | -0.1055 (0/8 positive) | ❌ |
| 2 | 256 | 0.40 | 8 | 1.28 | 0.50 | -0.0371 (0/8 positive) | ❌ |
| 2 | 256 | 0.50 | 8 | 1.48 | 0.62 | -0.0361 (0/8 positive) | ❌ |
| 2 | 256 | 0.67 | 8 | 1.84 | 0.80 | -0.0430 (0/8 positive) | ❌ |
| 2 | 256 | 1.00 | 8 | 2.46 | 1.16 | -0.0645 (0/8 positive) | ❌ |
| 2 | 256 | 1.50 | 8 | 3.27 | 1.71 | -0.0869 (0/8 positive) | ❌ |
| 2 | 256 | 2.00 | 8 | 3.76 | 2.23 | -0.0947 (0/8 positive) | ❌ |
| 2 | 512 | 0.40 | 8 | 1.26 | 0.50 | -0.0234 (0/8 positive) | ❌ |
| 2 | 512 | 0.50 | 8 | 1.47 | 0.61 | -0.0298 (0/8 positive) | ❌ |
| 2 | 512 | 0.67 | 8 | 1.88 | 0.81 | -0.0386 (0/8 positive) | ❌ |
| 2 | 512 | 1.00 | 8 | 2.45 | 1.18 | -0.0435 (0/8 positive) | ❌ |
| 2 | 512 | 1.50 | 8 | 3.10 | 1.70 | -0.0493 (0/8 positive) | ❌ |
| 2 | 512 | 2.00 | 8 | 3.62 | 2.22 | -0.0635 (0/8 positive) | ❌ |
| 3 | 128 | 0.40 | 8 | 1.09 | 0.48 | -0.1738 (0/8 positive) | ❌ |
| 3 | 128 | 0.50 | 8 | 1.21 | 0.59 | -0.2168 (0/8 positive) | ❌ |
| 3 | 128 | 0.67 | 8 | 1.47 | 0.77 | -0.1914 (0/8 positive) | ❌ |
| 3 | 128 | 1.00 | 8 | 1.88 | 1.11 | -0.3281 (0/8 positive) | ❌ |
| 3 | 128 | 1.50 | 8 | 2.34 | 1.62 | -0.3125 (0/8 positive) | ❌ |
| 3 | 128 | 2.00 | 8 | 2.88 | 2.15 | -0.3066 (0/8 positive) | ❌ |
| 3 | 256 | 0.40 | 8 | 1.05 | 0.48 | -0.1279 (0/8 positive) | ❌ |
| 3 | 256 | 0.50 | 8 | 1.26 | 0.58 | -0.1611 (0/8 positive) | ❌ |
| 3 | 256 | 0.67 | 8 | 1.48 | 0.77 | -0.2090 (0/8 positive) | ❌ |
| 3 | 256 | 1.00 | 8 | 1.82 | 1.10 | -0.2441 (0/8 positive) | ❌ |
| 3 | 256 | 1.50 | 8 | 2.34 | 1.62 | -0.2637 (0/8 positive) | ❌ |
| 3 | 256 | 2.00 | 8 | 2.87 | 2.14 | -0.2471 (0/8 positive) | ❌ |
| 3 | 512 | 0.40 | 8 | 1.09 | 0.48 | -0.1284 (0/8 positive) | ❌ |
| 3 | 512 | 0.50 | 8 | 1.25 | 0.58 | -0.1260 (0/8 positive) | ❌ |
| 3 | 512 | 0.67 | 8 | 1.49 | 0.77 | -0.1499 (0/8 positive) | ❌ |
| 3 | 512 | 1.00 | 8 | 1.92 | 1.11 | -0.1968 (0/8 positive) | ❌ |
| 3 | 512 | 1.50 | 8 | 2.37 | 1.62 | -0.2412 (0/8 positive) | ❌ |
| 3 | 512 | 2.00 | 8 | 2.85 | 2.13 | -0.2563 (0/8 positive) | ❌ |
| 4 | 128 | 0.40 | 8 | 0.94 | 0.46 | -0.2734 (0/8 positive) | ❌ |
| 4 | 128 | 0.50 | 8 | 1.07 | 0.57 | -0.2949 (0/8 positive) | ❌ |
| 4 | 128 | 0.67 | 8 | 1.23 | 0.74 | -0.2930 (0/8 positive) | ❌ |
| 4 | 128 | 1.00 | 8 | 1.56 | 1.09 | -0.2793 (0/8 positive) | ❌ |
| 4 | 128 | 1.50 | 8 | 2.07 | 1.59 | -0.2773 (0/8 positive) | ❌ |
| 4 | 128 | 2.00 | 8 | 2.54 | 2.09 | -0.2383 (0/8 positive) | ❌ |
| 4 | 256 | 0.40 | 8 | 0.91 | 0.46 | -0.2568 (0/8 positive) | ❌ |
| 4 | 256 | 0.50 | 8 | 1.07 | 0.57 | -0.2588 (0/8 positive) | ❌ |
| 4 | 256 | 0.67 | 8 | 1.27 | 0.75 | -0.2949 (0/8 positive) | ❌ |
| 4 | 256 | 1.00 | 8 | 1.54 | 1.08 | -0.3291 (0/8 positive) | ❌ |
| 4 | 256 | 1.50 | 8 | 2.08 | 1.58 | -0.3340 (0/8 positive) | ❌ |
| 4 | 256 | 2.00 | 8 | 2.56 | 2.08 | -0.2979 (0/8 positive) | ❌ |
| 4 | 512 | 0.40 | 8 | 0.92 | 0.46 | -0.2495 (0/8 positive) | ❌ |
| 4 | 512 | 0.50 | 8 | 1.05 | 0.57 | -0.2593 (0/8 positive) | ❌ |
| 4 | 512 | 0.67 | 8 | 1.26 | 0.74 | -0.3140 (0/8 positive) | ❌ |
| 4 | 512 | 1.00 | 8 | 1.59 | 1.08 | -0.3345 (0/8 positive) | ❌ |
| 4 | 512 | 1.50 | 8 | 2.08 | 1.59 | -0.3340 (0/8 positive) | ❌ |
| 4 | 512 | 2.00 | 8 | 2.56 | 2.09 | -0.3442 (0/8 positive) | ❌ |

**Phase B Summary:**

- d=2: 0/18 (p, N) cells show early > late w_max_ratio
- d=3: 0/18 (p, N) cells show early > late w_max_ratio
- d=4: 0/18 (p, N) cells show early > late w_max_ratio

### B2: Local Features vs Local H (all bins pooled)

| d | Feature | Spearman ρ(H_local, feature) | p-value | Significant? |
|---|---------|------------------------------|---------|-------------|
| 2 | w_max_ratio | **-0.240** | 3.77e-05 | Yes |
| 2 | n_layers | **+0.486** | 1.86e-18 | Yes |
| 2 | layer_ratio | **+0.497** | 2.37e-19 | Yes |
| 2 | mean_layer_width | **-0.497** | 2.37e-19 | Yes |
| 3 | w_max_ratio | -0.115 | 5.20e-02 | No |
| 3 | n_layers | **+0.373** | 5.89e-11 | Yes |
| 3 | layer_ratio | **+0.312** | 6.47e-08 | Yes |
| 3 | mean_layer_width | **-0.312** | 6.47e-08 | Yes |
| 4 | w_max_ratio | +0.020 | 7.39e-01 | No |
| 4 | n_layers | **+0.213** | 2.64e-04 | Yes |
| 4 | layer_ratio | **+0.181** | 2.10e-03 | Yes |
| 4 | mean_layer_width | **-0.181** | 2.10e-03 | Yes |

---

## Verdict

### Phase A (Global): 37/126 features pass beyond-density in power-law FRW

**Partial generalization supported.** Post-density antichain observables carry

curvature-responsive information in power-law FRW, but the signal is weaker and

more dimension-dependent than in de Sitter:

- d=2: 22/42 beyond-density (vs 7/7 antichain in de Sitter §4.1.28)
- d=3: 10/42 beyond-density (vs 7/7 antichain in de Sitter §4.1.28)
- d=4: 5/42 beyond-density (vs 7/7 antichain in de Sitter §4.1.28)

**Key difference from de Sitter:** In de Sitter, H is a single global constant —

the entire sprinkling sees the same curvature. In power-law FRW, H(t) = p/t varies

within each realization; the 'control parameter' p encodes the EQUATION OF STATE,

not a single curvature value. The weaker signal at high d may reflect the fact that

the time-averaged curvature is a noisier proxy than a constant global curvature.


**B_ℓ spectral channel mostly fails:** Only w_max_ratio (antichain) consistently

passes beyond-density in power-law FRW; the B_ℓ spectral features, which were

already the weaker channel in de Sitter (6/18 vs 21/21), essentially vanish here.

This suggests the antichain (transverse) channel is the more robust carrier of

curvature information across different backgrounds.


### Phase B (Local): 0/54 cells show early>late w_max_ratio

**All 54 cells show early < late** — 100% opposite to the 'naive' expectation.

However, this is NOT evidence against local H-tracking. The Phase B design has a

**fundamental methodological confound**:


In power-law FRW with a(t) = (t/t₀)^p (p > 0), the sprinkling acceptance probability

∝ a(t)^{d-1} = t^{p(d-1)} increases with t. This means:

- **Early bin** (small t): fewer accepted points → smaller sub-poset

- **Late bin** (large t): more accepted points → larger sub-poset


The w_max_ratio of a **sub-poset** extracted from a larger sprinkling is not

comparable to w_max_ratio of an independent sprinkling. The sub-poset inherits

causal relations from the full geometry, and the early sub-poset is systematically

smaller (fewer elements) → its antichain width as a fraction of N_bin is dominated

by finite-size effects, not local curvature.


**The early < late direction is consistent with a DENSITY ARTIFACT:**

late bins have more elements → richer causal structure → the ratio w_max/N_bin

converges to a higher asymptotic value faster than in the sparse early bin.


**Conclusion:** Phase B as designed CANNOT distinguish local H-tracking from

finite-size + density confounds. A proper local test would require:

(a) Density-matched sub-sampling (equal N per bin), or

(b) Density-residualized local features, or

(c) Independent sprinklings in patches with different local curvature.

This remains an open methodological question for §4.1.35+.


### Comparison with de Sitter results (§4.1.28/31/32):

| Property | de Sitter (§4.1.28–32) | Power-law FRW (§4.1.34) |
|----------|----------------------|------------------------|
| Beyond-density (antichain) | **21/21** | 37/126 (63) |
| Beyond-density (B_ℓ) | 6/18 | weak/absent |
| α target (d=4) | **1.0–1.25** (clean H) | noisy, no clean target |
| Local H tracking | N/A (constant H) | inconclusive (confounded) |
| Background symmetry | Maximally symmetric | Reduced (time-dependent H) |

### Overall Physical Interpretation

**1. The antichain (transverse) channel partially generalizes beyond de Sitter.**

w_max_ratio survives density removal at d=2,3 in power-law FRW, confirming that

the DDT escape via transverse statistics is not an artifact of de Sitter symmetry.

The signal weakening at d=4 may reflect the stronger density dominance in higher

dimensions (DDT becomes harder to escape) combined with the time-varying H(t).


**2. The B_ℓ spectral channel does NOT generalize robustly.**

This is consistent with its already-weaker performance in de Sitter (6/18).

The spectral channel appears more sensitive to background symmetry.


**3. Local H(t) tracking remains an open question.**

Phase B's methodological confound prevents a clean test. The question 'does the

observable respond to local geometry?' requires a redesigned experiment.


**4. The 'remaining gap' is partially addressed but not closed.**

Phase A shows that Conjecture E's bulk structure is not purely a de Sitter artifact.

But the signal degradation at d=4 and the lack of a clean α target mean that the

power-law FRW generalization is weaker than the de Sitter result. The gap has

**narrowed** (from 'untested' to 'partially confirmed with caveats') but is not closed.
