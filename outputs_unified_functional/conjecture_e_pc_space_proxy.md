# Conjecture E — PC-Space Curvature Proxy (§4.1.26)

Loaded 360 realizations, p_0 through p_4


## Part 1: Per-dimension optimal PC2+ proxy

Strategy: within each d, fit PCA on standardized p_k, then find the unit vector in PC2+ subspace that maximizes |Spearman(proxy, H²)|.


### d = 2

PC1 variance: 56.7%

**Pooled optimal direction** (across all N): |ρ| = 0.142

Weights on PC2..PC5: -0.144, -0.004, +0.005, +0.990

Back-projected to raw p_k coefficients (normalized):
  [+0.714, +1.000, +0.552, +0.859, +0.810]

Closed-form proxy ∝ Σ_k w_k · p_k where w = [+0.714, +1.000, +0.552, +0.859, +0.810]

#### Per-N slice results

| N | PC2-opt ρ | p | bdg_d2c ρ | p | KL ρ | p | PC2 alone ρ | p |
|---|----------|---|----------|---|------|---|------------|---|
| 128 | +0.236 | 1.4e-01 | +0.175 | 2.8e-01 | +0.899 | 3.4e-15 | -0.234 | 1.5e-01 |
| 256 | +0.006 | 9.7e-01 | +0.202 | 2.1e-01 | -0.891 | 1.3e-14 | -0.006 | 9.7e-01 |
| 512 | +0.224 | 1.7e-01 | +0.287 | 7.3e-02 | -0.980 | 2.7e-28 | -0.207 | 2.0e-01 |

#### Per-N locally-optimized proxy

| N | local-opt ρ | pooled-opt ρ | bdg_d2c ρ | KL ρ |
|---|-----------|-------------|----------|------|
| 128 | +0.292 | +0.236 | +0.175 | +0.899 |
| 256 | +0.156 | +0.006 | +0.202 | -0.891 |
| 512 | +0.262 | +0.224 | +0.287 | -0.980 |

### d = 3

PC1 variance: 87.9%

**Pooled optimal direction** (across all N): |ρ| = 0.388

Weights on PC2..PC5: -0.595, -0.671, +0.018, +0.443

Back-projected to raw p_k coefficients (normalized):
  [+0.173, +1.000, +0.924, -0.155, -0.738]

Closed-form proxy ∝ Σ_k w_k · p_k where w = [+0.173, +1.000, +0.924, -0.155, -0.738]

#### Per-N slice results

| N | PC2-opt ρ | p | bdg_d2c ρ | p | KL ρ | p | PC2 alone ρ | p |
|---|----------|---|----------|---|------|---|------------|---|
| 128 | -0.112 | 4.9e-01 | +0.815 | 1.6e-10 | +0.974 | 4.2e-26 | +0.070 | 6.7e-01 |
| 256 | +0.700 | 5.0e-07 | +0.894 | 7.6e-15 | +0.980 | 2.7e-28 | -0.583 | 7.8e-05 |
| 512 | +0.764 | 9.6e-09 | +0.925 | 1.5e-17 | +0.980 | 2.7e-28 | -0.781 | 2.8e-09 |

#### Per-N locally-optimized proxy

| N | local-opt ρ | pooled-opt ρ | bdg_d2c ρ | KL ρ |
|---|-----------|-------------|----------|------|
| 128 | +0.401 | -0.112 | +0.815 | +0.974 |
| 256 | +0.720 | +0.700 | +0.894 | +0.980 |
| 512 | +0.828 | +0.764 | +0.925 | +0.980 |

### d = 4

PC1 variance: 83.1%

**Pooled optimal direction** (across all N): |ρ| = 0.751

Weights on PC2..PC5: +0.999, -0.023, +0.021, -0.010

Back-projected to raw p_k coefficients (normalized):
  [+1.000, +0.959, +0.968, +0.790, +0.686]

Closed-form proxy ∝ Σ_k w_k · p_k where w = [+1.000, +0.959, +0.968, +0.790, +0.686]

#### Per-N slice results

| N | PC2-opt ρ | p | bdg_d2c ρ | p | KL ρ | p | PC2 alone ρ | p |
|---|----------|---|----------|---|------|---|------------|---|
| 128 | +0.609 | 3.1e-05 | +nan | nan | +0.972 | 1.3e-23 | +0.609 | 3.1e-05 |
| 256 | +0.970 | 6.7e-25 | +0.974 | 4.2e-26 | +0.981 | 1.7e-28 | +0.965 | 9.5e-24 |
| 512 | +0.979 | 1.1e-27 | +0.974 | 4.2e-26 | +0.980 | 2.7e-28 | +0.979 | 1.1e-27 |

#### Per-N locally-optimized proxy

| N | local-opt ρ | pooled-opt ρ | bdg_d2c ρ | KL ρ |
|---|-----------|-------------|----------|------|
| 128 | +0.719 | +0.609 | +nan | +0.972 |
| 256 | +0.974 | +0.970 | +0.974 | +0.981 |
| 512 | +0.980 | +0.979 | +0.974 | +0.980 |

## Part 2: Summary — PC-optimal vs bdg_d2c vs KL

### Significant positive slices (ρ > 0, p < 0.05)

| Method | Count | Slices |
|--------|-------|--------|
| PC2-opt (pooled) | 5/9 | d=3N=256, d=3N=512, d=4N=128, d=4N=256, d=4N=512 |
| PC2-opt (local) | 6/9 | d=3N=128, d=3N=256, d=3N=512, d=4N=128, d=4N=256, d=4N=512 |
| bdg_d2c | 5/9 | d=3N=128, d=3N=256, d=3N=512, d=4N=256, d=4N=512 |
| KL(p||p_flat) | 7/9 | d=2N=128, d=3N=128, d=3N=256, d=3N=512, d=4N=128, d=4N=256, d=4N=512 |
| PC2 alone | 3/9 | d=4N=128, d=4N=256, d=4N=512 |

### Mean |ρ| across all 9 slices

| Method | Mean |ρ| | Min ρ | Max ρ |
|--------|---------|-------|-------|
| PC2-opt (pooled) | 0.511 | -0.112 | +0.979 |
| PC2-opt (local) | 0.593 | +0.156 | +0.980 |
| bdg_d2c | 0.656 | +0.175 | +0.974 |
| KL(p||p_flat) | 0.960 | -0.980 | +0.981 |
| PC2 alone | 0.493 | -0.781 | +0.979 |

## Part 3: N-convergence of PC-optimal proxy

Does the shape-space curvature signal improve with N?

| d | N=128 ρ | N=256 ρ | N=512 ρ | Trend |
|---|---------|---------|---------|-------|
| 2 | +0.236 | +0.006 | +0.224 | ↓ degrading |
| 3 | -0.112 | +0.700 | +0.764 | ↑ improving |
| 4 | +0.609 | +0.970 | +0.979 | ↑ improving |

## Part 4: Interpretation

### Two-line architecture for pure-causal curvature encoding

**Line A — Wall / Admissibility (PC1 = density mode)**
- Encodes causal connectivity density
- Anti-monotone with curvature (raw C_k) or positive (normalized p_k)
- Already implemented as sigmoid wall in F7/F8a
- 88.6% of total variance → dominant, robust, cross-dimensional

**Line B — Bulk Curvature Proxy (PC2+ = shape residual)**
- Encodes shape distortion of C_k distribution
- Positive correlation with H² (when density mode is removed)
- Only 11.4% of variance → subleading, weaker signal
- bdg_d2c is an impure projection (71.7% shape, 28.3% density)
- PC-optimal proxy is a pure shape projection

### Key insight
> Wall succeeded first because it rides the dominant mode.
> Bulk proxy is harder because it must extract a subleading residual.
> This is not an engineering accident — it is dictated by the information hierarchy
> of causal interval distributions under de Sitter expansion.

