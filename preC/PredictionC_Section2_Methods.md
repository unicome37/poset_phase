# Paper Skeleton: Section 2 Methods

## 2. Ensemble, Observables, and Three-Tier Validation Design

### 2.1. Poset ensemble

The study uses the same finite causal poset ensemble as the companion papers [Prediction B, A], consisting of 8 families:

| Family | Description | HII regime | Primary role |
|--------|-------------|------------|-------------|
| `lorentzian_like_2d` | 2D causal diamond profile | High HII | Test subject |
| `lorentzian_like_3d` | 3D analogue, broader layers | Low HII | Comparison |
| `lorentzian_like_4d` | 4D analogue, flat layers | Lowest HII | Comparison |
| `KR_like` | 3-layer bipartite, random edges | Low HII (fixed layers) | Entropy benchmark |
| `multi_layer_random` | Multi-layer with random edges | Medium HII | Paired competitor |
| `transitive_percolation` | Random DAG via percolation | Low HII | Control |
| `absolute_layered` | Fixed-depth layered | Medium HII | Control |
| `interval_order` | Interval-representable posets | High HII | Structural test |

All families share the same generation framework (`generators.py`), with deterministic seeding: `seed = 1000 * N + sample_id`.

### 2.2. Exact entropy observable

The primary entropy observable is

$$
H(P) = \log |L(P)|,
$$

where $L(P)$ is the set of all linear extensions of poset $P$. Exact counting is performed via dynamic programming over antichains (`experiment.py`). For the confirmatory Tier 1 (N = 10–16), all 8 families are computed exactly. For Tier 2 (N = 30–48), exact counting is used for `lor2d`; for `mlr` at larger N, SIS approximation (2000 permutations, 3 independent runs, median) is used as fallback.

### 2.3. Hierarchy Integration Index (HII)

We define the Hierarchy Integration Index as a composite z-score:

$$
\text{HII}(P) = \frac{1}{5}\left(z_\ell + z_g + z_f - z_a - z_r\right)
$$

where the components are:

| Symbol | Observable | Physical meaning | Contribution |
|--------|-----------|-----------------|-------------|
| $z_\ell$ | `layer_count` | Number of layers in the longest chain decomposition | Positive: more layers → deeper hierarchy |
| $z_g$ | `mean_layer_gap` | Average gap between adjacent layers | Positive: larger gaps → more inter-layer structure |
| $z_f$ | `long_edge_fraction` | Fraction of edges spanning ≥2 layers | Positive: long-range causal connections |
| $z_a$ | `adjacent_edge_fraction` | Fraction of edges between adjacent layers only | Negative: local-only connections → shallow |
| $z_r$ | `reduction_edge_density` | Edge density in the transitive reduction | Negative: dense reduction → less constrained |

Each $z$-score is computed within the relevant analysis tier (Tier 1: across all 8 families at all N; Tier 2: across matched pairs).

**Important limitation**: As documented in Section 5, this five-component composite shows direction reversal within individual families. The cross-family negative correlation is robust, but within-family the composite HII can correlate positively with $\log H$. This motivates narrowing the predictor set to `layer_count` and `mean_layer_gap`.

### 2.4. Coarse-graining stability observable

The coarse-graining identity switch rate $\sigma_\text{CG}$ is defined as:

$$
\sigma_\text{CG}(P) = \frac{\text{number of CG trials where family identity changes}}{\text{total CG trials}}
$$

Coarse-graining proceeds by randomly merging pairs of causally adjacent elements and re-classifying the resulting poset using the full family scoring pipeline. This observable was computed in the companion pairwise compressibility duel framework as `cg_family_switch_rate`.

### 2.5. Control observables

Three structural observables serve as controls:

- `antichain_width` (aw): maximum antichain size, measuring "flatness"
- `comparable_fraction` (cf): fraction of element pairs connected by the partial order
- `geo_dim_eff` (geo_dim): effective geometric dimension from order-fraction estimator

These were the original Tier 1 controls and are retained for comparison, though the Simpson's Paradox analysis (Section 4) shows they are insufficient without explicit N control.

### 2.6. Three-tier validation design

The validation employs three tiers with complementary strengths:

| Tier | Dataset | N range | Analysis type | Confound control |
|------|---------|---------|--------------|-----------------|
| **1** | frozen_exact (320 samples, 8 families) | 10–16 | Partial correlation | Explicit N + family controls |
| **2** | Expanded matched pairs (34 pairs, lor2d vs mlr) | 30–48 | Δ correlation | N matched by pair design |
| **3** | CG stability duel (68 samples, lor2d vs mlr) | 30–48 | Predictor→outcome | N range shared |

**Design rationale**:
- Tier 1 provides the broadest family coverage but is susceptible to Simpson's Paradox without N control.
- Tier 2 uses matched pairs at identical N, eliminating the dominant confound by design.
- Tier 3 extends the chain from HII → log_H to HII → σ_CG, testing whether entropy compression translates into macroscopic stability.

### 2.7. Statistical methods

All correlation p-values are computed via permutation tests (2000 permutations, seed = 20260313), not parametric assumptions.

- **Partial correlations** (Tier 1): Pearson partial correlation with specified control sets; significance by 2000 permutation shuffles of the residuals.
- **Pairwise Δ** (Tier 2): For each matched pair (lor2d sample, mlr sample at same N), compute Δfeature = feature(lor2d) − feature(mlr) and Δlog_H similarly. Pearson r across pairs.
- **Single-predictor** (Tier 3): Pearson r between feature and cg_switch_rate across all 68 samples.

### 2.8. Data sources

| Tier | Input files | Output prefix |
|------|------------|--------------|
| 1 | `outputs_frozen_exact/raw_samples.csv` | `outputs_exploratory/prediction_c_comprehensive/tier1_*` |
| 2 | `outputs_exploratory/mlr_survivor_matched_lor2d_expanded/mlr_survivor_matched_pairs.csv` + residual CSV | `…/tier2_*` |
| 3 | `outputs_exploratory/pairwise_compressibility_duel_expanded/duel_results.csv` + submodel CSVs | `…/tier3_*` |

### 2.9. Near-wall sampling boundary

The matched-pair design (Tiers 2 and 3) requires finding MLR "survivors" — multi_layer_random samples that pass a structural similarity filter against lor2d. Empirically:

| N | MLR survivor hit rate |
|---|-----------------------|
| 30 | ~1/200 |
| 40 | ~1/2000 |
| 44 | ~1/5000 |
| 48 | ~1/15000 |
| 52 | ~1/60000 (1 survivor found) |
| 56 | 0/80000 |

These numbers define a hard sampling boundary. The three-tier results are therefore valid within their tested ranges but cannot be straightforwardly extended beyond N ≈ 50 without new methodologies.

### 2.10. Summary of Section 2

The methods section establishes:

1. HII is a clearly defined composite z-score with explicit physical motivation for each component.
2. Three tiers address different confound profiles, with complementary strengths.
3. All p-values are permutation-based, not parametric.
4. Near-wall sampling constrains the N-range of Tiers 2–3.
5. The composite HII has a known within-family limitation, motivating narrower predictors.

---

## Suggested Tables for Section 2

### Table 1
Family definitions, HII regime classification, and primary role in each tier.

### Table 2
HII component glossary: observable name, formula, physical meaning, sign in HII.

### Table 3
Three-tier design comparison: dataset, N range, analysis type, confound control mechanism, key output metric.

### Table 4
Near-wall sampling rates — MLR survivor hit rate by N.

---

## Writing Notes

- Introduce the Simpson's Paradox risk explicitly in this section (§2.6 design rationale) so Section 4 is not a surprise.
- Be precise about which tiers use exact vs SIS entropy — this is a credibility point inherited from Prediction B.
- The HII limitation should be flagged here, not first revealed in Section 5. This mirrors Prediction B's practice of front-loading methodological caveats.
- The near-wall boundary must appear in Methods, not just Discussion. Frame it as a design constraint, not an afterthought.
