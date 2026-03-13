# Prediction C — Manuscript Draft: Section 2

## 2. Ensemble, Observables, and Three-Tier Validation Design

### 2.1 Poset Ensemble

We study the same eight-family finite causal poset ensemble introduced in the companion papers [Prediction B, A]. Each family is parametrized by the poset size $N$ (number of elements) and a deterministic seed. The families and their structural signatures are summarized in Table 1.

**Table 1.** Poset family glossary.

| Family | Abbreviated | Construction | Layer regime | Primary role |
|--------|-------------|-------------|-------------|-------------|
| lorentzian_like_2d | Lor2D | Sprinkling in a 2D causal diamond | Deep (many layers, large gaps) | Test subject |
| lorentzian_like_3d | Lor3D | Sprinkling in a 3D causal diamond | Moderate | Dimensional comparison |
| lorentzian_like_4d | Lor4D | Sprinkling in a 4D causal diamond | Shallow (few layers, small gaps) | Dimensional comparison |
| KR_like | KR | 3-layer bipartite graph with random edges | Fixed 3 layers | Entropy benchmark |
| multi_layer_random | MLR | Multi-layer DAG with random inter-layer edges | Variable | Paired competitor (Tiers 2–3) |
| transitive_percolation | TP | Random DAG via Erdős–Rényi percolation | Low | Structural control |
| absolute_layered | AL | Uniform-width layered poset | Fixed $\lfloor N/4 \rfloor$ layers | Structural control |
| interval_order | IO | Interval-representable partial order | Deep | Structural test |

All families share the generation framework in `generators.py`, with seeding convention $\text{seed} = 1000 \times N + \text{sample\_id}$ for reproducibility. The ensemble spans $N = 10$ to $56$ across the three analysis tiers.

### 2.2 Combinatorial Entropy

The primary outcome variable is the combinatorial entropy

$$H(P) \;\equiv\; \log\bigl|L(P)\bigr|,$$

where $L(P)$ is the set of all linear extensions of the poset $P$. Each linear extension is a total order on the $N$ elements that is compatible with the partial order. This quantity controls the entropic weight of $P$ in a Boltzmann-weighted causal set ensemble [3].

**Exact computation.** For $N \le 16$ (Tier 1), $|L(P)|$ is computed exactly via dynamic programming over antichains [4], as implemented in `entropy_exact.py`. The algorithm enumerates all antichains and uses the inclusion–exclusion formula $|L(P)| = \sum_{A \in \mathcal{A}(P)} (-1)^{|A|+1} \cdot f(A)$, where $\mathcal{A}(P)$ is the set of non-empty antichains.

**Approximate computation.** For $N > 16$ (Tiers 2–3), Lor2D samples use exact DP where computationally feasible; MLR samples at $N \ge 30$ use sequential importance sampling (SIS) with 2000 permutations and 3 independent runs, taking the median. SIS estimates were validated against exact values at $N = 16$ with mean absolute error $< 0.5\%$.

### 2.3 Hierarchy Integration Index (HII)

We define the Hierarchy Integration Index as an equal-weighted composite of five z-scored structural observables:

$$\text{HII}(P) \;=\; \frac{1}{5}\!\left(\,z_\ell + z_g + z_f - z_a - z_r\,\right),$$

where the components are listed in Table 2.

**Table 2.** HII component definitions.

| Symbol | Observable | Definition | Sign in HII | Physical rationale |
|--------|-----------|------------|-------------|-------------------|
| $z_\ell$ | `layer_count` | Number of layers in the longest-chain decomposition | $+$ | More layers $\Rightarrow$ deeper temporal hierarchy |
| $z_g$ | `mean_layer_gap` | Mean element gap between consecutive layers | $+$ | Larger gaps $\Rightarrow$ stronger inter-layer separation |
| $z_f$ | `long_edge_fraction` | Fraction of covering edges spanning $\ge 2$ layers | $+$ | Long-range causal connections $\Rightarrow$ tighter ordering constraints |
| $z_a$ | `adjacent_edge_fraction` | Fraction of covering edges between consecutive layers | $-$ | Local-only connections $\Rightarrow$ shallow structure |
| $z_r$ | `reduction_edge_density` | Edge density in the transitive reduction | $-$ | Dense reduction $\Rightarrow$ weaker per-edge ordering constraint |

Each $z$-score is computed within the relevant analysis pool (e.g., all 320 samples in Tier 1, or the 46 matched pairs in Tier 2). The equal weighting $1/5$ is a pre-registered design choice: no post-hoc tuning of component weights was performed.

**Known limitation.** The composite HII shows direction reversal within individual families: at fixed family and $N$, higher HII can correlate *positively* with $\log H$ (see Section 5.4). The robust negative correlation holds only in cross-family comparisons at fixed $N$. This motivates the component-level analysis in Section 5, which identifies `layer_count` and `mean_layer_gap` as the primary within-family predictors.

### 2.4 Individual Structural Observables

Beyond the five HII components, several additional observables are used in the component decomposition (Section 5):

- **`cover_density`**: fraction of all possible edges present in the transitive reduction.
- **`layer_signature_redundancy`**: information-theoretic measure of how predictable the layer structure is from the element ordering (higher $\Rightarrow$ more regular layering).
- **`layer_size_std`**: standard deviation of layer sizes in the chain decomposition.

These observables are computed by `matched_residual_freedom_check.py::residual_metrics()`.

### 2.5 Coarse-Graining Stability Observable

The coarse-graining identity switch rate is

$$\sigma_{\text{CG}}(P) \;=\; \frac{\#\{\text{CG trials where family identity changes}\}}{\#\{\text{total CG trials}\}}.$$

Each CG trial randomly deletes $1 - \alpha$ of the elements (with $\alpha \in \{0.8, 0.6\}$, corresponding to 80% and 60% retention ratios) and reclassifies the resulting sub-poset using the full family scoring pipeline. Classification uses nearest-centroid assignment against pre-computed family centroids ($K = 12$ centroids per family per $N$, with 3 repeats per ratio).

A low $\sigma_{\text{CG}}$ indicates that the poset's structural identity is *robust* to random element removal under the current classification protocol. This is a proxy for structural robustness, not a direct physical observable; its interpretation depends on the fidelity of the classifier.

### 2.6 Three-Tier Validation Design

The core statistical question is:

> At fixed poset size $N$, is deeper hierarchy integration associated with lower combinatorial entropy?

No single analysis design can cleanly isolate this relationship, because (i) $N$ acts as a confound that can reverse the structural relationship (Section 4), and (ii) exact entropy computation is restricted to small $N$. We therefore employ three complementary tiers:

**Table 3.** Three-tier design summary.

| Tier | Dataset | $N$ range | Design | Confound control | Key metric |
|------|---------|-----------|--------|-----------------|-----------|
| **1** | 320 samples (8 families $\times$ 4 sizes $\times$ 10 samples) | 10–16 | Partial correlation | Explicit $N$ + family controls | $r_{\text{partial}}(\text{HII}, \log H \mid N)$ |
| **2** | 46 matched Lor2D–MLR pairs | 30–56 | $\Delta$-correlation | $N$ matched by construction | $r(\Delta\text{HII}, \Delta\log H)$ |
| **3** | 92 samples (Lor2D + MLR) | 30–56 | Feature $\to$ outcome | Shared $N$ range | $r(\text{layer\_count}, \sigma_\text{CG})$ |

**Design rationale:**

*Tier 1* provides the broadest family coverage (all 8 families) with exact entropy values. It tests the partial correlation $r(\text{HII}, \log H)$ under progressively richer control sets, culminating in explicit $N$ control. The Simpson's Paradox diagnosis (Section 4) emerges from this tier.

*Tier 2* eliminates $N$ as a confound by construction: each $\Delta$ is computed within a matched pair at identical $N$. The MLR "survivors" are filtered to fall within the Lor2D structural window (antichain width + comparable fraction), ensuring that the pair members are structurally comparable. Three filter stringencies (P10–P90, P5–P95, P0–P100) are tested as a sensitivity analysis.

*Tier 3* extends the association chain from HII $\to$ $\log H$ to structural features $\to$ $\sigma_\text{CG}$. If deeper hierarchy compresses entropy, and lower entropy promotes CG identity stability, then hierarchy features should directly predict classification stability under the current protocol.

### 2.7 MLR Survivor Matching Protocol

The matched-pair design (Tier 2) requires identifying MLR samples whose structural properties fall within the Lor2D reference window. The protocol is:

1. Generate a Lor2D reference sample ($n_\text{ref}$ samples at the target $N$).
2. Compute `antichain_width` and `comparable_fraction` for each reference sample.
3. Define the acceptance window as the $[q_L, q_U]$ quantile range of each observable.
4. Generate MLR samples sequentially; accept those where *both* observables fall within the window.

**Table 4.** MLR survivor acceptance rates by $N$ and filter stringency.

| $N$ | P10–P90 | P5–P95 | P0–P100 |
|-----|---------|--------|---------|
| 30 | 12/1,708 (0.70%) | — | — |
| 40 | 12/2,716 (0.44%) | — | — |
| 44 | 6 / 2,075 (0.29%) | — | — |
| 48 | 4 / 3,468 (0.12%) | — | — |
| **52** | **1/60,000 (0.002%)** | **6/12,601 (0.048%)** | 8/1,978 (0.40%) |
| **56** | **0/80,000 (0%)** | **6/45,121 (0.013%)** | 8/333 (2.40%) |

The acceptance rate drops super-exponentially with $N$, defining a "near-wall boundary" beyond which the matched-pair design becomes impractical under stringent filters. The P5–P95 moderate filter was introduced specifically to rescue $N = 52$ and $56$ data while retaining meaningful structural similarity constraints.

### 2.8 Statistical Methods

All significance tests use non-parametric permutation procedures:

- **Permutation count**: $n_\text{perm} = 2000$ with fixed seed (20260313) for reproducibility.
- **Test statistic**: Pearson correlation coefficient $r$.
- **$p$-value**: fraction of permuted $|r|$ values exceeding the observed $|r|$.
- **Partial correlation**: computed by regressing both the feature and the target on the control variables, then correlating the residuals. Permutation shuffles are applied to the residual of the feature variable.

No Bonferroni or multiple-testing corrections are applied in the confirmatory analysis (Tiers 1–3); the three-tier design provides replication across independent datasets rather than correction within a single test. The component decomposition (Section 5) is explicitly labeled as exploratory.

### 2.9 Confirmatory vs. Exploratory Separation

Following the framework established in Prediction B, we designate one **primary confirmatory endpoint** per tier before any results are examined:

- **Tier 1 primary**: $r_\text{partial}(\text{HII}, \log H \mid N)$ — sign and significance.
- **Tier 2 primary**: $r(\Delta\text{HII}, \Delta\log H)$ at P5–P95 — sign, magnitude, and $p$.
- **Tier 3 primary** (classifier-contingent): $r(\text{layer\_count}, \sigma_\text{CG})$ — sign and magnitude, interpreted as supporting evidence under the current nearest-centroid classifier, not as a co-equal mechanistic pillar.

All other analyses (alternative control sets in Tier 1, per-component correlations, within-family regressions, switch enhancement) are designated **exploratory** and reported without multiple-testing correction. The confirmatory–exploratory boundary is strict: no finding from the exploratory set is promoted to headline status.

| Status | Analysis | What it tests |
|--------|----------|--------------|
| **Confirmatory** | Tier 1 partial correlation with $N$ control | Pre-registered direction: HII ↑ → $\log H$ ↓ |
| **Confirmatory** | Tier 2 matched-pair $\Delta$ correlation | Same prediction, $N$-free design |
| **Confirmatory** | Tier 3 `layer_count` → $\sigma_\text{CG}$ | Association chain extension (classifier-contingent) |
| *Exploratory* | Component decomposition (Section 5) | Which HII components drive the signal |
| *Exploratory* | Three-stringency sensitivity | Robustness to filter relaxation |
| *Exploratory* | Within-family direction analysis | HII behavior inside individual families |

### 2.10 Summary

The methods establish:

1. A five-component HII with explicit physical motivation and a pre-registered equal-weighting scheme.
2. Three complementary tiers that address different confound profiles, spanning $N = 10$ to $56$.
3. Non-parametric permutation-based inference throughout.
4. A quantified near-wall sampling boundary that constrains the matched-pair design.
5. A known limitation of the composite HII (within-family direction reversal), front-loaded as a methodological caveat rather than deferred to discussion.
