# Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study

---

### Abstract

Companion papers (Predictions A and B) established that Lorentzian-like posets become competitive phases in action-weighted causal poset ensembles. What structural feature drives this advantage remains unexplained. We define a Hierarchy Integration Index (HII) — a pre-registered composite z-score of five structural observables — and test whether hierarchy depth predicts combinatorial entropy $\log H$ across three independent tiers: (i) all-family exact computation (8 families, $N = 10$–$16$, 320 samples), (ii) matched-pair $\Delta$-analysis (Lor2D vs MLR, 46 pairs, $N = 30$–$56$), and (iii) coarse-graining identity-stability linkage (92 samples). At fixed $N$, deeper hierarchy correlates with lower $\log H$: partial $r = -0.578$ controlling for $N$; matched-pair $r = -0.834$, stable across three filter stringencies (variation $< 0.005$). Component decomposition reveals that `layer_count` and `mean_layer_gap` carry most of the signal; the five-component HII composite never exceeds its best constituent. The mechanism chain extends to coarse-graining identity stability under the current classifier, with `layer_count` predicting switch rate at $r = -0.874$. A Simpson's Paradox in the raw data ($r = +0.336$ before $N$ control) is diagnosed: $N$ is the dominant sign-determining confound. Hard limitations include family-specific HII reversals (KR-like, absolute-layered), classifier-dependent Tier 3, and weak near-wall power at $N \geq 52$. These results provide correlational support — not causal demonstration — for a link between hierarchy depth and entropy suppression in finite causal posets.

**Keywords:** causal posets, hierarchy integration, combinatorial entropy, linear extensions, coarse-graining stability, Simpson's Paradox, causal set theory

---

## 1. Introduction

### 1.1 The "Why Does Lor2D Win?" Question

Two companion papers have established the following finite-size results for discrete causal orders:

- **Prediction B** [1]: Using the same eight-family ensemble and action functional $\mathcal{A} = -\beta H + \gamma I$, the companion study shows that a two-dimensional Lorentzian-like (Lor2D) family becomes competitive against KR-like high-entropy posets at critical coupling $\gamma_c$ remaining $O(1)$ over $N = 10$–$44$. Non-target-anchored observable replacements preserve the competitive window. The key fact needed here: *Lor2D achieves systematically lower $\log H$ than MLR at matched $N$, and this ordering is reproducible under action-score replacement.*

- **Prediction A** [2]: Replacing the target-anchored dimension penalty with a consistency constraint yields four-dimensional Lorentzian dominance in the majority of tested configurations, with a margin that grows with $N$ up to $N = 72$. The key fact needed here: *dimensional selection and entropy ordering are connected, but the structural mechanism is unspecified — which is the gap the present paper addresses.*

These results establish *that* Lorentzian-like structures can dominate, and *which* action terms drive the competition. They do not explain *why* — at the level of intrinsic structural properties — certain families achieve systematically lower entropy at fixed poset size.

This gap matters because without a structural mechanism, the competition results remain purely phenomenological: we know the winning family, but not the property that makes it win.

### 1.2 Qualitative Intuition

A simple argument suggests where to look. Combinatorial entropy $H = \log |\mathcal{L}(P)|$ counts the logarithm of the number of linear extensions — the total orders compatible with the partial order $P$. A poset with deeper causal hierarchy (more temporal layers, larger inter-layer gaps, more long-range ordering constraints) leaves fewer compatible total orders.

The analogy is crystallographic: a solid has lower entropy than a gas at the same energy because its internal structure constrains the accessible microstates. We hypothesise that deeper causal hierarchy plays a structurally analogous role in constraining the linear extensions of a finite poset.

### 1.3 From Intuition to Quantitative Test

Making this intuition precise requires three ingredients:

1. **A quantifiable measure.** We define the Hierarchy Integration Index (HII), a composite z-score combining five structural observables — layer count, mean layer gap, long-edge fraction, adjacent-edge fraction, and reduction-edge density — each computable from the poset's Hasse diagram alone (Section 2).

2. **A multi-tier validation design.** No single analysis cleanly separates scale effects from structural effects, as a Simpson's Paradox (Section 4) will demonstrate. We therefore employ three independent tiers: all-family partial correlations (Tier 1), matched-pair differences (Tier 2), and coarse-graining stability linkage (Tier 3), each addressing different potential confounds.

3. **Extension to coarse-graining stability.** If deeper hierarchy compresses entropy, it should also promote stability under coarse-graining: structures with fewer compatible orderings should be less sensitive to element merging. We test this linkage explicitly in Tier 3.

### 1.4 Key Contributions

The present paper provides a systematic, multi-tier correlational test of the hypothesis that hierarchy depth observables predict entropy variation across poset families at fixed $N$. The main findings are:

1. **Negative correlation at fixed $N$** (Tier 1): $r_\text{partial}(\text{HII}, \log H \mid N) = -0.578$ across 8 families, 320 samples.

2. **Strong matched-pair signal** (Tier 2): $r(\Delta\text{HII}, \Delta\log H) = -0.834$ across 46 Lor2D–MLR pairs ($N = 30$–$56$), stable at $-0.836$ (P10–P90) and $-0.839$ (P0–P100).

3. **Supporting confirmatory evidence via CG identity stability** (Tier 3, classifier-contingent): layer count predicts coarse-graining switch rate with $r = -0.874$ (92 samples, $N = 30$–$56$), under the current nearest-centroid classifier.

4. **Simpson's Paradox resolved**: the naïve positive correlation ($r = +0.336$) is a confound artifact driven solely by $N$.

5. **Component decomposition**: `layer_count` is the single strongest predictor; the five-component HII composite never exceeds its best constituent.

### 1.5 The Three-Prediction Tower

The three predictions form a layered structure:

| Level | Prediction | Claim | Status |
|-------|-----------|-------|--------|
| Base | A | Dimensional selection via bounded $\gamma_c$ | Supported [2] |
| Middle | B | Entropy–action ordering consistency | Supported [1] |
| Top | **C** | **Hierarchy–entropy mechanism** | **Correlational support** |

Each prediction is logically independent: falsifying any one does not invalidate the others. However, their combined significance exceeds the sum. If all three hold — as the current data support — the narrative becomes that existential screening not only selects dimension (A) and maintains self-consistency (B), but also selects a specific hierarchical structure as the mechanism for entropy suppression (C).

C provides the *why* behind A and B's *what*.

### 1.6 Epistemic Positioning

The results presented here are **correlational, not causal**. The three-tier design provides consistent directional support for the chain

$$\text{HII} \uparrow \;\to\; \log H \downarrow \;\to\; \sigma_\text{CG} \downarrow,$$

but strict causality would require either (i) interventional experiments — modifying hierarchy depth while holding all other properties fixed — or (ii) a counting-theoretic proof connecting layer structure to bounds on $|\mathcal{L}(P)|$.

This paper should be read as a structural hypothesis with strong correlational backing, positioned as a bridge between the phenomenological competition results (Predictions A and B) and a future theoretical account.

### 1.6a Hard Limitations (Preview)

To set expectations, we flag four material limitations up front; each is discussed in full in Sections 5 and 6.

1. **Family-specific reversals.** Within KR-like ($\text{layer\_count} \equiv 3$) and absolute-layered ($\text{layer\_count} \equiv \lfloor N/4 \rfloor$) families, HII variance is dominated by secondary components and the correlation can reverse sign (Section 5.5).
2. **Feature redundancy.** `adjacent_edge_fraction` and `long_edge_fraction` are algebraically complementary; `cover_density` aliases `reduction_edge_density`. The effective dimensionality of the five-component HII is at most three.
3. **Classifier-dependent Tier 3.** The coarse-graining identity switch rate ($\sigma_\text{CG}$) depends on a nearest-centroid classifier, not on a physics-derived observable. Tier 3 results should be interpreted as evidence under the current classification protocol, not as a physical robustness law.
4. **Near-wall power loss.** At $N \geq 52$, MLR survivor rates under P10–P90 drop below 0.01%, requiring relaxed filters (P5–P95) with corresponding loss of matching quality.

### 1.7 Paper Structure

The remainder of the paper is organised as follows.

- **Section 2** defines the ensemble, observables, HII formula, and the three-tier validation design.
- **Section 3** presents the main results across all three tiers, including a three-stringency sensitivity analysis.
- **Section 4** diagnoses and resolves the Simpson's Paradox in the naïve Tier 1 analysis.
- **Section 5** decomposes HII into its five components and identifies the primary structural drivers.
- **Section 6** discusses implications, limitations, the near-wall sampling boundary, and future directions.

---

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

---

## 3. Results

### 3.1 Scope and Roadmap

This section presents the quantitative results of the three-tier validation. The Simpson's Paradox diagnosis (Section 4) and the HII component decomposition (Section 5) each receive their own sections because they require extended treatment.

The central question across all three tiers is:

> At fixed poset size $N$, is deeper hierarchy integration associated with lower combinatorial entropy?

We present Tier 1 (all-family partial correlations with exact entropy), Tier 2 (matched-pair $\Delta$-correlations across three filter stringencies), Tier 3 (structural features predicting CG stability), and a cross-tier summary.

---

### 3.2 Tier 1: All-Family Partial Correlation

#### Dataset

Tier 1 uses 320 samples from `outputs_frozen_exact/raw_samples.csv`: 8 families $\times$ 4 sizes ($N = 10, 12, 14, 16$) $\times$ 10 samples per cell. All entropy values are exact (DP over antichains), eliminating approximation error.

#### Primary result

The partial correlation $r(\text{HII}, \log H)$ depends critically on which variables are controlled. Table 5 reports the result under five control sets.

**Table 5.** Partial correlation $r(\text{HII}, \log H)$ under different control variables ($n = 320$).

| Control set | $r_\text{partial}$ | $p$ (perm) | Direction |
|-------------|---------------------|-----------|-----------|
| aw, cf, geo_dim (original Prediction B controls) | $+0.336$ | $0.0005$ | **Positive — spurious** |
| **$N$ only** | $\mathbf{-0.578}$ | $\mathbf{< 0.001}$ | **Negative — true structural direction** |
| $N$ + family dummies | $-0.250$ | $< 0.001$ | Negative (attenuated) |
| Family dummies only | $+0.148$ | $0.043$ | Positive (N confound remains) |
| $N$ + aw + cf + geo_dim | $-0.434$ | $< 0.001$ | Negative |

The sign flip between the first and second rows — from $+0.336$ to $-0.578$ — is a textbook Simpson's Paradox. Adding $N$ as a control is the dominant factor needed to reveal the true structural direction. The diagnosis of this paradox is the subject of Section 4.

#### Fixed-$N$ cross-family correlations

To visualize the within-$N$ relationship without regression residuals, we compute bare Pearson $r(\text{HII}, \log H)$ within each $N$ slice.

**Table 6.** Cross-family Pearson $r(\text{HII}, \log H)$ at each fixed $N$ ($n = 80$ per slice, 8 families $\times$ 10 samples).

| $N$ | $r$ | $p$ | Direction |
|-----|-----|-----|-----------|
| 10 | $-0.86$ | $< 0.001$ | Negative |
| 12 | $-0.73$ | $< 0.001$ | Negative |
| 14 | $-0.74$ | $< 0.001$ | Negative |
| 16 | $-0.53$ | $< 0.001$ | Negative |

The correlation is consistently negative at all four $N$ values. The attenuation from $-0.86$ at $N = 10$ to $-0.53$ at $N = 16$ is consistent with increasing within-family entropy variance at larger $N$, which introduces noise that dilutes the cross-family signal.

#### Interim conclusion

Tier 1 confirms the predicted negative direction once $N$ is controlled. However, the naïve analysis without $N$ control yields the opposite sign. This methodological finding motivates the matched-pair design of Tier 2, which eliminates $N$ as a confound by construction.

---

### 3.3 Tier 2: Matched-Pair $\Delta$-Correlation

#### Dataset

Tier 2 uses 46 matched Lor2D–MLR pairs from the combined matched-pair pool:

| $N$ | Pairs | Filter |
|-----|-------|--------|
| 30 | 12 | P10–P90 |
| 40 | 12 | P10–P90 |
| 44 | 6 | P10–P90 |
| 48 | 4 | P10–P90 |
| 52 | 6 | P5–P95 |
| 56 | 6 | P5–P95 |
| **Total** | **46** | |

Each pair consists of one Lor2D sample and one MLR "survivor" at the *same* $N$. The $\Delta$ for each feature is defined as $\Delta f = f_\text{MLR} - f_\text{Lor2D}$, so a negative $\Delta\text{HII}$ (MLR has lower HII than Lor2D, as expected) co-occurring with a positive $\Delta\log H$ (MLR has higher entropy) would support the prediction.

#### Primary result

The composite HII$_\Delta$ achieves $r(\Delta\text{HII}, \Delta\log H) = -0.834$ ($p < 0.001$). All five individual HII components also show the predicted sign direction; `mean_layer_gap` ($r = -0.836$) and `layer_count` ($r = -0.816$) nearly match the composite. The full component-level breakdown is deferred to Table 10 (Section 5.2); here we note that `reduction_edge_density` is the weakest contributor ($r = +0.459$), consistent with its modest role in the composite.

The component hierarchy — mean_layer_gap $\approx$ HII $\approx$ layer_count $\gg$ long_edge_fraction $\gg$ reduction_edge_density — suggests that one or two components carry most of the signal (Section 5).

#### Sensitivity to filter stringency

A critical robustness test: does the correlation depend on how strictly the MLR survivors are filtered?

**Table 7.** Stability of $r(\Delta\text{HII}, \Delta\log H)$ across three filter windows.

| Filter | Quantile window | $N$ range | Pairs | $r$ |
|--------|----------------|-----------|-------|-----|
| Expanded | P10–P90 | 30–48 | 34 | $-0.836$ |
| **Moderate** | **P5–P95** | **30–56** | **46** | $\mathbf{-0.834}$ |
| Rescue | P0–P100 | 30–56 | 50 | $-0.839$ |

The correlation coefficient varies by less than 0.005 across the three stringency levels. This near-invariance is a strong robustness check: it indicates that the HII–$\log H$ relationship is not an artifact of the particular MLR-filtering protocol. Even under maximal relaxation (P0–P100, where any MLR sample within the Lor2D min–max range is accepted), the correlation is indistinguishable from the stringent version.

#### Interpretation

The matched-pair design eliminates $N$ as a confound by construction: each $\Delta$ is computed within a pair at the same $N$. The strong negative $r$ therefore reflects the within-$N$ structural relationship *without* the Simpson's Paradox contamination that complicates Tier 1.

That individual components (`mean_layer_gap`, `layer_count`) achieve correlations comparable to or exceeding the full HII composite indicates that the five-component formula adds little incremental prediction beyond its two dominant terms. This finding is expanded in Section 5.

---

### 3.4 Tier 3: Coarse-Graining Stability Linkage

#### Dataset

Tier 3 uses 92 samples from the combined duel CSVs (expanded + moderate), encompassing both Lor2D and MLR at $N = 30, 40, 44, 48, 52, 56$. The outcome variable is `cg_family_switch_rate` ($\sigma_\text{CG}$).

#### Primary result

The single strongest predictor of CG stability is `layer_count` ($r = -0.874$, $p < 0.001$), surpassing the composite HII ($r = -0.820$). `mean_layer_gap` ($r = -0.847$) and `long_edge_fraction` ($r = -0.803$) also achieve $|r| > 0.8$. The full component-level analysis, including `adjacent_edge_fraction` and `reduction_edge_density`, is presented in Table 11 (Section 5.3). This result extends the association chain: deeper hierarchy (more layers) not only correlates with lower entropy (Tier 2) but also with greater identity stability under the current classification protocol.

#### Physical interpretation

Posets with more layers impose more global ordering constraints. When elements are randomly removed (coarse-grained), these constraints create a "structural backbone" that preserves family identity. Posets with fewer layers (e.g., KR-like with only 3, or Lor4D with $\sim 2$) lose their identity more easily because their ordering constraints are sparser.

---

### 3.5 Cross-Tier Summary

**Table 8.** Key effect sizes across the three tiers.

| Tier | Design | Best metric | $r$ | $p$ | $N$ range | Samples |
|------|--------|------------|-----|-----|-----------|---------|
| 1 | All-family partial corr | HII vs $\log H$ (controlling $N$) | $-0.578$ | $< 0.001$ | 10–16 | 320 |
| 2 | Matched-pair $\Delta$ | HII$_\Delta$ vs $\Delta\log H$ | $-0.834$ | $< 0.001$ | 30–56 | 46 pairs |
| 3 | Feature $\to$ $\sigma_\text{CG}$ | layer_count vs $\sigma_\text{CG}$ | $-0.874$ | $< 0.001$ | 30–56 | 92 |

All three tiers show consistently large negative effect sizes ($|r| > 0.5$), with Tiers 2 and 3 exceeding $|r| = 0.8$. The direction is uniformly negative once the Simpson's Paradox is resolved.

#### Why three tiers matter

No single tier suffices:

- **Tier 1 alone** has the Simpson's Paradox pitfall: the naïve analysis yields a positive $r$, which could be mistaken for a falsification.
- **Tier 2 alone** covers only Lor2D vs. MLR at $N \ge 30$, missing the small-$N$ exact-entropy domain and six other families.
- **Tier 3 alone** does not directly test the HII–$\log H$ link; it only tests the downstream consequence ($\sigma_\text{CG}$).

Together, the tiers provide mutually reinforcing correlational support for an association chain:

$$\text{deeper hierarchy} \;\longrightarrow\; \text{lower entropy} \;\longrightarrow\; \text{greater CG identity stability.}$$

Each arrow is supported by $|r| > 0.8$ in at least one tier, under the current analysis protocol. The chain is correlational, not causal — establishing strict causality remains an open challenge (Section 6.3).

---

## 4. Simpson's Paradox: Poset Size as the Dominant Sign-Determining Confound

### 4.1 The Anomaly

The most counterintuitive result in Tier 1 is the *sign* of the naïve partial correlation. Controlling for antichain width, comparable fraction, and geometric dimension — the standard controls from Prediction B — gives

$$r_\text{partial}(\text{HII},\;\log H \mid \text{aw}, \text{cf}, \text{geo\_dim}) \;=\; +0.336 \qquad (p = 0.0005).$$

This is positive: higher HII appears to be associated with *higher* entropy, in direct contradiction to the prediction. If accepted at face value, this result would falsify Prediction C. This section shows that the positive sign is a Simpson's Paradox driven predominantly by the poset size $N$ — the dominant confound in Tier 1 — and that the true structural relationship is negative.

### 4.2 Stepwise Confound Identification

To isolate the responsible variable, we systematically vary the control set and observe the sign of the partial correlation.

**Table 9.** Stepwise diagnosis of the Simpson's Paradox ($n = 320$).

| Step | Control variables | $r_\text{partial}$ | Sign | Interpretation |
|------|------------------|---------------------|------|----------------|
| 0 | None (raw Pearson) | $+0.12$ | $+$ | Weak positive; $N$-scaling dominates |
| 1 | aw, cf, geo_dim | $+0.336$ | $+$ | **Simpson's Paradox** — controls absorb wrong variance |
| 2 | **$N$ only** | $\mathbf{-0.578}$ | $\mathbf{-}$ | **Sign flip** — $N$ resolves the paradox |
| 3 | $N$ + family dummies | $-0.250$ | $-$ | Same direction; attenuated by within-family noise |
| 4 | $N$ + aw + cf + geo_dim | $-0.434$ | $-$ | Full controls with $N$: still negative |

The diagnostic pattern is clear: $N$ is the dominant sign-determining confound. Adding $N$ at step 2 produces the largest absolute change in $r$ ($+0.336 \to -0.578$, a swing of $0.914$). No other single control achieves a comparable effect.

### 4.3 Why $N$ Is a Confound

$N$ simultaneously drives both HII and $\log H$, but with different magnitudes:

1. **$N \to \log H$**: The number of linear extensions grows combinatorially with $N$. Across all 8 families, $r(N, \log H) > 0.91$. This is a first-order scaling effect: more elements create exponentially more compatible orderings.

2. **$N \to \text{HII}$**: Larger posets generically admit deeper layer structures, so HII tends to increase with $N$ across families (e.g., Lor2D has $r(N, \text{HII}) \approx +0.65$). This is a second-order structural effect: more elements permit — but do not guarantee — deeper hierarchy.

3. **Cross-$N$ contamination**: Without controlling for $N$, the analysis conflates two effects:

   - A **scaling effect** (positive): larger $N$ $\Rightarrow$ both HII $\uparrow$ and $\log H \uparrow$.
   - A **structural effect** (negative): at fixed $N$, deeper hierarchy $\Rightarrow$ fewer compatible orderings $\Rightarrow$ $\log H \downarrow$.

   Because the scaling effect operates on a larger dynamic range than the structural effect (entropy grows super-linearly in $N$, while HII grows sub-linearly), the aggregate relationship is dominated by scaling and appears positive.

### 4.4 An Ecological Correlation

This is a textbook instance of an *ecological fallacy*: the aggregate-level trend (positive) is the reverse of the within-group trend (negative). The "groups" here are the fixed-$N$ slices.

The fixed-$N$ Pearson correlations (Table 6) confirm this directly: at every $N$ from 10 to 16, the cross-family $r(\text{HII}, \log H)$ is negative ($-0.86$ to $-0.53$). The aggregate positive trend exists *only* because families with higher HII (e.g., Lor2D, Interval Order) also have higher *baseline* entropy at every $N$, and larger $N$ amplifies both quantities.

### 4.5 Physical Analogy: Crystallization at Fixed Temperature

An intuitive analogy clarifies the structure of the paradox:

- **Temperature** (analogue of $N$): A hotter melt has more thermal kinetic energy (higher entropy baseline) and more molecular mobility (higher structural complexity baseline).
- **At fixed temperature**: Crystalline solids have lower entropy than liquids because their internal structure constrains the microstate count.
- **Across temperatures**: Comparing a cold crystal to a hot gas shows the crystal having *lower* entropy with *more* order — but this correlation disguises the fact that temperature (system size) is the true driver.

The Prediction C analogy:
- **$N$** plays the role of temperature.
- **HII** plays the role of crystalline order.
- **$\log H$** plays the role of entropy.

Studying the structure–entropy relationship requires *fixing the temperature* (controlling for $N$).

### 4.6 Verification: Family Controls Alone Are Insufficient

A natural question: can family dummies substitute for $N$ control? The answer is no:

$$r_\text{partial}(\text{HII},\;\log H \mid \text{family dummies only}) \;=\; +0.148 \qquad (p = 0.043).$$

The sign remains positive, because within each family, $N$ still acts as a confound. Family membership absorbs between-family variance but not between-$N$ variance. Only explicit $N$ control (or the matched-pair design of Tier 2) resolves the paradox.

### 4.7 Why the Original Controls Fail

The Prediction B controls (aw, cf, geo_dim) were designed for a different purpose: controlling for structural similarity in the action-score comparison. They fail here because:

1. **Collinearity with $N$**: `antichain_width` scales approximately as $N^{1/2}$; `comparable_fraction` weakly depends on $N$; `geo_dim_eff` has family-dependent $N$-scaling. These controls absorb some $N$-variance but not enough to overcome the dominant $N \to \log H$ pathway.

2. **Wrong regression target**: The Prediction B controls were not chosen to purge $N$-scaling from an entropy residual — they were chosen to equalize structural baselines across families. Re-purposing them for the HII–$\log H$ test introduces an inadequate adjustment that paradoxically amplifies the confound (from $r = +0.12$ to $r = +0.336$).

### 4.8 Methodological Implication

The Simpson's Paradox has a general lesson for cross-size poset studies:

> **Raw correlations between structural observables and entropy-like quantities should not be trusted without explicit $N$ control.** The super-linear growth of $\log H$ with $N$ creates confounds that can reverse structural relationships.

This conclusion extends beyond the present analysis: any study comparing poset families across different sizes must either (i) include $N$ as a regression control, (ii) stratify by $N$, or (iii) use a matched-pair design at identical $N$.

The matched-pair Tier 2 design was motivated precisely by this concern, and its results (Section 3.3) confirm the negative direction without any $N$-related ambiguity.

---

## 5. Component Decomposition and Mechanism Refinement

### 5.1 Motivation

The Hierarchy Integration Index is a composite of five z-scored structural observables:

$$\text{HII} = \frac{z_{\text{lc}} + z_{\text{mlg}} + z_{\text{lef}} - z_{\text{aef}} - z_{\text{red}}}{5}$$

where lc = layer_count, mlg = mean_layer_gap, lef = long_edge_fraction, aef = adjacent_edge_fraction, red = reduction_edge_density. The equal-weight design was specified *a priori*, as part of the confirmatory protocol.

Sections 3 and 4 established that HII correlates strongly and negatively with $\log H$ (Tier 2) and with coarse-graining instability $\sigma_\text{CG}$ (Tier 3) once the confound $N$ is properly controlled. A natural follow-up question is whether all five components contribute comparably, or whether a small subset carries most of the signal. This section decomposes the composite into its constituents across three independent analyses.

**Important caveat.** Within certain families the HII components can be structurally fixed (e.g., KR-like posets always have exactly 3 layers; absolute-layered always have $\lfloor N/4 \rfloor$ layers). The decomposition reported here is therefore a *cross-family, controlled-$N$* analysis, and its conclusions should not be projected to within-family variation without further evidence.

---

### 5.2 Component-Level $\Delta$-Correlations (Tier 2)

The matched-pair design eliminates $N$ as a confound by construction (Section 2.6). For each pair, we compute $\Delta_\text{component} = \text{value}_\text{MLR} - \text{value}_\text{Lor2D}$ alongside $\Delta_{\log H}$. Correlations are shown for the 46-pair dataset (P5–P95 filter stringency).

**Table 10.** Component-level Pearson $r$ with $\Delta_{\log H}$ (46 matched pairs).

| Component | $r(\Delta_\text{comp},\;\Delta_{\log H})$ | $p_\text{perm}$ | Expected sign |
|-----------|------------------------------------------|-----------------|---------------|
| mean_layer_gap | $-0.836$ | $< 0.001$ | $-$ ✓ |
| **HII (composite)** | $\mathbf{-0.834}$ | $< 0.001$ | $-$ ✓ |
| layer_count | $-0.816$ | $< 0.001$ | $-$ ✓ |
| long_edge_fraction | $-0.643$ | $< 0.001$ | $-$ ✓ |
| adjacent_edge_fraction | $+0.643$ | $< 0.001$ | $+$ ✓ |
| reduction_edge_density | $+0.459$ | $0.001$ | $+$ ✓ |

All five components have the predicted sign, and all are significant at the $\alpha = 0.01$ level. The top three — mean_layer_gap, HII, layer_count — are nearly indistinguishable in magnitude ($|r| \in [0.816, 0.836]$), while long_edge_fraction and adjacent_edge_fraction form a second tier ($|r| \approx 0.64$), and reduction_edge_density is the weakest ($|r| = 0.46$).

Note that adjacent_edge_fraction and long_edge_fraction are algebraically complementary (they sum to a near-constant for fixed $N$), so their $|r|$ values are identical by construction.

---

### 5.3 Component-Level CG Stability (Tier 3)

The coarse-graining stability test (Section 2.5) provides an entirely independent channel: the outcome variable is $\sigma_\text{CG}$ (switch rate under coarse-graining), not $\log H$. Component correlations from the 92-sample extended dataset (N = 30–56):

**Table 11.** Component-level Pearson $r$ with $\sigma_\text{CG}$ (92 samples, Tier 3 extended).

| Component | $r(\text{comp},\;\sigma_\text{CG})$ | $p_\text{perm}$ |
|-----------|-------------------------------------|-----------------|
| **layer_count** | $\mathbf{-0.874}$ | $< 0.001$ |
| mean_layer_gap | $-0.847$ | $< 0.001$ |
| **HII (composite)** | $-0.820$ | $< 0.001$ |
| long_edge_fraction | $-0.803$ | $< 0.001$ |
| adjacent_edge_fraction | $+0.803$ | $< 0.001$ |
| reduction_edge_density | $+0.571$ | $< 0.001$ |

The ranking shifts slightly compared to Tier 2: **layer_count now leads**, followed by mean_layer_gap, with HII in third place. The composite does not exceed its best individual component. This pattern is consistent across both filter stringencies and across the original (68-sample) and extended (92-sample) CG datasets.

---

### 5.4 Cross-Tier Comparison

Combining the three independent analyses — fixed-$N$ raw correlation, matched-pair $\Delta$-method, and coarse-graining stability — yields a coherent ranking.

**Table 12.** Component-level $|r|$ across three analyses.

| Component | Fixed $N = 14$ | Tier 2 $\Delta$ (46 pairs) | Tier 3 CG (92 samples) | Rank |
|-----------|---------------|---------------------------|------------------------|------|
| **layer_count** | $0.787$ | $0.816$ | $\mathbf{0.874}$ | **1st** |
| mean_layer_gap | $0.745$ | $\mathbf{0.836}$ | $0.847$ | 2nd |
| long_edge_fraction | $0.771$ | $0.643$ | $0.803$ | 3rd |
| adj_edge_fraction | $0.771$ | $0.643$ | $0.803$ | 3rd (mirror) |
| red_edge_density | $0.693$ | $0.459$ | $0.571$ | 5th |
| **HII (composite)** | $0.649$ | $0.834$ | $0.820$ | — |

Two patterns emerge:

1. **`layer_count` is the single strongest predictor in 2 of 3 analyses** (fixed-$N$ and Tier 3), and a close second in Tier 2. Together with its high correlation partner `mean_layer_gap`, the layer-depth pair accounts for the bulk of the signal.

2. **The five-component HII composite never exceeds its best constituent.** In the fixed-$N$ analysis, HII ($|r| = 0.649$) actually *underperforms* all individual components. In Tiers 2 and 3, HII lies between the top and bottom components. Equal weighting dilutes the strongest signal (layer_count) with weaker ones (reduction_edge_density).

---

### 5.5 Within-Family Behavior at $N = 14$

To isolate cross-family effects from within-family noise, we fix $N = 14$ and compute raw Pearson $r(\text{HII}, \log H)$ within each of the eight families ($n = 10$ per family).

**Table 13.** Within-family $r(\text{HII}, \log H)$ at $N = 14$.

| Family | $r$ | Sign | $p$ | HII variability |
|--------|-----|------|-----|-----------------|
| multi_layer_random | $-0.861$ | $-$ | $0.001$ | High |
| interval_order | $-0.746$ | $-$ | $0.013$ | Moderate |
| transitive_percolation | $-0.678$ | $-$ | $0.031$ | Moderate |
| lorentzian_like_3d | $-0.675$ | $-$ | $0.032$ | Moderate |
| lorentzian_like_4d | $-0.435$ | $-$ | $0.209$ | Low |
| lorentzian_like_2d | $-0.402$ | $-$ | $0.249$ | Moderate |
| absolute_layered | $+0.294$ | $+$ | $0.410$ | **Very low** |
| **KR_like** | $\mathbf{+0.782}$ | $+$ | $0.008$ | **Very low** |

Six of eight families show the predicted negative sign. The two positive outliers — KR-like and absolute-layered — have near-constant layer structure by construction: KR always has exactly 3 layers; absolute-layered always has $\lfloor N/4 \rfloor = 3$ layers at $N = 14$. Their HII variance is almost entirely noise from the secondary components. The positive sign in KR-like ($r = +0.782$, significant) likely reflects a minor counter-directional effect of reduction_edge_density within a structurally frozen hierarchy.

This pattern is *expected* rather than anomalous: when the primary driver (layer_count) is clamped, the composite HII loses its predictive content, and residual fluctuations from weaker components dominate.

---

### 5.6 Physical Interpretation: Layer Depth as the Core Mechanism

The component decomposition converges on a single structural feature: **the number of distinct temporal layers** in the poset.

Physically, each additional layer imposes a global ordering constraint. A poset with $k$ layers has at most

$$|\mathcal{L}(P)| \;\leq\; \prod_{i=1}^{k} |L_i|!$$

linear extensions, where $L_i$ is the $i$-th antichain-layer. For a fixed total $N = \sum_{i} |L_i|$, increasing $k$ (hence decreasing the average $|L_i|$) reduces the product combinatorially. More layers means smaller antichains per layer, hence fewer permutation degrees of freedom, hence lower entropy. The correlation $r(\text{layer\_count}, \log H) < 0$ is therefore not merely empirical but has a transparent combinatorial origin.

The secondary component `mean_layer_gap` is strongly correlated with `layer_count` ($r > 0.9$ across all $N$ values). It captures the average separation between adjacent layers — a geometric re-expression of the same layering depth. Their joint contribution adds no incremental signal beyond either alone (Table 12).

`long_edge_fraction` and `adjacent_edge_fraction` capture *edge-level* consequences of layering: deeper hierarchies tend to have more long edges (spanning non-adjacent layers) and fewer adjacent edges. These are downstream companions, not independent mechanisms.

`reduction_edge_density` is the weakest component, and the only one that is non-significant in the Tier 1 partial-correlation analysis ($|r| < 0.04$, $p > 0.5$). It measures the density of the transitive reduction, which depends on local connectivity rather than global hierarchy. Its inclusion in HII adds noise rather than signal.

---

### 5.7 Implications for Index Refinement

The data suggest that a parsimonious alternative:

$$\text{HII}_\text{narrow} = \frac{z_{\text{layer\_count}} + z_{\text{mean\_layer\_gap}}}{2}$$

would achieve comparable or better predictive power. However, we do *not* adopt this modification, for three reasons:

1. **A priori discipline.** The five-component HII was defined before any data analysis. Post-hoc narrowing, no matter how well-motivated by the current results, risks overfitting to a specific ensemble.

2. **Sufficiency.** The original HII achieves $|r| > 0.82$ across all controlled analyses (Tables 10–11). There is no practical need for refinement at this stage.

3. **Generalisability.** A dedicated study with held-out families or larger $N$ ranges should evaluate whether the dominance of layer_count persists, or whether reduction_edge_density and the edge-fraction pair gain importance in structurally richer regimes.

This observation is flagged as a direction for future work (Section 6).

---

### 5.8 Summary

The five-component HII composite is predictively valid but mechanistically reducible. `layer_count` and its geometric correlate `mean_layer_gap` account for the majority of the entropy–hierarchy relationship in cross-family comparisons. The composite's equal weighting dilutes but does not destroy the signal. The within-family analysis at a single $N$ slice ($N = 14$) provides limited additional support: 6 of 8 families show the predicted negative direction, but the analysis is confined to one small-$N$ point and should not be over-generalised. The two exceptions (KR-like, absolute-layered) are structurally expected: when the primary driver is clamped by construction, the composite reverts to noise.

---

## 6. Discussion

### 6.1 Summary of Findings

This study tested the prediction that deeper hierarchy integration in finite causal posets correlates negatively with combinatorial entropy. Three independent tiers of evidence support this correlation once the confound $N$ (poset size) is properly controlled:

- **Tier 1** (all-family partial correlations, $N = 10$–$16$, $n = 320$): $r_\text{partial}(\text{HII}, \log H \mid N) = -0.578$; consistent negative sign at each of the four $N$ values.
- **Tier 2** (matched Lor2D–MLR pairs, $N = 30$–$56$, 46 pairs at P5–P95): $r(\Delta\text{HII}, \Delta\log H) = -0.834$; stable at $-0.836$ (P10–P90, 34 pairs) and $-0.839$ (P0–P100, 50 pairs).
- **Tier 3** (coarse-graining stability, $N = 30$–$56$, 92 samples): $r(\text{layer\_count}, \sigma_\text{CG}) = -0.874$; $r(\text{HII}, \sigma_\text{CG}) = -0.820$.

A Simpson's Paradox in the naïve Tier 1 analysis — where controlling for structural covariates *without* $N$ produces a positive $r = +0.336$ — was fully diagnosed (Section 4). The component decomposition (Section 5) identified `layer_count` as the primary structural driver.

---

### 6.2 What Is Established

Three claims are supported by the present data.

**First**, a robust *correlational association* between hierarchy depth and combinatorial entropy in cross-family comparisons at fixed $N$. The negative correlation holds across three independent designs, two disjoint $N$ ranges, and three filter stringency levels. No tested confound structure changes the sign once $N$ is controlled, though family-specific reversals (KR-like, absolute-layered) and classifier dependence (Tier 3) constrain the generality of this finding.

**Second**, a *methodological principle*: cross-size poset studies must control for $N$ or use matched-pair designs. Failure to do so can reverse structural relationships (Table 9). This Simpson's Paradox is not a peculiarity of our particular composite; it arises from the systematic $N$-scaling of both HII components and $\log H$.

**Third**, *correlational support for an association chain*:

$$\text{deeper hierarchy} \;\to\; \text{lower entropy} \;\to\; \text{greater CG identity stability}.$$

The first link is supported by Tiers 1 and 2; the second by Tier 3, where higher-HII posets exhibit lower coarse-graining switch rates under the current classifier. The chain is consistent across analyses with $|r| > 0.8$ in all controlled tests, but remains correlational and classifier-dependent (Tier 3).

---

### 6.3 What Is *Not* Established

**Causality.** All three tiers provide correlational evidence only. We do not claim that hierarchy depth *causes* entropy reduction. The relationship may be mediated by unmeasured structural properties, or both HII and $\log H$ may be driven by a common underlying geometric constraint yet to be identified. A causal claim would require an intervention design — e.g., controlled injection of additional layers into a poset while holding all other properties fixed.

**Universality.** The results are limited to the eight families in our ensemble. Other causal set ensembles — particularly those with dynamic growth rules or non-layered topology — may behave differently.

**Continuum limit.** All analyses operate at finite $N \leq 56$. Whether the HII–$\log H$ correlation persists, strengthens, or weakens as $N \to \infty$ is unknown. The near-wall dead zone (§6.5) already signals that the accessible $N$ range for matched-pair tests is severely bounded.

---

### 6.4 Relation to Predictions A and B

Prediction C is the third in a series that examines structural selection in finite causal posets.

**Prediction A** established that a bounded path-length ratio $\gamma_c \sim 4$ selects for low-dimensional ($3\!+\!1$) manifold-like posets from a diverse ensemble. It answers *what* is selected.

**Prediction B** demonstrated that action-score ordering and entropy ordering are non-circularly consistent: replaceable observables reproduce the rank order predicted by the action functional. It answers *that* the selection is self-consistent.

**Prediction C** asks *why* certain families achieve lower entropy. The answer, within the correlational evidence available, is that they possess deeper hierarchy integration — more temporal layers, larger inter-layer gaps, fewer adjacent edges.

#### 6.4.1 Logical Independence

Each prediction is independently testable. Prediction C's HII–$\log H$ correlation can be verified without knowing whether action–entropy orderings are consistent (B) or whether $\gamma_c$ selects $3\!+\!1$ dimensions (A). Falsifying any one prediction does not logically falsify the others.

#### 6.4.2 Semantic Layering

Although logically independent, the three predictions are semantically layered:

- If B fails but C holds, the HII–$\log H$ correlation remains valid but loses its role as "the mechanism behind B's success."
- If B holds but C fails, then B's success arises from structural reasons other than hierarchy depth.
- If both B and C hold — as the current data support — then C *elevates* B from an interesting numerical coincidence to a regularity backed by structural mechanisms.

The three predictions together form a tower:

| Level | Prediction | Claim | Evidential status |
|-------|-----------|-------|-------------------|
| Base | A | Dimensional selection via $\gamma_c$ | Supported [2] |
| Middle | B | Entropy–action ordering consistency | Supported [1] (non-circular) |
| Top | C | Hierarchy–entropy correlation | **Correlational support** |

The top level is explicitly weaker than the bottom two: A and B make sharp, falsifiable predictions; C establishes a correlational pattern and identifies a mechanism *candidate*.

---

### 6.5 The Near-Wall Boundary

The matched-pair design (Tier 2) requires MLR "survivors" — multi_layer_random posets that pass structural similarity filters against the Lor2D reference window. At larger $N$, survivors become vanishingly rare.

**Table 14.** Near-wall survival rates across three filter stringencies.

| $N$ | P10–P90 | P5–P95 | P0–P100 |
|-----|---------|--------|---------|
| 48 | 4 / 3,468 (0.115%) | — | — |
| 52 | 1 / 60,000 (0.002%) | 6 / 12,601 (0.048%) | 8 / 1,978 (0.40%) |
| 56 | 0 / 80,000 (0%) | 6 / 45,121 (0.013%) | 8 / 333 (2.4%) |

The P10–P90 window is effectively dead at $N \geq 52$. The P5–P95 moderate window rescues workable samples, but at exponentially increasing computational cost. Nevertheless, the Tier 2 correlation is stable across all three levels: $r = -0.836$ (P10–P90), $-0.834$ (P5–P95), $-0.839$ (P0–P100). The less-than-0.005 variation demonstrates that the HII–$\log H$ relationship is insensitive to the matching tightness.

This boundary imposes a practical ceiling on the matched-pair approach. Extending Tier 2 beyond $N \approx 60$ will require either targeted generative models (constrained random posets matching specific Lor2D structural profiles) or a fundamentally different matching strategy.

---

### 6.6 Limitations

1. **Exact entropy is computationally bounded.** The DP algorithm for $\log H$ is exponential in antichain width. Tier 1 operates at $N \leq 16$; Tiers 2 and 3 use approximate methods. This limits both the $N$ range and the precision of entropy estimates.

2. **HII is a bespoke composite.** The five-component, equal-weight formula was defined for this study. While the a priori specification prevents data snooping, the component decomposition (Section 5) shows that the weighting is suboptimal. Future work should investigate data-driven or information-theoretic weighting — on held-out samples.

3. **Within-family HII variance is structurally suppressed.** KR-like and absolute-layered posets have near-constant layer structure by construction. Their within-family HII fluctuations are noise from secondary components. The cross-family analysis is therefore essential, but it introduces the $N$-confound challenge diagnosed in Section 4.

4. **No causal identification strategy.** The three-tier design provides convergent correlational evidence but does not support causal inference. An intervention design — e.g., surgically adding or removing elements from specific layers — would be needed to test the mechanism chain directly.

5. **Two-family matching only.** Tier 2 uses Lor2D–MLR pairs exclusively. A more comprehensive test would match across all ${8 \choose 2} = 28$ family pairs, but most combinations lack sufficient structural overlap for meaningful matching at large $N$.

6. **Matching quality degrades with $N$.** Even within the Lor2D–MLR pair, balance diagnostics (standardized mean differences on non-matched covariates) worsen as the acceptance window widens from P10–P90 to P5–P95. At $N = 52$–$56$, the moderate-filter MLR survivors are structurally less similar to Lor2D than at $N = 30$–$40$. The correlation stability across stringency levels (Table 7) mitigates this concern, but does not eliminate it: the $\Delta$-correlation captures family-mean differences, not individual-pair precision.

7. **Language boundaries respected throughout.** All findings are stated in correlational terms. Readers should not interpret the association chain (§6.2) as a demonstrated causal pathway.

---

### 6.7 Future Directions

1. **HII refinement on held-out data.** Generate new families or use alternative poset growth rules to test whether a two-component $\text{HII}_\text{narrow}$ (layer_count + mean_layer_gap) matches or exceeds the full five-component index.

2. **Near-wall generative models.** Replace exhaustive rejection sampling with constrained generative algorithms that produce MLR-like posets at $N > 60$, enabling Tier 2 extension.

3. **Causal intervention design.** Construct posets by inserting or removing layers while preserving $N$ and edge density, testing whether the entropy change matches the HII–$\log H$ correlation magnitude.

4. **Polynomial-time entropy estimation.** Develop approximate counting algorithms for linear extensions in the $N = 50$–$200$ range, enabling Tier 1-style all-family analysis at larger scales.

5. **Analytic derivation.** Seek a proof linking layer count to bounds on $|\mathcal{L}(P)|$ for specific poset families, providing a theoretical underpinning for the empirical correlation.

6. **Joint A–B–C modelling.** Integrate $\gamma_c$, action scores, and HII in a single statistical framework to test for mediation effects across the three-prediction tower.

---

## 7. Conclusion

We have shown that at fixed poset size $N$, hierarchy depth observables — principally `layer_count` and `mean_layer_gap` — correlate negatively with combinatorial entropy across three independent statistical designs spanning $N = 10$ to $56$. The pre-registered five-component HII composite achieves $|r| > 0.8$ in Tiers 2 and 3, but never exceeds its best constituent; the signal is carried by the depth pair. The correlation is stable across three levels of MLR-filter stringency ($r$ variation $< 0.005$). A Simpson's Paradox in the naïve analysis reveals that cross-size comparisons are unreliable without $N$ controls — a methodological finding of independent value. Hard limitations remain: family-specific reversals (KR-like, absolute-layered), classifier-dependent Tier 3, and weak near-wall power at $N \geq 52$. Within these boundaries, the results provide correlational support for a link between hierarchy depth and entropy suppression in finite causal posets, offering a structural account of *why* certain families achieve lower entropy as observed in Prediction B.

---

## References

[1] Zhang, G. Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy. *Entropy* **2026**, submitted. Code and data: https://github.com/unicome37/poset_phase (DOI: 10.5281/zenodo.18980657).

[2] Zhang, G. Dimensional Selection Without Dimensional Priors: 4D Lorentzian Dominance in Finite Causal Poset Ensembles Under Consistency-Based Actions. *Entropy* **2026**, submitted. Code and data: https://github.com/unicome37/poset_phase (DOI: 10.5281/zenodo.18980657).

[3] Surya, S. Evidence for the continuum in 2D causal set quantum gravity. *Class. Quantum Grav.* **2012**, *29*, 132001.

[4] Brightwell, G.; Winkler, P. Counting linear extensions. *Order* **1991**, *8*, 225–242.

[5] Pratt, J.W.; Gibbons, J.D. *Concepts of Nonparametric Theory*; Springer: New York, NY, USA, 1981.