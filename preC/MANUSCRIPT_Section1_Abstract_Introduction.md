# Prediction C — Manuscript Draft: Section 1

## Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study

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
