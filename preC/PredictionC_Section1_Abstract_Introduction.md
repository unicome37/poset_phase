# Paper Skeleton: Abstract and Section 1 Introduction

## Working Title Options

### Option A

**Hierarchy Integration Lowers Combinatorial Entropy: A Three-Tier Correlational Study in Finite Causal Posets**

### Option B

**Why Lorentzian-like Posets Win: Structural Hierarchy as an Entropy Compression Mechanism in Discrete Causal Orders**

### Option C

**Deeper Causal Layers, Fewer Linear Extensions: Correlational Evidence for Structural Entropy Suppression**

Recommended current choice:

> **Hierarchy Integration Lowers Combinatorial Entropy: A Three-Tier Correlational Study in Finite Causal Posets**

This title is precise, avoids causal overclaim ("lowers" states the correlation direction, not a proven causal mechanism), and makes the three-tier methodology immediately visible.

---

## Abstract (~200 words)

**Structure**: Context → Gap → Method → Key Results → Implication

- **Context**: Companion papers [Prediction B, A] established that Lorentzian-like posets can become competitive phases in action-weighted causal poset ensembles. What remains unexplained is the structural *mechanism* by which certain families achieve lower combinatorial entropy at fixed poset size.

- **Gap**: The structural feature that connects causal architecture to entropy compression has not been identified or systematically tested.

- **Method**: We define a Hierarchy Integration Index (HII), a composite z-score measure of layer depth, layer gap structure, and edge composition. We test the correlation between HII and combinatorial entropy $\log H$ across three independent analysis tiers: (i) all-family exact computation (8 families, N = 10–16, 320 samples), (ii) matched-pair $\Delta$ analysis (Lor2D vs MLR, 46 pairs, N = 30–56), and (iii) coarse-graining stability linkage (68 samples).

- **Key Results**: At fixed N, higher HII correlates strongly with lower $\log H$ (partial $r = -0.578$ controlling for N; matched-pair $r = -0.834$, stable across three filter stringencies). The mechanism chain extends to coarse-graining stability: layer count alone predicts CG identity switch rate with $r = -0.888$. A Simpson's Paradox in the raw data (positive $r$ before N control) is fully resolved by identifying N as a confound.

- **Implication**: The results provide correlational support for a specific mechanism chain — deeper causal hierarchy → fewer linear extensions → greater coarse-graining stability — that may underlie the finite-size competitive advantage of Lorentzian-like posets established in companion work. The composite HII shows within-family direction reversal, motivating a narrower focus on layer count and mean layer gap as primary predictors.

---

## Suggested Keywords

- causal posets
- hierarchy integration
- combinatorial entropy
- linear extensions
- coarse-graining stability
- Simpson's Paradox
- structural mechanism
- causal set theory

---

## 1. Introduction

### 1.1 The "Why Does Lor2D Win?" Question

Two companion papers have established the following finite-size results:

- **Prediction B** [ref]: A 2D Lorentzian-like family becomes competitive against KR-like high-entropy posets under an action $\mathcal{A} = -\beta H + \gamma I$, with the critical coupling $\gamma_c$ remaining $O(1)$ over $N = 10$–$44$. Non-target-anchored replacements for dimension-specific terms preserve the competitive window.

- **Prediction A** [ref]: Replacing the target-anchored dimension penalty with a consistency constraint yields unconditional 4D Lorentzian dominance (92/98 configurations), with growing margin up to $N = 72$.

These results establish *that* Lorentzian-like structures can become competitive, and *which* action terms drive the competition. They do not yet explain *why* — at the level of intrinsic structural properties — certain poset families achieve systematically lower entropy at fixed size.

This gap matters because without a structural mechanism, the competition results remain purely phenomenological: we know the winning family, but not the structural property that makes it win.

### 1.2 Qualitative Intuition

A simple qualitative argument suggests where to look. Combinatorial entropy $H = \log |L(P)|$ counts the logarithm of the number of linear extensions — the total orders compatible with the partial order. A poset with deeper causal hierarchy (more layers, larger inter-layer gaps, more long-range causal connections) imposes more ordering constraints, leaving fewer compatible total orders.

By analogy: in statistical physics, a crystalline solid has lower entropy than a gas at the same energy because its internal structure constrains the microstate count. The hypothesis of the present paper is that deeper causal hierarchy plays a structurally analogous role in constraining the linear extensions of a finite poset.

### 1.3 From Intuition to Quantitative Test

Making this intuition precise requires:

1. **A quantifiable measure of "depth of hierarchy"**: We define the Hierarchy Integration Index (HII), a composite z-score combining five structural observables.

2. **A multi-tier validation design**: Because no single analysis can cleanly separate scale effects from structural effects (as we will show through a Simpson's Paradox), we employ three independent tiers spanning different $N$ ranges and analysis designs.

3. **Extension to coarse-graining stability**: If deeper hierarchy compresses entropy, it should also promote stability under coarse-graining — structures with fewer compatible orderings should be less likely to change identity when elements are merged. We test this linkage explicitly.

### 1.4 Key Contribution

The present paper provides the first systematic, multi-tier correlational test of the hypothesis that causal hierarchy integration is the structural mechanism driving entropy compression in finite posets. The main findings are:

1. **Negative correlation at fixed N** (Tier 1): partial $r(\text{HII}, \log H \mid N) = -0.578$ across 8 families.

2. **Strong matched-pair signal** (Tier 2): $\Delta\text{HII}$ vs $\Delta\log H$ gives $r = -0.834$ across 46 Lor2D–MLR pairs (N = 30–56), stable at $r = -0.839$ under rescue filtering (50 pairs).

3. **Extended mechanism chain** (Tier 3): layer count predicts CG switch rate with $r = -0.888$.

4. **Simpson's Paradox resolved**: The raw positive correlation (+0.336 before N control) is a confound artifact; N is the sole driver.

5. **HII narrowing**: Within-family, the composite HII shows direction reversal. The primary within-family predictors are `layer_count` and `mean_layer_gap`, not the full five-component composite.

### 1.5 The Three-Prediction Tower

The three predictions form a layered structure:

- **A** (base): Existential selection picks 3+1 dimensions via bounded $\gamma_c$.
- **B** (middle): Self-consistency — action and entropy agree on family rankings.
- **C** (present work): Structural mechanism — hierarchy integration explains *why* certain families achieve lower entropy.

Each prediction is logically independent and can be falsified without affecting the others. However, their combined semantic weight is greater than the sum: if all three hold, the narrative becomes that *existence progressively narrows structural possibilities through layered constraints*. C provides the "why" behind A and B's "what."

### 1.6 Epistemic Positioning

We emphasize from the outset that the results presented here are **correlational**, not causal. The three-tier design provides consistent directional support for the mechanism chain HII $\uparrow$ → $\log H \downarrow$ → $\sigma_\text{CG} \downarrow$, but establishing strict causality would require either (i) interventional experiments (modifying HII while holding other properties constant) or (ii) a counting-theoretic proof connecting hierarchy depth to bounds on $|L(P)|$.

The present work should therefore be read as a structural hypothesis with strong correlational backing, positioned as a bridge between the phenomenological competition results (Predictions A and B) and a future theoretical account of why certain causal architectures suppress combinatorial entropy.

### 1.7 Paper Structure

- Section 2: Ensemble, observables, HII definition, and three-tier validation design.
- Section 3: Main results across all three tiers.
- Section 4: Simpson's Paradox — diagnosis, resolution, and physical interpretation.
- Section 5: HII component decomposition and the case for narrowing to core predictors.
- Section 6: Discussion, limitations, near-wall sampling boundary, and outlook.

---

## Writing Notes

- Mirror Prediction B's tone: technical, explicitly careful about epistemic scope.
- Avoid "causal chain" or "proves mechanism"; use "correlational support for mechanism chain" consistently.
- Emphasize the Simpson's Paradox finding as a methodological contribution — it demonstrates that naïve all-family correlations can be misleading without proper N control.
- The three-tier design should be framed as a strength: each tier addresses a different potential confound.
- Note the near-wall boundary early (Introduction §1.5 or Methods §2) so it's not a surprise in Discussion.
