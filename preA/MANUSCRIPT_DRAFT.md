# Why 3+1? Dimensional Selection via Structural Consistency in Causal Set Geometry

**Author**: Gang Zhang (Independent Researcher, unicome@gmail.com)

---

## Abstract

A companion paper established that a bounded geometric phase transition exists in finite causal posets under an action combining entropy and geometric penalties, and that the target-anchored dimension proxy (d = 2) can be replaced by a non-target-anchored consistency penalty without destroying the transition. Here we test a direct prediction of that replacement: if the action rewards dimensional *consistency* rather than a specific dimension, then higher-dimensional Lorentzian-like posets should dominate. We extend the ensemble to include a 4D Lorentzian-like family and compare four families (2D, 3D, 4D Lorentzian-like, and KR) across N = 20–72 under six action variants. Under the original target-anchored action, 4D posets and KR posets are locked in a near-perfect tie (36 vs 35 wins out of 98 configurations). Systematic ablation localizes the suppression to a single term (g_dim). Replacing it with the non-target-anchored g_con yields unconditional 4D dominance (92/98), with a margin that grows monotonically with N. Seed sensitivity at N = 72 confirms 100% win rate (21/21) under the consistency action, versus 43% under the original. These results provide the first finite-size evidence that dimensional consistency constraints, without dimension-specific priors, naturally favor higher-dimensional Lorentzian structures in discrete causal order ensembles.

---

## 1. Introduction

The question "why does spacetime have 3+1 dimensions?" is among the oldest in fundamental physics. In most approaches to quantum gravity, the dimensionality of spacetime is either postulated (as in loop quantum gravity [6]), compactified from a higher-dimensional starting point (as in string theory), or left implicit in the measure over geometries. Carlip [5] has argued that dimensional reduction at short distances is generic across quantum gravity approaches, but the question of why the *large-scale* dimension is 3+1 remains open. Causal set theory (CST) [2,3] offers a distinctive setting for this question: spacetime emerges as a statistical phase from a competition between combinatorial entropy and structural order in ensembles of locally finite partial orders (posets). In this framework, the dimensionality of the emergent geometry is, in principle, *selected* by the dynamics rather than imposed — making the transition from *a priori* dimensionality to *emergent* dimensionality concrete and numerically testable.

In a companion paper [1], we took a first step: we showed that under an action A₂ combining combinatorial entropy H = log|L(P)| with geometric penalty terms, the critical coupling γ_c at which Lorentzian-like 2D posets overtake entropically dominant Kleitman–Rothschild (KR) structures [4] remains O(1) across N = 10–44. This established a *bounded geometric phase transition* — Lorentzian geometry can compete as a phase. A key further finding was that the dimension-specific proxy g_dim (penalizing deviation from d = 2) can be replaced by a non-target-anchored consistency penalty g_con (penalizing local–global dimension *inconsistency* without specifying a target), with the phase transition preserved.

This replacement raises a natural second-order question: if the action rewards dimensional *consistency* rather than proximity to a specific d, then among Lorentzian-like structures of different dimensionality, which dimension does the action prefer? This is the "why 3+1?" question recast in a finite, discrete, numerically tractable setting.

We test this question directly. We extend the poset ensemble to include a 4D Lorentzian-like family and conduct systematic experiments across N = 20–72, comparing 2D, 3D, 4D Lorentzian-like, and KR posets under the original and consistency-based actions. Our approach follows a "diagnose → ablate → replace → confirm" sequence:

1. Under the original target-anchored action, the 4D family is suppressed by a single identifiable term (g_dim), resulting in a closely contested lor4d–KR tie.
2. Systematic ablation identifies g_dim as the dominant source of this low-dimensional bias.
3. Replacing g_dim with the non-target-anchored g_con yields unconditional 4D dominance with growing margin.
4. Seed sensitivity confirms robustness at N = 72.

This establishes that the dimension preference is not an intrinsic property of the geometric penalties as a whole, but a specific consequence of whether the dimension term is target-anchored or consistency-based — a nontrivial step toward turning "why 3+1?" from a metaphysical presupposition into a quantitatively testable structural question.

---

## 2. Ensemble and Methods

### 2.1 Extended Poset Ensemble

We retain three families from the companion paper and add one new family:

- **Lorentzian-like 2D** (`lor2d`): Layered posets mimicking 2D causal diamond structure, with layer sizes following a diamond profile and inter-layer edges determined by spatial proximity.
- **Lorentzian-like 3D** (`lor3d`): Three-dimensional analogue with broader layers and lower comparable fraction.
- **KR-like** (`kr`): Three-layer bipartite posets with random inter-layer edges, approximating the KR structure [4].
- **Lorentzian-like 4D** (`lor4d`): *New*. Four-dimensional analogue with layer sizes following a 4D causal diamond profile.

The 4D generator follows the same construction principle as the 2D and 3D generators: elements are distributed across layers whose sizes follow a d-dimensional diamond profile, and inter-layer edges are determined by spatial proximity in d−1 spatial dimensions. The key structural consequence of higher dimensionality: `lor4d` has (i) higher entropy H (broader layers, more linear extensions), (ii) higher g_dim penalty (since d_eff ≈ 4 ≠ 2), and (iii) *lower* g_con penalty (strong internal dimensional coherence).

For each family and each N ∈ {20, 24, 28, …, 72} (14 sizes), we generate K = 4 independent samples. Linear extension counts are computed exactly for `lor2d` (sub-second up to N = 72) and via SIS (Sequential Importance Sampling, 4096 runs) for other families at N > 32.

### 2.2 Action Variants

All actions share the form S(P) = −βH(P) + γ·I(P), with β = 1 and penalty I decomposed into neutral and geometric components. We test six variants:

| Variant | Dim term | Other geo | Purpose |
|---------|----------|-----------|---------|
| A₂ full | g_dim (target d=2) | all | Original baseline |
| A₂ no-dim | removed | all | Ablation |
| A₂ interval-only | removed | interval only | Minimal geometry |
| A₂ no-wh | g_dim | no g_wh | Ablation |
| A₂ consistency | g_con (w=0.6829) | all | Core replacement |
| A₂ multi-con. | g_mcon | all | Robustness check |

Consistency weight w_con = 0.6829 is calibrated to match the scale of g_dim across the original 7-family ensemble [1]. For each (N, γ, variant) configuration (γ ∈ {0.0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0}), we compute the mean action score per family and identify the **winner** (lowest score). The **margin of victory** is the score gap to the runner-up.

---

## 3. Ablation: Locating the Low-Dimensional Bias

### 3.1 The Negative Result

Under the original A₂ full action, the 4D family does *not* achieve stable dominance. Across all 98 (N, γ) configurations, the winner counts are:

| Lor4D | KR | Lor2D | Lor3D |
|-------|----|-------|-------|
| 36 | 35 | 19 | 8 |

The near-perfect 36:35 tie between `lor4d` and KR is itself informative: the target-anchored penalty creates a dimension-specific "glass ceiling" that prevents `lor4d` from expressing its structural advantages (higher entropy and lower consistency penalty).

### 3.2 Locating the Source of Bias

To identify which component suppresses `lor4d`, we perform systematic ablation:

1. **Remove g_dim only** (A₂ no-dim): `lor4d` win count rises sharply. The suppression is localized to this single term.
2. **Retain only interval terms** (A₂ interval-only): No stable dimension preference emerges. Interval geometry alone is insufficient.
3. **Remove g_wh** (A₂ no-wh): `lor4d` remains suppressed. The width-height term is not the cause.

This mirrors the companion paper's finding: in both the 2D-vs-KR phase transition (Prediction B) and the cross-dimension competition (Prediction A), g_dim is the single most consequential term.

---

## 4. Non-Target-Anchored Replacement

We replace g_dim with g_con [1], keeping all other penalty terms unchanged:

> g_con(P) = Var[d_local] + (d̄_local − d_global)²

where d_local is the order-fraction dimension proxy applied to sampled causal intervals with |I| ≥ 4, and d_global is the same estimator applied to the full poset.

| Variant | Lor4D | KR | Lor2D | Lor3D |
|---------|-------|----|-------|-------|
| A₂ full (original) | 36 | 35 | 19 | 8 |
| A₂ consistency (g_con) | 92 | — | — | 6 |
| A₂ multi-con. (g_mcon) | 89 | — | — | 9 |

Under the consistency replacement, `lor4d` wins 92 out of 98 configurations. KR and `lor2d` *never* win. The only non-`lor4d` wins go to `lor3d` at a handful of small-N, low-γ configurations. The multi-estimator variant (g_mcon) yields a nearly identical pattern (89/98), confirming estimator independence.

Against `lor2d` specifically, `lor4d` wins all 98 pairwise comparisons (100%). Against `lor3d`, it wins 92/98 (the 6 losses occur at N = 20–28, γ = 0.0, where entropy differences are small).

---

## 5. Finite-Size Scaling

### 5.1 Growing Margin of Victory

The margin of `lor4d` over both competitors grows monotonically with N under A₂ consistency:

| N | vs Lor2D (mean) | vs Lor2D (min) | vs Lor3D (mean) | vs Lor3D (min) |
|---|-----------------|----------------|-----------------|----------------|
| 20 | 7.41 | 1.70 | 2.27 | −0.17 |
| 28 | 14.21 | 7.77 | 2.59 | −2.07 |
| 36 | 19.30 | 11.41 | 4.20 | −1.72 |
| 44 | 27.12 | 19.42 | 6.02 | 0.36 |
| 52 | 36.22 | 30.40 | 11.29 | 7.27 |
| 60 | 45.41 | 38.88 | 14.50 | 10.02 |
| 68 | 50.82 | 41.51 | 16.41 | 9.44 |
| 72 | 56.72 | 47.15 | 16.53 | 10.06 |

Against `lor2d`, the mean margin grows nearly linearly from +7 (N = 20) to +57 (N = 72), and even the minimum margin (worst-case γ) is positive at every N. Against `lor3d`, the margin crosses into uniformly positive territory by N ≥ 44.

### 5.2 Seed Sensitivity at N = 72

We repeat the N = 72 experiment with three independent generator seed offsets (0, 50000, 100000):

| Variant | Seed 0 | Seed 50k | Seed 100k |
|---------|--------|----------|-----------|
| A₂ full | KR 4, Lor4d 3 | KR 4, Lor4d 3 | KR 4, Lor4d 3 |
| A₂ consistency | Lor4d 7/7 | Lor4d 7/7 | Lor4d 7/7 |
| A₂ multi-con. | Lor4d 7/7 | Lor4d 7/7 | Lor4d 7/7 |

The KR–`lor4d` split under A₂ full is identical (4:3) at every seed — a structural property of the action, not a sampling artifact. Under both consistency variants, `lor4d` achieves 100% win rate (21/21).

| Variant | vs | Mean margin | Min margin | Win rate |
|---------|----|-------------|------------|----------|
| A₂ full | Lor2d | 28.51 | −17.56 | 81.0% |
| A₂ full | Lor3d | 0.78 | −27.59 | 57.1% |
| A₂ consistency | Lor2d | 56.00 | 47.15 | 100% |
| A₂ consistency | Lor3d | 17.28 | 10.06 | 100% |
| A₂ multi-con. | Lor2d | 55.82 | 46.74 | 100% |
| A₂ multi-con. | Lor3d | 17.20 | 9.91 | 100% |

---

## 6. Discussion and Limitations

### 6.1 Interpretation

The results confirm the prediction: the consistency-based action naturally selects higher-dimensional geometric structures. The logic is straightforward. The 4D family has higher entropy (competitive at low γ), low consistency penalty (competitive at high γ), but high target-anchored penalty under g_dim (suppressed under the original action). The consistency replacement removes this artificial suppression while retaining the dimensional coherence reward.

Both Prediction A and the companion paper's Prediction B employ the same methodology — diagnosis, ablation, replacement, robustness check — and arrive at the same structural conclusion: the mechanism depends on dimensional *consistency*, not dimensional *identity*.

### 6.2 Relation to Prior Work

Carlip [5] argued that dimensional reduction is generic in quantum gravity approaches. Surya [6] and Loomis and Carlip [7] studied causal set dynamics numerically. Our work provides a complementary perspective: rather than asking whether spacetime dimensionality *flows*, we ask whether a specific dimension can be *selected* by consistency constraints in a discrete statistical ensemble. The finite-size evidence presented here suggests that such selection is possible, at least among the tested families.

### 6.3 Limitations

1. **Finite ensemble**: Only four families are tested; the claim is "4D wins among tested types," not global optimality.
2. **SIS estimation**: For N > 32 in KR, `lor3d`, and `lor4d`, entropy is estimated via SIS (4096 runs). Systematic biases cannot be fully excluded.
3. **Generator dependence**: The 4D generator's structural properties determine its penalty values; different 4D constructions might yield different results.
4. **Physical dimension selection**: "Higher dimension wins" does not directly explain why our universe has d = 3+1. Additional constraints would be needed.
5. **Consistency weight**: w_con = 0.6829 is calibrated on the companion paper's 7-family ensemble.

### 6.4 Outlook

Three directions are natural: (i) extending to 5D, 6D, and beyond to test whether dominance continues or saturates; (ii) developing faster exact algorithms to eliminate SIS uncertainty; (iii) incorporating dynamical or field-theoretic constraints to investigate whether a preferred dimension emerges.

---

## References

[1] G. Zhang, "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy," *Entropy* (2026), submitted.

[2] L. Bombelli, J. Lee, D. Meyer, R. D. Sorkin, "Space-time as a causal set," Phys. Rev. Lett. **59**, 521 (1987).

[3] R. D. Sorkin, "Causal sets: Discrete gravity," in *Lectures on Quantum Gravity* (Springer, 2005). arXiv:gr-qc/0309009.

[4] D. J. Kleitman and B. L. Rothschild, "Asymptotic enumeration of partial orders on a finite set," Trans. Amer. Math. Soc. **205**, 205–220 (1975).

[5] S. Carlip, "Dimension and dimensional reduction in quantum gravity," Class. Quantum Grav. **34**, 193001 (2017).

[6] S. Surya, "The causal set approach to quantum gravity," Living Rev. Relativ. **22**, 5 (2019).

[7] S. P. Loomis and S. Carlip, "Suppression of non-manifold-like sets in the causal set path integral," Class. Quantum Grav. **35**, 024002 (2018).

[8] G. Brightwell and P. Winkler, "Counting linear extensions," Order **8**, 225–242 (1991).

[9] W. J. Cunningham and S. Surya, "Dimensionally restricted causal set quantum gravity," Class. Quantum Grav. **37**, 054002 (2020).

---

## Acknowledgments

**AI Assistance Statement**: This work was conducted in collaboration with large language model assistants (Claude, Anthropic). The AI contributed to experimental design, Python code implementation, data analysis, theoretical interpretation, and manuscript drafting. All scientific decisions, theoretical motivations, and final interpretations were made by the human author, who takes full responsibility for the content.

## Code and Data Availability

All code, configuration files, and output data are available at: https://github.com/unicome37/poset_phase

Archived version with DOI: https://doi.org/10.5281/zenodo.18980657
