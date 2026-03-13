# Prediction C — Manuscript Draft: Section 6

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
