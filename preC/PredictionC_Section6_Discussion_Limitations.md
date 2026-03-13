# Paper Skeleton: Section 6 — Discussion, Limitations, and Outlook

## 6. Discussion

### 6.1. Summary of Findings

This study tested the prediction that deeper hierarchy integration (measured by HII) correlates negatively with combinatorial entropy (log_H) in finite causal posets. Three independent tiers of evidence consistently support this correlation once poset size N is controlled:

- **Tier 1** (all-family, N = 10–16): $r_{\text{partial}} = -0.578$ after N control; consistent negative sign at all four N values.
- **Tier 2** (matched lor2d–mlr pairs, N = 30–56): $r(\Delta\text{HII}, \Delta\log H) = -0.834$ (46 pairs, P5–P95); stable at $-0.839$ under P0–P100 rescue (50 pairs).
- **Tier 3** (CG stability, N = 30–56): $r(\text{layer\_count}, \sigma_{\text{CG}}) = -0.874$ (92 samples).

A Simpson's Paradox in the naïve Tier 1 analysis was fully diagnosed: N is the sole confound. The component decomposition and switch-enhancement analysis both identified `layer_count` as the primary structural driver, with the full five-component HII adding no incremental predictive power.

### 6.2. What Is Established

1. **A robust correlational association** between hierarchy integration and entropy across multiple designs, N ranges, and family subsets.
2. **A methodological lesson**: cross-size poset studies must control for N or use matched-pair designs. Failure to do so can reverse structural relationships.
3. **A mechanism chain receiving correlational support**: deeper hierarchy → lower entropy → greater CG stability. This chain is consistent across all three tiers with $|r| > 0.5$.

### 6.3. What Is NOT Established

1. **Causality**. All three tiers provide correlational evidence. We do not claim that hierarchy integration *causes* entropy reduction. The relationship may be mediated by unmeasured structural properties, or both HII and log_H may be driven by a common underlying constraint.

2. **Universality beyond the 8-family ensemble**. The findings are limited to the specific families studied (absolute_layered, KR_like, lorentzian_like_2d/3d/4d, transitive_percolation, interval_order, multi_layer_random). Other causal set ensembles may behave differently.

3. **Continuum limit relevance**. All results are at finite N (≤ 48). Whether the HII–log_H correlation persists, strengthens, or weakens as $N \to \infty$ remains unknown.

### 6.4. Relation to Predictions A and B

**Prediction A** established that bounded path-length ratio $\gamma_c \sim 4$ selects for low-dimensional (3+1) manifold-like posets. Prediction C extends this by asking *why* certain families achieve lower γ_c: deeper hierarchy integration is a structural correlate.

**Prediction B** demonstrated non-circular consistency between action-score ordering and entropy ordering (replaceable by new observables). Prediction C provides a mechanism-level account: the entropy ordering observed in Prediction B correlates with the hierarchy structure measured by HII.

Together, the three predictions form a tower:
- **A** (base): Dimensional selection via γ_c  
- **B** (middle): Entropy–action consistency  
- **C** (top): Hierarchy–entropy mechanism (correlational support)

#### 6.4.1. Logical Independence vs. Semantic Layering

Each prediction is independently testable. Prediction C does not depend on the truth of A or B for its correlational claims — the HII–log_H negative correlation can be verified without knowing whether action–entropy orderings are consistent (B) or whether γ_c selects 3+1 dimensions (A).

However, the **semantic significance** of C depends on B:
- If B fails (action–entropy ordering is inconsistent), C's HII–log_H correlation remains valid, but loses its role as "the mechanism explaining B's success."
- If B holds and C fails, then B's success arises from structural reasons other than hierarchy depth.
- If both B and C hold (the current data supports this), then C elevates B from an "interesting numerical coincidence" to a "regularity supported by structural mechanisms."

#### 6.4.2. Philosophical Implication

When the three-prediction tower is complete, the narrative becomes:

> Existential screening not only selects dimension (A) and maintains self-consistency (B), but also selects a specific hierarchical structure as the mechanism for achieving low entropy (C). This suggests that "existence" is not accidental, but is progressively narrowed through layers of structural constraint.

This is the deeper philosophical significance that motivates the present study: C provides the *why* behind A and B's *what*.

### 6.5. The Near-Wall Boundary

The matched-pair design (Tier 2) depends on the existence of MLR "survivors" — multi_layer_random posets that pass structural similarity filters against the Lor2D reference window. At larger N, these survivors become extremely rare under the standard P10–P90 quantile window:

| N | Window | MLR attempts | MLR accepted | Rate |
|---|--------|-------------|--------------|------|
| 48 | P10–P90 | 3,468 | 4 | 0.115% |
| **52** | P10–P90 | **60,000** | **1** | **0.002%** |
| **56** | P10–P90 | **80,000** | **0** | **0%** |

This "near-wall dead zone" has been addressed with three strategies:

1. **Rescue (P0–P100)**: Expanding to the full min–max range of the Lor2D reference. Recovery rates: N=52 → 0.4% (8/1978), N=56 → 2.4% (8/333). Data available but structural similarity constraint is eliminated.

2. **Moderate (P5–P95)**: A middle-ground configuration preserving partial structural similarity while widening window enough for survival. Currently running for N=52/56. If successful, these samples provide stronger evidence than rescue because the MLR survivors are structurally more Lor2D-like.

3. **Targeted generative models**: A future design that would produce constrained random posets matching specific Lor2D structural profiles (e.g., width, comparable_fraction, layer_count) without exhaustive rejection sampling.

The progressive widening sequence (P10–P90 → P5–P95 → P0–P100) also serves as a **sensitivity analysis**: if the HII–log_H correlation persists across all three filter stringencies, it is robust to the matching quality.

### 6.6. Limitations

1. **Exact entropy computation is exponential.** Log_H was computed exactly (DP over antichains) only up to N = 16 (Tier 1). Tier 2 uses an approximate method. This constrains the accessible N range.

2. **HII is a bespoke composite.** The five-component HII was defined for this study. While the a priori definition prevents data snooping, the choice of components and equal weighting may not be optimal. Future work should investigate information-theoretic component weighting.

3. **Within-family HII variance is low for structured families.** KR_like and absolute_layered have approximately constant layer structure, making within-family HII variation uninformative. The cross-family analysis is therefore essential but introduces the N-confound challenge addressed in Section 4.

4. **No causal identification strategy.** Observational correlation in a computational ensemble does not support causal claims without an intervention design (e.g., controlled injection of additional layers into a poset while holding other properties fixed).

5. **Language boundaries.** The mechanism chain "deeper hierarchy → lower entropy → greater CG stability" is supported by correlational evidence only. We deliberately avoid causal language throughout.

### 6.7. Future Directions

1. **HII refinement study** (held-out data): Test a narrowed two-component HII (layer_count + mean_layer_gap) on independently generated samples.

2. **Near-wall exploration**: Develop generative models capable of producing MLR-like posets at N > 48 to extend Tier 2.

3. **Causal intervention design**: Construct posets by inserting/removing layers while preserving N and edge density, testing whether the entropy change is consistent with the HII–log_H correlation.

4. **Approximate entropy methods**: Develop polynomial-time estimators of log_H for larger N, enabling Tier 1-style analysis at N > 16.

5. **Theoretical derivation**: Attempt an analytic proof linking layer count to linear extension count for specific poset families (e.g., layered posets with uniform antichain sizes).

6. **Integration with Pred A/B**: Jointly model γ_c, action scores, and HII in a single statistical framework to test for mediation effects.

---

## 7. Conclusion

[Writing note: The conclusion should be 1 paragraph, maximum 150 words. Restate the main finding in one sentence, acknowledge the correlational limitation, and point to the most promising future direction.]

Draft:

> We have shown that hierarchy integration, as measured by the five-component HII score, correlates negatively with combinatorial entropy across three independent statistical designs spanning N = 10 to 56. The correlation is large ($|r| > 0.5$ in all tiers, $|r| > 0.8$ in Tiers 2 and 3) and stable across three levels of MLR-filter stringency (P10–P90, P5–P95, P0–P100). A Simpson's Paradox in the naïve analysis reveals that cross-size comparisons are unreliable without N controls — a methodological finding of independent value. Component decomposition and switch-enhancement analysis both identify layer_count as the single strongest predictor. These results provide the strongest correlational support to date for a mechanism chain linking deep causal hierarchy to low combinatorial entropy to coarse-graining stability in finite posets.

---

## Suggested Figure

**Figure (final)**: "Three-prediction tower" schematic showing how Predictions A, B, and C relate. Pred A (base): dimensional selection. Pred B: consistency. Pred C: mechanism correlate. Arrows labeled with what each tier establishes and the remaining open questions.

---

## Writing Notes

- The discussion should be balaned: neither overclaiming ("we proved the mechanism") nor underselling ("these are just correlations").
- Frame the Simpson's Paradox as a contribution, not a mistake.
- The conclusion is deliberately short — all detail is in the body.
- Be explicit about "correlational, not causal" in at least three places: Summary (§6.1), What Is NOT Established (§6.3), and Limitations (§6.6).
- The relation to Predictions A and B (§6.4) is important for readers who have followed the series. It should be self-contained enough for new readers.
- The near-wall boundary (§6.5) is notable because it may limit future replication attempts. Flagging it proactively is better than having reviewers discover it.
- Future directions should be specific and actionable, not generic. Each should have a clear next step.
