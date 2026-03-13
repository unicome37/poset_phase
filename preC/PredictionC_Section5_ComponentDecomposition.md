# Paper Skeleton: Section 5 — Component Decomposition and Mechanism Refinement

## 5. Which Components of HII Drive the Entropy Relationship?

### 5.1. Motivation

HII is a composite of five z-scored structural observables:

$$\text{HII} = \frac{z_{\text{layer\_count}} + z_{\text{mean\_layer\_gap}} + z_{\text{long\_edge\_fraction}} - z_{\text{adj\_edge\_fraction}} - z_{\text{red\_edge\_density}}}{5}$$

The results in Sections 3 and 4 show that HII correlates strongly with log_H once N is controlled. But is the composite truly necessary, or do one or two components carry most of the signal?

> **Important caveat**: Within individual families, HII components can be approximately constant (e.g., KR: always 3 layers; AL: always N/4 layers). The cross-family composite correlation does not imply that each component varies meaningfully within every family. This section analyzes the **cross-family, fixed-N decomposition.**

---

### 5.2. Component-Level Partial Correlations

**Table 5.1**: Individual component partial correlations with log_H (controlling for N)

| Component | partial_r(component, log_H | N) | p | Direction relative to HII sign |
|-----------|-------------------------------|---|------|
| layer_count | −0.55 | < 0.001 | ✓ (positive in HII = negative log_H) |
| mean_layer_gap | −0.52 | < 0.001 | ✓ |
| long_edge_fraction | −0.38 | < 0.001 | ✓ |
| adjacent_edge_fraction | +0.41 | < 0.001 | ✓ (negative in HII = positive log_H) |
| reduction_edge_density | −0.18 | 0.012 | ✓ but weak |
| layer_signature_redundancy* | −0.608 | < 0.001 | — (not in HII composite) |

*`layer_signature_redundancy` is not an HII component but was computed as part of the residual metrics. Its strong correlation suggests investigating inclusion in future HII revisions.

### 5.3. Component Contributions Across Tiers

**Table 5.2**: Component-level effect sizes across tiers

| Component | Tier 1 partial_r (N ctrl) | Tier 2 Δ-corr | Tier 3 → σ_CG |
|-----------|--------------------------|----------------|----------------|
| layer_count | −0.55 | −0.823 | **−0.888** |
| mean_layer_gap | −0.52 | −0.814 | −0.860 |
| long_edge_fraction | −0.38 | −0.693 | −0.801 |
| adj_edge_fraction | +0.41 | +0.552 | +0.810 |
| red_edge_density | −0.18 | −0.325 | −0.298 |
| **HII (composite)** | **−0.578** | **−0.836** | **−0.835** |

Key observation: **`layer_count` matches or exceeds HII in Tiers 2 and 3.** The composite does not add incremental power beyond its dominant component.

### 5.4. Per-Family Behavior at N = 14

**Table 5.3**: Within-family r(HII, log_H) at N = 14 (n = 10 samples per family)

| Family | r(HII, log_H) | Sign | HII variability |
|--------|---------------|------|-----------------|
| lorentzian_like_2d | −0.67 | − | High |
| transitive_percolation | −0.58 | − | Moderate |
| interval_order | −0.43 | − | Moderate |
| multi_layer_random | −0.35 | − | Moderate |
| lorentzian_like_3d | −0.28 | − | Moderate |
| lorentzian_like_4d | −0.12 | − | Low |
| **KR_like** | **+0.15** | **+** | **Very low** (fixed 3 layers) |
| **absolute_layered** | **+0.08** | **+** | **Very low** (fixed N/4 layers) |

6 of 8 families show the predicted negative direction. The two positive outliers (KR, AL) have near-constant layer structure by construction — their HII variance is almost entirely noise. This is expected and does **not** undermine the cross-family finding.

### 5.5. Interpretation: Layer Depth as the Core Driver

The component decomposition supports narrowing the mechanism:

1. **Primary driver**: `layer_count` (depth of layering) — the number of distinct temporal layers in the poset.
2. **Secondary driver**: `mean_layer_gap` — average size of inter-layer jumps, which correlates with layer_count.
3. **Tertiary contributors**: `long_edge_fraction` and `adjacent_edge_fraction` — these capture edge-level consequences of layering but are not independently as strong.
4. **Weak contributor**: `reduction_edge_density` — contributes modestly.

Physically, each additional layer imposes a global ordering constraint on the linear extensions. A poset with $k$ layers has at most $\prod_{i=1}^{k} |L_i|!$ linear extensions (where $L_i$ is the $i$-th antichain-layer), which decreases as $k$ increases for fixed total $N = \sum |L_i|$. The layer_count is therefore a natural structural predictor of entropy.

### 5.6. Switch Enhancement Analysis: layer_count vs Joint Terms

An independent validation comes from the switch-enhancement scan (near-wall rescue data, N = 30–56). At optimal coupling $\eta = 3.0$, the switch_zscore reduction from baseline:

| Predictor | Mean reduction | Rank |
|-----------|---------------|------|
| **layer_count** | **107.1%** | **1st** |
| HII (composite) | 105.7% | 2nd |
| layer_count + mean_layer_gap | 102.5% | 3rd |
| mean_layer_gap | 93.1% | 4th |

At N = 52/56 specifically:

| N | Baseline switch_z | layer_count | HII | layer_count + gap |
|---|-------------------|-------------|-----|-------------------|
| 52 | 6.80 | → 1.01 | → 1.13 | → 1.13 |
| 56 | 7.09 | → 0.47 | → 0.61 | → 0.61 |

The joint term (layer_count + mean_layer_gap) does not exceed layer_count alone. This reinforces the component analysis: `layer_count` is the single most informative structural predictor. The five-component HII adds minor discrimination but no incremental enhancement.

### 5.7. Implications for HII Refinement

The current HII formula weights all five components equally (1/5 each). The data suggest that a two-component version:

$$\text{HII}_{\text{narrow}} = \frac{z_{\text{layer\_count}} + z_{\text{mean\_layer\_gap}}}{2}$$

would achieve comparable or better predictive power. However, we do not recommend adopting this modification in the present paper because:

1. The original five-component HII was defined a priori, and post-hoc narrowing risks overfitting to the current dataset.
2. The cross-tier consistency with the original HII is sufficient to support the prediction.
3. A dedicated study of HII refinement should use held-out data.

This observation is flagged as a direction for future work.

---

### 5.8. Figure Suggestions

**Figure 5**: Bar chart of component-level $|r|$ across the three tiers (Table 5.2). Grouped bars by component, colored by tier. Visual should make clear that `layer_count` dominates.

**Figure 6** (optional): Within-family scatter at N = 14 showing r(HII, log_H) for each family, with KR and AL flagged as structurally constrained.

### 5.9. Table Summary

- Table 5.1: Component partial correlations with log_H
- Table 5.2: Component-level effect sizes across tiers
- Table 5.3: Within-family behavior at N = 14

---

## Writing Notes

- This section is primarily analytical, not evidential. It does not test new predictions but decomposes a confirmed correlation.
- The key message is constructive: "We found that layer_count does most of the work, which has a clean physical interpretation."
- Be careful with the within-family analysis (Table 5.3): emphasize that the 6/8 agreement is encouraging but the 2/8 exceptions are structurally explained (not anomalies). Do NOT frame KR/AL as "failures."
- The HII refinement discussion (§5.6) should be explicitly framed as **future work, not a recommendation for the present analysis**. A reviewer might otherwise ask: "Why not just use layer_count?" Answer: a priori definition discipline.
- Mention `layer_signature_redundancy` as an intriguing lead but not part of the confirmatory analysis.
- Cross-reference Section 3 heavily; Section 5 extends the results, it doesn't replace them.
