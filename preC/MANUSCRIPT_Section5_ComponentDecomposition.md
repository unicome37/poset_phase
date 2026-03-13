# Prediction C — Manuscript Draft: Section 5

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
