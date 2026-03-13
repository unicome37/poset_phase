# Hierarchical Depth, Entropy Reduction, and the Mechanism Layer Behind Lorentzian-Like Competition in Finite Causal Posets

## Abstract

We study a mechanism-level question that sits downstream of finite Lorentzian phase competition: among already matched quasi-geometric competitors, does deeper causal hierarchy systematically reduce combinatorial entropy? The question is motivated by a finite causal-poset program in which Lorentzian-like structures emerge as genuine competitors against high-entropy alternatives, but where the strongest remaining rivals inside the filtered window are not generic random layered posets, but highly Lorentzianized near-miss structures.

Prediction C is formulated as a fixed-\(N\) structural claim. We test whether hierarchical depth observables such as `layer_count`, `mean_layer_gap`, and a broader `hierarchy_integration` index are negatively associated with \(\log H(G)\), where \(H(G)\) is the number of linear extensions of a finite causal poset \(G\). The central evidence comes from matched-pair comparisons between window-surviving `multi_layer_random` competitors and nearest-neighbor `lorentzian_like_2d` samples, supplemented by pooled and stratified regressions, cross-family robustness checks, and a second-filter enhancement analysis.

The main result is again two-sided. First, the pairwise signal is strong and stable: in the original matched window, the depth-entropic correlation is strongly negative; after extending the near-wall range through a moderate `P5~P95` window and comparing against both stricter and looser matching windows, the main correlation remains essentially unchanged. Second, the strongest global claim must be narrower than the strongest pairwise claim. Pooled analyses display a direct Simpson's paradox driven mainly by scale \(N\), and once `n + family` fixed effects are introduced, the broader `hierarchy_integration` composite weakens sharply while the direct depth observables `layer_count` and `mean_layer_gap` remain the most stable negative predictors.

These results do not derive Einstein dynamics, and they do not show that hierarchical depth alone replaces the coarse-graining stability filter. They do, however, support a more precise mechanism statement: at fixed finite size, deeper hierarchical organization is associated with lower combinatorial entropy inside the quasi-geometric competition window, and explicit depth terms can reduce the burden carried by the second filter. Prediction C should therefore be read as a mechanism layer beneath Prediction B: logically independent of Lorentzian phase competition, but theoretically most important when both results hold together.

## 1. Introduction

The broader structural program asks whether Lorentzian-like geometry can emerge as a competitive phase inside a finite ensemble of causal partial orders, rather than being inserted directly into the microscopic generating rule. That problem was framed in earlier work as a competition between combinatorial entropy and structural regularization. The present paper asks a narrower, downstream question:

> once a Lorentzian-like family has already survived finite competition, what deeper structural feature helps explain its remaining entropy advantage over quasi-geometric competitors?

This is the question addressed by Prediction C. It is not a claim that Einstein equations have already been derived from discrete causal structure. It is a more limited mechanism test. The working hypothesis is that the relevant difference between a genuine Lorentzian-like structure and its strongest near-miss competitors is not exhausted by surface geometric similarity. Instead, the decisive difference may lie in causal hierarchy itself: deeper, more integrated layered organization may reduce the number of compatible global orderings, and therefore lower combinatorial entropy.

This hypothesis arises from a specific empirical situation. In the filtered competition window defined in the preceding Prediction B work, the strongest remaining non-Lorentzian competitors are not generic `multi_layer_random` samples. They are rare survivors that already look strongly Lorentzian under width, comparability, and coarse interval geometry. Yet even after this filtering, they retain systematically higher \(\log H\) than matched `lorentzian_like_2d` samples. The natural next question is therefore not "why does geometry matter?" in the abstract, but:

> among already matched quasi-geometric competitors, what residual structural variable still tracks the entropy gap?

Prediction C addresses this question through three linked levels of analysis:

1. matched-pair tests inside the quasi-geometric competition window;
2. pooled and stratified regressions over a larger reconstructed dataset;
3. second-filter enhancement scans asking whether explicit depth terms reduce the weight needed on `switch_zscore`.

This paper should be read as the mechanism companion to Prediction B. The logical relation is subtle and important. Prediction C does **not** require Prediction B to be true in order to be testable. One can ask whether hierarchical depth correlates with lower \(\log H\) without already knowing whether Lorentzian-like structures win the action competition. But the theoretical meaning of Prediction C is additive rather than isolated. If Prediction B fails, then Prediction C may still describe a real structural regularity, but it no longer explains Lorentzian-like competitiveness. If Prediction B succeeds while Prediction C fails, then the Lorentzian advantage must be attributed to some other structural cause. Only when both hold together does Prediction C supply a mechanism layer beneath Prediction B.

In this sense, the present paper occupies the third position in a three-prediction tower:

1. Prediction A asks which effective dimension becomes favored;
2. Prediction B asks whether Lorentzian-like geometry becomes a competitive phase at all;
3. Prediction C asks what structural mechanism may help explain that competitiveness once the phase competition window is already in place.

The strongest current claim of Prediction C must be formulated carefully. Pairwise evidence is strong, but global pooled evidence is not uniform. In particular, pooled regressions display a strong Simpson's paradox driven mainly by finite size \(N\). Thus the main claim is not that every broad hierarchy composite universally lowers entropy across the full ensemble. The better claim is narrower:

> at fixed \(N\), direct depth observables, especially `layer_count` and `mean_layer_gap`, are negatively associated with \(\log H\) inside the quasi-geometric competition regime.

The rest of the paper is organized as follows. Section 2 defines the finite ensemble, the matched quasi-geometric competitors, and the depth observables. Section 3 introduces the pairwise design and explains why fixed-\(N\) comparisons are central. Section 4 presents the main pairwise results, including the three-window near-wall sensitivity analysis. Section 5 discusses cross-family robustness, pooled regression, and the Simpson problem. Section 6 studies the relation between depth terms and the second filter. Section 7 discusses current interpretation, limits, and the relation to Predictions B and A. Section 8 concludes.

## 2. Windowed Quasi-Geometric Competitors and Depth Observables

The microscopic objects are finite causal posets \(G=(V,\prec)\) with \(|V|=N\). The entropy observable is

\[
H(G) = \log N_{\mathrm{ext}}(G),
\]

where \(N_{\mathrm{ext}}(G)\) is the number of linear extensions of \(G\). As in the broader program, \(\log H\) is treated as a combinatorial freedom proxy.

Prediction C focuses on a narrower subset of the broader ensemble than Predictions A and B. The central comparison is not between generic family means, but between:

- window-surviving `multi_layer_random` samples,
- nearest-neighbor `lorentzian_like_2d` samples matched on quasi-geometric covariates.

The most important matching covariates are:

- `antichain_width`
- `comparable_fraction`
- `geo_dim_eff`
- `geo_interval_shape`

This restriction matters. The purpose of Prediction C is not to rediscover the obvious fact that radically different families differ in entropy. It is to identify what still explains entropy differences **after** the strongest quasi-geometric competitors have already been brought close together in observable structure.

The main depth-related observables are:

- `layer_count`
- `mean_layer_gap`
- `adjacent_edge_fraction`
- `long_edge_fraction`
- `reduction_edge_density`
- `cover_density`
- `layer_signature_redundancy`
- composite `hierarchy_integration_index`

At the current stage, the most important distinction is between direct depth measures and broader composites. The direct depth measures are:

- `layer_count`, the number of hierarchical layers;
- `mean_layer_gap`, the typical vertical separation between connected layers.

The broader `hierarchy_integration_index` is useful as a descriptive composite, but later sections show that it is less stable than the direct depth measures once fixed effects are introduced.

Conceptually, the working interpretation is straightforward. A structure can look quasi-Lorentzian in width, comparability, and interval geometry while still organizing its causal links in a shallower and more local way. If shallower organization leaves more admissible global orderings unconstrained, then deeper hierarchy should correspond to lower \(\log H\).

Prediction C therefore asks whether depth observables continue to matter **after** the main geometric screen has already been passed.

## 3. Pairwise Design and Why Fixed-\(N\) Comparisons Matter

Prediction C is most naturally formulated as a fixed-\(N\) question. This is not a cosmetic choice. The pooled dataset contains a strong scale effect:

\[
N \uparrow \quad \Longrightarrow \quad \log H \uparrow
\]

and many depth-related observables also scale upward with \(N\). As a result, pooled correlations can be pushed in the wrong direction even when within-\(N\) structure is negative. The pairwise design is meant to neutralize this problem as directly as possible.

The basic matched-pair object is a quasi-geometric difference:

\[
\Delta X = X_{\mathrm{MLR\ survivor}} - X_{\mathrm{Lor2D\ match}},
\]

with the corresponding entropy difference

\[
\Delta \log H = \log H_{\mathrm{MLR\ survivor}} - \log H_{\mathrm{Lor2D\ match}}.
\]

If shallower hierarchy is one reason why the surviving quasi-geometric `MLR` competitors remain entropically expensive, then one expects:

- negative correlations for `layer_count_delta` and `mean_layer_gap_delta`,
- a negative correlation for the broader `hierarchy_integration_delta_index`,
- weaker or secondary roles for more local edge-ratio observables.

This pairwise design should not be confused with a full causal identification strategy. It is still finite, observational, and model-dependent. But it is much closer to the actual mechanism question than broad family-level pooling, because it holds constant the main geometric filters that already define the competition window.

The near-wall extension makes this design even more important. At larger \(N\), the standard strict matching window becomes sample-starved. This creates a real methodological tradeoff: one can loosen the window and gain samples, but only at the cost of weakening structural similarity constraints. Section 4 therefore treats near-wall matching strictness itself as a sensitivity dimension rather than hiding it as a preprocessing detail.

## 4. Main Pairwise Results and Near-Wall Sensitivity

The strongest current evidence for Prediction C comes from matched-pair analysis inside the quasi-geometric competition window.

### 4.1 Expanded window result up to \(N=48\)

The first stabilized pairwise result used 34 matched pairs over \(N=30\) to \(48\). In this sample:

- `hierarchy_integration_delta_index -> log_H_delta` gives \(r \approx -0.836\);
- `layer_count_delta -> log_H_delta` gives \(r \approx -0.823\);
- `mean_layer_gap_delta -> log_H_delta` gives \(r \approx -0.814\).

The permutation result for the composite depth index is approximately \(p \approx 0.0002\). The key point is not merely that the signs are negative, but that the strongest individual signals come from the simplest direct depth measures. This already suggests that the decisive mechanism is better described as hierarchical depth than as a diffuse bundle of all layer-related observables.

### 4.2 Near-wall extension: three matching strictness levels

The original strict near-wall continuation faced severe sample depletion. Under the `P10~P90` window, usable `MLR` survivors at \(N=52\) and \(N=56\) were nearly exhausted. This led to a three-level sensitivity design:

- `P10~P90`: strict expanded window;
- `P5~P95`: moderate near-wall window;
- `P0~P100`: rescue window.

The resulting summary is:

| Window | Pair count | N range | `HII_delta r` | `layer_count_delta r` |
|---|---:|---|---:|---:|
| `P10~P90` | 34 | 30-48 | -0.836 | -0.823 |
| `P5~P95` | 46 | 30-56 | -0.834 | -0.816 |
| `P0~P100` | 50 | 30-56 | -0.839 | -0.827 |

This is one of the strongest current robustness results in the Prediction C line. The main pairwise correlation varies by less than 0.005 across all three matching strictness levels. In other words:

> the main depth-entropic signal is effectively insensitive to whether the near-wall window is strict, moderate, or rescue-level loose.

This result matters methodologically. The rescue window is still too loose to be treated as the ideal main specification. But it is no longer carrying the burden alone. The moderate `P5~P95` window shows that one can rescue usable near-wall sample size **without** materially changing the main pairwise conclusion.

### 4.3 What the pairwise result currently supports

The current pairwise evidence supports a narrower statement than some earlier formulations suggested. It does **not** yet justify the claim that every broad hierarchy composite universally lowers entropy. It does support the following:

> among already matched quasi-geometric competitors, deeper hierarchical organization, especially greater `layer_count` and larger `mean_layer_gap`, is associated with lower \(\log H\).

This is currently the most defensible core statement of Prediction C.

## 5. Cross-Family Directional Robustness, Pooled Regression, and the Simpson Problem

The pairwise result is strong, but it is not the whole story. To assess external robustness, a larger reconstructed cache was used to generate 544 unique posets and 396 cross-family matched pairs. These compare `KR_like`, `Lor3D`, and `Lor4D` against `Lor2D` after matching on key geometric covariates.

The resulting directional correlations are very strong:

- `KR vs Lor2D`: `HII_delta -> log_H_delta` about `-0.929`;
- `Lor3D vs Lor2D`: about `-0.773`;
- `Lor4D vs Lor2D`: about `-0.863`.

The corresponding direct depth measures show similar signs. As a directional statement, this is useful. It shows that "deeper -> lower \(\log H\)" is not confined to the narrow `MLR vs Lor2D` window.

But this robustness check has to be interpreted carefully. It is not a substitute for the windowed mechanism analysis, for three reasons:

1. the cross-family comparisons are structurally much broader than the quasi-geometric `MLR vs Lor2D` competition;
2. several groups, especially `KR_like`, are partially degenerate for some depth variables;
3. matching quality is not strong enough to justify a clean causal reading.

The pooled regression results make the scale problem unavoidable. On the full reconstructed dataset, controlling only for width, comparability, and geometric dimension proxies pushes all major depth observables **positive**:

- `hierarchy_integration_index`: about `+0.879`
- `layer_count`: about `+0.855`
- `mean_layer_gap`: about `+0.878`

This is a direct Simpson's paradox. The dominant source is scale \(N\). Once `n + family` fixed effects are introduced, the picture changes sharply:

- `layer_count`: about `-0.450`
- `mean_layer_gap`: about `-0.424`
- `hierarchy_integration_index`: about `+0.086`

This split is one of the most important current interpretation constraints. It shows that the global claim must be narrower than the pairwise claim. The stronger and more stable fixed-effect predictors are the direct depth measures, not the broader composite.

### 5.1 Current best global statement

The strongest current global statement is therefore:

> at fixed \(N\), direct hierarchical depth observables, especially `layer_count` and `mean_layer_gap`, remain negatively associated with \(\log H\).

This is why the pairwise results are so valuable. They absorb the strongest scale confound by construction.

### 5.2 Additional family-level and methodological limits

Several limits must remain explicit:

- `KR_like` has `layer_count = 3` identically, so some within-family depth interpretation collapses there.
- In `Lor4D`, `layer_count` is only marginal while `mean_layer_gap` remains stronger.
- Around \(N=52\) and \(N=56\), statistical power is still weak even after near-wall rescue.
- The broader composite `hierarchy_integration_index` is close to a significance knife-edge under fixed effects and should not be oversold.

In short, the direction of the signal survives, but the effect size and evidential status depend strongly on whether one is looking at matched within-\(N\) competition or broad pooled comparisons.

## 6. Depth Terms and the Second Filter

The Prediction B line identified a second filter inside the quasi-geometric competition window. In its current reduced form, this filter is best represented by `switch_zscore`, a standardized measure of coarse-grained family identity persistence. Prediction C asks whether hierarchical depth simply duplicates that filter, or whether it explains part of what the filter is doing.

To test this, the enhancement scan extends the effective comparison rule from

\[
\Delta \mathrm{score}_{A2} + \zeta \,\mathrm{switch\_zscore}
\]

to

\[
\Delta \mathrm{score}_{A2} + \zeta \,\mathrm{switch\_zscore} + \eta \,\mathrm{depth\ term}.
\]

The question is whether adding a depth term lowers the crossing weight \(\zeta_{\mathrm{cross}}\) required to flip the pairwise ranking.

At the moderate working point \(\eta = 1.0\), the answer is yes:

- `layer_count` lowers mean `zeta_cross` by about `38.3%`;
- `hierarchy_integration` lowers it by about `38.2%`;
- `mean_layer_gap` lowers it by about `34.6%`.

The `layer_count` path is the cleanest. Over the earlier four-scale window:

- \(N=30\): `2.36 -> 1.28`
- \(N=40\): `3.54 -> 2.21`
- \(N=44\): `4.53 -> 2.92`
- \(N=48\): `6.29 -> 4.14`

Earlier near-wall rescue scans also showed that a joint `layer_count + mean_layer_gap` term can work, but does not beat `layer_count` alone. This reinforces the same lesson seen in the fixed-effect regressions: the direct depth observables are currently stronger and cleaner than wider composites.

The correct interpretation is conservative:

> depth does not replace the second filter, but it does reduce part of the burden that the second filter would otherwise have to carry alone.

This is exactly the kind of result Prediction C was meant to test. It takes the earlier phase-competition picture one step deeper. The second filter is no longer just an empirical stabilizer. It now looks at least partly like a proxy for underlying hierarchical depth.

## 7. Current Interpretation, Limits, and Relation to Predictions B and A

Prediction C is strongest when written as a mechanism statement, not as a grand unification claim.

The current evidence supports three linked conclusions.

First, within matched quasi-geometric competition, deeper hierarchy is associated with lower combinatorial entropy. This is the cleanest result and the least dependent on broad modeling choices.

Second, the main signal is robust to near-wall matching strictness. The three-window comparison shows that the main pairwise result is not a fragile artifact of one specific quantile cutoff.

Third, the strongest global claim must be narrower than the strongest pairwise claim. Once scale and family are controlled, the best-supported variables are `layer_count` and `mean_layer_gap`, not the broader `hierarchy_integration` composite.

These conclusions also clarify the relation to Prediction B. The relation is:

- logically independent,
- semantically additive.

Prediction C can be true even if Prediction B eventually weakens. But its most important role is to explain part of **why** Lorentzian-like structures may be competitive once the phase window is already in place.

At the same time, several limits remain nontrivial.

1. The evidence is finite-size throughout. No asymptotic theorem has been established.
2. Some broader matched datasets suffer from poor covariate balance, so large pooled correlations should not be read naively.
3. `KR_like` degeneracy and near-wall sample sparsity continue to restrict family-level interpretation.
4. The second-filter enhancement remains exploratory; it does not yet yield a unique reduced law.
5. Nothing in the present paper derives Einstein equations. The present result is a mechanism precursor, not a continuum field-equation derivation.

The right interpretation is therefore intermediate:

> Prediction C has moved beyond a vague intuition, but it remains a mechanism-level finite-size result rather than a completed theory of gravitational emergence.

## 8. Conclusion

This paper asked whether the residual entropy gap between matched quasi-geometric competitors is linked to hierarchical depth. The answer is now provisionally yes.

The strongest current evidence is pairwise and fixed-\(N\): among matched quasi-geometric competitors, deeper hierarchy, especially more layers and larger mean layer gaps, is associated with lower \(\log H\). This conclusion survives extension from the original \(N \le 48\) window into the near-wall \(N \le 56\) regime, and it remains essentially unchanged across strict, moderate, and rescue matching windows.

At the same time, the present work also narrows its own claim. Broad pooled analyses are distorted by Simpson's paradox, and once `n + family` fixed effects are introduced, the composite hierarchy index weakens sharply while direct depth measures remain the most stable predictors. The proper current statement is therefore not "all hierarchy composites lower entropy everywhere," but:

> at fixed finite size, direct hierarchical depth observables are the most stable currently known mechanism-level predictors of lower combinatorial entropy inside the quasi-geometric competition regime.

Prediction C should therefore be read as a third-layer result beneath the earlier Lorentzian phase competition line. It does not prove geometric emergence in full. It does show that the Lorentzian-like advantage is increasingly difficult to describe as a surface geometric accident alone. Some of that advantage appears to be carried by deeper hierarchical organization.

## 9. Suggested Figures and Tables

The current unified draft is best supported by the following tables and figures:

1. A table summarizing the three near-wall matching windows:
   - pair count,
   - \(N\) range,
   - `HII_delta r`,
   - `layer_count_delta r`.
2. A pairwise scatter plot of `layer_count_delta` versus `log_H_delta`.
3. A paired bar or line figure showing how `zeta_cross` changes when a `layer_count` term is added.
4. A pooled-versus-fixed-effect comparison table showing the sign flip for:
   - `hierarchy_integration_index`,
   - `layer_count`,
   - `mean_layer_gap`.
5. A schematic "three-prediction tower" figure:
   - Prediction A: dimension selection,
   - Prediction B: Lorentzian phase competition,
   - Prediction C: hierarchical mechanism.

## References to Supporting Drafts

This unified draft compresses material developed in the following working documents:

- `结构存在论_物理假设与推论提取.md`
- `结构存在论评论分析.md`
- `结构存在论_数值实验与方法学讨论整合纪要.md`
- `PredictionA_MainPaper_Unified_Draft.md`
- `PredictionB_MainPaper_Unified_Draft.md`

The main numerical support currently comes from the pairwise validation, pooled regression, cross-family matched validation, and switch-enhancement outputs under `理论体系/poset_phase/outputs_exploratory/`.
