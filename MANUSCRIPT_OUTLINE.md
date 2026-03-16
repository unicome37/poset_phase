# Bounded Geometric Phase Transition in Finite Causal Posets: Exact Computation and Ablation Analysis

## Manuscript Outline (Letter format, targeting *Entropy* or *Foundations of Physics*)

---

## Title Options

1. "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy"
2. "Emergence of Lorentzian-like Order in Action-Weighted Poset Ensembles"
3. "Ablation-Verified Geometric Dominance in Discrete Causal Structures"

**Recommended**: Option 1 — most precise, makes the key result (bounded + exact) immediately visible.

---

## Abstract (~150 words)

**Structure**: Problem → Method → Key Result → Implication

- **Problem**: In causal set theory, Kleitman–Rothschild (KR) posets dominate the entropy landscape of finite partial orders, presenting an obstacle for geometric phase emergence. Whether action-weighted ensembles can produce a bounded phase transition favoring geometric structures over KR remains open.
- **Method**: We construct a 7-family poset ensemble with exact linear extension counts (N = 10–44) and a decomposable geometric penalty. Three action paths (neutral-only, neutral+geometric, geometric-only) are compared.
- **Key Result**: Under the combined action, the critical coupling γ_c remains O(1) for all tested N. Ablation identifies a minimal backbone of two structural constraints. A non-target-anchored replacement for the dimension proxy preserves γ_c, establishing non-circularity.
- **Implication**: Finite-size evidence for a geometric phase transition in discrete causal orders, driven by structural self-consistency rather than dimension-specific priors.

---

## 1. Introduction (~500 words)

### 1.1 The Entropy Problem in Discrete Causal Orders
- Causal set theory (CST) posits spacetime as a locally finite partial order
- Kleitman–Rothschild (1975): asymptotically almost all finite posets are 3-layered, maximally entropic, and non-geometric
- This creates a "KR entropy catastrophe": naive ensemble averages are dominated by non-geometric structures
- Carlip (2017, CQG): proposed midpoint-scaling arguments to suppress KR, but no explicit action-based phase transition demonstrated at finite N

### 1.2 Our Approach
- Construct explicit poset ensembles with 7 candidate families
- Compute exact linear extension counts (entropy) using DP, avoiding SIS approximation errors
- Define decomposable action with separable neutral and geometric penalties
- Ask: does there exist a bounded γ_c such that geometric structures become competitive?

### 1.3 Key Contribution
- First explicit finite-size computation of γ_c across N = 10–44 with exact entropy
- Systematic ablation of geometric sub-terms
- Non-circular replacement: dim_proxy → dim_consistency
- Evidence that geometric phase emergence is driven by structural self-consistency, not dimension-specific priors

---

## 2. Framework (~800 words)

### 2.1 Poset Ensemble
- Definition of Poset data structure (DAG adjacency, transitive reduction)
- 7 candidate families with generators (Table: family name, structure description, key parameters)
- Sample sizes and generation protocol

### 2.2 Entropy: Exact Linear Extensions
- DP algorithm over antichains for counting linear extensions
- Computational complexity: Lor2D remains sub-second to N=48, about 3.8 s at N=72, about 33.8 s at N=88, and about 127.4 s at N=104; KR becomes expensive at N≥40. Use [lor2d_exact_timing_frontier.png](/d:/Kiro/理论体系/poset_phase/outputs_exploratory/lor2d_exact_timing_frontier/lor2d_exact_timing_frontier.png) and [lor2d_exact_timing_frontier.csv](/d:/Kiro/理论体系/poset_phase/outputs_exploratory/lor2d_exact_timing_frontier/lor2d_exact_timing_frontier.csv) as the compact frontier summary.
- Manuscript figure set now includes [fig4_exact_timing_frontier.pdf](/d:/Kiro/理论体系/poset_phase/manuscript_figures/fig4_exact_timing_frontier.pdf) and [fig4_exact_timing_frontier.png](/d:/Kiro/理论体系/poset_phase/manuscript_figures/fig4_exact_timing_frontier.png), showing the Lor2D frontier against KR/TP confirmatory exact timings.
- Manuscript figure set now also includes [fig5_mixed_lor2d_vs_kr.pdf](/d:/Kiro/理论体系/poset_phase/manuscript_figures/fig5_mixed_lor2d_vs_kr.pdf), [fig5_mixed_lor2d_vs_kr.png](/d:/Kiro/理论体系/poset_phase/manuscript_figures/fig5_mixed_lor2d_vs_kr.png), and source table [fig5_mixed_lor2d_vs_kr.csv](/d:/Kiro/理论体系/poset_phase/manuscript_figures/fig5_mixed_lor2d_vs_kr.csv), showing the near-wall mixed A2_full comparison between Lor2D and KR at N=52 and 56.
- This asymmetry is a structural fact (antichain width), not an algorithmic artifact

### 2.3 Action Paths
- A1 = −βH + γ · I_neutral
- A2 = −βH + γ · (I_neutral + I_geometric)
- A3 = −βH + γ · I_geometric
- Neutral penalty I_neutral: comparable fraction, degree variance, layer entropy
- Geometric penalty I_geometric: 7 sub-terms (Table with names and physical motivation)

### 2.4 Phase Transition Criterion
- Define γ_c as the value where Score(Lor2D, γ) = Score(KR, γ)
- If γ_c exists and is O(1) as N grows → evidence for bounded phase transition
- If γ_c diverges → geometric penalty must grow without bound to overcome entropy

---

## 3. Results (~1200 words)

### 3.1 Confirmatory γ_c Curve (Main Result)

**Table 1**: γ_c(N) for N = 10, 12, 14, 16, 20, 24, 28, 32, 36, 40, 44

| N | γ_c (frozen exact) |
|---|---|
| 10 | 9.19* |
| 12 | 0.744 |
| 14 | 0.781 |
| 16 | 0.693 |
| 20 | 0.146 |
| 24 | 0.262 |
| 28 | 0.995 |
| 32 | 0.694 |
| 36 | 0.570 |
| 40 | 0.248 |
| 44 | 0.149 |

*N=10 outlier due to finite-size effects (small poset regime).

**Figure 1**: γ_c(N) plot showing bounded O(1) behavior.

**Key observation**: Under A1 (neutral only), no γ_c exists at any N — Lor2D never overtakes KR. The geometric penalty is necessary.

### 3.2 Geometric Sub-Term Ablation

**Table 2**: Effect of removing individual geometric sub-terms on γ_c existence

| Removed sub-term | γ_c survival (N=20–44) | Classification |
|---|---|---|
| None (A2 full) | All N | Baseline |
| geo_width_height | Fails from N≥36 | **Critical driver** |
| geo_dim_proxy_penalty | Fails at N=36,44 | **Critical driver** |
| geo_interval_shape | All N (shifted higher) | Enhancer |
| geo_comparability_window | All N | Non-critical |
| geo_cover_density | All N | Non-critical |
| geo_interval_profile | All N | Non-critical |
| geo_layer_smoothness | All N | Non-critical |
| A1 (all geometric removed) | No γ_c at any N | Baseline control |

**Combined ablation**: Removing both critical drivers simultaneously → γ_c vanishes at all N. Retaining only interval terms → γ_c also vanishes.

### 3.3 Non-Circularity: dim_consistency Replacement

**The circularity concern**: `dim_proxy_penalty` penalizes deviation from d=2, which could be circular when the test target is a 2D structure.

**Resolution**: Replace `dim_proxy_penalty` with `dim_consistency_penalty`, which penalizes local-global dimension inconsistency without anchoring to any target dimension.

**Table 3**: γ_c under different backbone configurations

| Configuration | N=20–40 | N=44 | Notes |
|---|---|---|---|
| A2 full (7 terms) | ✓ all | ✓ | Baseline |
| replace dim→consistency | ✓ all | ✓ | Non-circular, same magnitude |
| wh + consistency only | ✓ all | ✗ | Minimal non-circular backbone |
| wh + consistency (1.3× weight) | ✓ all | ✓ | Weight-adjusted |
| width_height only | N=24,28 only | ✗ | Insufficient alone |

**Key conclusion**: The phase transition window does not require a target-dimension prior. A structurally self-consistent dimension constraint suffices.

---

## 4. Discussion (~600 words)

### 4.1 Physical Interpretation
- The two essential constraints have independent physical motivation:
  - `width_height_balance`: causal chain depth non-degeneracy (any dimension)
  - `dim_consistency`: dimensional scale coherence (no target required)
- Lor2D wins not because the action "knows" d=2, but because Lor2D structures naturally satisfy both constraints
- This is emergence, not prior
- Evidence chain should be written explicitly:
  - confirmatory exact window `N=10–44` establishes bounded `γ_c`
  - one-sided exact frontier shows `Lor2D` remains tractable up to `N=104`
  - mixed near-wall `N=52/56` extends the competition without overclaiming a crossing

### 4.2 Relation to Prior Work
- Carlip (2017): argued for KR suppression via midpoint scaling; our work provides the first explicit finite-size γ_c computation
- Surya (2012): CST partition function formalism; our ensemble is a discrete realization
- Loomis & Carlip (2018): 2D causal set simulations; our approach uses exact enumeration rather than Markov chain sampling

### 4.3 Limitations and Outlook
- **Finite-size**: N = 10–44 is far from thermodynamic limit. The key open question is whether γ_c remains bounded as N → ∞.
- **Boundary wording**: `N=52/56` mixed scans should be presented as near-degeneracy evidence, not as confirmatory continuation of the exact bounded-window result.
- **Family coverage**: 7 families do not exhaust all posets. The claim is "competitive phase among tested families," not global dominance.
- **Weight sensitivity**: The consistency-only backbone at N=44 requires mild weight adjustment (1.3×), suggesting the minimal backbone is at criticality edge for larger N.
- **Exploratory extensions** (Supplemental): Two-sieve framework, MLR survivor analysis, and paired locality validation provide additional structural context but are not part of the confirmatory chain.

---

## Supplemental Material

### S1. Generator Definitions
Detailed construction for all 7 families.

### S2. Complete Ablation Tables
Full N-by-N γ_c values for all ablation variants.

### S3. Exact Timing Benchmark
Runtime scaling for Lor2D, KR, TP across N.

### S4. Two-Sieve Framework (Exploratory)
Compressibility window definition, switch_zscore, MLR survivor analysis.

### S5. Paired Locality Validation (Exploratory)
locality_delta ↔ ΔlogH: Pearson r=0.759, p≈0.0005.

### S6. Non-Cyclic Replacement Details
dim_consistency_penalty definition, weight scan, trend analysis.

---

## Figures List

1. **Fig. 1**: γ_c(N) for N=10–44 under A2 (main result)
2. **Fig. 2**: Ablation summary — bar chart showing γ_c survival pattern
3. **Fig. 3**: Non-circular replacement comparison — γ_c under A2_full vs replace_dim_consistency
4. **Fig. 4**: Exact timing frontier — single-sided exact tractability of `Lor2D` versus KR/TP timing wall
5. **Fig. 5**: Near-wall mixed extension — `Lor2D` vs `KR` at `N=52,56`, showing high-γ near-degeneracy without crossing

---

## Target Journals (prioritized)

1. **Entropy (MDPI)** — Best fit: interdisciplinary, accepts combinatorial/statistical physics, no physics credentials required, fast review (~6 weeks), IF ~2.1
2. **Foundations of Physics (Springer)** — Good fit: accepts foundational/non-mainstream rigorous work, CST papers published here
3. **Physical Review Research (APS)** — Open access, cross-disciplinary, lower barrier than PRL/PRD
4. **Journal of Physics A: Mathematical and Theoretical (IOP)** — Mathematical physics framing
5. **Classical and Quantum Gravity (IOP)** — Core CST journal, but may require arXiv cross-listing

---

## Estimated Length

- Main text: ~3500 words (4–5 journal pages)
- Supplemental: ~3000 words
- Figures: 3 main + 2–3 supplemental
- Tables: 3 main + 2–3 supplemental
