**Date:** March 13, 2026

**To:** The Academic Editor, *Entropy* (MDPI)

**Re:** Submission of manuscript entitled "Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study"

---

Dear Editor,

I am pleased to submit the enclosed manuscript for consideration in *Entropy*. This paper is the third in a series of three companion papers: Prediction B ("Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy") and Prediction A ("Dimensional Selection Without Dimensional Priors"), both concurrently submitted to *Entropy*. The present paper (Prediction C) addresses the structural mechanism question: *why* do certain poset families achieve systematically lower combinatorial entropy?

**Motivation and scope fit.** The central objects are combinatorial entropy $H(P) = \log |L(P)|$ of finite causal posets and a pre-registered Hierarchy Integration Index (HII) — a composite z-score of five structural depth observables. The interplay between structural order and combinatorial counting connects directly to *Entropy's* scope on information-theoretic and statistical-mechanical foundations of discrete systems.

**Summary of contributions.**

1. **Three-tier validation design.** We test the HII–entropy hypothesis across three independent statistical tiers: (i) all-family partial correlation (8 families, $N = 10$–$16$, 320 samples with exact entropy), (ii) matched-pair $\Delta$-analysis between Lor2D and MLR families (46 pairs, $N = 30$–$56$), and (iii) coarse-graining identity-stability linkage (92 samples). Each tier addresses different potential confounds.

2. **Strong, stable correlation.** At fixed $N$, deeper hierarchy — principally `layer_count` and `mean_layer_gap` — correlates negatively with $\log H$: matched-pair $r = -0.834$, stable across three levels of filter stringency (variation $< 0.005$). The association extends to coarse-graining stability ($r = -0.874$, classifier-contingent).

3. **Simpson's Paradox diagnosis.** The naïve cross-$N$ correlation is positive ($r = +0.336$) and reverses sign only after controlling for $N$ ($r = -0.578$). We provide a complete decomposition of this paradox — a methodological finding of independent value for researchers working with size-heterogeneous discrete structure data.

4. **Component decomposition.** The five-component HII composite never exceeds its best single constituent (`layer_count`), revealing that the signal is carried by a two-component depth pair.

5. **Hard limitations stated.** Family-specific HII reversals (KR-like, absolute-layered), classifier-dependent Tier 3, and weak near-wall statistical power at $N \geq 52$ are documented as explicit boundaries.

**Novelty.** No previous study of causal set ensembles has: (a) defined a pre-registered composite structural index for hierarchy depth; (b) validated it across three independent statistical designs spanning two disjoint $N$ ranges; (c) diagnosed a Simpson's Paradox in the structural–entropy relationship; or (d) extended the association chain from hierarchy to coarse-graining stability.

**Declarations.**

- This manuscript is original, has not been published elsewhere, and is not under consideration by another journal.
- Code and data: https://github.com/unicome37/poset_phase (Zenodo archive: https://doi.org/10.5281/zenodo.18980657).
- AI Assistance: Large language model assistants (Claude, Anthropic) contributed to code implementation, data analysis, and manuscript drafting. All scientific decisions were made by the human author. This is fully disclosed in the manuscript.
- No external funding was received.
- The author declares no conflicts of interest.

Thank you for considering this work. I look forward to your editorial feedback.

Sincerely,

Gang Zhang  
Independent Researcher  
unicome@gmail.com
