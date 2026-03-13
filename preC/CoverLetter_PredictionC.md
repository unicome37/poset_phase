# Cover Letter — MDPI *Entropy* Submission

> MDPI *Entropy* 投稿通过在线系统 SuSy 完成：https://susy.mdpi.com/
> 下方 Cover Letter 粘贴至投稿系统的 Cover Letter 栏。

---

**Date:** [INSERT DATE]

**To:** The Academic Editor, *Entropy* (MDPI)

**Re:** Submission of manuscript entitled "Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study"

---

Dear Editor,

We are pleased to submit the enclosed manuscript for consideration as a Research Article in *Entropy*.

**Motivation and scope fit.** This paper investigates which intrinsic structural properties of finite causal partial orders (posets) predict their combinatorial entropy — quantified by the exact logarithmic count of linear extensions, $H(P) = \log |\mathcal{L}(P)|$. The central object of study is a pre-registered Hierarchy Integration Index (HII), a composite z-score of five depth-related observables, tested against entropy across three independent statistical designs. The interplay between structural order and combinatorial counting connects directly to *Entropy's* scope on information-theoretic and statistical-mechanical foundations of discrete systems.

**Summary of contributions.**

1. **Three-tier validation design.** We test the HII–entropy hypothesis across three independent statistical tiers: (i) all-family partial correlation (8 families, $N = 10$–$16$, 320 samples), (ii) matched-pair $\Delta$-analysis between Lor2D and MLR families (46 pairs, $N = 30$–$56$), and (iii) coarse-graining identity-stability linkage (92 samples). Each tier addresses different potential confounds.

2. **Strong, stable correlation.** At fixed $N$, deeper hierarchy—principally `layer_count` and `mean_layer_gap`—correlates negatively with $\log H$: matched-pair $r = -0.834$, stable across three levels of filter stringency (variation $< 0.005$). The association extends to coarse-graining stability ($r = -0.874$, classifier-contingent).

3. **Simpson's Paradox diagnosis.** The naïve cross-$N$ correlation is positive ($r = +0.336$) and reverses sign only after controlling for $N$ ($r = -0.578$). We provide a complete decomposition of this paradox, establishing $N$ as the dominant sign-determining confound—a methodological finding of independent value for researchers working with size-heterogeneous discrete structure data.

4. **Component decomposition.** The five-component HII composite never exceeds its best single constituent (`layer_count`), revealing that the signal is carried by a two-component depth pair rather than distributed across all five observables.

5. **Hard limitations stated.** Family-specific HII reversals (KR-like, absolute-layered), classifier-dependent Tier 3, and weak near-wall statistical power at $N \geq 52$ are documented as explicit boundaries of the current evidence.

**Relation to companion papers.** This manuscript is the third in a series of three companion papers, all submitted to *Entropy*:

- **Prediction B** [1]: Establishes that a 2D Lorentzian-like family becomes competitive against KR-like high-entropy posets under an action functional, with a bounded critical coupling $\gamma_c = O(1)$ over $N = 10$–$44$.
- **Prediction A** [2]: Shows that a 4D Lorentzian-like family achieves systematic dominance under consistency-based actions over $N = 20$–$72$.
- **Prediction C** (this paper): Identifies hierarchy depth as the structural mechanism underlying the entropy ordering observed in Predictions A and B.

Each paper is self-contained and can be evaluated independently.

**Novelty.** While prior work on causal set theory has studied entropy scaling and phase transitions in discrete causal orders, no previous study has: (a) defined a pre-registered composite structural index for hierarchy depth; (b) validated it across three independent statistical designs; (c) diagnosed a Simpson's Paradox in the structural–entropy relationship; or (d) extended the association chain to coarse-graining stability.

**Declarations.**

- This manuscript is original, has not been published elsewhere, and is not under consideration by another journal.
- **Code and data**: https://github.com/unicome37/poset_phase (archived at Zenodo: DOI 10.5281/zenodo.18980657). All computational results are fully reproducible from the published code and configuration files.
- **AI Assistance**: Large language model assistants (Claude, Anthropic) contributed to code implementation, data analysis, and manuscript drafting. All scientific decisions, hypothesis formulation, and experimental design were made by the human author. This is fully disclosed in the manuscript.
- **Funding**: No external funding was received for this research.
- **Conflicts of Interest**: The author declares no conflicts of interest.

**Suggested reviewers** (optional):

1. Sumati Surya (Raman Research Institute) — expert on causal set quantum gravity and entropy in discrete spacetimes
2. Graham Brightwell (London School of Economics) — expert on linear extension counting and poset combinatorics
3. Seth Major (Hamilton College) — expert on discrete quantum gravity approaches

Thank you for considering this submission.

Sincerely,

Gang Zhang
Independent Researcher
