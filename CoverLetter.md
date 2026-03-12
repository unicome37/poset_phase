**Date:** March 12, 2026

**To:** The Academic Editor, *Entropy* (MDPI)

**Re:** Submission of manuscript entitled "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy"

---

Dear Editor,

We are pleased to submit the enclosed manuscript for consideration in *Entropy*.

**Motivation and scope fit.** This paper studies how combinatorial entropy — quantified by exact linear extension counts of finite partially ordered sets — competes with structural order constraints to determine which causal configurations emerge as dominant phases. The central mathematical object is $H(P) = \log |L(P)|$, the logarithm of the number of linear extensions of a poset $P$, which serves as the discrete analogue of thermodynamic entropy in this setting. The interplay between entropy maximization and geometric penalty terms directly mirrors the free energy competition familiar in statistical mechanics. We believe this falls squarely within *Entropy's* scope on information-theoretic and statistical-mechanical foundations.

**Summary of contributions.**

1. **Exact entropy computation.** We construct a 7-family poset ensemble and compute exact linear extension counts via antichain-based dynamic programming for posets with $N = 10$–$44$ elements. This exact approach avoids the systematic biases inherent in Monte Carlo approximations and reveals that computational cost itself varies dramatically across structural types — a structural fact about the posets, not an algorithmic artifact.

2. **Bounded phase transition threshold.** Under a decomposable action combining entropic and geometric terms, we show that the critical coupling $\gamma_c$ — the minimum geometric penalty strength needed for Lorentzian-like posets to overtake the entropically dominant Kleitman–Rothschild (KR) structures — remains bounded and $O(1)$ across all tested sizes. This constitutes the first finite-size evidence of such a bounded threshold.

3. **Systematic ablation.** Ablation of seven individual geometric sub-terms identifies a minimal two-component backbone (global layer-shape balance + dimensional scale information) necessary and sufficient for the phase transition.

4. **Non-circular replacement.** The dimension-dependent backbone component can be replaced by a non-target-anchored alternative — penalizing local–global dimension *inconsistency* rather than deviation from a specific target $d = 2$ — with the phase transition preserved at the same order of magnitude. This directly addresses the concern that geometric penalties might circularly "reward looking geometric."

**Novelty.** While prior work (Carlip 2017; Loomis & Carlip 2018) has argued qualitatively and via MCMC for KR suppression in causal set models, no previous study has provided an explicit $\gamma_c(N)$ curve with exact entropy across a systematic range of sizes, combined with full ablation and non-circular replacement analysis.

**Declarations.**

- This manuscript is original, has not been published elsewhere, and is not under consideration by another journal.
- Code and data: https://github.com/unicome37/poset_phase (Zenodo archive: https://doi.org/10.5281/zenodo.18963421).
- AI Assistance: Large language model assistants (Claude, Anthropic) contributed to code implementation, data analysis, and manuscript drafting. All scientific decisions were made by the human author. This is fully disclosed in the manuscript.
- No external funding was received.
- The author declares no conflicts of interest.

Thank you for considering this work. We look forward to your editorial feedback.

Sincerely,

Gang Zhang  
Independent Researcher  
unicome@gmail.com