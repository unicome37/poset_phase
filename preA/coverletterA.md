**Date:** March 12, 2026

**To:** The Academic Editor, *Entropy* (MDPI)

**Re:** Submission of manuscript entitled "Dimensional Selection Without Dimensional Priors: 4D Lorentzian Dominance in Finite Causal Poset Ensembles Under Consistency-Based Actions"

---

Dear Editor,

I am pleased to submit the enclosed manuscript for consideration in *Entropy*. This paper is a companion to a concurrently submitted work "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy" (Prediction B), which establishes the existence of a bounded geometric phase transition in finite causal posets and introduces a non-target-anchored consistency penalty as a replacement for dimension-specific priors.

**Motivation and scope fit.** This paper asks a second-order question arising from the companion work: among Lorentzian-like structures of different dimensionality, *which dimension does a consistency-based action prefer?* The central objects are combinatorial entropy $H(P) = \log |L(P)|$ of finite posets and the competition between entropy maximization and structural consistency penalties — a free energy competition that falls squarely within *Entropy's* scope on information-theoretic and statistical-mechanical foundations.

**Summary of contributions.** The paper follows a "diagnose → ablate → replace → confirm" narrative:

1. **Diagnosis.** Under the original target-anchored action ($A_2$ full), a 4D Lorentzian-like family and KR posets are locked in a near-perfect tie (36 vs 35 wins out of 98 configurations). The 4D family does *not* achieve stable dominance — the negative result is presented honestly.

2. **Ablation.** Systematic removal of action components identifies a single term ($g_{\text{dim}}$, the target-anchored dimension proxy) as the dominant source of low-dimensional bias. This parallels the companion paper's finding that the same term is the critical backbone for the 2D phase transition.

3. **Replacement.** Replacing $g_{\text{dim}}$ with the non-target-anchored consistency penalty $g_{\text{con}}$ yields systematic 4D dominance across the tested grid (92/98 configurations on a $14\,N \times 7\,\gamma$ grid), with the 4D signal surviving inside a full geometric action package.

4. **Scaling and robustness.** The margin of victory shows a strong upward trend from +7 ($N = 20$) to +57 ($N = 72$), with finite-size dips at $N = 44$ and $N = 64$. Three independent generator seeds at $N = 68$ and $N = 72$ confirm 100% win rate (42/42 configurations) under consistency variants, versus 43% under the original action.

**Novelty.** This is the first finite-size numerical evidence that dimensional consistency constraints, without dimension-specific priors, naturally favor higher-dimensional Lorentzian structures in discrete causal order ensembles — a nontrivial step toward turning "why 3+1?" from a metaphysical presupposition into a quantitatively testable structural question.

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