# Work Plan (2–3 Weeks) — poset_phase / Two-Layer Screening

Date: 2026-03-27
Goal: convert the current numerical synthesis into a submission-ready package with precise claims, strong robustness, and clean paper structure.

## Week 1 (Stabilize Claims + Close the Small-N Hole)

### 1.1 Lock “claim boundaries” (write once, reuse everywhere)
- Define canonical claim set:
  - `S_BD` is admissibility (necessary, not sufficient).
  - `S_MD` (Mahalanobis) is identity (sufficient for N ≥ N0).
  - Joint screening selects Lor4D uniquely in the tested library.
- Pick and enforce a safe N-threshold for Mahalanobis (default: `N ≥ 20`), unless the following work upgrades it.
- Propagate wording to: `SMD_OPERATOR_LETTER.md`, `INTRODUCTION_DRAFT.md`, `DISCUSSION_THEORY_IMPLICATIONS.md`, `MASTER_NARRATIVE.md`.

### 1.2 Decide a principled small-N treatment (choose one, implement, document)
Target: handle the `N=16` Lor5D “intruder” for pure Mahalanobis without introducing tunable knobs.

Options:
- (Preferred) Add an explicit “resolution limit” clause:
  - State that `N=16` does not encode enough information to separate 4D vs 5D reliably.
  - Keep Mahalanobis claim as `N ≥ 20`.
- (If you want “all N”) Add a *principled* prior:
  - Use `d*=4` as a physical prior at small N (hybrid metric: Mahalanobis + fixed (d-4)^2 term, with a rule-based crossover N0).
  - Alternatively, use shrinkage Σ→(1-α)Σ+α·I with α fixed by an information criterion or analytic prescription (not hand tuned).

Status note (local tests, 2026-03-27):
- Anchoring only the d-component mean to `d*=4` did not reliably eliminate the `N=16` intruder by itself.
- The `N=16` top-1 winner is sensitive to finite-sample estimation of `(mu,Sigma)` when the Lor4D reference ensemble is small (e.g. 20 samples); with larger Lor4D ensembles the baseline/full Mahalanobis stabilizes to Lor4D #1 in seed sweeps.
- Reference: `outputs_carlip/SMALL_N_TREATMENT_NOTE.md` and `outputs_carlip/mahalanobis_n16_stability_results.md`.

Deliverable:
- 1–2 paragraphs “Small-N limitation / treatment” inserted into `SMD_OPERATOR_LETTER.md` and the full paper.
- A short note summarizing the decision and rationale in `DISCUSSION_THEORY_IMPLICATIONS.md` (limitations subsection).

### 1.3 Re-run or extend one robustness slice (only if needed)
- Expand N-grid around the boundary (e.g., N=16,18,20,24,28) for the Mahalanobis ranking to show the “turn-on” of uniqueness.
- Keep this as a single compact table for the Letter.

## Week 2 (Turn the Theory Object Into a Clean Story)

### 2.1 Reference manifold: from “fit” to “object”
- Consolidate μ(N), Σ(N) definitions:
  - sampling protocol (how many reps, seeds),
  - finite-size ansatz,
  - what is derived vs what is fitted.
- Make the asymptotic pedigrees crisp:
  - `d_eff → 4` (Myrheim-Meyer),
  - `C1/C0` (Beta integral framework),
  - `width` scaling (diamond cross-section / N^{-1/d} behavior).

Deliverable:
- One self-contained “Methods: reference manifold construction” section for the full paper.
- One compact paragraph + table for the Letter.

### 2.2 Minimality and feature story (keep it honest)
- Use `feature_ablation_test.md` to define what is “minimal complete basis” vs “critical feature”.
- Adjust messaging so it matches the data:
  - current ablation shows `d_eff` is critical; others mainly increase margin.
- If you want the stronger “triple is minimal” claim, update ablation experiments to demonstrate a genuine failure when dropping each (or soften the claim).

Deliverable:
- Update the “Limitations / choice of features” section (Letter + full paper) so it is defensible.

### 2.3 Relegate gradient bridge properly
- Move gradient-bridge material to a short Discussion subsubsection or appendix-like note.
- Remove any remaining “cos=0.97 proves derivation” style language per `DISCUSSION_THEORY_IMPLICATIONS.md`.

## Week 3 (Submission Packaging)

### 3.1 Letter finalization (CQG Letters / PRL-style)
- Ensure the Letter is self-contained:
  - 1 figure (two-layer picture or margin vs N),
  - 1 table (25 families + N grid summary),
  - 1 limitations paragraph (small N, finite N max, feature choice).
- Tighten abstract language:
  - replace “all N” with “N ≥ N0” if needed,
  - keep “zero free parameters” only for the finalized definition.

Deliverable:
- `SMD_OPERATOR_LETTER.md` v2 ready for LaTeX conversion.

### 3.2 Full paper skeleton stabilization
- Freeze section ordering and what goes where:
  - Introduction: entropy catastrophe + BD limits + motivation
  - Methods: family library + features + reference manifold
  - Results: S_BD gate, S_MD identity, basin deepening, adversarial tests
  - Discussion: two-layer meaning + order-raising (aux) + limitations

Deliverable:
- A single “paper map” page (can be appended to `MASTER_NARRATIVE.md` or kept separate).

## Concrete Definition of “Done”

- Claims are consistent across all drafts.
- Mahalanobis small-N behavior is either (a) cleanly bounded by N0 with a resolution-limit explanation, or (b) fixed by a principled rule that does not reintroduce tunable knobs.
- Letter draft reads as a short, defensible numerical discovery with explicit limits and strong robustness.
