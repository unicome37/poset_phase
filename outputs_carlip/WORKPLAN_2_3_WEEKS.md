# Work Plan (2–3 Weeks) — poset_phase / Two-Layer Screening

Date: 2026-03-30
Goal: convert the current numerical synthesis into a submission-ready package with precise claims, strong robustness, and clean paper structure.

> **2026-03-30 supersession note**
>
> This work plan is now partially superseded by the F2 margin-aware refit and P0 metric-faithful PhaseA results:
> - the safe onset claim is now **`N ≥ 10`**, not `N ≥ 14`;
> - the curvature wording should now be **`background-dependent robustness`**;
> - the current immediate deliverables are: (i) full manuscript wording sync, (ii) a one-page FLRW `κ=1.0` defense brief, and (iii) a reproducibility runbook.

## Week 1 (Stabilize Claims + Close the Small-N Hole) — ✅ COMPLETED 2026-03-27

### 1.1 Lock "claim boundaries" — ✅ DONE
- **Turn-on boundary superseded: N ≥ 10 under fixed-reference F2 protocol**
- F2 protocol: 20 seeds × 120 reps, separate reference ensemble
- N=10–24: 20/20 Lor4D rank #1; manuscript-safe onset already holds at N=10
- Older `N=12/14/16` instability belongs to shared-reference contamination diagnostics, not the current safe claim boundary
- Evidence: `f2_turnon_margin_refit.md`, `falsify_c1_turnon_refit.md`, `P1_STRONGER_THAN_OLD_DIAGNOSIS_20260330.md`

### 1.2 Small-N treatment — ✅ DECIDED
- **Decision (updated)**: Claim `N ≥ 10` under the fixed-reference protocol.
- Earlier N=12–16 instability should be described as a reference/test contamination artifact in the shared-seed CV regime, not as the current safe physical boundary.
- No tunable knobs are introduced; the change is protocol-level (independent reference ensemble + larger reference reps).

### 1.3 N-boundary robustness  DONE
- Full N=12,14,16,18,20,24,28,32 grid with 80 runs each
- Turn-on table and 4-panel diagnostic figure ready for Letter

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

### 3.0 Current front-loaded priorities (2026-03-30)
- Produce `FLRW_KAPPA1_DEFENSE_BRIEF_20260330.md`
- Produce `REPRO_RUNBOOK_20260330.md`
- Sweep remaining docs for old `N ≥ 14/16` and “uniform mild-curvature robustness” wording

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
- Mahalanobis small-N behavior is stated using the current safe threshold (`N ≥ 10` under the fixed-reference protocol), with old CV/shared-reference instability clearly labeled as historical/diagnostic.
- Letter draft reads as a short, defensible numerical discovery with explicit limits and strong robustness.
