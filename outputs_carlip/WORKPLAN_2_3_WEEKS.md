# Work Plan (2–3 Weeks) — poset_phase / Two-Layer Screening

Date: 2026-03-27
Goal: convert the current numerical synthesis into a submission-ready package with precise claims, strong robustness, and clean paper structure.

## Week 1 (Stabilize Claims + Close the Small-N Hole) — ✅ COMPLETED 2026-03-27

### 1.1 Lock "claim boundaries" — ✅ DONE
- **Turn-on boundary confirmed: N ≥ 14** (previously expected N ≥ 20)
- REPS=80, 25 families, 10 independent seeds
- N=14–32: 10/10 Lor4D rank #1, margin +1.28 → +4.59
- N=12: 6/10 (physical resolution limit, intruders: Lor5D + KR_2layer)
- Evidence: `mahalanobis_n_boundary_turnon.md`, `fig_margin_vs_n.png`, `n12_failure_diagnosis.md`

### 1.2 Small-N treatment — ✅ DECIDED
- **Decision**: Claim N ≥ 14 with resolution-limit explanation for N=12
- N=12 failures are physical: d_eff overlap between Lor4D/Lor5D at 12 sprinkled points
- Not a statistical artifact (cond(Σ) < 50, not pathological)
- No tunable knobs introduced

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
