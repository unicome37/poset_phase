# Project Status (One Page) — poset_phase / outputs_carlip

Date: 2026-03-30
Scope: Post-Carlip rebuild → “Two-layer screening” (S_BD admissibility + S_MD identity)

Update note (2026-03-30):
- F2 margin-aware refit now supports the safe onset claim **`N ≥ 10`** under a fixed-reference protocol.
- FLRW `κ=1.0` should now be treated as a **boundary-sensitive** background case under **background-dependent robustness**, not as evidence of global failure.

## 1) Where We Are Now (Current Thesis)

The project has converged from the original F7/logH-based narrative to a two-layer screening architecture:

- Layer 1 (linear admissibility): Benincasa-Dowker action `S_BD = c^T ΔC` filters grossly non-physical curvature content but cannot uniquely select Lor4D.
- Layer 2 (quadratic identity): a minimum-distortion / Mahalanobis structural action
  `S_MD(P,N) = (I(P)-μ(N))^T Σ(N)^{-1} (I(P)-μ(N))`
  uniquely identifies Lor4D as the closest family to the Lor4D reference manifold in a 3-feature space.

Key design choice: replace logH with purely causal-geometric observables:
`I(P) = (d_eff, C1/C0, width_ratio)`.

Primary claim target (paper): the poset landscape is not governed by a single action minimum; Lor4D emerges as the unique survivor of a gate (S_BD) plus an identity basin (S_MD).

## 2) What Is Locked In (High-Confidence Results)

### A. F7/logH is retired (and why)
- On the expanded 17-family space, F7 fails at N ≥ 28 and Lor4D is no longer top-ranked.
- The failure mechanism is structural: logH grows with N and dominates; “wall” terms decay; KR_2layer is a decisive counterexample in the (d_eff, R) plane.
Sources: `诊断总结_F7_17家族.md`, `f7_17family_test.md`, `综合诊断_两轮实验.md`.

### B. “Well” is the correct architecture (non-monotone discrimination)
- Any monotone “more is better” LSD attempt is defeated by layered families and/or KR pathologies.
- A quadratic well around Lor4D centroids (triple well) succeeds across N and families.
Sources: `综合诊断_两轮实验.md`, `lsd_17family_test.md`, `lsd_well_17family_test.md`, `lsd_well_n_adapted.md`.

### C. Zero-parameter identity via Mahalanobis (core advance)
- Mahalanobis LSD replaces hand-tuned weights with Σ(N)^{-1} from Lor4D ensemble statistics.
- Under the fixed-reference F2 protocol (20 seeds × 120 reps, separate reference ensemble), Lor4D is ranked #1 from **N = 10** onward; margins then grow rapidly with N.
Sources: `mahalanobis_lsd.md`, `prediction_b_cross_validation.md`, `prediction_b_seed_reproducibility.md`, `large_n_extreme.md`.

### D. Basin deepening + reference manifold as “theory object”
- Empirically, the Lor4D reference manifold (μ(N), Σ(N)) sharpens with N (CLT-like σ^2 ~ N^{-1} behavior).
- Separation (margin/gap) diverges at large N in extreme-N tests up to N=1024.
Sources: `mu_trajectory.md`, `basin_deepening_results.txt`, `large_n_extreme.md`, `MASTER_NARRATIVE.md`.

### E. Adversarial robustness (expanded family library)
- A 25-family library (17 standard + 8 adversarial) has been tested.
- LSD-Well keeps Lor4D rank #1 at all tested N in that report.
- The old Mahalanobis `N=16` exception now belongs to the **historical CV/shared-reference regime**; the current safe onset statement is `N ≥ 10` under the fixed-reference F2 protocol.
Source: `expanded_family_results.txt`.

## 3) What Is Still “Open” (Main Risks / Gaps)

### Risk 0: Background-dependent curvature response (newly elevated)
- Split low-N F3 runs (N=256/512) now separate the C2 picture by background:
  - de Sitter: pass (top-2 maintained)
  - weak-field Schwarzschild: pass (top-2 maintained)
  - matter-FLRW: hard-fail triggered at $\kappa=1.0$ (repeatable top-2 loss under current threshold)
- Completed high-N runs (N=768/1024) show **partial recovery** for FLRW at $\kappa=1.0$ with failure ratio $0.3 < 0.5$.
- P0 metric-faithful PhaseA refines the same picture: metric branch `11/60 = 0.183`, proxy branch `0/60`.
- Implication: avoid any unified “mild-curvature robustness” wording; use **background-dependent robustness**.
- Immediate priority: package FLRW `κ=1.0` into a one-page reviewer-defense brief, then decide whether PhaseB (`N=1024/1536/2048`) is needed.

### Risk 1: Small-N exception for pure Mahalanobis
- The old N=16 Mahalanobis exception should now be treated as a **historical diagnostic** tied to CV/shared-reference contamination and small reference ensembles.
- Current safe claim: under the fixed-reference protocol, Mahalanobis uniquely selects Lor4D for **N ≥ 10**.
Source: `n16_mahalanobis_cv_rootcause.md`, plus `expanded_family_results.txt` (N=16 line).

Implication for claims:
- Avoid “all N” absolutes for Mahalanobis; state the current safe threshold **`N ≥ 10` under the fixed-reference protocol** and keep the older CV instability explicitly labeled as historical/diagnostic.

Principled small-N treatment (attempted):
- Anchoring only the dimension mean to d*=4 does not reliably eliminate the N=16 intruder.
- The N=16 winner is sensitive to how well (mu,Sigma) are estimated; increasing the Lor4D ensemble size at N=16 stabilizes baseline/full Mahalanobis to Lor4D #1 in local seed-sweeps.
Source: `SMALL_N_TREATMENT_NOTE.md`, `mahalanobis_n16_stability_results.md`.

### Risk 2: Gradient-bridge overclaim (must remain auxiliary)
- Magnitude-level alignment can be stable, but sign is not stable due to Jacobian pseudoinverse ill-conditioning and parameterization sensitivity.
- Keep “order-raising relatedness” language; do not claim pointwise gradient identity.
Sources: `DISCUSSION_THEORY_IMPLICATIONS.md`, `gradient_alignment_v2.md`, `gradient_signflip_diagnostic.txt`.

### Risk 3: “Zero parameter” messaging must be precise
- “Zero parameter” is true for S_MD once μ(N), Σ(N) are defined as Lor4D ensemble objects.
- But paper must be explicit about how μ(N), Σ(N) are estimated (finite reps, sampling protocol) and how results change under CV/seed changes.
Sources: `SMD_OPERATOR_LETTER.md`, `INTRODUCTION_DRAFT.md`, `mahalanobis_lsd.md`.

## 4) What We Should Say (Safe Paper-Grade Positioning)

- Core: two-layer architecture with different algebraic orders (linear gate vs quadratic identity), nearly orthogonal signals.
- Main evidence pillars:
  - admissibility/identity separation,
  - basin deepening scalings,
  - expanded/adversarial family robustness,
  - reference manifold collapse with N.
- Auxiliary only: gradient alignment → “order-raising transition” intuition.
Source synthesis: `MASTER_NARRATIVE.md`, `TWO_LAYER_SCREENING_THEORY.md`, `DISCUSSION_THEORY_IMPLICATIONS.md`.

## 5) Immediate Deliverables Already in Draft Form

- Letter-style short paper draft: `SMD_OPERATOR_LETTER.md`
- Full-paper entry points: `INTRODUCTION_DRAFT.md`, `DISCUSSION_THEORY_IMPLICATIONS.md`
- Unified “one sentence + tables”: `MASTER_NARRATIVE.md`
- FLRW one-page defense brief: `FLRW_KAPPA1_DEFENSE_BRIEF_20260330.md`
- Reproducibility guide: `REPRO_RUNBOOK_20260330.md`

## 6) Newly Added Status Artifacts (2026-03-28)

- `outputs_carlip/F1_ONEPAGE_SUMMARY_20260328.md`
- `outputs_carlip/F3_LOWN_SPLIT_SUMMARY_20260328.md`
- `outputs_carlip/CURRENT_REVIEW_RISK_STATUS_20260328.md`
