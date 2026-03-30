# Small-N Treatment Note (Historical N=16 Issue → superseded by F2)

Date: 2026-03-30
Context: this note records the **historical** `N=16` Mahalanobis intruder issue from earlier robustness runs. The current safe claim is no longer derived from this note, but from the fixed-reference F2 protocol.

## What We Tried (Principled, No Hand-Tuned Knobs)

1. **Anchor only the dimension component** of the reference mean:
   `mu_anchor(N) = (4.0, mu_c(N), mu_w(N))`, keeping `Sigma(N)` from Lor4D.
   - Outcome: does **not** reliably eliminate the intruder by itself.

2. **Diagonal (Fisher) metric** as a robust alternative:
   replace `Sigma^{-1}` by `diag(1/var_i)` where `var_i` are Lor4D within-class variances.
   - Outcome: improves stability at `N=16` in our tests, but is not uniformly best at all N.

3. **Key observation: the “intruder” is sample-size sensitive at N=16**.
   With too few Lor4D samples used to estimate `(mu, Sigma)`, finite-sample fluctuation can invert the top-2 ordering between Lor4D and Lor5D.

## Evidence (Reproducible Local Runs)

- Low reps (`REPS=20`) across 10 seed-bases: baseline/full Mahalanobis is **not** always Lor4D at N=16.
- Higher reps (`REPS=80`) across the same 10 seed-bases: baseline/full Mahalanobis becomes **stable**:
  Lor4D is #1 in 10/10 seeds for baseline/full; diagonal is also 10/10.

See: `d:\Kiro\理论体系\poset_phase\outputs_carlip\mahalanobis_n16_stability_results.md`.

## Practical Resolution for “All N”

If the paper wants to claim “all N including N=16”, the cleanest *principled* move is:
- **Define** `mu(N), Sigma(N)` as ensemble objects estimated to a fixed precision, and
- **Use a larger Lor4D ensemble at N=16** when constructing `mu(16), Sigma(16)` (computational budget, not a tunable model knob).

If one wants to quote the historical conservative paper-grade statement, it was:
- “Mahalanobis selects Lor4D uniquely for `N >= 20`”, and treat `N=16` as a resolution-limit regime.

This wording is now **obsolete** and should be retained only as a historical record of the pre-F2 stage.

## 2026-03-30 Update: F2 margin-aware refit resolved this

The F2 runner (`f2_turnon_margin_runner.py`) uses `reference_reps=120` with a **separate seed offset** (`seed_base + 100000`) to build the Lor4D reference ensemble, completely independent from the test families. With this protocol:

- N=10: 20/20 Lor4D rank#1, min_margin=0.198, ci95_lower=0.268 > 0, max_cond_σ=54.9 < 60 → **manuscript-safe**
- N=12..24: 20/20 rank#1, margins monotonically increasing

**The N=16/20 conservative boundary is now obsolete.** The correct statement is: "Mahalanobis selects Lor4D uniquely for N ≥ 10 under the fixed-reference protocol (reference_reps=120, separate seed)." The old instability at N=12–16 was entirely due to reference/test contamination with small reps, not a physical limitation.

