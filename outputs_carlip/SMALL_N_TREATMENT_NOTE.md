# Small-N Treatment Note (N=16) — Attempt Toward “All N”

Date: 2026-03-27
Context: Mahalanobis `S_MD(P,N) = (I(P)-mu(N))^T Sigma(N)^{-1} (I(P)-mu(N))` sometimes shows an `N=16` Lor5D intruder in earlier robustness runs.

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

If we prefer the most conservative paper-grade statement, keep:
- “Mahalanobis selects Lor4D uniquely for `N >= 20`”, and treat `N=16` as a resolution-limit regime.

