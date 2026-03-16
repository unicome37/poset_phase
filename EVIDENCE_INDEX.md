# Claim → Evidence → Boundary Index
## Date: 2026-03-16  |  Covers: Prediction B + A

---

## PREDICTION B: Bounded γ_c Window for Lorentzian-like Causal Sets

### B.1 Core Claim
Under the structure action A₂, `lorentzian_like_2d` exhibits a bounded phase-transition
coupling γ_c against `KR_like`, with γ_c remaining finite and non-divergent across N = 10–44.

| Claim | Evidence | Boundary |
|---|---|---|
| γ_c exists for N=10–16 (exact, frozen) | `outputs_confirmatory/frozen_exact/gamma_c_report.csv` | Small-N only; confirmed under frozen protocol |
| γ_c exists for N=20–44 (exact, medium) | `outputs_confirmatory/medium_exact_scan/gamma_c_report.csv` | N=44 is the exact confirmatory frontier |
| γ_c shows no divergence trend N=10–44 | Combined: frozen_exact + medium_exact_scan gamma_c reports | Not a proof of boundedness; only absence of divergence |
| N=52–72 mixed: competition still active | `outputs_exploratory/prediction_b_nearwall_mixed_frontier/` | Mixed only; no new crossing confirmed |
| Near-degeneracy at N=52,56 | `prediction_b_nearwall_gap_summary.csv` in above dir | Gap < 1 at γ=2.0; does NOT prove crossing restored |

### B.2 Tier-1/Tier-2 Discrimination
| Claim | Evidence | Boundary |
|---|---|---|
| Tier-1/Tier-2 pool definition | `TIER_DEFINITION_FROZEN.md` | Frozen 2026-03-16 |
| aw/√N separates Lor2D from layered | `outputs_exploratory/tier_scaling_analysis/aw_scaling_tests.csv` | All p < 10⁻¹⁰; per-N all p < 0.005 |
| Layered families beat Lor2D under current action | `outputs_exploratory/prediction_b_candidate_expansion_layered/` | Tier-2 stress test; does not invalidate core B claim |
| Action discriminative ceiling exists | Implied by Tier-2 results | Future work: embed aw/√N or interval statistics in action |

### B.3 Computational Feasibility
| Claim | Evidence | Boundary |
|---|---|---|
| Lor2D exact extends to N=104 | `outputs_exploratory/lor2d_exact_timing_frontier/` | N≈100 enters minute-level; still feasible |
| Lor3D exact extends to N=56 | `outputs_exploratory/exact_timing_lor_push/exact_timing_benchmark.csv` | 280s/288s/327s for N=48/52/56 |
| Lor4D infeasible at N≥48 | Same benchmark | Memory overflow confirmed |
| KR exact wall at N≈44 | `outputs_confirmatory/exact_timing/exact_timing_summary.csv` | 56s at N=44; infeasible beyond |

### B.4 Robustness
| Claim | Evidence | Boundary |
|---|---|---|
| Weight sensitivity scan | `outputs_exploratory/prediction_b_weight_sensitivity/` | γ_c stable under ±30% weight perturbation |
| Weight sensitivity focus | `outputs_exploratory/prediction_b_weight_sensitivity_focus/` | Focused grid on sensitive weights |
| Seed sensitivity (c_threshold) | `outputs_exploratory/seed_sensitivity_c_threshold_*/` (6 dirs) | Multiple seed offsets at N=44,52 |

---

## PREDICTION A: 3+1-Dimensional Emergence under Consistency Replacement

### A.1 Core Claim
When the targeted dimension-proxy penalty is replaced with non-target-anchored consistency
constraints, `lorentzian_like_4d` systematically wins across N = 20–72.

| Claim | Evidence | Boundary |
|---|---|---|
| A2_full does NOT support 4D | `outputs_exploratory/prediction_a_pilot_fast/` | Lor2D wins under full action |
| dim_proxy_penalty identified as bias source | `outputs_exploratory/prediction_a_geometric_ablation/` | Single-component ablation confirms |
| Consistency replacement → Lor4D wins N=20–48 | `outputs_exploratory/prediction_a_dim_replacement/` | 30/35 at sp=4, 36/42 at sp=8 |
| Extension to N=44,48 near-wall | `outputs_exploratory/prediction_a_dim_replacement_n44_n48/` | 14/14 (consistency), 13/14 (multi) |
| Extension to N=52 mixed | `outputs_exploratory/prediction_a_n52_mixed/` | 7/7 Lor4D under consistency |
| Extension to N=56–72 mixed | `prediction_a_n56_mixed/` through `prediction_a_n72_mixed/` | Lor4D dominant; mixed oral |
| Phase summary visualization | `outputs_exploratory/prediction_a_phase_summary/` | Side-by-side winner maps |
| Margin summary | `outputs_exploratory/prediction_a_margin_summary/` | Quantified advantage gaps |

### A.2 Robustness
| Claim | Evidence | Boundary |
|---|---|---|
| sp=8 replication | `outputs_exploratory/prediction_a_dim_replacement_sp8/` | Results consistent with sp=4 |
| Robust scan (wider grid) | `outputs_exploratory/prediction_a_robust/` | Basic robustness check |
| Seed sensitivity N=68 | `outputs_exploratory/prediction_a_seed_sensitivity_n68/` | Limited seed variation |
| Seed sensitivity N=72 | `outputs_exploratory/prediction_a_seed_sensitivity_n72/` | Limited seed variation |

### A.3 Boundaries (explicitly stated)
- **NOT a proof** that 3+1 dimensions emerge from the axioms
- A2_full still selects low-dimensional families → the replacement is necessary
- Finite-size scaling NOT yet performed
- No analytic explanation for why consistency replacement favors 4D
- All N > 48 results are mixed (SIS-based), not exact

---

## SHARED INFRASTRUCTURE

| Item | Path |
|---|---|
| Frozen config files | `configs/` directory |
| Family generators | `generators.py` |
| Observable definitions | `observables.py`, `observables_geo.py` |
| Exact counting engine | `entropy_exact.py` |
| Benchmark (fault-tolerant) | `benchmark_exact_timing_safe.py` |
| Experiment runner | `experiment.py` |
| Code DOI | 10.5281/zenodo.18980657 |
| Paper B DOI (preprint) | 10.5281/zenodo.19048146 |
| Paper A DOI (preprint) | 10.5281/zenodo.19048324 |
| Paper C DOI (preprint) | 10.5281/zenodo.19048405 |
