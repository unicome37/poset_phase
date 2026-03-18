"""Prediction D: Sample-Level Stratified Permutation Test for CG Perturbation.

Upgrades the family-level Test B (6 families, underpowered) to a sample-level
continuous quasi-intervention test with 144 paired observations and stratified
permutation inference.

Design:
  X = Δpenalty_cg  (treatment: how much CG stability changed due to perturbation)
  Y = Δscore_local (outcome: how much the local score changed — pure, no CG term)

  These are mechanistically independent: penalty_cg is computed from CG operator
  output (drift, family switch, rank shift), while score_local = -β·log_H + γ·S/N
  is computed from the poset's entropy and action before CG. The perturbation
  (cover removal) changes the poset structure, which affects both independently.

  Strata = (N, family): 18 strata × 8 samples each = 144 paired observations.
  Within each stratum, permute X to break the causal link while preserving any
  stratum-level confounders.

  Statistic: weighted-mean Spearman ρ across strata (weighted by n_samples).

  If the causal claim is true (structural stability → better CG behavior → better
  ranking), then Δpenalty_cg should predict Δscore_local even after stratification.

  Additional tests:
  - Partial correlation controlling for family (pooled)
  - Within-family Spearman (more granular)
  - Dose-response: larger perturbation → stronger Δpenalty ~ Δscore coupling?
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

OUT_DIR = Path("outputs_exploratory/prediction_d_perturbation")
N_PERM = 100_000
SEED = 42


# ===================================================================
# Statistics
# ===================================================================

def _tieaware_rank(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    n = len(x)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1, dtype=float)
    xs = x[order]
    i = 0
    while i < n:
        j = i + 1
        while j < n and xs[j] == xs[i]:
            j += 1
        if j - i > 1:
            avg = 0.5 * ((i + 1) + j)
            ranks[order[i:j]] = avg
        i = j
    return ranks


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    xc = x - x.mean()
    yc = y - y.mean()
    denom = math.sqrt(float((xc * xc).sum() * (yc * yc).sum()))
    return float((xc * yc).sum() / denom) if denom > 1e-12 else 0.0


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(_tieaware_rank(x), _tieaware_rank(y))


# ===================================================================
# Main
# ===================================================================

def main():
    t0 = time.time()

    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║   Prediction D: Stratified Permutation Quasi-Intervention Test     ║")
    print("║   X = Δpenalty_cg,  Y = Δscore_local  (mechanistically independent)║")
    print("╚══════════════════════════════════════════════════════════════════════╝")

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    sample_df = pd.read_csv(OUT_DIR / "perturbation_sample_cg.csv")
    baseline = sample_df[sample_df["perturb"] == "p00"].copy()
    merge_cols = ["n", "family", "sample_id"]

    all_results = []

    for perturb_label in ["p05", "p10", "p20"]:
        frac_val = {"p05": 0.05, "p10": 0.10, "p20": 0.20}[perturb_label]
        pert = sample_df[sample_df["perturb"] == perturb_label].copy()

        merged = pert.merge(
            baseline[merge_cols + ["mean_penalty_cg", "mean_score_local_orig"]],
            on=merge_cols,
            suffixes=("_pert", "_orig"),
        )
        merged["delta_penalty"] = merged["mean_penalty_cg_pert"] - merged["mean_penalty_cg_orig"]
        merged["delta_score_local"] = merged["mean_score_local_orig_pert"] - merged["mean_score_local_orig_orig"]

        n_total = len(merged)
        print(f"\n{'='*70}")
        print(f"Perturbation: {perturb_label} ({frac_val:.0%} covers removed), n={n_total}")
        print(f"{'='*70}")

        # ==============================================================
        # Test 1: Pooled Spearman (no stratification)
        # ==============================================================
        x_all = merged["delta_penalty"].to_numpy(float)
        y_all = merged["delta_score_local"].to_numpy(float)
        rho_pooled = spearman(x_all, y_all)
        r_pooled = pearson(x_all, y_all)

        rng = np.random.default_rng(SEED)
        count_pooled = 1
        for _ in range(N_PERM):
            if abs(spearman(rng.permutation(x_all), y_all)) >= abs(rho_pooled):
                count_pooled += 1
        p_pooled = count_pooled / (N_PERM + 1)

        print(f"\n  [Test 1] Pooled (n={n_total}):")
        print(f"    Pearson  r(Δpenalty, Δscore_local)  = {r_pooled:+.4f}")
        print(f"    Spearman ρ(Δpenalty, Δscore_local)  = {rho_pooled:+.4f}")
        print(f"    Permutation p (|ρ|) = {p_pooled:.6f}  (n_perm={N_PERM})")

        # ==============================================================
        # Test 2: Stratified Permutation Test
        # ==============================================================
        strata = []
        for (n_val, fam), sub in merged.groupby(["n", "family"], sort=True):
            x = sub["delta_penalty"].to_numpy(float)
            y = sub["delta_score_local"].to_numpy(float)
            if len(x) >= 4:
                rx = _tieaware_rank(x)
                ry = _tieaware_rank(y)
                rx_c = rx - rx.mean()
                ry_c = ry - ry.mean()
                denom = math.sqrt(float((rx_c * rx_c).sum() * (ry_c * ry_c).sum()))
                denom = denom if denom > 1e-12 else 1.0
                strata.append({
                    "key": (n_val, fam),
                    "n": len(x),
                    "rx_c": rx_c,
                    "ry_c": ry_c,
                    "denom": denom,
                    "rho": float(np.dot(rx_c, ry_c) / denom),
                })

        # Weighted mean observed statistic
        total_n = sum(s["n"] for s in strata)
        obs_wmean = sum(s["rho"] * s["n"] for s in strata) / total_n if total_n else 0.0

        # Per-stratum details
        print(f"\n  [Test 2] Stratified by (N, family): {len(strata)} strata")
        for s in strata:
            print(f"    {str(s['key']):>45s}  n={s['n']}  ρ={s['rho']:+.4f}")
        print(f"    Weighted mean ρ = {obs_wmean:+.4f}")

        # Stratified permutation
        rng2 = np.random.default_rng(SEED + 1)
        count_strat = 1
        for _ in range(N_PERM):
            perm_wmean = 0.0
            for s in strata:
                perm_idx = rng2.permutation(len(s["rx_c"]))
                perm_rho = float(np.dot(s["rx_c"][perm_idx], s["ry_c"]) / s["denom"])
                perm_wmean += perm_rho * s["n"]
            perm_wmean /= total_n
            if abs(perm_wmean) >= abs(obs_wmean):
                count_strat += 1
        p_strat = count_strat / (N_PERM + 1)

        print(f"    Stratified permutation p (|ρ̄|) = {p_strat:.6f}  (n_perm={N_PERM})")

        # ==============================================================
        # Test 3: Per-N Spearman (within-N pooled across families)
        # ==============================================================
        print(f"\n  [Test 3] Per-N pooled:")
        for n_val in sorted(merged["n"].unique()):
            sub = merged[merged["n"] == n_val]
            x_n = sub["delta_penalty"].to_numpy(float)
            y_n = sub["delta_score_local"].to_numpy(float)
            rho_n = spearman(x_n, y_n)
            r_n = pearson(x_n, y_n)

            rng3 = np.random.default_rng(SEED + n_val)
            count_n = 1
            n_perm_n = 50_000
            for _ in range(n_perm_n):
                if abs(spearman(rng3.permutation(x_n), y_n)) >= abs(rho_n):
                    count_n += 1
            p_n = count_n / (n_perm_n + 1)

            sig = "***" if p_n < 0.001 else ("**" if p_n < 0.01 else ("*" if p_n < 0.05 else "ns"))
            print(f"    N={n_val}: r={r_n:+.4f}, ρ={rho_n:+.4f}, p={p_n:.5f} {sig}  (n={len(sub)})")

        # ==============================================================
        # Test 4: Within-family Spearman (8 samples per family per N)
        # ==============================================================
        print(f"\n  [Test 4] Within-family (most granular):")
        within_rhos = []
        for (n_val, fam), sub in merged.groupby(["n", "family"], sort=True):
            x_wf = sub["delta_penalty"].to_numpy(float)
            y_wf = sub["delta_score_local"].to_numpy(float)
            if len(x_wf) >= 4:
                rho_wf = spearman(x_wf, y_wf)
                within_rhos.append(rho_wf)
                print(f"    N={n_val} {fam:30s}: ρ={rho_wf:+.4f}  (n={len(x_wf)})")
        mean_within = np.mean(within_rhos) if within_rhos else float("nan")
        pos_frac = sum(1 for r in within_rhos if r > 0) / len(within_rhos) if within_rhos else 0
        print(f"    Mean within-family ρ = {mean_within:+.4f}")
        print(f"    Fraction positive: {pos_frac:.2f} ({sum(1 for r in within_rhos if r > 0)}/{len(within_rhos)})")

        # ==============================================================
        # Test 5: Δpenalty component decomposition
        # ==============================================================
        # Which component of penalty_cg is driving the coupling?
        # penalty_cg = 10*drift + 1.5*switch + 0.5*rank_shift
        merged["delta_drift"] = pert.merge(
            baseline[merge_cols + ["mean_self_drift"]], on=merge_cols,
            suffixes=("_pert", "_orig"),
        )["mean_self_drift_pert"] - pert.merge(
            baseline[merge_cols + ["mean_self_drift"]], on=merge_cols,
            suffixes=("_pert", "_orig"),
        )["mean_self_drift_orig"]

        merged["delta_switch"] = pert.merge(
            baseline[merge_cols + ["mean_fam_switch"]], on=merge_cols,
            suffixes=("_pert", "_orig"),
        )["mean_fam_switch_pert"] - pert.merge(
            baseline[merge_cols + ["mean_fam_switch"]], on=merge_cols,
            suffixes=("_pert", "_orig"),
        )["mean_fam_switch_orig"]

        print(f"\n  [Test 5] Component decomposition:")
        for comp_name, comp_col in [("Δdrift", "delta_drift"), ("Δswitch", "delta_switch")]:
            c = merged[comp_col].to_numpy(float)
            rho_c = spearman(c, y_all)
            r_c = pearson(c, y_all)
            print(f"    {comp_name:>10s} ~ Δscore_local:  r={r_c:+.4f}, ρ={rho_c:+.4f}")

        # Save results
        all_results.append({
            "perturb": perturb_label,
            "perturb_frac": frac_val,
            "n_total": n_total,
            "n_strata": len(strata),
            "rho_pooled": rho_pooled,
            "r_pooled": r_pooled,
            "p_pooled": p_pooled,
            "rho_stratified_wmean": obs_wmean,
            "p_stratified": p_strat,
            "mean_within_family_rho": mean_within,
            "frac_positive_within": pos_frac,
        })

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT_DIR / "stratified_intervention_results.csv", index=False, encoding="utf-8-sig")

    print("\n" + "=" * 70)
    print("SUMMARY: Sample-Level Stratified Quasi-Intervention Test")
    print("X = Δpenalty_cg,  Y = Δscore_local  (mechanistically independent)")
    print("=" * 70)
    print(f"{'pert':>6s}  {'n':>4s}  {'ρ_pool':>8s}  {'p_pool':>10s}  {'ρ̄_strat':>8s}  {'p_strat':>10s}  {'within_ρ̄':>8s}  {'%pos':>5s}")
    for _, r in results_df.iterrows():
        sig_p = "***" if r["p_pooled"] < 0.001 else ("**" if r["p_pooled"] < 0.01 else ("*" if r["p_pooled"] < 0.05 else "ns"))
        sig_s = "***" if r["p_stratified"] < 0.001 else ("**" if r["p_stratified"] < 0.01 else ("*" if r["p_stratified"] < 0.05 else "ns"))
        print(f"  {r['perturb']:>4s}  {int(r['n_total']):>4d}  {r['rho_pooled']:>+8.4f}  {r['p_pooled']:>9.6f}{sig_p:>3s}  "
              f"{r['rho_stratified_wmean']:>+8.4f}  {r['p_stratified']:>9.6f}{sig_s:>3s}  "
              f"{r['mean_within_family_rho']:>+8.4f}  {r['frac_positive_within']:>5.0%}")

    elapsed = time.time() - t0
    print(f"\nRuntime: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
