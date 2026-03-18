"""Prediction D: High-Power Stratified Permutation Test (32 samples/family).

Re-runs the CG perturbation experiment with 4× more samples per stratum
to determine whether the within-stratum null (ρ̄ ≈ -0.06) seen at n=8
was real or a power artefact.

Previous result (n=8 per stratum):
  Pooled: p05 ρ=-0.21*, p20 ρ=-0.31***
  Stratified: all ns (p>0.11)

This test: n=32 per stratum → 576 total paired observations.
If true within-stratum ρ ≈ 0.3, power ~60% per stratum, and the
aggregated 18-stratum test should reach ~95%+ power.

Design unchanged from §6.5:
  X = Δpenalty_cg, Y = Δscore_local (mechanistically independent)
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Import core functions from existing scripts
from prediction_d_perturbation_intervention import (
    hasse_covers,
    transitive_closure,
    perturb_remove_covers,
    run_cg_pipeline_for_posets,
    FAMILY_NAMES,
    N_VALUES,
    KEEP_RATIO,
    GAMMA,
    BETA,
    ACTION_MODE,
    CG_REPEATS,
    SIS_RUNS,
    ZETA_VALUES,
    SEED_OFFSET,
)
from experiment import FAMILIES
from stability import (
    signature_dict,
    family_centroids,
    SIGNATURE_COLUMNS,
)
from runtime_utils import estimate_entropy
from action import action_value, get_action_penalty
from generators import Poset

# ===================================================================
# Configuration
# ===================================================================
SAMPLES_PER_FAMILY = 32  # 4× increase from 8
PERTURB_FRACS = (0.0, 0.05, 0.10, 0.20)
OUT_DIR = Path("outputs_exploratory/prediction_d_perturbation_n32")
N_PERM = 100_000
SEED = 42


# ===================================================================
# Statistics (copied from stratified test for self-containment)
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
# Stratified permutation test
# ===================================================================

def run_stratified_test(merged: pd.DataFrame, perturb_label: str, frac_val: float):
    """Run all stratified tests on a merged baseline-vs-perturbed dataframe."""
    n_total = len(merged)
    x_all = merged["delta_penalty"].to_numpy(float)
    y_all = merged["delta_score_local"].to_numpy(float)

    print(f"\n{'='*70}")
    print(f"Perturbation: {perturb_label} ({frac_val:.0%}), n={n_total}")
    print(f"{'='*70}")

    # Test 1: Pooled
    rho_pooled = spearman(x_all, y_all)
    r_pooled = pearson(x_all, y_all)
    rng = np.random.default_rng(SEED)
    count_p = 1
    for _ in range(N_PERM):
        if abs(spearman(rng.permutation(x_all), y_all)) >= abs(rho_pooled):
            count_p += 1
    p_pooled = count_p / (N_PERM + 1)
    sig_p = "***" if p_pooled < 0.001 else ("**" if p_pooled < 0.01 else ("*" if p_pooled < 0.05 else "ns"))
    print(f"\n  [Pooled] r={r_pooled:+.4f}, ρ={rho_pooled:+.4f}, p={p_pooled:.6f} {sig_p}")

    # Test 2: Stratified
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

    total_n = sum(s["n"] for s in strata)
    obs_wmean = sum(s["rho"] * s["n"] for s in strata) / total_n if total_n else 0.0

    print(f"\n  [Stratified] {len(strata)} strata, weighted-mean ρ = {obs_wmean:+.4f}")
    for s in strata:
        print(f"    {str(s['key']):>45s}  n={s['n']}  ρ={s['rho']:+.4f}")

    rng2 = np.random.default_rng(SEED + 1)
    count_s = 1
    for _ in range(N_PERM):
        perm_wmean = 0.0
        for s in strata:
            perm_idx = rng2.permutation(len(s["rx_c"]))
            perm_rho = float(np.dot(s["rx_c"][perm_idx], s["ry_c"]) / s["denom"])
            perm_wmean += perm_rho * s["n"]
        perm_wmean /= total_n
        if abs(perm_wmean) >= abs(obs_wmean):
            count_s += 1
    p_strat = count_s / (N_PERM + 1)
    sig_s = "***" if p_strat < 0.001 else ("**" if p_strat < 0.01 else ("*" if p_strat < 0.05 else "ns"))
    print(f"    Stratified permutation p = {p_strat:.6f} {sig_s}")

    # Test 3: Per-N pooled
    print(f"\n  [Per-N]")
    for n_val in sorted(merged["n"].unique()):
        sub = merged[merged["n"] == n_val]
        x_n = sub["delta_penalty"].to_numpy(float)
        y_n = sub["delta_score_local"].to_numpy(float)
        rho_n = spearman(x_n, y_n)
        rng3 = np.random.default_rng(SEED + n_val)
        count_n = 1
        n_perm_n = 20_000
        for _ in range(n_perm_n):
            if abs(spearman(rng3.permutation(x_n), y_n)) >= abs(rho_n):
                count_n += 1
        p_n = count_n / (n_perm_n + 1)
        sig_n = "***" if p_n < 0.001 else ("**" if p_n < 0.01 else ("*" if p_n < 0.05 else "ns"))
        print(f"    N={n_val}: ρ={rho_n:+.4f}, p={p_n:.5f} {sig_n}  (n={len(sub)})")

    # Test 4: Within-family
    within_rhos = []
    for (n_val, fam), sub in merged.groupby(["n", "family"], sort=True):
        x_wf = sub["delta_penalty"].to_numpy(float)
        y_wf = sub["delta_score_local"].to_numpy(float)
        if len(x_wf) >= 4:
            within_rhos.append(spearman(x_wf, y_wf))
    mean_within = np.mean(within_rhos) if within_rhos else float("nan")
    neg_frac = sum(1 for r in within_rhos if r < 0) / len(within_rhos) if within_rhos else 0
    print(f"\n  [Within-family] mean ρ = {mean_within:+.4f}, negative fraction: {neg_frac:.0%} ({sum(1 for r in within_rhos if r < 0)}/{len(within_rhos)})")

    return {
        "perturb": perturb_label,
        "frac": frac_val,
        "n_total": n_total,
        "n_per_stratum": SAMPLES_PER_FAMILY,
        "rho_pooled": rho_pooled,
        "p_pooled": p_pooled,
        "rho_strat_wmean": obs_wmean,
        "p_strat": p_strat,
        "mean_within_rho": mean_within,
        "neg_fraction": neg_frac,
    }


# ===================================================================
# Main
# ===================================================================

def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("╔════════════════════════════════════════════════════════════════════╗")
    print("║  Prediction D: High-Power Stratified Test (32 samples/family)   ║")
    print("║  X = Δpenalty_cg,  Y = Δscore_local  (independent variables)    ║")
    print("╚════════════════════════════════════════════════════════════════════╝")

    seed_base = SEED_OFFSET * 10_000_000

    # ------------------------------------------------------------------
    # Step 1: Generate posets (32 per family)
    # ------------------------------------------------------------------
    print(f"\nStep 1: Generating {SAMPLES_PER_FAMILY} posets per family...", flush=True)
    posets_orig: dict[tuple[int, str, int], Poset] = {}
    sig_rows = []
    for n in N_VALUES:
        for family in FAMILY_NAMES:
            generator = FAMILIES[family]
            for sid in range(SAMPLES_PER_FAMILY):
                seed = seed_base + 1000 * n + sid
                poset = generator(n=n, seed=seed)
                posets_orig[(n, family, sid)] = poset
                sig = signature_dict(poset)
                sig["n"] = n
                sig["family"] = family
                sig["sample_id"] = sid
                sig_rows.append(sig)

    sig_df = pd.DataFrame(sig_rows)
    centroids = family_centroids(sig_df[["n", "family"] + SIGNATURE_COLUMNS])
    print(f"  Generated {len(posets_orig)} posets. Time: {time.time() - t0:.1f}s")

    # Base family ranks
    base_ranks: dict[tuple[int, str], int] = {}
    for n in N_VALUES:
        fam_scores = {}
        for family in FAMILY_NAMES:
            scores = []
            for sid in range(SAMPLES_PER_FAMILY):
                poset = posets_orig[(n, family, sid)]
                log_h, _ = estimate_entropy(poset, sis_runs=SIS_RUNS, seed=n * 1000 + sid, exact_threshold=0)
                pen = get_action_penalty(poset, ACTION_MODE)
                sc = action_value(log_h, pen, BETA, GAMMA)
                scores.append(sc)
            fam_scores[family] = np.mean(scores)
        sorted_fams = sorted(fam_scores.items(), key=lambda x: x[1])
        for rank_i, (fam, _) in enumerate(sorted_fams, 1):
            base_ranks[(n, fam)] = rank_i

    print(f"  Base ranks computed. Time: {time.time() - t0:.1f}s")

    # ------------------------------------------------------------------
    # Step 2: Perturbation + CG pipeline
    # ------------------------------------------------------------------
    all_sample_frames = []

    for perturb_frac in PERTURB_FRACS:
        label = f"p{int(perturb_frac * 100):02d}"
        print(f"\nStep 2-{label}: {perturb_frac:.0%} cover removal...", flush=True)

        posets_pert: dict[tuple[int, str, int], Poset] = {}
        for (n, family, sid), poset_orig in sorted(posets_orig.items()):
            if perturb_frac == 0.0:
                posets_pert[(n, family, sid)] = poset_orig
            else:
                closure = poset_orig.closure.copy()
                new_closure = perturb_remove_covers(
                    closure, perturb_frac,
                    seed=SEED_OFFSET + n * 10000 + sid * 100 + int(perturb_frac * 100),
                )
                posets_pert[(n, family, sid)] = Poset(new_closure)

        _, sample_cg = run_cg_pipeline_for_posets(
            posets_pert, centroids, base_ranks,
            keep_ratio=KEEP_RATIO,
            cg_repeats=CG_REPEATS,
            sis_runs=SIS_RUNS,
            gamma=GAMMA,
            beta=BETA,
            action_mode=ACTION_MODE,
            zeta_values=ZETA_VALUES,
            perturb_label=label,
        )
        all_sample_frames.append(sample_cg)
        print(f"  Done {label}. Time: {time.time() - t0:.1f}s")

    # Merge all sample data
    sample_df = pd.concat(all_sample_frames, ignore_index=True)
    sample_df.to_csv(OUT_DIR / "perturbation_sample_cg_n32.csv", index=False, encoding="utf-8-sig")
    print(f"\nSaved {len(sample_df)} rows to perturbation_sample_cg_n32.csv")

    # ------------------------------------------------------------------
    # Step 3: Stratified permutation test
    # ------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("Step 3: Stratified Permutation Tests (high-power, n=32/stratum)")
    print(f"{'='*70}")

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

        result = run_stratified_test(merged, perturb_label, frac_val)
        all_results.append(result)

    # ------------------------------------------------------------------
    # Summary : n=8 vs n=32 comparison
    # ------------------------------------------------------------------
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT_DIR / "stratified_highpower_results.csv", index=False, encoding="utf-8-sig")

    print(f"\n{'='*70}")
    print("SUMMARY: n=32 vs n=8 per stratum")
    print(f"{'='*70}")
    print(f"{'pert':>6s}  {'n':>5s}  {'ρ_pool':>8s}  {'p_pool':>10s}  {'ρ̄_strat':>8s}  {'p_strat':>10s}  {'within_ρ̄':>8s}  {'%neg':>5s}")
    for _, r in results_df.iterrows():
        sig_p = "***" if r["p_pooled"] < 0.001 else ("**" if r["p_pooled"] < 0.01 else ("*" if r["p_pooled"] < 0.05 else "ns"))
        sig_s = "***" if r["p_strat"] < 0.001 else ("**" if r["p_strat"] < 0.01 else ("*" if r["p_strat"] < 0.05 else "ns"))
        print(f"  {r['perturb']:>4s}  {int(r['n_total']):>5d}  {r['rho_pooled']:>+8.4f}  {r['p_pooled']:>9.6f}{sig_p:>3s}  "
              f"{r['rho_strat_wmean']:>+8.4f}  {r['p_strat']:>9.6f}{sig_s:>3s}  "
              f"{r['mean_within_rho']:>+8.4f}  {r['neg_fraction']:>5.0%}")

    print(f"\nPrevious results (n=8/stratum) for comparison:")
    print(f"  p05: ρ_pool=-0.21 p=0.011*, ρ̄_strat=-0.14 p=0.113 ns")
    print(f"  p10: ρ_pool=-0.16 p=0.062 ns, ρ̄_strat=-0.05 p=0.613 ns")
    print(f"  p20: ρ_pool=-0.31 p=0.0002***, ρ̄_strat=-0.06 p=0.497 ns")

    elapsed = time.time() - t0
    print(f"\nRuntime: {elapsed:.1f}s ({elapsed / 60:.1f} min)")


if __name__ == "__main__":
    main()
