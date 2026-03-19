"""Prediction D Independence Verification

Goal: Prove that W1 (Wasserstein spectral distance) and ΔH_int are
independent of F7 — i.e., they carry NEW information beyond what F7 already
captures about dynamic stability.

Three tests:
  1. Partial correlation: P_basin ~ F7 + D_spec (does D_spec remain significant?)
  2. Matched-on-F7: stratify by F7 quartile, compare D_spec within strata
  3. Within-family: within each family, does low D_spec → higher P_basin?

Uses data from _d_recovery_experiment.py (outputs_d_recovery/).
Expands to N=16,20,28,36 for larger sample size.

Usage:
    python _d_independence_verification.py
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_kr_like,
    transitive_closure,
)
from unified_functional import (
    compute_log_H,
    compute_pi_geo,
    compute_sigma_hist,
    compute_xi_dim,
)
from observables import comparable_fraction
from coarse_grain import coarse_grain_delete_nodes
from prediction_a_mid_gap import interval_spectrum
from prediction_a_bd_bridge import count_intervals_fast

# ══════════════════════════════════════════════════════════════════════════
# F7 computation (same as _d_recovery_experiment.py)
# ══════════════════════════════════════════════════════════════════════════

def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_R(poset: Poset) -> float:
    counts = count_intervals_fast(poset)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    C0 = counts.get(0, 0)
    return 1.0 - C0 / total


def compute_F7(poset, alpha0=16.0, q=-0.5, lam=10.0, eta=0.6,
               Rc=0.25, w=0.015, N0=20.0, n_sis=64):
    N = poset.n
    log_H = compute_log_H(poset, n_runs=n_sis)
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim, _ = compute_xi_dim(poset)
    R = compute_R(poset)
    alpha_N = alpha0 * (N0 / N) ** abs(q)
    wall = alpha_N * sigmoid((R - Rc) / w)
    return log_H + 0.0004 * pi_geo - lam * sigma_hist + eta * xi_dim + wall


def compute_F7_fast(poset, **kw):
    return compute_F7(poset, n_sis=32, **kw)


# ══════════════════════════════════════════════════════════════════════════
# Hasse diagram + perturbation + dynamics (from _d_recovery_experiment.py)
# ══════════════════════════════════════════════════════════════════════════

def hasse_diagram(closure):
    c = closure.astype(np.uint8)
    has_intermediate = (c @ c).astype(bool)
    cover = closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return cover


def perturb_add_edge(poset, rng):
    closure = poset.closure
    comparable = closure | closure.T
    np.fill_diagonal(comparable, True)
    incomparable = ~comparable
    candidates = np.argwhere(np.triu(incomparable, k=1))
    if len(candidates) == 0:
        return None
    idx = rng.integers(len(candidates))
    i, j = candidates[idx]
    new_closure = closure.copy()
    if rng.random() < 0.5:
        new_closure[i, j] = True
    else:
        new_closure[j, i] = True
    new_closure = transitive_closure(new_closure)
    if np.any(new_closure & new_closure.T):
        return None
    return Poset(new_closure)


def perturb_delete_edge(poset, rng):
    cover = hasse_diagram(poset.closure)
    edges = np.argwhere(cover)
    if len(edges) == 0:
        return None
    idx = rng.integers(len(edges))
    i, j = edges[idx]
    new_cover = cover.copy()
    new_cover[i, j] = False
    new_closure = transitive_closure(new_cover)
    return Poset(new_closure)


def perturb_delete_node(poset, rng):
    if poset.n < 6:
        return None
    node = rng.integers(poset.n)
    keep = [i for i in range(poset.n) if i != node]
    return Poset(poset.closure[np.ix_(keep, keep)])


PERTURBATIONS = {
    "add_edge": perturb_add_edge,
    "del_edge": perturb_delete_edge,
    "del_node": perturb_delete_node,
}


def propose_swap_move(closure, rng):
    cover = hasse_diagram(closure)
    edges = np.argwhere(cover)
    if len(edges) == 0:
        return None
    rm_idx = rng.integers(len(edges))
    ri, rj = edges[rm_idx]
    new_cover = cover.copy()
    new_cover[ri, rj] = False
    intermediate = transitive_closure(new_cover)
    comparable = intermediate | intermediate.T
    np.fill_diagonal(comparable, True)
    incomparable = ~comparable
    candidates = np.argwhere(np.triu(incomparable, k=1))
    if len(candidates) == 0:
        return None
    add_idx = rng.integers(len(candidates))
    ai, aj = candidates[add_idx]
    if rng.random() < 0.5:
        intermediate[ai, aj] = True
    else:
        intermediate[aj, ai] = True
    new_closure = transitive_closure(intermediate)
    if np.any(new_closure & new_closure.T):
        return None
    return new_closure


def recovery_dynamics(perturbed, original_cf, n_steps=150, temperature=1.0, seed=42):
    rng = np.random.default_rng(seed)
    current = perturbed
    current_F = compute_F7_fast(current)
    cf_now = comparable_fraction(current)
    returned_basin = False
    return_step_basin = -1

    for step in range(1, n_steps + 1):
        new_closure = propose_swap_move(current.closure, rng)
        if new_closure is None:
            continue
        new_poset = Poset(new_closure)
        new_F = compute_F7_fast(new_poset)
        delta_F = new_F - current_F
        if delta_F < 0 or rng.random() < math.exp(-delta_F / max(temperature, 1e-10)):
            current = new_poset
            current_F = new_F

        cf_now = comparable_fraction(current)
        if not returned_basin and abs(cf_now - original_cf) < 0.05:
            returned_basin = True
            return_step_basin = step

    return {
        "basin_retention": abs(cf_now - original_cf) < 0.05,
        "returned_basin": returned_basin,
        "return_step_basin": return_step_basin,
        "final_F": current_F,
        "final_cf": cf_now,
    }


# ══════════════════════════════════════════════════════════════════════════
# Spectral fidelity computation
# ══════════════════════════════════════════════════════════════════════════

def compute_spectral_vars(poset, keep_ratio=0.7, n_cg=3, seed=42):
    spec_orig = interval_spectrum(poset)
    H_orig = spec_orig["spectral_entropy"]

    counts_orig = count_intervals_fast(poset)
    total_orig = sum(counts_orig.values())
    if total_orig == 0:
        return {"H_int_orig": 0.0, "delta_H_int": 0.0, "wasserstein": 0.0}

    max_k_orig = max(counts_orig.keys()) if counts_orig else 0

    H_cg_list, wass_list = [], []
    for i in range(n_cg):
        cg = coarse_grain_delete_nodes(poset, keep_ratio=keep_ratio, seed=seed + i)
        if cg.n < 4:
            continue
        spec_cg = interval_spectrum(cg)
        H_cg_list.append(spec_cg["spectral_entropy"])

        counts_cg = count_intervals_fast(cg)
        total_cg = sum(counts_cg.values())
        if total_cg == 0:
            continue
        max_k_cg = max(counts_cg.keys()) if counts_cg else 0
        max_k = max(max_k_orig, max_k_cg)
        dist_o = np.array([counts_orig.get(k, 0) / total_orig for k in range(max_k + 1)])
        dist_c = np.array([counts_cg.get(k, 0) / total_cg for k in range(max_k + 1)])
        cdf_o = np.cumsum(dist_o)
        cdf_c = np.cumsum(dist_c)
        wass_list.append(float(np.sum(np.abs(cdf_o - cdf_c))))

    H_cg_mean = float(np.mean(H_cg_list)) if H_cg_list else H_orig
    wass_mean = float(np.mean(wass_list)) if wass_list else 0.0

    return {
        "H_int_orig": H_orig,
        "delta_H_int": abs(H_cg_mean - H_orig),
        "wasserstein": wass_mean,
    }


# ══════════════════════════════════════════════════════════════════════════
# Main experiment
# ══════════════════════════════════════════════════════════════════════════

FAMILIES = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR": generate_kr_like,
}


def main():
    print("=" * 90)
    print("PREDICTION D: INDEPENDENCE VERIFICATION")
    print("Does W1 / ΔH_int carry information beyond F7?")
    print("=" * 90)

    N_VALUES = (16, 20, 28, 36)
    REPS = 5
    PERT_REPS = 2
    seed_base = 42

    # ── Phase 1: Generate all samples + compute F7 & spectral vars ──
    print(f"\n[1/4] Generating samples (N={N_VALUES}, {REPS} reps, 5 families)...")
    samples = []
    for N in N_VALUES:
        for fam_name, gen_func in FAMILIES.items():
            for rep in range(REPS):
                s = seed_base + rep * 1000 + N * 7
                poset = gen_func(N, seed=s)
                F7 = compute_F7(poset)
                cf = comparable_fraction(poset)
                sv = compute_spectral_vars(poset, seed=s)
                samples.append({
                    "family": fam_name, "N": N, "rep": rep,
                    "poset": poset, "F7": F7, "comp_frac": cf,
                    **sv,
                })

        done_n = sum(1 for s in samples if s["N"] == N)
        print(f"  N={N}: {done_n} samples generated")

    print(f"  Total: {len(samples)} samples")

    # ── Phase 2: Perturbation + recovery for P_basin ──
    print(f"\n[2/4] Perturbation + recovery ({len(samples)} × 3 perts × {PERT_REPS} reps)...")
    recovery_map = {}  # key=(family,N,rep) -> list of basin_retention bools

    trial_count = 0
    for si, sample in enumerate(samples):
        key = (sample["family"], sample["N"], sample["rep"])
        recovery_map[key] = []

        for pert_name, pert_func in PERTURBATIONS.items():
            for pr in range(PERT_REPS):
                pseed = abs(hash((key, pert_name, pr))) % (2**31)
                pert_rng = np.random.default_rng(pseed)
                perturbed = pert_func(sample["poset"], pert_rng)
                if perturbed is None:
                    continue

                rec = recovery_dynamics(
                    perturbed, sample["comp_frac"],
                    n_steps=150, temperature=1.0,
                    seed=(pseed + 7) % (2**31),
                )
                recovery_map[key].append(rec["basin_retention"])
                trial_count += 1

        if (si + 1) % 20 == 0:
            print(f"  {si+1}/{len(samples)} samples processed ({trial_count} trials)")

    print(f"  Total recovery trials: {trial_count}")

    # ── Merge into analysis table ──
    data = []
    for sample in samples:
        key = (sample["family"], sample["N"], sample["rep"])
        brs = recovery_map.get(key, [])
        if not brs:
            continue
        data.append({
            "family": sample["family"],
            "N": sample["N"],
            "rep": sample["rep"],
            "F7": sample["F7"],
            "comp_frac": sample["comp_frac"],
            "wasserstein": sample["wasserstein"],
            "delta_H_int": sample["delta_H_int"],
            "P_basin": float(np.mean(brs)),
            "n_trials": len(brs),
        })

    print(f"  Analysis table: {len(data)} rows")

    F7_arr = np.array([d["F7"] for d in data])
    W1_arr = np.array([d["wasserstein"] for d in data])
    dH_arr = np.array([d["delta_H_int"] for d in data])
    Pb_arr = np.array([d["P_basin"] for d in data])
    fam_arr = np.array([d["family"] for d in data])
    N_arr = np.array([d["N"] for d in data])

    # ══════════════════════════════════════════════════════════════════
    # TEST 1: Partial correlation (residual regression)
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 90)
    print("TEST 1: Partial correlation — P_basin ~ F7 + D_spec")
    print("=" * 90)

    def partial_corr(x, y, z):
        """Partial Spearman correlation of x and y controlling z."""
        rho_xz, _ = sp_stats.spearmanr(x, z)
        rho_yz, _ = sp_stats.spearmanr(y, z)
        res_x = sp_stats.rankdata(x) - rho_xz * sp_stats.rankdata(z)
        res_y = sp_stats.rankdata(y) - rho_yz * sp_stats.rankdata(z)
        rho_partial, p_partial = sp_stats.spearmanr(res_x, res_y)
        return rho_partial, p_partial

    # Raw correlations
    rho_F7_Pb, p_F7_Pb = sp_stats.spearmanr(F7_arr, Pb_arr)
    rho_W1_Pb, p_W1_Pb = sp_stats.spearmanr(W1_arr, Pb_arr)
    rho_dH_Pb, p_dH_Pb = sp_stats.spearmanr(dH_arr, Pb_arr)

    print(f"\n  Raw correlations with P_basin:")
    print(f"    F7:           ρ = {rho_F7_Pb:+.3f}, p = {p_F7_Pb:.4f}")
    print(f"    W1:           ρ = {rho_W1_Pb:+.3f}, p = {p_W1_Pb:.4f}")
    print(f"    ΔH_int:       ρ = {rho_dH_Pb:+.3f}, p = {p_dH_Pb:.4f}")

    # Partial correlations controlling F7
    rho_W1_partial, p_W1_partial = partial_corr(W1_arr, Pb_arr, F7_arr)
    rho_dH_partial, p_dH_partial = partial_corr(dH_arr, Pb_arr, F7_arr)

    print(f"\n  Partial correlations with P_basin (controlling F7):")
    print(f"    W1 | F7:      ρ = {rho_W1_partial:+.3f}, p = {p_W1_partial:.4f} "
          f"{'★★★' if p_W1_partial < 0.001 else ('★★' if p_W1_partial < 0.01 else ('★' if p_W1_partial < 0.05 else 'ns'))}")
    print(f"    ΔH_int | F7:  ρ = {rho_dH_partial:+.3f}, p = {p_dH_partial:.4f} "
          f"{'★★★' if p_dH_partial < 0.001 else ('★★' if p_dH_partial < 0.01 else ('★' if p_dH_partial < 0.05 else 'ns'))}")

    # Also control F7 + N jointly
    # Residualize on [F7, N] via rank regression
    def partial_corr_2(x, y, z1, z2):
        """Partial Spearman correlation of x,y controlling z1 AND z2."""
        # Residualize x on (z1, z2)
        from numpy.linalg import lstsq
        Z = np.column_stack([sp_stats.rankdata(z1), sp_stats.rankdata(z2)])
        rx = sp_stats.rankdata(x)
        ry = sp_stats.rankdata(y)
        coef_x, _, _, _ = lstsq(Z, rx, rcond=None)
        coef_y, _, _, _ = lstsq(Z, ry, rcond=None)
        res_x = rx - Z @ coef_x
        res_y = ry - Z @ coef_y
        return sp_stats.spearmanr(res_x, res_y)

    rho_W1_p2, p_W1_p2 = partial_corr_2(W1_arr, Pb_arr, F7_arr, N_arr)
    rho_dH_p2, p_dH_p2 = partial_corr_2(dH_arr, Pb_arr, F7_arr, N_arr)

    print(f"\n  Partial correlations with P_basin (controlling F7 + N):")
    print(f"    W1 | F7,N:    ρ = {rho_W1_p2:+.3f}, p = {p_W1_p2:.4f} "
          f"{'★★★' if p_W1_p2 < 0.001 else ('★★' if p_W1_p2 < 0.01 else ('★' if p_W1_p2 < 0.05 else 'ns'))}")
    print(f"    ΔH | F7,N:    ρ = {rho_dH_p2:+.3f}, p = {p_dH_p2:.4f} "
          f"{'★★★' if p_dH_p2 < 0.001 else ('★★' if p_dH_p2 < 0.01 else ('★' if p_dH_p2 < 0.05 else 'ns'))}")

    # ══════════════════════════════════════════════════════════════════
    # TEST 2: Matched-on-F7 stratified comparison
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 90)
    print("TEST 2: Matched-on-F7 — stratify by F7 quartile")
    print("=" * 90)

    quartiles = np.percentile(F7_arr, [25, 50, 75])
    q_labels = np.digitize(F7_arr, quartiles)  # 0,1,2,3

    print(f"\n  F7 quartile boundaries: {quartiles}")
    print(f"  {'Q':>3} {'n':>4} {'mean_F7':>8} {'mean_W1':>8} {'mean_dH':>8} "
          f"{'mean_Pb':>8} {'ρ(W1,Pb)':>10} {'ρ(dH,Pb)':>10}")

    for q in range(4):
        mask = q_labels == q
        if mask.sum() < 4:
            continue
        mF7 = np.mean(F7_arr[mask])
        mW1 = np.mean(W1_arr[mask])
        mdH = np.mean(dH_arr[mask])
        mPb = np.mean(Pb_arr[mask])

        if mask.sum() >= 5 and np.std(W1_arr[mask]) > 1e-10 and np.std(Pb_arr[mask]) > 1e-10:
            rho_w, _ = sp_stats.spearmanr(W1_arr[mask], Pb_arr[mask])
        else:
            rho_w = float('nan')
        if mask.sum() >= 5 and np.std(dH_arr[mask]) > 1e-10 and np.std(Pb_arr[mask]) > 1e-10:
            rho_d, _ = sp_stats.spearmanr(dH_arr[mask], Pb_arr[mask])
        else:
            rho_d = float('nan')

        print(f"  Q{q}: {mask.sum():4d} {mF7:8.2f} {mW1:8.4f} {mdH:8.3f} "
              f"{mPb:8.1%} {rho_w:+10.3f} {rho_d:+10.3f}")

    # Median split within each quartile
    print(f"\n  Median-split within F7 quartiles:")
    print(f"  {'Q':>3} {'var':>6} {'P_basin(low)':>12} {'P_basin(high)':>13} {'diff':>8} {'direction':>10}")

    for q in range(4):
        mask = q_labels == q
        if mask.sum() < 4:
            continue
        for var_name, var_arr in [("W1", W1_arr), ("ΔH", dH_arr)]:
            med = np.median(var_arr[mask])
            low = Pb_arr[mask & (var_arr <= med)]
            high = Pb_arr[mask & (var_arr > med)]
            if len(low) == 0 or len(high) == 0:
                continue
            ml = np.mean(low)
            mh = np.mean(high)
            direction = "low→stable" if ml > mh else "high→stable"
            print(f"  Q{q}: {var_name:>6} {ml:12.1%} {mh:13.1%} {ml-mh:+8.1%} {direction:>10}")

    # ══════════════════════════════════════════════════════════════════
    # TEST 3: Within-family tests
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 90)
    print("TEST 3: Within-family — does low D_spec → higher P_basin?")
    print("=" * 90)

    print(f"\n  {'family':6s} {'n':>4} {'ρ(W1,Pb)':>10} {'p':>8} {'ρ(dH,Pb)':>10} {'p':>8}")

    fam_rho_W1 = []
    fam_rho_dH = []
    for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR"]:
        mask = fam_arr == fam
        n_fam = mask.sum()
        if n_fam < 5:
            continue
        W1_f = W1_arr[mask]
        dH_f = dH_arr[mask]
        Pb_f = Pb_arr[mask]

        if np.std(W1_f) > 1e-10 and np.std(Pb_f) > 1e-10:
            rho_w, p_w = sp_stats.spearmanr(W1_f, Pb_f)
        else:
            rho_w, p_w = float('nan'), 1.0
        if np.std(dH_f) > 1e-10 and np.std(Pb_f) > 1e-10:
            rho_d, p_d = sp_stats.spearmanr(dH_f, Pb_f)
        else:
            rho_d, p_d = float('nan'), 1.0

        fam_rho_W1.append(rho_w)
        fam_rho_dH.append(rho_d)
        print(f"  {fam:6s} {n_fam:4d} {rho_w:+10.3f} {p_w:8.4f} {rho_d:+10.3f} {p_d:8.4f}")

    # Summary: sign consistency
    n_neg_W1 = sum(1 for r in fam_rho_W1 if r < 0 and not np.isnan(r))
    n_neg_dH = sum(1 for r in fam_rho_dH if r < 0 and not np.isnan(r))
    n_valid_W1 = sum(1 for r in fam_rho_W1 if not np.isnan(r))
    n_valid_dH = sum(1 for r in fam_rho_dH if not np.isnan(r))

    print(f"\n  Sign consistency (negative = low D_spec → high P_basin):")
    print(f"    W1:     {n_neg_W1}/{n_valid_W1} families negative")
    print(f"    ΔH_int: {n_neg_dH}/{n_valid_dH} families negative")

    # Within-family + within-N (strictest control)
    print(f"\n  Within (family, N) subgroups:")
    print(f"  {'family':6s} {'N':>3} {'n':>3} {'ρ(W1,Pb)':>10} {'ρ(dH,Pb)':>10}")

    subgroup_W1, subgroup_dH = [], []
    for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR"]:
        for N in N_VALUES:
            mask = (fam_arr == fam) & (N_arr == N)
            n_sub = mask.sum()
            if n_sub < 3:
                continue
            W1_s = W1_arr[mask]
            dH_s = dH_arr[mask]
            Pb_s = Pb_arr[mask]

            if np.std(W1_s) > 1e-10 and np.std(Pb_s) > 1e-10:
                rho_w, _ = sp_stats.spearmanr(W1_s, Pb_s)
            else:
                rho_w = float('nan')
            if np.std(dH_s) > 1e-10 and np.std(Pb_s) > 1e-10:
                rho_d, _ = sp_stats.spearmanr(dH_s, Pb_s)
            else:
                rho_d = float('nan')

            subgroup_W1.append(rho_w)
            subgroup_dH.append(rho_d)
            print(f"  {fam:6s} {N:3d} {n_sub:3d} {rho_w:+10.3f} {rho_d:+10.3f}")

    n_neg_sub_W1 = sum(1 for r in subgroup_W1 if r < 0 and not np.isnan(r))
    n_neg_sub_dH = sum(1 for r in subgroup_dH if r < 0 and not np.isnan(r))
    n_valid_sub_W1 = sum(1 for r in subgroup_W1 if not np.isnan(r))
    n_valid_sub_dH = sum(1 for r in subgroup_dH if not np.isnan(r))

    print(f"\n  Subgroup sign consistency:")
    print(f"    W1:     {n_neg_sub_W1}/{n_valid_sub_W1} subgroups negative")
    print(f"    ΔH_int: {n_neg_sub_dH}/{n_valid_sub_dH} subgroups negative")

    # ══════════════════════════════════════════════════════════════════
    # VERDICT
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 90)
    print("VERDICT SUMMARY")
    print("=" * 90)

    print(f"\n  Test 1 (partial corr | F7):     W1 ρ={rho_W1_partial:+.3f} p={p_W1_partial:.4f}, "
          f"ΔH ρ={rho_dH_partial:+.3f} p={p_dH_partial:.4f}")
    print(f"  Test 1 (partial corr | F7+N):   W1 ρ={rho_W1_p2:+.3f} p={p_W1_p2:.4f}, "
          f"ΔH ρ={rho_dH_p2:+.3f} p={p_dH_p2:.4f}")
    print(f"  Test 3 (within-family sign):    W1 {n_neg_W1}/{n_valid_W1} neg, "
          f"ΔH {n_neg_dH}/{n_valid_dH} neg")
    print(f"  Test 3 (within-fam×N sign):     W1 {n_neg_sub_W1}/{n_valid_sub_W1} neg, "
          f"ΔH {n_neg_sub_dH}/{n_valid_sub_dH} neg")

    W1_pass = p_W1_partial < 0.05
    dH_pass = p_dH_partial < 0.05
    W1_pass2 = p_W1_p2 < 0.05
    dH_pass2 = p_dH_p2 < 0.05

    print(f"\n  W1 independent of F7?     {'✓ YES' if W1_pass else '✗ NO'} (partial | F7)")
    print(f"  W1 independent of F7+N?   {'✓ YES' if W1_pass2 else '✗ NO'} (partial | F7,N)")
    print(f"  ΔH independent of F7?     {'✓ YES' if dH_pass else '✗ NO'} (partial | F7)")
    print(f"  ΔH independent of F7+N?   {'✓ YES' if dH_pass2 else '✗ NO'} (partial | F7,N)")

    if W1_pass or dH_pass:
        print("\n  ★ At least one spectral variable carries INDEPENDENT information "
              "about P_basin beyond F7.")
        if W1_pass and dH_pass:
            print("  → Both W1 and ΔH_int are independently significant.")
            print("  → D_main = W1 (strongest), D_clean = ΔH_int (most independent)")
        elif W1_pass:
            print("  → W1 is the primary D indicator.")
        else:
            print("  → ΔH_int is the primary D indicator (most independent).")
    else:
        print("\n  ✗ Neither variable is independent of F7 at p<0.05.")
        print("  → D remains an open problem; spectral variables may be F7-redundant.")

    # Save analysis table
    outdir = Path("outputs_d_recovery")
    outdir.mkdir(exist_ok=True)
    save_keys = [k for k in data[0].keys()]
    with open(outdir / "d_independence_analysis.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=save_keys)
        w.writeheader()
        w.writerows(data)
    print(f"\n  → {outdir / 'd_independence_analysis.csv'} ({len(data)} rows)")

    print("\n" + "=" * 90)
    print("INDEPENDENCE VERIFICATION COMPLETE")
    print("=" * 90)


if __name__ == "__main__":
    main()
