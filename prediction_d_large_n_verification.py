"""Prediction D Large-N Verification

Goal: Extend the dual-indicator (W1, ΔH_int) independence verification
from N ∈ {16,20,28,36} to N ∈ {36,52,72,100} to test whether the
coarse-graining stability signal scales to larger posets.

Three core questions:
  Q1. Do W1 and ΔH_int remain significant partial correlates of P_basin
      after controlling for F7, at large N?
  Q2. Does the family discrimination (Kruskal-Wallis) of W1/ΔH_int
      strengthen or weaken with N?
  Q3. Is there an N-scaling trend in the partial correlation strength?

Design:
  - N ∈ {36, 52, 72, 100}
  - 5 families: Lor2D, Lor3D, Lor4D, Lor5D, KR
  - 5 reps per (family, N) = 100 posets total
  - Per poset: compute F7, spectral vars (W1, ΔH_int), perturbation recovery
  - Adaptive: reduce SIS runs and dynamics steps at large N for feasibility
  - Output: CSV + Markdown report with N-scaling analysis

Usage:
    python prediction_d_large_n_verification.py [--ns 36 52 72 100] [--reps 5]
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from dataclasses import dataclass, field
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
# F7 computation (adaptive SIS runs based on N)
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


def compute_F7(poset: Poset, n_sis: int = 64,
               alpha0: float = 16.0, q: float = -0.5,
               lam: float = 10.0, eta: float = 0.6,
               Rc: float = 0.25, w: float = 0.015,
               N0: float = 20.0) -> float:
    N = poset.n
    log_H = compute_log_H(poset, n_runs=n_sis)
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim, _ = compute_xi_dim(poset)
    R = compute_R(poset)
    alpha_N = alpha0 * (N0 / N) ** abs(q)
    wall = alpha_N * sigmoid((R - Rc) / w)
    return log_H + 0.0004 * pi_geo - lam * sigma_hist + eta * xi_dim + wall


def compute_F7_fast(poset: Poset) -> float:
    return compute_F7(poset, n_sis=16)


def adaptive_sis_runs(N: int) -> int:
    """Reduce SIS runs at large N for feasibility."""
    if N <= 36:
        return 64
    elif N <= 52:
        return 48
    elif N <= 72:
        return 32
    else:
        return 24


def adaptive_recovery_steps(N: int) -> int:
    """Reduce recovery dynamics steps at large N."""
    if N <= 36:
        return 120
    elif N <= 52:
        return 80
    elif N <= 72:
        return 60
    else:
        return 40


# ══════════════════════════════════════════════════════════════════════════
# Hasse diagram + perturbation + dynamics
# ══════════════════════════════════════════════════════════════════════════

def hasse_diagram(closure: np.ndarray) -> np.ndarray:
    c = closure.astype(np.uint8)
    has_intermediate = (c @ c).astype(bool)
    cover = closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return cover


def perturb_add_edge(poset: Poset, rng: np.random.Generator) -> Poset | None:
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


def perturb_delete_edge(poset: Poset, rng: np.random.Generator) -> Poset | None:
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


def perturb_delete_node(poset: Poset, rng: np.random.Generator) -> Poset | None:
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


def propose_swap_move(closure: np.ndarray, rng: np.random.Generator) -> np.ndarray | None:
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


def recovery_dynamics(
    perturbed: Poset,
    original_cf: float,
    n_steps: int = 100,
    temperature: float = 1.0,
    seed: int = 42,
) -> dict:
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

def compute_spectral_vars(poset: Poset, keep_ratio: float = 0.7,
                          n_cg: int = 3, seed: int = 42) -> dict:
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
# Multi-step RG drift
# ══════════════════════════════════════════════════════════════════════════

def multistep_rg_drift(
    poset: Poset, n_steps: int = 3, keep_ratio: float = 0.7, seed: int = 42,
) -> dict:
    current = poset
    F_traj = [compute_F7(current, n_sis=adaptive_sis_runs(current.n))]
    phi_traj = [comparable_fraction(current)]

    for step in range(n_steps):
        if current.n < 6:
            break
        cg = coarse_grain_delete_nodes(current, keep_ratio=keep_ratio, seed=seed + step)
        if cg.n < 4:
            break
        current = cg
        F_traj.append(compute_F7(current, n_sis=adaptive_sis_runs(current.n)))
        phi_traj.append(comparable_fraction(current))

    cum_drift_F = sum(abs(F_traj[i + 1] - F_traj[i]) for i in range(len(F_traj) - 1))
    cum_drift_phi = sum(abs(phi_traj[i + 1] - phi_traj[i]) for i in range(len(phi_traj) - 1))

    step_drifts = [abs(F_traj[i + 1] - F_traj[i]) for i in range(len(F_traj) - 1)]
    converged = len(step_drifts) >= 2 and step_drifts[-1] < step_drifts[0]

    return {
        "rg_drift_F": cum_drift_F,
        "rg_drift_phi": cum_drift_phi,
        "rg_converged": converged,
        "rg_steps": len(F_traj) - 1,
    }


# ══════════════════════════════════════════════════════════════════════════
# Family generators
# ══════════════════════════════════════════════════════════════════════════

FAMILIES = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR": generate_kr_like,
}


# ══════════════════════════════════════════════════════════════════════════
# Report generation
# ══════════════════════════════════════════════════════════════════════════

def partial_corr(x, y, z):
    """Partial Spearman correlation of x and y controlling z."""
    rho_xz, _ = sp_stats.spearmanr(x, z)
    rho_yz, _ = sp_stats.spearmanr(y, z)
    res_x = sp_stats.rankdata(x) - rho_xz * sp_stats.rankdata(z)
    res_y = sp_stats.rankdata(y) - rho_yz * sp_stats.rankdata(z)
    rho_partial, p_partial = sp_stats.spearmanr(res_x, res_y)
    return rho_partial, p_partial


def partial_corr_2(x, y, z1, z2):
    """Partial Spearman correlation controlling z1 AND z2."""
    from numpy.linalg import lstsq
    Z = np.column_stack([sp_stats.rankdata(z1), sp_stats.rankdata(z2)])
    rx = sp_stats.rankdata(x)
    ry = sp_stats.rankdata(y)
    coef_x, _, _, _ = lstsq(Z, rx, rcond=None)
    coef_y, _, _, _ = lstsq(Z, ry, rcond=None)
    res_x = rx - Z @ coef_x
    res_y = ry - Z @ coef_y
    return sp_stats.spearmanr(res_x, res_y)


def generate_report(data: list[dict], outdir: Path) -> str:
    """Generate Markdown report with full N-scaling analysis."""
    lines = []
    lines.append("# Prediction D: Large-N Verification Report\n")
    lines.append(f"**Total samples**: {len(data)}\n")

    N_vals = sorted(set(d["N"] for d in data))
    lines.append(f"**N values**: {N_vals}\n")
    lines.append(f"**Families**: {sorted(set(d['family'] for d in data))}\n")

    F7_arr = np.array([d["F7"] for d in data])
    W1_arr = np.array([d["wasserstein"] for d in data])
    dH_arr = np.array([d["delta_H_int"] for d in data])
    Pb_arr = np.array([d["P_basin"] for d in data])
    N_arr = np.array([d["N"] for d in data])
    fam_arr = np.array([d["family"] for d in data])

    # ── Q1: Pooled partial correlations ──
    lines.append("\n---\n")
    lines.append("## Q1: Pooled Partial Correlations\n")

    rho_F7_Pb, p_F7_Pb = sp_stats.spearmanr(F7_arr, Pb_arr)
    rho_W1_Pb, p_W1_Pb = sp_stats.spearmanr(W1_arr, Pb_arr)
    rho_dH_Pb, p_dH_Pb = sp_stats.spearmanr(dH_arr, Pb_arr)

    lines.append("### Raw correlations with P_basin\n")
    lines.append("| Variable | ρ | p |\n|----------|---|---|\n")
    lines.append(f"| F7 | {rho_F7_Pb:+.3f} | {p_F7_Pb:.4f} |\n")
    lines.append(f"| W1 | {rho_W1_Pb:+.3f} | {p_W1_Pb:.4f} |\n")
    lines.append(f"| ΔH_int | {rho_dH_Pb:+.3f} | {p_dH_Pb:.4f} |\n")

    rho_W1_p, p_W1_p = partial_corr(W1_arr, Pb_arr, F7_arr)
    rho_dH_p, p_dH_p = partial_corr(dH_arr, Pb_arr, F7_arr)

    lines.append("\n### Partial correlations (controlling F7)\n")
    lines.append("| Variable | ρ_partial | p | sig |\n|----------|-----------|---|-----|\n")
    sig_W1 = "★★★" if p_W1_p < 0.001 else ("★★" if p_W1_p < 0.01 else ("★" if p_W1_p < 0.05 else "ns"))
    sig_dH = "★★★" if p_dH_p < 0.001 else ("★★" if p_dH_p < 0.01 else ("★" if p_dH_p < 0.05 else "ns"))
    lines.append(f"| W1 \\| F7 | {rho_W1_p:+.3f} | {p_W1_p:.4f} | {sig_W1} |\n")
    lines.append(f"| ΔH_int \\| F7 | {rho_dH_p:+.3f} | {p_dH_p:.4f} | {sig_dH} |\n")

    rho_W1_p2, p_W1_p2 = partial_corr_2(W1_arr, Pb_arr, F7_arr, N_arr)
    rho_dH_p2, p_dH_p2 = partial_corr_2(dH_arr, Pb_arr, F7_arr, N_arr)

    lines.append("\n### Partial correlations (controlling F7 + N)\n")
    lines.append("| Variable | ρ_partial | p | sig |\n|----------|-----------|---|-----|\n")
    sig_W1_2 = "★★★" if p_W1_p2 < 0.001 else ("★★" if p_W1_p2 < 0.01 else ("★" if p_W1_p2 < 0.05 else "ns"))
    sig_dH_2 = "★★★" if p_dH_p2 < 0.001 else ("★★" if p_dH_p2 < 0.01 else ("★" if p_dH_p2 < 0.05 else "ns"))
    lines.append(f"| W1 \\| F7,N | {rho_W1_p2:+.3f} | {p_W1_p2:.4f} | {sig_W1_2} |\n")
    lines.append(f"| ΔH_int \\| F7,N | {rho_dH_p2:+.3f} | {p_dH_p2:.4f} | {sig_dH_2} |\n")

    # ── Q2: Family discrimination per N (Kruskal-Wallis) ──
    lines.append("\n---\n")
    lines.append("## Q2: Family Discrimination per N (Kruskal-Wallis)\n")
    lines.append("| N | var | H_stat | p | sig |\n|---|-----|--------|---|-----|\n")

    kw_results = {}
    for N in N_vals:
        mask_N = N_arr == N
        for var_name, var_arr in [("W1", W1_arr), ("ΔH_int", dH_arr),
                                   ("rg_drift_F", np.array([d["rg_drift_F"] for d in data]))]:
            groups = []
            for fam in sorted(FAMILIES.keys()):
                vals = var_arr[(fam_arr == fam) & mask_N]
                if len(vals) >= 2:
                    groups.append(vals)
            if len(groups) >= 2:
                try:
                    H_stat, p_val = sp_stats.kruskal(*groups)
                except ValueError:
                    H_stat, p_val = 0.0, 1.0
                sig = "★★★" if p_val < 0.001 else ("★★" if p_val < 0.01 else ("★" if p_val < 0.05 else "ns"))
                lines.append(f"| {N} | {var_name} | {H_stat:.2f} | {p_val:.4f} | {sig} |\n")
                kw_results[(N, var_name)] = (H_stat, p_val)

    # ── Q3: N-scaling of partial correlation ──
    lines.append("\n---\n")
    lines.append("## Q3: N-Scaling of Partial Correlation Strength\n")
    lines.append("| N | n | ρ(W1,Pb\\|F7) | p | ρ(ΔH,Pb\\|F7) | p |\n")
    lines.append("|---|---|-------------|---|--------------|---|\n")

    n_scaling_W1 = []
    n_scaling_dH = []
    for N in N_vals:
        mask_N = N_arr == N
        n_N = mask_N.sum()
        if n_N < 8:
            lines.append(f"| {N} | {n_N} | — | — | — | — |\n")
            continue
        F7_N = F7_arr[mask_N]
        W1_N = W1_arr[mask_N]
        dH_N = dH_arr[mask_N]
        Pb_N = Pb_arr[mask_N]
        if np.std(W1_N) > 1e-10 and np.std(Pb_N) > 1e-10 and np.std(F7_N) > 1e-10:
            rho_w, p_w = partial_corr(W1_N, Pb_N, F7_N)
            n_scaling_W1.append((N, rho_w))
        else:
            rho_w, p_w = float("nan"), 1.0
        if np.std(dH_N) > 1e-10 and np.std(Pb_N) > 1e-10 and np.std(F7_N) > 1e-10:
            rho_d, p_d = partial_corr(dH_N, Pb_N, F7_N)
            n_scaling_dH.append((N, rho_d))
        else:
            rho_d, p_d = float("nan"), 1.0
        lines.append(f"| {N} | {n_N} | {rho_w:+.3f} | {p_w:.4f} | {rho_d:+.3f} | {p_d:.4f} |\n")

    # Spearman of |ρ_partial| vs N
    if len(n_scaling_W1) >= 3:
        ns_w = [x[0] for x in n_scaling_W1]
        rs_w = [abs(x[1]) for x in n_scaling_W1]
        rho_trend_W1, p_trend_W1 = sp_stats.spearmanr(ns_w, rs_w)
        lines.append(f"\n**W1 |ρ_partial| trend with N**: Spearman ρ = {rho_trend_W1:+.3f}, p = {p_trend_W1:.4f}\n")
    if len(n_scaling_dH) >= 3:
        ns_d = [x[0] for x in n_scaling_dH]
        rs_d = [abs(x[1]) for x in n_scaling_dH]
        rho_trend_dH, p_trend_dH = sp_stats.spearmanr(ns_d, rs_d)
        lines.append(f"**ΔH_int |ρ_partial| trend with N**: Spearman ρ = {rho_trend_dH:+.3f}, p = {p_trend_dH:.4f}\n")

    # ── Family means per N ──
    lines.append("\n---\n")
    lines.append("## Family Means per N\n")
    lines.append("| N | family | n | mean_F7 | mean_W1 | mean_ΔH | mean_Pb | mean_rg_dF | mean_cf |\n")
    lines.append("|---|--------|---|---------|---------|---------|---------|------------|--------|\n")
    for N in N_vals:
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR"]:
            fam_data = [d for d in data if d["N"] == N and d["family"] == fam]
            if not fam_data:
                continue
            mF7 = np.mean([d["F7"] for d in fam_data])
            mW1 = np.mean([d["wasserstein"] for d in fam_data])
            mdH = np.mean([d["delta_H_int"] for d in fam_data])
            mPb = np.mean([d["P_basin"] for d in fam_data])
            mRG = np.mean([d["rg_drift_F"] for d in fam_data])
            mCF = np.mean([d["comp_frac"] for d in fam_data])
            lines.append(f"| {N} | {fam} | {len(fam_data)} | {mF7:.2f} | {mW1:.4f} | "
                         f"{mdH:.4f} | {mPb:.1%} | {mRG:.2f} | {mCF:.3f} |\n")

    # ── Within-family tests ──
    lines.append("\n---\n")
    lines.append("## Within-Family Partial Correlations\n")
    lines.append("| family | n | ρ(W1,Pb) | p | ρ(ΔH,Pb) | p |\n")
    lines.append("|--------|---|----------|---|----------|---|\n")
    for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR"]:
        mask = fam_arr == fam
        n_f = mask.sum()
        if n_f < 5:
            continue
        if np.std(W1_arr[mask]) > 1e-10 and np.std(Pb_arr[mask]) > 1e-10:
            rw, pw = sp_stats.spearmanr(W1_arr[mask], Pb_arr[mask])
        else:
            rw, pw = float("nan"), 1.0
        if np.std(dH_arr[mask]) > 1e-10 and np.std(Pb_arr[mask]) > 1e-10:
            rd, pd = sp_stats.spearmanr(dH_arr[mask], Pb_arr[mask])
        else:
            rd, pd = float("nan"), 1.0
        lines.append(f"| {fam} | {n_f} | {rw:+.3f} | {pw:.4f} | {rd:+.3f} | {pd:.4f} |\n")

    # ── RG drift convergence ──
    lines.append("\n---\n")
    lines.append("## RG Drift N-Scaling\n")
    lines.append("| N | family | mean_rg_dF | mean_rg_dφ | converged_frac |\n")
    lines.append("|---|--------|-----------|-----------|---------------|\n")
    for N in N_vals:
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR"]:
            fd = [d for d in data if d["N"] == N and d["family"] == fam]
            if not fd:
                continue
            mdf = np.mean([d["rg_drift_F"] for d in fd])
            mdp = np.mean([d["rg_drift_phi"] for d in fd])
            conv = np.mean([d["rg_converged"] for d in fd])
            lines.append(f"| {N} | {fam} | {mdf:.2f} | {mdp:.4f} | {conv:.0%} |\n")

    # ── Verdict ──
    lines.append("\n---\n")
    lines.append("## Verdict\n")

    W1_indep = p_W1_p < 0.05
    dH_indep = p_dH_p < 0.05
    W1_indep2 = p_W1_p2 < 0.05
    dH_indep2 = p_dH_p2 < 0.05

    lines.append(f"- **W1 independent of F7**: {'✓ YES' if W1_indep else '✗ NO'} "
                 f"(ρ={rho_W1_p:+.3f}, p={p_W1_p:.4f})\n")
    lines.append(f"- **W1 independent of F7+N**: {'✓ YES' if W1_indep2 else '✗ NO'} "
                 f"(ρ={rho_W1_p2:+.3f}, p={p_W1_p2:.4f})\n")
    lines.append(f"- **ΔH_int independent of F7**: {'✓ YES' if dH_indep else '✗ NO'} "
                 f"(ρ={rho_dH_p:+.3f}, p={p_dH_p:.4f})\n")
    lines.append(f"- **ΔH_int independent of F7+N**: {'✓ YES' if dH_indep2 else '✗ NO'} "
                 f"(ρ={rho_dH_p2:+.3f}, p={p_dH_p2:.4f})\n")

    # N-scaling verdict
    if len(n_scaling_W1) >= 3:
        trend = "strengthening ↑" if rho_trend_W1 > 0.3 else ("weakening ↓" if rho_trend_W1 < -0.3 else "stable ↔")
        lines.append(f"- **W1 N-scaling trend**: {trend} (ρ={rho_trend_W1:+.3f})\n")
    if len(n_scaling_dH) >= 3:
        trend = "strengthening ↑" if rho_trend_dH > 0.3 else ("weakening ↓" if rho_trend_dH < -0.3 else "stable ↔")
        lines.append(f"- **ΔH_int N-scaling trend**: {trend} (ρ={rho_trend_dH:+.3f})\n")

    # Overall
    n_pass = sum([W1_indep, dH_indep, W1_indep2, dH_indep2])
    if n_pass >= 3:
        verdict = "**STRONG**: Both indicators survive large-N extension"
    elif n_pass >= 1:
        verdict = "**PARTIAL**: At least one indicator survives"
    else:
        verdict = "**WEAK**: Neither indicator survives large-N — D remains fragile"
    lines.append(f"\n### Overall: {verdict}\n")

    report = "".join(lines)
    return report


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Prediction D Large-N Verification")
    parser.add_argument("--ns", type=int, nargs="+", default=[36, 52, 72, 100])
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--pert-reps", type=int, default=2)
    parser.add_argument("--outdir", type=str, default="outputs_d_recovery")
    parser.add_argument("--suffix", type=str, default="",
                        help="Suffix for output filenames, e.g. '_r15'")
    args = parser.parse_args()

    N_VALUES = tuple(args.ns)
    REPS = args.reps
    PERT_REPS = args.pert_reps
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    total_start = time.time()

    print("=" * 90)
    print("PREDICTION D: LARGE-N VERIFICATION")
    print(f"N = {N_VALUES}, {REPS} reps/family, {len(FAMILIES)} families")
    print(f"Total posets: {len(N_VALUES) * len(FAMILIES) * REPS}")
    print("=" * 90)

    # ── Phase 1: Generate samples + compute F7 & spectral vars ──
    print(f"\n[1/4] Generating samples...")
    samples = []
    seed_base = 7777

    for N in N_VALUES:
        t0 = time.time()
        n_sis = adaptive_sis_runs(N)
        for fam_name, gen_func in FAMILIES.items():
            for rep in range(REPS):
                s = seed_base + rep * 1000 + N * 7
                poset = gen_func(N, seed=s)
                F7 = compute_F7(poset, n_sis=n_sis)
                cf = comparable_fraction(poset)
                sv = compute_spectral_vars(poset, seed=s)
                samples.append({
                    "family": fam_name, "N": N, "rep": rep,
                    "poset": poset, "F7": F7, "comp_frac": cf,
                    **sv,
                })
        dt = time.time() - t0
        done_n = sum(1 for s in samples if s["N"] == N)
        print(f"  N={N}: {done_n} samples ({dt:.1f}s, SIS={n_sis})")

    print(f"  Total: {len(samples)} samples ({time.time()-total_start:.1f}s)")

    # ── Phase 2: Perturbation + recovery for P_basin ──
    print(f"\n[2/4] Perturbation + recovery...")
    recovery_map = {}

    trial_count = 0
    for si, sample in enumerate(samples):
        key = (sample["family"], sample["N"], sample["rep"])
        recovery_map[key] = []

        n_steps = adaptive_recovery_steps(sample["N"])
        for pert_name, pert_func in PERTURBATIONS.items():
            for pr in range(PERT_REPS):
                pseed = abs(hash((key, pert_name, pr))) % (2**31)
                pert_rng = np.random.default_rng(pseed)
                perturbed = pert_func(sample["poset"], pert_rng)
                if perturbed is None:
                    continue

                rec = recovery_dynamics(
                    perturbed, sample["comp_frac"],
                    n_steps=n_steps, temperature=1.0,
                    seed=(pseed + 7) % (2**31),
                )
                recovery_map[key].append(rec["basin_retention"])
                trial_count += 1

        if (si + 1) % 10 == 0 or si + 1 == len(samples):
            elapsed = time.time() - total_start
            print(f"  {si+1}/{len(samples)} samples ({trial_count} trials, {elapsed:.0f}s)")

    print(f"  Total recovery trials: {trial_count}")

    # ── Phase 3: RG drift ──
    print(f"\n[3/4] Computing RG drift...")
    rg_results = {}
    for si, sample in enumerate(samples):
        key = (sample["family"], sample["N"], sample["rep"])
        rg = multistep_rg_drift(sample["poset"], n_steps=3, seed=seed_base + sample["rep"])
        rg_results[key] = rg
        if (si + 1) % 10 == 0 or si + 1 == len(samples):
            print(f"  {si+1}/{len(samples)}")

    # ── Merge into analysis table ──
    data = []
    for sample in samples:
        key = (sample["family"], sample["N"], sample["rep"])
        brs = recovery_map.get(key, [])
        rg = rg_results.get(key, {})
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
            "H_int_orig": sample["H_int_orig"],
            "P_basin": float(np.mean(brs)),
            "n_trials": len(brs),
            "rg_drift_F": rg.get("rg_drift_F", 0.0),
            "rg_drift_phi": rg.get("rg_drift_phi", 0.0),
            "rg_converged": rg.get("rg_converged", False),
        })

    print(f"\n  Analysis table: {len(data)} rows")

    # ── Phase 4: Analysis + Report ──
    print(f"\n[4/4] Analysis and report generation...")

    report = generate_report(data, outdir)

    # Save CSV
    csv_path = outdir / f"prediction_d_large_n{args.suffix}.csv"
    if data:
        keys = [k for k in data[0].keys()]
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(data)
        print(f"  → {csv_path} ({len(data)} rows)")

    # Save report
    md_path = outdir / f"prediction_d_large_n{args.suffix}.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"  → {md_path}")

    # Console summary
    print("\n" + "=" * 90)
    # Extract key numbers from data
    F7_arr = np.array([d["F7"] for d in data])
    W1_arr = np.array([d["wasserstein"] for d in data])
    dH_arr = np.array([d["delta_H_int"] for d in data])
    Pb_arr = np.array([d["P_basin"] for d in data])

    rho_W1_p, p_W1_p = partial_corr(W1_arr, Pb_arr, F7_arr)
    rho_dH_p, p_dH_p = partial_corr(dH_arr, Pb_arr, F7_arr)

    print(f"  Pooled partial corr (|F7):")
    print(f"    W1:     ρ = {rho_W1_p:+.3f}, p = {p_W1_p:.4f}")
    print(f"    ΔH_int: ρ = {rho_dH_p:+.3f}, p = {p_dH_p:.4f}")

    total_time = time.time() - total_start
    print(f"\n  Total time: {total_time:.0f}s ({total_time/60:.1f}min)")
    print("=" * 90)
    print("LARGE-N VERIFICATION COMPLETE")
    print("=" * 90)


if __name__ == "__main__":
    main()
