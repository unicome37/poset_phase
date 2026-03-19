"""Prediction D Recovery Experiment — Dynamic Stability under F7

New D hypothesis: "structures that survive in the Lorentzian window are not
merely low-F7, but dynamically stable — they return to the window after
perturbation."

This replaces the falsified Π_cg closure hypothesis (§5.8.6) with a
dynamical recovery framework.

Three perturbation types:
  1. Add a random cover edge
  2. Delete a random cover edge
  3. Small-scale coarse-grain (delete 1-2 nodes)

Four candidate D variables:
  (a) Return probability P_return — does the system return to the window?
  (b) Return time τ_return — how many steps to return?
  (c) Multi-step ΔF7 cumulative drift under iterated CG
  (d) Spectral fidelity — Wasserstein distance of interval spectrum pre/post CG

Usage:
    python _d_recovery_experiment.py
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from generators import Poset, transitive_closure
from unified_functional import (
    FunctionalWeights,
    compute_log_H,
    compute_pi_geo,
    compute_sigma_hist,
    compute_xi_dim,
)
from observables import comparable_fraction, layer_profile
from observables_geo import height_ratio, width_ratio
from coarse_grain import coarse_grain_delete_nodes
from prediction_a_mid_gap import interval_spectrum
from prediction_a_bd_bridge import count_intervals_fast

# ══════════════════════════════════════════════════════════════════════════
# F7 computation (no κΠ_cg, with sigmoid wall + N-scaling)
# ══════════════════════════════════════════════════════════════════════════

def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_R(poset: Poset) -> float:
    """Interval occupancy ratio R = 1 - f_link."""
    counts = count_intervals_fast(poset)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    C0 = counts.get(0, 0)
    return 1.0 - C0 / total


def compute_F7(
    poset: Poset,
    alpha0: float = 16.0,
    q: float = -0.5,
    lam: float = 10.0,
    eta: float = 0.6,
    Rc: float = 0.25,
    w: float = 0.015,
    N0: float = 20.0,
    n_sis: int = 64,
) -> float:
    """F7 = logH + γΠ_geo - λΣ_hist + ηΞ_d + α(N)σ((R-Rc)/w).  No κΠ_cg."""
    N = poset.n
    log_H = compute_log_H(poset, n_runs=n_sis)
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim, _ = compute_xi_dim(poset)
    R = compute_R(poset)
    alpha_N = alpha0 * (N0 / N) ** abs(q)
    wall = alpha_N * sigmoid((R - Rc) / w)
    return log_H + 0.0004 * pi_geo - lam * sigma_hist + eta * xi_dim + wall


def compute_F7_fast(poset: Poset, **kw) -> float:
    """Alias with fewer SIS samples for dynamics inner loop."""
    return compute_F7(poset, n_sis=32, **kw)


# ══════════════════════════════════════════════════════════════════════════
# Hasse diagram utilities
# ══════════════════════════════════════════════════════════════════════════

def hasse_diagram(closure: np.ndarray) -> np.ndarray:
    n = closure.shape[0]
    c = closure.astype(np.uint8)
    has_intermediate = (c @ c).astype(bool)
    cover = closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return cover


# ══════════════════════════════════════════════════════════════════════════
# Perturbation operators
# ══════════════════════════════════════════════════════════════════════════

def perturb_add_edge(poset: Poset, rng: np.random.Generator) -> Poset | None:
    """Add a random covering relation between two incomparable elements."""
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
    """Delete a random covering relation."""
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
    """Delete a random node (small-scale coarse-grain)."""
    if poset.n < 6:
        return None
    node = rng.integers(poset.n)
    keep = [i for i in range(poset.n) if i != node]
    new_closure = poset.closure[np.ix_(keep, keep)]
    return Poset(new_closure)


PERTURBATIONS = {
    "add_edge": perturb_add_edge,
    "del_edge": perturb_delete_edge,
    "del_node": perturb_delete_node,
}


# ══════════════════════════════════════════════════════════════════════════
# Swap-move dynamics (microcanonical, from unified_functional_metropolis)
# ══════════════════════════════════════════════════════════════════════════

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


def is_in_window(poset: Poset) -> bool:
    """Check if poset is in the Lorentzian-like window (comp_frac 0.30-0.55)."""
    cf = comparable_fraction(poset)
    return 0.30 <= cf <= 0.55


def recovery_dynamics(
    perturbed: Poset,
    original_F: float,
    original_cf: float,
    n_steps: int = 200,
    temperature: float = 1.0,
    seed: int = 42,
) -> dict:
    """Run swap-move dynamics from perturbed state, track recovery metrics.
    
    Returns:
      returned_to_window: bool
      return_step: int (-1 if never returned)
      F_trajectory: list of F values at each step
      cf_trajectory: list of comp_frac at each step
      final_F: float
      basin_retention: bool (final comp_frac within ±0.05 of original)
    """
    rng = np.random.default_rng(seed)
    current = perturbed
    current_F = compute_F7_fast(current)
    
    F_traj = [current_F]
    cf_traj = [comparable_fraction(current)]
    returned = False
    return_step = -1
    
    for step in range(1, n_steps + 1):
        new_closure = propose_swap_move(current.closure, rng)
        if new_closure is None:
            F_traj.append(current_F)
            cf_traj.append(cf_traj[-1])
            continue
        
        new_poset = Poset(new_closure)
        new_F = compute_F7_fast(new_poset)
        delta_F = new_F - current_F
        
        if delta_F < 0 or rng.random() < math.exp(-delta_F / max(temperature, 1e-10)):
            current = new_poset
            current_F = new_F
        
        cf = comparable_fraction(current)
        F_traj.append(current_F)
        cf_traj.append(cf)
        
        if not returned and is_in_window(current):
            returned = True
            return_step = step
    
    final_cf = cf_traj[-1]
    basin_retention = abs(final_cf - original_cf) < 0.05
    
    return {
        "returned_to_window": returned,
        "return_step": return_step,
        "final_F": current_F,
        "delta_F_recovery": current_F - original_F,
        "basin_retention": basin_retention,
        "F_trajectory": F_traj,
        "cf_trajectory": cf_traj,
    }


# ══════════════════════════════════════════════════════════════════════════
# Candidate D variable: Multi-step RG drift
# ══════════════════════════════════════════════════════════════════════════

def multistep_rg_drift(
    poset: Poset,
    n_steps: int = 3,
    keep_ratio: float = 0.7,
    seed: int = 42,
) -> dict:
    """Multi-step coarse-grain: iterate CG and track cumulative F7 drift.
    
    Returns:
      F_rg_trajectory: list of F7 at each RG step
      cumulative_drift: total |ΔF7| across steps
      phi_trajectory: list of comp_frac at each RG step
      phi_drift: total |Δφ| across steps (φ = r_comp = comp_frac)
      converged: bool (drift decreasing across steps)
    """
    rng = np.random.default_rng(seed)
    current = poset
    F_traj = [compute_F7(current)]
    phi_traj = [comparable_fraction(current)]
    
    for step in range(n_steps):
        if current.n < 6:
            break
        cg = coarse_grain_delete_nodes(current, keep_ratio=keep_ratio, seed=seed + step)
        if cg.n < 4:
            break
        current = cg
        F_traj.append(compute_F7(current))
        phi_traj.append(comparable_fraction(current))
    
    cum_drift_F = sum(abs(F_traj[i+1] - F_traj[i]) for i in range(len(F_traj)-1))
    cum_drift_phi = sum(abs(phi_traj[i+1] - phi_traj[i]) for i in range(len(phi_traj)-1))
    
    # Check if drift is decreasing (converging)
    step_drifts = [abs(F_traj[i+1] - F_traj[i]) for i in range(len(F_traj)-1)]
    converged = len(step_drifts) >= 2 and step_drifts[-1] < step_drifts[0]
    
    return {
        "F_rg_trajectory": F_traj,
        "cumulative_drift_F": cum_drift_F,
        "phi_trajectory": phi_traj,
        "cumulative_drift_phi": cum_drift_phi,
        "converged": converged,
        "n_rg_steps": len(F_traj) - 1,
    }


# ══════════════════════════════════════════════════════════════════════════
# Candidate D variable: Spectral fidelity
# ══════════════════════════════════════════════════════════════════════════

def spectral_fidelity(
    poset: Poset,
    keep_ratio: float = 0.7,
    n_cg_samples: int = 3,
    seed: int = 42,
) -> dict:
    """Compare interval spectrum before/after CG.
    
    Returns:
      spectral_entropy_orig: H_int of original
      spectral_entropy_cg_mean: mean H_int of CG'd versions
      delta_H_int: |ΔH_int|
      wasserstein_mean: mean Wasserstein-1 distance of C_k distributions
    """
    spec_orig = interval_spectrum(poset)
    H_orig = spec_orig["spectral_entropy"]
    
    # Get original C_k distribution
    counts_orig = count_intervals_fast(poset)
    total_orig = sum(counts_orig.values())
    if total_orig == 0:
        return {
            "spectral_entropy_orig": 0.0,
            "spectral_entropy_cg_mean": 0.0,
            "delta_H_int": 0.0,
            "wasserstein_mean": 0.0,
        }
    max_k_orig = max(counts_orig.keys()) if counts_orig else 0
    dist_orig = np.array([counts_orig.get(k, 0) / total_orig for k in range(max_k_orig + 1)])
    
    H_cg_list = []
    wass_list = []
    
    for i in range(n_cg_samples):
        cg = coarse_grain_delete_nodes(poset, keep_ratio=keep_ratio, seed=seed + i)
        if cg.n < 4:
            continue
        spec_cg = interval_spectrum(cg)
        H_cg_list.append(spec_cg["spectral_entropy"])
        
        # Wasserstein-1 distance of C_k distributions
        counts_cg = count_intervals_fast(cg)
        total_cg = sum(counts_cg.values())
        if total_cg == 0:
            continue
        max_k_cg = max(counts_cg.keys()) if counts_cg else 0
        max_k = max(max_k_orig, max_k_cg)
        dist_cg = np.array([counts_cg.get(k, 0) / total_cg for k in range(max_k + 1)])
        dist_o = np.array([counts_orig.get(k, 0) / total_orig for k in range(max_k + 1)])
        
        # Wasserstein-1 = sum of |CDF_diff|
        cdf_o = np.cumsum(dist_o)
        cdf_c = np.cumsum(dist_cg)
        wass = float(np.sum(np.abs(cdf_o - cdf_c)))
        wass_list.append(wass)
    
    H_cg_mean = float(np.mean(H_cg_list)) if H_cg_list else H_orig
    wass_mean = float(np.mean(wass_list)) if wass_list else 0.0
    
    return {
        "spectral_entropy_orig": H_orig,
        "spectral_entropy_cg_mean": H_cg_mean,
        "delta_H_int": abs(H_cg_mean - H_orig),
        "wasserstein_mean": wass_mean,
    }


# ══════════════════════════════════════════════════════════════════════════
# Sample generation
# ══════════════════════════════════════════════════════════════════════════

def generate_samples(N_values=(16, 20), reps_per_family=3, seed=42):
    """Generate representative posets for each family at each N."""
    from generators import (
        generate_lorentzian_like_2d,
        generate_lorentzian_like_3d,
        generate_lorentzian_like_4d,
        generate_lorentzian_like_5d,
        generate_kr_like,
    )
    
    families = {
        "Lor2D": generate_lorentzian_like_2d,
        "Lor3D": generate_lorentzian_like_3d,
        "Lor4D": generate_lorentzian_like_4d,
        "Lor5D": generate_lorentzian_like_5d,
        "KR": generate_kr_like,
    }
    
    samples = []
    for N in N_values:
        for fam_name, gen_func in families.items():
            for rep in range(reps_per_family):
                s = seed + rep * 1000
                poset = gen_func(N, seed=s)
                F = compute_F7(poset)
                cf = comparable_fraction(poset)
                in_win = is_in_window(poset)
                samples.append({
                    "family": fam_name,
                    "N": N,
                    "rep": rep,
                    "poset": poset,
                    "F7": F,
                    "comp_frac": cf,
                    "in_window": in_win,
                })
    return samples


# ══════════════════════════════════════════════════════════════════════════
# Main experiment
# ══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 90)
    print("PREDICTION D RECOVERY EXPERIMENT")
    print("Dynamic stability under F7 (replacing falsified Π_cg closure)")
    print("=" * 90)
    
    # Phase 1: Generate samples
    print("\n[1/4] Generating representative posets...")
    samples = generate_samples(N_values=(16, 20), reps_per_family=3, seed=42)
    
    n_window = sum(1 for s in samples if s["in_window"])
    print(f"  Total samples: {len(samples)}, in window: {n_window}")
    for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR"]:
        fam_samples = [s for s in samples if s["family"] == fam]
        n_w = sum(1 for s in fam_samples if s["in_window"])
        mean_F = np.mean([s["F7"] for s in fam_samples])
        mean_cf = np.mean([s["comp_frac"] for s in fam_samples])
        print(f"  {fam:6s}: {len(fam_samples)} samples, {n_w} in window, "
              f"mean F7={mean_F:.2f}, mean cf={mean_cf:.3f}")
    
    # Phase 2: Perturbation + recovery
    print("\n[2/4] Perturbation + recovery dynamics...")
    print(f"  {'family':6s} {'N':>3} {'rep':>3} {'perturb':>10} "
          f"{'returned':>8} {'τ_ret':>5} {'ΔF':>8} {'basin':>5}")
    print("  " + "-" * 60)
    
    recovery_results = []
    rng = np.random.default_rng(42)
    
    for sample in samples:
        for pert_name, pert_func in PERTURBATIONS.items():
            for pert_rep in range(3):
                pert_seed = 1000 * sample["rep"] + 100 * pert_rep + hash(pert_name) % 100
                pert_rng = np.random.default_rng(abs(pert_seed) % (2**31))
                
                perturbed = pert_func(sample["poset"], pert_rng)
                if perturbed is None:
                    continue
                
                rec = recovery_dynamics(
                    perturbed,
                    original_F=sample["F7"],
                    original_cf=sample["comp_frac"],
                    n_steps=150,
                    temperature=1.0,
                    seed=abs(pert_seed + 7) % (2**31),
                )
                
                row = {
                    "family": sample["family"],
                    "N": sample["N"],
                    "rep": sample["rep"],
                    "perturbation": pert_name,
                    "pert_rep": pert_rep,
                    "in_window_orig": sample["in_window"],
                    "F7_orig": sample["F7"],
                    "cf_orig": sample["comp_frac"],
                    "returned": rec["returned_to_window"],
                    "return_step": rec["return_step"],
                    "delta_F_recovery": rec["delta_F_recovery"],
                    "basin_retention": rec["basin_retention"],
                    "final_F": rec["final_F"],
                }
                recovery_results.append(row)
                
                if len(recovery_results) % 20 == 1:
                    print(f"  {row['family']:6s} {row['N']:3d} {row['rep']:3d} "
                          f"{row['perturbation']:>10s} "
                          f"{'✓' if row['returned'] else '✗':>8s} "
                          f"{row['return_step']:5d} "
                          f"{row['delta_F_recovery']:+8.2f} "
                          f"{'✓' if row['basin_retention'] else '✗':>5s}")
    
    print(f"\n  Total recovery trials: {len(recovery_results)}")
    
    # Phase 3: Candidate D variables
    print("\n[3/4] Computing candidate D variables...")
    
    d_candidates = []
    for sample in samples:
        # Multi-step RG drift
        rg = multistep_rg_drift(sample["poset"], n_steps=3, seed=42 + sample["rep"])
        
        # Spectral fidelity
        sf = spectral_fidelity(sample["poset"], seed=42 + sample["rep"])
        
        d_candidates.append({
            "family": sample["family"],
            "N": sample["N"],
            "rep": sample["rep"],
            "F7": sample["F7"],
            "comp_frac": sample["comp_frac"],
            "in_window": sample["in_window"],
            # RG drift
            "rg_drift_F": rg["cumulative_drift_F"],
            "rg_drift_phi": rg["cumulative_drift_phi"],
            "rg_converged": rg["converged"],
            "rg_steps": rg["n_rg_steps"],
            # Spectral
            "H_int_orig": sf["spectral_entropy_orig"],
            "H_int_cg": sf["spectral_entropy_cg_mean"],
            "delta_H_int": sf["delta_H_int"],
            "wasserstein": sf["wasserstein_mean"],
        })
        
        if len(d_candidates) % 10 == 1:
            dc = d_candidates[-1]
            print(f"  {dc['family']:6s} N={dc['N']:2d} rep={dc['rep']}: "
                  f"RG_drift_F={dc['rg_drift_F']:.2f}, "
                  f"RG_drift_φ={dc['rg_drift_phi']:.4f}, "
                  f"ΔH_int={dc['delta_H_int']:.3f}, "
                  f"W1={dc['wasserstein']:.4f}")
    
    # Phase 4: Analysis
    print("\n[4/4] Analysis...")
    print("=" * 90)
    
    # 4a: Recovery rates by family
    print("\n--- Recovery rates by family ---")
    print(f"  {'family':6s} {'n_trials':>8} {'P_return':>8} "
          f"{'mean_τ':>7} {'P_basin':>8}")
    for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR"]:
        fam_rec = [r for r in recovery_results if r["family"] == fam]
        if not fam_rec:
            continue
        n_ret = sum(1 for r in fam_rec if r["returned"])
        n_basin = sum(1 for r in fam_rec if r["basin_retention"])
        ret_steps = [r["return_step"] for r in fam_rec if r["returned"]]
        mean_tau = np.mean(ret_steps) if ret_steps else -1
        print(f"  {fam:6s} {len(fam_rec):8d} {n_ret/len(fam_rec):8.1%} "
              f"{mean_tau:7.1f} {n_basin/len(fam_rec):8.1%}")
    
    # 4b: Recovery rates by perturbation type
    print("\n--- Recovery rates by perturbation type ---")
    print(f"  {'perturb':>10} {'n_trials':>8} {'P_return':>8} {'mean_τ':>7}")
    for pert in ["add_edge", "del_edge", "del_node"]:
        pr = [r for r in recovery_results if r["perturbation"] == pert]
        if not pr:
            continue
        n_ret = sum(1 for r in pr if r["returned"])
        ret_steps = [r["return_step"] for r in pr if r["returned"]]
        mean_tau = np.mean(ret_steps) if ret_steps else -1
        print(f"  {pert:>10} {len(pr):8d} {n_ret/len(pr):8.1%} {mean_tau:7.1f}")
    
    # 4c: Candidate D variables — family discrimination
    print("\n--- Candidate D variables: family means ---")
    print(f"  {'family':6s} {'RG_dF':>8} {'RG_dφ':>8} {'ΔH_int':>8} {'W1':>8}")
    for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR"]:
        fam_dc = [d for d in d_candidates if d["family"] == fam]
        if not fam_dc:
            continue
        print(f"  {fam:6s} "
              f"{np.mean([d['rg_drift_F'] for d in fam_dc]):8.2f} "
              f"{np.mean([d['rg_drift_phi'] for d in fam_dc]):8.4f} "
              f"{np.mean([d['delta_H_int'] for d in fam_dc]):8.3f} "
              f"{np.mean([d['wasserstein'] for d in fam_dc]):8.4f}")
    
    # 4d: Kruskal-Wallis tests for family discrimination
    print("\n--- Kruskal-Wallis family discrimination ---")
    for var_name in ["rg_drift_F", "rg_drift_phi", "delta_H_int", "wasserstein"]:
        groups = []
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR"]:
            vals = [d[var_name] for d in d_candidates if d["family"] == fam]
            if vals:
                groups.append(vals)
        if len(groups) >= 2:
            H_stat, p_val = sp_stats.kruskal(*groups)
            sig = "★★★" if p_val < 0.001 else ("★★" if p_val < 0.01 else ("★" if p_val < 0.05 else "ns"))
            print(f"  {var_name:>15s}: H={H_stat:7.2f}, p={p_val:.4f} {sig}")
    
    # 4e: Correlation with F7 (independence check)
    print("\n--- Correlation with F7 (independence check) ---")
    F7_vals = [d["F7"] for d in d_candidates]
    for var_name in ["rg_drift_F", "rg_drift_phi", "delta_H_int", "wasserstein"]:
        var_vals = [d[var_name] for d in d_candidates]
        rho, p = sp_stats.spearmanr(F7_vals, var_vals)
        independent = "independent" if abs(rho) < 0.3 else "correlated"
        print(f"  {var_name:>15s}: ρ={rho:+.3f}, p={p:.4f} [{independent}]")
    
    # 4f: Recovery probability correlation with candidate D vars
    print("\n--- Recovery probability vs candidate D variables ---")
    # Merge recovery results with D candidates
    # For each sample, compute mean recovery rate
    sample_recovery = {}
    for r in recovery_results:
        key = (r["family"], r["N"], r["rep"])
        if key not in sample_recovery:
            sample_recovery[key] = {"returned": [], "basin": []}
        sample_recovery[key]["returned"].append(r["returned"])
        sample_recovery[key]["basin"].append(r["basin_retention"])
    
    merged = []
    for dc in d_candidates:
        key = (dc["family"], dc["N"], dc["rep"])
        if key in sample_recovery:
            sr = sample_recovery[key]
            merged.append({
                **dc,
                "P_return": np.mean(sr["returned"]),
                "P_basin": np.mean(sr["basin"]),
            })
    
    if merged:
        for var_name in ["rg_drift_F", "rg_drift_phi", "delta_H_int", "wasserstein"]:
            var_vals = [m[var_name] for m in merged]
            p_ret = [m["P_return"] for m in merged]
            p_basin = [m["P_basin"] for m in merged]
            
            rho_ret, p_ret_p = sp_stats.spearmanr(var_vals, p_ret)
            rho_bas, p_bas_p = sp_stats.spearmanr(var_vals, p_basin)
            
            print(f"  {var_name:>15s}: ρ(P_return)={rho_ret:+.3f} (p={p_ret_p:.3f}), "
                  f"ρ(P_basin)={rho_bas:+.3f} (p={p_bas_p:.3f})")
    
    # Save results
    outdir = Path("outputs_d_recovery")
    outdir.mkdir(exist_ok=True)
    
    # Save recovery results
    if recovery_results:
        keys = [k for k in recovery_results[0].keys()]
        with open(outdir / "recovery_results.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(recovery_results)
        print(f"\n  → {outdir / 'recovery_results.csv'} ({len(recovery_results)} rows)")
    
    # Save D candidates
    if d_candidates:
        keys = [k for k in d_candidates[0].keys()]
        with open(outdir / "d_candidates.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(d_candidates)
        print(f"  → {outdir / 'd_candidates.csv'} ({len(d_candidates)} rows)")
    
    print("\n" + "=" * 90)
    print("EXPERIMENT COMPLETE")
    print("=" * 90)


if __name__ == "__main__":
    main()
