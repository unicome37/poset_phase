"""Recalibrate unified functional weights WITH the R (interval richness) term.

Current F5 = β·log_H + γ·Π_geo - λ·Σ_hist + η·Ξ_d + κ·Π_cg

New F6 = F5 - α_R · R

where R = 1 - C0/ΣCk (non-link fraction = interval occupancy).

The sign is NEGATIVE because:
- Lower F is "better" (more Lorentzian-like)
- Higher R means richer interval structure (more mediated causality)
- For 4D, R is moderate; for 2D, R is very high
- So -α_R·R penalizes low-d (high R) posets → pushes 4D to be minimum

This script:
1. Computes R for all 160 posets in raw_features.csv
2. Runs MH-MCMC for F6 weights (β≡1, γ, λ, η, κ, α_R)
3. Compares F5 vs F6 hit rates for Prediction A
"""
from __future__ import annotations

import csv
import math
import numpy as np
from generators import (
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_kr_like,
    Poset,
)
from prediction_a_bd_bridge import (
    count_intervals_fast,
    regenerate_poset,
    BAYESIAN_WEIGHTS,
    CALIBRATED_WEIGHTS,
)


# ── Compute R for a poset ──

def compute_R(poset: Poset) -> float:
    """R = 1 - C0/ΣCk (interval occupancy / non-link fraction)."""
    counts = count_intervals_fast(poset)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    C0 = counts.get(0, 0)
    return 1.0 - C0 / total


# ── Load data + compute R ──

print("Loading raw_features.csv and computing R for all posets...")
data = []
with open("outputs_unified_functional/raw_features.csv", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        data.append(row)

families_arr = np.array([d["family"] for d in data])
N_arr = np.array([int(d["N"]) for d in data], dtype=float)
log_H = np.array([float(d["log_H"]) for d in data])
pi_geo = np.array([float(d["pi_geo"]) for d in data])
sigma_hist = np.array([float(d["sigma_hist"]) for d in data])
xi_dim = np.array([float(d["xi_dim"]) for d in data])
pi_cg = np.array([float(d["pi_cg"]) for d in data])

n_data = len(data)
unique_N = sorted(set(N_arr))

# Compute R for each poset
generators_map = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}

R_arr = np.zeros(n_data)
lor_set = {"Lor2D", "Lor3D", "Lor4D", "Lor5D"}
computed = 0
for i, row in enumerate(data):
    fam = row["family"]
    if fam not in lor_set:
        R_arr[i] = 0.0
        continue
    n = int(row["N"])
    rep = int(row["rep"])
    seed = 42 + rep * 1000 + n * 100
    gen = generators_map[fam]
    poset = gen(n, seed=seed)
    # Fast R computation via matrix ops
    c = poset.closure.astype(np.int32)
    interval_sizes = c @ c
    mask = poset.closure
    C0 = 0; total = 0
    for ii in range(n):
        for jj in range(n):
            if mask[ii, jj]:
                total += 1
                if interval_sizes[ii, jj] == 0:
                    C0 += 1
    R_arr[i] = 1.0 - C0 / total if total > 0 else 0.0
    computed += 1

print(f"  Done. R computed for {n_data} posets.")

# Show R by family
print("\n  R by family:")
for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
    idx = families_arr == fam
    print(f"    {fam}: mean={np.mean(R_arr[idx]):.4f}, std={np.std(R_arr[idx]):.4f}")


# ── F5 and F6 computation ──

# Features: [γ, -λ, η, κ]  → F5 = log_H + γ·pi_geo - λ·sigma_hist + η·xi_dim + κ·pi_cg
# For F6: [γ, -λ, η, κ, -α_R] → F6 = F5 - α_R·R
f5_features = np.column_stack([pi_geo, sigma_hist, xi_dim, pi_cg])
f6_features = np.column_stack([pi_geo, sigma_hist, xi_dim, pi_cg, R_arr])


def compute_F5(w):
    """w = [γ, -λ, η, κ]"""
    return log_H + f5_features @ w


def compute_F6(w):
    """w = [γ, -λ, η, κ, -α_R]"""
    return log_H + f6_features @ w


# ── Precompute pair indices ──

B_pairs = []
for n_val in unique_N:
    lor_idx = np.where((N_arr == n_val) & (families_arr == "Lor2D"))[0]
    kr_idx = np.where((N_arr == n_val) & (families_arr == "KR_like"))[0]
    for i in lor_idx:
        for j in kr_idx:
            B_pairs.append((i, j))

A_pairs_strong = []  # 4D vs 5D (strong signal)
A_pairs_weak = []    # 4D vs 2D/3D (currently weak, R should help)
for n_val in unique_N:
    if n_val < 28:
        continue
    lor4_idx = np.where((N_arr == n_val) & (families_arr == "Lor4D"))[0]
    for other in ["Lor5D"]:
        other_idx = np.where((N_arr == n_val) & (families_arr == other))[0]
        for i in lor4_idx:
            for j in other_idx:
                A_pairs_strong.append((i, j))
    for other in ["Lor2D", "Lor3D"]:
        other_idx = np.where((N_arr == n_val) & (families_arr == other))[0]
        for i in lor4_idx:
            for j in other_idx:
                A_pairs_weak.append((i, j))

C_groups = {}
for n_val in unique_N:
    C_groups[n_val] = np.where(N_arr == n_val)[0]

print(f"\nPair counts: B={len(B_pairs)}, A_strong(4v5)={len(A_pairs_strong)}, "
      f"A_weak(4v2,3)={len(A_pairs_weak)}, C_groups={len(C_groups)}")


# ── Log-likelihood for F6 ──

def log_likelihood_f6(w, scale_B=2.0, scale_A_strong=2.0, scale_A_weak=1.0, scale_C=2.0):
    F = compute_F6(w)
    ll = 0.0
    
    # B: Lor2D < KR (Lorentzian-like should be lower)
    for i, j in B_pairs:
        diff = (F[j] - F[i]) / scale_B
        diff = np.clip(diff, -30, 30)
        ll += -np.log1p(np.exp(-diff))
    
    # A_strong: Lor4D < Lor5D
    for i, j in A_pairs_strong:
        diff = (F[j] - F[i]) / scale_B
        diff = np.clip(diff, -30, 30)
        ll += scale_A_strong * (-np.log1p(np.exp(-diff)))
    
    # A_weak: Lor4D < Lor2D/3D (this is the key improvement target)
    for i, j in A_pairs_weak:
        diff = (F[j] - F[i]) / scale_B
        diff = np.clip(diff, -30, 30)
        ll += scale_A_weak * (-np.log1p(np.exp(-diff)))
    
    # C: within-N correlation of Σ_hist with F should be negative
    for n_val, idx in C_groups.items():
        sh = sigma_hist[idx]
        f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            r = np.corrcoef(sh, f)[0, 1]
            ll += scale_C * (-r)
    
    return ll


def log_prior_f6(w):
    """Prior on [γ, -λ, η, κ, -α_R]."""
    gamma, neg_lam, eta, kappa, neg_alpha_R = w
    lam = -neg_lam
    alpha_R = -neg_alpha_R
    
    if gamma < 0 or lam < 0 or eta < 0 or kappa < 0 or alpha_R < 0:
        return -np.inf
    if gamma > 10 or lam > 10 or eta > 5 or kappa > 5 or alpha_R > 20:
        return -np.inf
    
    lp = 0.0
    lp += -0.5 * ((gamma - 0.25) / 0.5)**2
    lp += -0.5 * ((lam - 0.75) / 1.0)**2
    lp += -0.5 * ((eta - 0.05) / 0.3)**2
    lp += -0.5 * ((kappa - 0.025) / 0.2)**2
    # Weak prior on α_R: centered at 2, wide
    lp += -0.5 * ((alpha_R - 2.0) / 3.0)**2
    return lp


def log_posterior_f6(w):
    lp = log_prior_f6(w)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood_f6(w)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


# ── MCMC for F6 ──

def run_mcmc_f6(n_steps=200000, burn_in=40000, seed=42):
    rng = np.random.default_rng(seed)
    
    # Start: [γ, -λ, η, κ, -α_R]
    w0 = np.array([0.25, -0.75, 0.05, 0.025, -2.0])
    proposal_scale = np.array([0.04, 0.06, 0.015, 0.01, 0.3])
    
    chain = np.zeros((n_steps, 5))
    log_posts = np.zeros(n_steps)
    w_current = w0.copy()
    lp_current = log_posterior_f6(w_current)
    n_accept = 0
    
    for step in range(n_steps):
        w_prop = w_current + rng.normal(0, proposal_scale)
        lp_prop = log_posterior_f6(w_prop)
        
        if np.log(rng.random()) < (lp_prop - lp_current):
            w_current = w_prop
            lp_current = lp_prop
            n_accept += 1
        
        chain[step] = w_current
        log_posts[step] = lp_current
        
        if step < burn_in and step > 0 and step % 1000 == 0:
            rate = n_accept / (step + 1)
            if rate < 0.15:
                proposal_scale *= 0.8
            elif rate > 0.40:
                proposal_scale *= 1.2
    
    return chain[burn_in:], n_accept / n_steps


# ── Prediction A hit rate computation ──

def compute_hit_rates(F_values, label=""):
    """Compute pairwise hit rates for Prediction A."""
    lor_families = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]
    
    results = {}
    for n_val in unique_N:
        for fam_other in ["Lor2D", "Lor3D", "Lor5D"]:
            idx_4d = np.where((N_arr == n_val) & (families_arr == "Lor4D"))[0]
            idx_other = np.where((N_arr == n_val) & (families_arr == fam_other))[0]
            
            wins = 0
            total = 0
            for i in idx_4d:
                for j in idx_other:
                    if F_values[i] < F_values[j]:
                        wins += 1
                    total += 1
            
            key = f"4D_vs_{fam_other}_N{int(n_val)}"
            results[key] = (wins, total) if total > 0 else (0, 0)
    
    # Aggregate
    total_wins = 0
    total_pairs = 0
    by_comparison = {}
    for fam_other in ["Lor2D", "Lor3D", "Lor5D"]:
        comp_wins = 0
        comp_total = 0
        for n_val in unique_N:
            key = f"4D_vs_{fam_other}_N{int(n_val)}"
            w, t = results.get(key, (0, 0))
            comp_wins += w
            comp_total += t
        by_comparison[fam_other] = (comp_wins, comp_total)
        total_wins += comp_wins
        total_pairs += comp_total
    
    return total_wins, total_pairs, by_comparison


# ── Run ──

print("\n" + "=" * 70)
print("PHASE 1: F5 baseline (no R term)")
print("=" * 70)

# F5 with Bayesian weights
w5_bayes = np.array([0.0004, -0.888, 0.637, 0.068])
F5_bayes = compute_F5(w5_bayes)
wins5, total5, by5 = compute_hit_rates(F5_bayes, "F5_Bayesian")

print(f"\n  F5 (Bayesian weights): {wins5}/{total5} = {100*wins5/total5:.1f}%")
for fam, (w, t) in by5.items():
    pct = 100 * w / t if t > 0 else 0
    print(f"    4D vs {fam}: {w}/{t} = {pct:.1f}%")

# F5 with calibrated weights
w5_cal = np.array([0.5, -1.5, 0.1, 0.05])
F5_cal = compute_F5(w5_cal)
wins5c, total5c, by5c = compute_hit_rates(F5_cal, "F5_Calibrated")

print(f"\n  F5 (Calibrated weights): {wins5c}/{total5c} = {100*wins5c/total5c:.1f}%")
for fam, (w, t) in by5c.items():
    pct = 100 * w / t if t > 0 else 0
    print(f"    4D vs {fam}: {w}/{t} = {pct:.1f}%")

print("\n" + "=" * 70)
print("PHASE 2: Quick α_R scan (fix F5 weights, vary α_R)")
print("=" * 70)

# Before full MCMC, do a grid scan of α_R with fixed F5 weights
print(f"\n  F5 weights fixed at Bayesian: γ={0.0004}, λ={0.888}, η={0.637}, κ={0.068}")
print(f"  Scanning α_R from 0 to 20...")

print(f"\n  {'α_R':>6} | {'Total':>8} | {'vs 2D':>8} | {'vs 3D':>8} | {'vs 5D':>8}")
print("  " + "─" * 55)

best_alpha = 0
best_rate = 0

for alpha_R in np.arange(0, 21, 1):
    F6_scan = F5_bayes - alpha_R * R_arr
    wins, total, by_comp = compute_hit_rates(F6_scan)
    rate = wins / total if total > 0 else 0
    
    vs2 = by_comp["Lor2D"]
    vs3 = by_comp["Lor3D"]
    vs5 = by_comp["Lor5D"]
    
    print(f"  {alpha_R:6.1f} | {wins:3d}/{total:3d} ({100*rate:5.1f}%) | "
          f"{vs2[0]:3d}/{vs2[1]:3d} ({100*vs2[0]/vs2[1] if vs2[1] else 0:5.1f}%) | "
          f"{vs3[0]:3d}/{vs3[1]:3d} ({100*vs3[0]/vs3[1] if vs3[1] else 0:5.1f}%) | "
          f"{vs5[0]:3d}/{vs5[1]:3d} ({100*vs5[0]/vs5[1] if vs5[1] else 0:5.1f}%)")
    
    if rate > best_rate:
        best_rate = rate
        best_alpha = alpha_R

print(f"\n  Best α_R = {best_alpha:.1f} (hit rate = {100*best_rate:.1f}%)")

print("\n" + "=" * 70)
print("PHASE 3: Full MH-MCMC for F6 = F5 - α_R·R")
print("=" * 70)

print("\nRunning MCMC (200000 steps, 40000 burn-in)...")
samples_f6, acc_rate = run_mcmc_f6(n_steps=200000, burn_in=40000)
print(f"Acceptance rate: {acc_rate:.3f}")

# Thin
samples_f6 = samples_f6[::5]
print(f"Posterior samples (after thin-5): {len(samples_f6)}")

# Convert to display form
display = samples_f6.copy()
display[:, 1] = -display[:, 1]  # -λ → λ
display[:, 4] = -display[:, 4]  # -α_R → α_R

param_names = ["γ", "λ", "η", "κ", "α_R"]

print(f"\n  {'Param':>6s} {'Mean':>8s} {'Median':>8s} {'Std':>8s} {'2.5%':>8s} {'97.5%':>8s}")
print(f"  {'-'*6} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")

for i, name in enumerate(param_names):
    s = display[:, i]
    print(f"  {name:>6s} {np.mean(s):8.4f} {np.median(s):8.4f} "
          f"{np.std(s):8.4f} {np.percentile(s, 2.5):8.4f} {np.percentile(s, 97.5):8.4f}")

# Compute F6 with posterior median
w6_median = np.median(samples_f6, axis=0)
F6_median = compute_F6(w6_median)
wins6, total6, by6 = compute_hit_rates(F6_median, "F6_posterior_median")

print(f"\n  F6 (posterior median): {wins6}/{total6} = {100*wins6/total6:.1f}%")
for fam, (w, t) in by6.items():
    pct = 100 * w / t if t > 0 else 0
    print(f"    4D vs {fam}: {w}/{t} = {pct:.1f}%")

# Also compute with posterior mean
w6_mean = np.mean(samples_f6, axis=0)
F6_mean = compute_F6(w6_mean)
wins6m, total6m, by6m = compute_hit_rates(F6_mean, "F6_posterior_mean")

print(f"\n  F6 (posterior mean): {wins6m}/{total6m} = {100*wins6m/total6m:.1f}%")
for fam, (w, t) in by6m.items():
    pct = 100 * w / t if t > 0 else 0
    print(f"    4D vs {fam}: {w}/{t} = {pct:.1f}%")

print("\n" + "=" * 70)
print("COMPARISON: F5 vs F6")
print("=" * 70)

print(f"\n  {'Model':>20} | {'Total':>8} | {'vs 2D':>8} | {'vs 3D':>8} | {'vs 5D':>8}")
print("  " + "─" * 62)

for label, (w, t, by) in [
    ("F5 (Bayesian)", (wins5, total5, by5)),
    ("F5 (Calibrated)", (wins5c, total5c, by5c)),
    (f"F6 (α_R={best_alpha:.0f}, scan)", 
     compute_hit_rates(F5_bayes - best_alpha * R_arr) + ({},)),
    ("F6 (MCMC median)", (wins6, total6, by6)),
    ("F6 (MCMC mean)", (wins6m, total6m, by6m)),
]:
    if len(by) == 0:
        w_total, t_total, by_comp = compute_hit_rates(F5_bayes - best_alpha * R_arr)
    else:
        w_total, t_total, by_comp = w, t, by
    
    pct_total = 100 * w_total / t_total if t_total else 0
    vs2 = by_comp.get("Lor2D", (0, 0))
    vs3 = by_comp.get("Lor3D", (0, 0))
    vs5 = by_comp.get("Lor5D", (0, 0))
    
    p2 = 100 * vs2[0] / vs2[1] if vs2[1] else 0
    p3 = 100 * vs3[0] / vs3[1] if vs3[1] else 0
    p5 = 100 * vs5[0] / vs5[1] if vs5[1] else 0
    
    print(f"  {label:>20} | {w_total:3d}/{t_total:3d} ({pct_total:5.1f}%) | "
          f"({p2:5.1f}%) | ({p3:5.1f}%) | ({p5:5.1f}%)")

# Also check Prediction B (Lor2D < KR)
print("\n  Prediction B check (Lor2D < KR):")
for label, F_vals in [("F5_Bayesian", F5_bayes), ("F6_median", F6_median)]:
    wins_b = 0
    total_b = 0
    for i, j in B_pairs:
        if F_vals[i] < F_vals[j]:
            wins_b += 1
        total_b += 1
    print(f"    {label}: {wins_b}/{total_b} = {100*wins_b/total_b:.1f}%")


if __name__ == "__main__":
    pass
