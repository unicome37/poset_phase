"""Recalibrate F5 → F6 with R term. Self-contained, no heavy imports."""
from __future__ import annotations
import csv, sys, time
import numpy as np
from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)

# ── Load data ──
print("Loading raw_features.csv...")
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

# ── Compute R for Lorentzian posets ──
print(f"Computing R for {n_data} posets...")
t0 = time.time()
gens = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
}
R_arr = np.zeros(n_data)
for i, row in enumerate(data):
    fam = row["family"]
    if fam not in gens:
        R_arr[i] = 0.0
        continue
    n, rep = int(row["N"]), int(row["rep"])
    seed = 42 + rep * 1000 + n * 100
    poset = gens[fam](n, seed=seed)
    c = poset.closure.astype(np.int32)
    ks = c @ c
    mask = poset.closure
    C0 = int(np.sum(mask & (ks == 0)))
    total = int(np.sum(mask))
    R_arr[i] = 1.0 - C0 / total if total > 0 else 0.0

print(f"  Done in {time.time()-t0:.1f}s")
print("\n  R by family:")
for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
    idx = families_arr == fam
    vals = R_arr[idx]
    print(f"    {fam}: mean={np.mean(vals):.4f}, std={np.std(vals):.4f}, "
          f"min={np.min(vals):.4f}, max={np.max(vals):.4f}")

# ── F5 computation ──
f5_features = np.column_stack([pi_geo, sigma_hist, xi_dim, pi_cg])

def compute_F5(w):
    return log_H + f5_features @ w

# ── Pair indices for Prediction A ──
B_pairs = []
for n_val in unique_N:
    lor_idx = np.where((N_arr == n_val) & (families_arr == "Lor2D"))[0]
    kr_idx = np.where((N_arr == n_val) & (families_arr == "KR_like"))[0]
    for i in lor_idx:
        for j in kr_idx:
            B_pairs.append((i, j))

A_all_pairs = []
for n_val in unique_N:
    lor4_idx = np.where((N_arr == n_val) & (families_arr == "Lor4D"))[0]
    for other in ["Lor2D", "Lor3D", "Lor5D"]:
        other_idx = np.where((N_arr == n_val) & (families_arr == other))[0]
        for i in lor4_idx:
            for j in other_idx:
                A_all_pairs.append((i, j))

def compute_hit_rates(F_values):
    results = {}
    for fam_other in ["Lor2D", "Lor3D", "Lor5D"]:
        comp_wins = 0; comp_total = 0
        for n_val in unique_N:
            idx_4d = np.where((N_arr == n_val) & (families_arr == "Lor4D"))[0]
            idx_o = np.where((N_arr == n_val) & (families_arr == fam_other))[0]
            for i in idx_4d:
                for j in idx_o:
                    if F_values[i] < F_values[j]:
                        comp_wins += 1
                    comp_total += 1
        results[fam_other] = (comp_wins, comp_total)
    total_w = sum(v[0] for v in results.values())
    total_t = sum(v[1] for v in results.values())
    return total_w, total_t, results

def compute_B_rate(F_values):
    wins = sum(1 for i, j in B_pairs if F_values[i] < F_values[j])
    return wins, len(B_pairs)

# ── Phase 1: F5 baselines ──
print("\n" + "=" * 70)
print("PHASE 1: F5 baselines")
print("=" * 70)

# Bayesian weights: γ=0.0004, λ=0.888, η=0.637, κ=0.068
w5_bay = np.array([0.0004, -0.888, 0.637, 0.068])
F5_bay = compute_F5(w5_bay)
w5a, t5a, by5a = compute_hit_rates(F5_bay)
bw, bt = compute_B_rate(F5_bay)
print(f"\n  F5 (Bayesian): Pred A = {w5a}/{t5a} ({100*w5a/t5a:.1f}%), Pred B = {bw}/{bt} ({100*bw/bt:.1f}%)")
for fam, (w, t) in by5a.items():
    print(f"    4D vs {fam}: {w}/{t} ({100*w/t:.1f}%)")

# Calibrated weights: γ=0.5, λ=1.5, η=0.1, κ=0.05
w5_cal = np.array([0.5, -1.5, 0.1, 0.05])
F5_cal = compute_F5(w5_cal)
w5c, t5c, by5c = compute_hit_rates(F5_cal)
bwc, btc = compute_B_rate(F5_cal)
print(f"\n  F5 (Calibrated): Pred A = {w5c}/{t5c} ({100*w5c/t5c:.1f}%), Pred B = {bwc}/{btc} ({100*bwc/btc:.1f}%)")
for fam, (w, t) in by5c.items():
    print(f"    4D vs {fam}: {w}/{t} ({100*w/t:.1f}%)")

# ── Phase 2: α_R grid scan ──
print("\n" + "=" * 70)
print("PHASE 2: α_R grid scan (F5_Bayesian + α_R·R)")
print("=" * 70)

print(f"\n  {'α_R':>6} | {'Pred A':>10} | {'vs 2D':>10} | {'vs 3D':>10} | {'vs 5D':>10} | {'Pred B':>10}")
print("  " + "─" * 75)

best_alpha = 0; best_rate = 0
scan_results = []

for alpha_R in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 80, 100]:
    F6 = F5_bay + alpha_R * R_arr
    wa, ta, bya = compute_hit_rates(F6)
    wb, tb = compute_B_rate(F6)
    rate_a = wa / ta if ta else 0
    rate_b = wb / tb if tb else 0
    
    vs2 = bya["Lor2D"]; vs3 = bya["Lor3D"]; vs5 = bya["Lor5D"]
    print(f"  {alpha_R:6.1f} | {wa:3d}/{ta:3d} ({100*rate_a:5.1f}%) | "
          f"{vs2[0]:3d}/{vs2[1]:3d} ({100*vs2[0]/vs2[1]:5.1f}%) | "
          f"{vs3[0]:3d}/{vs3[1]:3d} ({100*vs3[0]/vs3[1]:5.1f}%) | "
          f"{vs5[0]:3d}/{vs5[1]:3d} ({100*vs5[0]/vs5[1]:5.1f}%) | "
          f"{wb:3d}/{tb:3d} ({100*rate_b:5.1f}%)")
    
    scan_results.append((alpha_R, rate_a, rate_b))
    if rate_a > best_rate:
        best_rate = rate_a; best_alpha = alpha_R

print(f"\n  Best α_R = {best_alpha:.1f} → Pred A = {100*best_rate:.1f}%")

# ── Phase 3: Fine scan around best ──
print("\n" + "=" * 70)
print(f"PHASE 3: Fine scan around α_R = {best_alpha}")
print("=" * 70)

fine_range = np.arange(max(0, best_alpha - 15), best_alpha + 16, 1)
print(f"\n  {'α_R':>6} | {'Pred A':>10} | {'vs 2D':>10} | {'vs 3D':>10} | {'vs 5D':>10} | {'Pred B':>10}")
print("  " + "─" * 75)

best_alpha_fine = best_alpha; best_rate_fine = 0

for alpha_R in fine_range:
    F6 = F5_bay + alpha_R * R_arr
    wa, ta, bya = compute_hit_rates(F6)
    wb, tb = compute_B_rate(F6)
    rate_a = wa / ta if ta else 0
    rate_b = wb / tb if tb else 0
    
    vs2 = bya["Lor2D"]; vs3 = bya["Lor3D"]; vs5 = bya["Lor5D"]
    print(f"  {alpha_R:6.1f} | {wa:3d}/{ta:3d} ({100*rate_a:5.1f}%) | "
          f"{vs2[0]:3d}/{vs2[1]:3d} ({100*vs2[0]/vs2[1]:5.1f}%) | "
          f"{vs3[0]:3d}/{vs3[1]:3d} ({100*vs3[0]/vs3[1]:5.1f}%) | "
          f"{vs5[0]:3d}/{vs5[1]:3d} ({100*vs5[0]/vs5[1]:5.1f}%) | "
          f"{wb:3d}/{tb:3d} ({100*rate_b:5.1f}%)")
    
    if rate_a > best_rate_fine:
        best_rate_fine = rate_a; best_alpha_fine = alpha_R

print(f"\n  Best α_R = {best_alpha_fine:.1f} → Pred A = {100*best_rate_fine:.1f}%")

# ── Phase 4: Full MCMC with R term ──
print("\n" + "=" * 70)
print("PHASE 4: MH-MCMC for F6 = log_H + γ·Π_geo - λ·Σ_hist + η·Ξ_d + κ·Π_cg + α_R·R")
print("=" * 70)

f6_features = np.column_stack([pi_geo, sigma_hist, xi_dim, pi_cg, R_arr])

def compute_F6_mcmc(w):
    """w = [γ, -λ, η, κ, +α_R] → F6 = log_H + γ·Π_geo - λ·Σ_hist + η·Ξ_d + κ·Π_cg + α_R·R"""
    return log_H + f6_features @ w

C_groups = {}
for n_val in unique_N:
    C_groups[n_val] = np.where(N_arr == n_val)[0]

def log_likelihood_f6(w, scale_B=2.0, scale_A=1.5, scale_C=2.0):
    F = compute_F6_mcmc(w)
    ll = 0.0
    for i, j in B_pairs:
        diff = np.clip((F[j] - F[i]) / scale_B, -30, 30)
        ll += -np.log1p(np.exp(-diff))
    for i, j in A_all_pairs:
        diff = np.clip((F[j] - F[i]) / scale_B, -30, 30)
        ll += scale_A * (-np.log1p(np.exp(-diff)))
    for n_val, idx in C_groups.items():
        sh = sigma_hist[idx]; f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            ll += scale_C * (-np.corrcoef(sh, f)[0, 1])
    return ll

def log_prior_f6(w):
    gamma, neg_lam, eta, kappa, alpha_R = w
    lam = -neg_lam
    if gamma < 0 or lam < 0 or eta < 0 or kappa < 0 or alpha_R < 0: return -np.inf
    if gamma > 10 or lam > 10 or eta > 5 or kappa > 5 or alpha_R > 100: return -np.inf
    lp = -0.5*((gamma-0.25)/0.5)**2 - 0.5*((lam-0.75)/1.0)**2
    lp += -0.5*((eta-0.05)/0.3)**2 - 0.5*((kappa-0.025)/0.2)**2
    lp += -0.5*((alpha_R - 30.0)/20.0)**2
    return lp

def log_post_f6(w):
    lp = log_prior_f6(w)
    if not np.isfinite(lp): return -np.inf
    ll = log_likelihood_f6(w)
    return lp + ll if np.isfinite(ll) else -np.inf

print("\nRunning MCMC (80000 steps, 20000 burn-in)...")
rng = np.random.default_rng(42)
n_steps, burn_in = 80000, 20000
w0 = np.array([0.25, -0.75, 0.05, 0.025, 30.0])
prop_scale = np.array([0.04, 0.06, 0.015, 0.01, 3.0])

chain = np.zeros((n_steps, 5))
w_cur = w0.copy(); lp_cur = log_post_f6(w_cur); n_acc = 0

for step in range(n_steps):
    w_prop = w_cur + rng.normal(0, prop_scale)
    lp_prop = log_post_f6(w_prop)
    if np.log(rng.random()) < (lp_prop - lp_cur):
        w_cur = w_prop; lp_cur = lp_prop; n_acc += 1
    chain[step] = w_cur
    if step < burn_in and step > 0 and step % 1000 == 0:
        rate = n_acc / (step + 1)
        if rate < 0.15: prop_scale *= 0.8
        elif rate > 0.40: prop_scale *= 1.2

samples = chain[burn_in::5]
print(f"Acceptance rate: {n_acc/n_steps:.3f}, samples: {len(samples)}")

disp = samples.copy()
disp[:, 1] = -disp[:, 1]  # -λ → λ
# α_R is already positive, no sign flip needed
names = ["γ", "λ", "η", "κ", "α_R"]

print(f"\n  {'Param':>6s} {'Mean':>8s} {'Median':>8s} {'Std':>8s} {'2.5%':>8s} {'97.5%':>8s}")
print(f"  {'-'*6} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
for i, name in enumerate(names):
    s = disp[:, i]
    print(f"  {name:>6s} {np.mean(s):8.4f} {np.median(s):8.4f} "
          f"{np.std(s):8.4f} {np.percentile(s, 2.5):8.4f} {np.percentile(s, 97.5):8.4f}")

# F6 with posterior median
w6_med = np.median(samples, axis=0)
F6_med = compute_F6_mcmc(w6_med)
w6a, t6a, by6a = compute_hit_rates(F6_med)
bw6, bt6 = compute_B_rate(F6_med)
print(f"\n  F6 (MCMC median): Pred A = {w6a}/{t6a} ({100*w6a/t6a:.1f}%), Pred B = {bw6}/{bt6} ({100*bw6/bt6:.1f}%)")
for fam, (w, t) in by6a.items():
    print(f"    4D vs {fam}: {w}/{t} ({100*w/t:.1f}%)")

# ── Phase 5: Summary ──
print("\n" + "=" * 70)
print("SUMMARY: F5 vs F6")
print("=" * 70)

# Best grid-scan F6
F6_best = F5_bay + best_alpha_fine * R_arr
w6b, t6b, by6b = compute_hit_rates(F6_best)
bwb, btb = compute_B_rate(F6_best)

print(f"\n  {'Model':>25} | {'Pred A':>12} | {'vs2D':>8} {'vs3D':>8} {'vs5D':>8} | {'Pred B':>10}")
print("  " + "─" * 85)

for label, F_vals in [
    ("F5 (Bayesian)", F5_bay),
    ("F5 (Calibrated)", F5_cal),
    (f"F6 (scan aR={best_alpha_fine:.0f})", F6_best),
    ("F6 (MCMC median)", F6_med),
]:
    wa, ta, bya = compute_hit_rates(F_vals)
    wb, tb = compute_B_rate(F_vals)
    vs2=bya["Lor2D"]; vs3=bya["Lor3D"]; vs5=bya["Lor5D"]
    print(f"  {label:>25} | {wa:3d}/{ta:3d} ({100*wa/ta:5.1f}%) | "
          f"{100*vs2[0]/vs2[1]:5.1f}% {100*vs3[0]/vs3[1]:5.1f}% {100*vs5[0]/vs5[1]:5.1f}% | "
          f"{wb:3d}/{tb:3d} ({100*wb/tb:5.1f}%)")

# N-breakdown for best model
print(f"\n  N-breakdown for F6 (MCMC median):")
for n_val in unique_N:
    for fam_other in ["Lor2D", "Lor3D", "Lor5D"]:
        idx_4d = np.where((N_arr == n_val) & (families_arr == "Lor4D"))[0]
        idx_o = np.where((N_arr == n_val) & (families_arr == fam_other))[0]
        wins = sum(1 for i in idx_4d for j in idx_o if F6_med[i] < F6_med[j])
        total = len(idx_4d) * len(idx_o)
        pct = 100 * wins / total if total else 0
        print(f"    N={int(n_val):2d}, 4D vs {fam_other}: {wins}/{total} ({pct:.0f}%)")
