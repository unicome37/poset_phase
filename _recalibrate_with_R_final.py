"""Final recalibration: F6 = F5 + α_R·R with R computed for ALL families.

Key fix from diagnostic:
- Sign must be POSITIVE: F6 = F5 + α_R·R (penalize high-R low-d posets)
- R must be computed for KR_like too (R≈0.30, not zero)
- MCMC likelihood vectorized for speed
"""
from __future__ import annotations
import csv, time
import numpy as np
from generators import (
    generate_lorentzian_like_2d, generate_lorentzian_like_3d,
    generate_lorentzian_like_4d, generate_lorentzian_like_5d,
    generate_kr_like,
)

# ── Load data ──
data = list(csv.DictReader(open("outputs_unified_functional/raw_features.csv")))
families_arr = np.array([d["family"] for d in data])
N_arr = np.array([int(d["N"]) for d in data], dtype=float)
log_H = np.array([float(d["log_H"]) for d in data])
pi_geo = np.array([float(d["pi_geo"]) for d in data])
sigma_hist = np.array([float(d["sigma_hist"]) for d in data])
xi_dim = np.array([float(d["xi_dim"]) for d in data])
pi_cg = np.array([float(d["pi_cg"]) for d in data])
n_data = len(data)
unique_N = sorted(set(N_arr))

# ── Compute R for ALL posets ──
print(f"Computing R for all {n_data} posets (including KR_like)...")
t0 = time.time()
gens = {
    "Lor2D": generate_lorentzian_like_2d, "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d, "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}
R_arr = np.zeros(n_data)
for i, row in enumerate(data):
    fam = row["family"]
    n, rep = int(row["N"]), int(row["rep"])
    seed = 42 + rep * 1000 + n * 100
    p = gens[fam](n, seed=seed)
    c = p.closure.astype(np.int32)
    ks = c @ c; mask = p.closure
    C0 = int(np.sum(mask & (ks == 0))); total = int(np.sum(mask))
    R_arr[i] = 1.0 - C0 / total if total > 0 else 0.0
print(f"  Done in {time.time()-t0:.1f}s")

print("\n  R by family:")
for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
    idx = families_arr == fam
    v = R_arr[idx]
    print(f"    {fam:>7}: mean={np.mean(v):.4f} ± {np.std(v):.4f}  [{np.min(v):.4f}, {np.max(v):.4f}]")

# ── Build pair index arrays (vectorized) ──
def build_pairs(fam_i, fam_j, min_N=0):
    pairs = []
    for nv in unique_N:
        if nv < min_N: continue
        ii = np.where((N_arr == nv) & (families_arr == fam_i))[0]
        jj = np.where((N_arr == nv) & (families_arr == fam_j))[0]
        for a in ii:
            for b in jj:
                pairs.append((a, b))
    return np.array(pairs, dtype=int) if pairs else np.zeros((0, 2), dtype=int)

# Pred A: 4D < other Lor families
A_pairs_2d = build_pairs("Lor4D", "Lor2D")
A_pairs_3d = build_pairs("Lor4D", "Lor3D")
A_pairs_5d = build_pairs("Lor4D", "Lor5D")
A_all = np.vstack([A_pairs_2d, A_pairs_3d, A_pairs_5d]) if len(A_pairs_2d) else np.zeros((0,2),dtype=int)

# Pred B: Lor2D < KR_like
B_pairs = build_pairs("Lor2D", "KR_like")

print(f"\n  Pairs: A_2d={len(A_pairs_2d)}, A_3d={len(A_pairs_3d)}, A_5d={len(A_pairs_5d)}, B={len(B_pairs)}")

# ── F5 features ──
f5_feat = np.column_stack([pi_geo, sigma_hist, xi_dim, pi_cg])

def compute_F5(w4):
    return log_H + f5_feat @ w4

def hit_rate(F, pairs):
    if len(pairs) == 0: return 0, 0
    wins = int(np.sum(F[pairs[:, 0]] < F[pairs[:, 1]]))
    return wins, len(pairs)

def all_hit_rates(F):
    w2, t2 = hit_rate(F, A_pairs_2d)
    w3, t3 = hit_rate(F, A_pairs_3d)
    w5, t5 = hit_rate(F, A_pairs_5d)
    wb, tb = hit_rate(F, B_pairs)
    wa = w2 + w3 + w5; ta = t2 + t3 + t5
    return wa, ta, w2, t2, w3, t3, w5, t5, wb, tb

def print_row(label, F):
    wa, ta, w2, t2, w3, t3, w5, t5, wb, tb = all_hit_rates(F)
    pa = 100*wa/ta if ta else 0
    p2 = 100*w2/t2 if t2 else 0
    p3 = 100*w3/t3 if t3 else 0
    p5 = 100*w5/t5 if t5 else 0
    pb = 100*wb/tb if tb else 0
    print(f"  {label:>28} | {wa:3d}/{ta:3d} ({pa:5.1f}%) | {p2:5.1f}% {p3:5.1f}% {p5:5.1f}% | {wb:3d}/{tb:3d} ({pb:5.1f}%)")

# ── Phase 1: Baselines ──
print("\n" + "=" * 80)
print("PHASE 1: Baselines")
print("=" * 80)
print(f"\n  {'Model':>28} | {'Pred A':>12} | {'vs2D':>6} {'vs3D':>6} {'vs5D':>6} | {'Pred B':>12}")
print("  " + "─" * 80)

w5_bay = np.array([0.0004, -0.888, 0.637, 0.068])
F5_bay = compute_F5(w5_bay)
print_row("F5 (Bayesian)", F5_bay)

w5_cal = np.array([0.5, -1.5, 0.1, 0.05])
F5_cal = compute_F5(w5_cal)
print_row("F5 (Calibrated)", F5_cal)

# ── Phase 2: α_R grid scan (with KR_like R) ──
print("\n" + "=" * 80)
print("PHASE 2: Grid scan F6 = F5_Bay + α_R·R (R computed for ALL families)")
print("=" * 80)
print(f"\n  {'α_R':>6} | {'Pred A':>12} | {'vs2D':>6} {'vs3D':>6} {'vs5D':>6} | {'Pred B':>12}")
print("  " + "─" * 70)

best_a = 0; best_rate = 0
for aR in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 80, 100]:
    F6 = F5_bay + aR * R_arr
    wa, ta, w2, t2, w3, t3, w5, t5, wb, tb = all_hit_rates(F6)
    pa = 100*wa/ta if ta else 0
    p2 = 100*w2/t2 if t2 else 0; p3 = 100*w3/t3 if t3 else 0; p5 = 100*w5/t5 if t5 else 0
    pb = 100*wb/tb if tb else 0
    print(f"  {aR:6.0f} | {wa:3d}/{ta:3d} ({pa:5.1f}%) | {p2:5.1f}% {p3:5.1f}% {p5:5.1f}% | {wb:3d}/{tb:3d} ({pb:5.1f}%)")
    if pa > best_rate: best_rate = pa; best_a = aR

print(f"\n  Coarse best: α_R = {best_a} → Pred A = {best_rate:.1f}%")

# Fine scan
print(f"\n  Fine scan [{max(0,best_a-15)}, {best_a+15}]:")
best_af = best_a; best_rf = 0
for aR in range(max(0, best_a - 15), best_a + 16):
    F6 = F5_bay + aR * R_arr
    wa, ta, w2, t2, w3, t3, w5, t5, wb, tb = all_hit_rates(F6)
    pa = 100*wa/ta if ta else 0
    p2 = 100*w2/t2 if t2 else 0; p3 = 100*w3/t3 if t3 else 0; p5 = 100*w5/t5 if t5 else 0
    pb = 100*wb/tb if tb else 0
    tag = " ***" if pa > best_rf else ""
    print(f"    {aR:4d} | {wa:3d}/{ta:3d} ({pa:5.1f}%) | {p2:5.1f}% {p3:5.1f}% {p5:5.1f}% | {wb:3d}/{tb:3d} ({pb:5.1f}%){tag}")
    if pa > best_rf: best_rf = pa; best_af = aR

print(f"\n  Fine best: α_R = {best_af} → Pred A = {best_rf:.1f}%")

# ── Phase 3: Vectorized MCMC ──
print("\n" + "=" * 80)
print("PHASE 3: MH-MCMC for F6 = log_H + γ·Π_geo - λ·Σ_hist + η·Ξ_d + κ·Π_cg + α_R·R")
print("=" * 80)

f6_feat = np.column_stack([pi_geo, sigma_hist, xi_dim, pi_cg, R_arr])

C_groups_idx = {}
for nv in unique_N:
    C_groups_idx[nv] = np.where(N_arr == nv)[0]

def log_lik(w5, scale_B=2.0, scale_A=1.5, scale_C=2.0):
    F = log_H + f6_feat @ w5
    ll = 0.0
    # B: Lor2D < KR (vectorized)
    if len(B_pairs):
        diff_b = (F[B_pairs[:, 1]] - F[B_pairs[:, 0]]) / scale_B
        diff_b = np.clip(diff_b, -30, 30)
        ll += np.sum(-np.log1p(np.exp(-diff_b)))
    # A: Lor4D < others (vectorized)
    if len(A_all):
        diff_a = (F[A_all[:, 1]] - F[A_all[:, 0]]) / scale_B
        diff_a = np.clip(diff_a, -30, 30)
        ll += scale_A * np.sum(-np.log1p(np.exp(-diff_a)))
    # C: sigma_hist correlation
    for nv, idx in C_groups_idx.items():
        sh = sigma_hist[idx]; f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            ll += scale_C * (-np.corrcoef(sh, f)[0, 1])
    return ll

def log_pri(w5):
    gamma, neg_lam, eta, kappa, alpha_R = w5
    lam = -neg_lam
    if gamma < 0 or lam < 0 or eta < 0 or kappa < 0 or alpha_R < 0: return -np.inf
    if gamma > 10 or lam > 10 or eta > 5 or kappa > 5 or alpha_R > 100: return -np.inf
    lp = -0.5*((gamma-0.25)/0.5)**2 - 0.5*((lam-0.75)/1.0)**2
    lp += -0.5*((eta-0.05)/0.3)**2 - 0.5*((kappa-0.025)/0.2)**2
    lp += -0.5*((alpha_R - 30.0)/20.0)**2
    return lp

def log_post(w5):
    lp = log_pri(w5)
    if not np.isfinite(lp): return -np.inf
    ll = log_lik(w5)
    return lp + ll if np.isfinite(ll) else -np.inf

n_steps, burn_in = 100000, 20000
rng = np.random.default_rng(42)
w0 = np.array([0.25, -0.75, 0.05, 0.025, float(best_af)])
prop = np.array([0.04, 0.06, 0.015, 0.01, 3.0])
chain = np.zeros((n_steps, 5))
w_cur = w0.copy(); lp_cur = log_post(w_cur); n_acc = 0

print(f"\nRunning MCMC ({n_steps} steps, {burn_in} burn-in)...")
t0 = time.time()
for step in range(n_steps):
    w_p = w_cur + rng.normal(0, prop)
    lp_p = log_post(w_p)
    if np.log(rng.random()) < (lp_p - lp_cur):
        w_cur = w_p; lp_cur = lp_p; n_acc += 1
    chain[step] = w_cur
    if step < burn_in and step > 0 and step % 2000 == 0:
        rate = n_acc / (step + 1)
        if rate < 0.15: prop *= 0.8
        elif rate > 0.40: prop *= 1.2
    if step % 20000 == 0 and step > 0:
        print(f"  step {step}, acc={n_acc/(step+1):.3f}, lp={lp_cur:.1f}")

dt = time.time() - t0
samples = chain[burn_in::5]
print(f"\n  Done in {dt:.1f}s. Acceptance: {n_acc/n_steps:.3f}, samples: {len(samples)}")

# Display
disp = samples.copy()
disp[:, 1] = -disp[:, 1]  # -λ → λ
names = ["γ", "λ", "η", "κ", "α_R"]

print(f"\n  {'Param':>6s} {'Mean':>8s} {'Median':>8s} {'Std':>8s} {'2.5%':>8s} {'97.5%':>8s}")
print(f"  {'-'*6} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
for i, name in enumerate(names):
    s = disp[:, i]
    print(f"  {name:>6s} {np.mean(s):8.4f} {np.median(s):8.4f} "
          f"{np.std(s):8.4f} {np.percentile(s, 2.5):8.4f} {np.percentile(s, 97.5):8.4f}")

# Hit rates
w6_med = np.median(samples, axis=0)
F6_med = log_H + f6_feat @ w6_med

print("\n" + "=" * 80)
print("FINAL COMPARISON")
print("=" * 80)
print(f"\n  {'Model':>28} | {'Pred A':>12} | {'vs2D':>6} {'vs3D':>6} {'vs5D':>6} | {'Pred B':>12}")
print("  " + "─" * 80)

print_row("F5 (Bayesian)", F5_bay)
print_row("F5 (Calibrated)", F5_cal)
print_row(f"F6 (scan aR={best_af})", F5_bay + best_af * R_arr)
print_row("F6 (MCMC median)", F6_med)

# N-breakdown for MCMC median
print(f"\n  N-breakdown for F6 (MCMC median):")
for nv in unique_N:
    for fo in ["Lor2D", "Lor3D", "Lor5D"]:
        i4 = np.where((N_arr == nv) & (families_arr == "Lor4D"))[0]
        io = np.where((N_arr == nv) & (families_arr == fo))[0]
        if len(i4) == 0 or len(io) == 0: continue
        wins = int(np.sum(F6_med[i4, None] < F6_med[None, io]))
        total = len(i4) * len(io)
        print(f"    N={int(nv):2d}, 4D vs {fo:>5}: {wins:2d}/{total:2d} ({100*wins/total:5.1f}%)")

# Pred B breakdown
print(f"\n  Pred B (Lor2D < KR_like) N-breakdown:")
for nv in unique_N:
    i2 = np.where((N_arr == nv) & (families_arr == "Lor2D"))[0]
    ik = np.where((N_arr == nv) & (families_arr == "KR_like"))[0]
    if len(i2) == 0 or len(ik) == 0: continue
    wins = int(np.sum(F6_med[i2, None] < F6_med[None, ik]))
    total = len(i2) * len(ik)
    print(f"    N={int(nv):2d}: {wins:2d}/{total:2d} ({100*wins/total:5.1f}%)")

# ── Phase 4: Optimal α_R that maximizes BOTH A and B ──
print("\n" + "=" * 80)
print("PHASE 4: Pareto search — maximize Pred A subject to Pred B ≥ 90%")
print("=" * 80)
print(f"\n  {'α_R':>6} | {'Pred A':>12} | {'Pred B':>12} | {'A+B':>6}")
print("  " + "─" * 50)

best_pareto = 0; best_pareto_a = 0
for aR in range(0, 101):
    F6 = F5_bay + aR * R_arr
    wa, ta, _, _, _, _, _, _, wb, tb = all_hit_rates(F6)
    pa = 100*wa/ta if ta else 0; pb = 100*wb/tb if tb else 0
    if pb >= 90 and pa > best_pareto_a:
        best_pareto = aR; best_pareto_a = pa
    if aR % 5 == 0:
        tag = " ***" if pb >= 90 else ""
        print(f"  {aR:6d} | {wa:3d}/{ta:3d} ({pa:5.1f}%) | {wb:3d}/{tb:3d} ({pb:5.1f}%) | {pa+pb:6.1f}{tag}")

print(f"\n  Pareto optimal: α_R = {best_pareto} → Pred A = {best_pareto_a:.1f}% with Pred B ≥ 90%")
if best_pareto > 0:
    print_row(f"F6 (Pareto aR={best_pareto})", F5_bay + best_pareto * R_arr)
