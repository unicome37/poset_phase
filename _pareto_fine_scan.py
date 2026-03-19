"""Fine Pareto scan around balanced region + full 5D MH-MCMC with multi-objective.

Key insight from grid scan: balanced region is α_R ≈ 30-33, η ≈ 3.5-4.0.
A–B tension is structural: R_2D > R_KR → α_R helps A but hurts B.
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

# ── Compute R ──
print("Computing R...")
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

# ── Build pair arrays ──
def build_pairs(fam_i, fam_j):
    pairs = []
    for nv in unique_N:
        ii = np.where((N_arr == nv) & (families_arr == fam_i))[0]
        jj = np.where((N_arr == nv) & (families_arr == fam_j))[0]
        for a in ii:
            for b in jj:
                pairs.append((a, b))
    return np.array(pairs, dtype=int) if pairs else np.zeros((0, 2), dtype=int)

A_pairs_2d = build_pairs("Lor4D", "Lor2D")
A_pairs_3d = build_pairs("Lor4D", "Lor3D")
A_pairs_5d = build_pairs("Lor4D", "Lor5D")
A_all = np.vstack([A_pairs_2d, A_pairs_3d, A_pairs_5d])
B_pairs = build_pairs("Lor2D", "KR_like")

feat_6 = np.column_stack([pi_geo, sigma_hist, xi_dim, pi_cg, R_arr])

def compute_F(w5):
    """w5 = [γ, -λ, η, κ, α_R]"""
    return log_H + feat_6 @ w5

def compute_scores(F):
    w2 = int(np.sum(F[A_pairs_2d[:, 0]] < F[A_pairs_2d[:, 1]]))
    w3 = int(np.sum(F[A_pairs_3d[:, 0]] < F[A_pairs_3d[:, 1]]))
    w5 = int(np.sum(F[A_pairs_5d[:, 0]] < F[A_pairs_5d[:, 1]]))
    t_a = len(A_pairs_2d) + len(A_pairs_3d) + len(A_pairs_5d)
    A = (w2 + w3 + w5) / t_a
    A2 = w2 / len(A_pairs_2d); A3 = w3 / len(A_pairs_3d); A5 = w5 / len(A_pairs_5d)
    wb = int(np.sum(F[B_pairs[:, 0]] < F[B_pairs[:, 1]]))
    B = wb / len(B_pairs)
    # C: fraction of N-groups with corr(Σ_hist, F) < 0
    neg = 0
    for nv in unique_N:
        idx = np.where(N_arr == nv)[0]
        sh = sigma_hist[idx]; f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            if np.corrcoef(sh, f)[0, 1] < 0: neg += 1
    C = neg / len(unique_N)
    D = abs(np.corrcoef(pi_cg, F)[0, 1])
    return A, A2, A3, A5, B, C, D

# ── Phase 1: Fine 5D grid around balanced region ──
print("\n" + "=" * 80)
print("PHASE 1: Fine 5D grid around balanced region")
print("  Center: γ=0.066, λ=3.27, η=3.5, κ=0.075, α_R=30")
print("=" * 80)

gamma_vals = [0.0, 0.05, 0.10]
lam_vals = [1.0, 2.0, 3.0, 4.0, 5.0]
eta_vals = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
kappa_vals = [0.0, 0.075, 0.15]
alpha_R_vals = [18, 21, 24, 27, 30, 33, 36, 39]

total = len(gamma_vals) * len(lam_vals) * len(eta_vals) * len(kappa_vals) * len(alpha_R_vals)
print(f"\n  Scanning {total} points...")
t0 = time.time()

results = []
for gam in gamma_vals:
    for lam in lam_vals:
        for eta in eta_vals:
            for kap in kappa_vals:
                for aR in alpha_R_vals:
                    w = np.array([gam, -lam, eta, kap, aR], dtype=float)
                    F = compute_F(w)
                    A, A2, A3, A5, B, C, D = compute_scores(F)
                    results.append((gam, lam, eta, kap, aR, A, A2, A3, A5, B, C, D))

dt = time.time() - t0
print(f"  Done: {len(results)} points in {dt:.1f}s")

# ── Classify results ──
# A-max: top 10 by A
# B-preserve: top 10 by A with B >= 0.95
# Balanced: top 10 by A + B

print("\n" + "=" * 80)
print("A-MAX: Top 10 by A (unconstrained)")
print("=" * 80)
results.sort(key=lambda x: -x[5])
print(f"  {'γ':>5} {'λ':>5} {'η':>5} {'κ':>5} {'α_R':>5} | {'A':>6} {'v2D':>5} {'v3D':>5} {'v5D':>5} | {'B':>6} {'C':>5} {'D':>5}")
print("  " + "-" * 75)
for r in results[:10]:
    print(f"  {r[0]:5.3f} {r[1]:5.2f} {r[2]:5.2f} {r[3]:5.3f} {r[4]:5.0f} | "
          f"{100*r[5]:5.1f}% {100*r[6]:4.0f}% {100*r[7]:4.0f}% {100*r[8]:4.0f}% | "
          f"{100*r[9]:5.1f}% {r[10]:5.2f} {r[11]:5.3f}")

print("\n" + "=" * 80)
print("B-PRESERVE: Top 10 by A with B ≥ 95%")
print("=" * 80)
b95 = [r for r in results if r[9] >= 0.95]
b95.sort(key=lambda x: -x[5])
print(f"  {'γ':>5} {'λ':>5} {'η':>5} {'κ':>5} {'α_R':>5} | {'A':>6} {'v2D':>5} {'v3D':>5} {'v5D':>5} | {'B':>6} {'C':>5} {'D':>5}")
print("  " + "-" * 75)
for r in b95[:10]:
    print(f"  {r[0]:5.3f} {r[1]:5.2f} {r[2]:5.2f} {r[3]:5.3f} {r[4]:5.0f} | "
          f"{100*r[5]:5.1f}% {100*r[6]:4.0f}% {100*r[7]:4.0f}% {100*r[8]:4.0f}% | "
          f"{100*r[9]:5.1f}% {r[10]:5.2f} {r[11]:5.3f}")

print("\n" + "=" * 80)
print("BALANCED: Top 10 by A + B")
print("=" * 80)
results.sort(key=lambda x: -(x[5] + x[9]))
print(f"  {'γ':>5} {'λ':>5} {'η':>5} {'κ':>5} {'α_R':>5} | {'A':>6} {'v2D':>5} {'v3D':>5} {'v5D':>5} | {'B':>6} {'C':>5} {'D':>5} | {'A+B':>5}")
print("  " + "-" * 80)
for r in results[:10]:
    print(f"  {r[0]:5.3f} {r[1]:5.2f} {r[2]:5.2f} {r[3]:5.3f} {r[4]:5.0f} | "
          f"{100*r[5]:5.1f}% {100*r[6]:4.0f}% {100*r[7]:4.0f}% {100*r[8]:4.0f}% | "
          f"{100*r[9]:5.1f}% {r[10]:5.2f} {r[11]:5.3f} | {100*(r[5]+r[9]):5.1f}")

print("\n" + "=" * 80)
print("BALANCED with C≥0.75: Top 10 by A + B with C ≥ 0.75")
print("=" * 80)
bal_c = [r for r in results if r[10] >= 0.75]
bal_c.sort(key=lambda x: -(x[5] + x[9]))
print(f"  {'γ':>5} {'λ':>5} {'η':>5} {'κ':>5} {'α_R':>5} | {'A':>6} {'v2D':>5} {'v3D':>5} {'v5D':>5} | {'B':>6} {'C':>5} {'D':>5} | {'A+B':>5}")
print("  " + "-" * 80)
for r in bal_c[:10]:
    print(f"  {r[0]:5.3f} {r[1]:5.2f} {r[2]:5.2f} {r[3]:5.3f} {r[4]:5.0f} | "
          f"{100*r[5]:5.1f}% {100*r[6]:4.0f}% {100*r[7]:4.0f}% {100*r[8]:4.0f}% | "
          f"{100*r[9]:5.1f}% {r[10]:5.2f} {r[11]:5.3f} | {100*(r[5]+r[9]):5.1f}")

# ── Phase 2: Multi-objective MH-MCMC ──
print("\n" + "=" * 80)
print("PHASE 2: Multi-objective MH-MCMC")
print("  Objective: maximize A_score + B_score + C_penalty")
print("  Using Bradley-Terry likelihood with equal A/B weighting")
print("=" * 80)

def log_lik_multi(w5, wt_A=1.0, wt_B=1.0, scale=2.0):
    F = compute_F(w5)
    ll = 0.0
    # A pairs
    diff_a = (F[A_all[:, 1]] - F[A_all[:, 0]]) / scale
    diff_a = np.clip(diff_a, -30, 30)
    ll += wt_A * np.sum(-np.log1p(np.exp(-diff_a)))
    # B pairs
    diff_b = (F[B_pairs[:, 1]] - F[B_pairs[:, 0]]) / scale
    diff_b = np.clip(diff_b, -30, 30)
    ll += wt_B * np.sum(-np.log1p(np.exp(-diff_b)))
    # C: penalize positive correlations
    for nv in unique_N:
        idx = np.where(N_arr == nv)[0]
        sh = sigma_hist[idx]; f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            r = np.corrcoef(sh, f)[0, 1]
            ll += 2.0 * (-r)  # reward negative correlation
    return ll

def log_pri_multi(w5):
    gam, neg_lam, eta, kap, aR = w5
    lam = -neg_lam
    if gam < 0 or lam < 0 or eta < 0 or kap < 0 or aR < 0: return -np.inf
    if gam > 5 or lam > 10 or eta > 8 or kap > 2 or aR > 80: return -np.inf
    lp = -0.5*((gam-0.05)/0.3)**2 - 0.5*((lam-3.0)/2.0)**2
    lp += -0.5*((eta-3.0)/2.0)**2 - 0.5*((kap-0.05)/0.2)**2
    lp += -0.5*((aR-25.0)/15.0)**2  # center prior at 25, broad
    return lp

def log_post_multi(w5):
    lp = log_pri_multi(w5)
    if not np.isfinite(lp): return -np.inf
    ll = log_lik_multi(w5)
    return lp + ll if np.isfinite(ll) else -np.inf

# Run MCMC
n_steps, burn_in = 100000, 20000
rng = np.random.default_rng(99)
w0 = np.array([0.05, -3.0, 3.5, 0.05, 25.0])
prop = np.array([0.03, 0.08, 0.15, 0.01, 2.0])
chain = np.zeros((n_steps, 5))
w_cur = w0.copy(); lp_cur = log_post_multi(w_cur); n_acc = 0

print(f"\nRunning MCMC ({n_steps} steps, equal A/B weighting)...")
t0 = time.time()
for step in range(n_steps):
    w_p = w_cur + rng.normal(0, prop)
    lp_p = log_post_multi(w_p)
    if np.log(rng.random()) < (lp_p - lp_cur):
        w_cur = w_p; lp_cur = lp_p; n_acc += 1
    chain[step] = w_cur
    if step < burn_in and step > 0 and step % 2000 == 0:
        rate = n_acc / (step + 1)
        if rate < 0.15: prop *= 0.8
        elif rate > 0.40: prop *= 1.2
    if step % 20000 == 0 and step > 0:
        print(f"  step {step}, acc={n_acc/(step+1):.3f}")

dt = time.time() - t0
samples = chain[burn_in::5]
print(f"  Done in {dt:.1f}s. Acceptance: {n_acc/n_steps:.3f}, samples: {len(samples)}")

disp = samples.copy()
disp[:, 1] = -disp[:, 1]
names = ["γ", "λ", "η", "κ", "α_R"]
print(f"\n  {'Param':>6s} {'Mean':>8s} {'Median':>8s} {'Std':>8s} {'2.5%':>8s} {'97.5%':>8s}")
print(f"  {'-'*6} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
for i, name in enumerate(names):
    s = disp[:, i]
    print(f"  {name:>6s} {np.mean(s):8.4f} {np.median(s):8.4f} "
          f"{np.std(s):8.4f} {np.percentile(s, 2.5):8.4f} {np.percentile(s, 97.5):8.4f}")

# Evaluate
w_med = np.median(samples, axis=0)
F_med = compute_F(w_med)
A, A2, A3, A5, B, C, D = compute_scores(F_med)
print(f"\n  MCMC median scores: A={100*A:.1f}% (v2D={100*A2:.1f}%, v3D={100*A3:.1f}%, v5D={100*A5:.1f}%), "
      f"B={100*B:.1f}%, C={C:.2f}, D={D:.3f}")

# Also evaluate mean
w_mean = np.mean(samples, axis=0)
F_mean = compute_F(w_mean)
Am, A2m, A3m, A5m, Bm, Cm, Dm = compute_scores(F_mean)
print(f"  MCMC mean scores:   A={100*Am:.1f}% (v2D={100*A2m:.1f}%, v3D={100*A3m:.1f}%, v5D={100*A5m:.1f}%), "
      f"B={100*Bm:.1f}%, C={Cm:.2f}, D={Dm:.3f}")

# ── Phase 3: B-heavy MCMC ──
print("\n" + "=" * 80)
print("PHASE 3: B-heavy MH-MCMC (wt_B=3.0 vs wt_A=1.0)")
print("=" * 80)

def log_post_b_heavy(w5):
    lp = log_pri_multi(w5)
    if not np.isfinite(lp): return -np.inf
    ll = log_lik_multi(w5, wt_A=1.0, wt_B=3.0)
    return lp + ll if np.isfinite(ll) else -np.inf

chain2 = np.zeros((n_steps, 5))
w_cur = w0.copy(); lp_cur = log_post_b_heavy(w_cur); n_acc = 0
prop2 = np.array([0.03, 0.08, 0.15, 0.01, 2.0])

print(f"\nRunning MCMC ({n_steps} steps, B-heavy)...")
t0 = time.time()
for step in range(n_steps):
    w_p = w_cur + rng.normal(0, prop2)
    lp_p = log_post_b_heavy(w_p)
    if np.log(rng.random()) < (lp_p - lp_cur):
        w_cur = w_p; lp_cur = lp_p; n_acc += 1
    chain2[step] = w_cur
    if step < burn_in and step > 0 and step % 2000 == 0:
        rate = n_acc / (step + 1)
        if rate < 0.15: prop2 *= 0.8
        elif rate > 0.40: prop2 *= 1.2

dt = time.time() - t0
samples2 = chain2[burn_in::5]
print(f"  Done in {dt:.1f}s. Acceptance: {n_acc/n_steps:.3f}")

disp2 = samples2.copy()
disp2[:, 1] = -disp2[:, 1]
print(f"\n  {'Param':>6s} {'Mean':>8s} {'Median':>8s} {'Std':>8s}")
print(f"  {'-'*6} {'-'*8} {'-'*8} {'-'*8}")
for i, name in enumerate(names):
    s = disp2[:, i]
    print(f"  {name:>6s} {np.mean(s):8.4f} {np.median(s):8.4f} {np.std(s):8.4f}")

w_med2 = np.median(samples2, axis=0)
F_med2 = compute_F(w_med2)
Ab, A2b, A3b, A5b, Bb, Cb, Db = compute_scores(F_med2)
print(f"\n  B-heavy MCMC median: A={100*Ab:.1f}% (v2D={100*A2b:.1f}%, v3D={100*A3b:.1f}%, v5D={100*A5b:.1f}%), "
      f"B={100*Bb:.1f}%, C={Cb:.2f}, D={Db:.3f}")

# ── FINAL COMPARISON ──
print("\n" + "=" * 80)
print("FINAL COMPARISON TABLE")
print("=" * 80)

def print_model(label, F):
    A, A2, A3, A5, B, C, D = compute_scores(F)
    print(f"  {label:>30} | {100*A:5.1f}% {100*A2:5.1f}% {100*A3:5.1f}% {100*A5:5.1f}% | {100*B:5.1f}% | {C:4.2f} {D:5.3f}")

print(f"  {'Model':>30} | {'A':>6} {'v2D':>6} {'v3D':>6} {'v5D':>6} | {'B':>6} | {'C':>4} {'D':>5}")
print("  " + "-" * 80)

# F5 baselines
F5_bay = log_H + 0.0004*pi_geo - 0.888*sigma_hist + 0.637*xi_dim + 0.068*pi_cg
print_model("F5 (Bayesian)", F5_bay)

F5_cal = log_H + 0.5*pi_geo - 1.5*sigma_hist + 0.1*xi_dim + 0.05*pi_cg
print_model("F5 (Calibrated)", F5_cal)

# Grid scan best balanced with C
if bal_c:
    r = bal_c[0]
    F_grid = log_H + r[0]*pi_geo - r[1]*sigma_hist + r[2]*xi_dim + r[3]*pi_cg + r[4]*R_arr
    print_model(f"Grid balanced (aR={r[4]:.0f})", F_grid)

# MCMC equal weight
print_model("MCMC equal A/B", F_med)

# MCMC B-heavy
print_model("MCMC B-heavy (3:1)", F_med2)

# N-breakdown for best model (MCMC equal)
print(f"\n  N-breakdown (MCMC equal A/B):")
for nv in unique_N:
    for fo in ["Lor2D", "Lor3D", "Lor5D"]:
        i4 = np.where((N_arr == nv) & (families_arr == "Lor4D"))[0]
        io = np.where((N_arr == nv) & (families_arr == fo))[0]
        wins = int(np.sum(F_med[i4, None] < F_med[None, io]))
        total = len(i4) * len(io)
        print(f"    N={int(nv):2d}, 4D vs {fo:>5}: {wins:2d}/{total:2d} ({100*wins/total:5.1f}%)")
    # B at this N
    i2 = np.where((N_arr == nv) & (families_arr == "Lor2D"))[0]
    ik = np.where((N_arr == nv) & (families_arr == "KR_like"))[0]
    wb = int(np.sum(F_med[i2, None] < F_med[None, ik]))
    tb = len(i2) * len(ik)
    print(f"    N={int(nv):2d}, 2D vs KR   : {wb:2d}/{tb:2d} ({100*wb/tb:5.1f}%)")
