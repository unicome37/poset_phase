"""Pareto grid scan: (α_R, η) 2D grid with A/B/C/D scores.

Fixed: β=1, γ=0.066, λ=3.27, κ=0.075 (MCMC medians from F6 run).
Scan: α_R ∈ [0, 60], η ∈ [0, 5].
Scores:
  A = 4D < 2D/3D/5D pairwise hit rate
  B = Lor2D < KR_like pairwise hit rate
  C = fraction of N-groups where corr(Σ_hist, F) < 0
  D = |corr(Π_cg, F)| (should remain nonzero)
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
print("Computing R for all posets...")
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
B_pairs = build_pairs("Lor2D", "KR_like")

# ── Decomposition: why B breaks ──
print("\n" + "=" * 80)
print("DECOMPOSITION: ΔF contributions at MCMC median weights")
print("  γ=0.066, λ=3.27, η=2.49, κ=0.075, α_R=47.5")
print("=" * 80)

# Compute individual contributions per family
gamma_fix, lam_fix, eta_fix, kappa_fix = 0.066, 3.27, 2.49, 0.075

print(f"\n  {'Fam':>7} {'N':>3} | {'logH':>7} {'γΠg':>7} {'-λΣh':>7} {'ηΞd':>7} {'κΠcg':>7} | {'F5':>7} | {'αR·R':>7} | {'F6':>7}")
print("  " + "-" * 85)
for fam in ["Lor2D", "Lor4D", "KR_like", "Lor3D", "Lor5D"]:
    for n in [16, 20, 28, 36]:
        idx = np.where((N_arr == n) & (families_arr == fam))[0]
        lh = np.mean(log_H[idx])
        gpg = gamma_fix * np.mean(pi_geo[idx])
        lsh = -lam_fix * np.mean(sigma_hist[idx])
        exd = eta_fix * np.mean(xi_dim[idx])
        kpc = kappa_fix * np.mean(pi_cg[idx])
        f5 = lh + gpg + lsh + exd + kpc
        ar = 47.5 * np.mean(R_arr[idx])
        f6 = f5 + ar
        print(f"  {fam:>7} {n:3d} | {lh:7.2f} {gpg:7.3f} {lsh:7.2f} {exd:7.2f} {kpc:7.4f} | {f5:7.2f} | {ar:7.2f} | {f6:7.2f}")
    print()

# Gap analysis
print("\nGAP ANALYSIS: F6(Lor2D) - F6(KR_like) per N")
print("  (Positive means Lor2D > KR → B fails; Negative means B passes)")
print(f"  {'N':>3} | {'ΔlogH':>7} {'ΔγΠg':>7} {'Δ-λΣh':>7} {'ΔηΞd':>7} {'ΔκΠcg':>7} | {'ΔF5':>7} | {'ΔαR·R':>7} | {'ΔF6':>7}")
print("  " + "-" * 80)
for n in [16, 20, 28, 36]:
    i2 = np.where((N_arr == n) & (families_arr == "Lor2D"))[0]
    ik = np.where((N_arr == n) & (families_arr == "KR_like"))[0]
    d_lh = np.mean(log_H[i2]) - np.mean(log_H[ik])
    d_gpg = gamma_fix * (np.mean(pi_geo[i2]) - np.mean(pi_geo[ik]))
    d_lsh = -lam_fix * (np.mean(sigma_hist[i2]) - np.mean(sigma_hist[ik]))
    d_exd = eta_fix * (np.mean(xi_dim[i2]) - np.mean(xi_dim[ik]))
    d_kpc = kappa_fix * (np.mean(pi_cg[i2]) - np.mean(pi_cg[ik]))
    d_f5 = d_lh + d_gpg + d_lsh + d_exd + d_kpc
    d_ar = 47.5 * (np.mean(R_arr[i2]) - np.mean(R_arr[ik]))
    d_f6 = d_f5 + d_ar
    print(f"  {n:3d} | {d_lh:+7.2f} {d_gpg:+7.3f} {d_lsh:+7.2f} {d_exd:+7.2f} {d_kpc:+7.4f} | {d_f5:+7.2f} | {d_ar:+7.2f} | {d_f6:+7.2f}")

print("\nGAP ANALYSIS: F6(Lor4D) - F6(Lor2D) per N")
print("  (Negative means 4D < 2D → A passes)")
for n in [16, 20, 28, 36]:
    i4 = np.where((N_arr == n) & (families_arr == "Lor4D"))[0]
    i2 = np.where((N_arr == n) & (families_arr == "Lor2D"))[0]
    d_f5 = np.mean(log_H[i4] + gamma_fix*pi_geo[i4] - lam_fix*sigma_hist[i4] + eta_fix*xi_dim[i4] + kappa_fix*pi_cg[i4]) - \
           np.mean(log_H[i2] + gamma_fix*pi_geo[i2] - lam_fix*sigma_hist[i2] + eta_fix*xi_dim[i2] + kappa_fix*pi_cg[i2])
    d_ar = 47.5 * (np.mean(R_arr[i4]) - np.mean(R_arr[i2]))
    print(f"  N={n:2d}: ΔF5={d_f5:+7.2f}, Δα_R·R={d_ar:+7.2f}, ΔF6={d_f5+d_ar:+7.2f}")

# ── 2D Grid scan ──
print("\n" + "=" * 80)
print("2D GRID SCAN: (α_R, η), fixed γ=0.066, λ=3.27, κ=0.075")
print("=" * 80)

def compute_all_scores(F):
    # A: 4D < others
    w2 = int(np.sum(F[A_pairs_2d[:, 0]] < F[A_pairs_2d[:, 1]])) if len(A_pairs_2d) else 0
    w3 = int(np.sum(F[A_pairs_3d[:, 0]] < F[A_pairs_3d[:, 1]])) if len(A_pairs_3d) else 0
    w5 = int(np.sum(F[A_pairs_5d[:, 0]] < F[A_pairs_5d[:, 1]])) if len(A_pairs_5d) else 0
    t_a = len(A_pairs_2d) + len(A_pairs_3d) + len(A_pairs_5d)
    A_score = (w2 + w3 + w5) / t_a if t_a else 0
    A_2d = w2 / len(A_pairs_2d) if len(A_pairs_2d) else 0
    A_3d = w3 / len(A_pairs_3d) if len(A_pairs_3d) else 0
    A_5d = w5 / len(A_pairs_5d) if len(A_pairs_5d) else 0
    # B: Lor2D < KR
    wb = int(np.sum(F[B_pairs[:, 0]] < F[B_pairs[:, 1]])) if len(B_pairs) else 0
    B_score = wb / len(B_pairs) if len(B_pairs) else 0
    # C: fraction of N-groups with corr(Σ_hist, F) < 0
    neg_count = 0
    for nv in unique_N:
        idx = np.where(N_arr == nv)[0]
        sh = sigma_hist[idx]; f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            if np.corrcoef(sh, f)[0, 1] < 0:
                neg_count += 1
    C_score = neg_count / len(unique_N)
    # D: |corr(Π_cg, F)| overall
    D_score = abs(np.corrcoef(pi_cg, F)[0, 1]) if np.std(pi_cg) > 1e-10 else 0
    return A_score, A_2d, A_3d, A_5d, B_score, C_score, D_score

# Grid
alpha_R_grid = np.arange(0, 61, 3)  # 0..60 step 3
eta_grid = np.arange(0, 5.1, 0.25)  # 0..5 step 0.25

results = []
print(f"\n  Scanning {len(alpha_R_grid)} x {len(eta_grid)} = {len(alpha_R_grid)*len(eta_grid)} points...")
t0 = time.time()

for aR in alpha_R_grid:
    for eta in eta_grid:
        F = log_H + gamma_fix * pi_geo - lam_fix * sigma_hist + eta * xi_dim + kappa_fix * pi_cg + aR * R_arr
        A, A2, A3, A5, B, C, D = compute_all_scores(F)
        results.append((aR, eta, A, A2, A3, A5, B, C, D))

dt = time.time() - t0
print(f"  Done in {dt:.1f}s")

# Filter feasible island: A ≥ 0.65, B ≥ 0.95, C ≥ 0.75 (3/4 N-groups)
print("\n" + "=" * 80)
print("FEASIBLE ISLAND: A ≥ 0.65, B ≥ 0.95, C ≥ 0.75")
print("=" * 80)
feasible = [(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]) 
            for r in results if r[2] >= 0.65 and r[6] >= 0.95 and r[7] >= 0.75]

if feasible:
    print(f"\n  Found {len(feasible)} feasible points!")
    print(f"  {'α_R':>5} {'η':>5} | {'A':>6} {'vs2D':>6} {'vs3D':>6} {'vs5D':>6} | {'B':>6} {'C':>5} {'D':>5}")
    print("  " + "-" * 65)
    for r in sorted(feasible, key=lambda x: -x[2]):
        print(f"  {r[0]:5.0f} {r[1]:5.2f} | {100*r[2]:5.1f}% {100*r[3]:5.1f}% {100*r[4]:5.1f}% {100*r[5]:5.1f}% | "
              f"{100*r[6]:5.1f}% {r[7]:5.2f} {r[8]:5.3f}")
else:
    print("\n  No points satisfy A≥0.65 AND B≥0.95!")
    # Relax to B ≥ 0.90
    print("\n  Relaxing to B ≥ 0.90:")
    feasible90 = [(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8])
                  for r in results if r[2] >= 0.65 and r[6] >= 0.90 and r[7] >= 0.75]
    if feasible90:
        print(f"  Found {len(feasible90)} points with B≥90%:")
        print(f"  {'α_R':>5} {'η':>5} | {'A':>6} {'vs2D':>6} {'vs3D':>6} {'vs5D':>6} | {'B':>6} {'C':>5} {'D':>5}")
        print("  " + "-" * 65)
        for r in sorted(feasible90, key=lambda x: -x[2])[:30]:
            print(f"  {r[0]:5.0f} {r[1]:5.2f} | {100*r[2]:5.1f}% {100*r[3]:5.1f}% {100*r[4]:5.1f}% {100*r[5]:5.1f}% | "
                  f"{100*r[6]:5.1f}% {r[7]:5.2f} {r[8]:5.3f}")
    else:
        print("  No points even with B≥90%!")

# Best A with B >= 0.95
print("\n" + "-" * 60)
b95 = [r for r in results if r[6] >= 0.95]
if b95:
    best = max(b95, key=lambda x: x[2])
    print(f"  Best A with B≥95%: α_R={best[0]:.0f}, η={best[1]:.2f}, A={100*best[2]:.1f}%, B={100*best[6]:.1f}%, C={best[7]:.2f}, D={best[8]:.3f}")

# Best A with B >= 0.90
b90 = [r for r in results if r[6] >= 0.90]
if b90:
    best = max(b90, key=lambda x: x[2])
    print(f"  Best A with B≥90%: α_R={best[0]:.0f}, η={best[1]:.2f}, A={100*best[2]:.1f}%, B={100*best[6]:.1f}%, C={best[7]:.2f}, D={best[8]:.3f}")

# Best balanced: A + B
best_sum = max(results, key=lambda x: x[2] + x[6])
print(f"  Best A+B sum: α_R={best_sum[0]:.0f}, η={best_sum[1]:.2f}, A={100*best_sum[2]:.1f}%, B={100*best_sum[6]:.1f}%, C={best_sum[7]:.2f}, D={best_sum[8]:.3f}")

# Best A (unconstrained)
best_a = max(results, key=lambda x: x[2])
print(f"  Best A (unconstrained): α_R={best_a[0]:.0f}, η={best_a[1]:.2f}, A={100*best_a[2]:.1f}%, B={100*best_a[6]:.1f}%, C={best_a[7]:.2f}, D={best_a[8]:.3f}")

# ── Pareto front: non-dominated points on (A, B) ──
print("\n" + "=" * 80)
print("PARETO FRONT: A vs B (non-dominated)")
print("=" * 80)

# For each unique (A, B) pair, keep the one with best C
results_arr = np.array([(r[2], r[6]) for r in results])
pareto = []
for i, (a, b) in enumerate(results_arr):
    dominated = False
    for j, (a2, b2) in enumerate(results_arr):
        if i != j and a2 >= a and b2 >= b and (a2 > a or b2 > b):
            dominated = True
            break
    if not dominated:
        pareto.append(results[i])

pareto.sort(key=lambda x: -x[2])
print(f"\n  {len(pareto)} Pareto-optimal points:")
print(f"  {'α_R':>5} {'η':>5} | {'A':>6} {'vs2D':>6} {'vs3D':>6} {'vs5D':>6} | {'B':>6} {'C':>5} {'D':>5}")
print("  " + "-" * 65)
for r in pareto:
    tag = " ★" if r[2] >= 0.65 and r[6] >= 0.95 else (" ●" if r[2] >= 0.65 and r[6] >= 0.90 else "")
    print(f"  {r[0]:5.0f} {r[1]:5.2f} | {100*r[2]:5.1f}% {100*r[3]:5.1f}% {100*r[4]:5.1f}% {100*r[5]:5.1f}% | "
          f"{100*r[6]:5.1f}% {r[7]:5.2f} {r[8]:5.3f}{tag}")

# ── Also scan λ dimension (maybe λ can help B) ──
print("\n" + "=" * 80)
print("3D TARGETED SCAN: (α_R, η, λ) around best balanced region")
print("=" * 80)

# Use the balanced point as center, scan λ too
alpha_R_fine = np.arange(15, 36, 3)
eta_fine = np.arange(1.5, 4.1, 0.5)
lam_fine = np.arange(1.0, 5.1, 0.5)

best_3d = None; best_3d_score = -1
results_3d = []

for aR in alpha_R_fine:
    for eta in eta_fine:
        for lam in lam_fine:
            F = log_H + gamma_fix * pi_geo - lam * sigma_hist + eta * xi_dim + kappa_fix * pi_cg + aR * R_arr
            A, A2, A3, A5, B, C, D = compute_all_scores(F)
            results_3d.append((aR, eta, lam, A, A2, A3, A5, B, C, D))
            score = A + B  # simple balanced metric
            if B >= 0.95 and A > best_3d_score:
                best_3d_score = A; best_3d = (aR, eta, lam, A, A2, A3, A5, B, C, D)

print(f"\n  Scanned {len(results_3d)} points")
if best_3d:
    r = best_3d
    print(f"\n  Best A with B≥95%:")
    print(f"    α_R={r[0]:.0f}, η={r[1]:.2f}, λ={r[2]:.2f}")
    print(f"    A={100*r[3]:.1f}% (vs2D={100*r[4]:.1f}%, vs3D={100*r[5]:.1f}%, vs5D={100*r[6]:.1f}%)")
    print(f"    B={100*r[7]:.1f}%, C={r[8]:.2f}, D={r[9]:.3f}")

# Top 10 by A with B ≥ 0.95
print("\n  Top 10 (A with B≥95%):")
feas = [r for r in results_3d if r[7] >= 0.95]
feas.sort(key=lambda x: -x[3])
print(f"  {'α_R':>5} {'η':>5} {'λ':>5} | {'A':>6} {'vs2D':>6} {'vs3D':>6} {'vs5D':>6} | {'B':>6} {'C':>5} {'D':>5}")
print("  " + "-" * 70)
for r in feas[:10]:
    print(f"  {r[0]:5.0f} {r[1]:5.2f} {r[2]:5.2f} | {100*r[3]:5.1f}% {100*r[4]:5.1f}% {100*r[5]:5.1f}% {100*r[6]:5.1f}% | "
          f"{100*r[7]:5.1f}% {r[8]:5.2f} {r[9]:5.3f}")

# Top 10 by A with B ≥ 0.90
print("\n  Top 10 (A with B≥90%):")
feas90 = [r for r in results_3d if r[7] >= 0.90]
feas90.sort(key=lambda x: -x[3])
for r in feas90[:10]:
    print(f"  {r[0]:5.0f} {r[1]:5.2f} {r[2]:5.2f} | {100*r[3]:5.1f}% {100*r[4]:5.1f}% {100*r[5]:5.1f}% {100*r[6]:5.1f}% | "
          f"{100*r[7]:5.1f}% {r[8]:5.2f} {r[9]:5.3f}")
