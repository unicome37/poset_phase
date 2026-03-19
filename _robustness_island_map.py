"""Robustness island maps for sigmoid barrier F7.

Produces two ASCII heatmaps:
  1. (Rc, w) plane at fixed α=26, η=1.0, λ=5.0
  2. (α, η) plane at fixed Rc=0.25, w=0.015, λ=5.0

Each cell shows a symbol: ★=A≥65%+B≥95%+C≥75%, ●=A≥65%+B≥95%, ◆=A≥65%+B≥90%, ·=else
Plus prints exact A/B/C for the main model point and its neighborhood.
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
B_pairs = build_pairs("Lor2D", "KR_like")

gam0, kap0 = 0.0004, 0.068

def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))

def compute_scores(F):
    w2 = int(np.sum(F[A_pairs_2d[:, 0]] < F[A_pairs_2d[:, 1]]))
    w3 = int(np.sum(F[A_pairs_3d[:, 0]] < F[A_pairs_3d[:, 1]]))
    w5 = int(np.sum(F[A_pairs_5d[:, 0]] < F[A_pairs_5d[:, 1]]))
    t_a = len(A_pairs_2d) + len(A_pairs_3d) + len(A_pairs_5d)
    A = (w2 + w3 + w5) / t_a
    A2 = w2 / len(A_pairs_2d); A3 = w3 / len(A_pairs_3d); A5 = w5 / len(A_pairs_5d)
    wb = int(np.sum(F[B_pairs[:, 0]] < F[B_pairs[:, 1]]))
    B = wb / len(B_pairs)
    neg = 0
    for nv in unique_N:
        idx = np.where(N_arr == nv)[0]
        sh = sigma_hist[idx]; f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            if np.corrcoef(sh, f)[0, 1] < 0: neg += 1
    C = neg / len(unique_N)
    D = abs(np.corrcoef(pi_cg, F)[0, 1])
    return A, A2, A3, A5, B, C, D

def grade(A, B, C):
    if A >= 0.65 and B >= 0.95 and C >= 0.75: return "★"
    if A >= 0.65 and B >= 0.95: return "●"
    if A >= 0.65 and B >= 0.90: return "◆"
    if A >= 0.55 and B >= 0.90: return "○"
    return "·"

# ══════════════════════════════════════════════════════════════════════════
# MAP 1: (Rc, w) plane, fixed α=26, η=1.0, λ=5.0
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("MAP 1: (Rc, w) robustness island — fixed α=26, η=1.0, λ=5.0")
print("  ★=A≥65%+B≥95%+C≥75%  ●=A≥65%+B≥95%  ◆=A≥65%+B≥90%  ○=A≥55%+B≥90%  ·=else")
print("=" * 90)

Rc_vals = np.arange(0.10, 0.46, 0.02)
w_vals = np.array([0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.08, 0.10, 0.15])
alpha_fix, eta_fix, lam_fix = 26, 1.0, 5.0

# Header
print(f"  {'Rc':>5}", end="")
for w in w_vals:
    print(f" {w:.3f}", end="")
print()
print("  " + "-" * (6 + len(w_vals) * 6))

map1_data = {}
for Rc in Rc_vals:
    row = ""
    for w in w_vals:
        wall = alpha_fix * sigmoid((R_arr - Rc) / w)
        F = log_H + gam0*pi_geo - lam_fix*sigma_hist + eta_fix*xi_dim + kap0*pi_cg + wall
        sc = compute_scores(F)
        g = grade(sc[0], sc[4], sc[5])
        row += f"    {g} "
        map1_data[(Rc, w)] = sc
    print(f"  {Rc:5.2f}{row}")

# Print details for the main model neighborhood
print(f"\n  Details around main point (Rc=0.25, w=0.015):")
print(f"  {'Rc':>5} {'w':>6} | {'A':>6} {'v2D':>5} {'v3D':>5} {'v5D':>5} | {'B':>6} {'C':>4} {'D':>5}")
print("  " + "-" * 60)
for Rc in [0.22, 0.24, 0.25, 0.26, 0.28]:
    for w in [0.010, 0.015, 0.020]:
        wall = alpha_fix * sigmoid((R_arr - Rc) / w)
        F = log_H + gam0*pi_geo - lam_fix*sigma_hist + eta_fix*xi_dim + kap0*pi_cg + wall
        A, A2, A3, A5, B, C, D = compute_scores(F)
        g = grade(A, B, C)
        print(f"  {Rc:5.2f} {w:6.3f} | {100*A:5.1f}% {100*A2:4.0f}% {100*A3:4.0f}% {100*A5:4.0f}% | {100*B:5.1f}% {C:4.2f} {D:5.3f} {g}")

# ══════════════════════════════════════════════════════════════════════════
# MAP 2: (α, η) plane, fixed Rc=0.25, w=0.015, λ=5.0
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("MAP 2: (α, η) robustness island — fixed Rc=0.25, w=0.015, λ=5.0")
print("=" * 90)

alpha_vals = np.arange(5, 51, 3)
eta_vals = np.arange(0.0, 4.1, 0.25)
Rc_fix, w_fix = 0.25, 0.015

# Header
print(f"  {'α\\η':>5}", end="")
for eta in eta_vals:
    print(f" {eta:4.1f}", end="")
print()
print("  " + "-" * (6 + len(eta_vals) * 5))

for alpha in alpha_vals:
    row = ""
    for eta in eta_vals:
        wall = alpha * sigmoid((R_arr - Rc_fix) / w_fix)
        F = log_H + gam0*pi_geo - lam_fix*sigma_hist + eta*xi_dim + kap0*pi_cg + wall
        sc = compute_scores(F)
        g = grade(sc[0], sc[4], sc[5])
        row += f"    {g}"
    print(f"  {alpha:5.0f}{row}")

# Details around main point
print(f"\n  Details around main point (α=26, η=1.0):")
print(f"  {'α':>5} {'η':>5} | {'A':>6} {'v2D':>5} {'v3D':>5} {'v5D':>5} | {'B':>6} {'C':>4} {'D':>5}")
print("  " + "-" * 60)
for alpha in [20, 23, 26, 29, 32]:
    for eta in [0.5, 0.75, 1.0, 1.25, 1.5]:
        wall = alpha * sigmoid((R_arr - Rc_fix) / w_fix)
        F = log_H + gam0*pi_geo - lam_fix*sigma_hist + eta*xi_dim + kap0*pi_cg + wall
        A, A2, A3, A5, B, C, D = compute_scores(F)
        g = grade(A, B, C)
        print(f"  {alpha:5.0f} {eta:5.2f} | {100*A:5.1f}% {100*A2:4.0f}% {100*A3:4.0f}% {100*A5:4.0f}% | {100*B:5.1f}% {C:4.2f} {D:5.3f} {g}")

# ══════════════════════════════════════════════════════════════════════════
# MAP 3: (α, λ) plane, fixed Rc=0.25, w=0.015, η=1.0
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("MAP 3: (α, λ) robustness island — fixed Rc=0.25, w=0.015, η=1.0")
print("=" * 90)

lam_vals = np.arange(0.5, 8.1, 0.5)

print(f"  {'α\\λ':>5}", end="")
for lam in lam_vals:
    print(f" {lam:4.1f}", end="")
print()
print("  " + "-" * (6 + len(lam_vals) * 5))

for alpha in alpha_vals:
    row = ""
    for lam in lam_vals:
        wall = alpha * sigmoid((R_arr - Rc_fix) / w_fix)
        F = log_H + gam0*pi_geo - lam*sigma_hist + eta_fix*xi_dim + kap0*pi_cg + wall
        sc = compute_scores(F)
        g = grade(sc[0], sc[4], sc[5])
        row += f"    {g}"
    print(f"  {alpha:5.0f}{row}")

# ══════════════════════════════════════════════════════════════════════════
# ISLAND STATISTICS
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("ISLAND STATISTICS")
print("=" * 90)

# Count cells in each grade across all three maps
total = 0; n_star = 0; n_dot = 0; n_dia = 0; n_circ = 0
for Rc in Rc_vals:
    for w in w_vals:
        for eta in eta_vals:
            for lam in [3.0, 4.0, 5.0, 6.0]:
                for alpha in [20, 23, 26, 29, 32]:
                    wall = alpha * sigmoid((R_arr - Rc) / w)
                    F = log_H + gam0*pi_geo - lam*sigma_hist + eta*xi_dim + kap0*pi_cg + wall
                    sc = compute_scores(F)
                    g = grade(sc[0], sc[4], sc[5])
                    total += 1
                    if g == "★": n_star += 1
                    elif g == "●": n_dot += 1
                    elif g == "◆": n_dia += 1
                    elif g == "○": n_circ += 1

print(f"  Comprehensive scan: {total} points")
print(f"    ★ (A≥65%+B≥95%+C≥75%): {n_star} ({100*n_star/total:.1f}%)")
print(f"    ● (A≥65%+B≥95%):        {n_dot} ({100*n_dot/total:.1f}%)")
print(f"    ◆ (A≥65%+B≥90%):        {n_dia} ({100*n_dia/total:.1f}%)")
print(f"    ○ (A≥55%+B≥90%):        {n_circ} ({100*n_circ/total:.1f}%)")
print(f"    · (below threshold):     {total-n_star-n_dot-n_dia-n_circ} ({100*(total-n_star-n_dot-n_dia-n_circ)/total:.1f}%)")

# ── N-breakdown for main model ──
print("\n" + "=" * 90)
print("N-BREAKDOWN: Main model (Rc=0.25, w=0.015, α=26, η=1.0, λ=5.0)")
print("=" * 90)

wall_main = alpha_fix * sigmoid((R_arr - 0.25) / 0.015)
F_main = log_H + gam0*pi_geo - lam_fix*sigma_hist + eta_fix*xi_dim + kap0*pi_cg + wall_main

for nv in unique_N:
    print(f"\n  N = {int(nv)}:")
    for fo in ["Lor2D", "Lor3D", "Lor5D"]:
        i4 = np.where((N_arr == nv) & (families_arr == "Lor4D"))[0]
        io = np.where((N_arr == nv) & (families_arr == fo))[0]
        wins = int(np.sum(F_main[i4, None] < F_main[None, io]))
        total_p = len(i4) * len(io)
        print(f"    4D vs {fo:>5}: {wins:2d}/{total_p:2d} ({100*wins/total_p:5.1f}%)")
    i2 = np.where((N_arr == nv) & (families_arr == "Lor2D"))[0]
    ik = np.where((N_arr == nv) & (families_arr == "KR_like"))[0]
    wb = int(np.sum(F_main[i2, None] < F_main[None, ik]))
    tb = len(i2) * len(ik)
    print(f"    2D vs   KR: {wb:2d}/{tb:2d} ({100*wb/tb:5.1f}%)")
    # C at this N
    idx = np.where(N_arr == nv)[0]
    r = np.corrcoef(sigma_hist[idx], F_main[idx])[0, 1]
    print(f"    corr(Σ_hist, F) = {r:+.3f} {'✓' if r < 0 else '✗'}")
