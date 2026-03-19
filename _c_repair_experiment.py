"""C-repair experiment: fix sigmoid wall shape, tune (λ, α) to restore C≥0.75.

Phase 1: Constant λ scan — (α, λ) at fixed Rc=0.25, w=0.015, η=1.0
Phase 2: Weak N-scaling — λ(N) = λ0·(N/N0)^p, scan (λ0, p, α)

Success: A≥0.70, B≥0.95, B16≥0.95, C≥0.75
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
Rc_fix, w_fix = 0.25, 0.015

def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))

wall_base = sigmoid((R_arr - Rc_fix) / w_fix)  # precompute sigmoid activation

def compute_full_scores(F):
    """Return A, A2, A3, A5, B, C, D, per-N C and B details."""
    w2 = int(np.sum(F[A_pairs_2d[:, 0]] < F[A_pairs_2d[:, 1]]))
    w3 = int(np.sum(F[A_pairs_3d[:, 0]] < F[A_pairs_3d[:, 1]]))
    w5 = int(np.sum(F[A_pairs_5d[:, 0]] < F[A_pairs_5d[:, 1]]))
    t_a = len(A_pairs_2d) + len(A_pairs_3d) + len(A_pairs_5d)
    A = (w2 + w3 + w5) / t_a
    A2 = w2 / len(A_pairs_2d); A3 = w3 / len(A_pairs_3d); A5 = w5 / len(A_pairs_5d)
    wb = int(np.sum(F[B_pairs[:, 0]] < F[B_pairs[:, 1]]))
    B = wb / len(B_pairs)
    
    # Per-N details
    c_per_n = {}
    b_per_n = {}
    corr_per_n = {}
    for nv in unique_N:
        idx = np.where(N_arr == nv)[0]
        sh = sigma_hist[idx]; f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            r = np.corrcoef(sh, f)[0, 1]
            c_per_n[int(nv)] = 1 if r < 0 else 0
            corr_per_n[int(nv)] = r
        else:
            c_per_n[int(nv)] = 0
            corr_per_n[int(nv)] = 0.0
        # B per N
        i2 = np.where((N_arr == nv) & (families_arr == "Lor2D"))[0]
        ik = np.where((N_arr == nv) & (families_arr == "KR_like"))[0]
        if len(i2) and len(ik):
            b_per_n[int(nv)] = int(np.sum(F[i2, None] < F[None, ik])) / (len(i2) * len(ik))
    
    C = sum(c_per_n.values()) / len(unique_N)
    D = abs(np.corrcoef(pi_cg, F)[0, 1])
    return A, A2, A3, A5, B, C, D, c_per_n, b_per_n, corr_per_n

def print_header():
    print(f"  {'Params':>22} | {'A':>5} {'v2D':>4} {'v3D':>4} {'v5D':>4} | {'B':>5} {'B16':>4} | {'C':>4} {'C16':>3} {'C20':>3} {'C28':>3} {'C36':>3} | {'r16':>6} {'r20':>6}")

def print_row(label, sc):
    A, A2, A3, A5, B, C, D, cn, bn, rn = sc
    b16 = bn.get(16, 0)
    c16 = cn.get(16, 0); c20 = cn.get(20, 0); c28 = cn.get(28, 0); c36 = cn.get(36, 0)
    r16 = rn.get(16, 0); r20 = rn.get(20, 0)
    tag = ""
    if A >= 0.70 and B >= 0.95 and b16 >= 0.95 and C >= 0.75:
        tag = " ★★★"
    elif A >= 0.65 and B >= 0.95 and C >= 0.75:
        tag = " ★"
    elif A >= 0.65 and B >= 0.95:
        tag = " ●"
    print(f"  {label:>22} | {100*A:4.1f}% {100*A2:3.0f}% {100*A3:3.0f}% {100*A5:3.0f}% | {100*B:4.1f}% {100*b16:3.0f}% | {C:4.2f}  {c16}   {c20}   {c28}   {c36}  | {r16:+.3f} {r20:+.3f}{tag}")

# ══════════════════════════════════════════════════════════════════════════
# DIAGNOSTIC: Why does C fail at N=16,20?
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 100)
print("DIAGNOSTIC: Component magnitudes by N (at main model α=26, η=1.0, λ=5.0)")
print("=" * 100)

for nv in unique_N:
    idx = np.where(N_arr == nv)[0]
    print(f"\n  N = {int(nv)}:")
    print(f"    logH:       mean={np.mean(log_H[idx]):7.2f} ± {np.std(log_H[idx]):5.2f}")
    print(f"    -λ·Σ_hist:  mean={np.mean(-5.0*sigma_hist[idx]):7.2f} ± {np.std(-5.0*sigma_hist[idx]):5.2f}")
    print(f"    η·Ξ_d:      mean={np.mean(1.0*xi_dim[idx]):7.2f} ± {np.std(1.0*xi_dim[idx]):5.2f}")
    print(f"    α·σ(wall):   mean={np.mean(26*wall_base[idx]):7.2f} ± {np.std(26*wall_base[idx]):5.2f}")
    print(f"    Total F:     mean={np.mean(log_H[idx] - 5.0*sigma_hist[idx] + 1.0*xi_dim[idx] + 26*wall_base[idx]):7.2f}")
    
    # What fraction of variance does wall contribute?
    F_nowall = log_H[idx] - 5.0*sigma_hist[idx] + 1.0*xi_dim[idx]
    F_full = F_nowall + 26*wall_base[idx]
    var_nowall = np.var(F_nowall)
    var_wall = np.var(26*wall_base[idx])
    var_full = np.var(F_full)
    print(f"    Var(F_nowall)={var_nowall:.2f}, Var(wall)={var_wall:.2f}, Var(F)={var_full:.2f}")
    print(f"    Wall variance fraction: {var_wall/var_full:.2%}")
    
    # Correlation of each component with Σ_hist within this N
    sh = sigma_hist[idx]
    if np.std(sh) > 1e-10:
        r_lh = np.corrcoef(sh, log_H[idx])[0,1]
        r_wall = np.corrcoef(sh, wall_base[idx])[0,1]
        r_xi = np.corrcoef(sh, xi_dim[idx])[0,1]
        r_f = np.corrcoef(sh, F_full)[0,1]
        print(f"    corr(Σ_hist, logH)={r_lh:+.3f}, corr(Σ_hist, wall)={r_wall:+.3f}, corr(Σ_hist, Ξ_d)={r_xi:+.3f}")
        print(f"    corr(Σ_hist, F)={r_f:+.3f} {'✓' if r_f < 0 else '✗'}")

# ══════════════════════════════════════════════════════════════════════════
# PHASE 1: Constant λ — scan (α, λ) at fixed Rc=0.25, w=0.015, η=1.0
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 100)
print("PHASE 1: Constant λ scan — (α, λ) with fixed Rc=0.25, w=0.015, η=1.0")
print("  Target: A≥0.70, B≥0.95, B16≥0.95, C≥0.75")
print("=" * 100)
print_header()

results1 = []
for alpha in range(6, 20):
    for lam in np.arange(3.0, 12.1, 0.5):
        for eta in [0.637, 1.0, 1.5]:
            F = log_H + gam0*pi_geo - lam*sigma_hist + eta*xi_dim + kap0*pi_cg + alpha*wall_base
            sc = compute_full_scores(F)
            results1.append((alpha, lam, eta, sc))

# Show all ★★★ (full success)
stars3 = [r for r in results1 if r[3][0] >= 0.70 and r[3][4] >= 0.95 and r[3][8].get(16,0) >= 0.95 and r[3][5] >= 0.75]
if stars3:
    stars3.sort(key=lambda x: -x[3][0])
    print(f"\n  ★★★ FULL SUCCESS: A≥70%+B≥95%+B16≥95%+C≥75% ({len(stars3)} points)")
    for r in stars3[:15]:
        print_row(f"α={r[0]:2d} λ={r[1]:4.1f} η={r[2]:.1f}", r[3])

# Show ★ (A≥65%+B≥95%+C≥75%)
stars1 = [r for r in results1 if r[3][0] >= 0.65 and r[3][4] >= 0.95 and r[3][5] >= 0.75]
if stars1:
    stars1.sort(key=lambda x: -x[3][0])
    print(f"\n  ★ A≥65%+B≥95%+C≥75% ({len(stars1)} points)")
    for r in stars1[:15]:
        print_row(f"α={r[0]:2d} λ={r[1]:4.1f} η={r[2]:.1f}", r[3])

# Best balanced A+B+C
best_bal = sorted(results1, key=lambda x: -(x[3][0] + x[3][4] + 0.5*x[3][5]))
print(f"\n  Top 10 balanced (A+B+0.5C):")
print_header()
for r in best_bal[:10]:
    print_row(f"α={r[0]:2d} λ={r[1]:4.1f} η={r[2]:.1f}", r[3])

# ══════════════════════════════════════════════════════════════════════════
# PHASE 2: Weak N-scaling λ(N) = λ0·(N/20)^p
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 100)
print("PHASE 2: N-scaling λ(N) = λ0·(N/20)^p — scan (λ0, p, α)")
print("  Idea: if p>0, small N gets weaker λ (less Σ_hist damping)")
print("        if p<0, small N gets stronger λ (more Σ_hist signal)")
print("=" * 100)

N0 = 20.0
results2 = []
for alpha in range(8, 18):
    for lam0 in np.arange(3.0, 12.1, 1.0):
        for p in [-0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5, 0.8, 1.0]:
            for eta in [0.637, 1.0, 1.5]:
                lam_n = lam0 * (N_arr / N0) ** p
                F = log_H + gam0*pi_geo - lam_n*sigma_hist + eta*xi_dim + kap0*pi_cg + alpha*wall_base
                sc = compute_full_scores(F)
                results2.append((alpha, lam0, p, eta, sc))

# Full success
stars3_2 = [r for r in results2 if r[4][0] >= 0.70 and r[4][4] >= 0.95 and r[4][8].get(16,0) >= 0.95 and r[4][5] >= 0.75]
if stars3_2:
    stars3_2.sort(key=lambda x: -x[4][0])
    print(f"\n  ★★★ FULL SUCCESS ({len(stars3_2)} points)")
    print_header()
    for r in stars3_2[:20]:
        print_row(f"α={r[0]:2d} λ0={r[1]:4.1f} p={r[2]:+.1f} η={r[3]:.1f}", r[4])
else:
    print("\n  No ★★★ full success with N-scaling")

# A≥65%+B≥95%+C≥75%
stars1_2 = [r for r in results2 if r[4][0] >= 0.65 and r[4][4] >= 0.95 and r[4][5] >= 0.75]
if stars1_2:
    stars1_2.sort(key=lambda x: -x[4][0])
    print(f"\n  ★ A≥65%+B≥95%+C≥75% ({len(stars1_2)} points)")
    print_header()
    for r in stars1_2[:20]:
        print_row(f"α={r[0]:2d} λ0={r[1]:4.1f} p={r[2]:+.1f} η={r[3]:.1f}", r[4])
else:
    print("\n  No ★ island with N-scaling either")

# Best balanced
best_bal2 = sorted(results2, key=lambda x: -(x[4][0] + x[4][4] + 0.5*x[4][5]))
print(f"\n  Top 10 balanced:")
print_header()
for r in best_bal2[:10]:
    print_row(f"α={r[0]:2d} λ0={r[1]:4.1f} p={r[2]:+.1f} η={r[3]:.1f}", r[4])

# ══════════════════════════════════════════════════════════════════════════
# PHASE 3: N-scaling wall α(N) = α0·(N0/N)^q  (weaker wall at small N)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 100)
print("PHASE 3: N-scaling wall α(N) = α0·(20/N)^q — scan (α0, q, λ)")
print("  Idea: q>0 → weaker wall at large N, q<0 → weaker wall at small N")
print("=" * 100)

results3 = []
for alpha0 in [10, 12, 14, 16, 18, 20, 24]:
    for q in [-0.5, -0.3, 0.0, 0.3, 0.5, 0.8, 1.0]:
        for lam in np.arange(4.0, 10.1, 1.0):
            for eta in [0.637, 1.0, 1.5]:
                alpha_n = alpha0 * (N0 / N_arr) ** q
                F = log_H + gam0*pi_geo - lam*sigma_hist + eta*xi_dim + kap0*pi_cg + alpha_n*wall_base
                sc = compute_full_scores(F)
                results3.append((alpha0, q, lam, eta, sc))

stars3_3 = [r for r in results3 if r[4][0] >= 0.70 and r[4][4] >= 0.95 and r[4][8].get(16,0) >= 0.95 and r[4][5] >= 0.75]
if stars3_3:
    stars3_3.sort(key=lambda x: -x[4][0])
    print(f"\n  ★★★ FULL SUCCESS ({len(stars3_3)} points)")
    print_header()
    for r in stars3_3[:20]:
        print_row(f"α0={r[0]:2d} q={r[1]:+.1f} λ={r[2]:4.1f} η={r[3]:.1f}", r[4])
else:
    print("\n  No ★★★ full success with N-scaling wall")

stars1_3 = [r for r in results3 if r[4][0] >= 0.65 and r[4][4] >= 0.95 and r[4][5] >= 0.75]
if stars1_3:
    stars1_3.sort(key=lambda x: -x[4][0])
    print(f"\n  ★ A≥65%+B≥95%+C≥75% ({len(stars1_3)} points)")
    print_header()
    for r in stars1_3[:20]:
        print_row(f"α0={r[0]:2d} q={r[1]:+.1f} λ={r[2]:4.1f} η={r[3]:.1f}", r[4])

# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 100)
print("FINAL SUMMARY")
print("=" * 100)

all_stars = []
for r in results1:
    if r[3][0] >= 0.65 and r[3][4] >= 0.95 and r[3][5] >= 0.75:
        all_stars.append(("const-λ", f"α={r[0]} λ={r[1]:.1f} η={r[2]:.1f}", r[3]))
for r in results2:
    if r[4][0] >= 0.65 and r[4][4] >= 0.95 and r[4][5] >= 0.75:
        all_stars.append(("λ-scale", f"α={r[0]} λ0={r[1]:.1f} p={r[2]:+.1f} η={r[3]:.1f}", r[4]))
for r in results3:
    if r[4][0] >= 0.65 and r[4][4] >= 0.95 and r[4][5] >= 0.75:
        all_stars.append(("α-scale", f"α0={r[0]} q={r[1]:+.1f} λ={r[2]:.1f} η={r[3]:.1f}", r[4]))

print(f"\n  Total ★ island points (A≥65%+B≥95%+C≥75%): {len(all_stars)}")
if all_stars:
    # Sort by A
    all_stars.sort(key=lambda x: -x[2][0])
    print(f"\n  Top 15 by A:")
    print_header()
    for typ, label, sc in all_stars[:15]:
        print_row(f"[{typ}] {label}", sc)
    
    # Full success
    full = [x for x in all_stars if x[2][0] >= 0.70 and x[2][8].get(16,0) >= 0.95]
    if full:
        print(f"\n  ★★★ Full success (A≥70%+B≥95%+B16≥95%+C≥75%): {len(full)}")
        print_header()
        for typ, label, sc in full[:10]:
            print_row(f"[{typ}] {label}", sc)
    else:
        print(f"\n  No ★★★ full success across all phases")
        # Show best C with A≥65%+B≥95%
        high_c = sorted(all_stars, key=lambda x: -(x[2][5] + 0.3*x[2][0]))
        print(f"\n  Best C with A≥65%+B≥95%:")
        print_header()
        for typ, label, sc in high_c[:10]:
            print_row(f"[{typ}] {label}", sc)
else:
    print("\n  NO ★ island found across all three phases!")
    print("  C=0.75 is fundamentally incompatible with sigmoid wall at these α levels")
