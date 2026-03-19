"""A–B Decoupling Experiment: test 3 nonlinear R coupling forms.

Goal: find a coupling form where A≥0.65, B≥0.95, C≥0.75 simultaneously.

Form 1: Threshold wall    F = F5 + α·max(0, R - R_c)²
Form 2: Sigmoid wall      F = F5 + α·σ((R - R_c)/w)
Form 3: Entropy modulation F = β_eff(R)·logH + γΠ - λΣ + ηΞ + κΠcg
         where β_eff(R) = β0 + μ·R  or  β0·(1 + μ·R)

Key insight: linear R penalizes Lor2D(0.65) and KR(0.32) proportionally.
Nonlinear forms can create a "wall" that blocks R>R_c (hitting 2D hard)
while leaving KR (R≈0.32) relatively unaffected.
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
    # B at N=16 specifically
    i2_16 = np.where((N_arr == 16) & (families_arr == "Lor2D"))[0]
    ik_16 = np.where((N_arr == 16) & (families_arr == "KR_like"))[0]
    b16 = int(np.sum(F[i2_16, None] < F[None, ik_16])) / (len(i2_16) * len(ik_16)) if len(i2_16) and len(ik_16) else 0
    return A, A2, A3, A5, B, C, D, b16

def print_header():
    print(f"  {'Params':>30} | {'A':>6} {'v2D':>5} {'v3D':>5} {'v5D':>5} | {'B':>6} {'B16':>5} {'C':>4} {'D':>5}")
    print("  " + "-" * 80)

def print_row(label, scores):
    A, A2, A3, A5, B, C, D, b16 = scores
    tag = " ★" if A >= 0.65 and B >= 0.95 and C >= 0.75 else (" ●" if A >= 0.65 and B >= 0.95 else (" ◆" if A >= 0.65 and B >= 0.90 else ""))
    print(f"  {label:>30} | {100*A:5.1f}% {100*A2:4.0f}% {100*A3:4.0f}% {100*A5:4.0f}% | {100*B:5.1f}% {100*b16:4.0f}% {C:4.2f} {D:5.3f}{tag}")

# Fixed F5 base weights (Bayesian)
gam0, lam0, eta0, kap0 = 0.0004, 0.888, 0.637, 0.068
F5_base = log_H + gam0 * pi_geo - lam0 * sigma_hist + eta0 * xi_dim + kap0 * pi_cg

# ── Baseline ──
print("\n" + "=" * 90)
print("BASELINE")
print("=" * 90)
print_header()
print_row("F5 (Bayesian)", compute_scores(F5_base))
print_row("F5 + 47.5·R (linear)", compute_scores(F5_base + 47.5 * R_arr))

# ══════════════════════════════════════════════════════════════════════════
# FORM 1: Threshold Wall  F = F5 + α·max(0, R - R_c)²
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("FORM 1: THRESHOLD WALL  F = F5 + α·max(0, R - R_c)²")
print("  R values: 2D≈0.65, 3D≈0.34, KR≈0.32, 4D≈0.12, 5D≈0.05")
print("  Idea: R_c between KR(0.32) and 3D(0.34) → only penalize 3D/2D, not KR")
print("=" * 90)
print_header()

best1 = None; best1_score = -1
for R_c in [0.20, 0.25, 0.30, 0.33, 0.35, 0.38, 0.40, 0.45, 0.50]:
    for alpha in [50, 100, 200, 500, 1000, 2000, 5000]:
        wall = alpha * np.maximum(0, R_arr - R_c)**2
        F = F5_base + wall
        sc = compute_scores(F)
        A, A2, A3, A5, B, C, D, b16 = sc
        if A >= 0.65 and B >= 0.95 and C >= 0.75:
            if A > best1_score: best1 = (R_c, alpha, sc); best1_score = A
            print_row(f"Rc={R_c:.2f} α={alpha}", sc)
        elif A >= 0.65 and B >= 0.95:
            if best1 is None or A > best1_score: best1 = (R_c, alpha, sc); best1_score = A
            print_row(f"Rc={R_c:.2f} α={alpha}", sc)

if best1 is None:
    # Show best A with B >= 0.95
    best_b95 = None
    for R_c in [0.20, 0.25, 0.30, 0.33, 0.35, 0.38, 0.40, 0.45, 0.50]:
        for alpha in [50, 100, 200, 500, 1000, 2000, 5000]:
            wall = alpha * np.maximum(0, R_arr - R_c)**2
            F = F5_base + wall
            sc = compute_scores(F)
            if sc[4] >= 0.95:
                if best_b95 is None or sc[0] > best_b95[2][0]:
                    best_b95 = (R_c, alpha, sc)
    if best_b95:
        print(f"\n  No A≥65% + B≥95% island. Best A with B≥95%:")
        print_row(f"Rc={best_b95[0]:.2f} α={best_b95[1]}", best_b95[2])
    # Best balanced
    best_bal = None
    for R_c in [0.20, 0.25, 0.30, 0.33, 0.35, 0.38, 0.40, 0.45, 0.50]:
        for alpha in [50, 100, 200, 500, 1000, 2000, 5000]:
            wall = alpha * np.maximum(0, R_arr - R_c)**2
            F = F5_base + wall
            sc = compute_scores(F)
            score = sc[0] + sc[4]
            if best_bal is None or score > best_bal[3]:
                best_bal = (R_c, alpha, sc, score)
    if best_bal:
        print_row(f"Best A+B: Rc={best_bal[0]:.2f} α={best_bal[1]}", best_bal[2])
else:
    print(f"\n  ★ ISLAND FOUND! Rc={best1[0]:.2f}, α={best1[1]}")

# Now do finer scan around best region
print("\n  Fine scan around promising R_c values:")
print_header()
for R_c in np.arange(0.30, 0.42, 0.02):
    for alpha in [100, 200, 400, 800, 1500, 3000, 6000, 10000]:
        wall = alpha * np.maximum(0, R_arr - R_c)**2
        F = F5_base + wall
        sc = compute_scores(F)
        A, A2, A3, A5, B, C, D, b16 = sc
        if A >= 0.55 and B >= 0.93:
            print_row(f"Rc={R_c:.2f} α={alpha}", sc)

# ══════════════════════════════════════════════════════════════════════════
# FORM 2: Sigmoid Wall  F = F5 + α·σ((R - R_c)/w)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("FORM 2: SIGMOID WALL  F = F5 + α·σ((R - R_c)/w)")
print("  σ(x) = 1/(1+exp(-x)). Wall activates at R > R_c.")
print("=" * 90)
print_header()

def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))

best2 = None; best2_score = -1
for R_c in [0.25, 0.30, 0.33, 0.35, 0.38, 0.40, 0.45]:
    for w in [0.02, 0.05, 0.08, 0.10, 0.15]:
        for alpha in [5, 10, 20, 50, 100, 200]:
            wall = alpha * sigmoid((R_arr - R_c) / w)
            F = F5_base + wall
            sc = compute_scores(F)
            A, A2, A3, A5, B, C, D, b16 = sc
            if A >= 0.65 and B >= 0.95 and C >= 0.75:
                if A > best2_score: best2 = (R_c, w, alpha, sc); best2_score = A
            elif A >= 0.65 and B >= 0.95:
                if best2 is None or (best2[3][0] < 0.65 and A > best2[3][0]):
                    best2 = (R_c, w, alpha, sc)

if best2 and best2[3][0] >= 0.65 and best2[3][4] >= 0.95 and best2[3][5] >= 0.75:
    print(f"\n  ★ ISLAND FOUND! Rc={best2[0]:.2f}, w={best2[1]:.2f}, α={best2[2]}")
    print_row(f"Rc={best2[0]:.2f} w={best2[1]:.2f} α={best2[2]}", best2[3])

# Broad report
print("\n  Top results (A≥55% and B≥90%):")
results2 = []
for R_c in [0.25, 0.30, 0.33, 0.35, 0.38, 0.40, 0.45]:
    for w in [0.02, 0.05, 0.08, 0.10, 0.15]:
        for alpha in [5, 10, 20, 50, 100, 200]:
            wall = alpha * sigmoid((R_arr - R_c) / w)
            F = F5_base + wall
            sc = compute_scores(F)
            if sc[0] >= 0.55 and sc[4] >= 0.90:
                results2.append((R_c, w, alpha, sc))

results2.sort(key=lambda x: -(x[3][0] + x[3][4]))
print_header()
for r in results2[:15]:
    print_row(f"Rc={r[0]:.2f} w={r[1]:.2f} α={r[2]}", r[3])

# ══════════════════════════════════════════════════════════════════════════
# FORM 3: Entropy Modulation  F = β_eff(R)·logH + γΠ - λΣ + ηΞ + κΠcg
#   β_eff = β0 + μ·R  (additive)
#   β_eff = β0·(1 + μ·R)  (multiplicative)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("FORM 3: ENTROPY MODULATION  F = (β0 + μ·R)·logH + rest")
print("  Idea: R modulates how much logH counts → different slope per family")
print("=" * 90)

# Form 3a: additive β_eff = β0 + μ·R
print("\n--- Form 3a: β_eff = β0 + μ·R ---")
print_header()
results3a = []
for mu in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0]:
    for eta in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
        beta_eff = 1.0 + mu * R_arr
        F = beta_eff * log_H + gam0 * pi_geo - lam0 * sigma_hist + eta * xi_dim + kap0 * pi_cg
        sc = compute_scores(F)
        results3a.append((mu, eta, sc))

# Filter and sort
feas3a = [(r[0], r[1], r[2]) for r in results3a if r[2][0] >= 0.55 and r[2][4] >= 0.90]
feas3a.sort(key=lambda x: -(x[2][0] + x[2][4]))
for r in feas3a[:15]:
    print_row(f"μ={r[0]:.1f} η={r[1]:.1f}", r[2])

# Check for islands
islands3a = [(r[0], r[1], r[2]) for r in results3a if r[2][0] >= 0.65 and r[2][4] >= 0.95 and r[2][5] >= 0.75]
if islands3a:
    print(f"\n  ★ ISLAND FOUND in Form 3a! {len(islands3a)} points")
    for r in sorted(islands3a, key=lambda x: -x[2][0])[:5]:
        print_row(f"μ={r[0]:.1f} η={r[1]:.1f}", r[2])

# Form 3b: multiplicative β_eff = β0·(1 + μ·R)
print("\n--- Form 3b: β_eff = β0·(1 + μ·R) ---")
print_header()
results3b = []
for mu in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0]:
    for eta in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
        beta_eff = 1.0 * (1.0 + mu * R_arr)
        F = beta_eff * log_H + gam0 * pi_geo - lam0 * sigma_hist + eta * xi_dim + kap0 * pi_cg
        sc = compute_scores(F)
        results3b.append((mu, eta, sc))

feas3b = [(r[0], r[1], r[2]) for r in results3b if r[2][0] >= 0.55 and r[2][4] >= 0.90]
feas3b.sort(key=lambda x: -(x[2][0] + x[2][4]))
for r in feas3b[:15]:
    print_row(f"μ={r[0]:.1f} η={r[1]:.1f}", r[2])

islands3b = [(r[0], r[1], r[2]) for r in results3b if r[2][0] >= 0.65 and r[2][4] >= 0.95 and r[2][5] >= 0.75]
if islands3b:
    print(f"\n  ★ ISLAND FOUND in Form 3b! {len(islands3b)} points")
    for r in sorted(islands3b, key=lambda x: -x[2][0])[:5]:
        print_row(f"μ={r[0]:.1f} η={r[1]:.1f}", r[2])

# Form 3c: combined modulation + separate λ scan
print("\n--- Form 3c: β_eff = 1 + μ·R, vary (μ, η, λ) ---")
print_header()
results3c = []
for mu in [1.0, 2.0, 3.0, 5.0, 8.0, 12.0]:
    for eta in [1.0, 2.0, 3.0, 4.0, 5.0]:
        for lam in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0]:
            beta_eff = 1.0 + mu * R_arr
            F = beta_eff * log_H + gam0 * pi_geo - lam * sigma_hist + eta * xi_dim + kap0 * pi_cg
            sc = compute_scores(F)
            results3c.append((mu, eta, lam, sc))

islands3c = [(r[0], r[1], r[2], r[3]) for r in results3c if r[3][0] >= 0.65 and r[3][4] >= 0.95 and r[3][5] >= 0.75]
if islands3c:
    print(f"\n  ★ ISLAND FOUND in Form 3c! {len(islands3c)} points")
    islands3c.sort(key=lambda x: -x[3][0])
    for r in islands3c[:10]:
        print_row(f"μ={r[0]:.1f} η={r[1]:.1f} λ={r[2]:.1f}", r[3])
else:
    print("  No A≥65%+B≥95%+C≥75% island in Form 3c")
    # Relax: A≥65% + B≥95%
    feas3c_95 = [r for r in results3c if r[3][0] >= 0.65 and r[3][4] >= 0.95]
    if feas3c_95:
        feas3c_95.sort(key=lambda x: -x[3][0])
        print(f"  Best A with B≥95% (no C constraint): {len(feas3c_95)} points")
        for r in feas3c_95[:5]:
            print_row(f"μ={r[0]:.1f} η={r[1]:.1f} λ={r[2]:.1f}", r[3])
    # Relax: A≥55% + B≥93%
    feas3c_bal = [r for r in results3c if r[3][0] >= 0.55 and r[3][4] >= 0.93]
    feas3c_bal.sort(key=lambda x: -(x[3][0] + x[3][4]))
    print(f"\n  Top balanced (A≥55% + B≥93%): {len(feas3c_bal)} points")
    for r in feas3c_bal[:10]:
        print_row(f"μ={r[0]:.1f} η={r[1]:.1f} λ={r[2]:.1f}", r[3])

# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("FINAL SUMMARY: Did any form achieve the decoupling island?")
print("  Target: A ≥ 0.65, B ≥ 0.95, C ≥ 0.75")
print("=" * 90)

any_found = False
for name, islands in [("Threshold wall", [(best1[0], best1[1], best1[2])] if best1 and best1[2][0] >= 0.65 and best1[2][4] >= 0.95 and best1[2][5] >= 0.75 else []),
                       ("Sigmoid wall", [(best2[0], best2[2], best2[3])] if best2 and best2[3][0] >= 0.65 and best2[3][4] >= 0.95 and best2[3][5] >= 0.75 else []),
                       ("Entropy mod (a)", islands3a),
                       ("Entropy mod (b)", islands3b),
                       ("Entropy mod (c)", islands3c)]:
    if islands:
        print(f"  ★ {name}: {len(islands)} island point(s) found!")
        any_found = True
    else:
        print(f"  ✗ {name}: no island")

if not any_found:
    print("\n  CONCLUSION: A–B tension is NOT just a coupling-form issue.")
    print("  All three nonlinear forms fail to achieve the decoupling island.")
    print("  The tension is mechanism-level: R_2D > R_KR is a geometric fact,")
    print("  and any R-based penalty that lifts 2D above 4D also lifts 2D closer to KR.")
    
    # Show best achievable across all forms
    print("\n  Best achievable across all forms:")
    print_header()
    all_best = []
    # Form 1
    for R_c in np.arange(0.20, 0.55, 0.02):
        for alpha in [50, 100, 200, 500, 1000, 2000, 5000, 10000]:
            wall = alpha * np.maximum(0, R_arr - R_c)**2
            F = F5_base + wall
            sc = compute_scores(F)
            all_best.append(("wall", f"Rc={R_c:.2f} α={alpha}", sc))
    # Form 2
    for R_c in [0.25, 0.30, 0.33, 0.35, 0.38, 0.40, 0.45]:
        for w in [0.02, 0.05, 0.08, 0.10, 0.15]:
            for alpha in [5, 10, 20, 50, 100, 200]:
                wall = alpha * sigmoid((R_arr - R_c) / w)
                F = F5_base + wall
                sc = compute_scores(F)
                all_best.append(("sigm", f"Rc={R_c:.2f} w={w:.2f} α={alpha}", sc))
    # Form 3
    for r in results3c:
        all_best.append(("emod", f"μ={r[0]:.1f} η={r[1]:.1f} λ={r[2]:.1f}", r[3]))
    
    # Sort by A+B
    all_best.sort(key=lambda x: -(x[2][0] + x[2][4]))
    print("  Top 5 by A+B:")
    for typ, label, sc in all_best[:5]:
        print_row(f"[{typ}] {label}", sc)
    
    # Best A with B≥95%
    b95 = [x for x in all_best if x[2][4] >= 0.95]
    if b95:
        b95.sort(key=lambda x: -x[2][0])
        print("\n  Best A with B≥95%:")
        for typ, label, sc in b95[:5]:
            print_row(f"[{typ}] {label}", sc)
    
    # Best with A≥65%, B≥90%, C≥0.75
    relaxed = [x for x in all_best if x[2][0] >= 0.65 and x[2][4] >= 0.90 and x[2][5] >= 0.75]
    if relaxed:
        relaxed.sort(key=lambda x: -(x[2][0] + x[2][4]))
        print(f"\n  Best with A≥65% + B≥90% + C≥75% ({len(relaxed)} points):")
        for typ, label, sc in relaxed[:5]:
            print_row(f"[{typ}] {label}", sc)
    else:
        print("\n  No point achieves A≥65% + B≥90% + C≥75% across ANY form!")
else:
    print("\n  SUCCESS: At least one nonlinear form achieves the decoupling island!")
    print("  → A–B tension is a coupling-form issue, not a mechanism-level constraint.")
