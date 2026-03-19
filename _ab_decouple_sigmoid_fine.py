"""Fine scan around the breakthrough sigmoid wall at R_c≈0.25, w≈0.02.

Key discovery: sigmoid wall with sharp cutoff (w=0.02) at R_c=0.25 achieves:
  α=20: A=81%, B=100% (!!), C=0.50
  α=50: A=89.6%, B=95.3%, C=0.00
Now: fine scan (R_c, w, α, η, λ) to find the decoupling island.
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

# R distribution
print("\n  R distribution by family:")
for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
    idx = families_arr == fam
    v = R_arr[idx]
    print(f"    {fam:>7}: mean={np.mean(v):.4f} ± {np.std(v):.4f}  [{np.min(v):.4f}, {np.max(v):.4f}]")

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
    # B at each N
    b_n = {}
    for nv in unique_N:
        i2 = np.where((N_arr == nv) & (families_arr == "Lor2D"))[0]
        ik = np.where((N_arr == nv) & (families_arr == "KR_like"))[0]
        if len(i2) and len(ik):
            b_n[int(nv)] = int(np.sum(F[i2, None] < F[None, ik])) / (len(i2) * len(ik))
    return A, A2, A3, A5, B, C, D, b_n

def print_header():
    print(f"  {'Params':>35} | {'A':>6} {'v2D':>5} {'v3D':>5} {'v5D':>5} | {'B':>6} {'B16':>4} {'B20':>4} {'C':>4} {'D':>5}")
    print("  " + "-" * 90)

def print_row(label, sc):
    A, A2, A3, A5, B, C, D, bn = sc
    b16 = bn.get(16, 0); b20 = bn.get(20, 0)
    tag = " ★" if A >= 0.65 and B >= 0.95 and C >= 0.75 else (" ●" if A >= 0.65 and B >= 0.95 else (" ◆" if A >= 0.65 and B >= 0.90 else ""))
    print(f"  {label:>35} | {100*A:5.1f}% {100*A2:4.0f}% {100*A3:4.0f}% {100*A5:4.0f}% | {100*B:5.1f}% {100*b16:3.0f}% {100*b20:3.0f}% {C:4.2f} {D:5.3f}{tag}")

# ── Understand why sigmoid at R_c=0.25, w=0.02 works ──
print("\n" + "=" * 90)
print("WHY SIGMOID WORKS: sigmoid((R - 0.25)/0.02) activation per family")
print("=" * 90)
for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
    idx = families_arr == fam
    v = R_arr[idx]
    act = sigmoid((v - 0.25) / 0.02)
    print(f"  {fam:>7}: R={np.mean(v):.3f}, σ(R)={np.mean(act):.4f} [{np.min(act):.4f}, {np.max(act):.4f}]")

print("\n  Key: R_c=0.25 is BELOW KR(0.32), so KR gets sigmoid≈1 too!")
print("  But it still works because the ABSOLUTE increase is the same for 2D and KR,")
print("  while 2D already has higher F5 (B was already F5(2D) < F5(KR)).")
print("  The sigmoid adds a CONSTANT to both → doesn't change their relative order!")

# ── Fine scan: vary (R_c, w, α) with base F5 weights ──
print("\n" + "=" * 90)
print("FINE SCAN 1: sigmoid with F5 Bayesian base (γ=0.0004, λ=0.888, η=0.637, κ=0.068)")
print("=" * 90)

gam0, lam0, eta0, kap0 = 0.0004, 0.888, 0.637, 0.068
F5_bay = log_H + gam0*pi_geo - lam0*sigma_hist + eta0*xi_dim + kap0*pi_cg

print_header()
results = []
for Rc in np.arange(0.15, 0.40, 0.02):
    for w in [0.01, 0.015, 0.02, 0.03, 0.05, 0.08]:
        for alpha in [10, 15, 20, 25, 30, 35, 40, 50, 60, 80]:
            wall = alpha * sigmoid((R_arr - Rc) / w)
            F = F5_bay + wall
            sc = compute_scores(F)
            results.append((Rc, w, alpha, sc))

# Filter islands
islands = [r for r in results if r[3][0] >= 0.65 and r[3][4] >= 0.95 and r[3][5] >= 0.75]
islands.sort(key=lambda x: -x[3][0])

if islands:
    print(f"\n  ★★★ DECOUPLING ISLAND FOUND! {len(islands)} points ★★★")
    for r in islands[:20]:
        print_row(f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]}", r[3])
else:
    print("\n  No strict island (A≥65%+B≥95%+C≥75%)")

# A≥65% + B≥95%
feas95 = [r for r in results if r[3][0] >= 0.65 and r[3][4] >= 0.95]
feas95.sort(key=lambda x: -x[3][0])
print(f"\n  Points with A≥65% + B≥95% ({len(feas95)}):")
for r in feas95[:20]:
    print_row(f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]}", r[3])

# ── Fine scan 2: sigmoid + vary η to boost 5D ──
print("\n" + "=" * 90)
print("FINE SCAN 2: sigmoid + varied η (boost 4D vs 5D)")
print("=" * 90)
print_header()

results2 = []
for Rc in [0.21, 0.23, 0.25, 0.27, 0.29]:
    for w in [0.01, 0.015, 0.02, 0.03]:
        for alpha in [15, 20, 25, 30, 35, 40]:
            for eta in [0.637, 1.0, 1.5, 2.0, 2.5, 3.0]:
                wall = alpha * sigmoid((R_arr - Rc) / w)
                F = log_H + gam0*pi_geo - lam0*sigma_hist + eta*xi_dim + kap0*pi_cg + wall
                sc = compute_scores(F)
                results2.append((Rc, w, alpha, eta, sc))

islands2 = [r for r in results2 if r[4][0] >= 0.65 and r[4][4] >= 0.95 and r[4][5] >= 0.75]
islands2.sort(key=lambda x: -x[4][0])

if islands2:
    print(f"\n  ★★★ ISLAND WITH η VARIATION! {len(islands2)} points ★★★")
    for r in islands2[:20]:
        print_row(f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]} η={r[3]:.1f}", r[4])
else:
    print("\n  No strict island with η variation")

# Relax C to 0.50
islands2_c50 = [r for r in results2 if r[4][0] >= 0.65 and r[4][4] >= 0.95 and r[4][5] >= 0.50]
if islands2_c50:
    islands2_c50.sort(key=lambda x: -x[4][0])
    print(f"\n  With C≥0.50 ({len(islands2_c50)} points):")
    for r in islands2_c50[:15]:
        print_row(f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]} η={r[3]:.1f}", r[4])

# A≥65% + B≥95% (no C)
feas2_95 = [r for r in results2 if r[4][0] >= 0.65 and r[4][4] >= 0.95]
feas2_95.sort(key=lambda x: -x[4][0])
print(f"\n  Best A with B≥95% (any C) ({len(feas2_95)} points):")
for r in feas2_95[:15]:
    print_row(f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]} η={r[3]:.1f}", r[4])

# ── Fine scan 3: sigmoid + vary (η, λ) ──
print("\n" + "=" * 90)
print("FINE SCAN 3: sigmoid + varied (η, λ) for best overall")
print("=" * 90)
print_header()

results3 = []
for Rc in [0.23, 0.25, 0.27]:
    for w in [0.015, 0.02, 0.03]:
        for alpha in [18, 22, 26, 30, 35]:
            for eta in [0.637, 1.0, 1.5, 2.0, 3.0]:
                for lam in [0.888, 1.5, 2.0, 3.0, 5.0]:
                    wall = alpha * sigmoid((R_arr - Rc) / w)
                    F = log_H + gam0*pi_geo - lam*sigma_hist + eta*xi_dim + kap0*pi_cg + wall
                    sc = compute_scores(F)
                    results3.append((Rc, w, alpha, eta, lam, sc))

islands3 = [r for r in results3 if r[5][0] >= 0.65 and r[5][4] >= 0.95 and r[5][5] >= 0.75]
if islands3:
    islands3.sort(key=lambda x: -x[5][0])
    print(f"\n  ★★★ FULL ISLAND! {len(islands3)} points ★★★")
    for r in islands3[:20]:
        print_row(f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]} η={r[3]:.1f} λ={r[4]:.1f}", r[5])
else:
    print("\n  No strict island with (η, λ) variation either")

# Best with A≥65% + B≥95%
feas3 = [r for r in results3 if r[5][0] >= 0.65 and r[5][4] >= 0.95]
feas3.sort(key=lambda x: -x[5][0])
print(f"\n  Best A with B≥95% ({len(feas3)} points):")
for r in feas3[:15]:
    print_row(f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]} η={r[3]:.1f} λ={r[4]:.1f}", r[5])

# ── FINAL ANSWER ──
print("\n" + "=" * 90)
print("FINAL ANSWER")
print("=" * 90)

# Collect all islands across all scans
all_islands = []
for r in results:
    if r[3][0] >= 0.65 and r[3][4] >= 0.95:
        all_islands.append(("base", f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]}", r[3]))
for r in results2:
    if r[4][0] >= 0.65 and r[4][4] >= 0.95:
        all_islands.append(("η-var", f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]} η={r[3]:.1f}", r[4]))
for r in results3:
    if r[5][0] >= 0.65 and r[5][4] >= 0.95:
        all_islands.append(("η,λ", f"Rc={r[0]:.2f} w={r[1]:.3f} α={r[2]} η={r[3]:.1f} λ={r[4]:.1f}", r[5]))

if all_islands:
    # Sort by A
    all_islands.sort(key=lambda x: -x[2][0])
    print(f"\n  Total A≥65% + B≥95% points: {len(all_islands)}")
    
    # With C≥0.75
    with_c = [x for x in all_islands if x[2][5] >= 0.75]
    print(f"  Of these, C≥0.75: {len(with_c)}")
    
    # With C≥0.50
    with_c50 = [x for x in all_islands if x[2][5] >= 0.50]
    print(f"  Of these, C≥0.50: {len(with_c50)}")
    
    print(f"\n  Top 10 overall (A≥65% + B≥95%):")
    print_header()
    for typ, label, sc in all_islands[:10]:
        print_row(f"[{typ}] {label}", sc)
    
    if with_c:
        print(f"\n  ★ Top 10 with C≥0.75:")
        print_header()
        for typ, label, sc in sorted(with_c, key=lambda x: -x[2][0])[:10]:
            print_row(f"[{typ}] {label}", sc)
    
    if with_c50:
        print(f"\n  ● Top 10 with C≥0.50:")
        print_header()
        for typ, label, sc in sorted(with_c50, key=lambda x: -x[2][0])[:10]:
            print_row(f"[{typ}] {label}", sc)
    
    # Best balanced: maximize A + B + 0.5*C
    all_islands.sort(key=lambda x: -(x[2][0] + x[2][4] + 0.5*x[2][5]))
    print(f"\n  Best balanced (A + B + 0.5·C):")
    print_header()
    for typ, label, sc in all_islands[:5]:
        print_row(f"[{typ}] {label}", sc)
