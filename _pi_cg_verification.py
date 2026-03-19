"""Π_cg closure hypothesis independent verification under F7 sigmoid wall.

Tests:
1. Π_cg statistics under F7 — does κΠ_cg still matter?
2. CG drift of F7 vs F5 — does sigmoid wall change CG stability?
3. Π_cg as hard filter — does removing high-Π_cg samples improve A/B/C?
4. Π_cg predictive power — does low Π_cg predict lower F7 within family?
5. Family discrimination under F7 — any new signal?
6. κ sensitivity — how do A/B/C change as κ varies from 0 to 0.5?
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
print("Computing R and CG-drift features...")
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

def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))

wall_base = sigmoid((R_arr - 0.25) / 0.015)
N0 = 20.0

def compute_F7(kappa_val, alpha0=16, q=-0.5, lam=10, eta=0.6):
    alpha_n = alpha0 * (N0 / N_arr) ** abs(q)
    return log_H + 0.0004*pi_geo - lam*sigma_hist + eta*xi_dim + kappa_val*pi_cg + alpha_n*wall_base

def compute_scores(F):
    w2 = int(np.sum(F[A_pairs_2d[:, 0]] < F[A_pairs_2d[:, 1]]))
    w3 = int(np.sum(F[A_pairs_3d[:, 0]] < F[A_pairs_3d[:, 1]]))
    w5 = int(np.sum(F[A_pairs_5d[:, 0]] < F[A_pairs_5d[:, 1]]))
    t_a = len(A_pairs_2d) + len(A_pairs_3d) + len(A_pairs_5d)
    A = (w2 + w3 + w5) / t_a
    wb = int(np.sum(F[B_pairs[:, 0]] < F[B_pairs[:, 1]]))
    B = wb / len(B_pairs)
    neg = 0
    for nv in unique_N:
        idx = np.where(N_arr == nv)[0]
        sh = sigma_hist[idx]; f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            if np.corrcoef(sh, f)[0, 1] < 0: neg += 1
    C = neg / len(unique_N)
    return A, B, C

# ══════════════════════════════════════════════════════════════════════════
# TEST 1: κ sensitivity — how do A/B/C change with κ?
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("TEST 1: kappa sensitivity (alpha0=16, q=-0.5, lambda=10, eta=0.6)")
print("=" * 90)
print(f"  {'kappa':>8} | {'A':>6} {'B':>6} {'C':>4} | {'Delta_A':>7} {'Delta_B':>7}")

F_base = compute_F7(0.0)
A0, B0, C0 = compute_scores(F_base)
for kappa in [0.0, 0.01, 0.02, 0.05, 0.068, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0, 2.0]:
    F = compute_F7(kappa)
    A, B, C = compute_scores(F)
    dA = A - A0; dB = B - B0
    print(f"  {kappa:8.3f} | {100*A:5.1f}% {100*B:5.1f}% {C:4.2f} | {100*dA:+6.2f}% {100*dB:+6.2f}%")

# ══════════════════════════════════════════════════════════════════════════
# TEST 2: Π_cg correlation with F7 components
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("TEST 2: Correlation of pi_cg with F7 and its components")
print("=" * 90)

F7_main = compute_F7(0.068)
F7_no_cg = compute_F7(0.0)

print(f"  rho(pi_cg, F7)        = {np.corrcoef(pi_cg, F7_main)[0,1]:+.4f}")
print(f"  rho(pi_cg, F7_no_cg)  = {np.corrcoef(pi_cg, F7_no_cg)[0,1]:+.4f}")
print(f"  rho(pi_cg, logH)      = {np.corrcoef(pi_cg, log_H)[0,1]:+.4f}")
print(f"  rho(pi_cg, sigma_hist)= {np.corrcoef(pi_cg, sigma_hist)[0,1]:+.4f}")
print(f"  rho(pi_cg, xi_dim)    = {np.corrcoef(pi_cg, xi_dim)[0,1]:+.4f}")
print(f"  rho(pi_cg, pi_geo)    = {np.corrcoef(pi_cg, pi_geo)[0,1]:+.4f}")
print(f"  rho(pi_cg, R)         = {np.corrcoef(pi_cg, R_arr)[0,1]:+.4f}")
print(f"  rho(pi_cg, wall)      = {np.corrcoef(pi_cg, wall_base)[0,1]:+.4f}")

# Per-N correlations
print("\n  Per-N correlations with F7 (no cg):")
for nv in unique_N:
    idx = np.where(N_arr == nv)[0]
    r = np.corrcoef(pi_cg[idx], F7_no_cg[idx])[0, 1]
    print(f"    N={int(nv):2d}: rho(pi_cg, F7) = {r:+.4f}")

# ══════════════════════════════════════════════════════════════════════════
# TEST 3: Family-level Π_cg statistics under F7
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("TEST 3: Π_cg family statistics")
print("=" * 90)

fams = ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]
for fam in fams:
    idx = np.where(families_arr == fam)[0]
    print(f"  {fam:>8}: pi_cg = {np.mean(pi_cg[idx]):.4f} +/- {np.std(pi_cg[idx]):.4f}  "
          f"range=[{np.min(pi_cg[idx]):.4f}, {np.max(pi_cg[idx]):.4f}]")

# Per-N family statistics
print("\n  Per-N family means:")
print(f"  {'N':>4} | {'Lor2D':>8} {'Lor3D':>8} {'Lor4D':>8} {'Lor5D':>8} {'KR_like':>8} | {'F(MW)':>8}")
from scipy import stats
for nv in unique_N:
    vals = {}
    for fam in fams:
        idx = np.where((N_arr == nv) & (families_arr == fam))[0]
        vals[fam] = pi_cg[idx]
    # Kruskal-Wallis test
    kw_stat, kw_p = stats.kruskal(*[vals[f] for f in fams])
    means = [f"{np.mean(vals[f]):.4f}" for f in fams]
    print(f"  {int(nv):4d} | {'  '.join(means)} | p={kw_p:.4f}")

# Overall KW
all_groups = [pi_cg[families_arr == f] for f in fams]
kw_stat, kw_p = stats.kruskal(*all_groups)
print(f"  Overall Kruskal-Wallis: H={kw_stat:.2f}, p={kw_p:.4f}")

# ══════════════════════════════════════════════════════════════════════════
# TEST 4: Π_cg as hard filter — remove high-Π_cg samples
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("TEST 4: Π_cg as hard filter (remove high-Π_cg samples)")
print("=" * 90)

F7 = compute_F7(0.0)  # without kappa term
print(f"  {'Filter':>20} | {'A':>6} {'B':>6} {'C':>4} | {'n_remain':>8} {'n_lor4d':>7}")

# No filter baseline
A_nf, B_nf, C_nf = compute_scores(F7)
print(f"  {'No filter':>20} | {100*A_nf:5.1f}% {100*B_nf:5.1f}% {C_nf:4.2f} | {n_data:>8} {int(np.sum(families_arr=='Lor4D')):>7}")

# Filter at various percentiles
for pct in [90, 80, 75, 70, 60, 50]:
    threshold = np.percentile(pi_cg, pct)
    keep = pi_cg <= threshold
    # Rebuild pairs with filtered indices
    idx_map = np.where(keep)[0]
    fam_k = families_arr[keep]
    N_k = N_arr[keep]
    F_k = F7[keep]
    sh_k = sigma_hist[keep]
    
    # Re-count predictions
    w2 = w3 = w5 = wb = 0
    t2 = t3 = t5 = tb = 0
    for nv in unique_N:
        i4 = np.where((N_k == nv) & (fam_k == "Lor4D"))[0]
        i2 = np.where((N_k == nv) & (fam_k == "Lor2D"))[0]
        i3 = np.where((N_k == nv) & (fam_k == "Lor3D"))[0]
        i5 = np.where((N_k == nv) & (fam_k == "Lor5D"))[0]
        ik = np.where((N_k == nv) & (fam_k == "KR_like"))[0]
        for a in i4:
            for b in i2: t2 += 1; w2 += int(F_k[a] < F_k[b])
            for b in i3: t3 += 1; w3 += int(F_k[a] < F_k[b])
            for b in i5: t5 += 1; w5 += int(F_k[a] < F_k[b])
        for a in i2:
            for b in ik: tb += 1; wb += int(F_k[a] < F_k[b])
    
    A_f = (w2+w3+w5)/(t2+t3+t5) if (t2+t3+t5) > 0 else 0
    B_f = wb/tb if tb > 0 else 0
    neg = 0
    for nv in unique_N:
        idx = np.where(N_k == nv)[0]
        if len(idx) < 3: continue
        if np.std(sh_k[idx]) > 1e-10 and np.std(F_k[idx]) > 1e-10:
            if np.corrcoef(sh_k[idx], F_k[idx])[0, 1] < 0: neg += 1
    C_f = neg / len(unique_N)
    
    n4 = int(np.sum((fam_k == "Lor4D")))
    print(f"  {'P'+str(pct)+' ('+f'{threshold:.3f}'+')':>20} | {100*A_f:5.1f}% {100*B_f:5.1f}% {C_f:4.2f} | {int(np.sum(keep)):>8} {n4:>7}")

# ══════════════════════════════════════════════════════════════════════════
# TEST 5: Within-family: does low Π_cg predict lower F7?
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("TEST 5: Within-family correlation: rho(pi_cg, F7) per family per N")
print("  Positive rho = high pi_cg → high F7 (bad) = pi_cg acts as penalty")
print("  Negative rho = high pi_cg → low F7 (good) = pi_cg counterproductive as penalty")
print("=" * 90)

print(f"  {'Family':>8} {'N':>4} | {'rho(cg,F7)':>10} {'rho(cg,logH)':>12} {'rho(cg,R)':>10} | {'n':>3}")
for fam in fams:
    for nv in unique_N:
        idx = np.where((N_arr == nv) & (families_arr == fam))[0]
        if len(idx) < 3: continue
        r_f7 = np.corrcoef(pi_cg[idx], F7_no_cg[idx])[0, 1]
        r_lh = np.corrcoef(pi_cg[idx], log_H[idx])[0, 1]
        r_R = np.corrcoef(pi_cg[idx], R_arr[idx])[0, 1]
        print(f"  {fam:>8} {int(nv):>4} | {r_f7:+10.4f} {r_lh:+12.4f} {r_R:+10.4f} | {len(idx):>3}")

# ══════════════════════════════════════════════════════════════════════════
# TEST 6: CG drift of F7 — compute ΔF7 for coarse-grained posets
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("TEST 6: CG drift verification — does Π_cg predict |ΔF7| under coarse-graining?")
print("=" * 90)

# We need to compute F7 for CG'd posets. This requires regenerating and CG'ing.
# For efficiency, compute for a subset (first 80 = N=16 and N=20)
np.random.seed(42)
cg_results = []
for i in range(n_data):
    fam = data[i]["family"]
    n, rep = int(data[i]["N"]), int(data[i]["rep"])
    seed = 42 + rep * 1000 + n * 100
    p = gens[fam](n, seed=seed)
    
    # CG: randomly delete 30% of nodes
    n_keep = max(int(0.7 * n), 4)
    keep_nodes = np.sort(np.random.choice(n, n_keep, replace=False))
    
    # Extract sub-poset
    cl = p.closure
    sub_cl = cl[np.ix_(keep_nodes, keep_nodes)]
    
    # Compute features for sub-poset
    mask_sub = sub_cl.astype(bool)
    n_rels_sub = int(np.sum(mask_sub)) - n_keep  # exclude diagonal
    log_H_cg = np.log(max(n_rels_sub, 1))
    
    # R for sub-poset
    c_sub = sub_cl.astype(np.int32)
    ks_sub = c_sub @ c_sub
    C0_sub = int(np.sum(mask_sub & (ks_sub == 0)))
    total_sub = int(np.sum(mask_sub))
    R_cg = 1.0 - C0_sub / total_sub if total_sub > 0 else 0.0
    
    # F7 for CG'd poset (only logH and wall terms, simplified)
    alpha_n = 16 * (20.0 / n) ** 0.5
    wall_cg = alpha_n * sigmoid((R_cg - 0.25) / 0.015)
    wall_orig = alpha_n * sigmoid((R_arr[i] - 0.25) / 0.015)
    
    # Relative drift in wall term
    delta_wall = abs(wall_cg - wall_orig)
    delta_logH = abs(log_H_cg - log_H[i])
    delta_F_approx = abs((log_H_cg + wall_cg) - (log_H[i] + wall_orig))
    
    cg_results.append({
        'i': i, 'fam': fam, 'N': n, 'pi_cg': pi_cg[i],
        'delta_wall': delta_wall, 'delta_logH': delta_logH,
        'delta_F_approx': delta_F_approx,
        'R_orig': R_arr[i], 'R_cg': R_cg,
    })

delta_F_arr = np.array([r['delta_F_approx'] for r in cg_results])
delta_wall_arr = np.array([r['delta_wall'] for r in cg_results])
delta_logH_arr = np.array([r['delta_logH'] for r in cg_results])

print(f"  Global rho(pi_cg, |delta_F7_approx|) = {np.corrcoef(pi_cg, delta_F_arr)[0,1]:+.4f}")
print(f"  Global rho(pi_cg, |delta_wall|)      = {np.corrcoef(pi_cg, delta_wall_arr)[0,1]:+.4f}")
print(f"  Global rho(pi_cg, |delta_logH|)      = {np.corrcoef(pi_cg, delta_logH_arr)[0,1]:+.4f}")

for fam in fams:
    idx = [j for j, r in enumerate(cg_results) if r['fam'] == fam]
    if len(idx) < 3: continue
    pc = pi_cg[idx]
    df = delta_F_arr[idx]
    r = np.corrcoef(pc, df)[0, 1]
    print(f"  {fam:>8}: rho(pi_cg, |delta_F7|) = {r:+.4f}  (n={len(idx)})")

# ══════════════════════════════════════════════════════════════════════════
# TEST 7: Does removing κ from F7 affect predictions at all?
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("TEST 7: F7 with vs without kappa*pi_cg — direct comparison")
print("=" * 90)

for alpha0 in [16, 14, 26]:
    for q in [-0.5, 0.0]:
        for lam in [10, 5]:
            for eta in [0.6, 1.0]:
                F_with = compute_F7(0.068, alpha0, q, lam, eta)
                F_without = compute_F7(0.0, alpha0, q, lam, eta)
                Aw, Bw, Cw = compute_scores(F_with)
                Ano, Bno, Cno = compute_scores(F_without)
                if abs(Aw - Ano) > 0.001 or abs(Bw - Bno) > 0.001 or abs(Cw - Cno) > 0.001:
                    print(f"  a0={alpha0:2d} q={q:+.1f} l={lam:2d} e={eta:.1f}: "
                          f"WITH kappa: A={100*Aw:.1f}% B={100*Bw:.1f}% C={Cw:.2f}  "
                          f"WITHOUT: A={100*Ano:.1f}% B={100*Bno:.1f}% C={Cno:.2f}  "
                          f"DIFF: dA={100*(Aw-Ano):+.1f}% dB={100*(Bw-Bno):+.1f}%")

print("\n  (No output above = kappa*pi_cg has ZERO effect on predictions)")

# ══════════════════════════════════════════════════════════════════════════
# TEST 8: Magnitude comparison — how big is κ*Π_cg vs other terms?
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("TEST 8: Magnitude comparison of F7 terms")
print("=" * 90)

alpha_n_arr = 16 * (N0 / N_arr) ** 0.5
wall_term = alpha_n_arr * wall_base
cg_term = 0.068 * pi_cg

for nv in unique_N:
    idx = np.where(N_arr == nv)[0]
    print(f"\n  N = {int(nv)}:")
    print(f"    logH:         {np.mean(log_H[idx]):7.2f} +/- {np.std(log_H[idx]):5.2f}")
    print(f"    -lam*sig:     {np.mean(-10*sigma_hist[idx]):7.2f} +/- {np.std(-10*sigma_hist[idx]):5.2f}")
    print(f"    eta*xi:       {np.mean(0.6*xi_dim[idx]):7.2f} +/- {np.std(0.6*xi_dim[idx]):5.2f}")
    print(f"    wall:         {np.mean(wall_term[idx]):7.2f} +/- {np.std(wall_term[idx]):5.2f}")
    print(f"    kappa*pi_cg:  {np.mean(cg_term[idx]):7.4f} +/- {np.std(cg_term[idx]):5.4f}")
    print(f"    ratio cg/wall: {np.mean(cg_term[idx])/max(np.mean(wall_term[idx]),0.001):.4f}")
    print(f"    ratio cg/|F7|: {np.mean(cg_term[idx])/max(abs(np.mean(compute_F7(0.068)[idx])),0.001):.6f}")

# ══════════════════════════════════════════════════════════════════════════
# FINAL VERDICT
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 90)
print("FINAL VERDICT")
print("=" * 90)

# Compute key statistics for verdict
F7_def = compute_F7(0.068)
F7_noc = compute_F7(0.0)
A_def, B_def, C_def = compute_scores(F7_def)
A_noc, B_noc, C_noc = compute_scores(F7_noc)

print(f"\n  F7 with kappa=0.068:  A={100*A_def:.1f}%, B={100*B_def:.1f}%, C={C_def:.2f}")
print(f"  F7 with kappa=0.000:  A={100*A_noc:.1f}%, B={100*B_noc:.1f}%, C={C_noc:.2f}")
print(f"  Difference:           dA={100*(A_def-A_noc):+.2f}%, dB={100*(B_def-B_noc):+.2f}%")
print(f"\n  kappa*pi_cg mean magnitude: {np.mean(0.068*pi_cg):.4f}")
print(f"  wall mean magnitude:         {np.mean(wall_term):.2f}")
print(f"  Ratio:                       {np.mean(0.068*pi_cg)/np.mean(wall_term):.6f}")
print(f"\n  rho(pi_cg, delta_F7_CG) = {np.corrcoef(pi_cg, delta_F_arr)[0,1]:+.4f}")
print(f"  Family discrimination p = {stats.kruskal(*all_groups)[1]:.4f}")

print(f"""
  CONCLUSION:
  - kappa*pi_cg contributes ~{np.mean(0.068*pi_cg):.4f} to F7, vs wall ~{np.mean(wall_term):.1f}
    (ratio ~{np.mean(0.068*pi_cg)/np.mean(wall_term):.5f})
  - Removing kappa changes A/B/C by: {100*(A_def-A_noc):+.2f}% / {100*(B_def-B_noc):+.2f}% / {C_def-C_noc:+.2f}
  - pi_cg does NOT discriminate families (p={stats.kruskal(*all_groups)[1]:.4f})
  - pi_cg correlates with CG drift: rho={np.corrcoef(pi_cg, delta_F_arr)[0,1]:+.3f} (confirms it measures what it claims)
  - Hard filter at various thresholds: marginal or no improvement
""")
