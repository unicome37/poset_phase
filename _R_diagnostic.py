"""Diagnose why R term doesn't flip 4D vs 2D ordering in F5."""
import csv, numpy as np
from generators import (
    generate_lorentzian_like_2d, generate_lorentzian_like_3d,
    generate_lorentzian_like_4d, generate_lorentzian_like_5d,
)

data = list(csv.DictReader(open("outputs_unified_functional/raw_features.csv")))
gens = {
    "Lor2D": generate_lorentzian_like_2d, "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d, "Lor5D": generate_lorentzian_like_5d,
}

R_by = {}; F5_by = {}; logH_by = {}
for row in data:
    fam = row["family"]
    if fam not in gens:
        continue
    n, rep = int(row["N"]), int(row["rep"])
    seed = 42 + rep * 1000 + n * 100
    p = gens[fam](n, seed=seed)
    c = p.closure.astype(np.int32)
    ks = c @ c; mask = p.closure
    C0 = int(np.sum(mask & (ks == 0))); total = int(np.sum(mask))
    R = 1.0 - C0 / total if total > 0 else 0.0

    lH = float(row["log_H"])
    f5 = (lH + 0.0004 * float(row["pi_geo"])
          - 0.888 * float(row["sigma_hist"])
          + 0.637 * float(row["xi_dim"])
          + 0.068 * float(row["pi_cg"]))

    key = (fam, n)
    R_by.setdefault(key, []).append(R)
    F5_by.setdefault(key, []).append(f5)
    logH_by.setdefault(key, []).append(lH)

print("=" * 80)
print("DIAGNOSIS: Why R term doesn't flip 4D vs 2D")
print("=" * 80)

print(f"\n  {'fam':>6} {'N':>3} | {'R':>7} | {'log_H':>7} | {'F5':>8} | {'F5-5R':>8} | {'F5-20R':>8} | {'F5-50R':>8}")
print("  " + "-" * 75)
for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]:
    for n in [16, 20, 28, 36]:
        key = (fam, n)
        rm = np.mean(R_by[key])
        lm = np.mean(logH_by[key])
        fm = np.mean(F5_by[key])
        print(f"  {fam:>6} {n:3d} | {rm:7.4f} | {lm:7.2f} | {fm:8.2f} | "
              f"{fm - 5 * rm:8.2f} | {fm - 20 * rm:8.2f} | {fm - 50 * rm:8.2f}")

# Key: what is the GAP between 4D and 2D for each N?
print("\n\nGAP ANALYSIS: F(4D) - F(2D) for each N")
print("-" * 60)
print(f"  {'N':>3} | {'F5 gap':>8} | {'R gap':>8} | {'need α_R':>10} | {'F5-5R gap':>10}")
print("  " + "-" * 55)
for n in [16, 20, 28, 36]:
    f5_4d = np.mean(F5_by[("Lor4D", n)])
    f5_2d = np.mean(F5_by[("Lor2D", n)])
    r_4d = np.mean(R_by[("Lor4D", n)])
    r_2d = np.mean(R_by[("Lor2D", n)])
    gap_f5 = f5_4d - f5_2d
    gap_r = r_4d - r_2d
    # For F6(4D) < F6(2D): F5(4D) - α_R·R(4D) < F5(2D) - α_R·R(2D)
    # → F5(4D) - F5(2D) < α_R·(R(4D) - R(2D))
    # → gap_f5 < α_R · gap_r
    # Since R(4D) < R(2D), gap_r < 0, so need gap_f5 / gap_r > α_R
    # But gap_f5 > 0 and gap_r < 0, so gap_f5/gap_r < 0 → impossible with α_R > 0!
    # Unless we flip the sign: F6 = F5 + α_R·R (positive)
    # Then F5(4D) + α_R·R(4D) < F5(2D) + α_R·R(2D)
    # → gap_f5 < α_R·(R(2D) - R(4D)) = α_R·(-gap_r)
    # → α_R > gap_f5 / (-gap_r) = gap_f5 / (R(2D) - R(4D))
    need_alpha = gap_f5 / (r_2d - r_4d) if (r_2d - r_4d) > 0 else float('inf')
    gap_f6_5 = (f5_4d - 5 * r_4d) - (f5_2d - 5 * r_2d)
    print(f"  {n:3d} | {gap_f5:+8.2f} | {gap_r:+8.4f} | {need_alpha:10.1f} | {gap_f6_5:+10.2f}")

print("\n\nCRITICAL INSIGHT:")
print("  The sign in F6 = F5 - α_R·R is WRONG!")
print("  R is HIGH for 2D and LOW for 4D.")
print("  So -α_R·R makes 2D LOWER (better) and 4D HIGHER (worse) → WRONG direction!")
print("  We need F6 = F5 + α_R·R (POSITIVE sign) to penalize high-R (2D) posets.")
print("  Or equivalently, F6 = F5 - α_R·(1-R) = F5 + α_R·R - α_R")

# Now test with POSITIVE sign
print("\n\n" + "=" * 80)
print("CORRECTED: F6 = F5 + α_R·R (penalize high-R low-d posets)")
print("=" * 80)

families_arr = np.array([d["family"] for d in data])
N_arr = np.array([int(d["N"]) for d in data], dtype=float)
unique_N = sorted(set(N_arr))

# Recompute F5 and R for all
F5_all = np.zeros(len(data))
R_all = np.zeros(len(data))
for i, row in enumerate(data):
    fam = row["family"]
    lH = float(row["log_H"])
    f5 = (lH + 0.0004 * float(row["pi_geo"])
          - 0.888 * float(row["sigma_hist"])
          + 0.637 * float(row["xi_dim"])
          + 0.068 * float(row["pi_cg"]))
    F5_all[i] = f5
    if fam in gens:
        n, rep = int(row["N"]), int(row["rep"])
        seed = 42 + rep * 1000 + n * 100
        p = gens[fam](n, seed=seed)
        c = p.closure.astype(np.int32)
        ks = c @ c; mask = p.closure
        C0 = int(np.sum(mask & (ks == 0))); total = int(np.sum(mask))
        R_all[i] = 1.0 - C0 / total if total > 0 else 0.0
    else:
        R_all[i] = 0.0

def hit_rates(F_vals):
    tw = 0; tt = 0; by = {}
    for fo in ["Lor2D", "Lor3D", "Lor5D"]:
        w = 0; t = 0
        for nv in unique_N:
            i4 = np.where((N_arr == nv) & (families_arr == "Lor4D"))[0]
            io = np.where((N_arr == nv) & (families_arr == fo))[0]
            for a in i4:
                for b in io:
                    if F_vals[a] < F_vals[b]: w += 1
                    t += 1
        by[fo] = (w, t); tw += w; tt += t
    return tw, tt, by

print(f"\n  {'alpha_R':>8} | {'Total':>12} | {'vs2D':>12} | {'vs3D':>12} | {'vs5D':>12}")
print("  " + "-" * 70)

best_a = 0; best_r = 0
for aR in [0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100, 200]:
    F6 = F5_all + aR * R_all  # POSITIVE sign
    tw, tt, by = hit_rates(F6)
    rate = tw / tt if tt else 0
    v2 = by["Lor2D"]; v3 = by["Lor3D"]; v5 = by["Lor5D"]
    print(f"  {aR:8.0f} | {tw:3d}/{tt:3d} ({100*rate:5.1f}%) | "
          f"{v2[0]:3d}/{v2[1]:3d} ({100*v2[0]/v2[1]:5.1f}%) | "
          f"{v3[0]:3d}/{v3[1]:3d} ({100*v3[0]/v3[1]:5.1f}%) | "
          f"{v5[0]:3d}/{v5[1]:3d} ({100*v5[0]/v5[1]:5.1f}%)")
    if rate > best_r: best_r = rate; best_a = aR

print(f"\n  Best α_R = {best_a} → hit rate = {100*best_r:.1f}%")

# Fine scan
print(f"\n  Fine scan around {best_a}:")
for aR in np.arange(max(0, best_a - 10), best_a + 11, 1):
    F6 = F5_all + aR * R_all
    tw, tt, by = hit_rates(F6)
    rate = tw / tt if tt else 0
    v2 = by["Lor2D"]; v3 = by["Lor3D"]; v5 = by["Lor5D"]
    print(f"  {aR:8.0f} | {tw:3d}/{tt:3d} ({100*rate:5.1f}%) | "
          f"{v2[0]:3d}/{v2[1]:3d} ({100*v2[0]/v2[1]:5.1f}%) | "
          f"{v3[0]:3d}/{v3[1]:3d} ({100*v3[0]/v3[1]:5.1f}%) | "
          f"{v5[0]:3d}/{v5[1]:3d} ({100*v5[0]/v5[1]:5.1f}%)")
