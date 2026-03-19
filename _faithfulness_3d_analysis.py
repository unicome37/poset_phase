"""Deeper analysis of 3D vs 4D under faithfulness filter.

Questions:
1. Does the link_fraction gap between 3D and 4D widen with N?
2. Can a combined filter (link + interval) eliminate 3D partially?
3. After filter, what's the F5 gap between 3D and 4D survivors?
4. Can a MILD S_BD energy term close the remaining 3D-4D gap post-filter?
"""
import csv
import numpy as np
from prediction_a_bd_bridge import (
    regenerate_poset, compute_s_bd_ratio, compute_f5, CALIBRATED_WEIGHTS
)
from prediction_a_faithfulness import faithfulness_score

rows = list(csv.DictReader(open("outputs_unified_functional/raw_features.csv")))
lor_families = {"Lor2D", "Lor3D", "Lor4D", "Lor5D"}
rows = [r for r in rows if r["family"] in lor_families]
families_order = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]

# Compute everything
data = []
for r in rows:
    fam, n, rep = r["family"], int(r["N"]), int(r["rep"])
    poset = regenerate_poset(fam, n, rep, 42)
    f5 = compute_f5(r, CALIBRATED_WEIGHTS)
    sbd = compute_s_bd_ratio(poset)
    scores = faithfulness_score(poset)
    data.append({"family": fam, "N": n, "rep": rep, "f5": f5, "sbd": sbd, **scores})

n_values = sorted(set(d["N"] for d in data))

# --- 1. Link fraction gap 3D vs 4D by N ---
print("=" * 70)
print("1. LINK FRACTION GAP: 3D vs 4D by N")
print("=" * 70)
for n_val in n_values:
    lf_3d = [d["int_link_fraction"] for d in data if d["family"] == "Lor3D" and d["N"] == n_val]
    lf_4d = [d["int_link_fraction"] for d in data if d["family"] == "Lor4D" and d["N"] == n_val]
    gap = np.mean(lf_4d) - np.mean(lf_3d)
    # Overlap: what fraction of 3D samples have link_frac > min(4D)?
    min_4d = np.min(lf_4d)
    overlap = sum(1 for v in lf_3d if v >= min_4d) / len(lf_3d) * 100
    print(f"  N={n_val}: 3D={np.mean(lf_3d):.3f}+/-{np.std(lf_3d):.3f}  "
          f"4D={np.mean(lf_4d):.3f}+/-{np.std(lf_4d):.3f}  "
          f"gap={gap:.3f}  overlap={overlap:.0f}%")

print("\n  --> Does gap widen with N? (Needed for large-N extrapolation)")

# --- 2. Two-stage: filter then mild S_BD ---
print("\n" + "=" * 70)
print("2. TWO-STAGE: Faithfulness Filter + Mild S_BD")
print("=" * 70)
print("\nAfter link_fraction >= 0.50 filter (eliminates all 2D):")
print("Then apply F6 = F5 + alpha * S_BD_ratio:")

filtered = [d for d in data if d["int_link_fraction"] >= 0.50]

# Among survivors, what alpha is needed to flip 3D/4D?
f5_3d = np.mean([d["f5"] for d in filtered if d["family"] == "Lor3D"])
f5_4d = np.mean([d["f5"] for d in filtered if d["family"] == "Lor4D"])
sbd_3d = np.mean([d["sbd"] for d in filtered if d["family"] == "Lor3D"])
sbd_4d = np.mean([d["sbd"] for d in filtered if d["family"] == "Lor4D"])
f5_5d = np.mean([d["f5"] for d in filtered if d["family"] == "Lor5D"])
sbd_5d = np.mean([d["sbd"] for d in filtered if d["family"] == "Lor5D"])

print(f"\n  Survivors: 3D={sum(1 for d in filtered if d['family']=='Lor3D')} "
      f"4D={sum(1 for d in filtered if d['family']=='Lor4D')} "
      f"5D={sum(1 for d in filtered if d['family']=='Lor5D')}")
print(f"  F5:  3D={f5_3d:.2f}  4D={f5_4d:.2f}  5D={f5_5d:.2f}")
print(f"  SBD: 3D={sbd_3d:.4f}  4D={sbd_4d:.4f}  5D={sbd_5d:.4f}")

denom_34 = sbd_3d - sbd_4d
denom_45 = sbd_5d - sbd_4d  # This is negative (5D has lower SBD)
alpha_cross_34 = (f5_4d - f5_3d) / denom_34 if abs(denom_34) > 1e-10 else float('inf')
# For 5D: need F5(4D) + alpha*SBD(4D) < F5(5D) + alpha*SBD(5D)
# alpha * (SBD(4D) - SBD(5D)) < F5(5D) - F5(4D)
# alpha < (F5(5D) - F5(4D)) / (SBD(4D) - SBD(5D))
alpha_cross_45 = (f5_5d - f5_4d) / (sbd_4d - sbd_5d) if abs(sbd_4d - sbd_5d) > 1e-10 else float('inf')

print(f"\n  Post-filter crossing points:")
print(f"    4D crosses 3D at alpha = {alpha_cross_34:.2f} (dSBD = {denom_34:.4f})")
print(f"    4D crosses 5D at alpha = {alpha_cross_45:.2f} (dSBD = {sbd_4d-sbd_5d:.4f})")

if alpha_cross_34 < alpha_cross_45:
    print(f"\n  *** WINDOW EXISTS: alpha in [{alpha_cross_34:.1f}, {alpha_cross_45:.1f}] ***")
    alpha_sweet = (alpha_cross_34 + alpha_cross_45) / 2
    print(f"  Sweet spot alpha = {alpha_sweet:.1f}")

    # Check relative contribution at sweet spot
    for fam in ["Lor3D", "Lor4D", "Lor5D"]:
        fam_data = [d for d in filtered if d["family"] == fam]
        f5_m = np.mean([d["f5"] for d in fam_data])
        sbd_m = np.mean([d["sbd"] for d in fam_data])
        contrib = alpha_sweet * sbd_m
        ratio = contrib / f5_m * 100
        f6 = f5_m + alpha_sweet * sbd_m
        print(f"    {fam}: F5={f5_m:.1f} + alpha*SBD={contrib:.1f} ({ratio:.0f}%) = F6={f6:.1f}")

    # Per-N pairwise
    print(f"\n  Per-N pairwise at alpha={alpha_sweet:.1f} (post-filter):")
    for n_val in n_values:
        n_filt = [d for d in filtered if d["N"] == n_val]
        fam_f6 = {}
        for fam in ["Lor3D", "Lor4D", "Lor5D"]:
            fam_f6[fam] = [d["f5"] + alpha_sweet * d["sbd"]
                           for d in n_filt if d["family"] == fam]

        if not fam_f6["Lor3D"] or not fam_f6["Lor4D"]:
            continue

        wins_3d = sum(1 for a, b in zip(fam_f6["Lor4D"], fam_f6["Lor3D"]) if a < b)
        wins_5d = sum(1 for a, b in zip(fam_f6["Lor4D"], fam_f6["Lor5D"]) if a < b)
        n3, n4, n5 = len(fam_f6["Lor3D"]), len(fam_f6["Lor4D"]), len(fam_f6["Lor5D"])

        f6_means = {fam: np.mean(vals) for fam, vals in fam_f6.items() if vals}
        rank = sorted(f6_means.keys(), key=lambda f: f6_means[f])
        print(f"    N={n_val}: 4D<3D={wins_3d}/{min(n4,n3)} 4D<5D={wins_5d}/{min(n4,n5)} | {' < '.join(rank)}")

else:
    print(f"  No window: 5D crosses before 3D")

# --- 3. Pure filter approach: can tighter filter alone separate 3D from 4D? ---
print("\n" + "=" * 70)
print("3. PURE FILTER: Can tighter thresholds separate 3D from 4D?")
print("=" * 70)

# What if we also filter on mean_interval?
for int_thresh in [0.5, 0.4, 0.3, 0.25, 0.2]:
    print(f"\n  link>=0.50 AND mean_interval<={int_thresh}:")
    for fam in families_order:
        fam_data = [d for d in data if d["family"] == fam]
        passed = [d for d in fam_data
                  if d["int_link_fraction"] >= 0.50
                  and d["int_mean_interval"] <= int_thresh]
        if fam_data:
            f5_surv = np.mean([d["f5"] for d in passed]) if passed else float('nan')
            print(f"    {fam}: {len(passed)}/{len(fam_data)} survive, F5={f5_surv:.2f}")

# --- 4. N-scaling of link fraction ---
print("\n" + "=" * 70)
print("4. N-SCALING: Does 2D link fraction trend toward 0?")
print("=" * 70)
for fam in families_order:
    print(f"\n  {fam}:")
    for n_val in n_values:
        fam_n = [d for d in data if d["family"] == fam and d["N"] == n_val]
        lf = [d["int_link_fraction"] for d in fam_n]
        mi = [d["int_mean_interval"] for d in fam_n]
        print(f"    N={n_val}: link_frac={np.mean(lf):.4f}  mean_interval={np.mean(mi):.4f}")

# --- 5. Summary: the two-mechanism picture ---
print("\n" + "=" * 70)
print("5. SUMMARY: Two-Mechanism Dimensional Selection")
print("=" * 70)
print("""
  Mechanism 1: HARD FILTER (manifold faithfulness)
    Diagnostic: link_fraction >= 0.50
    Effect: eliminates Lor2D (0% pass) while preserving 3D/4D/5D
    Physics: 2D sprinklings at moderate N have too-dense causal structure
             (link fraction ~0.25-0.45 vs expected >0.50 for manifold-like)
    N-scaling: 2D link fraction DECREASES with N (more non-links at larger N)
              -> filter becomes STRONGER at large N

  Mechanism 2: UPPER BARRIER (existing F5 functional)
    Already in F5: Lor4D < Lor5D with ~95% pairwise win rate
    Physics: Ξ_d barrier cost increases sharply above d=4

  REMAINING GAP: Lor3D < Lor4D in F5 (~23 gap)
    Options:
    (a) Mild S_BD energy term post-filter (much smaller α needed)
    (b) Tighter faithfulness filter (mean_interval <= 0.3)
    (c) Accept: F5 alone doesn't select 4D over 3D; need BD-type physics
""")
