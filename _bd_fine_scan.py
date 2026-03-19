"""Fine-grain α_BD scan + per-N pairwise win rates."""
import csv
import numpy as np
from prediction_a_bd_bridge import regenerate_poset, compute_s_bd_ratio, compute_f5, CALIBRATED_WEIGHTS

rows = list(csv.DictReader(open("outputs_unified_functional/raw_features.csv")))
lor_families = {"Lor2D", "Lor3D", "Lor4D", "Lor5D"}
rows = [r for r in rows if r["family"] in lor_families]

# Compute S_BD_ratio for each row
data = []
for r in rows:
    fam, n, rep = r["family"], int(r["N"]), int(r["rep"])
    poset = regenerate_poset(fam, n, rep, 42)
    sbd = compute_s_bd_ratio(poset)
    f5 = compute_f5(r, CALIBRATED_WEIGHTS)
    data.append({"family": fam, "N": n, "rep": rep, "f5": f5, "sbd": sbd})

families_order = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]

# Fine-grain scan
print("FINE-GRAIN SCAN (S_BD_ratio, +sign)")
print(f"{'alpha':>5s} | {'F6(2D)':>8s} {'F6(3D)':>8s} {'F6(4D)':>8s} {'F6(5D)':>8s} | 4<5 4<2 4<3 | Ranking")
print("-" * 90)

for alpha in range(45, 71):
    means = {}
    for fam in families_order:
        vals = [d["f5"] + alpha * d["sbd"] for d in data if d["family"] == fam]
        means[fam] = np.mean(vals)
    c1 = means["Lor4D"] < means["Lor5D"]
    c2 = means["Lor4D"] < means["Lor2D"]
    c3 = means["Lor4D"] < means["Lor3D"]
    rank = sorted(families_order, key=lambda f: means[f])
    mark = " ***" if c1 and c2 and c3 else ""
    c1s = "Y" if c1 else "N"
    c2s = "Y" if c2 else "N"
    c3s = "Y" if c3 else "N"
    print(f"{alpha:5d} | {means['Lor2D']:8.1f} {means['Lor3D']:8.1f} {means['Lor4D']:8.1f} {means['Lor5D']:8.1f} | {c1s:>3s} {c2s:>3s} {c3s:>3s} | {' < '.join(rank)}{mark}")

# Per-N pairwise win rates at multiple alpha values
for alpha in [53, 55, 60, 65]:
    print()
    print("=" * 70)
    print(f"PER-N PAIRWISE WIN RATES at alpha_BD = {alpha}")
    print("=" * 70)

    n_values = sorted(set(d["N"] for d in data))
    for n_val in n_values:
        n_data = [d for d in data if d["N"] == n_val]
        fam_f6 = {}
        for fam in families_order:
            fam_f6[fam] = [d["f5"] + alpha * d["sbd"] for d in n_data if d["family"] == fam]

        print(f"\n  N = {n_val}:")
        for fam in families_order:
            vals = fam_f6[fam]
            print(f"    {fam}: F6 = {np.mean(vals):.2f} +/- {np.std(vals):.2f}")

        # Pairwise: 4D vs each
        for other in ["Lor2D", "Lor3D", "Lor5D"]:
            wins = sum(1 for a, b in zip(fam_f6["Lor4D"], fam_f6[other]) if a < b)
            total = min(len(fam_f6["Lor4D"]), len(fam_f6[other]))
            pct = 100 * wins / total if total > 0 else 0
            print(f"    4D < {other}: {wins}/{total} = {pct:.0f}%")

    # Overall ranking
    print(f"\n  Overall (all N):")
    overall = {}
    for fam in families_order:
        overall[fam] = np.mean([d["f5"] + alpha * d["sbd"] for d in data if d["family"] == fam])
    rank = sorted(families_order, key=lambda f: overall[f])
    print(f"    Ranking: {' < '.join(rank)}")
    for fam in rank:
        print(f"    {fam}: {overall[fam]:.2f}")

# Crossing point analysis
print()
print("=" * 70)
print("CROSSING POINT ANALYSIS")
print("=" * 70)

# Compute exact crossing alphas using per-family mean S_BD and F5
for fam in families_order:
    fam_data = [d for d in data if d["family"] == fam]
    mean_f5 = np.mean([d["f5"] for d in fam_data])
    mean_sbd = np.mean([d["sbd"] for d in fam_data])
    print(f"  {fam}: mean_F5 = {mean_f5:.4f}, mean_S_BD_ratio = {mean_sbd:.4f}")

print()
f5_4d = np.mean([d["f5"] for d in data if d["family"] == "Lor4D"])
sbd_4d = np.mean([d["sbd"] for d in data if d["family"] == "Lor4D"])
for other in ["Lor2D", "Lor3D", "Lor5D"]:
    f5_o = np.mean([d["f5"] for d in data if d["family"] == other])
    sbd_o = np.mean([d["sbd"] for d in data if d["family"] == other])
    denom = sbd_o - sbd_4d
    if abs(denom) > 1e-10:
        alpha_cross = (f5_4d - f5_o) / denom
        print(f"  4D crosses {other} at alpha = {alpha_cross:.2f} (dS_BD = {denom:.4f})")
    else:
        print(f"  4D vs {other}: no crossing (dS_BD ~ 0)")
