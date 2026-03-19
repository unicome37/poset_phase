"""Test nonlinear S_BD transforms: log, sqrt, and power-law to find
a more balanced sixth term that doesn't sledgehammer 2D."""
import csv
import numpy as np
from prediction_a_bd_bridge import regenerate_poset, compute_s_bd_ratio, compute_f5, CALIBRATED_WEIGHTS

rows = list(csv.DictReader(open("outputs_unified_functional/raw_features.csv")))
lor_families = {"Lor2D", "Lor3D", "Lor4D", "Lor5D"}
rows = [r for r in rows if r["family"] in lor_families]
families_order = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]

# Compute raw data
data = []
for r in rows:
    fam, n, rep = r["family"], int(r["N"]), int(r["rep"])
    poset = regenerate_poset(fam, n, rep, 42)
    sbd = compute_s_bd_ratio(poset)
    f5 = compute_f5(r, CALIBRATED_WEIGHTS)
    data.append({"family": fam, "N": n, "rep": rep, "f5": f5, "sbd_raw": sbd})

# Define nonlinear transforms
def log_transform(x):
    return np.log1p(x)  # log(1+x), monotonic, compresses large values

def sqrt_transform(x):
    return np.sqrt(x)

def power_transform(x, p=0.3):
    return x ** p  # x^0.3, strong compression

# Apply transforms
for d in data:
    d["sbd_log"] = log_transform(d["sbd_raw"])
    d["sbd_sqrt"] = sqrt_transform(d["sbd_raw"])
    d["sbd_p03"] = power_transform(d["sbd_raw"], 0.3)

# Report transformed means
print("=" * 80)
print("TRANSFORMED S_BD STATISTICS")
print("=" * 80)
for key, label in [("sbd_raw", "raw"), ("sbd_log", "log(1+x)"),
                    ("sbd_sqrt", "sqrt(x)"), ("sbd_p03", "x^0.3")]:
    print(f"\n  {label}:")
    print(f"    {'Family':<8s} {'Mean':>8s} {'Std':>8s}")
    for fam in families_order:
        vals = [d[key] for d in data if d["family"] == fam]
        print(f"    {fam:<8s} {np.mean(vals):8.4f} {np.std(vals):8.4f}")
    # Ratio: how much larger is 2D vs 4D?
    m2d = np.mean([d[key] for d in data if d["family"] == "Lor2D"])
    m4d = np.mean([d[key] for d in data if d["family"] == "Lor4D"])
    m3d = np.mean([d[key] for d in data if d["family"] == "Lor3D"])
    print(f"    2D/4D ratio: {m2d/m4d:.1f}x | 3D/4D ratio: {m3d/m4d:.1f}x")

# Scan each transform
print("\n" + "=" * 80)
print("CROSSING SCAN — ALL TRANSFORMS")
print("=" * 80)

for key, label in [("sbd_raw", "raw"), ("sbd_log", "log(1+x)"),
                    ("sbd_sqrt", "sqrt(x)"), ("sbd_p03", "x^0.3")]:
    # Compute exact crossing points
    f5_means = {}
    sbd_means = {}
    for fam in families_order:
        f5_means[fam] = np.mean([d["f5"] for d in data if d["family"] == fam])
        sbd_means[fam] = np.mean([d[key] for d in data if d["family"] == fam])

    print(f"\n--- Transform: {label} ---")
    crossings = {}
    for other in ["Lor2D", "Lor3D", "Lor5D"]:
        denom = sbd_means[other] - sbd_means["Lor4D"]
        if abs(denom) > 1e-10:
            alpha_cross = (f5_means["Lor4D"] - f5_means[other]) / denom
            crossings[other] = alpha_cross
            print(f"  4D crosses {other} at alpha = {alpha_cross:.2f}")
        else:
            print(f"  4D vs {other}: no crossing")

    # Check if window exists (4D crosses 2D AND 3D before crossing 5D)
    if "Lor2D" in crossings and "Lor3D" in crossings:
        alpha_min = max(crossings["Lor2D"], crossings["Lor3D"])
        alpha_max = crossings.get("Lor5D", float("inf"))
        if alpha_min < alpha_max:
            print(f"  ** WINDOW: alpha in [{alpha_min:.2f}, {alpha_max:.2f}]")

            # Check relative contribution at window entry
            alpha_test = alpha_min * 1.05  # just inside window
            for fam in families_order:
                sbd_contrib = alpha_test * sbd_means[fam]
                ratio = sbd_contrib / f5_means[fam] * 100
                print(f"     {fam}: F5={f5_means[fam]:.1f}, alpha*S_BD={sbd_contrib:.1f} ({ratio:.0f}%)")

            # Per-N pairwise at sweet spot
            alpha_sweet = alpha_min * 1.2
            print(f"\n  Per-N pairwise at alpha={alpha_sweet:.1f}:")
            for n_val in sorted(set(d["N"] for d in data)):
                n_data = [d for d in data if d["N"] == n_val]
                fam_f6 = {}
                for fam in families_order:
                    fam_f6[fam] = [d["f5"] + alpha_sweet * d[key]
                                   for d in n_data if d["family"] == fam]

                wins_2d = sum(1 for a, b in zip(fam_f6["Lor4D"], fam_f6["Lor2D"]) if a < b)
                wins_3d = sum(1 for a, b in zip(fam_f6["Lor4D"], fam_f6["Lor3D"]) if a < b)
                wins_5d = sum(1 for a, b in zip(fam_f6["Lor4D"], fam_f6["Lor5D"]) if a < b)
                total = len(fam_f6["Lor4D"])

                f6_means = {fam: np.mean(vals) for fam, vals in fam_f6.items()}
                rank = sorted(families_order, key=lambda f: f6_means[f])
                print(f"    N={n_val}: 4D<2D={wins_2d}/{total} 4D<3D={wins_3d}/{total} "
                      f"4D<5D={wins_5d}/{total}  |  {' < '.join(rank)}")
        else:
            print(f"  No window: 4D crosses 5D ({alpha_max:.2f}) before crossing 3D ({alpha_min:.2f})")
    else:
        print(f"  Insufficient crossing data")
