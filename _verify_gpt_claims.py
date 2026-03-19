"""Verify GPT's claims about Prediction A and other statistics."""
import csv
import numpy as np

with open("outputs_unified_functional/raw_features.csv") as f:
    rows = list(csv.DictReader(f))

# Compute F from components: F = β·log_H + γ·Π_geo - λ·Σ_hist + η·Ξ_d + κ·Π_cg
# CALIBRATED_WEIGHTS: β=2.0, γ=0.5, λ=1.5, η=0.1, κ=0.05
for r in rows:
    r["F"] = (2.0 * float(r["log_H"]) + 0.5 * float(r["pi_geo"])
              - 1.5 * float(r["sigma_hist"]) + 0.1 * float(r["xi_dim"])
              + 0.05 * float(r["pi_cg"]))

print(f"Total rows: {len(rows)}")

# ── Prediction A: Lor4D vs Lor2D/3D/5D at N>=28 ──
print("\n=== Prediction A: Lor4D vs others (F lower = win) ===")
for n_val in [28, 36]:
    lor4 = [r for r in rows if int(r["N"]) == n_val and r["family"] == "Lor4D"]
    for other_fam in ["Lor2D", "Lor3D", "Lor5D"]:
        others = [r for r in rows if int(r["N"]) == n_val and r["family"] == other_fam]
        wins = sum(1 for i in lor4 for j in others if float(i["F"]) < float(j["F"]))
        total = len(lor4) * len(others)
        print(f"  N={n_val} Lor4D vs {other_fam}: {wins}/{total} = {wins/max(total,1):.1%}")

# Overall A
total_wins = 0
total_pairs = 0
for n_val in [28, 36]:
    lor4 = [r for r in rows if int(r["N"]) == n_val and r["family"] == "Lor4D"]
    for other_fam in ["Lor2D", "Lor3D", "Lor5D"]:
        others = [r for r in rows if int(r["N"]) == n_val and r["family"] == other_fam]
        for i in lor4:
            for j in others:
                total_pairs += 1
                if float(i["F"]) < float(j["F"]):
                    total_wins += 1
print(f"\nOverall A: {total_wins}/{total_pairs} = {total_wins/max(total_pairs,1):.1%}")

# ── Mean F by family and N ──
print("\n=== Mean F by family and N ===")
for n_val in [16, 20, 28, 36]:
    print(f"\n  N={n_val}:")
    for fam in ["KR_like", "Lor2D", "Lor3D", "Lor4D", "Lor5D"]:
        subset = [float(r["F"]) for r in rows if int(r["N"]) == n_val and r["family"] == fam]
        if subset:
            print(f"    {fam:>8s}: F = {np.mean(subset):7.2f} +/- {np.std(subset):5.2f}")

# ── Pi_cg by family ──
print("\n=== Pi_cg by family ===")
for fam in ["KR_like", "Lor2D", "Lor3D", "Lor4D", "Lor5D"]:
    vals = [float(r["pi_cg"]) for r in rows if r["family"] == fam]
    print(f"  {fam:>8s}: pi_cg = {np.mean(vals):.4f} +/- {np.std(vals):.4f}")

lor2d_pcg = [float(r["pi_cg"]) for r in rows if r["family"] == "Lor2D"]
others_pcg = [float(r["pi_cg"]) for r in rows if r["family"] != "Lor2D"]
from scipy import stats
t, p = stats.ttest_ind(lor2d_pcg, others_pcg)
print(f"\n  Lor2D vs Others: t={t:.3f}, p={p:.4f}")
print(f"  Lor2D mean={np.mean(lor2d_pcg):.4f}, Others mean={np.mean(others_pcg):.4f}")

# ── Check which A wins come from which comparison ──
print("\n=== Breakdown of A wins by comparison ===")
for n_val in [28, 36]:
    lor4 = [r for r in rows if int(r["N"]) == n_val and r["family"] == "Lor4D"]
    for other_fam in ["Lor2D", "Lor3D", "Lor5D"]:
        others = [r for r in rows if int(r["N"]) == n_val and r["family"] == other_fam]
        wins = sum(1 for i in lor4 for j in others if float(i["F"]) < float(j["F"]))
        total = len(lor4) * len(others)
        f4_mean = np.mean([float(r["F"]) for r in lor4])
        fo_mean = np.mean([float(r["F"]) for r in others])
        print(f"  N={n_val} Lor4D(F={f4_mean:.1f}) vs {other_fam}(F={fo_mean:.1f}): {wins}/{total}")

# ── rho(Pi_cg, F_4) ──
print("\n=== Correlation Pi_cg with F (all 5 terms) ===")
pcg = np.array([float(r["pi_cg"]) for r in rows])
F_vals = np.array([float(r["F"]) for r in rows])
# F_4 = F - kappa*pi_cg (kappa=0.05)
kappa = 0.05
F4 = F_vals - kappa * pcg
print(f"  rho(Pi_cg, F_5): {np.corrcoef(pcg, F_vals)[0,1]:.4f}")
print(f"  rho(Pi_cg, F_4): {np.corrcoef(pcg, F4)[0,1]:.4f}")

# Window-only correlations
comp = np.array([float(r['comp_frac']) for r in rows])
win_mask = (comp >= 0.30) & (comp <= 0.55)
print(f"  Window rho(Pi_cg, F_4): {np.corrcoef(pcg[win_mask], F4[win_mask])[0,1]:.4f}")
print(f"  Window n={win_mask.sum()}")

# ── t=-517.5 check: where is it still cited? ──
print("\n=== Done. Check manuscript for t=-517.5 re-citations ===")
