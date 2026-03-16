"""Tier-1 / Tier-2 scaling analysis: aw/sqrt(N) as structural discriminant."""
import numpy as np
from scipy import stats
from observables import antichain_width, comparable_fraction
from experiment import FAMILIES
import csv, os
import pandas as pd

# === Tier definition ===
TIER1 = ['lorentzian_like_2d', 'lorentzian_like_3d', 'lorentzian_like_4d',
         'KR_like', 'transitive_percolation']
TIER2 = ['random_layered_k4_uniform', 'random_layered_k6_uniform',
         'random_layered_k6_tapered', 'random_layered_k6_middle_heavy',
         'random_layered_k6_longjump', 'random_layered_k8_uniform',
         'multi_layer_random']

ALL_FAMILIES = TIER1 + TIER2
Ns = [20, 28, 36, 44, 52, 56]
SEEDS = list(range(42, 52))  # 10 seeds

def main():
    rows = []
    total = len(ALL_FAMILIES) * len(Ns) * len(SEEDS)
    done = 0
    print(f"Computing aw/sqrt(N) for {len(ALL_FAMILIES)} families x {len(Ns)} N x {len(SEEDS)} seeds = {total} samples...")

    for fam in ALL_FAMILIES:
        gen = FAMILIES[fam]
        tier = "Tier-1" if fam in TIER1 else "Tier-2"
        for N in Ns:
            for s in SEEDS:
                p = gen(n=N, seed=s)
                aw = antichain_width(p)
                cf = comparable_fraction(p)
                ratio = aw / np.sqrt(N)
                rows.append({
                    "family": fam, "tier": tier, "N": N, "seed": s,
                    "aw": aw, "cf": cf, "aw_over_sqrtN": ratio
                })
                done += 1
                if done % 100 == 0:
                    print(f"  {done}/{total}...", flush=True)

    # Save raw data
    outdir = "outputs_exploratory/tier_scaling_analysis"
    os.makedirs(outdir, exist_ok=True)
    with open(f"{outdir}/aw_scaling_raw.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["family", "tier", "N", "seed", "aw", "cf", "aw_over_sqrtN"])
        w.writeheader()
        w.writerows(rows)

    df = pd.DataFrame(rows)

    # === Summary table ===
    print("\n=== aw/sqrt(N) summary by family (mean +/- std across all N,seed) ===")
    summary = df.groupby(["family", "tier"])["aw_over_sqrtN"].agg(["mean", "std", "count"]).reset_index()
    summary = summary.sort_values("mean")
    summary_rows = []
    for _, r in summary.iterrows():
        print(f'  {r["tier"]:7s}  {r["family"]:35s}  {r["mean"]:.3f} +/- {r["std"]:.3f}  (n={int(r["count"])})')
        summary_rows.append({
            "tier": r["tier"], "family": r["family"],
            "mean_aw_sqrtN": round(r["mean"], 4),
            "std_aw_sqrtN": round(r["std"], 4),
            "n_samples": int(r["count"])
        })
    with open(f"{outdir}/aw_scaling_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["tier", "family", "mean_aw_sqrtN", "std_aw_sqrtN", "n_samples"])
        w.writeheader()
        w.writerows(summary_rows)

    # === Mann-Whitney U tests ===
    print("\n=== Mann-Whitney U: Lor2D vs each Tier-2 family (aw/sqrt(N), all N pooled) ===")
    lor2d_vals = df[df["family"] == "lorentzian_like_2d"]["aw_over_sqrtN"].values
    test_rows = []
    for fam in TIER2:
        fam_vals = df[df["family"] == fam]["aw_over_sqrtN"].values
        u_stat, p_val = stats.mannwhitneyu(lor2d_vals, fam_vals, alternative="less")
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        print(f"  Lor2D vs {fam:35s}  U={u_stat:.0f}  p={p_val:.2e}  {sig}")
        test_rows.append({
            "comparison": f"Lor2D_vs_{fam}",
            "U_statistic": round(u_stat, 1),
            "p_value": f"{p_val:.4e}",
            "significance": sig,
            "lor2d_median": round(np.median(lor2d_vals), 4),
            "other_median": round(np.median(fam_vals), 4)
        })

    print("\n=== Mann-Whitney U: Lor2D vs Tier-1 families ===")
    for fam in ["lorentzian_like_3d", "lorentzian_like_4d", "KR_like", "transitive_percolation"]:
        fam_vals = df[df["family"] == fam]["aw_over_sqrtN"].values
        u_stat, p_val = stats.mannwhitneyu(lor2d_vals, fam_vals, alternative="less")
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        print(f"  Lor2D vs {fam:35s}  U={u_stat:.0f}  p={p_val:.2e}  {sig}")
        test_rows.append({
            "comparison": f"Lor2D_vs_{fam}",
            "U_statistic": round(u_stat, 1),
            "p_value": f"{p_val:.4e}",
            "significance": sig,
            "lor2d_median": round(np.median(lor2d_vals), 4),
            "other_median": round(np.median(fam_vals), 4)
        })

    with open(f"{outdir}/aw_scaling_tests.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["comparison", "U_statistic", "p_value", "significance", "lor2d_median", "other_median"])
        w.writeheader()
        w.writerows(test_rows)

    # === Per-N breakdown ===
    print("\n=== Per-N aw/sqrt(N): Lor2D vs key competitors ===")
    perN_rows = []
    for N in Ns:
        lor2d_n = df[(df["family"] == "lorentzian_like_2d") & (df["N"] == N)]["aw_over_sqrtN"].values
        for fam in ["random_layered_k4_uniform", "random_layered_k6_middle_heavy", "KR_like"]:
            other_n = df[(df["family"] == fam) & (df["N"] == N)]["aw_over_sqrtN"].values
            gap = np.mean(other_n) - np.mean(lor2d_n)
            if len(lor2d_n) >= 3 and len(other_n) >= 3:
                _, p = stats.mannwhitneyu(lor2d_n, other_n, alternative="less")
            else:
                p = float("nan")
            print(f"  N={N:3d}  Lor2D={np.mean(lor2d_n):.3f}  {fam[:20]:20s}={np.mean(other_n):.3f}  gap={gap:+.3f}  p={p:.3e}")
            perN_rows.append({
                "N": N, "lor2d_mean": round(np.mean(lor2d_n), 4),
                "competitor": fam,
                "competitor_mean": round(np.mean(other_n), 4),
                "gap": round(gap, 4), "p_value": f"{p:.4e}"
            })

    with open(f"{outdir}/aw_scaling_perN.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["N", "lor2d_mean", "competitor", "competitor_mean", "gap", "p_value"])
        w.writeheader()
        w.writerows(perN_rows)

    # === Tier-1 vs Tier-2 aggregate test ===
    print("\n=== Aggregate Tier-1 (geometric) vs Tier-2 (layered null-model) ===")
    # Only Lor2D from Tier-1 for the cleanest comparison
    t1_vals = df[df["family"] == "lorentzian_like_2d"]["aw_over_sqrtN"].values
    t2_vals = df[df["tier"] == "Tier-2"]["aw_over_sqrtN"].values
    u, p = stats.mannwhitneyu(t1_vals, t2_vals, alternative="less")
    print(f"  Lor2D (n={len(t1_vals)}) vs all Tier-2 (n={len(t2_vals)}): U={u:.0f}, p={p:.2e}")
    print(f"  Lor2D median={np.median(t1_vals):.4f}, Tier-2 median={np.median(t2_vals):.4f}")

    print(f"\nAll outputs saved to {outdir}/")
    print("DONE")

if __name__ == "__main__":
    main()
