"""
Prediction A — F10 4D<5D Weakness Diagnosis
=============================================
At N=20, 4D<5D win rate is only 62%. Diagnose why and explore fixes.

Hypothesis: d_eff(5D)≈4.4 is too close to d*=4.1.
  - (d_eff(4D)−d*)² = (4.0−4.1)² = 0.01
  - (d_eff(5D)−d*)² = (4.4−4.1)² = 0.09
  - Δ = 0.08 per element → 0.08×20 = 1.6 at N=20
  - But ΔlogH(5D−4D) ≈ 2.6 at N=20 → d_eff penalty insufficient

Possible fixes:
  1. Higher γ (but may hurt other pairs)
  2. d* closer to 4.0 (increases 4D<5D gap but may hurt 4D<3D)
  3. Asymmetric well: steeper on d>d* side
  4. Combine with ξ_d or wall terms
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats


def sigmoid(x):
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def load_csv(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            d = {}
            for k, v in row.items():
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
            d["N"] = int(d["N"])
            d["rep"] = int(d["rep"])
            rows.append(d)
    return rows


def compute_F10(row, d_star=4.1, gamma=1.0, lam=10.0, eta=0.0,
                alpha0=16.0, q=0.5, Rc=0.25, w=0.015, N0=20.0,
                gamma_hi=None):
    """F10 with optional asymmetric well (gamma_hi for d_eff > d*)."""
    N = row["N"]
    d_eff = row["d_eff"]
    
    if gamma_hi is not None and d_eff > d_star:
        dim_well = gamma_hi * N * (d_eff - d_star) ** 2
    else:
        dim_well = gamma * N * (d_eff - d_star) ** 2
    
    alpha_N = alpha0 * (N0 / max(N, 1)) ** abs(q)
    wall = alpha_N * sigmoid((row["R"] - Rc) / w)
    
    return (row["log_H"] + dim_well - lam * row["sigma_hist"] 
            + eta * row["xi_dim"] + wall)


def main():
    csv_path = "outputs_d_recovery/prediction_c_f7_large_n.csv"
    rows = load_csv(csv_path)
    n_values = sorted(set(r["N"] for r in rows))
    all_families = sorted(set(r["family"] for r in rows))
    lor_families = [f for f in all_families if f.startswith("Lor")]

    by_nf = defaultdict(list)
    for r in rows:
        by_nf[(r["N"], r["family"])].append(r)

    report = []
    report.append("# F10 4D<5D Weakness Diagnosis\n")

    # ── 1. Per-sample d_eff distributions ──
    report.append("## 1. d_eff Distributions at N=20\n")
    for f in ["Lor4D", "Lor5D"]:
        vals = by_nf.get((20, f), [])
        deff = [r["d_eff"] for r in vals]
        report.append(f"**{f}**: mean={np.mean(deff):.3f}, std={np.std(deff):.3f}, "
                      f"min={np.min(deff):.3f}, max={np.max(deff):.3f}")

    # d_eff overlap analysis
    d4 = [r["d_eff"] for r in by_nf.get((20, "Lor4D"), [])]
    d5 = [r["d_eff"] for r in by_nf.get((20, "Lor5D"), [])]
    report.append(f"\nd_eff(4D) range: [{min(d4):.2f}, {max(d4):.2f}]")
    report.append(f"d_eff(5D) range: [{min(d5):.2f}, {max(d5):.2f}]")
    overlap = max(0, min(max(d4), max(d5)) - max(min(d4), min(d5)))
    report.append(f"Overlap: {overlap:.2f}\n")

    # ── 2. Component budget at N=20 ──
    report.append("## 2. Component Budget (4D vs 5D)\n")
    report.append("| N | ΔlogH(5D−4D) | Δ_dim_well(d*=4.1) | Δ_wall | Δ_Σ_hist | Δ_F10 | 4D<5D win% |")
    report.append("|---|-------------|-------------------|--------|----------|-------|------------|")

    for N in n_values:
        v4 = by_nf.get((N, "Lor4D"), [])
        v5 = by_nf.get((N, "Lor5D"), [])
        if not v4 or not v5:
            continue

        dlh = np.mean([r["log_H"] for r in v5]) - np.mean([r["log_H"] for r in v4])
        
        dw4 = 1.0 * N * np.mean([(r["d_eff"] - 4.1)**2 for r in v4])
        dw5 = 1.0 * N * np.mean([(r["d_eff"] - 4.1)**2 for r in v5])
        d_dim = dw5 - dw4
        
        alpha_N = 16.0 * (20.0/max(N,1))**0.5
        dwall = (np.mean([alpha_N*sigmoid((r["R"]-0.25)/0.015) for r in v5]) - 
                 np.mean([alpha_N*sigmoid((r["R"]-0.25)/0.015) for r in v4]))
        
        dsh = -10 * (np.mean([r["sigma_hist"] for r in v5]) - np.mean([r["sigma_hist"] for r in v4]))
        
        df10 = dlh + d_dim + dwall + dsh
        
        # Win rate
        n_p = min(len(v4), len(v5))
        f10_4 = [compute_F10(v4[i], d_star=4.1, gamma=1.0, lam=10.0) for i in range(n_p)]
        f10_5 = [compute_F10(v5[i], d_star=4.1, gamma=1.0, lam=10.0) for i in range(n_p)]
        wins = sum(1 for i in range(n_p) if f10_4[i] < f10_5[i])
        wr = wins / n_p
        
        report.append(f"| {N} | {dlh:+.1f} | {d_dim:+.1f} | {dwall:+.2f} | {dsh:+.2f} | {df10:+.1f} | {wr:.0%} |")

    # ── 3. Asymmetric well: steeper on high-d side ──
    report.append("\n## 3. Asymmetric Well Exploration\n")
    report.append("Use γ_hi for d_eff > d* (penalize high-d more aggressively).\n")

    report.append("| d* | γ_lo | γ_hi | N=20 4D<5D | N=36 | N=52 | N=72 | N=100 | N=20 4D<3D | min(3D<KR) |")
    report.append("|-----|------|------|------------|------|------|------|-------|------------|------------|")

    best_asym = None
    best_asym_score = -1

    for d_star in [3.95, 4.00, 4.05, 4.10, 4.15]:
        for gamma_lo in [0.5, 1.0, 1.5, 2.0]:
            for gamma_hi in [1.5, 2.0, 3.0, 5.0, 8.0]:
                per_n = {}
                for N in n_values:
                    f10_vals = {}
                    for f in all_families:
                        vals = by_nf.get((N, f), [])
                        if vals:
                            f10_vals[f] = [compute_F10(r, d_star=d_star, gamma=gamma_lo, 
                                                        lam=10.0, gamma_hi=gamma_hi) for r in vals]
                    
                    def wr(left, right):
                        lv, rv = f10_vals.get(left,[]), f10_vals.get(right,[])
                        if not lv or not rv: return 0.0
                        n_p = min(len(lv), len(rv))
                        return sum(1 for i in range(n_p) if lv[i]<rv[i]) / n_p
                    
                    per_n[N] = {
                        "4D<5D": wr("Lor4D","Lor5D"),
                        "4D<3D": wr("Lor4D","Lor3D"),
                        "4D<2D": wr("Lor4D","Lor2D"),
                        "3D<KR": wr("Lor3D","KR_like"),
                        "4D<KR": wr("Lor4D","KR_like"),
                    }

                min_45 = min(per_n[N]["4D<5D"] for N in n_values)
                min_43 = min(per_n[N]["4D<3D"] for N in n_values)
                min_42 = min(per_n[N]["4D<2D"] for N in n_values)
                min_3kr = min(per_n[N]["3D<KR"] for N in n_values)
                min_4kr = min(per_n[N]["4D<KR"] for N in n_values)

                # Score: want all A-pairs ≥ 70% and B preserved
                score = (min_45 * 40 + min_43 * 30 + min_42 * 10 
                         + min(min_3kr, min_4kr) * 20)

                if score > best_asym_score:
                    best_asym_score = score
                    best_asym = (d_star, gamma_lo, gamma_hi, per_n)

                # Show promising ones
                if min_45 >= 0.7 and min_43 >= 0.8 and min_3kr >= 0.8:
                    n20_45 = per_n[20]["4D<5D"]
                    n36_45 = per_n[36]["4D<5D"]
                    n52_45 = per_n[52]["4D<5D"]
                    n72_45 = per_n[72]["4D<5D"]
                    n100_45 = per_n[100]["4D<5D"]
                    n20_43 = per_n[20]["4D<3D"]
                    report.append(f"| {d_star} | {gamma_lo} | {gamma_hi} | {n20_45:.0%} | "
                                  f"{n36_45:.0%} | {n52_45:.0%} | {n72_45:.0%} | {n100_45:.0%} | "
                                  f"{n20_43:.0%} | {min_3kr:.0%} |")

    # ── Show best asymmetric result ──
    if best_asym:
        d_star, gamma_lo, gamma_hi, per_n = best_asym
        report.append(f"\n### Best Asymmetric: d*={d_star}, γ_lo={gamma_lo}, γ_hi={gamma_hi}\n")
        report.append("| N | 4D<2D | 4D<3D | 4D<5D | 3D<KR | 4D<KR |")
        report.append("|---|-------|-------|-------|-------|-------|")
        for N in n_values:
            pn = per_n[N]
            report.append(f"| {N} | {pn['4D<2D']:.0%} | {pn['4D<3D']:.0%} | "
                          f"{pn['4D<5D']:.0%} | {pn['3D<KR']:.0%} | {pn['4D<KR']:.0%} |")

    # ── 4. Summary ──
    report.append("\n## 4. Summary\n")
    report.append("The 4D<5D weakness at N=20 is driven by:")
    report.append("1. ΔlogH(5D−4D) ≈ +2.6 at N=20 (5D has more entropy)")
    report.append("2. d_eff(5D)≈4.4 is close to d*=4.1 → small penalty gap (0.09−0.01=0.08 per element)")
    report.append("3. At N=20, total gap = 0.08×20 = 1.6 < 2.6 (logH dominates)\n")
    report.append("Asymmetric well with γ_hi > γ_lo can partially fix this.")

    out = "\n".join(report)
    outdir = Path("outputs_d_recovery")
    md_path = outdir / "prediction_a_f10_45_diagnosis.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
