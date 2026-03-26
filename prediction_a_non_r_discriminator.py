"""
Prediction A — Non-R Dimension Discriminator Search
=====================================================
The A-B tradeoff shows that R-based O(N) terms conflate KR with Lor4D.
We need a quantity φ(d) that:
  1. φ(4D) < φ(2D/3D/5D) — V-shaped with 4D minimum
  2. φ(KR) ≠ φ(Lor4D) — separates KR from Lor4D
  3. |φ| grows O(N) or can be amplified to O(N)

Candidates from existing data:
  - d_eff (Myrheim-Meyer dimension estimator)
  - xi_dim (spectral dimension functional)
  - sigma_hist (historical deposition)
  - Combinations: d_eff × N, (d_eff - d_target)² × N

New idea: use d_eff to build a dimension-aware penalty that
grows O(N) but distinguishes KR from Lor families by their
different d_eff values.
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
    report.append("# Prediction A — Non-R Dimension Discriminator Search\n")

    # ── 1. d_eff landscape ──
    report.append("## 1. d_eff (Myrheim-Meyer) Landscape\n")
    report.append("| N | d_eff(2D) | d_eff(3D) | d_eff(4D) | d_eff(5D) | d_eff(KR) | 4D vs KR separated? |")
    report.append("|---|-----------|-----------|-----------|-----------|-----------|---------------------|")

    for N in n_values:
        deff = {}
        for f in all_families:
            vals = by_nf.get((N, f), [])
            if vals:
                deff[f] = np.mean([r["d_eff"] for r in vals])
        if len(deff) >= 5:
            sep = abs(deff["Lor4D"] - deff["KR_like"]) > 0.5
            report.append(f"| {N} | {deff['Lor2D']:.2f} | {deff['Lor3D']:.2f} | "
                          f"{deff['Lor4D']:.2f} | {deff['Lor5D']:.2f} | {deff['KR_like']:.2f} | "
                          f"{'✅' if sep else '❌'} (Δ={deff['Lor4D']-deff['KR_like']:+.2f}) |")

    # ── 2. xi_dim landscape ──
    report.append("\n## 2. xi_dim Landscape\n")
    report.append("| N | ξ(2D) | ξ(3D) | ξ(4D) | ξ(5D) | ξ(KR) | 4D vs KR? |")
    report.append("|---|-------|-------|-------|-------|-------|-----------|")

    for N in n_values:
        xis = {}
        for f in all_families:
            vals = by_nf.get((N, f), [])
            if vals:
                xis[f] = np.mean([r["xi_dim"] for r in vals])
        if len(xis) >= 5:
            report.append(f"| {N} | {xis['Lor2D']:.3f} | {xis['Lor3D']:.3f} | "
                          f"{xis['Lor4D']:.3f} | {xis['Lor5D']:.3f} | {xis['KR_like']:.3f} | "
                          f"Δ={xis['Lor4D']-xis['KR_like']:+.3f} |")

    # ── 3. d_eff-based O(N) discriminators ──
    report.append("\n## 3. d_eff-Based O(N) Discriminator Candidates\n")

    # C1: N·(d_eff - d*)² — quadratic well in d_eff space
    report.append("### C1: logH + γ·N·(d_eff − d*)²\n")
    report.append("d* = target effective dimension. If d*(4D)≈4.0, this penalizes non-4D.\n")
    report.append("Key question: does KR have d_eff ≠ 4D?\n")

    for d_star in [2.5, 3.0, 3.5, 4.0]:
        for gamma in [0.5, 1.0, 2.0, 5.0]:
            ranks_all_n = []
            for N in n_values:
                means = {}
                for f in lor_families:
                    vals = by_nf.get((N, f), [])
                    if vals:
                        means[f] = np.mean([r["log_H"] + gamma * N * (r["d_eff"] - d_star) ** 2
                                            for r in vals])
                if len(means) == 4:
                    rank_sorted = sorted(means.items(), key=lambda x: x[1])
                    rank_4d = next(i + 1 for i, (f, _) in enumerate(rank_sorted) if f == "Lor4D")
                    ranks_all_n.append(rank_4d)

            if any(r == 1 for r in ranks_all_n):
                r_str = "/".join(str(r) for r in ranks_all_n)
                n_min = sum(1 for r in ranks_all_n if r == 1)
                report.append(f"  d*={d_star}, γ={gamma}: ranks={r_str}, 4D=min: {n_min}/5")

    # ── 4. Hybrid: F7 + d_eff-well (preserving Prediction B) ──
    report.append("\n## 4. F10 Hybrid: F7 + δ·N·(d_eff − d*)²\n")
    report.append("Add a d_eff-based O(N) term to F7 (which already satisfies B).\n")
    report.append("Key: d_eff(KR)≈2.7 ≠ d_eff(4D)≈4.0, so this should NOT destroy B!\n")

    best_f10 = None
    best_f10_score = -1

    for d_star in np.arange(2.0, 5.0, 0.25):
        for delta in np.arange(0.1, 3.0, 0.1):
            # Check Prediction A (dimension ordering)
            a_ranks = []
            for N in n_values:
                means = {}
                for f in lor_families:
                    vals = by_nf.get((N, f), [])
                    if vals:
                        means[f] = np.mean([
                            r["F7"] + delta * N * (r["d_eff"] - d_star) ** 2
                            for r in vals
                        ])
                if len(means) == 4:
                    rank_sorted = sorted(means.items(), key=lambda x: x[1])
                    rank_4d = next(i + 1 for i, (f, _) in enumerate(rank_sorted) if f == "Lor4D")
                    a_ranks.append(rank_4d)

            # Check Prediction B (Lor2D/3D < KR)
            b_ok = True
            for N in n_values:
                for lor in ["Lor2D", "Lor3D"]:
                    lv = by_nf.get((N, lor), [])
                    kv = by_nf.get((N, "KR_like"), [])
                    if not lv or not kv:
                        continue
                    f10_lor = [r["F7"] + delta * N * (r["d_eff"] - d_star) ** 2 for r in lv]
                    f10_kr = [r["F7"] + delta * N * (r["d_eff"] - d_star) ** 2 for r in kv]
                    n_p = min(len(f10_lor), len(f10_kr))
                    wins = sum(1 for i in range(n_p) if f10_lor[i] < f10_kr[i])
                    if wins / n_p < 0.5:
                        b_ok = False
                        break
                if not b_ok:
                    break

            a_score = sum(1 for r in a_ranks if r == 1)
            score = a_score * 10 + (50 if b_ok else 0)

            if score > best_f10_score:
                best_f10_score = score
                best_f10 = (d_star, delta, a_ranks, b_ok)

    if best_f10:
        d_star, delta, a_ranks, b_ok = best_f10
        report.append(f"\n**Best F10**: F7 + {delta:.1f}·N·(d_eff − {d_star:.2f})²")
        report.append(f"A ranks: {[f'{r}/4' for r in a_ranks]}, 4D=min: {sum(1 for r in a_ranks if r==1)}/5")
        report.append(f"B preserved (Lor2D/3D < KR at all N): {'✅' if b_ok else '❌'}\n")

        # Detailed per-N analysis
        report.append("### F10 Per-N Details\n")
        report.append("| N | F10(2D) | F10(3D) | F10(4D) | F10(5D) | F10(KR) | ordering (Lor only) | 4D_rank |")
        report.append("|---|---------|---------|---------|---------|---------|---------------------|---------|")

        for N in n_values:
            means = {}
            for f in all_families:
                vals = by_nf.get((N, f), [])
                if vals:
                    means[f] = np.mean([r["F7"] + delta * N * (r["d_eff"] - d_star) ** 2 for r in vals])
            if len(means) >= 5:
                lor_means = {f: means[f] for f in lor_families}
                rank_sorted = sorted(lor_means.items(), key=lambda x: x[1])
                rank_4d = next(i + 1 for i, (f, _) in enumerate(rank_sorted) if f == "Lor4D")
                order_str = " < ".join(f.replace("Lor", "") for f, _ in rank_sorted)
                report.append(f"| {N} | {means.get('Lor2D',0):.1f} | {means.get('Lor3D',0):.1f} | "
                              f"{means.get('Lor4D',0):.1f} | {means.get('Lor5D',0):.1f} | "
                              f"{means.get('KR_like',0):.1f} | {order_str} | {rank_4d}/4 |")

        # Per-N pairwise: A pairs + B pairs
        report.append("\n### F10 Pairwise Win Rates\n")

        all_pairs = [
            ("Lor4D", "Lor2D", "A: 4D<2D"),
            ("Lor4D", "Lor3D", "A: 4D<3D"),
            ("Lor4D", "Lor5D", "A: 4D<5D"),
            ("Lor2D", "KR_like", "B: 2D<KR"),
            ("Lor3D", "KR_like", "B: 3D<KR"),
            ("Lor4D", "KR_like", "B: 4D<KR"),
        ]

        for left, right, label in all_pairs:
            report.append(f"\n**{label}**")
            for N in n_values:
                lv = by_nf.get((N, left), [])
                rv = by_nf.get((N, right), [])
                if not lv or not rv:
                    continue
                n_p = min(len(lv), len(rv))
                f10_l = [lv[i]["F7"] + delta * N * (lv[i]["d_eff"] - d_star) ** 2 for i in range(n_p)]
                f10_r = [rv[i]["F7"] + delta * N * (rv[i]["d_eff"] - d_star) ** 2 for i in range(n_p)]
                wins = sum(1 for i in range(n_p) if f10_l[i] < f10_r[i])
                pct = wins / n_p
                u, p = sp_stats.mannwhitneyu(f10_l, f10_r, alternative="less")
                sig = "★★★" if p < 0.001 else "★★" if p < 0.01 else "★" if p < 0.05 else ""
                report.append(f"  N={N}: {pct:.0%} {sig}")

    out = "\n".join(report)
    outdir = Path("outputs_d_recovery")
    md_path = outdir / "prediction_a_non_r_discriminator.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
