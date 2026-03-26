"""
Prediction A — O(N) Wall Exploration
======================================
Goal: Find a term that grows O(N) and penalizes low-d MORE than 4D,
thereby restoring F(4D) as the global minimum.

Key insight from the data:
  - logH: 2D < 3D < 4D < 5D (monotone increasing with d) — drives wrong ordering
  - wall: 2D ≈ 3D >> 4D > 5D at small N, but ALL → 0 at large N
  - xi_dim: 2D ≈ 3D ≈ 0, 4D ≈ 0–0.7, 5D ≈ 0–0.1 (not enough)

What we need: a quantity Q(d) such that:
  Q(2D) > Q(3D) > Q(4D) < Q(5D)  (V-shaped with minimum at 4D)
  AND |Q| grows O(N) to compete with logH

Candidates to explore:
  C1: β·N·f(R)  — N-linear wall using R (occupancy)
  C2: logH/N·g(R)  — entropy density × admissibility
  C3: N·(R - R_target(d))²  — quadratic penalty around d=4 optimal R
  C4: logH - μ·N·R  — logH corrected by density (R absorbs dimension info)
  C5: logH·(1 - R)  — entropy × sparsity (low R = high d, high logH)
  C6: Direct dimension proxy via R — R is a monotone function of d
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats


def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def load_csv(path: str) -> list[dict]:
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


# ── Step 1: Component landscape analysis ──

def analyze_landscape(rows, n_values, lor_families):
    by_nf = defaultdict(list)
    for r in rows:
        by_nf[(r["N"], r["family"])].append(r)

    report = []
    report.append("# Prediction A — O(N) Wall Exploration\n")

    # R(d) scaling
    report.append("## 1. R(d) Scaling — Does R encode dimension?\n")
    report.append("| N | R(2D) | R(3D) | R(4D) | R(5D) | R monotone? |")
    report.append("|---|-------|-------|-------|-------|-------------|")

    for N in n_values:
        Rs = {}
        for f in lor_families:
            vals = by_nf.get((N, f), [])
            if vals:
                Rs[f] = np.mean([r["R"] for r in vals])
        if len(Rs) == 4:
            mono = Rs["Lor2D"] > Rs["Lor3D"] > Rs["Lor4D"] > Rs["Lor5D"]
            report.append(f"| {N} | {Rs['Lor2D']:.3f} | {Rs['Lor3D']:.3f} | "
                          f"{Rs['Lor4D']:.3f} | {Rs['Lor5D']:.3f} | {'✅ ↓' if mono else '❌'} |")

    # logH/N scaling (entropy density)
    report.append("\n## 2. Entropy Density logH/N\n")
    report.append("| N | h(2D) | h(3D) | h(4D) | h(5D) | h monotone? |")
    report.append("|---|-------|-------|-------|-------|-------------|")

    for N in n_values:
        hs = {}
        for f in lor_families:
            vals = by_nf.get((N, f), [])
            if vals:
                hs[f] = np.mean([r["log_H"] for r in vals]) / N
        if len(hs) == 4:
            mono = hs["Lor2D"] < hs["Lor3D"] < hs["Lor4D"] < hs["Lor5D"]
            report.append(f"| {N} | {hs['Lor2D']:.3f} | {hs['Lor3D']:.3f} | "
                          f"{hs['Lor4D']:.3f} | {hs['Lor5D']:.3f} | {'✅ ↑' if mono else '❌'} |")

    # N·R scaling (does N·R grow with dimension?)
    report.append("\n## 3. N·R Product\n")
    report.append("| N | NR(2D) | NR(3D) | NR(4D) | NR(5D) |")
    report.append("|---|--------|--------|--------|--------|")

    for N in n_values:
        for f in lor_families:
            pass
        nrs = {}
        for f in lor_families:
            vals = by_nf.get((N, f), [])
            if vals:
                nrs[f] = N * np.mean([r["R"] for r in vals])
        if len(nrs) == 4:
            report.append(f"| {N} | {nrs['Lor2D']:.1f} | {nrs['Lor3D']:.1f} | "
                          f"{nrs['Lor4D']:.1f} | {nrs['Lor5D']:.1f} |")

    # logH - μ·N·R for various μ
    report.append("\n## 4. Candidate C4: logH − μ·N·R\n")
    report.append("If R encodes dimension, then N·R correction could flatten logH's dimension dependence.\n")

    for mu in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
        report.append(f"\n### μ = {mu}\n")
        report.append("| N | C4(2D) | C4(3D) | C4(4D) | C4(5D) | 4D_rank | 4D=min? |")
        report.append("|---|--------|--------|--------|--------|---------|---------|")

        for N in n_values:
            c4s = {}
            for f in lor_families:
                vals = by_nf.get((N, f), [])
                if vals:
                    c4s[f] = np.mean([r["log_H"] - mu * N * r["R"] for r in vals])
            if len(c4s) == 4:
                rank_sorted = sorted(c4s.items(), key=lambda x: x[1])
                rank_4d = next(i+1 for i, (f, _) in enumerate(rank_sorted) if f == "Lor4D")
                is_min = rank_sorted[0][0] == "Lor4D"
                report.append(f"| {N} | {c4s['Lor2D']:.1f} | {c4s['Lor3D']:.1f} | "
                              f"{c4s['Lor4D']:.1f} | {c4s['Lor5D']:.1f} | {rank_4d}/4 | "
                              f"{'✅' if is_min else '❌'} |")

    # Candidate C5: logH·(1-R) — entropy × sparsity
    report.append("\n## 5. Candidate C5: logH·(1 − R)\n")
    report.append("Penalizes high entropy + high connectivity (low d = high R → low penalty)\n")
    report.append("| N | C5(2D) | C5(3D) | C5(4D) | C5(5D) | 4D_rank | 4D=min? |")
    report.append("|---|--------|--------|--------|--------|---------|---------|")

    for N in n_values:
        c5s = {}
        for f in lor_families:
            vals = by_nf.get((N, f), [])
            if vals:
                c5s[f] = np.mean([r["log_H"] * (1 - r["R"]) for r in vals])
        if len(c5s) == 4:
            rank_sorted = sorted(c5s.items(), key=lambda x: x[1])
            rank_4d = next(i+1 for i, (f, _) in enumerate(rank_sorted) if f == "Lor4D")
            is_min = rank_sorted[0][0] == "Lor4D"
            report.append(f"| {N} | {c5s['Lor2D']:.1f} | {c5s['Lor3D']:.1f} | "
                          f"{c5s['Lor4D']:.1f} | {c5s['Lor5D']:.1f} | {rank_4d}/4 | "
                          f"{'✅' if is_min else '❌'} |")

    # Candidate C3: N·(R - R*)² quadratic well
    report.append("\n## 6. Candidate C3: logH + γ·N·(R − R*)²\n")
    report.append("Quadratic penalty around optimal R=R*. If R*(4D) ≈ 0.25, this might create a V-shaped well.\n")

    for R_star in [0.15, 0.20, 0.25, 0.30, 0.35]:
        for gamma in [5.0, 10.0, 20.0]:
            # Quick check: what's the 4D rank at N=100?
            c3s = {}
            for f in lor_families:
                vals = by_nf.get((100, f), [])
                if vals:
                    c3s[f] = np.mean([r["log_H"] + gamma * N * (r["R"] - R_star)**2 for r in vals])
            if len(c3s) == 4:
                rank_sorted = sorted(c3s.items(), key=lambda x: x[1])
                is_min = rank_sorted[0][0] == "Lor4D"
                if is_min:
                    report.append(f"  ✅ R*={R_star}, γ={gamma}: 4D=min at N=100! "
                                  f"Order: {' < '.join(f.replace('Lor','') for f,_ in rank_sorted)}")

    # Full scan of C3 for promising parameters
    report.append("\n### Full C3 scan (best parameters)\n")
    report.append("| R* | γ | N=20 4D_rank | N=36 | N=52 | N=72 | N=100 | all_min? |")
    report.append("|-----|---|-------------|------|------|------|-------|----------|")

    for R_star in [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
        for gamma in [2.0, 5.0, 10.0, 20.0, 50.0]:
            ranks = []
            for N in n_values:
                c3s = {}
                for f in lor_families:
                    vals = by_nf.get((N, f), [])
                    if vals:
                        c3s[f] = np.mean([r["log_H"] + gamma * N * (r["R"] - R_star)**2 for r in vals])
                if len(c3s) == 4:
                    rank_sorted = sorted(c3s.items(), key=lambda x: x[1])
                    rank_4d = next(i+1 for i, (f, _) in enumerate(rank_sorted) if f == "Lor4D")
                    ranks.append(rank_4d)
                else:
                    ranks.append(0)

            all_min = all(r == 1 for r in ranks)
            if any(r == 1 for r in ranks):  # only show if 4D is min at least once
                report.append(f"| {R_star} | {gamma} | {ranks[0]}/4 | {ranks[1]}/4 | "
                              f"{ranks[2]}/4 | {ranks[3]}/4 | {ranks[4]}/4 | "
                              f"{'✅✅' if all_min else '❌'} |")

    # Candidate C6: logH − μ·N·R + quadratic correction
    report.append("\n## 7. Candidate C6 (hybrid): logH − μ·N·R + γ·N·(R−R*)²\n")

    best_combo = None
    best_score = -1

    for mu in np.arange(0.5, 4.0, 0.5):
        for R_star in np.arange(0.10, 0.45, 0.05):
            for gamma in [2.0, 5.0, 10.0, 20.0]:
                ranks = []
                for N in n_values:
                    c6s = {}
                    for f in lor_families:
                        vals = by_nf.get((N, f), [])
                        if vals:
                            c6s[f] = np.mean([
                                r["log_H"] - mu * N * r["R"] + gamma * N * (r["R"] - R_star)**2
                                for r in vals
                            ])
                    if len(c6s) == 4:
                        rank_sorted = sorted(c6s.items(), key=lambda x: x[1])
                        rank_4d = next(i+1 for i, (f, _) in enumerate(rank_sorted) if f == "Lor4D")
                        ranks.append(rank_4d)

                if len(ranks) == 5:
                    score = sum(1 for r in ranks if r == 1)
                    if score > best_score:
                        best_score = score
                        best_combo = (mu, R_star, gamma, ranks)

    if best_combo:
        mu, Rs, gam, ranks = best_combo
        report.append(f"\n**Best combo**: μ={mu:.1f}, R*={Rs:.2f}, γ={gam:.0f}")
        report.append(f"4D ranks: {[f'{r}/4' for r in ranks]}")
        report.append(f"4D=min count: {sum(1 for r in ranks if r == 1)}/5\n")

        # Show the full ordering for this combo
        report.append("| N | F9(2D) | F9(3D) | F9(4D) | F9(5D) | ordering |")
        report.append("|---|--------|--------|--------|--------|----------|")

        for N in n_values:
            c6s = {}
            for f in lor_families:
                vals = by_nf.get((N, f), [])
                if vals:
                    c6s[f] = np.mean([
                        r["log_H"] - mu * N * r["R"] + gam * N * (r["R"] - Rs)**2
                        for r in vals
                    ])
            if len(c6s) == 4:
                rank_sorted = sorted(c6s.items(), key=lambda x: x[1])
                order_str = " < ".join(f.replace("Lor", "") for f, _ in rank_sorted)
                report.append(f"| {N} | {c6s['Lor2D']:.1f} | {c6s['Lor3D']:.1f} | "
                              f"{c6s['Lor4D']:.1f} | {c6s['Lor5D']:.1f} | {order_str} |")

    # ── Pairwise win rates for best F9 ──
    if best_combo:
        mu, Rs, gam, _ = best_combo
        report.append(f"\n## 8. F9 Pairwise Win Rates (μ={mu}, R*={Rs:.2f}, γ={gam})\n")

        key_pairs = [("Lor4D", "Lor2D"), ("Lor4D", "Lor3D"), ("Lor4D", "Lor5D")]
        pair_labels = ["4D<2D", "4D<3D", "4D<5D"]

        report.append("| N | " + " | ".join(pair_labels) + " |")
        report.append("|---|" + "|".join(["------"] * 3) + "|")

        for N in n_values:
            cells = []
            for left, right in key_pairs:
                lv = by_nf.get((N, left), [])
                rv = by_nf.get((N, right), [])
                if not lv or not rv:
                    cells.append("—")
                    continue
                n_p = min(len(lv), len(rv))
                wins = 0
                for i in range(n_p):
                    fl = lv[i]["log_H"] - mu * N * lv[i]["R"] + gam * N * (lv[i]["R"] - Rs)**2
                    fr = rv[i]["log_H"] - mu * N * rv[i]["R"] + gam * N * (rv[i]["R"] - Rs)**2
                    if fl < fr:
                        wins += 1
                pct = wins / n_p
                sig = "★★★" if pct >= 0.8 else "★" if pct >= 0.6 else ""
                cells.append(f"{pct:.0%}{sig}")
            report.append(f"| {N} | {' | '.join(cells)} |")

    return "\n".join(report)


def main():
    csv_path = "outputs_d_recovery/prediction_c_f7_large_n.csv"
    print(f"Loading: {csv_path}")
    rows = load_csv(csv_path)
    n_values = sorted(set(r["N"] for r in rows))
    lor_families = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]
    print(f"Loaded {len(rows)} rows")

    report = analyze_landscape(rows, n_values, lor_families)

    outdir = Path("outputs_d_recovery")
    md_path = outdir / "prediction_a_on_wall_exploration.md"
    md_path.write_text(report, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + report)


if __name__ == "__main__":
    main()
