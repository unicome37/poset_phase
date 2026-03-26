"""
Prediction A — F7 Dimension Selection at Large N
=================================================
Core question: Does the definitive F7 select 4D as the optimal dimension
across N = 20–100?

Prediction A states: 3+1 dimensional spacetime is a statistical stability window.
This means F7(Lor4D) should be LOWER than F7(Lor2D), F7(Lor3D), and F7(Lor5D).

Two mechanisms:
  Mechanism I (low-d rejection): R-wall excludes 2D, H_int energy flip excludes 3D
  Mechanism II (high-d barrier): Ξ_d + logH convexity creates 4D→5D barrier

Design:
  - N ∈ {20, 36, 52, 72, 100}
  - 4 families: Lor2D, Lor3D, Lor4D, Lor5D  (no KR — that's Pred B)
  - 10 reps per (family, N) = 200 posets total
  - Per poset: compute F7 and all components
  - Analysis:
    Q1. At each N, does Lor4D have the lowest mean F7?
    Q2. Pairwise win rates: Lor4D vs each other dimension
    Q3. F7 component decomposition: which terms drive 4D selection?
    Q4. N-scaling: does the 4D advantage grow or shrink with N?
    Q5. Margin of victory and statistical significance

Usage:
    python prediction_a_f7_dimension_selection.py [--reps 10] [--quick]
"""
from __future__ import annotations

import argparse
import csv
import math
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)
from unified_functional import (
    compute_log_H,
    compute_pi_geo,
    compute_sigma_hist,
    compute_xi_dim,
)
from prediction_a_bd_bridge import count_intervals_fast


# ══════════════════════════════════════════════════════════════════════════
# F7 computation (definitive §5.10.7 model)
# ══════════════════════════════════════════════════════════════════════════

def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_R(poset: Poset) -> float:
    counts = count_intervals_fast(poset)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    C0 = counts.get(0, 0)
    return 1.0 - C0 / total


def adaptive_sis_runs(N: int) -> int:
    if N <= 36:
        return 512
    elif N <= 52:
        return 256
    elif N <= 72:
        return 128
    else:
        return 64


def compute_F7(poset: Poset,
               alpha0: float = 16.0, q: float = -0.5,
               lam: float = 10.0, eta: float = 0.6,
               Rc: float = 0.25, w: float = 0.015,
               N0: float = 20.0) -> dict:
    N = poset.n
    sis = adaptive_sis_runs(N)
    log_H = compute_log_H(poset, n_runs=sis)
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim_val, d_eff = compute_xi_dim(poset)
    R = compute_R(poset)

    alpha_N = alpha0 * (N0 / max(N, 1)) ** abs(q)
    wall = alpha_N * sigmoid((R - Rc) / w)
    F7 = log_H + 0.0004 * pi_geo - lam * sigma_hist + eta * xi_dim_val + wall

    return {
        "F7": F7,
        "log_H": log_H,
        "pi_geo": pi_geo,
        "sigma_hist": sigma_hist,
        "xi_dim": xi_dim_val,
        "d_eff": d_eff,
        "R": R,
        "wall": wall,
        "alpha_N": alpha_N,
        "term_logH": log_H,
        "term_pi": 0.0004 * pi_geo,
        "term_sigma": -lam * sigma_hist,
        "term_xi": eta * xi_dim_val,
        "term_wall": wall,
    }


FAMILIES = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
}


# ══════════════════════════════════════════════════════════════════════════
# Data collection
# ══════════════════════════════════════════════════════════════════════════

def collect_data(n_values, reps, seed_base=42):
    rows = []
    total = len(n_values) * len(FAMILIES) * reps
    count = 0
    t0 = time.time()

    for N in n_values:
        for fam_name, gen in FAMILIES.items():
            for rep in range(reps):
                count += 1
                seed = seed_base + N * 1000 + rep
                poset = gen(N, seed=seed)
                comp = compute_F7(poset)
                row = {
                    "family": fam_name,
                    "N": N,
                    "rep": rep,
                    "seed": seed,
                    **comp,
                }
                rows.append(row)

                if count % 20 == 0 or count == total:
                    elapsed = time.time() - t0
                    eta_s = (elapsed / count) * (total - count)
                    print(f"  [{count}/{total}] {fam_name} N={N} rep={rep}: "
                          f"F7={comp['F7']:.2f}  (elapsed {elapsed:.0f}s, ETA {eta_s:.0f}s)")

    return rows


# ══════════════════════════════════════════════════════════════════════════
# Analysis
# ══════════════════════════════════════════════════════════════════════════

def generate_report(rows, n_values):
    by_nf = defaultdict(list)
    for r in rows:
        by_nf[(r["N"], r["family"])].append(r)

    families = sorted(set(r["family"] for r in rows))
    report = []
    report.append("# Prediction A — F7 Dimension Selection at Large N\n")
    report.append(f"**Data**: {len(rows)} samples, N ∈ {{{', '.join(str(n) for n in n_values)}}}\n")
    report.append(f"**F7 model**: §5.10.7 definitive (α₀=16, q=−0.5, λ=10, η=0.6, Rc=0.25, w=0.015)\n")

    # ── Q1: Mean F7 per family per N ──
    report.append("\n## Q1: Mean F7 by Family and N\n")
    report.append("| N | Lor2D | Lor3D | Lor4D | Lor5D | Winner | 4D_rank |")
    report.append("|---|-------|-------|-------|-------|--------|---------|")

    for N in n_values:
        means = {}
        for fam in families:
            vals = [r["F7"] for r in by_nf.get((N, fam), [])]
            means[fam] = np.mean(vals) if vals else float('inf')

        winner = min(means, key=means.get)
        sorted_fams = sorted(families, key=lambda f: means[f])
        rank_4d = sorted_fams.index("Lor4D") + 1

        report.append(f"| {N} | {means.get('Lor2D', 0):.2f} | {means.get('Lor3D', 0):.2f} | "
                      f"{means.get('Lor4D', 0):.2f} | {means.get('Lor5D', 0):.2f} | "
                      f"**{winner}** | {rank_4d} |")

    # ── Q2: Pairwise win rates ──
    report.append("\n## Q2: Pairwise Win Rates (Lor4D vs Others)\n")
    comparisons = [("Lor4D", "Lor2D"), ("Lor4D", "Lor3D"), ("Lor4D", "Lor5D")]

    for left, right in comparisons:
        report.append(f"\n### {left} vs {right}\n")
        report.append("| N | wins | total | win% | mean_ΔF7 | p_value | sig |")
        report.append("|---|------|-------|------|----------|---------|-----|")

        total_wins = 0
        total_pairs = 0
        for N in n_values:
            l_rows = by_nf.get((N, left), [])
            r_rows = by_nf.get((N, right), [])
            n_pairs = min(len(l_rows), len(r_rows))

            wins = 0
            deltas = []
            for i in range(n_pairs):
                d = l_rows[i]["F7"] - r_rows[i]["F7"]
                deltas.append(d)
                if d < 0:
                    wins += 1

            total_wins += wins
            total_pairs += n_pairs
            pct = 100 * wins / n_pairs if n_pairs > 0 else 0
            mean_d = np.mean(deltas) if deltas else 0

            # Wilcoxon signed-rank test
            if len(deltas) >= 5:
                try:
                    stat, p_val = sp_stats.wilcoxon(deltas, alternative='less')
                except ValueError:
                    p_val = 1.0
            else:
                p_val = float('nan')

            sig = "★★★" if p_val < 0.001 else "★★" if p_val < 0.01 else "★" if p_val < 0.05 else ""
            report.append(f"| {N} | {wins} | {n_pairs} | {pct:.0f}% | {mean_d:+.2f} | "
                          f"{p_val:.4f} | {sig} |")

        overall_pct = 100 * total_wins / total_pairs if total_pairs > 0 else 0
        report.append(f"\n**Overall**: {total_wins}/{total_pairs} ({overall_pct:.1f}%)")

    # ── Q3: Component decomposition ──
    report.append("\n\n## Q3: F7 Component Decomposition (Mean per Family per N)\n")
    report.append("| N | family | logH | −λΣ_hist | ηΞ_d | wall | F7 |")
    report.append("|---|--------|------|----------|------|------|----|")

    for N in n_values:
        for fam in families:
            vals = by_nf.get((N, fam), [])
            if not vals:
                continue
            mH = np.mean([r["term_logH"] for r in vals])
            mS = np.mean([r["term_sigma"] for r in vals])
            mX = np.mean([r["term_xi"] for r in vals])
            mW = np.mean([r["term_wall"] for r in vals])
            mF = np.mean([r["F7"] for r in vals])
            report.append(f"| {N} | {fam} | {mH:.2f} | {mS:+.2f} | {mX:+.2f} | {mW:.2f} | {mF:.2f} |")

    # ── Q4: N-scaling of 4D advantage ──
    report.append("\n## Q4: N-Scaling of 4D Margin\n")
    report.append("Margin = mean F7(competitor) − mean F7(Lor4D)\n")
    report.append("| N | Δ(2D−4D) | Δ(3D−4D) | Δ(5D−4D) | min_margin |")
    report.append("|---|----------|----------|----------|------------|")

    margins_vs_2d = []
    margins_vs_3d = []
    margins_vs_5d = []

    for N in n_values:
        m4 = np.mean([r["F7"] for r in by_nf.get((N, "Lor4D"), [])])
        m2 = np.mean([r["F7"] for r in by_nf.get((N, "Lor2D"), [])])
        m3 = np.mean([r["F7"] for r in by_nf.get((N, "Lor3D"), [])])
        m5 = np.mean([r["F7"] for r in by_nf.get((N, "Lor5D"), [])])

        d2 = m2 - m4
        d3 = m3 - m4
        d5 = m5 - m4
        margins_vs_2d.append(d2)
        margins_vs_3d.append(d3)
        margins_vs_5d.append(d5)
        min_m = min(d2, d3, d5)
        report.append(f"| {N} | {d2:+.2f} | {d3:+.2f} | {d5:+.2f} | {min_m:+.2f} |")

    # Trend analysis
    for label, margins in [("vs 2D", margins_vs_2d), ("vs 3D", margins_vs_3d), ("vs 5D", margins_vs_5d)]:
        rho, p = sp_stats.spearmanr(n_values, margins)
        trend = "growing" if rho > 0.5 else "shrinking" if rho < -0.5 else "stable"
        report.append(f"- **{label}** N-trend: ρ={rho:+.3f} (p={p:.3f}) — {trend}")

    # ── Q5: Overall verdict ──
    report.append("\n## Q5: Verdict\n")

    # Check if 4D wins at every N
    all_winner = True
    for N in n_values:
        means = {fam: np.mean([r["F7"] for r in by_nf.get((N, fam), [])]) for fam in families}
        winner = min(means, key=means.get)
        if winner != "Lor4D":
            all_winner = False

    # Check minimum margin
    min_margins = [min(margins_vs_2d[i], margins_vs_3d[i], margins_vs_5d[i]) for i in range(len(n_values))]
    all_positive = all(m > 0 for m in min_margins)

    if all_winner and all_positive:
        report.append("**STRONG**: Lor4D has lowest F7 at ALL N values with positive margins.")
    elif all_winner:
        report.append("**PASS**: Lor4D wins at all N, but some margins are negative (individual rep level).")
    else:
        failing_ns = [n_values[i] for i in range(len(n_values)) if min_margins[i] <= 0]
        report.append(f"**PARTIAL**: Lor4D does NOT win at N={failing_ns}.")

    # Component attribution
    report.append("\n### Component Attribution\n")
    report.append("Which F7 terms favor 4D?\n")
    for N in n_values:
        v4 = by_nf.get((N, "Lor4D"), [])
        v5 = by_nf.get((N, "Lor5D"), [])
        if not v4 or not v5:
            continue
        for term in ["term_logH", "term_sigma", "term_xi", "term_wall"]:
            m4t = np.mean([r[term] for r in v4])
            m5t = np.mean([r[term] for r in v5])
            delta = m4t - m5t
            favor = "4D" if delta < 0 else "5D"
            report.append(f"- N={N}: {term} → Δ(4D−5D) = {delta:+.2f} (favors {favor})")
        report.append("")

    # Mechanism I check (low-d rejection)
    report.append("### Mechanism I: Low-Dimension Rejection\n")
    for N in n_values:
        v2 = by_nf.get((N, "Lor2D"), [])
        v4 = by_nf.get((N, "Lor4D"), [])
        if v2 and v4:
            wall_2d = np.mean([r["term_wall"] for r in v2])
            wall_4d = np.mean([r["term_wall"] for r in v4])
            R_2d = np.mean([r["R"] for r in v2])
            R_4d = np.mean([r["R"] for r in v4])
            report.append(f"- N={N}: wall(2D)={wall_2d:.2f} vs wall(4D)={wall_4d:.2f}, "
                          f"R(2D)={R_2d:.3f} vs R(4D)={R_4d:.3f}")

    # Mechanism II check (high-d barrier)
    report.append("\n### Mechanism II: High-Dimension Barrier\n")
    for N in n_values:
        v4 = by_nf.get((N, "Lor4D"), [])
        v5 = by_nf.get((N, "Lor5D"), [])
        if v4 and v5:
            logH_4 = np.mean([r["term_logH"] for r in v4])
            logH_5 = np.mean([r["term_logH"] for r in v5])
            xi_4 = np.mean([r["term_xi"] for r in v4])
            xi_5 = np.mean([r["term_xi"] for r in v5])
            sig_4 = np.mean([r["term_sigma"] for r in v4])
            sig_5 = np.mean([r["term_sigma"] for r in v5])
            report.append(f"- N={N}: ΔlogH(5D−4D)={logH_5-logH_4:+.2f}, "
                          f"ΔΞ(5D−4D)={xi_5-xi_4:+.2f}, "
                          f"ΔΣ(5D−4D)={sig_5-sig_4:+.2f}")

    return "\n".join(report)


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=10)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--ns", type=int, nargs="+", default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--suffix", type=str, default="")
    args = parser.parse_args()

    if args.quick:
        n_values = [20, 36, 52]
        reps = min(args.reps, 5)
    else:
        n_values = args.ns or [20, 36, 52, 72, 100]
        reps = args.reps

    print(f"=== Prediction A — F7 Dimension Selection ===")
    print(f"N values: {n_values}")
    print(f"Families: {list(FAMILIES.keys())}")
    print(f"Reps: {reps}")
    print()

    t0 = time.time()
    rows = collect_data(n_values, reps, seed_base=args.seed)
    elapsed = time.time() - t0
    print(f"\nData collection: {elapsed:.0f}s ({len(rows)} samples)")

    # Save CSV
    outdir = Path("outputs_d_recovery")
    outdir.mkdir(exist_ok=True)
    suffix = args.suffix or ""

    csv_path = outdir / f"prediction_a_f7_dimension{suffix}.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {csv_path}")

    # Generate report
    report = generate_report(rows, n_values)
    md_path = outdir / f"prediction_a_f7_dimension{suffix}.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"Saved: {md_path}")
    print("\n" + report)


if __name__ == "__main__":
    main()
