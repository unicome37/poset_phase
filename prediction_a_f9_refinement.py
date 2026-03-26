"""
Prediction A — F9 Refinement: Fine-grained parameter search
=============================================================
The C6 exploration found that logH − μ·N·R + γ·N·(R−R*)² can make 4D the
global minimum at large N. Now refine parameters to maximize coverage.

Also test adding back Σ_hist and ξ_d terms from F7.
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


def compute_f9(row, mu, R_star, gamma, lam=10.0, eta=0.6):
    """F9 = logH − μ·N·R + γ·N·(R−R*)² − λ·Σ_hist + η·ξ_d"""
    N = row["N"]
    R = row["R"]
    return (row["log_H"]
            - mu * N * R
            + gamma * N * (R - R_star) ** 2
            - lam * row["sigma_hist"]
            + eta * row["xi_dim"])


def evaluate_params(by_nf, n_values, lor_families, mu, R_star, gamma, lam=10.0, eta=0.6):
    """Evaluate a parameter set. Returns per-N 4D rank and pairwise win rates."""
    key_pairs = [("Lor4D", "Lor2D"), ("Lor4D", "Lor3D"), ("Lor4D", "Lor5D")]
    ranks = []
    win_rates = {p: [] for p in key_pairs}

    for N in n_values:
        means = {}
        f9_vals = {}
        for f in lor_families:
            vals = by_nf.get((N, f), [])
            if vals:
                f9s = [compute_f9(r, mu, R_star, gamma, lam, eta) for r in vals]
                means[f] = np.mean(f9s)
                f9_vals[f] = f9s

        if len(means) == 4:
            rank_sorted = sorted(means.items(), key=lambda x: x[1])
            rank_4d = next(i + 1 for i, (f, _) in enumerate(rank_sorted) if f == "Lor4D")
            ranks.append(rank_4d)

        for left, right in key_pairs:
            lv = f9_vals.get(left, [])
            rv = f9_vals.get(right, [])
            if lv and rv:
                n_p = min(len(lv), len(rv))
                wins = sum(1 for i in range(n_p) if lv[i] < rv[i])
                win_rates[(left, right)].append(wins / n_p)

    return ranks, win_rates


def main():
    csv_path = "outputs_d_recovery/prediction_c_f7_large_n.csv"
    rows = load_csv(csv_path)
    n_values = sorted(set(r["N"] for r in rows))
    lor_families = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]

    by_nf = defaultdict(list)
    for r in rows:
        by_nf[(r["N"], r["family"])].append(r)

    report = []
    report.append("# Prediction A — F9 Refinement\n")
    report.append("F9 = logH − μ·N·R + γ·N·(R−R*)² − λ·Σ_hist + η·ξ_d\n")

    # ── Phase 1: Fine grid search ──
    report.append("## Phase 1: Fine Grid Search\n")
    report.append("Goal: find (μ, R*, γ, λ, η) where 4D is rank 1 at ALL N and ALL pairwise >80%.\n")

    best_overall = None
    best_score = -1

    candidates = []

    for mu in np.arange(0.0, 2.5, 0.25):
        for R_star in np.arange(0.10, 0.40, 0.025):
            for gamma in [2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]:
                for lam in [0.0, 5.0, 10.0]:
                    for eta in [0.0, 0.3, 0.6]:
                        ranks, wr = evaluate_params(by_nf, n_values, lor_families,
                                                    mu, R_star, gamma, lam, eta)
                        if len(ranks) != 5:
                            continue

                        # Score: weighted by rank and win rates
                        rank_score = sum(1 for r in ranks if r == 1)
                        min_wr = {}
                        for pair, rates in wr.items():
                            if rates:
                                min_wr[pair] = min(rates)
                            else:
                                min_wr[pair] = 0

                        # Composite score: all 3 pairs minimum win rate
                        min_all = min(min_wr.values()) if min_wr else 0
                        score = rank_score * 10 + min_all * 100

                        if score > best_score:
                            best_score = score
                            best_overall = (mu, R_star, gamma, lam, eta, ranks, wr)

                        # Also collect good candidates
                        if rank_score >= 4 and min_all >= 0.5:
                            candidates.append((score, mu, R_star, gamma, lam, eta, ranks, wr))

    # Sort and show top 20
    candidates.sort(key=lambda x: -x[0])
    report.append("### Top 20 Parameter Sets\n")
    report.append("| # | μ | R* | γ | λ | η | ranks | min(4D<2D) | min(4D<3D) | min(4D<5D) | score |")
    report.append("|---|---|-----|---|---|---|-------|------------|------------|------------|-------|")

    for i, (sc, mu, Rs, gam, lam, eta, ranks, wr) in enumerate(candidates[:20]):
        key_pairs = [("Lor4D", "Lor2D"), ("Lor4D", "Lor3D"), ("Lor4D", "Lor5D")]
        min_wrs = []
        for p in key_pairs:
            rates = wr.get(p, [])
            min_wrs.append(min(rates) if rates else 0)

        rank_str = "/".join(str(r) for r in ranks)
        report.append(f"| {i+1} | {mu:.2f} | {Rs:.3f} | {gam:.0f} | {lam:.0f} | {eta:.1f} | "
                      f"{rank_str} | {min_wrs[0]:.0%} | {min_wrs[1]:.0%} | {min_wrs[2]:.0%} | {sc:.0f} |")

    # ── Phase 2: Detailed analysis of the best ──
    if best_overall:
        mu, Rs, gam, lam, eta, ranks, wr = best_overall
        report.append(f"\n## Phase 2: Best Parameter Set Analysis\n")
        report.append(f"**F9** = logH − {mu:.2f}·N·R + {gam:.0f}·N·(R − {Rs:.3f})² "
                      f"− {lam:.0f}·Σ_hist + {eta:.1f}·ξ_d\n")

        # Per-N details
        report.append("### Per-N Dimension Ordering\n")
        report.append("| N | F9(2D) | F9(3D) | F9(4D) | F9(5D) | ordering | 4D_rank |")
        report.append("|---|--------|--------|--------|--------|----------|---------|")

        for N in n_values:
            means = {}
            for f in lor_families:
                vals = by_nf.get((N, f), [])
                if vals:
                    f9s = [compute_f9(r, mu, Rs, gam, lam, eta) for r in vals]
                    means[f] = np.mean(f9s)
            if len(means) == 4:
                rank_sorted = sorted(means.items(), key=lambda x: x[1])
                rank_4d = next(i + 1 for i, (f, _) in enumerate(rank_sorted) if f == "Lor4D")
                order_str = " < ".join(f.replace("Lor", "") for f, _ in rank_sorted)
                report.append(f"| {N} | {means['Lor2D']:.1f} | {means['Lor3D']:.1f} | "
                              f"{means['Lor4D']:.1f} | {means['Lor5D']:.1f} | {order_str} | {rank_4d}/4 |")

        # Per-N pairwise win rates with MWU
        report.append("\n### Per-N Pairwise Win Rates\n")
        key_pairs = [("Lor4D", "Lor2D"), ("Lor4D", "Lor3D"), ("Lor4D", "Lor5D")]
        pair_labels = ["4D<2D", "4D<3D", "4D<5D"]

        report.append("| N | " + " | ".join(pair_labels) + " |")
        report.append("|---|" + "|".join(["----------"] * 3) + "|")

        for N in n_values:
            cells = []
            for left, right in key_pairs:
                lv = by_nf.get((N, left), [])
                rv = by_nf.get((N, right), [])
                if not lv or not rv:
                    cells.append("—")
                    continue
                n_p = min(len(lv), len(rv))
                f9_l = [compute_f9(lv[i], mu, Rs, gam, lam, eta) for i in range(n_p)]
                f9_r = [compute_f9(rv[i], mu, Rs, gam, lam, eta) for i in range(n_p)]
                wins = sum(1 for i in range(n_p) if f9_l[i] < f9_r[i])
                pct = wins / n_p
                u, p = sp_stats.mannwhitneyu(f9_l, f9_r, alternative="less")
                sig = "★★★" if p < 0.001 else "★★" if p < 0.01 else "★" if p < 0.05 else ""
                cells.append(f"{pct:.0%} (p={p:.4f}) {sig}")
            report.append(f"| {N} | {' | '.join(cells)} |")

        # Component decomposition: what does F9 look like?
        report.append("\n### Component Decomposition\n")
        report.append("| N | family | logH | −μNR | γN(R−R*)² | −λΣ | ηξ | F9 |")
        report.append("|---|--------|------|------|-----------|-----|-----|-----|")

        for N in n_values:
            for f in lor_families:
                vals = by_nf.get((N, f), [])
                if not vals:
                    continue
                m_lh = np.mean([r["log_H"] for r in vals])
                m_muNR = -mu * N * np.mean([r["R"] for r in vals])
                m_quad = gam * N * np.mean([(r["R"] - Rs) ** 2 for r in vals])
                m_lam = -lam * np.mean([r["sigma_hist"] for r in vals])
                m_eta = eta * np.mean([r["xi_dim"] for r in vals])
                m_f9 = m_lh + m_muNR + m_quad + m_lam + m_eta
                report.append(f"| {N} | {f} | {m_lh:.1f} | {m_muNR:.1f} | {m_quad:.1f} | "
                              f"{m_lam:.1f} | {m_eta:.1f} | {m_f9:.1f} |")

    # ── Phase 3: Physical interpretation ──
    report.append("\n## Phase 3: Physical Interpretation\n")
    report.append("The quadratic R-well γ·N·(R−R*)² has a clear physical meaning:\n")
    report.append("- R (occupancy/density) is a **dimension proxy**: R(2D)≈0.87 > R(3D)≈0.63 > R(4D)≈0.36 > R(5D)≈0.15")
    report.append("- The quadratic well creates a **preferred density band** around R*")
    report.append("- Deviations from R* are penalized proportionally to N (grows O(N))")
    report.append("- If R* ≈ R(4D), then 4D sits at the well minimum while 2D/3D (too dense) and 5D (too sparse) are penalized")
    report.append("- The −μ·N·R linear term shifts the well center and slope")
    report.append("\n**Physical correspondence**:")
    report.append("- In causal set theory, R = 1 − f_link relates to the causal connectivity density")
    report.append("- R(d) = 1 − E[exp(−ρ·c_d·τ^d)] where c_d is a dimension-dependent volume factor")
    report.append("- The quadratic well R(d)=R* effectively selects d* by matching causal connectivity")
    report.append("- This is analogous to the Einstein-Hilbert action selecting d=4 through the balance")
    report.append("  between R·√g (curvature prefers low d) and the path integral measure (entropy prefers high d)")

    # ── Also check: does F9 still satisfy Prediction B? ──
    if best_overall:
        mu, Rs, gam, lam, eta, _, _ = best_overall
        report.append(f"\n## Phase 4: F9 Cross-Check with Prediction B\n")
        report.append("Does F9(Lor) < F9(KR) hold?\n")
        report.append("| N | pair | win% | sig |")
        report.append("|---|------|------|-----|")

        for N in n_values:
            kr_vals = by_nf.get((N, "KR_like"), [])
            if not kr_vals:
                continue
            f9_kr = [compute_f9(r, mu, Rs, gam, lam, eta) for r in kr_vals]

            for lor in lor_families:
                lv = by_nf.get((N, lor), [])
                if not lv:
                    continue
                f9_lor = [compute_f9(r, mu, Rs, gam, lam, eta) for r in lv]
                n_p = min(len(f9_lor), len(f9_kr))
                wins = sum(1 for i in range(n_p) if f9_lor[i] < f9_kr[i])
                pct = wins / n_p
                u, p = sp_stats.mannwhitneyu(f9_lor[:n_p], f9_kr[:n_p], alternative="less")
                sig = "★★★" if p < 0.001 else "★★" if p < 0.01 else "★" if p < 0.05 else ""
                report.append(f"| {N} | {lor}<KR | {pct:.0%} | {sig} |")

    out = "\n".join(report)
    outdir = Path("outputs_d_recovery")
    md_path = outdir / "prediction_a_f9_refinement.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
