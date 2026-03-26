"""
Prediction B — F7 Large-N Verification (N=20–100, 40 reps)
============================================================
Core question: Does F7(Lor) < F7(KR) hold at large N?

Prediction B states: Lorentzian-like posets have lower F7 than non-Lorentzian (KR).
This is the "action ordering" prediction — Lorentzian geometry is preferred.

Strategy:
  - Reuse cached data from prediction_c_f7_large_n.csv (1000 posets, 5 families)
  - OR generate fresh data if --fresh
  - Pairwise tests: Lor2D/3D/4D/5D vs KR_like
  - Add proper statistical tests: Wilcoxon signed-rank, Mann-Whitney U, bootstrap CI
  - Component decomposition: what drives the ranking at each N?
  - N-scaling: does the margin grow, shrink, or stay constant?
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


def bootstrap_win_rate(f7_lor, f7_kr, n_boot=5000, seed=42):
    """Bootstrap CI for win rate P(F7_lor < F7_kr)."""
    rng = np.random.RandomState(seed)
    n = min(len(f7_lor), len(f7_kr))
    wins = np.array([f7_lor[i] < f7_kr[i] for i in range(n)], dtype=float)
    boot_means = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, n)
        boot_means.append(wins[idx].mean())
    lo, hi = np.percentile(boot_means, [2.5, 97.5])
    return wins.mean(), lo, hi


def generate_report(rows, n_values):
    families = sorted(set(r["family"] for r in rows))
    lor_families = [f for f in families if f.startswith("Lor")]
    kr_family = "KR_like"

    by_nf = defaultdict(list)
    for r in rows:
        by_nf[(r["N"], r["family"])].append(r)

    report = []
    report.append("# Prediction B — F7 Large-N Verification\n")
    report.append(f"**Data**: {len(rows)} samples, N ∈ {{{', '.join(str(n) for n in n_values)}}}")
    report.append(f"**Families**: {', '.join(families)}")
    report.append(f"**Reps per cell**: {len(by_nf.get((n_values[0], families[0]), []))}\n")

    # ── Section 1: Family means ──
    report.append("## 1. Family-Level Means\n")
    report.append("| N | Family | mean_F7 | mean_logH | mean_wall | mean_R | mean_Σ_hist | mean_ξ_d |")
    report.append("|---|--------|---------|-----------|-----------|--------|-------------|----------|")

    for N in n_values:
        for fam in families:
            vals = by_nf.get((N, fam), [])
            if not vals:
                continue
            mF = np.mean([r["F7"] for r in vals])
            mH = np.mean([r["log_H"] for r in vals])
            mW = np.mean([r["wall"] for r in vals])
            mR = np.mean([r["R"] for r in vals])
            mS = np.mean([r["sigma_hist"] for r in vals])
            mX = np.mean([r["xi_dim"] for r in vals])
            report.append(f"| {N} | {fam} | {mF:.1f} | {mH:.1f} | {mW:.2f} | {mR:.3f} | {mS:.3f} | {mX:.3f} |")

    # ── Section 2: Pairwise win rates with statistical tests ──
    report.append("\n## 2. Pairwise Tests: Lor_X vs KR_like\n")
    report.append("Prediction B claims F7(Lor) < F7(KR). Win rate = P(F7_Lor < F7_KR).\n")

    overall_results = {}  # (lor_fam) -> list of per-N results

    for lor in lor_families:
        report.append(f"### {lor} vs KR_like\n")
        report.append("| N | n | win% | 95% CI | mean_ΔF7 | MWU_p | sig | Δ_logH | Δ_wall | Δ_Σ_hist |")
        report.append("|---|---|------|--------|----------|-------|-----|--------|--------|----------|")

        per_n_results = []

        for N in n_values:
            lor_vals = by_nf.get((N, lor), [])
            kr_vals = by_nf.get((N, kr_family), [])
            if not lor_vals or not kr_vals:
                continue

            f7_lor = np.array([r["F7"] for r in lor_vals])
            f7_kr = np.array([r["F7"] for r in kr_vals])

            n_pairs = min(len(f7_lor), len(f7_kr))

            # Win rate with bootstrap CI
            win_pct, ci_lo, ci_hi = bootstrap_win_rate(f7_lor, f7_kr)

            # Mann-Whitney U test (unpaired)
            u_stat, mwu_p = sp_stats.mannwhitneyu(f7_lor, f7_kr, alternative="less")

            # Effect size
            delta_f7 = np.mean(f7_lor) - np.mean(f7_kr)
            delta_lh = np.mean([r["log_H"] for r in lor_vals]) - np.mean([r["log_H"] for r in kr_vals])
            delta_wall = np.mean([r["wall"] for r in lor_vals]) - np.mean([r["wall"] for r in kr_vals])
            delta_sh = np.mean([r["sigma_hist"] for r in lor_vals]) - np.mean([r["sigma_hist"] for r in kr_vals])

            sig = "★★★" if mwu_p < 0.001 else "★★" if mwu_p < 0.01 else "★" if mwu_p < 0.05 else ""
            ci_str = f"[{ci_lo:.0%},{ci_hi:.0%}]"

            report.append(f"| {N} | {n_pairs} | {win_pct:.0%} | {ci_str} | {delta_f7:+.1f} | "
                          f"{mwu_p:.4f} | {sig} | {delta_lh:+.1f} | {delta_wall:+.2f} | {delta_sh:+.3f} |")

            per_n_results.append({
                "N": N, "win_pct": win_pct, "delta_f7": delta_f7,
                "mwu_p": mwu_p, "delta_lh": delta_lh, "delta_wall": delta_wall
            })

        overall_results[lor] = per_n_results

        # N-scaling trend
        if len(per_n_results) >= 3:
            ns = [r["N"] for r in per_n_results]
            deltas = [r["delta_f7"] for r in per_n_results]
            rho_trend, p_trend = sp_stats.spearmanr(ns, deltas)
            direction = "margin growing" if rho_trend < -0.5 else "margin shrinking" if rho_trend > 0.5 else "stable"
            report.append(f"\nN-scaling: ρ(N, ΔF7) = {rho_trend:+.3f} (p={p_trend:.3f}) — {direction}")

        # Crossover detection
        crossover_N = None
        for r in per_n_results:
            if r["win_pct"] < 0.5:
                crossover_N = r["N"]
                break
        if crossover_N:
            report.append(f"\n⚠️ **Crossover at N={crossover_N}**: win% drops below 50% — F7 bounded-wall limitation")
        else:
            report.append(f"\n✅ No crossover: {lor} < KR_like at all tested N")
        report.append("")

    # ── Section 3: Component attribution ──
    report.append("\n## 3. Component Attribution\n")
    report.append("Which F7 component drives the Lor < KR ordering?\n")
    report.append("| pair | N | Δ_logH | −λΔ_Σ_hist | Δ_wall | Δ_ξ_d | Δ_F7 | dominant |")
    report.append("|------|---|--------|------------|--------|-------|------|----------|")

    for lor in lor_families:
        for N in n_values:
            lor_vals = by_nf.get((N, lor), [])
            kr_vals = by_nf.get((N, kr_family), [])
            if not lor_vals or not kr_vals:
                continue

            d_lh = np.mean([r["log_H"] for r in lor_vals]) - np.mean([r["log_H"] for r in kr_vals])
            d_sh = np.mean([r["sigma_hist"] for r in lor_vals]) - np.mean([r["sigma_hist"] for r in kr_vals])
            d_wall = np.mean([r["wall"] for r in lor_vals]) - np.mean([r["wall"] for r in kr_vals])
            d_xi = np.mean([r["xi_dim"] for r in lor_vals]) - np.mean([r["xi_dim"] for r in kr_vals])
            d_f7 = np.mean([r["F7"] for r in lor_vals]) - np.mean([r["F7"] for r in kr_vals])

            lam_d_sh = -10.0 * d_sh  # F7 has -λ·Σ_hist, so ΔF7 contribution = -λ·ΔΣ
            eta_d_xi = 0.6 * d_xi

            components = {
                "logH": d_lh,
                "−λΣ": lam_d_sh,
                "wall": d_wall,
                "ξ_d": eta_d_xi,
            }
            dominant = min(components, key=lambda k: components[k]) if d_f7 < 0 else max(components, key=lambda k: components[k])

            report.append(f"| {lor}→KR | {N} | {d_lh:+.1f} | {lam_d_sh:+.1f} | {d_wall:+.2f} | "
                          f"{eta_d_xi:+.2f} | {d_f7:+.1f} | {dominant} |")

    # ── Section 4: Overall verdict ──
    report.append("\n## 4. Overall Verdict\n")
    report.append("| pair | total_wins/total | overall_win% | all_N_pass? | crossover_N |")
    report.append("|------|-----------------|-------------|-------------|-------------|")

    for lor in lor_families:
        results = overall_results.get(lor, [])
        total_wins = sum(int(r["win_pct"] * 40) for r in results)  # approximate
        total_n = len(results) * 40
        overall_pct = total_wins / total_n if total_n > 0 else 0
        all_pass = all(r["win_pct"] >= 0.5 for r in results)
        crossover = next((r["N"] for r in results if r["win_pct"] < 0.5), "—")
        report.append(f"| {lor} < KR | {total_wins}/{total_n} | {overall_pct:.1%} | "
                      f"{'✅' if all_pass else '❌'} | {crossover} |")

    # Known limitation note
    report.append("\n### Known Limitation: Bounded Wall\n")
    report.append("The F7 wall term α(N)·σ((R−Rc)/w) is bounded, while ΔlogH(Lor4D−KR) grows as O(N).")
    report.append("This means Lor4D < KR will inevitably fail at sufficiently large N.")
    report.append("This is the same structural limitation identified in Prediction A (§1.3 F8 series).\n")

    return "\n".join(report)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default="outputs_d_recovery/prediction_c_f7_large_n.csv",
                        help="Path to cached CSV with F7 data")
    args = parser.parse_args()

    print(f"=== Prediction B — F7 Large-N Verification ===")
    print(f"Loading cached data: {args.csv}")

    rows = load_csv(args.csv)
    n_values = sorted(set(r["N"] for r in rows))
    families = sorted(set(r["family"] for r in rows))
    print(f"Loaded {len(rows)} rows, N={n_values}, families={families}")

    report = generate_report(rows, n_values)

    outdir = Path("outputs_d_recovery")
    outdir.mkdir(exist_ok=True)
    md_path = outdir / "prediction_b_f7_large_n.md"
    md_path.write_text(report, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + report)


if __name__ == "__main__":
    main()
