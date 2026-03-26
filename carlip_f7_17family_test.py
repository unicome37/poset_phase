"""
Carlip Response: F7 Functional on All 17 Families
===================================================
Tests whether F7 (the current unified functional) correctly ranks
Lorentzian-like structures above 2-layer, 4-layer, and random layered
structures — addressing Carlip's critique about cherry-picking families.

Key question: Does F7's sigmoid wall + Σ_hist + Ξ_d machinery
discriminate Lorentzian structure from the combinatorially dominant
layered structures that Carlip points out?
"""
from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from bd_action import count_intervals_fast
from generators import (
    Poset,
    generate_kr_like,
    generate_kr_2layer,
    generate_kr_4layer,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_transitive_percolation,
    generate_interval_order,
    generate_absolute_layered,
    generate_multi_layer_random,
    generate_random_layered_k4_uniform,
    generate_random_layered_k6_uniform,
    generate_random_layered_k8_uniform,
    generate_random_layered_k6_tapered,
    generate_random_layered_k6_middle_heavy,
    generate_random_layered_k6_longjump,
)
from unified_functional import (
    compute_log_H,
    compute_pi_geo,
    compute_sigma_hist,
    compute_xi_dim,
)


def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_R(poset: Poset) -> float:
    counts = count_intervals_fast(poset, k_max=3)
    if counts.total_relations <= 0:
        return 0.0
    return 1.0 - float(counts.get(0)) / float(counts.total_relations)


def compute_F7(poset, N, sis_runs=512):
    alpha0, q, lam, eta = 16.0, -0.5, 10.0, 0.6
    Rc, w, N0 = 0.25, 0.015, 20.0

    log_H = compute_log_H(poset, n_runs=sis_runs)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim_val, d_eff = compute_xi_dim(poset)
    R = compute_R(poset)

    alpha_N = alpha0 * (N0 / max(N, 1)) ** abs(q)
    wall = alpha_N * sigmoid((R - Rc) / w)
    F7 = log_H - lam * sigma_hist + eta * xi_dim_val + wall

    return {
        "F7": F7, "log_H": log_H, "sigma_hist": sigma_hist,
        "xi_dim": xi_dim_val, "d_eff": d_eff, "R": R,
        "wall": wall, "alpha_N": alpha_N,
    }


# All 17 families
FAMILIES = {
    # Lorentzian (from manifold sprinklings)
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    # KR and variants (combinatorially dominant)
    "KR_like": generate_kr_like,
    "KR_2layer": generate_kr_2layer,
    "KR_4layer": generate_kr_4layer,
    # Random layered (Carlip's "next most common" structures)
    "AbsLayer": generate_absolute_layered,
    "MLR": generate_multi_layer_random,
    "RLk4": generate_random_layered_k4_uniform,
    "RLk6": generate_random_layered_k6_uniform,
    "RLk8": generate_random_layered_k8_uniform,
    "RLk6_tap": generate_random_layered_k6_tapered,
    "RLk6_mid": generate_random_layered_k6_middle_heavy,
    "RLk6_lj": generate_random_layered_k6_longjump,
    # Other
    "TransPerc": generate_transitive_percolation,
    "IntOrder": generate_interval_order,
}

CATEGORY = {}
for f in ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]:
    CATEGORY[f] = "Lorentzian"
for f in ["KR_like", "KR_2layer", "KR_4layer"]:
    CATEGORY[f] = "KR-family"
for f in ["AbsLayer", "MLR", "RLk4", "RLk6", "RLk8", "RLk6_tap", "RLk6_mid", "RLk6_lj"]:
    CATEGORY[f] = "Layered"
for f in ["TransPerc", "IntOrder"]:
    CATEGORY[f] = "Other"


def main():
    N_VALUES = [16, 20, 28, 36]
    REPS = 8
    SIS_RUNS = 256
    SEED_BASE = 42

    print("=" * 80)
    print("CARLIP RESPONSE: F7 on All 17 Families")
    print("=" * 80)
    print(f"N = {N_VALUES}, reps = {REPS}, SIS runs = {SIS_RUNS}")
    print(f"Families: {len(FAMILIES)}")
    print()

    all_rows = []
    total = len(FAMILIES) * len(N_VALUES) * REPS
    done = 0

    for fam_name, gen_fn in FAMILIES.items():
        for N in N_VALUES:
            for rep in range(REPS):
                seed = SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep
                seed = seed % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    result = compute_F7(poset, N, sis_runs=SIS_RUNS)
                    result["family"] = fam_name
                    result["category"] = CATEGORY[fam_name]
                    result["N"] = N
                    result["rep"] = rep
                    all_rows.append(result)
                except Exception as e:
                    print(f"  ERROR: {fam_name} N={N} rep={rep}: {e}")

                done += 1
                if done % 50 == 0:
                    print(f"  [{done}/{total}] {fam_name} N={N}", flush=True)

    # ═══ Analysis ═══
    report = []
    report.append("# Carlip Response: F7 on All 17 Families\n")
    report.append(f"N = {N_VALUES}, reps = {REPS}, SIS = {SIS_RUNS}\n")

    # Group by (N, family) → mean F7
    by_nf = defaultdict(list)
    for r in all_rows:
        by_nf[(r["N"], r["family"])].append(r)

    for N in N_VALUES:
        report.append(f"\n## N = {N}\n")
        report.append("| Rank | Family | Category | F7 mean | logH | Σ_hist | d_eff | R | wall |")
        report.append("|------|--------|----------|---------|------|--------|-------|---|------|")

        means = {}
        details = {}
        for fam in FAMILIES:
            rows = by_nf.get((N, fam), [])
            if rows:
                f7s = [r["F7"] for r in rows]
                means[fam] = np.mean(f7s)
                details[fam] = {
                    "F7": np.mean(f7s),
                    "log_H": np.mean([r["log_H"] for r in rows]),
                    "sigma_hist": np.mean([r["sigma_hist"] for r in rows]),
                    "d_eff": np.mean([r["d_eff"] for r in rows]),
                    "R": np.mean([r["R"] for r in rows]),
                    "wall": np.mean([r["wall"] for r in rows]),
                }

        ranked = sorted(means, key=means.get)
        for rank, fam in enumerate(ranked, 1):
            d = details[fam]
            cat = CATEGORY[fam]
            tag = ""
            if cat == "Lorentzian":
                tag = " ◆"
            elif cat == "KR-family":
                tag = " ●"
            report.append(
                f"| {rank} | {fam}{tag} | {cat} | {d['F7']:+.2f} | "
                f"{d['log_H']:.1f} | {d['sigma_hist']:.3f} | "
                f"{d['d_eff']:.2f} | {d['R']:.3f} | {d['wall']:.2f} |"
            )

    # ═══ Pairwise: Lor vs each non-Lor ═══
    report.append("\n## Pairwise Win Rates: Lor4D vs Each Family\n")
    report.append("| N | Opponent | F7(Lor4D)<F7(opp) % | ΔF7 mean | interpretation |")
    report.append("|---|----------|--------------------:|----------|----------------|")

    for N in N_VALUES:
        lor4d_rows = by_nf.get((N, "Lor4D"), [])
        if not lor4d_rows:
            continue
        lor4d_f7s = [r["F7"] for r in lor4d_rows]
        lor4d_mean = np.mean(lor4d_f7s)

        for opp in sorted(FAMILIES.keys()):
            if opp == "Lor4D":
                continue
            opp_rows = by_nf.get((N, opp), [])
            if not opp_rows:
                continue
            opp_f7s = [r["F7"] for r in opp_rows]
            opp_mean = np.mean(opp_f7s)

            # Pairwise comparison (all pairs)
            wins = 0
            total_pairs = 0
            for l in lor4d_f7s:
                for o in opp_f7s:
                    total_pairs += 1
                    if l < o:
                        wins += 1
            win_rate = wins / total_pairs if total_pairs > 0 else 0
            delta = lor4d_mean - opp_mean
            interp = "✅ Lor4D wins" if win_rate > 0.6 else ("❌ Lor4D loses" if win_rate < 0.4 else "≈ tied")
            report.append(f"| {N} | {opp} | {win_rate:.1%} | {delta:+.2f} | {interp} |")

    # ═══ Key diagnostic: What drives the ranking? ═══
    report.append("\n## Diagnostic: Component Decomposition at N=20\n")
    report.append("| Family | Category | logH | -λΣ_hist | wall | F7 | logH rank | F7 rank |")
    report.append("|--------|----------|------|---------|------|-----|-----------|---------|")

    N_diag = 20
    fam_details_20 = {}
    for fam in FAMILIES:
        rows = by_nf.get((N_diag, fam), [])
        if rows:
            fam_details_20[fam] = {
                "log_H": np.mean([r["log_H"] for r in rows]),
                "neg_lam_sig": -10.0 * np.mean([r["sigma_hist"] for r in rows]),
                "wall": np.mean([r["wall"] for r in rows]),
                "F7": np.mean([r["F7"] for r in rows]),
            }

    logh_ranked = sorted(fam_details_20, key=lambda f: fam_details_20[f]["log_H"])
    f7_ranked = sorted(fam_details_20, key=lambda f: fam_details_20[f]["F7"])

    for fam in f7_ranked:
        d = fam_details_20[fam]
        lr = logh_ranked.index(fam) + 1
        fr = f7_ranked.index(fam) + 1
        report.append(
            f"| {fam} | {CATEGORY[fam]} | {d['log_H']:.1f} | "
            f"{d['neg_lam_sig']:.1f} | {d['wall']:.2f} | "
            f"{d['F7']:.1f} | {lr} | {fr} |"
        )

    # ═══ Conclusion ═══
    report.append("\n## Conclusion\n")
    report.append("### Does F7 address Carlip's critique?\n")
    report.append("- If Lor4D ranks below random layered/KR variants → F7 still cherry-picks")
    report.append("- If Lor4D ranks above ALL non-Lorentzian → F7 genuinely selects geometry")
    report.append("- If mixed → F7 partially selects but has blind spots\n")

    # Check: does any Lor family beat all non-Lor families?
    for N in N_VALUES:
        lor_fams = [f for f in FAMILIES if CATEGORY[f] == "Lorentzian"]
        non_lor_fams = [f for f in FAMILIES if CATEGORY[f] != "Lorentzian"]

        for lf in lor_fams:
            lr = by_nf.get((N, lf), [])
            if not lr:
                continue
            lf_mean = np.mean([r["F7"] for r in lr])
            beats_all = True
            for nf in non_lor_fams:
                nr = by_nf.get((N, nf), [])
                if nr:
                    nf_mean = np.mean([r["F7"] for r in nr])
                    if lf_mean >= nf_mean:
                        beats_all = False
                        break
            if beats_all:
                report.append(f"- N={N}: **{lf} beats ALL non-Lorentzian** ✅")
            else:
                # Find who beats it
                losers = []
                for nf in non_lor_fams:
                    nr = by_nf.get((N, nf), [])
                    if nr:
                        nf_mean = np.mean([r["F7"] for r in nr])
                        if lf_mean >= nf_mean:
                            losers.append(f"{nf}({nf_mean:.1f})")
                report.append(f"- N={N}: {lf}(F7={lf_mean:.1f}) LOSES to: {', '.join(losers)}")

    out = "\n".join(report)
    outdir = Path("outputs_carlip")
    outdir.mkdir(exist_ok=True)
    md_path = outdir / "f7_17family_test.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"\nSaved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
