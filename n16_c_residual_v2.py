"""N=16 C residual check v2: uses ORIGINAL C definition (all families at each N).

The original _c_repair_experiment.py computed C across ALL families at each N,
not just within Lor2D. This is a harder test because the wall term varies
wildly across families (different R values).

This script tests both definitions under F7 and F8a_v3.
"""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from bd_action import count_intervals_fast
from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_kr_like,
)
from unified_functional import compute_pi_geo, compute_sigma_hist, compute_xi_dim
from entropy_sis import estimate_log_linear_extensions_sis
from observables_geo import estimate_dimension_proxy_from_order_fraction

OUT_DIR = Path("outputs_unified_functional")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_R(poset: Poset) -> float:
    counts = count_intervals_fast(poset, k_max=3)
    if counts.total_relations <= 0:
        return 0.0
    return 1.0 - float(counts.get(0)) / float(counts.total_relations)


def compute_d_eff(poset: Poset) -> float:
    c = poset.closure
    n = poset.n
    pairs = n * (n - 1) // 2
    if pairs == 0:
        return 0.0
    comp = int(c.sum()) - n
    R_val = comp / (2 * pairs)
    return estimate_dimension_proxy_from_order_fraction(R_val)


def compute_log_H_sis(poset: Poset, seed: int = 0) -> float:
    mean_est, _ = estimate_log_linear_extensions_sis(poset, n_runs=2048, seed=seed)
    return mean_est


GENERATORS = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}


def F_model(row: dict, model: str) -> float:
    N = row["N"]
    log_H = row["log_H"]
    R = row["R"]
    pi_geo = row["pi_geo"]
    sigma_hist = row["sigma_hist"]
    xi_dim = row["xi_dim"]
    d_eff = row["d_eff"]

    if model == "F7":
        alpha_N = 16.0 * (20.0 / max(N, 1)) ** 0.5
        wall = alpha_N * sigmoid((R - 0.25) / 0.015)
        return log_H + 0.0004 * pi_geo - 10.0 * sigma_hist + 0.6 * xi_dim + wall

    elif model == "F8a_v3":
        nlogn = N * math.log(N) if N > 1 else 1.0
        h = log_H / nlogn
        alpha_N = 16.0 * (20.0 / max(N, 1)) ** 0.5
        wall_up = alpha_N * sigmoid((R - 0.25) / 0.015)
        alpha_Nf = 12.0 * (20.0 / max(N, 1)) ** 0.5
        wall_lo = alpha_Nf * sigmoid((0.18 - R) / 0.02)
        xi_d_plus = 120.0 * max(0.0, d_eff - 4.0) ** 3
        return (h + 0.0004 * pi_geo - 10.0 * sigma_hist
                + 0.6 * xi_dim + xi_d_plus + wall_up + wall_lo)

    raise ValueError(f"Unknown model: {model}")


def main():
    N_VALUES = [16, 20, 28, 36]
    REPS = 10  # match original scan's rep count
    MODELS = ["F7", "F8a_v3"]
    SEED_BASE = 42  # match original scan's seed pattern

    print("=" * 70)
    print("N=16 C Residual Check v2: Original Definition (ALL families)")
    print("=" * 70)

    # Generate data matching original seed convention
    all_data: list[dict] = []
    for N in N_VALUES:
        print(f"\nGenerating N={N}...")
        for fam_name, gen_fn in GENERATORS.items():
            for rep in range(REPS):
                seed = 42 + rep * 1000 + N * 100
                poset = gen_fn(N, seed=seed)
                log_H = compute_log_H_sis(poset, seed=seed * 7 + 13)
                R = compute_R(poset)
                pi_geo = compute_pi_geo(poset)
                sigma_hist = compute_sigma_hist(poset)
                xi_dim, _ = compute_xi_dim(poset)
                d_eff = compute_d_eff(poset)

                all_data.append({
                    "N": N, "family": fam_name, "rep": rep,
                    "log_H": log_H, "R": R,
                    "pi_geo": pi_geo, "sigma_hist": sigma_hist,
                    "xi_dim": xi_dim, "d_eff": d_eff,
                })

    print(f"\nTotal data: {len(all_data)} rows")

    lines = ["# N=16 C Residual Check v2\n"]
    lines.append("**Date:** 2026-03-23\n")
    lines.append("**Purpose:** Test C with TWO definitions under F7 and F8a_v3\n\n")

    for model in MODELS:
        lines.append(f"## {model}\n")

        # ── Definition 1: ALL families at each N (original) ──
        lines.append("### Definition 1: ALL families at each N (original weight scan)\n")
        lines.append("| N | n | ρ(Σ_hist, F) | p-value | direction |")
        lines.append("|---|---|-------------|---------|-----------|")

        c1_pass = 0
        c1_total = 0
        for N in N_VALUES:
            grp = [r for r in all_data if r["N"] == N]
            sh = np.array([r["sigma_hist"] for r in grp])
            fs = np.array([F_model(r, model) for r in grp])
            if sh.std() < 1e-10:
                lines.append(f"| {N} | {len(grp)} | – | – | no variance |")
                continue
            rho, pval = sp_stats.spearmanr(sh, fs)
            ok = "✅" if rho < 0 else "❌"
            lines.append(f"| {N} | {len(grp)} | {rho:+.4f} | {pval:.3e} | {ok} |")
            c1_total += 1
            if rho < 0:
                c1_pass += 1
            print(f"  {model} Def1 N={N}: rho={rho:+.4f} p={pval:.3e} {'✓' if rho < 0 else '✗'}")

        lines.append(f"\n**C (Def1) = {c1_pass}/{c1_total} = {c1_pass/max(c1_total,1):.2f}**\n")

        # ── Definition 2: Lor2D only at each N ──
        lines.append("### Definition 2: Lor2D only at each N (f8a_acd style)\n")
        lines.append("| N | n | ρ(Σ_hist, F) | p-value | direction |")
        lines.append("|---|---|-------------|---------|-----------|")

        c2_pass = 0
        c2_total = 0
        for N in N_VALUES:
            grp = [r for r in all_data if r["N"] == N and r["family"] == "Lor2D"]
            sh = np.array([r["sigma_hist"] for r in grp])
            fs = np.array([F_model(r, model) for r in grp])
            if sh.std() < 1e-10 or len(grp) < 5:
                lines.append(f"| {N} | {len(grp)} | – | – | insufficient |")
                continue
            rho, pval = sp_stats.spearmanr(sh, fs)
            ok = "✅" if rho < 0 else "❌"
            lines.append(f"| {N} | {len(grp)} | {rho:+.4f} | {pval:.3e} | {ok} |")
            c2_total += 1
            if rho < 0:
                c2_pass += 1
            print(f"  {model} Def2 N={N}: rho={rho:+.4f} p={pval:.3e} {'✓' if rho < 0 else '✗'}")

        lines.append(f"\n**C (Def2) = {c2_pass}/{c2_total} = {c2_pass/max(c2_total,1):.2f}**\n")

        # ── Component decomposition at N=16 for Def1 ──
        lines.append(f"### {model}: N=16 Component Analysis (all families)\n")
        grp = [r for r in all_data if r["N"] == 16]
        sh = np.array([r["sigma_hist"] for r in grp])
        fams = np.array([r["family"] for r in grp])
        R_vals = np.array([r["R"] for r in grp])
        lh = np.array([r["log_H"] for r in grp])

        lines.append("| Family | mean(R) | mean(logH) | mean(Σ_hist) | mean(F) |")
        lines.append("|--------|---------|-----------|-------------|---------|")
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
            idx = [i for i, r in enumerate(grp) if r["family"] == fam]
            if not idx:
                continue
            sub = [grp[i] for i in idx]
            mR = np.mean([r["R"] for r in sub])
            mH = np.mean([r["log_H"] for r in sub])
            mS = np.mean([r["sigma_hist"] for r in sub])
            mF = np.mean([F_model(r, model) for r in sub])
            lines.append(f"| {fam} | {mR:.4f} | {mH:.2f} | {mS:.4f} | {mF:.4f} |")

        lines.append("")

        # Per-family corr(Σ_hist, wall)
        if model == "F7":
            alpha_16 = 16.0 * (20.0 / 16.0) ** 0.5
            wall_vals = np.array([alpha_16 * sigmoid((r["R"] - 0.25) / 0.015) for r in grp])
            rho_wall, _ = sp_stats.spearmanr(sh, wall_vals)
            rho_lh, _ = sp_stats.spearmanr(sh, lh)
            lines.append(f"**corr(Σ_hist, wall) = {rho_wall:+.4f}**, corr(Σ_hist, logH) = {rho_lh:+.4f}\n")
            lines.append(f"**std(wall) = {wall_vals.std():.4f}**, std(logH) = {lh.std():.4f}\n")

    # ── Also test with 50 reps to match v1 script's power ──
    lines.append("\n## Higher Power Test (50 reps, all families)\n")
    print("\nGenerating 50-rep data for higher-power test...")
    all_data_50: list[dict] = []
    for N in N_VALUES:
        for fam_name, gen_fn in GENERATORS.items():
            for rep in range(50):
                seed = 20260323 + N * 1000 + rep + hash(fam_name) % 10000
                poset = gen_fn(N, seed=seed)
                log_H = compute_log_H_sis(poset, seed=seed * 7 + 13)
                R = compute_R(poset)
                pi_geo = compute_pi_geo(poset)
                sigma_hist = compute_sigma_hist(poset)
                xi_dim, _ = compute_xi_dim(poset)
                d_eff = compute_d_eff(poset)
                all_data_50.append({
                    "N": N, "family": fam_name, "rep": rep,
                    "log_H": log_H, "R": R,
                    "pi_geo": pi_geo, "sigma_hist": sigma_hist,
                    "xi_dim": xi_dim, "d_eff": d_eff,
                })

    for model in MODELS:
        lines.append(f"### {model} (50 reps, Def1 = all families)\n")
        lines.append("| N | n | ρ(Σ_hist, F) | p-value | direction |")
        lines.append("|---|---|-------------|---------|-----------|")
        c_pass = 0
        c_total = 0
        for N in N_VALUES:
            grp = [r for r in all_data_50 if r["N"] == N]
            sh = np.array([r["sigma_hist"] for r in grp])
            fs = np.array([F_model(r, model) for r in grp])
            if sh.std() < 1e-10:
                continue
            rho, pval = sp_stats.spearmanr(sh, fs)
            ok = "✅" if rho < 0 else "❌"
            lines.append(f"| {N} | {len(grp)} | {rho:+.4f} | {pval:.3e} | {ok} |")
            c_total += 1
            if rho < 0:
                c_pass += 1
            print(f"  {model} 50rep Def1 N={N}: rho={rho:+.4f} {'✓' if rho < 0 else '✗'}")
        lines.append(f"\n**C = {c_pass}/{c_total} = {c_pass/max(c_total,1):.2f}**\n")

    # Conclusion
    lines.append("\n## Conclusion\n")
    lines.append("(Filled after results)\n")

    report = "\n".join(lines)
    out_path = OUT_DIR / "n16_c_residual_v2.md"
    out_path.write_text(report, encoding="utf-8")
    print(f"\nReport saved to {out_path}")
    print("\n" + "=" * 70)
    print("DONE")


if __name__ == "__main__":
    main()
