"""Conjecture E — Layer 3: F7 ↔ S_BD Coefficient Fitting.

Tests the core Layer 3 question: does F7 decompose as
    F7 = a · S_BD + (family/N corrections) + residual
and does the coefficient `a` stabilize as N → ∞?

For each N separately, we fit:
    F7_i = a(N) · S_BD_i + Σ_j β_j · 1[family_i = j] + ε_i

Then track a(N) across N = 16, 20, 28, 36, 48, 56, 64.

Multiple S_BD variants are tested:
  - bd_ratio (interval richness metric)
  - bdg_d4_standard (BDG d=4 action)
  - bdg_d2_corrected (BDG d=2 with C1 correction)
  - R (interval occupancy = 1 - f_link, same quantity in sigmoid wall)
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from bd_action import (
    bd_ratio_metric,
    bdg_action_d2_corrected,
    bdg_action_d4_standard,
    count_intervals_fast,
)
from conjecture_e_f7_bridge import compute_F7_from_components, compute_R_from_counts
from generators import (
    Poset,
    generate_kr_like,
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


FAMILY_GENS = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}


@dataclass(frozen=True)
class Row:
    family: str
    N: int
    rep: int
    seed: int
    log_H: float
    pi_geo: float
    sigma_hist: float
    xi_dim: float
    R: float
    wall: float
    F7: float
    bd_ratio: float
    bdg_d2c: float
    bdg_d4s: float


def make_seed(family: str, n: int, rep: int, base: int) -> int:
    fam_idx = list(FAMILY_GENS.keys()).index(family)
    return base + fam_idx * 100000 + n * 100 + rep


def run_single(family: str, N: int, rep: int, base_seed: int) -> Row:
    s = make_seed(family, N, rep, base_seed)
    poset = FAMILY_GENS[family](N, seed=s)

    log_H = compute_log_H(poset, n_runs=512)
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim_val, _ = compute_xi_dim(poset)

    counts = count_intervals_fast(poset)
    R = compute_R_from_counts(counts)
    br = bd_ratio_metric(counts)
    bdg2c = bdg_action_d2_corrected(counts, N, normalized=True)
    bdg4s = bdg_action_d4_standard(counts, N, normalized=True)

    F7, wall_val = compute_F7_from_components(
        log_H, pi_geo, sigma_hist, xi_dim_val, R, N
    )

    return Row(
        family=family, N=N, rep=rep, seed=s,
        log_H=log_H, pi_geo=pi_geo, sigma_hist=sigma_hist,
        xi_dim=xi_dim_val, R=R, wall=wall_val, F7=F7,
        bd_ratio=br, bdg_d2c=bdg2c, bdg_d4s=bdg4s,
    )


def ols_with_family(
    y: np.ndarray,
    x_bd: np.ndarray,
    families: list[str],
) -> dict:
    """OLS: y = a*x_bd + Σ β_j · 1[family=j] + intercept.
    Returns dict with a, r2, se_a, family betas."""
    unique_fams = sorted(set(families))
    # Design: [intercept, x_bd, family dummies (drop last)]
    n = len(y)
    ncols = 2 + len(unique_fams) - 1
    X = np.zeros((n, ncols))
    X[:, 0] = 1.0  # intercept
    X[:, 1] = x_bd  # BD metric
    for i, fam in enumerate(families):
        idx = unique_fams.index(fam)
        if idx < len(unique_fams) - 1:
            X[i, 2 + idx] = 1.0

    try:
        beta, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    except np.linalg.LinAlgError:
        return {"a": float("nan"), "r2": float("nan"), "se_a": float("nan")}

    pred = X @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # Standard error of a (coefficient index 1)
    if n > ncols and ss_res > 0:
        mse = ss_res / (n - ncols)
        try:
            cov = mse * np.linalg.inv(X.T @ X)
            se_a = float(np.sqrt(cov[1, 1]))
        except np.linalg.LinAlgError:
            se_a = float("nan")
    else:
        se_a = float("nan")

    return {
        "a": float(beta[1]),
        "intercept": float(beta[0]),
        "r2": r2,
        "se_a": se_a,
        "family_betas": {unique_fams[i]: float(beta[2 + i])
                         for i in range(len(unique_fams) - 1)},
        "reference_family": unique_fams[-1],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Layer 3: F7↔S_BD coefficient fitting")
    ap.add_argument("--ns", nargs="*", type=int, default=[16, 20, 28, 36, 48, 56, 64])
    ap.add_argument("--families", nargs="*", default=["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"])
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default="outputs_unified_functional/conjecture_e_layer3_coefficient.csv")
    ap.add_argument("--report", default="outputs_unified_functional/conjecture_e_layer3_coefficient.md")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[Row] = []
    total = len(args.families) * len(args.ns) * args.reps
    done = 0

    for N in args.ns:
        for fam in args.families:
            for rep in range(args.reps):
                row = run_single(fam, N, rep, args.seed)
                rows.append(row)
                done += 1
                if done % 25 == 0 or done == total:
                    print(f"  [{done}/{total}] {fam} N={N} rep={rep} "
                          f"F7={row.F7:.1f} bd_ratio={row.bd_ratio:.3f}")

    # Save CSV
    fn = [f.name for f in fields(Row)]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fn)
        w.writeheader()
        for r in rows:
            w.writerow({k: getattr(r, k) for k in fn})
    print(f"\nSaved: {out_path}")

    # ── Analysis: coefficient fitting per N ──
    BD_METRICS = {
        "bd_ratio": lambda r: r.bd_ratio,
        "bdg_d4s": lambda r: r.bdg_d4s,
        "bdg_d2c": lambda r: r.bdg_d2c,
        "R": lambda r: r.R,
    }

    report_lines: list[str] = []
    report_lines.append("# Conjecture E — Layer 3: F7↔S_BD Coefficient Fitting\n")
    report_lines.append(f"Total: {len(rows)} realizations\n")

    # For each BD metric, fit per-N and track a(N)
    all_results: dict[str, list[dict]] = {m: [] for m in BD_METRICS}

    for metric_name, metric_fn in BD_METRICS.items():
        report_lines.append(f"\n## Metric: {metric_name}\n")
        report_lines.append("| N | a(N) | SE(a) | R² | n |")
        report_lines.append("|---|------|-------|-----|---|")

        for N in args.ns:
            subset = [r for r in rows if r.N == N]
            if len(subset) < 5:
                continue
            y = np.array([r.F7 for r in subset])
            x = np.array([metric_fn(r) for r in subset])
            fams = [r.family for r in subset]
            result = ols_with_family(y, x, fams)
            result["N"] = N
            result["n"] = len(subset)
            all_results[metric_name].append(result)

            report_lines.append(
                f"| {N} | {result['a']:+.3f} | {result['se_a']:.3f} | "
                f"{result['r2']:.4f} | {result['n']} |"
            )

    # ── Coefficient stability analysis ──
    report_lines.append("\n## Coefficient Stability: a(N) convergence\n")

    for metric_name, results in all_results.items():
        if len(results) < 3:
            continue
        ns_arr = np.array([r["N"] for r in results], dtype=float)
        a_arr = np.array([r["a"] for r in results])
        se_arr = np.array([r["se_a"] for r in results])

        report_lines.append(f"\n### {metric_name}\n")
        report_lines.append(f"a(N) trajectory: " +
                           " → ".join(f"{a:+.3f}" for a in a_arr))

        # Test 1: Is a(N) bounded? (range / mean < threshold)
        a_range = float(np.max(a_arr) - np.min(a_arr))
        a_mean = float(np.mean(a_arr))
        cv = a_range / abs(a_mean) if abs(a_mean) > 1e-12 else float("inf")
        report_lines.append(f"- Range: {a_range:.3f}, Mean: {a_mean:+.3f}, CV(range/|mean|): {cv:.3f}")

        # Test 2: Does a(N) trend toward a limit? (last 3 vs first 3)
        if len(a_arr) >= 6:
            first3 = a_arr[:3]
            last3 = a_arr[-3:]
            diff_first = float(np.std(first3))
            diff_last = float(np.std(last3))
            report_lines.append(f"- std(first 3): {diff_first:.3f}, std(last 3): {diff_last:.3f}")
            if diff_last < diff_first:
                report_lines.append(f"- → **Stabilizing** (std shrinking)")
            else:
                report_lines.append(f"- → Not yet stable")

        # Test 3: Spearman(a, 1/N) — if a converges, it should decorrelate from N
        rho_aN, p_aN = sp_stats.spearmanr(ns_arr, a_arr)
        report_lines.append(f"- Spearman(N, a): ρ={rho_aN:+.3f}, p={p_aN:.3e}")
        if abs(rho_aN) < 0.5:
            report_lines.append(f"- → No significant N-trend (good: a is N-independent)")
        elif rho_aN > 0.5:
            report_lines.append(f"- → a INCREASES with N")
        else:
            report_lines.append(f"- → a DECREASES with N")

    # ── Pooled fit (all N together) with N as covariate ──
    report_lines.append("\n## Pooled Fit: F7 ~ a·S_BD + b·N + family dummies\n")
    report_lines.append("| metric | a | b(N) | R² | ρ(a_per_N, N) |")
    report_lines.append("|--------|---|------|-----|---------------|")

    for metric_name, metric_fn in BD_METRICS.items():
        y = np.array([r.F7 for r in rows])
        x_bd = np.array([metric_fn(r) for r in rows])
        x_N = np.array([r.N for r in rows], dtype=float)
        fams_all = [r.family for r in rows]

        # Extended design: [intercept, x_bd, N, family dummies]
        unique_fams = sorted(set(fams_all))
        n_obs = len(rows)
        ncols = 3 + len(unique_fams) - 1
        X = np.zeros((n_obs, ncols))
        X[:, 0] = 1.0
        X[:, 1] = x_bd
        X[:, 2] = x_N
        for i, fam in enumerate(fams_all):
            idx = unique_fams.index(fam)
            if idx < len(unique_fams) - 1:
                X[i, 3 + idx] = 1.0

        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        pred = X @ beta
        ss_res = float(np.sum((y - pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

        per_n_results = all_results[metric_name]
        if len(per_n_results) >= 3:
            rho_trend = sp_stats.spearmanr(
                [r["N"] for r in per_n_results],
                [r["a"] for r in per_n_results]
            )[0]
        else:
            rho_trend = float("nan")

        report_lines.append(
            f"| {metric_name} | {float(beta[1]):+.3f} | {float(beta[2]):+.4f} | "
            f"{r2:.4f} | {rho_trend:+.3f} |"
        )

    # ── Summary ──
    report_lines.append("\n## Summary\n")

    best_metric = None
    best_cv = float("inf")
    for metric_name, results in all_results.items():
        if len(results) < 3:
            continue
        a_arr = np.array([r["a"] for r in results])
        a_range = float(np.max(a_arr) - np.min(a_arr))
        a_mean = float(np.mean(a_arr))
        cv = a_range / abs(a_mean) if abs(a_mean) > 1e-12 else float("inf")
        if cv < best_cv:
            best_cv = cv
            best_metric = metric_name

    if best_metric:
        results = all_results[best_metric]
        a_arr = np.array([r["a"] for r in results])
        report_lines.append(f"Most stable metric: **{best_metric}** (CV={best_cv:.3f})")
        report_lines.append(f"a(N) range: [{np.min(a_arr):+.3f}, {np.max(a_arr):+.3f}]")
        report_lines.append(f"a(N) mean ± std: {np.mean(a_arr):+.3f} ± {np.std(a_arr):.3f}")

    report_text = "\n".join(report_lines) + "\n"
    report_path = Path(args.report)
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved: {report_path}")

    # Console summary
    print("\n" + "=" * 60)
    print("LAYER 3: COEFFICIENT STABILITY SUMMARY")
    print("=" * 60)
    for metric_name, results in all_results.items():
        if not results:
            continue
        a_vals = [r["a"] for r in results]
        print(f"\n  {metric_name}:")
        for r in results:
            print(f"    N={r['N']:4d}: a={r['a']:+.3f} ± {r['se_a']:.3f}  R²={r['r2']:.3f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
