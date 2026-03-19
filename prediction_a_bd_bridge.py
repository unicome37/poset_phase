"""Prediction A — Benincasa-Dowker Bridge Experiment

Minimal discriminant test: can adding a BD-inspired sixth term to the
unified functional flip the Lor2D < Lor4D ordering?

F6[X] = F5[X] + α_BD · S_BD(X)

where S_BD is computed from the causal interval distribution of the poset.

The Benincasa-Dowker action for a d-dimensional causal set:
  S_BD^(d) depends on counts C_k = #{pairs (i≺j) with |I(i,j)| = k}
  For d=4:  S_BD = ε N - (2/√6) C_1 + (8/3√6) C_2 - ...
  where ε = ±1 and higher-order terms involve alternating coefficients.

Key physics:
  - 2D sprinklings have dense short intervals → low C_2/C_1 ratio
  - 4D sprinklings have sparser, richer interval structure → higher C_2/C_1
  - S_BD should naturally penalize 2D more than 4D

Execution:
  1. Regenerate posets with same seeds as raw_features.csv
  2. Compute S_BD for each poset
  3. Scan α_BD to find crossing window
  4. Output ranking table per α_BD

Usage:
  python prediction_a_bd_bridge.py [--n_list 16 20 28 36] [--reps 8] [--seed 42]
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np

from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)


# ---------------------------------------------------------------------------
# Benincasa-Dowker action computation
# ---------------------------------------------------------------------------

def count_intervals(poset: Poset) -> dict[int, int]:
    """Count causal intervals: C_k = #{pairs (i≺j) with exactly k elements
    strictly between i and j}.

    Returns dict mapping k -> count.
    """
    c = poset.closure  # c[i,j] = True iff i ≺ j (strict)
    n = poset.n
    counts: dict[int, int] = {}

    for i in range(n):
        for j in range(n):
            if not c[i, j]:
                continue
            # Count elements strictly between i and j
            # k is between i and j iff i≺k and k≺j
            between = 0
            for k in range(n):
                if k == i or k == j:
                    continue
                if c[i, k] and c[k, j]:
                    between += 1
            counts[between] = counts.get(between, 0) + 1

    return counts


def count_intervals_fast(poset: Poset) -> dict[int, int]:
    """Vectorized interval counting using matrix operations."""
    c = poset.closure.astype(np.int32)
    n = poset.n

    # For each pair (i,j) where i≺j, count #{k : i≺k ∧ k≺j}
    # This is (c @ c)[i,j] = Σ_k c[i,k] * c[k,j]
    interval_sizes = c @ c  # interval_sizes[i,j] = number of elements between i and j

    # Mask to only consider pairs where i≺j
    mask = poset.closure

    counts: dict[int, int] = {}
    for i in range(n):
        for j in range(n):
            if mask[i, j]:
                k = int(interval_sizes[i, j])
                counts[k] = counts.get(k, 0) + 1

    return counts


def compute_s_bd_4d(poset: Poset) -> float:
    """Compute the Benincasa-Dowker action for d=4.

    S_BD^(4) = N - (2/√6) C_0 + (8/(3√6)) C_1 - ...

    where C_k = number of pairs with exactly k-element intervals.

    We use the first 3 terms (sufficient for discrimination).
    The action is normalized by N to make it intensive.
    """
    n = poset.n
    if n < 4:
        return 0.0

    counts = count_intervals_fast(poset)

    C_0 = counts.get(0, 0)  # links (no element between)
    C_1 = counts.get(1, 0)  # 1-element intervals
    C_2 = counts.get(2, 0)  # 2-element intervals

    sqrt6 = math.sqrt(6.0)

    # BD action coefficients for d=4 (Benincasa & Dowker 2010)
    # S = ε(N - 2/√6 · C_0 + 8/(3√6) · C_1 - 4/(√6) · C_2 + ...)
    # We take ε = +1 and truncate at C_2
    s_bd = (
        n
        - (2.0 / sqrt6) * C_0
        + (8.0 / (3.0 * sqrt6)) * C_1
        - (4.0 / sqrt6) * C_2
    )

    # Normalize by N to make intensive
    return s_bd / n


def compute_s_bd_ratio(poset: Poset) -> float:
    """Alternative BD-inspired metric: interval richness ratio.

    R = (C_1 + C_2) / (C_0 + 1)

    Higher-dimensional sprinklings produce richer interval structures
    (more elements between causally related pairs).
    Low-d sprinklings are dominated by direct links (C_0).

    Normalized by total number of relations.
    """
    n = poset.n
    if n < 4:
        return 0.0

    counts = count_intervals_fast(poset)
    C_0 = counts.get(0, 0)
    C_1 = counts.get(1, 0)
    C_2 = counts.get(2, 0)
    C_3 = counts.get(3, 0)

    total_relations = sum(counts.values())
    if total_relations == 0:
        return 0.0

    # Fraction of relations that are non-links (have at least 1 element between)
    non_link_frac = 1.0 - C_0 / total_relations

    # Weighted interval depth: average interval size
    weighted_sum = sum(k * v for k, v in counts.items())
    mean_interval = weighted_sum / total_relations if total_relations > 0 else 0.0

    # BD-inspired: combine non-link fraction with mean interval depth
    # This should be larger for higher-d (richer interval structure)
    # and smaller for lower-d (dominated by direct links)
    return non_link_frac * (1.0 + mean_interval)


# ---------------------------------------------------------------------------
# F5 computation from raw_features.csv
# ---------------------------------------------------------------------------

CALIBRATED_WEIGHTS = {
    "beta": 2.0,
    "gamma": 0.5,
    "lam": 1.5,
    "eta": 0.1,
    "kappa": 0.05,
}

BAYESIAN_WEIGHTS = {
    "beta": 1.0,
    "gamma": 0.0004,
    "lam": 0.888,
    "eta": 0.637,
    "kappa": 0.068,
}


def compute_f5(row: dict, weights: dict | None = None) -> float:
    """Compute F5 from a raw_features.csv row."""
    if weights is None:
        weights = CALIBRATED_WEIGHTS
    return (
        weights["beta"] * float(row["log_H"])
        + weights["gamma"] * float(row["pi_geo"])
        - weights["lam"] * float(row["sigma_hist"])
        + weights["eta"] * float(row["xi_dim"])
        + weights["kappa"] * float(row["pi_cg"])
    )


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def load_raw_features(path: str) -> list[dict]:
    """Load raw_features.csv."""
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def regenerate_poset(family: str, n: int, rep: int, base_seed: int = 42) -> Poset:
    """Regenerate a poset with the same seed as used in raw_features extraction."""
    generators = {
        "Lor2D": generate_lorentzian_like_2d,
        "Lor3D": generate_lorentzian_like_3d,
        "Lor4D": generate_lorentzian_like_4d,
        "Lor5D": generate_lorentzian_like_5d,
    }
    gen = generators[family]
    seed = base_seed + rep * 1000 + n * 100
    return gen(n, seed=seed)


def run_experiment(
    raw_features_path: str,
    alpha_values: list[float],
    weight_set: str = "calibrated",
    base_seed: int = 42,
) -> list[dict]:
    """Run the α_BD crossing experiment.

    For each α_BD value, compute F6 = F5 + α_BD · S_BD for all 4 Lorentzian
    families, and check ranking.
    """
    rows = load_raw_features(raw_features_path)

    # Filter to Lorentzian families only
    lor_families = {"Lor2D", "Lor3D", "Lor4D", "Lor5D"}
    rows = [r for r in rows if r["family"] in lor_families]

    weights = CALIBRATED_WEIGHTS if weight_set == "calibrated" else BAYESIAN_WEIGHTS

    print(f"Weight set: {weight_set}")
    print(f"Weights: {weights}")
    print(f"Families: {lor_families}")
    print(f"Rows: {len(rows)}")
    print()

    # --- Step 1: Compute S_BD for each poset ---
    print("Step 1: Computing S_BD for all posets...")

    bd_data = []  # list of (row_dict, s_bd_4d, s_bd_ratio, f5)
    for i, r in enumerate(rows):
        family = r["family"]
        n = int(r["N"])
        rep = int(r["rep"])

        poset = regenerate_poset(family, n, rep, base_seed)
        s_bd_4d = compute_s_bd_4d(poset)
        s_bd_ratio = compute_s_bd_ratio(poset)
        f5 = compute_f5(r, weights)

        bd_data.append({
            "family": family,
            "N": n,
            "rep": rep,
            "f5": f5,
            "s_bd_4d": s_bd_4d,
            "s_bd_ratio": s_bd_ratio,
            "comp_frac": float(r["comp_frac"]),
        })

        if (i + 1) % 16 == 0:
            print(f"  [{i+1}/{len(rows)}] computed")

    print(f"  Total: {len(bd_data)} posets\n")

    # --- Step 2: Report S_BD statistics by family ---
    print("=" * 70)
    print("S_BD STATISTICS BY FAMILY")
    print("=" * 70)

    families_order = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]
    for metric_name, metric_key in [("S_BD(4d)", "s_bd_4d"), ("S_BD(ratio)", "s_bd_ratio")]:
        print(f"\n{metric_name}:")
        print(f"  {'Family':<8s} {'Mean':>10s} {'Std':>10s} {'Min':>10s} {'Max':>10s}")
        print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
        for fam in families_order:
            vals = [d[metric_key] for d in bd_data if d["family"] == fam]
            if vals:
                print(f"  {fam:<8s} {np.mean(vals):>10.4f} {np.std(vals):>10.4f} "
                      f"{np.min(vals):>10.4f} {np.max(vals):>10.4f}")

    # --- Step 3: α_BD crossing scan ---
    print("\n" + "=" * 70)
    print("α_BD CROSSING SCAN")
    print("=" * 70)

    # From the data:
    #   S_BD_ratio: HIGHER for 2D (~1.3-3.1), LOWER for 4D/5D (~0.07-0.18)
    #   S_BD_4d: varies with N, sign depends on BD coefficient balance
    #
    # To penalize 2D (push it up), we want:
    #   +α * S_BD_ratio  (penalize rich interval structure = low-d signature)
    #   -α * S_BD_4d     (if S_BD_4d is lower for 2D)
    #
    # We test all sign combinations to find any crossing.

    results = []

    for alpha_bd in alpha_values:
        for metric_name, metric_key, sign in [
            ("S_BD_4d(+)", "s_bd_4d", +1),
            ("S_BD_4d(-)", "s_bd_4d", -1),
            ("S_BD_ratio(+)", "s_bd_ratio", +1),  # penalize high ratio (= low-d)
            ("S_BD_ratio(-)", "s_bd_ratio", -1),  # reward high ratio
        ]:
            # F6 = F5 + sign * alpha_bd * metric
            family_means = {}
            for fam in families_order:
                f6_vals = [
                    d["f5"] + sign * alpha_bd * d[metric_key]
                    for d in bd_data if d["family"] == fam
                ]
                family_means[fam] = np.mean(f6_vals)

            # Check three conditions
            cond1 = family_means["Lor4D"] < family_means["Lor5D"]  # preserve 4D < 5D
            cond2 = family_means["Lor4D"] < family_means["Lor2D"]  # NEW: 4D < 2D
            cond3 = family_means["Lor4D"] < family_means["Lor3D"]  # bonus: 4D < 3D

            # Ranking
            ranking = sorted(families_order, key=lambda f: family_means[f])

            results.append({
                "alpha_BD": alpha_bd,
                "metric": metric_name,
                "sign": sign,
                "F6_Lor2D": family_means["Lor2D"],
                "F6_Lor3D": family_means["Lor3D"],
                "F6_Lor4D": family_means["Lor4D"],
                "F6_Lor5D": family_means["Lor5D"],
                "cond1_4D<5D": cond1,
                "cond2_4D<2D": cond2,
                "cond3_4D<3D": cond3,
                "all_three": cond1 and cond2 and cond3,
                "ranking": " < ".join(ranking),
            })

    # --- Step 4: Print scan table ---
    for metric_name in ["S_BD_4d(+)", "S_BD_4d(-)", "S_BD_ratio(+)", "S_BD_ratio(-)"]:
        metric_results = [r for r in results if r["metric"] == metric_name]
        if not metric_results:
            continue
        sign = metric_results[0]["sign"]
        sign_str = "+" if sign > 0 else "-"

        print(f"\n--- F6 = F5 {sign_str} α·{metric_name} ---")
        print(f"  {'α_BD':>8s} | {'F6(2D)':>10s} {'F6(3D)':>10s} {'F6(4D)':>10s} "
              f"{'F6(5D)':>10s} | {'4D<5D':>5s} {'4D<2D':>5s} {'4D<3D':>5s} | Ranking")
        print(f"  {'-'*8} | {'-'*10} {'-'*10} {'-'*10} {'-'*10} | {'-'*5} {'-'*5} {'-'*5} | {'-'*30}")

        for r in metric_results:
            c1 = "✓" if r["cond1_4D<5D"] else "✗"
            c2 = "✓" if r["cond2_4D<2D"] else "✗"
            c3 = "✓" if r["cond3_4D<3D"] else "✗"
            marker = " ★" if r["all_three"] else ""
            print(f"  {r['alpha_BD']:8.2f} | {r['F6_Lor2D']:10.4f} {r['F6_Lor3D']:10.4f} "
                  f"{r['F6_Lor4D']:10.4f} {r['F6_Lor5D']:10.4f} | "
                  f"{c1:>5s} {c2:>5s} {c3:>5s} | {r['ranking']}{marker}")

    # --- Step 5: Find crossing windows ---
    print("\n" + "=" * 70)
    print("CROSSING WINDOW ANALYSIS")
    print("=" * 70)

    for metric_name in ["S_BD_4d(+)", "S_BD_4d(-)", "S_BD_ratio(+)", "S_BD_ratio(-)"]:
        metric_results = [r for r in results if r["metric"] == metric_name]
        if not metric_results:
            continue
        windows = [r for r in metric_results if r["all_three"]]

        if windows:
            alpha_min = min(r["alpha_BD"] for r in windows)
            alpha_max = max(r["alpha_BD"] for r in windows)
            print(f"\n  {metric_name}: ★ WINDOW FOUND ★")
            print(f"    α_BD ∈ [{alpha_min:.2f}, {alpha_max:.2f}]")
            print(f"    {len(windows)}/{len(metric_results)} scan points satisfy all 3 conditions")

            # Show best point (4D most below 2D)
            best = min(windows, key=lambda r: r["F6_Lor4D"] - r["F6_Lor2D"])
            gap = best["F6_Lor2D"] - best["F6_Lor4D"]
            print(f"    Best α_BD = {best['alpha_BD']:.2f}: F6(2D)-F6(4D) = {gap:.4f}")
            print(f"    Ranking: {best['ranking']}")
        else:
            print(f"\n  {metric_name}: No crossing window found in scanned range.")
            # Report closest miss
            partial = [r for r in metric_results if r["cond1_4D<5D"] and r["cond2_4D<2D"]]
            if partial:
                print(f"    But {len(partial)} points satisfy cond1+cond2 (4D<5D AND 4D<2D)")
                best = min(partial, key=lambda r: r["F6_Lor4D"] - r["F6_Lor3D"])
                gap = best["F6_Lor3D"] - best["F6_Lor4D"]
                print(f"    Closest to cond3: α_BD={best['alpha_BD']:.2f}, F6(3D)-F6(4D)={gap:.4f}")

    # --- Step 6: Per-N breakdown for best metric ---
    print("\n" + "=" * 70)
    print("PER-N BREAKDOWN (at selected α_BD values)")
    print("=" * 70)

    # Find best α_BD from either metric
    all_windows = [r for r in results if r["all_three"]]
    if all_windows:
        best_result = min(all_windows, key=lambda r: r["F6_Lor4D"] - r["F6_Lor2D"])
        best_alpha = best_result["alpha_BD"]
        best_metric_key = "s_bd_4d" if best_result["metric"] == "S_BD_4d" else "s_bd_ratio"
        best_sign = best_result["sign"]
    else:
        # Use moderate α for diagnostic
        best_alpha = alpha_values[len(alpha_values) // 2]
        best_metric_key = "s_bd_ratio"
        best_sign = -1

    print(f"\n  Using α_BD = {best_alpha:.2f}, metric = {best_metric_key}, sign = {best_sign:+d}")

    n_values = sorted(set(d["N"] for d in bd_data))
    for n_val in n_values:
        print(f"\n  N = {n_val}:")
        print(f"    {'Family':<8s} {'F5':>10s} {'S_BD':>10s} {'F6':>10s}")
        print(f"    {'-'*8} {'-'*10} {'-'*10} {'-'*10}")

        n_data = [d for d in bd_data if d["N"] == n_val]
        for fam in families_order:
            fam_data = [d for d in n_data if d["family"] == fam]
            if fam_data:
                f5_mean = np.mean([d["f5"] for d in fam_data])
                sbd_mean = np.mean([d[best_metric_key] for d in fam_data])
                f6_mean = f5_mean + best_sign * best_alpha * sbd_mean
                print(f"    {fam:<8s} {f5_mean:>10.4f} {sbd_mean:>10.4f} {f6_mean:>10.4f}")

        # Ranking for this N
        fam_f6 = {}
        for fam in families_order:
            fam_data = [d for d in n_data if d["family"] == fam]
            if fam_data:
                f5_mean = np.mean([d["f5"] for d in fam_data])
                sbd_mean = np.mean([d[best_metric_key] for d in fam_data])
                fam_f6[fam] = f5_mean + best_sign * best_alpha * sbd_mean
        ranking = sorted(fam_f6.keys(), key=lambda f: fam_f6[f])
        print(f"    Ranking: {' < '.join(ranking)}")

    # --- Step 7: Faithfulness diagnostic (external, not in functional) ---
    print("\n" + "=" * 70)
    print("FAITHFULNESS DIAGNOSTIC (external indicator)")
    print("=" * 70)
    print("\nInterval distribution statistics:")
    print(f"  {'Family':<8s} {'C_0/total':>10s} {'C_1/total':>10s} {'C_2/total':>10s} "
          f"{'C_3+/total':>10s} {'mean_interval':>14s}")
    print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*14}")

    for fam in families_order:
        fam_data_all = [d for d in bd_data if d["family"] == fam]
        # Need to recompute interval distributions
        c0_fracs, c1_fracs, c2_fracs, c3p_fracs, mean_ints = [], [], [], [], []
        for d in fam_data_all:
            poset = regenerate_poset(d["family"], d["N"], d["rep"], base_seed)
            counts = count_intervals_fast(poset)
            total = sum(counts.values())
            if total == 0:
                continue
            c0_fracs.append(counts.get(0, 0) / total)
            c1_fracs.append(counts.get(1, 0) / total)
            c2_fracs.append(counts.get(2, 0) / total)
            c3p = sum(v for k, v in counts.items() if k >= 3) / total
            c3p_fracs.append(c3p)
            mean_int = sum(k * v for k, v in counts.items()) / total
            mean_ints.append(mean_int)

        print(f"  {fam:<8s} {np.mean(c0_fracs):>10.4f} {np.mean(c1_fracs):>10.4f} "
              f"{np.mean(c2_fracs):>10.4f} {np.mean(c3p_fracs):>10.4f} "
              f"{np.mean(mean_ints):>14.4f}")

    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prediction A — Benincasa-Dowker Bridge Experiment"
    )
    parser.add_argument("--raw_features", type=str,
                        default="outputs_unified_functional/raw_features.csv",
                        help="Path to raw_features.csv")
    parser.add_argument("--weights", type=str, default="calibrated",
                        choices=["calibrated", "bayesian"],
                        help="Weight set to use for F5")
    parser.add_argument("--alpha_min", type=float, default=0.0,
                        help="Minimum α_BD to scan")
    parser.add_argument("--alpha_max", type=float, default=20.0,
                        help="Maximum α_BD to scan")
    parser.add_argument("--alpha_steps", type=int, default=21,
                        help="Number of α_BD values to scan")
    parser.add_argument("--seed", type=int, default=42,
                        help="Base random seed (must match raw_features generation)")
    args = parser.parse_args()

    alpha_values = np.linspace(args.alpha_min, args.alpha_max, args.alpha_steps).tolist()

    print("=" * 70)
    print("PREDICTION A — BENINCASA-DOWKER BRIDGE EXPERIMENT")
    print("=" * 70)
    print(f"  α_BD range: [{args.alpha_min}, {args.alpha_max}], {args.alpha_steps} points")
    print(f"  Weight set: {args.weights}")
    print(f"  Seed: {args.seed}")
    print()

    results = run_experiment(
        raw_features_path=args.raw_features,
        alpha_values=alpha_values,
        weight_set=args.weights,
        base_seed=args.seed,
    )

    # Save results
    outdir = Path("outputs_unified_functional")
    outdir.mkdir(parents=True, exist_ok=True)

    outpath = str(outdir / "prediction_a_bd_scan.csv")
    with open(outpath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "alpha_BD", "metric", "sign",
            "F6_Lor2D", "F6_Lor3D", "F6_Lor4D", "F6_Lor5D",
            "cond1_4D<5D", "cond2_4D<2D", "cond3_4D<3D", "all_three",
            "ranking",
        ])
        writer.writeheader()
        writer.writerows(results)
    print(f"\n→ Saved scan results to {outpath}")


if __name__ == "__main__":
    main()
