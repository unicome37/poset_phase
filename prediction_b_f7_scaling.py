"""Prediction B — F7 N-scaling test: does F7(Lor) < F7(KR) hold at large N?

Extends the core Prediction B test to N=20–72 using the definitive F7 model
(§5.10.7) with sigmoid wall.  Also tests Lor4D vs KR (the key physical pair).

F7 = logH + 0.0004·Π_geo − 10·Σ_hist + 0.6·Ξ_d + α(N)·σ((R−Rc)/w)

Output:
  - outputs_unified_functional/prediction_b_f7_scaling.csv   (raw per-rep data)
  - outputs_unified_functional/prediction_b_f7_scaling.md    (summary report)
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from bd_action import count_intervals_fast
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


# ── F7 computation (copied from conjecture_e_f7_bridge.py) ──────────────

def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_R(poset: Poset) -> float:
    """Interval occupancy ratio R = 1 - C0/total."""
    counts = count_intervals_fast(poset, k_max=3)
    if counts.total_relations <= 0:
        return 0.0
    return 1.0 - float(counts.get(0)) / float(counts.total_relations)


def compute_F7(
    poset: Poset,
    N: int,
    *,
    alpha0: float = 16.0,
    q: float = -0.5,
    lam: float = 10.0,
    eta: float = 0.6,
    Rc: float = 0.25,
    w: float = 0.015,
    N0: float = 20.0,
    sis_runs: int = 512,
) -> dict:
    """Compute F7 and all components for a single poset."""
    log_H = compute_log_H(poset, n_runs=sis_runs)
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
    }


# ── Generators ───────────────────────────────────────────────────────────

FAMILIES = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}

# Key pairwise tests for Prediction B
PAIRS = [
    ("Lor2D", "KR_like"),
    ("Lor4D", "KR_like"),
    ("Lor3D", "KR_like"),
]


# ── Main experiment ──────────────────────────────────────────────────────

def run_scaling(
    n_values: list[int],
    families: list[str],
    reps: int,
    base_seed: int,
    sis_runs: int,
) -> list[dict]:
    rows = []
    total = len(n_values) * len(families) * reps
    count = 0

    for N in n_values:
        for fam in families:
            gen = FAMILIES[fam]
            for rep in range(reps):
                count += 1
                seed = base_seed + rep * 1000 + N * 100
                poset = gen(N, seed=seed)
                comp = compute_F7(poset, N, sis_runs=sis_runs)
                row = {
                    "family": fam,
                    "N": N,
                    "rep": rep,
                    "seed": seed,
                    **comp,
                }
                rows.append(row)

                if count % 50 == 0 or count == total:
                    print(f"  [{count}/{total}] {fam} N={N} rep={rep}: "
                          f"F7={comp['F7']:.2f} wall={comp['wall']:.2f} R={comp['R']:.3f}")

    return rows


def analyze(rows: list[dict], pairs: list[tuple[str, str]]) -> str:
    """Generate summary markdown report."""
    lines = ["# Prediction B — F7 N-Scaling Report\n"]
    lines.append(f"Total data points: {len(rows)}\n")

    # Group by (N, family)
    from collections import defaultdict
    by_nf: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for r in rows:
        by_nf[(r["N"], r["family"])].append(r)

    n_values = sorted(set(r["N"] for r in rows))

    # Family-level means table
    lines.append("## Family-level means\n")
    lines.append("| N | Family | mean F7 | mean logH | mean wall | mean R | mean d_eff |")
    lines.append("|---|--------|---------|-----------|-----------|--------|------------|")
    for N in n_values:
        for fam in sorted(set(r["family"] for r in rows)):
            key = (N, fam)
            if key not in by_nf:
                continue
            grp = by_nf[key]
            mF = np.mean([r["F7"] for r in grp])
            mH = np.mean([r["log_H"] for r in grp])
            mW = np.mean([r["wall"] for r in grp])
            mR = np.mean([r["R"] for r in grp])
            mD = np.mean([r["d_eff"] for r in grp])
            lines.append(f"| {N} | {fam} | {mF:.2f} | {mH:.2f} | {mW:.2f} | {mR:.3f} | {mD:.2f} |")

    # Pairwise win-rate tables
    lines.append("\n## Pairwise Prediction B tests\n")

    for left, right in pairs:
        lines.append(f"### {left} vs {right}\n")
        lines.append("| N | left_wins | total | win% | mean_ΔF7 | mean_Δwall | mean_ΔR |")
        lines.append("|---|-----------|-------|------|----------|------------|---------|")

        all_pass = True
        for N in n_values:
            left_rows = by_nf.get((N, left), [])
            right_rows = by_nf.get((N, right), [])
            if not left_rows or not right_rows:
                continue

            wins = 0
            delta_F7s = []
            delta_walls = []
            delta_Rs = []
            n_pairs = min(len(left_rows), len(right_rows))
            for i in range(n_pairs):
                dF = left_rows[i]["F7"] - right_rows[i]["F7"]
                dW = left_rows[i]["wall"] - right_rows[i]["wall"]
                dR = left_rows[i]["R"] - right_rows[i]["R"]
                delta_F7s.append(dF)
                delta_walls.append(dW)
                delta_Rs.append(dR)
                if dF < 0:
                    wins += 1

            pct = 100 * wins / n_pairs if n_pairs > 0 else 0
            if pct < 50:
                all_pass = False
            mdf = np.mean(delta_F7s)
            mdw = np.mean(delta_walls)
            mdr = np.mean(delta_Rs)
            lines.append(f"| {N} | {wins} | {n_pairs} | {pct:.0f}% | {mdf:+.2f} | {mdw:+.2f} | {mdr:+.3f} |")

        verdict = "**PASS**" if all_pass else "**PARTIAL**"
        lines.append(f"\nVerdict: {verdict}\n")

    # Summary
    lines.append("\n## Summary\n")
    for left, right in pairs:
        all_wins = []
        for N in n_values:
            left_rows = by_nf.get((N, left), [])
            right_rows = by_nf.get((N, right), [])
            n_pairs = min(len(left_rows), len(right_rows))
            for i in range(n_pairs):
                all_wins.append(left_rows[i]["F7"] < right_rows[i]["F7"])
        if all_wins:
            total_win = sum(all_wins)
            pct = 100 * total_win / len(all_wins)
            lines.append(f"- **{left} < {right}**: {total_win}/{len(all_wins)} ({pct:.1f}%)")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Prediction B — F7 N-scaling")
    parser.add_argument("--quick", action="store_true", help="Quick test with fewer N/reps")
    parser.add_argument("--reps", type=int, default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sis-runs", type=int, default=512)
    args = parser.parse_args()

    if args.quick:
        n_values = [20, 36, 48]
        families = ["Lor2D", "Lor4D", "KR_like"]
        reps = args.reps or 5
    else:
        n_values = [20, 28, 36, 44, 52, 56, 60, 64, 72]
        families = ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]
        reps = args.reps or 10

    print(f"=== Prediction B — F7 N-Scaling ===")
    print(f"N values: {n_values}")
    print(f"Families: {families}")
    print(f"Reps: {reps}, SIS runs: {args.sis_runs}")
    print()

    rows = run_scaling(n_values, families, reps, args.seed, args.sis_runs)

    # Save raw CSV
    outdir = Path("outputs_unified_functional")
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = outdir / "prediction_b_f7_scaling.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nSaved: {csv_path}")

    # Determine which pairs are testable
    fam_set = set(r["family"] for r in rows)
    active_pairs = [(l, r) for l, r in PAIRS if l in fam_set and r in fam_set]

    # Generate report
    report = analyze(rows, active_pairs)
    md_path = outdir / "prediction_b_f7_scaling.md"
    md_path.write_text(report, encoding="utf-8")
    print(f"Saved: {md_path}")

    # Print summary to console
    print("\n" + "=" * 70)
    print(report)


if __name__ == "__main__":
    main()
