"""Curvature Calibration Analysis — Layer 2 Follow-up.

Computes the calibration ratio c(d, N) = R_hat / R_dS for de Sitter sprinklings
and tests whether c stabilizes as N → ∞.

For d-dimensional de Sitter with Hubble parameter H:
  R_dS = d(d-1) H²   (Ricci scalar)
  d=2: R=2H², d=3: R=6H², d=4: R=12H²

We exclude H=0 (R_dS=0) and compute per-(d, N, H) the mean ratio.
Then track c(d, N) across N to test convergence.

Uses existing CSV output from Layer 2 experiments.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats


def load_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main() -> int:
    base = Path("outputs_unified_functional")
    files = [
        base / "curvature_layer2_recovery_full.csv",
        base / "curvature_layer2_recovery_n2048.csv",
    ]

    all_rows = []
    for fp in files:
        if fp.exists():
            rows = load_csv(fp)
            all_rows.extend(rows)
            print(f"Loaded {len(rows)} rows from {fp}")
        else:
            print(f"WARN: {fp} not found")

    if not all_rows:
        print("No data found!")
        return 1

    # Parse and filter H > 0
    data = []
    for r in all_rows:
        d = int(r["d"])
        N = int(r["N"])
        H = float(r["hubble"])
        R_hat = float(r["R_hat"])
        R_dS = float(r["R_dS"])
        n_pairs = int(r["n_causal_pairs"])

        if H <= 0 or R_dS == 0:
            continue

        ratio = R_hat / R_dS
        data.append({
            "d": d, "N": N, "H": H, "R_hat": R_hat,
            "R_dS": R_dS, "ratio": ratio, "n_pairs": n_pairs,
        })

    print(f"\nTotal non-flat entries: {len(data)}")

    # ── Per-(d, N) calibration ratio ──
    # Group by (d, N)
    groups: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for r in data:
        groups[(r["d"], r["N"])].append(r)

    lines: list[str] = []
    lines.append("# Curvature Calibration: R_hat / R_dS Convergence\n")
    lines.append(f"Total entries (H>0): {len(data)}\n")

    # ── Table 1: Per-(d, N) mean ratio ──
    lines.append("\n## 1. Mean Calibration Ratio c(d, N) = mean(R_hat / R_dS)\n")
    lines.append("| d | N | mean(ratio) | std(ratio) | median | n | mean(n_pairs) |")
    lines.append("|---|---|-------------|------------|--------|---|---------------|")

    ds = sorted(set(r["d"] for r in data))
    ns = sorted(set(r["N"] for r in data))

    ratio_table: dict[tuple[int, int], dict] = {}
    for d in ds:
        for N in ns:
            key = (d, N)
            if key not in groups:
                continue
            g = groups[key]
            ratios = np.array([r["ratio"] for r in g])
            pairs = np.array([r["n_pairs"] for r in g])
            result = {
                "mean": float(np.mean(ratios)),
                "std": float(np.std(ratios)),
                "median": float(np.median(ratios)),
                "n": len(g),
                "mean_pairs": float(np.mean(pairs)),
            }
            ratio_table[key] = result
            lines.append(
                f"| {d} | {N} | {result['mean']:+.2f} | {result['std']:.2f} | "
                f"{result['median']:+.2f} | {result['n']} | {result['mean_pairs']:.0f} |"
            )

    # ── Table 2: Per-(d, N, H) breakdown ──
    lines.append("\n## 2. Per-(d, N, H) Ratio Breakdown\n")
    lines.append("| d | N | H | mean(ratio) | std | n | mean(n_pairs) |")
    lines.append("|---|---|---|-------------|-----|---|---------------|")

    hs = sorted(set(r["H"] for r in data))
    for d in ds:
        for N in ns:
            for H in hs:
                subset = [r for r in data if r["d"] == d and r["N"] == N and r["H"] == H]
                if not subset:
                    continue
                ratios = np.array([r["ratio"] for r in subset])
                pairs = np.array([r["n_pairs"] for r in subset])
                lines.append(
                    f"| {d} | {N} | {H} | {np.mean(ratios):+.2f} | "
                    f"{np.std(ratios):.2f} | {len(subset)} | {np.mean(pairs):.0f} |"
                )

    # ── Table 3: N-convergence of c(d) ──
    lines.append("\n## 3. N-Convergence of c(d)\n")
    lines.append("Does c(d, N) stabilize as N → ∞?\n")

    for d in ds:
        d_ns = sorted([N for N in ns if (d, N) in ratio_table])
        if len(d_ns) < 2:
            continue

        c_vals = [ratio_table[(d, N)]["mean"] for N in d_ns]
        c_stds = [ratio_table[(d, N)]["std"] for N in d_ns]

        lines.append(f"\n### d={d}\n")
        lines.append(f"c(N) trajectory: " +
                     " → ".join(f"{c:+.2f}" for c in c_vals))
        lines.append(f"std trajectory:  " +
                     " → ".join(f"{s:.2f}" for s in c_stds))

        # Range and coefficient of variation
        c_arr = np.array(c_vals)
        c_range = float(np.max(c_arr) - np.min(c_arr))
        c_mean = float(np.mean(c_arr))
        cv = c_range / abs(c_mean) if abs(c_mean) > 1e-12 else float("inf")
        lines.append(f"- Range: {c_range:.2f}, Mean: {c_mean:+.2f}, CV: {cv:.3f}")

        # Spearman trend
        if len(d_ns) >= 3:
            rho, p = sp_stats.spearmanr(d_ns, c_vals)
            lines.append(f"- Spearman(N, c): ρ={rho:+.3f} (p={p:.3e})")

        # std shrinking?
        s_arr = np.array(c_stds)
        if len(s_arr) >= 3:
            rho_s, p_s = sp_stats.spearmanr(d_ns, c_stds)
            lines.append(f"- Spearman(N, std): ρ={rho_s:+.3f} (p={p_s:.3e})")
            if rho_s < -0.5:
                lines.append(f"- → **std shrinking** (calibration becoming more precise)")

    # ── Table 4: H-dependence check ──
    lines.append("\n## 4. H-Independence Check\n")
    lines.append("If R_hat ∝ R_dS, then ratio should be H-independent.\n")

    for d in ds:
        for N in ns:
            subset = [r for r in data if r["d"] == d and r["N"] == N]
            if len(subset) < 5:
                continue
            h_vals = np.array([r["H"] for r in subset])
            r_vals = np.array([r["ratio"] for r in subset])
            if len(set(h_vals)) < 2:
                continue
            rho, p = sp_stats.spearmanr(h_vals, r_vals)
            lines.append(f"- d={d}, N={N}: ρ(H, ratio)={rho:+.3f} (p={p:.3e})")

    # ── Summary ──
    lines.append("\n## Summary\n")
    for d in ds:
        d_ns = sorted([N for N in ns if (d, N) in ratio_table])
        if d_ns:
            last_c = ratio_table[(d, d_ns[-1])]["mean"]
            last_s = ratio_table[(d, d_ns[-1])]["std"]
            lines.append(f"- **d={d}**: c(N_max={d_ns[-1]}) = {last_c:+.2f} ± {last_s:.2f}")

    report_text = "\n".join(lines) + "\n"
    report_path = base / "curvature_calibration_analysis.md"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"\nSaved: {report_path}")

    # Console summary
    print("\n" + "=" * 60)
    print("CALIBRATION RATIO c(d, N) = R_hat / R_dS")
    print("=" * 60)
    for d in ds:
        d_ns = sorted([N for N in ns if (d, N) in ratio_table])
        print(f"\n  d={d}:")
        for N in d_ns:
            r = ratio_table[(d, N)]
            print(f"    N={N:5d}: c = {r['mean']:+.2f} ± {r['std']:.2f}"
                  f"  (median={r['median']:+.2f}, n={r['n']})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
