"""Curvature Calibration — Power-Law Fitting.

The naive ratio c = R_hat / R_dS is H-dependent (grows with H),
indicating R_hat is NOT simply proportional to R_dS.

This script fits the actual relationship:
  R_hat = c(d,N) · H^α(d,N)

If α → 2.0 as N → ∞, then R_hat ∝ H² ∝ R_dS (correct scaling).
If α ≠ 2, the discrete curvature proxy has anomalous scaling.

Also fits: R_hat = a · R_dS^β (power law in R_dS directly).

Uses log-log regression on H > 0 data.
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

    # Parse — keep H > 0 only
    data = []
    for r in all_rows:
        d = int(r["d"])
        N = int(r["N"])
        H = float(r["hubble"])
        R_hat = float(r["R_hat"])
        R_dS = float(r["R_dS"])
        n_pairs = int(r["n_causal_pairs"])
        if H <= 0:
            continue
        data.append({
            "d": d, "N": N, "H": H, "R_hat": R_hat,
            "R_dS": R_dS, "n_pairs": n_pairs,
        })

    # Also get H=0 bias for subtraction
    flat_data = []
    for r in all_rows:
        d = int(r["d"])
        N = int(r["N"])
        H = float(r["hubble"])
        R_hat = float(r["R_hat"])
        if H == 0:
            flat_data.append({"d": d, "N": N, "R_hat_flat": R_hat})

    # Compute mean flat bias per (d, N)
    flat_bias: dict[tuple[int, int], float] = {}
    groups_flat: dict[tuple[int, int], list[float]] = defaultdict(list)
    for r in flat_data:
        groups_flat[(r["d"], r["N"])].append(r["R_hat_flat"])
    for key, vals in groups_flat.items():
        flat_bias[key] = float(np.mean(vals))

    print(f"\nFlat bias per (d, N):")
    for (d, N), bias in sorted(flat_bias.items()):
        print(f"  d={d}, N={N}: R_hat(H=0) = {bias:+.2f}")

    ds = sorted(set(r["d"] for r in data))
    ns = sorted(set(r["N"] for r in data))
    hs = sorted(set(r["H"] for r in data))

    lines: list[str] = []
    lines.append("# Curvature Calibration: Power-Law Fitting\n")
    lines.append("## Motivation\n")
    lines.append("The naive ratio c = R_hat/R_dS grows with H, indicating R_hat is NOT")
    lines.append("simply proportional to H². We fit R_hat = c·H^α and check if α→2.\n")

    # ── Strategy A: Fit R_hat(H) = c · H^α using median per (d,N,H) ──
    lines.append("\n## Strategy A: Fit mean(R_hat) vs H per (d, N)\n")
    lines.append("Using bias-subtracted R_hat: R_hat_corr = R_hat - R_hat(H=0)\n")
    lines.append("| d | N | α (exponent) | log(c) | R² | α_raw (no bias sub) |")
    lines.append("|---|---|-------------|--------|-----|---------------------|")

    alpha_table: dict[tuple[int, int], dict] = {}

    for d in ds:
        for N in ns:
            subset = [r for r in data if r["d"] == d and r["N"] == N]
            if len(subset) < 5:
                continue

            # Compute mean R_hat per H
            h_means: dict[float, float] = {}
            h_means_raw: dict[float, float] = {}
            for H in hs:
                h_sub = [r for r in subset if r["H"] == H]
                if h_sub:
                    mean_rhat = float(np.mean([r["R_hat"] for r in h_sub]))
                    h_means_raw[H] = mean_rhat
                    # Bias-subtracted
                    bias = flat_bias.get((d, N), 0.0)
                    h_means[H] = mean_rhat - bias

            # Log-log fit: log(R_hat_corr) = α·log(H) + log(c)
            # Only use H values where R_hat_corr > 0
            log_H_arr = []
            log_R_arr = []
            log_R_raw_arr = []
            for H in sorted(h_means.keys()):
                if h_means[H] > 0:
                    log_H_arr.append(math.log(H))
                    log_R_arr.append(math.log(h_means[H]))
                if h_means_raw.get(H, 0) > 0:
                    log_R_raw_arr.append(math.log(h_means_raw[H]))

            if len(log_H_arr) < 2:
                continue

            log_H = np.array(log_H_arr)
            log_R = np.array(log_R_arr)

            slope, intercept, r_val, p_val, se = sp_stats.linregress(log_H, log_R)

            # Raw (no bias subtraction)
            alpha_raw = float("nan")
            if len(log_R_raw_arr) == len(log_H_arr):
                log_R_raw = np.array(log_R_raw_arr)
                slope_raw, _, _, _, _ = sp_stats.linregress(log_H, log_R_raw)
                alpha_raw = slope_raw

            result = {
                "alpha": float(slope),
                "log_c": float(intercept),
                "r2": float(r_val ** 2),
                "alpha_raw": alpha_raw,
                "se": float(se),
            }
            alpha_table[(d, N)] = result

            lines.append(
                f"| {d} | {N} | **{result['alpha']:.3f}** ± {result['se']:.3f} | "
                f"{result['log_c']:.2f} | {result['r2']:.4f} | {alpha_raw:.3f} |"
            )

    # ── N-convergence of α(d) ──
    lines.append("\n## α(d, N) Convergence\n")
    lines.append("Does α → 2.0 (linear in R_dS) as N → ∞?\n")

    for d in ds:
        d_ns = sorted([N for N in ns if (d, N) in alpha_table])
        if len(d_ns) < 2:
            continue

        alphas = [alpha_table[(d, N)]["alpha"] for N in d_ns]
        lines.append(f"\n### d={d}\n")
        lines.append(f"α(N): " + " → ".join(f"{a:.3f}" for a in alphas))

        a_arr = np.array(alphas)
        lines.append(f"- Mean α: {np.mean(a_arr):.3f} ± {np.std(a_arr):.3f}")
        lines.append(f"- Target: α = 2.0 (proportional to H²)")
        lines.append(f"- Deviation from 2: {np.mean(a_arr) - 2.0:+.3f}")

        if len(d_ns) >= 3:
            rho, p = sp_stats.spearmanr(d_ns, alphas)
            lines.append(f"- Spearman(N, α): ρ={rho:+.3f} (p={p:.3e})")
            if rho > 0.5 and np.mean(a_arr) < 2.0:
                lines.append(f"- → α **increasing toward 2** (converging)")
            elif abs(rho) < 0.5:
                lines.append(f"- → α stable (no N trend)")

    # ── Strategy B: slope dR_hat/dH² — already from Layer 2 ──
    lines.append("\n## Strategy B: Linear Slope in H² (recap from Layer 2)\n")
    lines.append("If R_hat = c·H² + bias, the slope is the calibration constant.\n")

    for d in ds:
        for N in ns:
            subset = [r for r in data if r["d"] == d and r["N"] == N]
            if len(subset) < 5:
                continue
            H_arr = np.array([r["H"] for r in subset])
            R_arr = np.array([r["R_hat"] for r in subset])
            H2_arr = H_arr ** 2

            slope, intercept, r_val, p_val, se = sp_stats.linregress(H2_arr, R_arr)
            R_factor = d * (d - 1)
            c_eff = slope / R_factor if R_factor > 0 else float("nan")

            lines.append(
                f"- d={d}, N={N}: slope(R_hat vs H²) = {slope:.1f}, "
                f"R² = {r_val**2:.3f}, c_eff = slope/{R_factor} = **{c_eff:.1f}**"
            )

    # ── Strategy C: Per-H calibration constants ──
    lines.append("\n## Strategy C: Per-H Calibration at N=2048\n")
    lines.append("Best estimate using largest N.\n")
    lines.append("| d | H | R_dS | mean(R_hat) | bias-corrected | c = corr/R_dS |")
    lines.append("|---|---|------|-------------|----------------|---------------|")

    for d in ds:
        bias = flat_bias.get((d, 2048), 0.0)
        for H in hs:
            subset = [r for r in data if r["d"] == d and r["N"] == 2048 and r["H"] == H]
            if not subset:
                continue
            R_dS = d * (d - 1) * H ** 2
            mean_rhat = float(np.mean([r["R_hat"] for r in subset]))
            corrected = mean_rhat - bias
            c = corrected / R_dS if R_dS > 0 else float("nan")
            lines.append(
                f"| {d} | {H} | {R_dS:.2f} | {mean_rhat:+.1f} | "
                f"{corrected:+.1f} | {c:+.1f} |"
            )

    # ── Summary ──
    lines.append("\n## Summary\n")
    for d in ds:
        d_ns = sorted([N for N in ns if (d, N) in alpha_table])
        if d_ns:
            last = alpha_table[(d, d_ns[-1])]
            lines.append(
                f"- **d={d}**: α(N={d_ns[-1]}) = **{last['alpha']:.3f}** "
                f"(target=2.0, deviation={last['alpha']-2.0:+.3f}, R²={last['r2']:.3f})"
            )

    lines.append("")
    lines.append("**Interpretation**: If α > 2, R_hat responds super-linearly to curvature —")
    lines.append("the discrete proxy amplifies strong curvature. If α < 2, it under-responds.")
    lines.append("Convergence α → 2 as N → ∞ would confirm correct continuum scaling.")

    report_text = "\n".join(lines) + "\n"
    report_path = base / "curvature_calibration_powerlaw.md"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"\nSaved: {report_path}")

    # Console summary
    print("\n" + "=" * 60)
    print("POWER-LAW CALIBRATION: R_hat = c · H^α")
    print("=" * 60)
    for d in ds:
        d_ns = sorted([N for N in ns if (d, N) in alpha_table])
        print(f"\n  d={d}:")
        for N in d_ns:
            r = alpha_table[(d, N)]
            print(f"    N={N:5d}: α = {r['alpha']:.3f} ± {r['se']:.3f}  R²={r['r2']:.3f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
