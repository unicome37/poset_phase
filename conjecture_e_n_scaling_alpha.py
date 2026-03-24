"""Conjecture E — T5: N-Scaling of α_eff (EH Bridge Convergence Test).

Tests the central hypothesis: does α_eff(N, d=4) → 2 as N → ∞?

If discrete observables track H^α with α(N) → 2, then the H→R bridge
is a continuum-limit phenomenon and the BDG theorem is confirmed
constructively.

Design:
  - d=4 only (critical dimension where α≈1 at N=512)
  - N = 128, 256, 512, 1024
  - H = 0, 0.25, 0.5, 1.0, 2.0
  - 16 reps per (N, H) cell → 320 total realizations
  - Features: w_max_ratio (antichain), b1_std (spectral), mean_layer_width
  - NO eigenvalue computation (skipped for speed at large N)
  - Per-N pooled group-mean α scan via R² grid search
  - Fit α_eff(N) = 2 - c·N^(-γ) to extrapolate

References:
  - §4.1.32: α≈1.25 (w_max_ratio), α≈1.00 (b1_std) at d=4, N≤512
  - §4.1.33: No algebraic bridge works (0/9)
  - T4: No spectral ratio achieves α≈2 (confirmed §4.1.33)
"""

from __future__ import annotations

import argparse
import csv
import time
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats
from scipy.optimize import curve_fit

from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import build_causal_matrix_de_sitter
from conjecture_e_sorkin_dalembertian import build_dalembertian_matrix, compute_b1_features
from conjecture_e_antichain_structure import compute_antichain_features


# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ExpRow:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    H2: float
    n_causal_pairs: int
    w_max_ratio: float
    n_layers: float
    layer_ratio: float
    mean_layer_width: float
    layer_width_std: float
    b1_mean: float
    b1_std: float


def run_single(d: int, N: int, hubble: float, rep: int, seed: int) -> ExpRow:
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 1000),
    )
    causal = build_causal_matrix_de_sitter(points, hubble)
    n_causal_pairs = int(np.sum(causal))
    R_dS = d * (d - 1) * hubble ** 2
    H2 = hubble ** 2

    # Antichain features (skip Dilworth for N>600)
    ac_feats = compute_antichain_features(causal, compute_dilworth=(N <= 600))

    # B_ℓ spectral: only b1_mean/b1_std (no eigenvalues — too expensive)
    B = build_dalembertian_matrix(causal, d, rho=1.0)
    ones = np.ones(N)
    b1 = B @ ones
    b1_mean = float(np.mean(b1))
    b1_std = float(np.std(b1))

    return ExpRow(
        d=d, N=N, hubble=hubble, rep=rep,
        R_dS=R_dS, H2=H2,
        n_causal_pairs=n_causal_pairs,
        w_max_ratio=ac_feats.get("w_max_ratio", float("nan")),
        n_layers=ac_feats.get("n_layers", float("nan")),
        layer_ratio=ac_feats.get("layer_ratio", float("nan")),
        mean_layer_width=ac_feats.get("mean_layer_width", float("nan")),
        layer_width_std=ac_feats.get("layer_width_std", float("nan")),
        b1_mean=b1_mean,
        b1_std=b1_std,
    )


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------
def ols_residualize(feat: np.ndarray, dens: np.ndarray) -> np.ndarray:
    coeffs = np.polyfit(dens, feat, 1)
    return feat - np.polyval(coeffs, dens)


def ols_r2(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return 0.0
    coeffs = np.polyfit(x, y, 1)
    pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return max(0.0, 1.0 - ss_res / ss_tot) if ss_tot > 1e-15 else 0.0


def group_means(H_arr: np.ndarray, Y_arr: np.ndarray):
    H_levels = sorted(set(np.round(H_arr, 6)))
    H_out, means, stds, counts = [], [], [], []
    for h in H_levels:
        mask = np.abs(H_arr - h) < 1e-6
        if mask.sum() > 0:
            H_out.append(h)
            means.append(np.mean(Y_arr[mask]))
            stds.append(np.std(Y_arr[mask]))
            counts.append(int(mask.sum()))
    return np.array(H_out), np.array(means), np.array(stds), np.array(counts)


def alpha_grid(H_arr: np.ndarray, Y_arr: np.ndarray) -> dict:
    alpha_range = np.arange(0.25, 8.05, 0.25)
    results = {}
    for alpha in alpha_range:
        x = np.power(np.abs(H_arr), alpha)
        results[round(alpha, 2)] = ols_r2(x, Y_arr)
    return results


def best_alpha(ag: dict) -> tuple[float, float]:
    sorted_ag = sorted(ag.items(), key=lambda x: -x[1])
    return sorted_ag[0]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260324)
    parser.add_argument("--skip-1024", action="store_true",
                        help="Skip N=1024 for faster testing")
    args = parser.parse_args()

    d = 4
    ns = [128, 256, 512] if args.skip_1024 else [128, 256, 512, 1024]
    hubbles = [0.0, 0.25, 0.5, 1.0, 2.0]
    reps = args.reps
    total = len(ns) * len(hubbles) * reps

    out_dir = Path(__file__).parent / "outputs_unified_functional"
    out_dir.mkdir(exist_ok=True)

    print(f"T5: N-scaling of α_eff (d={d})")
    print(f"  N values: {ns}")
    print(f"  H values: {hubbles}")
    print(f"  Reps: {reps}")
    print(f"  Total: {total} realizations")
    print()

    rows: list[ExpRow] = []
    t0 = time.time()
    count = 0

    for N in ns:
        t_n = time.time()
        for h in hubbles:
            for rep in range(reps):
                row = run_single(d, N, h, rep, args.seed)
                rows.append(row)
                count += 1
                if count % 20 == 0:
                    elapsed = time.time() - t0
                    eta = elapsed / count * (total - count)
                    print(f"  [{count}/{total}] N={N} H={h} rep={rep} "
                          f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")
        dt_n = time.time() - t_n
        print(f"  N={N} complete ({dt_n:.1f}s)")

    total_time = time.time() - t0
    print(f"\nAll realizations complete ({total_time:.1f}s)")

    # Save CSV
    csv_path = out_dir / "conjecture_e_n_scaling_alpha.csv"
    field_names = [f.name for f in fields(ExpRow)]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=field_names)
        writer.writeheader()
        for r in rows:
            writer.writerow({fn: getattr(r, fn) for fn in field_names})
    print(f"CSV saved: {csv_path}")

    # ---------------------------------------------------------------------------
    # Analysis
    # ---------------------------------------------------------------------------
    features = ["w_max_ratio", "mean_layer_width", "layer_width_std", "b1_std"]
    lines: list[str] = []

    lines.append("# T5: N-Scaling of α_eff — EH Bridge Convergence Test\n")
    lines.append("## Design\n")
    lines.append(f"- d = {d}")
    lines.append(f"- N = {ns}")
    lines.append(f"- H = {hubbles}")
    lines.append(f"- Reps = {reps}")
    lines.append(f"- Total = {len(rows)} realizations")
    lines.append(f"- Runtime = {total_time:.0f}s\n")
    lines.append("## Hypothesis\n")
    lines.append("If α_eff(N) → 2 as N → ∞, the H→R bridge is a continuum-limit")
    lines.append("phenomenon and the BDG theorem is constructively confirmed.\n")

    # ================================================================
    # Section 1: Per-N density-residual correlation with H²
    # ================================================================
    lines.append("## 1. Density-Residual ρ(residual, H²) per N\n")
    lines.append("| N | n_rows | " + " | ".join(features) + " |")
    lines.append("|---|--------|" + "|".join(["------"] * len(features)) + "|")

    # Store for α analysis
    per_n_data: dict[int, dict[str, tuple[np.ndarray, np.ndarray]]] = {}

    for N in ns:
        subset = [r for r in rows if r.N == N]
        H_arr = np.array([r.hubble for r in subset])
        h2_arr = H_arr ** 2
        dens_arr = np.array([float(r.n_causal_pairs) for r in subset])

        per_n_data[N] = {}
        cells = []

        for feat in features:
            fa = np.array([getattr(r, feat) for r in subset])
            mask = ~np.isnan(fa)
            if mask.sum() < 15:
                cells.append("N/A")
                continue
            resid = ols_residualize(fa[mask], dens_arr[mask])
            rho_r, p_r = sp_stats.spearmanr(h2_arr[mask], resid)
            per_n_data[N][feat] = (H_arr[mask], resid)
            sig = "**" if abs(rho_r) > 0.3 and p_r < 0.01 else ""
            cells.append(f"{sig}{rho_r:+.3f}{sig}")

        lines.append(f"| {N} | {len(subset)} | " + " | ".join(cells) + " |")

    # ================================================================
    # Section 2: Per-N α grid scan (KEY RESULT)
    # ================================================================
    lines.append("\n## 2. Per-N α Grid Scan (KEY RESULT)\n")
    lines.append("Group-mean R² power-law scan: which H^α best fits the density residual?\n")

    lines.append("| N | Feature | α_best | R²_best | R²(α=1) | R²(α=2) | R²(2)/R²(1) |")
    lines.append("|---|---------|--------|---------|---------|---------|-------------|")

    alpha_results: dict[str, list[tuple[int, float, float]]] = {f: [] for f in features}

    for N in ns:
        if N not in per_n_data:
            continue
        for feat in features:
            if feat not in per_n_data[N]:
                lines.append(f"| {N} | {feat} | N/A | N/A | N/A | N/A | N/A |")
                continue

            H_arr, resid = per_n_data[N][feat]
            nonzero = H_arr > 1e-10
            if nonzero.sum() < 10:
                lines.append(f"| {N} | {feat} | N/A | N/A | N/A | N/A | N/A |")
                continue

            H_nz = H_arr[nonzero]
            R_nz = resid[nonzero]

            # Group means
            H_gm, R_gm, R_std, R_cnt = group_means(H_nz, R_nz)

            if len(H_gm) < 3:
                lines.append(f"| {N} | {feat} | N/A | N/A | N/A | N/A | N/A |")
                continue

            ag = alpha_grid(H_gm, R_gm)
            ab, r2b = best_alpha(ag)
            r2_1 = ag.get(1.0, 0.0)
            r2_2 = ag.get(2.0, 0.0)
            ratio = r2_2 / r2_1 if r2_1 > 1e-10 else float("inf")

            alpha_results[feat].append((N, ab, r2b))

            lines.append(
                f"| {N} | {feat} | **{ab:.2f}** | {r2b:.4f} "
                f"| {r2_1:.4f} | {r2_2:.4f} | {ratio:.3f} |"
            )

    # ================================================================
    # Section 3: α_eff(N) trend (THE KEY TABLE)
    # ================================================================
    lines.append("\n## 3. α_eff(N) Trend — Does α → 2 as N → ∞?\n")

    lines.append("| Feature | " + " | ".join(f"α(N={N})" for N in ns) + " | ρ(N, α) | Trend |")
    lines.append("|---------|" + "|".join(["------"] * len(ns)) + "|---------|-------|")

    for feat in features:
        ar = alpha_results[feat]
        if len(ar) < 2:
            continue
        ns_vals = [a[0] for a in ar]
        al_vals = [a[1] for a in ar]

        cells = []
        for N in ns:
            match = [a for a in ar if a[0] == N]
            if match:
                cells.append(f"**{match[0][1]:.2f}**")
            else:
                cells.append("N/A")

        ns_a = np.array(ns_vals, dtype=float)
        al_a = np.array(al_vals)
        if np.std(al_a) < 1e-10:
            rho_trend = 0.0
        else:
            rho_trend, _ = sp_stats.spearmanr(ns_a, al_a)

        if rho_trend > 0.5 and al_vals[-1] > al_vals[0]:
            trend = "↑ CONVERGING toward 2"
        elif abs(rho_trend) < 0.3:
            trend = "→ flat"
        else:
            trend = "↓ decreasing or unstable"

        lines.append(
            f"| {feat} | " + " | ".join(cells) +
            f" | {rho_trend:+.2f} | {trend} |"
        )

    # ================================================================
    # Section 4: R²(α=2)/R²(α=1) ratio trend
    # ================================================================
    lines.append("\n## 4. R²(α=2)/R²(α=1) Ratio Trend\n")
    lines.append("If α_eff → 2, this ratio should increase with N.\n")

    lines.append("| Feature | " + " | ".join(f"ratio(N={N})" for N in ns) + " | ρ(N, ratio) | Trend |")
    lines.append("|---------|" + "|".join(["------"] * len(ns)) + "|-------------|-------|")

    for feat in features:
        ratio_list = []
        cells = []
        for N in ns:
            if N not in per_n_data or feat not in per_n_data[N]:
                cells.append("N/A")
                continue
            H_arr, resid = per_n_data[N][feat]
            nonzero = H_arr > 1e-10
            if nonzero.sum() < 10:
                cells.append("N/A")
                continue
            H_nz = H_arr[nonzero]
            R_nz = resid[nonzero]
            H_gm, R_gm, _, _ = group_means(H_nz, R_nz)
            if len(H_gm) < 3:
                cells.append("N/A")
                continue
            ag = alpha_grid(H_gm, R_gm)
            r2_1 = ag.get(1.0, 0.0)
            r2_2 = ag.get(2.0, 0.0)
            ratio = r2_2 / r2_1 if r2_1 > 1e-10 else float("inf")
            ratio_list.append((N, ratio))
            cells.append(f"{ratio:.3f}")

        if len(ratio_list) >= 2:
            ns_a = np.array([r[0] for r in ratio_list], dtype=float)
            rt_a = np.array([r[1] for r in ratio_list])
            rho_r, _ = sp_stats.spearmanr(ns_a, rt_a)
            trend = "↑ R preference growing" if rho_r > 0.5 else "→ flat" if abs(rho_r) < 0.3 else "↓ H preference growing"
        else:
            rho_r = float("nan")
            trend = "N/A"

        lines.append(
            f"| {feat} | " + " | ".join(cells) +
            f" | {rho_r:+.2f} | {trend} |"
        )

    # ================================================================
    # Section 5: Group-mean curves (for visual inspection)
    # ================================================================
    lines.append("\n## 5. Group-Mean Curves per N\n")
    lines.append("Mean ± std of density-residualized feature at each H level.\n")

    for feat in ["w_max_ratio", "b1_std"]:
        lines.append(f"### {feat}\n")
        lines.append("| H | " + " | ".join(f"N={N}" for N in ns) + " |")
        lines.append("|---|" + "|".join(["------"] * len(ns)) + "|")

        for h in hubbles:
            cells = []
            for N in ns:
                if N not in per_n_data or feat not in per_n_data[N]:
                    cells.append("N/A")
                    continue
                H_arr, resid = per_n_data[N][feat]
                mask = np.abs(H_arr - h) < 1e-6
                if mask.sum() == 0:
                    cells.append("N/A")
                    continue
                m = np.mean(resid[mask])
                s = np.std(resid[mask])
                cells.append(f"{m:+.4f}±{s:.4f}")
            lines.append(f"| {h} | " + " | ".join(cells) + " |")
        lines.append("")

    # ================================================================
    # Section 6: Verdict
    # ================================================================
    lines.append("## 6. Verdict\n")

    # Check if any feature shows α increasing with N
    converging = []
    for feat in features:
        ar = alpha_results[feat]
        if len(ar) >= 3:
            al_vals = [a[1] for a in ar]
            ns_vals = [a[0] for a in ar]
            rho_t, _ = sp_stats.spearmanr(ns_vals, al_vals)
            if rho_t > 0.5 and al_vals[-1] > al_vals[0]:
                converging.append((feat, al_vals, rho_t))

    if converging:
        lines.append("### ✅ Convergence signal detected\n")
        for feat, al_vals, rho_t in converging:
            lines.append(f"- **{feat}**: α = {' → '.join(f'{a:.2f}' for a in al_vals)}, ρ(N,α) = {rho_t:+.2f}")
        lines.append("")
        lines.append("At least one feature shows α_eff increasing with N toward 2.")
        lines.append("This supports the hypothesis that the H→R bridge is a continuum-limit phenomenon.\n")

        # Attempt extrapolation
        lines.append("### Extrapolation: α_eff(N) = 2 - c·N^(-γ)\n")
        for feat, al_vals, rho_t in converging:
            ar = alpha_results[feat]
            ns_a = np.array([a[0] for a in ar], dtype=float)
            al_a = np.array([a[1] for a in ar])

            # Only attempt fit if α values are reasonable (not saturated at 8)
            if max(al_a) > 3.0:
                lines.append(f"- **{feat}**: α values too scattered for reliable extrapolation")
                continue

            try:
                def model(N, c, gamma):
                    return 2.0 - c * np.power(N, -gamma)

                popt, pcov = curve_fit(model, ns_a, al_a, p0=[10, 0.5],
                                       bounds=([0, 0], [1000, 5]))
                c_fit, gamma_fit = popt
                perr = np.sqrt(np.diag(pcov))

                # Predict N needed for α=1.9 (90% of target)
                if c_fit > 0 and gamma_fit > 0:
                    N_90 = (c_fit / 0.1) ** (1.0 / gamma_fit)
                else:
                    N_90 = float("inf")

                lines.append(f"- **{feat}**: c = {c_fit:.2f}±{perr[0]:.2f}, "
                             f"γ = {gamma_fit:.3f}±{perr[1]:.3f}")
                lines.append(f"  - Predicted N for α=1.9: **{N_90:.0f}**")
                lines.append(f"  - Current α at max N: {al_a[-1]:.2f}")
            except Exception as e:
                lines.append(f"- **{feat}**: fit failed ({e})")
        lines.append("")
    else:
        lines.append("### ❌ No clear convergence signal\n")
        lines.append("α_eff does not systematically increase with N at d=4.")
        lines.append("Possible explanations:")
        lines.append("1. The convergence is too slow to detect at N ≤ 1024")
        lines.append("2. The α scan methodology is unreliable (step-function artifact)")
        lines.append("3. The discrete observables genuinely track H, not R, at all N\n")

    # Final summary
    lines.append("### Summary\n")
    lines.append("| Question | Answer |")
    lines.append("|----------|--------|")
    if converging:
        lines.append("| Does α_eff increase with N? | **Yes** (at least partial) |")
        lines.append("| Is the bridge a continuum-limit phenomenon? | **Supported** |")
    else:
        lines.append("| Does α_eff increase with N? | **No clear signal** |")
        lines.append("| Is the bridge a continuum-limit phenomenon? | **Inconclusive** |")
    lines.append("| E-bulk-second-order confidence change? | See above |")

    lines.append(f"\n---\n")
    lines.append(f"*Generated by conjecture_e_n_scaling_alpha.py*")
    lines.append(f"*{len(rows)} realizations, d={d}, runtime={total_time:.0f}s*")

    # Write report
    report_path = out_dir / "conjecture_e_n_scaling_alpha.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Report saved: {report_path}")


if __name__ == "__main__":
    main()
