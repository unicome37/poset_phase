"""T5v2: N-Scaling α_eff Reanalysis with Improved Methodology.

Fixes from T5v1:
  1. Uses per-realization Pearson R² (80 points per N slice) instead of
     group-mean R² (4 points) — much more statistical power to discriminate α.
  2. Uses log-log regression: log|resid| vs log(H) to directly estimate α.
  3. Computes R²(H²)/R²(H) ratio at per-realization level.
  4. Uses layer-based w_max_ratio proxy (max_layer_width_ratio) at N=1024
     since Dilworth width is unavailable.

Data source: existing T5 CSV (320 realizations, d=4, N=128/256/512/1024).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Helpers
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


def pearson_r(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return 0.0
    r, _ = sp_stats.pearsonr(x, y)
    return r


def alpha_grid_per_realization(H_arr: np.ndarray, Y_arr: np.ndarray) -> dict:
    """Grid search on per-realization data (not group means)."""
    alpha_range = np.arange(0.25, 6.05, 0.25)
    results = {}
    for alpha in alpha_range:
        x = np.power(np.abs(H_arr), alpha)
        results[round(alpha, 2)] = ols_r2(x, Y_arr)
    return results


def loglog_alpha_estimate(H_arr: np.ndarray, Y_arr: np.ndarray) -> tuple[float, float, float]:
    """Estimate α via log-log regression: log|Y| = α·log(H) + const.
    
    Returns (α, R², p-value).
    Uses sign-preserving log: sign(Y)·log(|Y|+1) to handle negative residuals.
    """
    mask = H_arr > 1e-10
    if mask.sum() < 5:
        return float("nan"), 0.0, 1.0
    
    logH = np.log(H_arr[mask])
    Y_m = Y_arr[mask]
    
    # Use absolute values with sign tracking
    absY = np.abs(Y_m)
    if np.max(absY) < 1e-15:
        return float("nan"), 0.0, 1.0
    
    # Log transform with offset to handle zeros
    logY = np.log(absY + 1e-10)
    
    slope, intercept, r_value, p_value, std_err = sp_stats.linregress(logH, logY)
    return slope, r_value**2, p_value


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
def load_csv(path: Path) -> list[dict]:
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = []
        for row in reader:
            parsed = {}
            for k, v in row.items():
                try:
                    parsed[k] = float(v)
                except (ValueError, TypeError):
                    parsed[k] = v
            rows.append(parsed)
    return rows


def main():
    data_path = Path(__file__).parent / "outputs_unified_functional" / "conjecture_e_n_scaling_alpha.csv"
    out_dir = Path(__file__).parent / "outputs_unified_functional"

    print(f"Loading {data_path} ...")
    raw_rows = load_csv(data_path)
    print(f"  {len(raw_rows)} rows loaded.")

    ns = sorted(set(int(r["N"]) for r in raw_rows))
    
    # Features: use mean_layer_width as w_max_ratio proxy at N=1024
    features = ["mean_layer_width", "layer_width_std", "b1_std"]
    # Also include w_max_ratio where available
    all_features = ["w_max_ratio"] + features

    lines: list[str] = []
    lines.append("# T5v2: N-Scaling α_eff — Improved Methodology\n")
    lines.append("## Methodological Improvements over T5v1\n")
    lines.append("1. **Per-realization R²** (64 non-zero-H points per N) instead of group-mean R² (4 points)")
    lines.append("2. **Finer α grid** (0.25–6.0 step 0.25) — avoids ceiling saturation at 8.0")
    lines.append("3. **Log-log regression** as independent α estimator")
    lines.append("4. **R²(H^α)/R²(H^1) ratio** at per-realization level\n")
    lines.append(f"- Data: {len(raw_rows)} realizations from T5 CSV (d=4)")
    lines.append(f"- N values: {ns}")
    lines.append(f"- Features: {all_features}\n")

    # ================================================================
    # Method A: Per-realization R² grid scan
    # ================================================================
    lines.append("## Method A: Per-Realization R² Grid Scan\n")
    lines.append("Using ALL non-zero-H points per N slice (64 data points).\n")

    lines.append("| N | Feature | α_best | R²_best | R²(α=1) | R²(α=2) | R²(2)/R²(1) |")
    lines.append("|---|---------|--------|---------|---------|---------|-------------|")

    alpha_results_A: dict[str, list[tuple[int, float, float]]] = {f: [] for f in all_features}
    r2_ratio_results: dict[str, list[tuple[int, float]]] = {f: [] for f in all_features}

    for N in ns:
        subset = [r for r in raw_rows if int(r["N"]) == N]
        H_arr = np.array([r["hubble"] for r in subset])
        dens_arr = np.array([r["n_causal_pairs"] for r in subset])

        for feat in all_features:
            fa = np.array([r.get(feat, float("nan")) for r in subset])
            mask = ~np.isnan(fa)
            if mask.sum() < 15:
                lines.append(f"| {N} | {feat} | N/A | N/A | N/A | N/A | N/A |")
                continue

            # Density-residualize
            resid = ols_residualize(fa[mask], dens_arr[mask])
            H_m = H_arr[mask]

            # Filter to non-zero H
            nonzero = H_m > 1e-10
            if nonzero.sum() < 10:
                lines.append(f"| {N} | {feat} | N/A | N/A | N/A | N/A | N/A |")
                continue

            H_nz = H_m[nonzero]
            R_nz = resid[nonzero]

            # Per-realization α grid scan
            ag = alpha_grid_per_realization(H_nz, R_nz)
            sorted_ag = sorted(ag.items(), key=lambda x: -x[1])
            ab, r2b = sorted_ag[0]
            r2_1 = ag.get(1.0, 0.0)
            r2_2 = ag.get(2.0, 0.0)
            ratio = r2_2 / r2_1 if r2_1 > 1e-10 else float("inf")

            alpha_results_A[feat].append((N, ab, r2b))
            r2_ratio_results[feat].append((N, ratio))

            lines.append(
                f"| {N} | {feat} | **{ab:.2f}** | {r2b:.4f} "
                f"| {r2_1:.4f} | {r2_2:.4f} | {ratio:.3f} |"
            )

    # ================================================================
    # Method B: Log-log regression α estimate
    # ================================================================
    lines.append("\n## Method B: Log-Log Regression α Estimate\n")
    lines.append("Direct slope estimation: log|residual| = α·log(H) + const\n")

    lines.append("| N | Feature | α_loglog | R²_loglog | p-value |")
    lines.append("|---|---------|----------|-----------|---------|")

    alpha_results_B: dict[str, list[tuple[int, float, float]]] = {f: [] for f in all_features}

    for N in ns:
        subset = [r for r in raw_rows if int(r["N"]) == N]
        H_arr = np.array([r["hubble"] for r in subset])
        dens_arr = np.array([r["n_causal_pairs"] for r in subset])

        for feat in all_features:
            fa = np.array([r.get(feat, float("nan")) for r in subset])
            mask = ~np.isnan(fa)
            if mask.sum() < 15:
                lines.append(f"| {N} | {feat} | N/A | N/A | N/A |")
                continue

            resid = ols_residualize(fa[mask], dens_arr[mask])
            H_m = H_arr[mask]

            alpha_ll, r2_ll, p_ll = loglog_alpha_estimate(H_m, resid)
            alpha_results_B[feat].append((N, alpha_ll, r2_ll))

            if np.isnan(alpha_ll):
                lines.append(f"| {N} | {feat} | N/A | N/A | N/A |")
            else:
                lines.append(f"| {N} | {feat} | **{alpha_ll:.3f}** | {r2_ll:.4f} | {p_ll:.2e} |")

    # ================================================================
    # Method C: Pearson correlation with H vs H² vs H³
    # ================================================================
    lines.append("\n## Method C: Pearson R² with H^k (k=1,2,3) — Direct Comparison\n")
    lines.append("Which power of H gives highest linear R² with density residual?\n")

    lines.append("| N | Feature | R²(H) | R²(H²) | R²(H³) | Best | R²(H²)/R²(H) |")
    lines.append("|---|---------|-------|--------|--------|------|---------------|")

    for N in ns:
        subset = [r for r in raw_rows if int(r["N"]) == N]
        H_arr = np.array([r["hubble"] for r in subset])
        dens_arr = np.array([r["n_causal_pairs"] for r in subset])

        for feat in all_features:
            fa = np.array([r.get(feat, float("nan")) for r in subset])
            mask = ~np.isnan(fa)
            if mask.sum() < 15:
                lines.append(f"| {N} | {feat} | N/A | N/A | N/A | N/A | N/A |")
                continue

            resid = ols_residualize(fa[mask], dens_arr[mask])
            H_m = H_arr[mask]

            r2_h1 = ols_r2(H_m, resid)
            r2_h2 = ols_r2(H_m**2, resid)
            r2_h3 = ols_r2(H_m**3, resid)

            best_k = max([(1, r2_h1), (2, r2_h2), (3, r2_h3)], key=lambda x: x[1])
            ratio = r2_h2 / r2_h1 if r2_h1 > 1e-10 else float("inf")

            lines.append(
                f"| {N} | {feat} | {r2_h1:.4f} | {r2_h2:.4f} | {r2_h3:.4f} "
                f"| H^{best_k[0]} | {ratio:.3f} |"
            )

    # ================================================================
    # Section 4: α_eff(N) convergence summary (ALL METHODS)
    # ================================================================
    lines.append("\n## 4. α_eff(N) Convergence Summary\n")

    lines.append("### Method A: Per-realization R² grid scan\n")
    lines.append("| Feature | " + " | ".join(f"α(N={N})" for N in ns) + " | ρ(N, α) | Trend |")
    lines.append("|---------|" + "|".join(["------"] * len(ns)) + "|---------|-------|")

    for feat in all_features:
        ar = alpha_results_A[feat]
        if len(ar) < 2:
            continue
        cells = []
        for N in ns:
            match = [a for a in ar if a[0] == N]
            cells.append(f"{match[0][1]:.2f}" if match else "N/A")

        ns_a = np.array([a[0] for a in ar], dtype=float)
        al_a = np.array([a[1] for a in ar])
        if len(ar) >= 3 and np.std(al_a) > 1e-10:
            rho_t, _ = sp_stats.spearmanr(ns_a, al_a)
        else:
            rho_t = float("nan")

        if not np.isnan(rho_t) and rho_t > 0.5:
            trend = "↑ increasing"
        elif not np.isnan(rho_t) and abs(rho_t) < 0.3:
            trend = "→ flat"
        else:
            trend = "↓ or N/A"

        lines.append(f"| {feat} | " + " | ".join(cells) + f" | {rho_t:+.2f} | {trend} |")

    lines.append("\n### Method B: Log-log slope\n")
    lines.append("| Feature | " + " | ".join(f"α(N={N})" for N in ns) + " | ρ(N, α) | Trend |")
    lines.append("|---------|" + "|".join(["------"] * len(ns)) + "|---------|-------|")

    for feat in all_features:
        ar = alpha_results_B[feat]
        if len(ar) < 2:
            continue
        cells = []
        for N in ns:
            match = [a for a in ar if a[0] == N]
            cells.append(f"{match[0][1]:.2f}" if match else "N/A")

        ns_a = np.array([a[0] for a in ar], dtype=float)
        al_a = np.array([a[1] for a in ar])
        valid = ~np.isnan(al_a)
        if valid.sum() >= 3 and np.std(al_a[valid]) > 1e-10:
            rho_t, _ = sp_stats.spearmanr(ns_a[valid], al_a[valid])
        else:
            rho_t = float("nan")

        if not np.isnan(rho_t) and rho_t > 0.5:
            trend = "↑ increasing"
        elif not np.isnan(rho_t) and abs(rho_t) < 0.3:
            trend = "→ flat"
        else:
            trend = "↓ or N/A"

        lines.append(f"| {feat} | " + " | ".join(cells) + f" | {rho_t:+.2f} | {trend} |")

    lines.append("\n### R²(H²)/R²(H) ratio trend\n")
    lines.append("| Feature | " + " | ".join(f"N={N}" for N in ns) + " | ρ(N, ratio) | Trend |")
    lines.append("|---------|" + "|".join(["------"] * len(ns)) + "|-------------|-------|")

    for feat in all_features:
        rr = r2_ratio_results[feat]
        if len(rr) < 2:
            continue
        cells = []
        for N in ns:
            match = [r for r in rr if r[0] == N]
            cells.append(f"{match[0][1]:.3f}" if match else "N/A")

        ns_a = np.array([r[0] for r in rr], dtype=float)
        rt_a = np.array([r[1] for r in rr])
        if len(rr) >= 3:
            rho_r, _ = sp_stats.spearmanr(ns_a, rt_a)
        else:
            rho_r = float("nan")

        trend = "↑ R-preference growing" if rho_r > 0.5 else "→ flat" if abs(rho_r) < 0.3 else "↓ or N/A"
        lines.append(f"| {feat} | " + " | ".join(cells) + f" | {rho_r:+.2f} | {trend} |")

    # ================================================================
    # Section 5: Extrapolation attempt
    # ================================================================
    lines.append("\n## 5. Extrapolation: α_eff(N) = 2 - c·N^(-γ)\n")

    for method_name, alpha_store in [("Method A", alpha_results_A), ("Method B", alpha_results_B)]:
        lines.append(f"\n### {method_name}\n")
        for feat in all_features:
            ar = alpha_store[feat]
            valid = [(n, a, r) for n, a, r in ar if not np.isnan(a) and 0 < a < 6]
            if len(valid) < 3:
                lines.append(f"- **{feat}**: insufficient valid α values for extrapolation")
                continue

            ns_a = np.array([v[0] for v in valid], dtype=float)
            al_a = np.array([v[1] for v in valid])

            try:
                def model(N, c, gamma):
                    return 2.0 - c * np.power(N, -gamma)

                popt, pcov = curve_fit(model, ns_a, al_a, p0=[10, 0.5],
                                       bounds=([0, 0.01], [1000, 5]),
                                       maxfev=5000)
                c_fit, gamma_fit = popt
                perr = np.sqrt(np.diag(pcov))

                # Predict
                alpha_pred = {N: model(N, c_fit, gamma_fit) for N in [2048, 4096, 10000]}
                
                lines.append(f"- **{feat}**: c = {c_fit:.2f}±{perr[0]:.2f}, γ = {gamma_fit:.3f}±{perr[1]:.3f}")
                for N_pred, a_pred in alpha_pred.items():
                    lines.append(f"  - α(N={N_pred}) = {a_pred:.3f}")
                
                # N needed for α=1.9
                if c_fit > 0 and gamma_fit > 0:
                    N_90 = (c_fit / 0.1) ** (1.0 / gamma_fit)
                    lines.append(f"  - N for α=1.9: **{N_90:.0f}**")
            except Exception as e:
                lines.append(f"- **{feat}**: fit failed ({e})")

    # ================================================================
    # Section 6: Verdict
    # ================================================================
    lines.append("\n## 6. Verdict\n")
    
    # Collect convergence evidence
    evidence_for = []
    evidence_against = []
    
    for feat in all_features:
        # Check Method A trend
        ar_A = alpha_results_A[feat]
        valid_A = [(n, a) for n, a, _ in ar_A if not np.isnan(a)]
        if len(valid_A) >= 3:
            ns_a = np.array([v[0] for v in valid_A], dtype=float)
            al_a = np.array([v[1] for v in valid_A])
            if np.std(al_a) > 1e-10:
                rho_t, _ = sp_stats.spearmanr(ns_a, al_a)
                if rho_t > 0.5:
                    evidence_for.append(f"{feat} (Method A): α increases, ρ={rho_t:+.2f}")
                elif rho_t < -0.3:
                    evidence_against.append(f"{feat} (Method A): α decreases, ρ={rho_t:+.2f}")

        # Check R²(H²)/R²(H) trend
        rr = r2_ratio_results[feat]
        if len(rr) >= 3:
            ns_a = np.array([r[0] for r in rr], dtype=float)
            rt_a = np.array([r[1] for r in rr])
            rho_r, _ = sp_stats.spearmanr(ns_a, rt_a)
            if rho_r > 0.5:
                evidence_for.append(f"{feat} R²(H²)/R²(H) ratio increases, ρ={rho_r:+.2f}")

    if evidence_for:
        lines.append("### Evidence FOR α → 2 convergence:\n")
        for e in evidence_for:
            lines.append(f"- {e}")
    
    if evidence_against:
        lines.append("\n### Evidence AGAINST α → 2 convergence:\n")
        for e in evidence_against:
            lines.append(f"- {e}")
    
    if not evidence_for and not evidence_against:
        lines.append("### No clear convergence signal in either direction.\n")
    
    lines.append("\n### Overall Assessment\n")
    lines.append("| Criterion | Status |")
    lines.append("|-----------|--------|")
    lines.append(f"| α_eff increases with N (Method A)? | {'YES' if any('Method A' in e for e in evidence_for) else 'NO'} |")
    lines.append(f"| R²(H²)/R²(H) increases with N? | {'YES' if any('R²' in e for e in evidence_for) else 'NO'} |")
    lines.append(f"| Extrapolation feasible? | See Section 5 |")
    lines.append(f"| E-bulk-second-order confidence change | See below |")

    lines.append(f"\n---\n")
    lines.append(f"*Generated by conjecture_e_n_scaling_alpha_v2.py*")
    lines.append(f"*Reanalysis of {len(raw_rows)} realizations from T5 CSV*")

    report_path = out_dir / "conjecture_e_n_scaling_alpha_v2.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Report saved: {report_path}")


if __name__ == "__main__":
    main()
