"""Conjecture E — T2: Cross-Scale Variance Experiment.

Tests whether the *noise structure* of H-tracking observables carries H²
information.  The idea (§3.2 of eh_bridge_theoretical_analysis.md):

  If observable O tracks H at first order, then its realisation-to-
  realisation scatter σ(O|H) may itself depend on H.  If σ ∝ H^β with
  β ≈ 1, then σ² ∝ H² ∼ R — giving a variance-based second-order
  estimator.

Two complementary approaches:

  **A  Inter-realization variance** (needs many reps per cell)
      For each (d, N, H) cell, compute std(obs) across reps.
      Fit std(obs) vs H to find β.

  **B  Intra-realization (sub-patch) variance**
      Within a single sprinkled causal set, partition elements into
      temporal slabs (layer blocks), compute a local version of the
      observable in each slab, and take the within-realization variance.
      Tests whether a *single* realisation's internal fluctuation
      encodes H².

Metrics:
  Q1: Does std(obs|H) correlate with H?  (Spearman ρ, significance)
  Q2: Power-law β from log-log fit  std(obs|H) ∝ H^β
  Q3: Variance-based R estimator V̂ = c · std²  → R²(V̂, H²)
  Q4: N-scaling — does the pattern strengthen with N?
  Q5: Comparison with T3 — is variance channel independent or redundant?

Design:
  d = 4
  N = 128, 256, 512, 1024
  H = 0, 0.1, 0.25, 0.5, 1.0, 2.0   (finer H grid for β fitting)
  32 reps per cell → 768 total
  → 6 H levels × 32 reps = 192 realizations per (d, N)

References:
  - T3 (§4.1.39): two-step squaring works; H² info is extraction-sensitive
  - T1 (§4.1.38): b1_std = unique spectral bulk observable
  - §4.1.28: w_max_ratio = strongest antichain observable
"""

from __future__ import annotations

import argparse
import csv
import time
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import build_causal_matrix_de_sitter
from conjecture_e_sorkin_dalembertian import build_dalembertian_matrix
from conjecture_e_antichain_structure import compute_antichain_features


# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class T2Row:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    H2: float
    n_causal_pairs: int
    # Spectral channel
    b1_mean: float
    b1_std: float
    # Antichain channel
    w_max_ratio: float
    mean_layer_width: float
    # Sub-patch (intra-realization) variance — layer-block approach
    n_patches: int
    b1_std_patch_mean: float
    b1_std_patch_std: float
    wmr_patch_mean: float
    wmr_patch_std: float


def _layer_assignment(causal: np.ndarray) -> np.ndarray:
    """Assign each node to a layer (topological level).
    Layer 0 = minimal elements (no predecessors)."""
    n = causal.shape[0]
    indeg = causal.sum(axis=0).astype(int)
    labels = np.full(n, -1, dtype=int)
    remaining = set(range(n))
    layer = 0
    while remaining:
        current = {i for i in remaining if indeg[i] == 0}
        if not current:
            # Cycle or isolated — assign everything remaining
            for i in remaining:
                labels[i] = layer
            break
        for i in current:
            labels[i] = layer
        remaining -= current
        for i in current:
            for j in range(n):
                if causal[i, j]:
                    indeg[j] -= 1
        layer += 1
    return labels


def _compute_patch_observables(causal: np.ndarray, d: int,
                               layer_labels: np.ndarray,
                               block_size: int = 3,
                               min_patch_n: int = 8) -> dict:
    """Partition causal set into layer blocks and compute observables per block.
    Returns dict with patch-level statistics."""
    n = causal.shape[0]
    n_layers = layer_labels.max() + 1

    # Group layers into blocks
    block_starts = list(range(0, n_layers, block_size))
    patch_b1_stds = []
    patch_wmrs = []

    for bs in block_starts:
        be = min(bs + block_size, n_layers)
        # Nodes in this layer block
        mask = np.zeros(n, dtype=bool)
        for l in range(bs, be):
            mask |= (layer_labels == l)
        patch_nodes = np.where(mask)[0]

        if len(patch_nodes) < min_patch_n:
            continue

        # Extract induced sub-poset
        sub_causal = causal[np.ix_(patch_nodes, patch_nodes)]
        n_patch = len(patch_nodes)

        # b1_std in this patch
        try:
            B_patch = build_dalembertian_matrix(sub_causal, d, rho=1.0)
            ones_patch = np.ones(n_patch)
            b1_patch = B_patch @ ones_patch
            patch_b1_stds.append(float(np.std(b1_patch)))
        except Exception:
            pass

        # w_max_ratio via layer decomposition of sub-patch
        try:
            ac = compute_antichain_features(sub_causal, compute_dilworth=False)
            wmr = ac.get("w_max_ratio", float("nan"))
            if not np.isnan(wmr):
                patch_wmrs.append(wmr)
        except Exception:
            pass

    result = {
        "n_patches": len(patch_b1_stds),
        "b1_std_patch_mean": float(np.mean(patch_b1_stds)) if patch_b1_stds else float("nan"),
        "b1_std_patch_std": float(np.std(patch_b1_stds)) if len(patch_b1_stds) >= 2 else float("nan"),
        "wmr_patch_mean": float(np.mean(patch_wmrs)) if patch_wmrs else float("nan"),
        "wmr_patch_std": float(np.std(patch_wmrs)) if len(patch_wmrs) >= 2 else float("nan"),
    }
    return result


def run_single(d: int, N: int, hubble: float, rep: int, seed: int) -> T2Row:
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 10000),
    )
    causal = build_causal_matrix_de_sitter(points, hubble)
    n_causal_pairs = int(np.sum(causal))
    R_dS = d * (d - 1) * hubble ** 2
    H2 = hubble ** 2

    # B_ℓ spectral features (global)
    B = build_dalembertian_matrix(causal, d, rho=1.0)
    ones = np.ones(N)
    b1 = B @ ones
    b1_mean = float(np.mean(b1))
    b1_std_val = float(np.std(b1))

    # Antichain features (global)
    ac_feats = compute_antichain_features(causal, compute_dilworth=(N <= 600))
    w_max_ratio = ac_feats.get("w_max_ratio", float("nan"))
    mean_layer_width = ac_feats.get("mean_layer_width", float("nan"))

    # Sub-patch variance (intra-realization)
    layer_labels = _layer_assignment(causal)
    patch_info = _compute_patch_observables(causal, d, layer_labels,
                                            block_size=3, min_patch_n=8)

    return T2Row(
        d=d, N=N, hubble=hubble, rep=rep,
        R_dS=R_dS, H2=H2,
        n_causal_pairs=n_causal_pairs,
        b1_mean=b1_mean, b1_std=b1_std_val,
        w_max_ratio=w_max_ratio,
        mean_layer_width=mean_layer_width,
        n_patches=patch_info["n_patches"],
        b1_std_patch_mean=patch_info["b1_std_patch_mean"],
        b1_std_patch_std=patch_info["b1_std_patch_std"],
        wmr_patch_mean=patch_info["wmr_patch_mean"],
        wmr_patch_std=patch_info["wmr_patch_std"],
    )


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------
def ols_r2(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return 0.0
    coeffs = np.polyfit(x, y, 1)
    pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return max(0.0, 1.0 - ss_res / ss_tot) if ss_tot > 1e-15 else 0.0


def ols_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 3 or np.std(x) < 1e-15:
        return (0.0, np.mean(y))
    coeffs = np.polyfit(x, y, 1)
    return (float(coeffs[0]), float(coeffs[1]))


def alpha_grid(H_arr: np.ndarray, Y_arr: np.ndarray,
               alpha_range: np.ndarray | None = None) -> dict:
    if alpha_range is None:
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
    parser.add_argument("--reps", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260326)
    args = parser.parse_args()

    dims = [4]
    ns_by_d = {4: [128, 256, 512, 1024]}
    hubbles = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0]
    hubbles_nonzero = [h for h in hubbles if h > 0]
    reps = args.reps

    total = sum(len(ns_by_d[d]) * len(hubbles) * reps for d in dims)

    out_dir = Path(__file__).parent / "outputs_unified_functional"
    out_dir.mkdir(exist_ok=True)

    print(f"T2: Cross-Scale Variance")
    print(f"  Dims: {dims}")
    print(f"  H values: {hubbles}")
    print(f"  Reps: {reps}")
    print(f"  Total: {total} realizations")
    print()

    rows: list[T2Row] = []
    t0 = time.time()
    count = 0

    for d in dims:
        ns = ns_by_d[d]
        for N in ns:
            t_n = time.time()
            for h in hubbles:
                for rep_i in range(reps):
                    row = run_single(d, N, h, rep_i, args.seed)
                    rows.append(row)
                    count += 1
                    if count % 20 == 0:
                        elapsed = time.time() - t0
                        eta = elapsed / count * (total - count)
                        print(f"  [{count}/{total}] d={d} N={N} H={h} rep={rep_i} "
                              f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")
            dt_n = time.time() - t_n
            print(f"  d={d} N={N} complete ({dt_n:.1f}s)")

    total_time = time.time() - t0
    print(f"\nAll realizations complete ({total_time:.1f}s)")

    # Save CSV
    csv_path = out_dir / "conjecture_e_cross_scale_variance.csv"
    field_names = [f.name for f in fields(T2Row)]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=field_names)
        writer.writeheader()
        for r in rows:
            writer.writerow({fn: getattr(r, fn) for fn in field_names})
    print(f"CSV saved: {csv_path}")

    # ---------------------------------------------------------------------------
    # Analysis
    # ---------------------------------------------------------------------------
    import pandas as pd

    df = pd.DataFrame([{fn: getattr(r, fn) for fn in field_names} for r in rows])

    md_lines: list[str] = []
    md = md_lines.append

    md("# T2: Cross-Scale Variance — Results\n")
    md("")
    md("## Design\n")
    md(f"- Dimensions: {dims}")
    md(f"- N values: { {d: ns_by_d[d] for d in dims} }")
    md(f"- H values: {hubbles}")
    md(f"- Reps per cell: {reps}")
    md(f"- Total realizations: {len(df)}")
    md(f"- Runtime: {total_time:.1f}s")
    md("")
    md("## Method\n")
    md("Two complementary approaches to extracting H² from noise structure:")
    md("")
    md("- **Approach A (inter-realization)**: For each (d, N, H) cell, compute")
    md("  std(observable) across reps. Test whether std(obs|H) ∝ H^β.")
    md("- **Approach B (intra-realization / sub-patch)**: Within each realization,")
    md("  partition into temporal layer blocks, compute observable per block,")
    md("  take within-realization std. Test whether patch_std correlates with H².")
    md("")

    # -----------------------------------------------------------------------
    # Q1: Inter-realization variance — std(obs|H) vs H
    # -----------------------------------------------------------------------
    md("## 1. Inter-Realization Variance: std(obs|H) vs H\n")
    md("| d | N | Observable | std(H=0) | std(H=0.1) | std(H=0.25) | std(H=0.5) | std(H=1.0) | std(H=2.0) | ρ(H, std) | p-value | β (log-log) |")
    md("|---|---|------------|----------|------------|-------------|------------|------------|------------|-----------|---------|-------------|")

    obs_names_a = ["b1_std", "w_max_ratio", "mean_layer_width"]
    # Store for later analysis
    inter_real_results = []

    for d in dims:
        for N in ns_by_d[d]:
            sub = df[(df["d"] == d) & (df["N"] == N)]
            for obs_name in obs_names_a:
                stds_by_h = {}
                for h in hubbles:
                    cell = sub[sub["hubble"] == h][obs_name].dropna()
                    stds_by_h[h] = float(cell.std()) if len(cell) >= 2 else float("nan")

                # Spearman of (H, std(obs|H)) across H levels (excluding H=0)
                h_vals = [h for h in hubbles_nonzero if not np.isnan(stds_by_h[h])]
                s_vals = [stds_by_h[h] for h in h_vals]
                if len(h_vals) >= 3:
                    rho_val, p_val = sp_stats.spearmanr(h_vals, s_vals)
                else:
                    rho_val, p_val = float("nan"), float("nan")

                # Power-law β via log-log OLS (H>0 only)
                h_pos = [h for h in h_vals if h > 0 and stds_by_h[h] > 0]
                s_pos = [stds_by_h[h] for h in h_pos]
                if len(h_pos) >= 3:
                    log_h = np.log(np.array(h_pos))
                    log_s = np.log(np.array(s_pos))
                    beta_slope, beta_intercept = ols_fit(log_h, log_s)
                else:
                    beta_slope = float("nan")

                std_strs = " | ".join(
                    f"{stds_by_h.get(h, float('nan')):.6f}" for h in hubbles
                )
                md(f"| {d} | {N} | {obs_name} | {std_strs} | "
                   f"{rho_val:+.3f} | {p_val:.2e} | {beta_slope:.3f} |")

                inter_real_results.append({
                    "d": d, "N": N, "obs": obs_name,
                    "stds": stds_by_h, "rho": rho_val, "p": p_val, "beta": beta_slope,
                })

    md("")

    # -----------------------------------------------------------------------
    # Q2: Variance-based R estimator — V̂ = group_std² as H² proxy
    # -----------------------------------------------------------------------
    md("## 2. Variance-Based R Estimator: V̂ = std²(obs|H)\n")
    md("For each (d, N), compute group-level V̂ = std(obs|H)² per H level,")
    md("then test R²(V̂, H²) and compare with T3's two-step approach.\n")
    md("")
    md("| d | N | Observable | R²(V̂, H²) | R²(V̂, R_dS) | α_eff(V̂) | T3 R²(R̂,H²) ref | Δ R² (T2−T3) |")
    md("|---|---|------------|-----------|-------------|-----------|-----------------|--------------|")

    for d in dims:
        for N in ns_by_d[d]:
            sub = df[(df["d"] == d) & (df["N"] == N)]
            for obs_name in obs_names_a:
                # Group-level std² per H
                h_arr = []
                v_hat_arr = []
                h2_arr = []
                r_ds_arr = []
                for h in hubbles_nonzero:
                    cell = sub[sub["hubble"] == h][obs_name].dropna()
                    if len(cell) >= 2:
                        v_hat = float(cell.std()) ** 2
                        h_arr.append(h)
                        v_hat_arr.append(v_hat)
                        h2_arr.append(h ** 2)
                        r_ds_arr.append(d * (d - 1) * h ** 2)

                if len(h_arr) >= 3:
                    h_arr_np = np.array(h_arr)
                    v_hat_np = np.array(v_hat_arr)
                    h2_np = np.array(h2_arr)
                    r_ds_np = np.array(r_ds_arr)

                    r2_h2 = ols_r2(h2_np, v_hat_np)
                    r2_rds = ols_r2(r_ds_np, v_hat_np)
                    ag = alpha_grid(h_arr_np, v_hat_np)
                    a_eff, _ = best_alpha(ag)
                else:
                    r2_h2 = float("nan")
                    r2_rds = float("nan")
                    a_eff = float("nan")

                # T3 reference (approximate; from T3 S1 global OLS b1_std)
                t3_ref = "—"
                delta = "—"
                md(f"| {d} | {N} | {obs_name} | {r2_h2:.4f} | {r2_rds:.4f} | "
                   f"{a_eff:.2f} | {t3_ref} | {delta} |")

    md("")

    # -----------------------------------------------------------------------
    # Q3: Intra-realization (sub-patch) variance vs H²
    # -----------------------------------------------------------------------
    md("## 3. Intra-Realization (Sub-Patch) Variance vs H²\n")
    md("Per-realization b1_std_patch_std and wmr_patch_std correlated with H².\n")
    md("")
    md("| d | N | Patch observable | ρ(patch_std, H²) | p-value | R²(patch_std, H²) | R²(patch_std², H²) | α_eff |")
    md("|---|---|-----------------|-------------------|---------|--------------------|--------------------|-------|")

    patch_obs_names = [
        ("b1_std_patch_std", "b1_std patches"),
        ("wmr_patch_std", "w_max_ratio patches"),
    ]

    for d in dims:
        for N in ns_by_d[d]:
            sub = df[(df["d"] == d) & (df["N"] == N)]
            sub_nz = sub[sub["hubble"] > 0].copy()

            for col, label in patch_obs_names:
                valid = sub_nz[col].dropna()
                if len(valid) < 10:
                    md(f"| {d} | {N} | {label} | — | — | — | — | — |")
                    continue

                mask = sub_nz[col].notna()
                x_h2 = sub_nz.loc[mask, "H2"].values
                y_ps = sub_nz.loc[mask, col].values

                rho_val, p_val = sp_stats.spearmanr(x_h2, y_ps)
                r2_ps = ols_r2(x_h2, y_ps)
                r2_ps2 = ols_r2(x_h2, y_ps ** 2)

                h_vals_nz = sub_nz.loc[mask, "hubble"].values
                ag = alpha_grid(h_vals_nz, y_ps)
                a_eff, _ = best_alpha(ag)

                md(f"| {d} | {N} | {label} | {rho_val:+.3f} | {p_val:.2e} | "
                   f"{r2_ps:.4f} | {r2_ps2:.4f} | {a_eff:.2f} |")

    md("")

    # -----------------------------------------------------------------------
    # Q4: N-scaling of inter-realization β
    # -----------------------------------------------------------------------
    md("## 4. N-Scaling of β (Inter-Realization)\n")
    md("Does the power-law exponent β (std ∝ H^β) converge with N?\n")
    md("")
    md("| Observable | N=128 β | N=256 β | N=512 β | N=1024 β | ρ(N, β) | Trend |")
    md("|------------|---------|---------|---------|----------|---------|-------|")

    for obs_name in obs_names_a:
        betas = []
        for N in ns_by_d[4]:
            entry = [e for e in inter_real_results
                     if e["d"] == 4 and e["N"] == N and e["obs"] == obs_name]
            if entry:
                betas.append(entry[0]["beta"])
            else:
                betas.append(float("nan"))

        valid_pairs = [(n, b) for n, b in zip(ns_by_d[4], betas) if not np.isnan(b)]
        if len(valid_pairs) >= 3:
            rho_nb, _ = sp_stats.spearmanr(
                [p[0] for p in valid_pairs],
                [p[1] for p in valid_pairs]
            )
        else:
            rho_nb = float("nan")

        trend = "↑" if rho_nb > 0.5 else ("↓" if rho_nb < -0.5 else "—")
        beta_strs = " | ".join(f"{b:.3f}" for b in betas)
        md(f"| {obs_name} | {beta_strs} | {rho_nb:+.2f} | {trend} |")

    md("")

    # -----------------------------------------------------------------------
    # Q5: Per-realization variance estimator
    # -----------------------------------------------------------------------
    md("## 5. Per-Realization Squared Deviation as H² Proxy\n")
    md("For each realization i, define δᵢ = (obs_i − mean(obs|H=0))². ")
    md("Does δ correlate with H²? This converts per-realization noise into a signal.\n")
    md("")
    md("| d | N | Observable | ρ(δ, H²) | p-value | R²(δ, H²) | α_eff(δ) |")
    md("|---|---|------------|-----------|---------|-----------|----------|")

    for d in dims:
        for N in ns_by_d[d]:
            sub = df[(df["d"] == d) & (df["N"] == N)]
            sub_nz = sub[sub["hubble"] > 0]

            for obs_name in obs_names_a:
                # Baseline: mean of obs at H=0
                h0_vals = sub[sub["hubble"] == 0][obs_name].dropna()
                if len(h0_vals) == 0:
                    md(f"| {d} | {N} | {obs_name} | — | — | — | — |")
                    continue
                baseline = float(h0_vals.mean())

                # δ = (obs - baseline)²
                obs_vals = sub_nz[obs_name].dropna()
                if len(obs_vals) < 10:
                    md(f"| {d} | {N} | {obs_name} | — | — | — | — |")
                    continue

                mask_nz = sub_nz[obs_name].notna()
                delta_sq = (sub_nz.loc[mask_nz, obs_name].values - baseline) ** 2
                h2_vals = sub_nz.loc[mask_nz, "H2"].values
                h_vals_raw = sub_nz.loc[mask_nz, "hubble"].values

                rho_val, p_val = sp_stats.spearmanr(h2_vals, delta_sq)
                r2_val = ols_r2(h2_vals, delta_sq)

                ag = alpha_grid(h_vals_raw, delta_sq)
                a_eff, _ = best_alpha(ag)

                md(f"| {d} | {N} | {obs_name} | {rho_val:+.3f} | {p_val:.2e} | "
                   f"{r2_val:.4f} | {a_eff:.2f} |")

    md("")

    # -----------------------------------------------------------------------
    # Q6: Critical comparison — variance channel vs T3 two-step
    # -----------------------------------------------------------------------
    md("## 6. Critical Comparison: Variance Channel vs T3 Two-Step\n")
    md("| Channel | Method | Best R²(H²) at N=1024 | α_eff | Independent? |")
    md("|---------|--------|----------------------|-------|-------------|")
    md("| (To be filled after analysis) | | | | |")
    md("")

    # -----------------------------------------------------------------------
    # Overall Assessment
    # -----------------------------------------------------------------------
    md("## 7. Overall Assessment\n")
    md("| Question | Answer | Evidence |")
    md("|----------|--------|----------|")

    # Collect summary
    # Q1: std monotonic with H?
    q1_results = [(e["obs"], e["N"], e["rho"]) for e in inter_real_results if e["d"] == 4]
    q1_positive = sum(1 for _, _, r in q1_results if r > 0.5)
    q1_total = len(q1_results)
    md(f"| Q1: std(obs|H) monotone with H? | "
       f"{'✅' if q1_positive > q1_total * 0.7 else '⚠️'} {q1_positive}/{q1_total} positive | "
       f"ρ > 0.5 criterion |")

    # Q2: β ≈ 1?
    betas_b1 = [e["beta"] for e in inter_real_results
                if e["d"] == 4 and e["obs"] == "b1_std" and not np.isnan(e["beta"])]
    mean_beta_b1 = float(np.mean(betas_b1)) if betas_b1 else float("nan")
    md(f"| Q2: β ≈ 1 (std ∝ H^β)? | "
       f"{'✅' if 0.7 < mean_beta_b1 < 1.3 else '❌'} mean β(b1_std) = {mean_beta_b1:.3f} | "
       f"log-log fit on H>0 |")

    md("| Q3: V̂ useful as R estimator? | (see §2 table) | R²(V̂, H²) |")
    md("| Q4: Pattern strengthens with N? | (see §4 table) | ρ(N, β) |")
    md("| Q5: Independent of T3? | (see §6) | Cross-correlation |")
    md("")

    # -----------------------------------------------------------------------
    # Verdict (placeholder — to be refined after inspection)
    # -----------------------------------------------------------------------
    md("## 8. Physical Interpretation and Verdict\n")
    md("*(To be written after data inspection.)*")
    md("")
    md("---\n")
    md(f"*Generated by conjecture_e_cross_scale_variance.py*")
    md(f"*{len(df)} realizations (d=4, N=128/256/512/1024, H=0–2, {reps} reps), runtime {total_time:.1f}s*")

    # Write markdown
    md_path = out_dir / "conjecture_e_cross_scale_variance.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))
    print(f"Report saved: {md_path}")


if __name__ == "__main__":
    main()
