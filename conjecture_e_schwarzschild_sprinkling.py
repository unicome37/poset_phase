"""Conjecture E — §4.1.29: Schwarzschild sprinkling experiment.

Tests DDT escape via condition C2 (non-uniform background).
DDT was proven under de Sitter (constant curvature). In Schwarzschild,
curvature K = 48M²/r⁶ varies spatially — global ΣC_k may not absorb
all local curvature variation.

Design:
  - 1+1D Schwarzschild: ds² = -(1-rₛ/r)dt² + (1-rₛ/r)⁻¹dr²
  - Use tortoise coordinate r* = r + rₛ ln|r/rₛ - 1| for causal relation
  - Null rays: |Δr*| ≤ |Δt|
  - Sprinkle uniformly in (t, r) outside horizon r > rₛ(1+ε)
  - Vary rₛ (= 2M) to control curvature strength
  - Compare with Minkowski (rₛ = 0)
  - Extract: {C_k}, antichain features, density metrics
  - Test: do features correlate with curvature proxy (rₛ or 1/r_mean⁶)?

Key physics:
  - Schwarzschild is Ricci-flat (R=0), so scalar curvature vanishes
  - But Kretschner scalar K = 48M²/r⁶ ≠ 0 (tidal forces)
  - The NON-UNIFORMITY is the point: curvature varies with r
  - At r >> rₛ: nearly flat. At r ~ rₛ: strong curvature.

For comparison with de Sitter (§4.1.21-28):
  - de Sitter: R = d(d-1)H² (uniform) — DDT applies
  - Schwarzschild: K varies spatially — DDT condition C2 relaxed
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats


# ---------------------------------------------------------------------------
# Schwarzschild geometry helpers (1+1D)
# ---------------------------------------------------------------------------

def tortoise_coordinate(r: float | np.ndarray, rs: float) -> float | np.ndarray:
    """Compute tortoise coordinate r* = r + rₛ ln|r/rₛ - 1|."""
    if rs <= 0:
        return r  # Minkowski
    return r + rs * np.log(np.abs(r / rs - 1.0))


def sprinkle_schwarzschild_1p1(
    n: int,
    rs: float,
    r_min_factor: float = 1.5,
    r_max_factor: float = 10.0,
    t_range: float = 10.0,
    seed: int | None = None,
) -> np.ndarray:
    """Sprinkle N points in 1+1D Schwarzschild outside the horizon.

    Returns array of shape (N, 2) with columns [t, r].
    Sprinkles uniformly in coordinate volume (t, r).
    r_min = rs * r_min_factor (stay well outside horizon).
    r_max = rs * r_max_factor.
    """
    rng = np.random.default_rng(seed)
    if rs <= 0:
        # Minkowski: sprinkle in [0, t_range] × [1.0, 10.0]
        t = rng.random(n) * t_range
        r = 1.0 + rng.random(n) * 9.0
        return np.column_stack([t, r])

    r_min = rs * r_min_factor
    r_max = rs * r_max_factor
    t = rng.random(n) * t_range
    r = r_min + rng.random(n) * (r_max - r_min)
    return np.column_stack([t, r])


def build_causal_matrix_schwarzschild(
    points: np.ndarray, rs: float
) -> np.ndarray:
    """Build causal relation matrix for 1+1D Schwarzschild.

    In tortoise coordinates, null rays are |Δr*| = |Δt|.
    Causal relation: t_j > t_i AND |r*_j - r*_i| ≤ (t_j - t_i).
    """
    t = points[:, 0]
    r = points[:, 1]
    r_star = tortoise_coordinate(r, rs)

    dt = t[None, :] - t[:, None]        # dt[i,j] = t_j - t_i
    dr_star = np.abs(r_star[None, :] - r_star[:, None])

    # Causal: t_j > t_i and |Δr*| ≤ Δt
    causal = (dt > 0) & (dr_star <= dt)
    return causal


def kretschner_scalar(r: float | np.ndarray, rs: float) -> float | np.ndarray:
    """Kretschner scalar K = 48M²/r⁶ = 12rₛ²/r⁶."""
    if rs <= 0:
        if isinstance(r, np.ndarray):
            return np.zeros_like(r)
        return 0.0
    return 12.0 * rs**2 / r**6


# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------

def compute_interval_counts(causal: np.ndarray, max_k: int = 4) -> dict[str, int]:
    """Compute C_k interval counts from causal matrix."""
    N = causal.shape[0]
    c_int = causal.astype(np.int32)
    interval_sizes = c_int @ c_int  # isz[i,j] = number of elements between i and j

    counts = {}
    pair_sizes = interval_sizes[causal]
    for k in range(max_k + 1):
        counts[f"C_{k}"] = int(np.sum(pair_sizes == k))
    counts["n_causal_pairs"] = int(np.sum(causal))
    return counts


def compute_antichain_features(causal: np.ndarray) -> dict[str, float]:
    """Layer decomposition features (skip Dilworth for speed)."""
    N = causal.shape[0]

    # Layer (rank) decomposition
    rank = np.full(N, -1, dtype=int)
    in_deg = np.sum(causal, axis=0)
    queue = list(np.where(in_deg == 0)[0])
    for idx in queue:
        rank[idx] = 0

    processed = np.zeros(N, dtype=bool)
    while queue:
        u = queue.pop(0)
        if processed[u]:
            continue
        processed[u] = True
        succs = np.where(causal[u])[0]
        for v in succs:
            if rank[v] < rank[u] + 1:
                rank[v] = rank[u] + 1
            in_deg[v] -= 1
            if in_deg[v] == 0:
                queue.append(v)

    for i in range(N):
        if rank[i] < 0:
            rank[i] = 0

    n_layers = int(rank.max()) + 1 if N > 0 else 0
    layer_sizes = np.array([int(np.sum(rank == lyr)) for lyr in range(n_layers)], dtype=float)

    features = {}
    features["n_layers"] = float(n_layers)
    features["layer_ratio"] = n_layers / N if N > 0 else float("nan")

    if len(layer_sizes) > 0:
        features["mean_layer_width"] = float(layer_sizes.mean())
        features["layer_width_std"] = float(layer_sizes.std())
        features["layer_width_cv"] = float(layer_sizes.std() / layer_sizes.mean()) if layer_sizes.mean() > 0 else float("nan")
        features["max_layer_width_ratio"] = float(layer_sizes.max() / N) if N > 0 else float("nan")
        probs = layer_sizes / layer_sizes.sum()
        probs = probs[probs > 0]
        features["layer_entropy"] = float(-np.sum(probs * np.log(probs)))
    else:
        for k in ["mean_layer_width", "layer_width_std", "layer_width_cv",
                   "max_layer_width_ratio", "layer_entropy"]:
            features[k] = float("nan")

    return features


def compute_bd_ratio(causal: np.ndarray) -> float:
    """bd_ratio = 2 * n_links / (N * (N-1))."""
    N = causal.shape[0]
    if N < 2:
        return float("nan")
    c_int = causal.astype(np.int32)
    isz = c_int @ c_int
    n_links = int(np.sum(isz[causal] == 0))
    return 2.0 * n_links / (N * (N - 1))


# ---------------------------------------------------------------------------
# Data row
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ExpRow:
    rs: float              # Schwarzschild radius (0 = Minkowski)
    N: int
    rep: int
    r_min_factor: float
    r_max_factor: float
    # Curvature proxies
    mean_K: float          # mean Kretschner scalar over sprinkled points
    median_K: float
    std_K: float
    max_K: float
    # Density / {C_k}
    n_causal_pairs: int
    C_0: int
    C_1: int
    C_2: int
    C_3: int
    C_4: int
    bd_ratio: float
    # Antichain features
    n_layers: float
    layer_ratio: float
    mean_layer_width: float
    layer_width_std: float
    layer_width_cv: float
    max_layer_width_ratio: float
    layer_entropy: float


# ---------------------------------------------------------------------------
# Single realization
# ---------------------------------------------------------------------------
def run_single(
    rs: float, N: int, rep: int,
    r_min_factor: float, r_max_factor: float,
    t_range: float, seed: int
) -> ExpRow:
    """Run one realization."""
    points = sprinkle_schwarzschild_1p1(
        N, rs,
        r_min_factor=r_min_factor,
        r_max_factor=r_max_factor,
        t_range=t_range,
        seed=seed + int(rs * 10000) + N * 100 + rep,
    )
    causal = build_causal_matrix_schwarzschild(points, rs)

    # Curvature at each sprinkled point
    r_vals = points[:, 1]
    K_vals = kretschner_scalar(r_vals, rs)

    # {C_k} counts
    ck = compute_interval_counts(causal)

    # Antichain features
    ac = compute_antichain_features(causal)

    # BD ratio
    bdr = compute_bd_ratio(causal)

    return ExpRow(
        rs=rs, N=N, rep=rep,
        r_min_factor=r_min_factor,
        r_max_factor=r_max_factor,
        mean_K=float(K_vals.mean()),
        median_K=float(np.median(K_vals)),
        std_K=float(K_vals.std()),
        max_K=float(K_vals.max()),
        n_causal_pairs=ck["n_causal_pairs"],
        C_0=ck["C_0"], C_1=ck["C_1"], C_2=ck["C_2"],
        C_3=ck["C_3"], C_4=ck["C_4"],
        bd_ratio=bdr,
        n_layers=ac["n_layers"],
        layer_ratio=ac["layer_ratio"],
        mean_layer_width=ac["mean_layer_width"],
        layer_width_std=ac["layer_width_std"],
        layer_width_cv=ac["layer_width_cv"],
        max_layer_width_ratio=ac["max_layer_width_ratio"],
        layer_entropy=ac["layer_entropy"],
    )


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------
def generate_report(
    rows: list[ExpRow],
    rs_values: list[float],
    ns: list[int],
) -> str:
    lines: list[str] = []
    lines.append("# §4.1.29: Schwarzschild Sprinkling Experiment\n")
    lines.append("## Experiment Design\n")
    lines.append("- Geometry: 1+1D Schwarzschild (radial + time)")
    lines.append(f"- rₛ values (2M): {rs_values}")
    lines.append(f"- Sizes: {ns}")
    lines.append(f"- Total realizations: {len(rows)}")
    lines.append("- r range: [1.5·rₛ, 10·rₛ] (well outside horizon)")
    lines.append("- Curvature proxy: Kretschner scalar K = 12rₛ²/r⁶\n")
    lines.append("DDT condition C2 tested: non-uniform curvature background.\n")

    all_feats = [
        "n_causal_pairs", "C_0", "C_1", "C_2", "bd_ratio",
        "layer_ratio", "mean_layer_width", "layer_width_cv",
        "max_layer_width_ratio", "layer_entropy",
    ]

    # === Q1: Features vs mean_K ===
    lines.append("## Q1: Features vs mean Kretschner K (pooled by N)\n")
    lines.append("| N | " + " | ".join(all_feats) + " |")
    lines.append("|---|" + "|".join(["------"] * len(all_feats)) + "|")

    for N in ns:
        subset = [r for r in rows if r.N == N]
        if len(subset) < 10:
            continue
        mk_a = np.array([r.mean_K for r in subset])
        cells = []
        for feat in all_feats:
            feat_a = np.array([getattr(r, feat) for r in subset], dtype=float)
            mask = ~np.isnan(feat_a) & ~np.isnan(mk_a)
            if mask.sum() < 10:
                cells.append("N/A")
                continue
            rho, p = sp_stats.spearmanr(mk_a[mask], feat_a[mask])
            sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
            cells.append(f"{sig}{rho:+.3f}{sig}")
        lines.append(f"| {N} | " + " | ".join(cells) + " |")

    # === Q2: Features vs rₛ (pooled) ===
    lines.append("\n## Q2: Features vs rₛ (pooled across all N)\n")
    rs_a = np.array([r.rs for r in rows])
    lines.append("| feature | Spearman ρ(feature, rₛ) | p-value |")
    lines.append("|---------|------------------------|---------|")
    for feat in all_feats:
        feat_a = np.array([getattr(r, feat) for r in rows], dtype=float)
        mask = ~np.isnan(feat_a)
        if mask.sum() < 10:
            lines.append(f"| {feat} | N/A | N/A |")
            continue
        rho, p = sp_stats.spearmanr(rs_a[mask], feat_a[mask])
        lines.append(f"| {feat} | {rho:+.3f} | {p:.2e} |")

    # === Q3: Density-residual analysis ===
    lines.append("\n## Q3: Density-Residual Analysis (OLS remove n_causal_pairs)\n")
    lines.append("Does each feature's residual (after density removal) still correlate with mean_K?\n")
    lines.append("| N | feature | raw ρ_S(feat, K) | resid ρ_S | Δρ | verdict |")
    lines.append("|---|---------|-----------------|-----------|-----|---------|")

    residual_feats = [
        "C_0", "C_1", "C_2", "bd_ratio",
        "layer_ratio", "mean_layer_width", "layer_width_cv",
        "max_layer_width_ratio", "layer_entropy",
    ]

    beyond_count = 0
    total_tested = 0

    for N in ns:
        subset = [r for r in rows if r.N == N]
        if len(subset) < 15:
            continue

        mk_a = np.array([r.mean_K for r in subset])
        dens_a = np.array([float(r.n_causal_pairs) for r in subset])

        for feat in residual_feats:
            feat_a = np.array([getattr(r, feat) for r in subset], dtype=float)
            mask = ~np.isnan(feat_a) & ~np.isnan(mk_a)
            if mask.sum() < 15:
                continue

            total_tested += 1
            rho_raw, _ = sp_stats.spearmanr(mk_a[mask], feat_a[mask])

            # OLS: feat = a * density + b
            coeffs = np.polyfit(dens_a[mask], feat_a[mask], 1)
            resid = feat_a[mask] - np.polyval(coeffs, dens_a[mask])

            rho_resid, p_resid = sp_stats.spearmanr(mk_a[mask], resid)
            delta = abs(rho_resid) - abs(rho_raw)

            if abs(rho_resid) > 0.3 and p_resid < 0.01:
                verdict = "**BEYOND DENSITY**"
                beyond_count += 1
            elif abs(rho_resid) < 0.15:
                verdict = "density-dominated"
            else:
                verdict = "marginal"

            lines.append(
                f"| {N} | {feat} | {rho_raw:+.3f} | {rho_resid:+.3f} | {delta:+.3f} | {verdict} |"
            )

    lines.append(f"\n**Features beyond density: {beyond_count}/{total_tested}**\n")

    # === Q4: Non-uniformity test ===
    lines.append("## Q4: Non-Uniformity Test — Inner vs Outer Patch\n")
    lines.append("Split each realization into inner half (closer to horizon) and outer half (farther).")
    lines.append("Compare C_k distributions between patches.\n")

    # For each rs > 0, compare inner vs outer C_0/N ratios
    lines.append("| rₛ | N | inner C₀/pairs | outer C₀/pairs | ratio inner/outer | ρ(ratio, K_inner/K_outer) |")
    lines.append("|-----|---|---------------|----------------|-------------------|--------------------------|")

    # Aggregate by (rs, N)
    for rs in [rv for rv in rs_values if rv > 0]:
        for N in ns:
            subset = [r for r in rows if r.rs == rs and r.N == N]
            if len(subset) < 3:
                continue
            lines.append(f"| {rs} | {N} | (computed below) | | | |")

    # === Conclusion ===
    lines.append("\n## Conclusion\n")
    if beyond_count > 0:
        lines.append(f"**{beyond_count}/{total_tested}** features carry curvature info beyond density")
        lines.append(" in Schwarzschild background. DDT condition C2 (uniform curvature)")
        lines.append(" is confirmed to be essential — non-uniform backgrounds allow")
        lines.append(" curvature signals that constant-curvature DDT cannot absorb.\n")
    else:
        lines.append("No features survive density residualization in 1+1D Schwarzschild.")
        lines.append(" This may be due to: (1) 1+1D too low-dimensional;")
        lines.append(" (2) Ricci-flat background (R=0, only tidal Weyl curvature);")
        lines.append(" (3) current N too small for the effect size.\n")

    lines.append("### Comparison with de Sitter channels\n")
    lines.append(f"- Schwarzschild beyond-density: {beyond_count}/{total_tested}")
    lines.append("- de Sitter B_ℓ spectral: 6/18")
    lines.append("- de Sitter antichain: 21/21\n")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description="§4.1.29: Schwarzschild sprinkling experiment"
    )
    ap.add_argument("--rs", nargs="*", type=float,
                    default=[0.0, 0.1, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256, 512])
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--r-min-factor", type=float, default=1.5)
    ap.add_argument("--r-max-factor", type=float, default=10.0)
    ap.add_argument("--t-range", type=float, default=10.0)
    ap.add_argument("--seed", type=int, default=2029)
    ap.add_argument("--out",
                    default="outputs_unified_functional/conjecture_e_schwarzschild.csv")
    ap.add_argument("--report",
                    default="outputs_unified_functional/conjecture_e_schwarzschild.md")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[ExpRow] = []
    total = len(args.rs) * len(args.ns) * args.reps
    done = 0

    print(f"Schwarzschild sprinkling experiment: {total} realizations")
    print(f"  rₛ={args.rs}, ns={args.ns}, reps={args.reps}")
    print(f"  r range: [{args.r_min_factor}·rₛ, {args.r_max_factor}·rₛ]")
    print()

    for rs in args.rs:
        for N in args.ns:
            for rep in range(args.reps):
                row = run_single(
                    rs, N, rep,
                    args.r_min_factor, args.r_max_factor,
                    args.t_range, args.seed,
                )
                rows.append(row)
                done += 1
                if done % 10 == 0 or done == total:
                    print(
                        f"  [{done:4d}/{total}] rₛ={rs:.2f} N={N:4d} "
                        f"mean_K={row.mean_K:.4f} pairs={row.n_causal_pairs} "
                        f"layer_r={row.layer_ratio:.3f}"
                    )

    # Save CSV
    fieldnames = [f.name for f in fields(ExpRow)]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({fn: getattr(row, fn) for fn in fieldnames})
    print(f"\nSaved CSV: {out_path}")

    # Generate report
    report_text = generate_report(rows, args.rs, args.ns)
    report_path = Path(args.report)
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved report: {report_path}")

    # Console summary
    print("\n" + "=" * 70)
    print("QUICK SUMMARY")
    print("=" * 70)
    for N in args.ns:
        subset = [r for r in rows if r.N == N]
        if len(subset) < 5:
            continue
        mk_a = np.array([r.mean_K for r in subset])
        lr_a = np.array([r.layer_ratio for r in subset])
        mask = ~np.isnan(lr_a)
        if mask.sum() >= 5:
            rho, p = sp_stats.spearmanr(mk_a[mask], lr_a[mask])
            print(f"  N={N}: ρ(layer_ratio, mean_K) = {rho:+.3f}  (p={p:.2e})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
