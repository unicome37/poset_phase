"""Conjecture E — §4.1.28: Antichain structure experiment.

Tests whether antichain (transverse) statistics of causal sets carry
curvature information beyond the interval-count family {C_k}.

Antichain statistics measure the "width" of the poset — how many elements
are mutually unrelated (spacelike-separated). This is orthogonal to the
"longitudinal" interval counts that DDT has shown to be density-degenerate.

Metrics:
  1. w_max / N  — maximum antichain width ratio (Dilworth's theorem)
  2. n_layers   — number of layers in the longest chain decomposition
  3. layer_ratio = n_layers / N
  4. mean_layer_width — mean elements per layer
  5. layer_width_std  — std of elements per layer
  6. layer_width_cv   — coefficient of variation of layer widths
  7. max_layer_width / N — widest single layer normalized

All metrics are computed from the causal matrix of de Sitter sprinklings.
Density-residual analysis (Q5) tests whether signals survive after
OLS removal of n_causal_pairs.

References:
  - Dilworth's theorem: width = N - maximum matching in comparability graph
  - Layer decomposition: greedily assign elements to antichains (rank function)
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import build_causal_matrix_de_sitter


# ---------------------------------------------------------------------------
# Antichain and layer metrics
# ---------------------------------------------------------------------------

def compute_transitive_closure(causal: np.ndarray) -> np.ndarray:
    """Warshall-style boolean transitive closure."""
    c = causal.copy()
    n = c.shape[0]
    for k in range(n):
        c |= c[:, [k]] & c[[k], :]
    np.fill_diagonal(c, False)
    return c


def antichain_width_from_closure(closure: np.ndarray) -> int:
    """Maximum antichain size via Dilworth reduction to bipartite matching.

    width = N - max_matching_size in the comparability bipartite graph.
    Uses Hopcroft-Karp style augmenting paths.
    """
    n = closure.shape[0]
    # Build adjacency list: i -> {j : i < j}
    adj = [list(np.where(closure[i])[0]) for i in range(n)]
    match_r = [-1] * n

    def try_augment(u: int, seen: list[bool]) -> bool:
        for v in adj[u]:
            if seen[v]:
                continue
            seen[v] = True
            if match_r[v] == -1 or try_augment(match_r[v], seen):
                match_r[v] = u
                return True
        return False

    match_size = 0
    for u in range(n):
        seen = [False] * n
        if try_augment(u, seen):
            match_size += 1
    return int(n - match_size)


def compute_layer_decomposition(causal: np.ndarray) -> list[int]:
    """Compute layer (rank) decomposition of a causal set.

    Each element gets a rank = 1 + max rank of its predecessors.
    Elements with no predecessors get rank 0.
    Returns list of layer sizes.
    """
    n = causal.shape[0]
    # in-degree from causal[j, i] = True means j ≺ i
    # predecessors of i: all j where causal[j, i] is True
    rank = np.full(n, -1, dtype=int)

    # Topological sort via in-degree
    in_deg = np.sum(causal, axis=0)  # causal[j,i] = j≺i, so col sum = in-degree
    queue = list(np.where(in_deg == 0)[0])
    for idx in queue:
        rank[idx] = 0

    # BFS-like topological processing
    processed = np.zeros(n, dtype=bool)
    order = []

    while queue:
        u = queue.pop(0)
        if processed[u]:
            continue
        processed[u] = True
        order.append(u)
        # u's successors: causal[u, v] = True
        succs = np.where(causal[u])[0]
        for v in succs:
            if rank[v] < rank[u] + 1:
                rank[v] = rank[u] + 1
            in_deg[v] -= 1
            if in_deg[v] == 0:
                queue.append(v)

    # Handle any unprocessed (shouldn't happen for a DAG)
    for i in range(n):
        if rank[i] < 0:
            rank[i] = 0

    # Count elements per layer
    n_layers = int(rank.max()) + 1 if n > 0 else 0
    layer_sizes = []
    for lyr in range(n_layers):
        layer_sizes.append(int(np.sum(rank == lyr)))

    return layer_sizes


def compute_antichain_features(
    causal: np.ndarray,
    compute_dilworth: bool = True,
) -> dict[str, float]:
    """Extract antichain/layer features from a causal matrix."""
    n = causal.shape[0]
    features: dict[str, float] = {}

    # Layer decomposition (fast — O(N²) via topological sort)
    layer_sizes = compute_layer_decomposition(causal)
    n_layers = len(layer_sizes)

    features["n_layers"] = float(n_layers)
    features["layer_ratio"] = n_layers / n if n > 0 else float("nan")

    if layer_sizes:
        ls = np.array(layer_sizes, dtype=float)
        features["mean_layer_width"] = float(ls.mean())
        features["layer_width_std"] = float(ls.std())
        features["layer_width_cv"] = float(ls.std() / ls.mean()) if ls.mean() > 0 else float("nan")
        features["max_layer_width"] = float(ls.max())
        features["max_layer_width_ratio"] = float(ls.max() / n) if n > 0 else float("nan")
        features["min_layer_width"] = float(ls.min())
        # Entropy of layer width distribution
        probs = ls / ls.sum()
        probs = probs[probs > 0]
        features["layer_entropy"] = float(-np.sum(probs * np.log(probs)))
    else:
        for k in ["mean_layer_width", "layer_width_std", "layer_width_cv",
                   "max_layer_width", "max_layer_width_ratio", "min_layer_width",
                   "layer_entropy"]:
            features[k] = float("nan")

    # Dilworth width (expensive — O(N^2.5) matching)
    if compute_dilworth and n <= 600:
        closure = compute_transitive_closure(causal)
        w = antichain_width_from_closure(closure)
        features["w_max"] = float(w)
        features["w_max_ratio"] = float(w / n) if n > 0 else float("nan")
    else:
        features["w_max"] = float("nan")
        features["w_max_ratio"] = float("nan")

    return features


# ---------------------------------------------------------------------------
# Data row
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
    # Antichain features
    w_max: float
    w_max_ratio: float
    n_layers: float
    layer_ratio: float
    mean_layer_width: float
    layer_width_std: float
    layer_width_cv: float
    max_layer_width: float
    max_layer_width_ratio: float
    min_layer_width: float
    layer_entropy: float


# ---------------------------------------------------------------------------
# Single realization
# ---------------------------------------------------------------------------
def run_single(
    d: int, N: int, hubble: float, rep: int, seed: int
) -> ExpRow:
    """Run one realization: sprinkle, build causal matrix, extract antichain features."""
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 1000),
    )
    causal = build_causal_matrix_de_sitter(points, hubble)

    n_causal_pairs = int(np.sum(causal))
    R_dS = d * (d - 1) * hubble ** 2
    H2 = hubble ** 2

    feats = compute_antichain_features(causal, compute_dilworth=(N <= 512))

    return ExpRow(
        d=d, N=N, hubble=hubble, rep=rep,
        R_dS=R_dS, H2=H2,
        n_causal_pairs=n_causal_pairs,
        w_max=feats["w_max"],
        w_max_ratio=feats["w_max_ratio"],
        n_layers=feats["n_layers"],
        layer_ratio=feats["layer_ratio"],
        mean_layer_width=feats["mean_layer_width"],
        layer_width_std=feats["layer_width_std"],
        layer_width_cv=feats["layer_width_cv"],
        max_layer_width=feats["max_layer_width"],
        max_layer_width_ratio=feats["max_layer_width_ratio"],
        min_layer_width=feats["min_layer_width"],
        layer_entropy=feats["layer_entropy"],
    )


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------
def generate_report(
    rows: list[ExpRow],
    dims: list[int],
    ns: list[int],
    hubbles: list[float],
) -> str:
    lines: list[str] = []
    lines.append("# §4.1.28: Antichain Structure Experiment\n")
    lines.append("## Experiment Design\n")
    lines.append(f"- Dimensions: {dims}")
    lines.append(f"- Sizes: {ns}")
    lines.append(f"- Hubble values: {hubbles}")
    lines.append(f"- Total realizations: {len(rows)}\n")
    lines.append("Antichain (transverse) statistics measure the 'width' of the poset,")
    lines.append("orthogonal to the 'longitudinal' interval counts {C_k} closed by DDT.\n")

    # === Q1: Raw correlations of antichain features with H² ===
    lines.append("## Q1: Antichain Features vs H² (Spearman, pooled by d)\n")
    ac_feats = [
        "w_max_ratio", "layer_ratio", "mean_layer_width",
        "layer_width_std", "layer_width_cv", "max_layer_width_ratio",
        "layer_entropy",
    ]
    lines.append("| d | " + " | ".join(ac_feats) + " |")
    lines.append("|---|" + "|".join(["------"] * len(ac_feats)) + "|")

    for d in dims:
        vals = {}
        h2_arr = []
        for feat in ac_feats:
            vals[feat] = []
        for r in rows:
            if r.d != d:
                continue
            h2_arr.append(r.H2)
            for feat in ac_feats:
                vals[feat].append(getattr(r, feat))
        if len(h2_arr) < 10:
            continue
        h2_a = np.array(h2_arr)
        cells = []
        for feat in ac_feats:
            feat_a = np.array(vals[feat])
            mask = ~np.isnan(feat_a)
            if mask.sum() < 10:
                cells.append("N/A")
                continue
            rho, p = sp_stats.spearmanr(h2_a[mask], feat_a[mask])
            sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
            cells.append(f"{sig}{rho:+.3f}{sig}")
        lines.append(f"| {d} | " + " | ".join(cells) + " |")

    # === Q2: Per-slice correlations ===
    lines.append("\n## Q2: Per-Slice Spearman (d, N) for Key Features\n")
    key_feats = ["w_max_ratio", "layer_ratio", "layer_width_cv", "max_layer_width_ratio", "layer_entropy"]
    lines.append("| d | N | " + " | ".join(key_feats) + " |")
    lines.append("|---|---|" + "|".join(["------"] * len(key_feats)) + "|")

    sig_count = 0
    total_slices = 0
    for d in dims:
        for N in ns:
            subset = [r for r in rows if r.d == d and r.N == N]
            if len(subset) < 10:
                continue
            total_slices += 1
            h2_a = np.array([r.H2 for r in subset])
            cells = []
            for feat in key_feats:
                feat_a = np.array([getattr(r, feat) for r in subset])
                mask = ~np.isnan(feat_a)
                if mask.sum() < 10:
                    cells.append("N/A")
                    continue
                rho, p = sp_stats.spearmanr(h2_a[mask], feat_a[mask])
                sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
                cells.append(f"{sig}{rho:+.3f}{sig}")
            lines.append(f"| {d} | {N} | " + " | ".join(cells) + " |")

    # === Q3: Density correlation check ===
    lines.append("\n## Q3: Antichain Features vs n_causal_pairs (Density)\n")
    lines.append("If |ρ| ≈ 1.0, the feature is density-dominated.\n")
    lines.append("| d | " + " | ".join(ac_feats) + " |")
    lines.append("|---|" + "|".join(["------"] * len(ac_feats)) + "|")

    for d in dims:
        h2_arr, dens_arr = [], []
        vals = {f: [] for f in ac_feats}
        for r in rows:
            if r.d != d:
                continue
            dens_arr.append(float(r.n_causal_pairs))
            for feat in ac_feats:
                vals[feat].append(getattr(r, feat))
        if len(dens_arr) < 10:
            continue
        dens_a = np.array(dens_arr)
        cells = []
        for feat in ac_feats:
            feat_a = np.array(vals[feat])
            mask = ~np.isnan(feat_a)
            if mask.sum() < 10:
                cells.append("N/A")
                continue
            rho, _ = sp_stats.spearmanr(dens_a[mask], feat_a[mask])
            cells.append(f"{rho:+.3f}")
        lines.append(f"| {d} | " + " | ".join(cells) + " |")

    # === Q5: Density-Residual Analysis ===
    lines.append("\n## Q5: Density-Residual Analysis\n")
    lines.append("After OLS-removing n_causal_pairs (density proxy ≈ ΣC_k),")
    lines.append("does each antichain feature's residual still correlate with H²?\n")
    lines.append("| d | feature | raw ρ_S | residual ρ_S | Δρ | verdict |")
    lines.append("|---|---------|--------|-------------|-----|---------|")

    density_residual_results: dict[tuple[int, str], tuple[float, float]] = {}

    for d in dims:
        for feat in ac_feats:
            h2_arr, feat_arr, dens_arr = [], [], []
            for r in rows:
                if r.d != d:
                    continue
                v = getattr(r, feat, None)
                if v is not None and not math.isnan(v) and not math.isnan(r.n_causal_pairs):
                    h2_arr.append(r.H2)
                    feat_arr.append(v)
                    dens_arr.append(float(r.n_causal_pairs))
            if len(h2_arr) < 15:
                continue

            h2_a = np.array(h2_arr)
            feat_a = np.array(feat_arr)
            dens_a = np.array(dens_arr)

            rho_raw, _ = sp_stats.spearmanr(h2_a, feat_a)

            # OLS: feat = a * density + b → residual
            coeffs_ols = np.polyfit(dens_a, feat_a, 1)
            resid = feat_a - np.polyval(coeffs_ols, dens_a)

            rho_resid, p_resid = sp_stats.spearmanr(h2_a, resid)
            delta = abs(rho_resid) - abs(rho_raw)

            if abs(rho_resid) > 0.3 and p_resid < 0.01:
                verdict = "**BEYOND DENSITY**"
            elif abs(rho_resid) < 0.15:
                verdict = "density-dominated"
            else:
                verdict = "marginal"

            density_residual_results[(d, feat)] = (rho_raw, rho_resid)
            lines.append(
                f"| {d} | {feat} | {rho_raw:+.3f} | {rho_resid:+.3f} | {delta:+.3f} | {verdict} |"
            )

    beyond_count = sum(1 for v in density_residual_results.values()
                       if abs(v[1]) > 0.3)
    total_tested = len(density_residual_results)
    lines.append(f"\n**Features with |residual ρ| > 0.3: {beyond_count}/{total_tested}**\n")

    if beyond_count > 0:
        lines.append("### Interpretation\n")
        lines.append("Antichain features that survive density removal carry ")
        lines.append("information **beyond** the {C_k} family. ")
        lines.append("This would validate DDT escape via condition C1 ")
        lines.append("(transverse structure ≠ interval counts).\n")
    else:
        lines.append("### Interpretation\n")
        lines.append("All antichain features collapse to density after residualization. ")
        lines.append("Transverse structure at these N values is still dominated by ")
        lines.append("the same density signal. This does NOT rule out antichain signals ")
        lines.append("at larger N or in non-uniform backgrounds.\n")

    # === Conclusion ===
    lines.append("\n## Conclusion\n")

    # Collect beyond-density summary
    spectral_sig = {}
    for d in dims:
        for feat in ac_feats:
            key = (d, feat)
            if key in density_residual_results:
                _, rho_r = density_residual_results[key]
                if abs(rho_r) > 0.3:
                    spectral_sig.setdefault(d, []).append((feat, rho_r))

    for d in dims:
        feats_found = spectral_sig.get(d, [])
        if feats_found:
            feat_str = ", ".join(f"{f}(ρ_resid={r:+.3f})" for f, r in feats_found)
            lines.append(f"- **d={d}**: {len(feats_found)} features BEYOND DENSITY: {feat_str}")
        else:
            lines.append(f"- d={d}: no antichain features survive density removal")

    lines.append("")
    if beyond_count >= 3:
        lines.append(f"**VERDICT: Antichain structure encodes curvature BEYOND {{C_k}}.**")
        lines.append(f" {beyond_count}/{total_tested} features survive density residualization.")
        lines.append(" DDT condition C1 escaped via transverse (antichain) observables.\n")
    elif beyond_count > 0:
        lines.append(f"**Partial signal**: {beyond_count}/{total_tested} features beyond density.")
        lines.append(" Some transverse features carry independent curvature info.\n")
    else:
        lines.append(f"**Antichain features density-dominated** at current N.")
        lines.append(" Transverse structure does not provide independent curvature signal")
        lines.append(" beyond what total density (ΣC_k) already captures.\n")

    # Comparison with B_ℓ spectral (§4.1.27)
    lines.append("### Comparison with B_ℓ Spectral (§4.1.27)\n")
    lines.append(f"- Antichain beyond-density: {beyond_count}/{total_tested}")
    lines.append("- B_ℓ spectral beyond-density: 6/18 (d=4: 4/6)")
    if beyond_count > 0:
        lines.append("\nBoth transverse and spectral channels provide independent curvature info.")
    else:
        lines.append("\nAntichain channel weaker than spectral — the operator-level")
        lines.append("structure of B_ℓ encodes more curvature than raw transverse widths.\n")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description="§4.1.28: Antichain structure experiment"
    )
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256, 512])
    ap.add_argument("--hubbles", nargs="*", type=float,
                    default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2028)
    ap.add_argument("--out",
                    default="outputs_unified_functional/conjecture_e_antichain_structure.csv")
    ap.add_argument("--report",
                    default="outputs_unified_functional/conjecture_e_antichain_structure.md")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[ExpRow] = []
    total = len(args.dims) * len(args.ns) * len(args.hubbles) * args.reps
    done = 0

    print(f"Antichain structure experiment: {total} realizations")
    print(f"  dims={args.dims}, ns={args.ns}, hubbles={args.hubbles}, reps={args.reps}")
    print()

    for d in args.dims:
        for N in args.ns:
            for hubble in args.hubbles:
                for rep in range(args.reps):
                    row = run_single(d, N, hubble, rep, args.seed)
                    rows.append(row)
                    done += 1
                    if done % 10 == 0 or done == total:
                        print(
                            f"  [{done:4d}/{total}] d={d} N={N:4d} H={hubble:.2f} "
                            f"w_max_r={row.w_max_ratio:.3f} "
                            f"layer_r={row.layer_ratio:.3f} "
                            f"pairs={row.n_causal_pairs}"
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
    report_text = generate_report(rows, args.dims, args.ns, args.hubbles)
    report_path = Path(args.report)
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved report: {report_path}")

    # Console summary
    print("\n" + "=" * 70)
    print("QUICK SUMMARY")
    print("=" * 70)
    for d in args.dims:
        h2_all, wr_all = [], []
        for r in rows:
            if r.d == d and not math.isnan(r.w_max_ratio):
                h2_all.append(r.H2)
                wr_all.append(r.w_max_ratio)
        if len(h2_all) >= 5:
            rho, p = sp_stats.spearmanr(h2_all, wr_all)
            print(f"  d={d}: ρ(w_max_ratio, H²) = {rho:+.3f}  (p={p:.2e})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
