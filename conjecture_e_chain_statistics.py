"""Conjecture E — §4.1.30: Causal chain statistics experiment.

Tests DDT escape via path I.4: k-chain counts N_k.

A k-chain is a totally ordered subset x_1 ≺ x_2 ≺ ... ≺ x_k.
Chain counts N_k encode "longitudinal path structure" — how many
directed paths of length k exist in the causal set.

Key relation to {C_k}:
  - C_k = # pairs with exactly k elements between them (interval richness)
  - N_k = total # of k-chains = 1^T · causal^{k-1} · 1
  - N_2 = n_causal_pairs (= ΣC_k + links)
  - N_3 counts 3-element chains — encodes path correlations beyond pairwise

Brightwell-Gregory connection:
  Interval counts and chain counts are related via inclusion-exclusion,
  but chain counts encode higher-order path correlations that {C_k} alone
  may not capture.

DDT prediction:
  If chain counts are merely functions of {C_k}, they should be
  density-degenerate. If they encode independent structure, they
  should survive density residualization.

Design:
  - de Sitter sprinklings: d=2,3,4; N=128,256,512; H=0–2; 8 reps
  - Extract: N_2, N_3, N_4, N_5, longest chain, chain ratios, chain entropy
  - Q5: OLS remove n_causal_pairs, test residual correlation with H²
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
# Chain statistics
# ---------------------------------------------------------------------------

def count_chains(causal: np.ndarray, max_k: int = 5) -> dict[str, float]:
    """Count k-chains via matrix powers of the causal matrix.

    N_k = 1^T · C^{k-1} · 1 = sum of all entries of C^{k-1}.
    C is the strict causal order (boolean, but compute as int for matmul).
    """
    N = causal.shape[0]
    c = causal.astype(np.float64)
    ones = np.ones(N, dtype=np.float64)

    counts: dict[str, float] = {}

    # N_1 = N (trivial 1-chains = individual elements)
    counts["N_1"] = float(N)

    # N_2 = sum(C) = number of comparable pairs
    counts["N_2"] = float(c.sum())

    # Higher chains via matrix powers
    power = c.copy()
    for k in range(3, max_k + 1):
        power = power @ c
        counts[f"N_{k}"] = float(power.sum())

    return counts


def longest_chain_length(causal: np.ndarray) -> int:
    """Longest chain (= height + 1) via longest path in DAG.

    Uses topological ordering + dynamic programming.
    """
    N = causal.shape[0]
    if N == 0:
        return 0

    # In-degree for topological sort
    in_deg = causal.sum(axis=0).astype(int)
    queue = list(np.where(in_deg == 0)[0])

    dp = np.ones(N, dtype=int)  # dp[i] = longest chain ending at i
    topo_order: list[int] = []

    while queue:
        u = queue.pop(0)
        topo_order.append(u)
        succs = np.where(causal[u])[0]
        for v in succs:
            dp[v] = max(dp[v], dp[u] + 1)
            in_deg[v] -= 1
            if in_deg[v] == 0:
                queue.append(v)

    return int(dp.max())


def compute_chain_features(causal: np.ndarray) -> dict[str, float]:
    """Extract all chain-based features."""
    N = causal.shape[0]
    features: dict[str, float] = {}

    # Chain counts
    ck = count_chains(causal, max_k=5)
    for key, val in ck.items():
        features[key] = val

    # Longest chain
    lc = longest_chain_length(causal)
    features["longest_chain"] = float(lc)
    features["longest_chain_ratio"] = lc / N if N > 0 else float("nan")

    # Chain ratios (higher-order vs lower-order)
    n2 = ck["N_2"]
    n3 = ck["N_3"]
    n4 = ck["N_4"]
    n5 = ck["N_5"]

    features["N3_over_N2"] = n3 / n2 if n2 > 0 else float("nan")
    features["N4_over_N3"] = n4 / n3 if n3 > 0 else float("nan")
    features["N5_over_N4"] = n5 / n4 if n4 > 0 else float("nan")
    features["N4_over_N2sq"] = n4 / (n2 ** 2) if n2 > 0 else float("nan")

    # Chain "entropy" — distribution of chain counts
    chain_vec = np.array([ck[f"N_{k}"] for k in range(2, 6)], dtype=float)
    chain_vec = chain_vec[chain_vec > 0]
    if len(chain_vec) > 0:
        p = chain_vec / chain_vec.sum()
        features["chain_entropy"] = float(-np.sum(p * np.log(p)))
    else:
        features["chain_entropy"] = float("nan")

    # Normalized chain counts
    features["N2_per_N"] = n2 / N if N > 0 else float("nan")
    features["N3_per_N"] = n3 / N if N > 0 else float("nan")
    features["N4_per_N"] = n4 / N if N > 0 else float("nan")

    return features


def compute_interval_counts(causal: np.ndarray) -> dict[str, int]:
    """Compute C_k interval counts from causal matrix."""
    c_int = causal.astype(np.int32)
    interval_sizes = c_int @ c_int
    pair_sizes = interval_sizes[causal]

    counts: dict[str, int] = {}
    for k in range(5):
        counts[f"C_{k}"] = int(np.sum(pair_sizes == k))
    counts["n_causal_pairs"] = int(np.sum(causal))
    return counts


# ---------------------------------------------------------------------------
# Data row
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ExpRow:
    d: int
    N: int
    H: float
    H_sq: float
    rep: int
    # Density
    n_causal_pairs: int
    C_0: int
    C_1: int
    C_2: int
    # Chain counts
    N_2: float
    N_3: float
    N_4: float
    N_5: float
    # Chain features
    longest_chain: float
    longest_chain_ratio: float
    N3_over_N2: float
    N4_over_N3: float
    N5_over_N4: float
    N4_over_N2sq: float
    chain_entropy: float
    N2_per_N: float
    N3_per_N: float
    N4_per_N: float


# ---------------------------------------------------------------------------
# Single realization
# ---------------------------------------------------------------------------
def run_single(
    d: int, N: int, H: float, rep: int, seed: int
) -> ExpRow:
    d_spatial = d - 1
    points = sprinkle_de_sitter_like_diamond(N, d_spatial, H, seed=seed)
    causal = build_causal_matrix_de_sitter(points, H)

    # Interval counts
    ic = compute_interval_counts(causal)

    # Chain features
    cf = compute_chain_features(causal)

    return ExpRow(
        d=d, N=N, H=H, H_sq=H * H, rep=rep,
        n_causal_pairs=ic["n_causal_pairs"],
        C_0=ic["C_0"], C_1=ic["C_1"], C_2=ic["C_2"],
        N_2=cf["N_2"], N_3=cf["N_3"], N_4=cf["N_4"], N_5=cf["N_5"],
        longest_chain=cf["longest_chain"],
        longest_chain_ratio=cf["longest_chain_ratio"],
        N3_over_N2=cf["N3_over_N2"],
        N4_over_N3=cf["N4_over_N3"],
        N5_over_N4=cf["N5_over_N4"],
        N4_over_N2sq=cf["N4_over_N2sq"],
        chain_entropy=cf["chain_entropy"],
        N2_per_N=cf["N2_per_N"],
        N3_per_N=cf["N3_per_N"],
        N4_per_N=cf["N4_per_N"],
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
    lines.append("# §4.1.30: Causal Chain Statistics Experiment\n")
    lines.append("## Experiment Design\n")
    lines.append(f"- Dimensions: {dims}")
    lines.append(f"- Sizes: {ns}")
    lines.append(f"- Hubble values: {hubbles}")
    lines.append(f"- Total realizations: {len(rows)}")
    lines.append("- Chain features: N_2..N_5, longest chain, chain ratios, chain entropy")
    lines.append("- DDT path I.4: do k-chain counts carry curvature info beyond {C_k}?\n")

    chain_feats = [
        "N_3", "N_4", "N_5",
        "longest_chain", "longest_chain_ratio",
        "N3_over_N2", "N4_over_N3", "N5_over_N4",
        "N4_over_N2sq", "chain_entropy",
        "N3_per_N", "N4_per_N",
    ]
    # NOTE: N2_per_N excluded — it equals n_causal_pairs/N, trivially
    # proportional to density within each (d,N) slice. OLS residuals
    # are numerical noise from perfect collinearity, not real signal.

    # === Q1: Raw correlation with H² ===
    lines.append("## Q1: Chain Features vs H² (pooled by d)\n")
    lines.append("| d | " + " | ".join(chain_feats) + " |")
    lines.append("|---|" + "|".join(["------"] * len(chain_feats)) + "|")

    for d in dims:
        subset = [r for r in rows if r.d == d]
        if len(subset) < 10:
            continue
        h2_a = np.array([r.H_sq for r in subset])
        cells = []
        for feat in chain_feats:
            feat_a = np.array([getattr(r, feat) for r in subset], dtype=float)
            mask = ~np.isnan(feat_a) & ~np.isnan(h2_a)
            if mask.sum() < 10:
                cells.append("N/A")
                continue
            rho, p = sp_stats.spearmanr(h2_a[mask], feat_a[mask])
            sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
            cells.append(f"{sig}{rho:+.3f}{sig}")
        lines.append(f"| {d} | " + " | ".join(cells) + " |")

    # === Q2: Correlation with density ===
    lines.append("\n## Q2: Chain Features vs Density (n_causal_pairs)\n")
    lines.append("| d | " + " | ".join(chain_feats) + " |")
    lines.append("|---|" + "|".join(["------"] * len(chain_feats)) + "|")

    for d in dims:
        subset = [r for r in rows if r.d == d]
        if len(subset) < 10:
            continue
        dens_a = np.array([float(r.n_causal_pairs) for r in subset])
        cells = []
        for feat in chain_feats:
            feat_a = np.array([getattr(r, feat) for r in subset], dtype=float)
            mask = ~np.isnan(feat_a)
            if mask.sum() < 10:
                cells.append("N/A")
                continue
            rho, _ = sp_stats.spearmanr(dens_a[mask], feat_a[mask])
            cells.append(f"{rho:+.3f}")
        lines.append(f"| {d} | " + " | ".join(cells) + " |")

    # === Q5: Density-residual analysis ===
    lines.append("\n## Q5: Density-Residual Analysis (OLS remove n_causal_pairs)\n")
    lines.append("After removing density, does the residual still correlate with H²?\n")
    lines.append("| d | N | feature | raw ρ(feat,H²) | ρ(feat,density) | resid ρ | p_resid | verdict |")
    lines.append("|---|---|---------|----------------|-----------------|---------|---------|---------|")

    beyond_count = 0
    total_tested = 0

    for d in dims:
        for N in ns:
            subset = [r for r in rows if r.d == d and r.N == N]
            if len(subset) < 15:
                continue

            h2_a = np.array([r.H_sq for r in subset])
            dens_a = np.array([float(r.n_causal_pairs) for r in subset])

            for feat in chain_feats:
                feat_a = np.array([getattr(r, feat) for r in subset], dtype=float)
                mask = ~np.isnan(feat_a) & ~np.isnan(h2_a)
                if mask.sum() < 15:
                    continue

                total_tested += 1
                rho_raw, _ = sp_stats.spearmanr(h2_a[mask], feat_a[mask])
                rho_dens, _ = sp_stats.spearmanr(dens_a[mask], feat_a[mask])

                # OLS: feat = a * density + b
                coeffs = np.polyfit(dens_a[mask], feat_a[mask], 1)
                resid = feat_a[mask] - np.polyval(coeffs, dens_a[mask])

                rho_resid, p_resid = sp_stats.spearmanr(h2_a[mask], resid)

                if abs(rho_resid) > 0.3 and p_resid < 0.01:
                    verdict = "**BEYOND DENSITY**"
                    beyond_count += 1
                elif abs(rho_resid) < 0.15:
                    verdict = "density-dominated"
                else:
                    verdict = "marginal"

                lines.append(
                    f"| {d} | {N} | {feat} | {rho_raw:+.3f} | {rho_dens:+.3f} "
                    f"| {rho_resid:+.3f} | {p_resid:.2e} | {verdict} |"
                )

    lines.append(f"\n**Total: {beyond_count}/{total_tested} features BEYOND DENSITY**\n")

    # === Conclusion ===
    lines.append("## Conclusion\n")
    if beyond_count > 0:
        pct = beyond_count / total_tested * 100 if total_tested > 0 else 0
        lines.append(f"**{beyond_count}/{total_tested}** ({pct:.0f}%) chain features carry")
        lines.append(" curvature information beyond density.\n")
    else:
        lines.append("**No chain features survive density residualization.**")
        lines.append(" Chain counts are functions of {C_k} interval structure")
        lines.append(" and do not provide an independent curvature channel.\n")

    lines.append("### Comparison with other DDT escape channels\n")
    lines.append("| Channel | Path | Beyond density |")
    lines.append("|---------|------|---------------|")
    lines.append("| B_ℓ spectral | I.1 | 6/18 |")
    lines.append("| Antichain transverse | I.2 | **21/21** |")
    lines.append("| Schwarzschild | I.3 | 2/27 |")
    lines.append(f"| **Chain statistics** | **I.4** | **{beyond_count}/{total_tested}** |")

    # Physical interpretation
    lines.append("\n### Physical Interpretation\n")
    if beyond_count == 0:
        lines.append("k-chain counts N_k are determined by the causal matrix adjacency structure,")
        lines.append(" which for sprinklings is governed by the same density that controls {C_k}.")
        lines.append(" The chain ratios N_{k+1}/N_k measure the average 'branching factor' of")
        lines.append(" causal paths, but this branching is itself a function of density.")
        lines.append(" Unlike antichains (transverse/spacelike) which measure an orthogonal")
        lines.append(" geometric dimension, chains are longitudinal (timelike) like {C_k}.")
        lines.append(" Both live in the same 'longitudinal' sector that DDT constrains.\n")
    else:
        lines.append("Some chain features carry curvature signal beyond density.")
        lines.append(" The surviving features likely encode path-correlation structure")
        lines.append(" that goes beyond pairwise interval counting.\n")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description="§4.1.30: Causal chain statistics experiment"
    )
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256, 512])
    ap.add_argument("--hubbles", nargs="*", type=float,
                    default=[0.0, 0.25, 0.5, 1.0, 1.5, 2.0])
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2030)
    ap.add_argument("--out",
                    default="outputs_unified_functional/conjecture_e_chain_statistics.csv")
    ap.add_argument("--report",
                    default="outputs_unified_functional/conjecture_e_chain_statistics.md")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[ExpRow] = []
    total = len(args.dims) * len(args.ns) * len(args.hubbles) * args.reps
    done = 0

    print(f"Chain statistics experiment: {total} realizations")
    print(f"  dims={args.dims}, ns={args.ns}, H={args.hubbles}, reps={args.reps}")
    print()

    for d in args.dims:
        for N in args.ns:
            for H in args.hubbles:
                for rep in range(args.reps):
                    seed = args.seed + d * 100000 + N * 1000 + int(H * 100) + rep
                    row = run_single(d, N, H, rep, seed)
                    rows.append(row)
                    done += 1
                    if done % 20 == 0 or done == total:
                        print(
                            f"  [{done:4d}/{total}] d={d} N={N:4d} H={H:.2f} "
                            f"N3={row.N_3:.0f} N4={row.N_4:.0f} "
                            f"lc={row.longest_chain:.0f} "
                            f"N3/N2={row.N3_over_N2:.3f}"
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
        subset = [r for r in rows if r.d == d]
        if len(subset) < 10:
            continue
        h2_a = np.array([r.H_sq for r in subset])
        for feat in ["N3_over_N2", "N4_over_N3", "longest_chain_ratio", "chain_entropy"]:
            feat_a = np.array([getattr(r, feat) for r in subset], dtype=float)
            mask = ~np.isnan(feat_a)
            if mask.sum() >= 10:
                rho, p = sp_stats.spearmanr(h2_a[mask], feat_a[mask])
                print(f"  d={d}: ρ({feat}, H²) = {rho:+.3f}  (p={p:.2e})")
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
