"""
Prediction A — Unification Test: Link density ↔ Geometric consistency
======================================================================

Goal: Show that the link action (C0/N) and the BDG geometric consistency
penalty are measuring the same underlying structural signal — effective
dimension.  If both observables correlate strongly with d_eff across
families, then the link-line and consistency-line are *unified*, closing
the second gap toward Prediction A confirmation.

We re-generate the same posets (same seeds as large-N script) for N≤80
and collect FULL geometric_components() output to compare with link_density.
"""

import sys, pathlib, warnings
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)
from observables_geo import geometric_components
from runtime_utils import estimate_entropy_by_family

OUT_DIR = pathlib.Path("outputs_exploratory/prediction_a_unification")
OUT_DIR.mkdir(parents=True, exist_ok=True)

warnings.filterwarnings("ignore")

# ===== Configuration =====
N_VALUES = [20, 36, 52, 68]   # N≤80, but skip 80 for speed (geo is slow)
SAMPLES = 4
SIS_RUNS = 4096
SEED_BASE = 980000

GENERATORS = {
    "2d": generate_lorentzian_like_2d,
    "3d": generate_lorentzian_like_3d,
    "4d": generate_lorentzian_like_4d,
    "5d": generate_lorentzian_like_5d,
}

FAMILY_EXACT_THRESHOLDS = {
    "lorentzian_like_2d": 68,
    "lorentzian_like_3d": 24,
    "lorentzian_like_4d": 24,
    "lorentzian_like_5d": 24,
}


def hasse_links(poset: Poset) -> int:
    c = poset.closure.astype(np.uint8, copy=False)
    has_intermediate = (c @ c).astype(bool, copy=False)
    cover = poset.closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return int(cover.sum())


def collect_data() -> pd.DataFrame:
    rows = []
    for n in N_VALUES:
        print(f"\n  N={n}", end="", flush=True)
        for dim_label, gen in GENERATORS.items():
            fam_name = f"lorentzian_like_{dim_label}"
            exact_thr = FAMILY_EXACT_THRESHOLDS[fam_name]

            for sid in range(SAMPLES):
                seed = SEED_BASE + 1000 * n + sid
                poset = gen(n=n, seed=seed)

                # Entropy
                log_h, method = estimate_entropy_by_family(
                    poset, family=fam_name, sis_runs=SIS_RUNS, seed=seed,
                    default_exact_threshold=24,
                    family_exact_thresholds={fam_name: exact_thr},
                )

                # Link action
                n_links = hasse_links(poset)
                link_density = n_links / (n * (n - 1) / 2)
                C0_per_N = n_links / n

                # Geometric components — the full dict
                try:
                    geo = geometric_components(poset)
                except Exception:
                    geo = {}

                row = {
                    "n": n, "dim": dim_label,
                    "sample_id": sid, "seed": seed,
                    "log_H": float(log_h),
                    "C0_links": n_links,
                    "link_density": link_density,
                    "C0_per_N": C0_per_N,
                    "S_link_d2_norm": (n - 2 * n_links) / n,
                }
                # Add all geo_ keys
                for k, v in geo.items():
                    row[k] = float(v) if isinstance(v, (int, float, np.floating)) else v

                rows.append(row)
            print(f" [{dim_label}]", end="", flush=True)
        print()
    return pd.DataFrame(rows)


def analyze(df: pd.DataFrame) -> str:
    """Analyze correlation between link and geo observables."""
    lines = []
    lines.append("=" * 72)
    lines.append("UNIFICATION TEST: Link Action ↔ Geometric Components")
    lines.append("=" * 72)

    # Identify geo columns
    geo_cols = [c for c in df.columns if c.startswith("geo_")]
    lines.append(f"\n  Geo columns available: {geo_cols}")

    # Aggregate by (n, dim)
    agg_cols = ["link_density", "C0_per_N", "S_link_d2_norm", "log_H"] + geo_cols
    agg_dict = {c: "mean" for c in agg_cols}
    agg = df.groupby(["n", "dim"]).agg(agg_dict).reset_index()
    
    dim_num = {"2d": 2, "3d": 3, "4d": 4, "5d": 5}
    agg["d_num"] = agg["dim"].map(dim_num)

    # ------ Part 1: Per-N correlations ------
    lines.append("\n\n--- Part 1: Per-N Pearson r (link_density vs geo measures) ---")
    target_geos = ["geo_d_order", "geo_d_chain", "geo_dim_eff",
                    "geo_dim_multi_consistency", "geo_dim_proxy_penalty", "geo_total"]
    target_geos = [g for g in target_geos if g in agg.columns]
    
    hdr = f"\n  {'N':>4}"
    for g in target_geos:
        short = g.replace("geo_", "")
        hdr += f"  {short:>20}"
    lines.append(hdr)

    for n in sorted(agg["n"].unique()):
        sub = agg[agg["n"] == n]
        if len(sub) < 3:
            continue
        line = f"  {n:>4}"
        x = sub["link_density"].values
        for g in target_geos:
            y = sub[g].values
            if np.std(x) < 1e-12 or np.std(y) < 1e-12:
                line += f"  {'--':>20}"
            else:
                r, p = pearsonr(x, y)
                line += f"  {r:>7.3f} (p={p:.3f})  "
        lines.append(line)

    # ------ Part 2: Overall correlation (pooled across N) ------
    lines.append("\n\n--- Part 2: Pooled Pearson r across ALL (N, dim) rows ---")
    x_all = agg["link_density"].values
    for g in target_geos:
        y_all = agg[g].values
        mask = ~(np.isnan(x_all) | np.isnan(y_all))
        if mask.sum() < 3 or np.std(x_all[mask]) < 1e-12 or np.std(y_all[mask]) < 1e-12:
            lines.append(f"  link_density vs {g}: not computable")
        else:
            r, p = pearsonr(x_all[mask], y_all[mask])
            rs, ps = spearmanr(x_all[mask], y_all[mask])
            lines.append(f"  link_density vs {g:>30}: r={r:.4f} (p={p:.4f}), rho={rs:.4f} (p={ps:.4f})")

    # ------ Part 3: Profile table ------
    lines.append("\n\n--- Part 3: Observable profiles across dimensions ---")
    for n in sorted(agg["n"].unique()):
        sub = agg[agg["n"] == n].sort_values("d_num")
        lines.append(f"\n  N={n}:")
        header = f"    {'dim':>4}  {'link_dens':>10}  {'C0/N':>8}  {'d_order':>8}"
        if "geo_dim_eff" in sub.columns:
            header += f"  {'d_eff':>8}"
        if "geo_dim_multi_consistency" in sub.columns:
            header += f"  {'multi_cons':>12}"
        if "geo_dim_proxy_penalty" in sub.columns:
            header += f"  {'proxy_pen':>12}"
        if "geo_total" in sub.columns:
            header += f"  {'geo_total':>10}"
        lines.append(header)

        for _, row in sub.iterrows():
            line = f"    {row['dim']:>4}  {row['link_density']:>10.4f}  {row['C0_per_N']:>8.4f}"
            if "geo_d_order" in row.index:
                line += f"  {row['geo_d_order']:>8.3f}"
            if "geo_dim_eff" in row.index:
                line += f"  {row['geo_dim_eff']:>8.3f}"
            if "geo_dim_multi_consistency" in row.index:
                line += f"  {row['geo_dim_multi_consistency']:>12.4f}"
            if "geo_dim_proxy_penalty" in row.index:
                line += f"  {row['geo_dim_proxy_penalty']:>12.4f}"
            if "geo_total" in row.index:
                line += f"  {row['geo_total']:>10.3f}"
            lines.append(line)

    # ------ Part 4: Link density gap structure ------
    lines.append("\n\n--- Part 4: Δ(link_density) between adjacent dimensions ---")
    lines.append(f"  {'N':>4}  {'Δ(2→3)':>10}  {'Δ(3→4)':>10}  {'Δ(4→5)':>10}  {'gap_ratio':>12}")
    for n in sorted(agg["n"].unique()):
        sub = agg[agg["n"] == n].set_index("dim")
        if all(d in sub.index for d in ["2d", "3d", "4d", "5d"]):
            d23 = sub.loc["3d", "link_density"] - sub.loc["2d", "link_density"]
            d34 = sub.loc["4d", "link_density"] - sub.loc["3d", "link_density"]
            d45 = sub.loc["5d", "link_density"] - sub.loc["4d", "link_density"]
            ratio = abs(d45) / abs(d34) if abs(d34) > 1e-12 else float("inf")
            lines.append(f"  {n:>4}  {d23:>10.4f}  {d34:>10.4f}  {d45:>10.4f}  {ratio:>12.2f}")

    # ------ Part 5: Key unification metric ------
    lines.append("\n\n--- Part 5: Unification Verdict ---")
    # Test: does d_order track link_density across dims?
    if "geo_d_order" in agg.columns:
        r_link_dorder, p_link_dorder = pearsonr(agg["link_density"].values, agg["geo_d_order"].values)
        lines.append(f"\n  r(link_density, d_order) = {r_link_dorder:.4f} (p={p_link_dorder:.6f})")
        if abs(r_link_dorder) > 0.7:
            lines.append("  → Link density and d_order are STRONGLY correlated ✓")
            lines.append("  → The link action is effectively measuring causal dimension")
            lines.append("  → BDG consistency term and link action select by the SAME signal")
            lines.append("\n  ★ UNIFICATION: Link-line and consistency-line converge ★")
        else:
            lines.append("  → Correlation is moderate or weak — partial unification only")

    # Additional: does the link action separate 4D from 5D more than 3D from 4D?
    # This is quantified by the gap ratio
    gap_ratios = []
    for n in sorted(agg["n"].unique()):
        sub = agg[agg["n"] == n].set_index("dim")
        if all(d in sub.index for d in ["3d", "4d", "5d"]):
            d34 = abs(sub.loc["4d", "link_density"] - sub.loc["3d", "link_density"])
            d45 = abs(sub.loc["5d", "link_density"] - sub.loc["4d", "link_density"])
            if d34 > 1e-12 and d45 > 1e-12:
                gap_ratios.append(d45 / d34)

    if gap_ratios:
        mean_ratio = np.mean(gap_ratios)
        lines.append(f"\n  Mean |Δ(4→5)| / |Δ(3→4)| = {mean_ratio:.2f}")
        if mean_ratio > 1.5:
            lines.append("  → The 4→5 link density gap is LARGER than 3→4")
            lines.append("  → This asymmetry explains why Ξ₄→₅ > Ξ₃→₄: structural sparsity jumps at d=5")
        else:
            lines.append("  → Gap ratio near 1 — link density gaps are roughly symmetric")

    return "\n".join(lines)


def main():
    print("=" * 72)
    print("Prediction A — Unification: Link Action ↔ Geometric Components")
    print("=" * 72)

    print("\n[1/2] Collecting data (N ≤ 68, 4 dims, 4 samples)...")
    df = collect_data()
    df.to_csv(OUT_DIR / "unification_data.csv", index=False)
    print(f"\n  Saved {len(df)} observations")

    print("\n[2/2] Analyzing correlations...")
    report = analyze(df)
    print(report)

    with open(OUT_DIR / "unification_report.txt", "w", encoding="utf-8") as f:
        f.write(report)

    print(f"\n  Report saved to {OUT_DIR}/unification_report.txt")


if __name__ == "__main__":
    main()
