"""
Prediction D: I_cg dose-response analysis
==========================================
Analogous to Prediction C's k-family dose-response (§4.9):
  C showed that slope(log_H vs k_layers) steepens with N.
  D asks: does slope(improve_rank vs I_cg) steepen with N?

Three complementary analyses:
  1. Within-block OLS slope of improve_rank on I_cg, aggregated per (N, rep)
  2. I_cg tertile binning: mean improve_rank by I_cg bin, checking monotonicity
  3. N-scaling: does the I_cg → improve_rank slope strengthen with N?

Uses confirmatory v9 data from rep3/4/5/6 (4 independent replications).
Focuses on gamma=0.2 (frozen window), zeta>0 (where CG actually shifts ranks),
and icg_variant ∈ {full, switch, no_switch}.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(__file__).parent
OUT_DIR = BASE / "outputs_exploratory" / "prediction_d_dose_response"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BLOCK_COLS = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode", "zeta", "icg_variant"]


def load_all_reps() -> pd.DataFrame:
    """Load v9 confirmatory data from rep3-6, add improve_rank."""
    frames = []
    for repn in [3, 4, 5, 6]:
        path = BASE / f"outputs_confirmatory/prediction_d_dynamic_v9_confirm_rep{repn}/cg_zeta_scan_rankings_variants.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path)
        df["rep"] = repn
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df["improve_rank"] = -(df["rank_eval"].astype(float) - df["rank_local_eval"].astype(float))
    return df


def _ols_slope(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Return (slope, r, p) from OLS regression of y on x."""
    if len(x) < 3 or np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return (np.nan, np.nan, np.nan)
    slope, intercept, r, p, se = stats.linregress(x, y)
    return (slope, r, p)


# ═══════════════════════════════════════════════════════════════
# Analysis 1: Within-block OLS slope, aggregated by (N, variant)
# ═══════════════════════════════════════════════════════════════
def analysis_block_slopes(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each block, compute OLS slope of improve_rank on mean_I_cg.
    Then aggregate per (N, icg_variant, rep) and per (N, icg_variant).
    """
    print("=" * 70)
    print("Analysis 1: Within-Block OLS Slope of improve_rank on I_cg")
    print("=" * 70)

    rows = []
    for keys, block in df.groupby(BLOCK_COLS):
        if len(block) < 4:
            continue
        x = block["mean_I_cg"].to_numpy(float)
        y = block["improve_rank"].to_numpy(float)
        slope, r, p = _ols_slope(x, y)
        key_dict = dict(zip(BLOCK_COLS, keys))
        rows.append({
            **key_dict,
            "n_families": len(block),
            "slope": slope,
            "r": r,
            "p": p,
            "rep": int(block["rep"].iloc[0]),
        })

    block_df = pd.DataFrame(rows)
    block_df.to_csv(OUT_DIR / "block_slopes.csv", index=False)

    # Aggregate by (N, icg_variant, rep)
    for variant in ["full", "switch", "no_switch"]:
        vsub = block_df[block_df["icg_variant"] == variant]
        if vsub.empty:
            continue
        print(f"\n  [{variant}] Per-rep mean slope:")
        agg = vsub.groupby(["n", "rep"]).agg(
            mean_slope=("slope", "mean"),
            mean_r=("r", "mean"),
            n_blocks=("slope", "count"),
        ).reset_index()
        for _, row in agg.sort_values(["n", "rep"]).iterrows():
            print(f"    N={int(row['n']):3d}  rep={int(row['rep'])}  "
                  f"mean_slope={row['mean_slope']:+.4f}  mean_r={row['mean_r']:+.4f}  "
                  f"n_blocks={int(row['n_blocks'])}")

    # Grand aggregate by (N, icg_variant) across all reps
    print("\n  Grand aggregate (all reps):")
    grand = block_df.groupby(["n", "icg_variant"]).agg(
        mean_slope=("slope", "mean"),
        sd_slope=("slope", "std"),
        mean_r=("r", "mean"),
        n_blocks=("slope", "count"),
    ).reset_index()
    grand.to_csv(OUT_DIR / "slope_by_n_variant.csv", index=False)

    for _, row in grand.sort_values(["icg_variant", "n"]).iterrows():
        print(f"    {row['icg_variant']:12s}  N={int(row['n']):3d}  "
              f"slope={row['mean_slope']:+.4f} ± {row['sd_slope']:.4f}  "
              f"r={row['mean_r']:+.4f}  (n_blocks={int(row['n_blocks'])})")

    return grand


# ═══════════════════════════════════════════════════════════════
# Analysis 2: I_cg tertile binning — monotonicity check
# ═══════════════════════════════════════════════════════════════
def analysis_tertile_binning(df: pd.DataFrame) -> pd.DataFrame:
    """
    Within each (N, icg_variant), bin families into I_cg tertiles
    (low / mid / high) and compute mean improve_rank per bin.
    Check for monotonic increase.
    """
    print("\n" + "=" * 70)
    print("Analysis 2: I_cg Tertile Binning → improve_rank Monotonicity")
    print("=" * 70)

    rows = []
    for (n_val, variant), group in df.groupby(["n", "icg_variant"]):
        x = group["mean_I_cg"].values
        y = group["improve_rank"].values

        # Compute tertile boundaries
        t1 = np.percentile(x, 33.3)
        t2 = np.percentile(x, 66.7)

        bins = np.where(x <= t1, "low", np.where(x <= t2, "mid", "high"))
        for bname in ["low", "mid", "high"]:
            mask = bins == bname
            if mask.sum() == 0:
                continue
            rows.append({
                "n": int(n_val),
                "icg_variant": str(variant),
                "I_cg_bin": bname,
                "n_obs": int(mask.sum()),
                "mean_I_cg": float(np.mean(x[mask])),
                "mean_improve_rank": float(np.mean(y[mask])),
                "sd_improve_rank": float(np.std(y[mask], ddof=1)) if mask.sum() > 1 else 0.0,
            })

    bin_df = pd.DataFrame(rows)
    bin_df.to_csv(OUT_DIR / "tertile_binning.csv", index=False)

    # Check monotonicity
    for variant in ["full", "switch", "no_switch"]:
        vsub = bin_df[bin_df["icg_variant"] == variant]
        if vsub.empty:
            continue
        print(f"\n  [{variant}]")
        for n_val in sorted(vsub["n"].unique()):
            nsub = vsub[vsub["n"] == n_val].set_index("I_cg_bin")
            vals = [nsub.loc[b, "mean_improve_rank"] if b in nsub.index else np.nan
                    for b in ["low", "mid", "high"]]
            icg_means = [nsub.loc[b, "mean_I_cg"] if b in nsub.index else np.nan
                         for b in ["low", "mid", "high"]]
            mono = "✓ monotonic" if vals[0] <= vals[1] <= vals[2] else "✗ NOT monotonic"
            # Jonckheere-Terpstra style trend test (simple version: correlation of bin order with Y)
            print(f"    N={int(n_val):3d}  "
                  f"low({icg_means[0]:+.3f})→{vals[0]:+.3f}  "
                  f"mid({icg_means[1]:+.3f})→{vals[1]:+.3f}  "
                  f"high({icg_means[2]:+.3f})→{vals[2]:+.3f}  {mono}")

    return bin_df


# ═══════════════════════════════════════════════════════════════
# Analysis 3: N-scaling of the dose-response slope
# ═══════════════════════════════════════════════════════════════
def analysis_n_scaling(grand: pd.DataFrame):
    """
    Check whether the I_cg → improve_rank slope strengthens with N.
    If it does, this parallels C's k-family dose-response steepening.
    """
    print("\n" + "=" * 70)
    print("Analysis 3: N-Scaling of the Dose-Response Slope")
    print("=" * 70)

    rows = []
    for variant in ["full", "switch", "no_switch"]:
        vsub = grand[grand["icg_variant"] == variant].sort_values("n")
        if len(vsub) < 2:
            continue
        ns = vsub["n"].values
        slopes = vsub["mean_slope"].values
        abs_slopes = np.abs(slopes)

        # Linear trend in |slope| vs N
        if len(ns) >= 3:
            trend_slope, trend_r, trend_p = _ols_slope(ns.astype(float), abs_slopes)
        else:
            trend_slope, trend_r, trend_p = np.nan, np.nan, np.nan

        # Is |slope| monotonically increasing with N?
        diffs = np.diff(abs_slopes)
        mono = bool(np.all(diffs >= 0))

        rows.append({
            "icg_variant": variant,
            "n_points": len(ns),
            "slopes": str(dict(zip(ns.tolist(), [f"{s:+.4f}" for s in slopes]))),
            "abs_slopes": str(dict(zip(ns.tolist(), [f"{s:.4f}" for s in abs_slopes]))),
            "monotonic_increasing": mono,
            "trend_slope_per_N": trend_slope,
            "trend_r": trend_r,
            "trend_p": trend_p,
        })

        direction = "↑ steepening" if mono else "— NOT steepening"
        print(f"\n  [{variant}]")
        for n_val, s in zip(ns, slopes):
            print(f"    N={int(n_val):3d}  slope = {s:+.4f}")
        print(f"    Monotonically steepening: {mono}  {direction}")
        if np.isfinite(trend_p):
            print(f"    Linear trend: d(|slope|)/dN = {trend_slope:.5f}, r={trend_r:.4f}, p={trend_p:.4f}")

    trend_df = pd.DataFrame(rows)
    trend_df.to_csv(OUT_DIR / "n_scaling.csv", index=False)
    return trend_df


# ═══════════════════════════════════════════════════════════════
# Analysis 4: Pooled dose-response curve (all reps, per N)
# ═══════════════════════════════════════════════════════════════
def analysis_pooled_dose_response(df: pd.DataFrame):
    """
    Pool all blocks at each N and fit I_cg → improve_rank OLS.
    Report slope, r, p, n for each N and variant.
    This is the direct analogue of C's Table in §4.9.
    """
    print("\n" + "=" * 70)
    print("Analysis 4: Pooled Dose-Response Table (Analogue of C §4.9)")
    print("=" * 70)

    rows = []
    for variant in ["full", "switch", "no_switch"]:
        for n_val in sorted(df["n"].unique()):
            sub = df[(df["n"] == n_val) & (df["icg_variant"] == variant)]
            if len(sub) < 5:
                continue
            x = sub["mean_I_cg"].to_numpy(float)
            y = sub["improve_rank"].to_numpy(float)
            slope, r, p = _ols_slope(x, y)

            # Also compute Spearman
            rho_s, p_s = stats.spearmanr(x, y)

            rows.append({
                "icg_variant": variant,
                "n": int(n_val),
                "n_obs": len(sub),
                "ols_slope": slope,
                "pearson_r": r,
                "pearson_p": p,
                "spearman_rho": rho_s,
                "spearman_p": p_s,
            })
            print(f"  {variant:12s}  N={int(n_val):3d}  n={len(sub):4d}  "
                  f"slope={slope:+.4f}  r={r:+.4f}  p={p:.2e}  "
                  f"rho_s={rho_s:+.4f}  p_s={p_s:.2e}")

    pooled_df = pd.DataFrame(rows)
    pooled_df.to_csv(OUT_DIR / "pooled_dose_response.csv", index=False)
    return pooled_df


# ═══════════════════════════════════════════════════════════════
# Analysis 5: Component decomposition of I_cg
# ═══════════════════════════════════════════════════════════════
def analysis_component_decomposition(df: pd.DataFrame):
    """
    Decompose: which components of the I_cg signal drive the
    improve_rank association?
    - mean_retain_identity (identity preservation)
    - mean_penalty_cg (CG penalty stability)
    - mean_sig_* (signature stability metrics)
    """
    print("\n" + "=" * 70)
    print("Analysis 5: I_cg Component Decomposition")
    print("=" * 70)

    components = [
        "mean_retain_identity",
        "mean_penalty_cg",
        "mean_sig_comp",
        "mean_sig_d_eff",
        "mean_sig_height_ratio",
        "mean_sig_width_ratio",
        "mean_sig_degree_var",
    ]

    available = [c for c in components if c in df.columns]

    rows = []
    for variant in ["full"]:  # focus on full for decomposition
        for n_val in sorted(df["n"].unique()):
            sub = df[(df["n"] == n_val) & (df["icg_variant"] == variant)]
            y = sub["improve_rank"].to_numpy(float)
            print(f"\n  [{variant}, N={int(n_val)}]")
            for comp in available:
                x = sub[comp].to_numpy(float)
                if np.std(x) < 1e-12:
                    continue
                rho, p = stats.spearmanr(x, y)
                rows.append({
                    "icg_variant": variant,
                    "n": int(n_val),
                    "component": comp,
                    "spearman_rho": rho,
                    "spearman_p": p,
                    "n_obs": len(sub),
                })
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                print(f"    {comp:30s}  rho={rho:+.4f}  p={p:.3e}  {sig}")

    comp_df = pd.DataFrame(rows)
    comp_df.to_csv(OUT_DIR / "component_decomposition.csv", index=False)
    return comp_df


# ═══════════════════════════════════════════════════════════════
# SYNTHESIS
# ═══════════════════════════════════════════════════════════════
def synthesis(grand, pooled_df):
    print("\n" + "=" * 70)
    print("SYNTHESIS: Prediction D Dose-Response Assessment")
    print("=" * 70)

    # Check if pooled slopes are all positive
    full = pooled_df[pooled_df["icg_variant"] == "full"]
    all_positive = bool((full["ols_slope"] > 0).all())
    all_sig = bool((full["pearson_p"] < 0.001).all())

    # Check monotonicity across N
    full_grand = grand[grand["icg_variant"] == "full"].sort_values("n")
    slopes = full_grand["mean_slope"].values
    abs_slopes = np.abs(slopes)
    steepening = bool(np.all(np.diff(abs_slopes) >= 0)) if len(abs_slopes) >= 2 else False

    print(f"""
  Dose-response summary (icg_variant = full):
    All slopes positive: {all_positive}
    All p < 0.001: {all_sig}
    |slope| steepening with N: {steepening}

  Interpretation:
    {'✓' if all_positive and all_sig else '✗'} Direction: higher I_cg → higher improve_rank (consistently positive)
    {'✓' if steepening else '—'} N-scaling: {'slope strengthens with N (parallels C dose-response)' if steepening else 'slope does NOT consistently strengthen with N'}

  Comparison with Prediction C dose-response:
    C: k-family slope = -0.30 (N=14), -0.53 (N=16), -0.70 (N=18) — steepens
    D: see table above

  Assessment:
    {'DOSE-RESPONSE SUPPORTED' if all_positive and all_sig else 'DOSE-RESPONSE PARTIALLY SUPPORTED'}
    {'with N-steepening' if steepening else 'WITHOUT N-steepening (dose-response exists but does not strengthen)'}
""")


if __name__ == "__main__":
    print()
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║   Prediction D: I_cg Dose-Response Analysis                       ║")
    print("║   Does improve_rank increase monotonically with I_cg dose?        ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    print()

    df = load_all_reps()
    # Focus on gamma=0.2 (frozen window), zeta>0 (where CG matters)
    df = df[(df["gamma"] == 0.2) & (df["zeta"] > 0)].copy()
    df = df[df["icg_variant"].isin(["full", "switch", "no_switch"])].copy()
    print(f"Data: {len(df)} rows, N ∈ {sorted(df['n'].unique())}, "
          f"reps = {sorted(df['rep'].unique())}, variants = {sorted(df['icg_variant'].unique())}")
    print()

    grand = analysis_block_slopes(df)
    analysis_tertile_binning(df)
    trend = analysis_n_scaling(grand)
    pooled_df = analysis_pooled_dose_response(df)
    analysis_component_decomposition(df)
    synthesis(grand, pooled_df)
