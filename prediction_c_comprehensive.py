"""Prediction C comprehensive validation: deep hierarchy integration → lower log_H.

This script performs a large-scale, all-family analysis testing the core
hypothesis that deeper causal hierarchy integration systematically lowers
combinatorial entropy (log_H) and improves coarse-graining stability.

Three analysis tiers:
  Tier 1 — All-family partial correlation (frozen_exact, N=10..16, 8 families)
  Tier 2 — Pairwise delta analysis on expanded matched pairs (N=30..48)
  Tier 3 — Layer-count dominance decomposition & CG-stability linkage
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from experiment import FAMILIES
from matched_residual_freedom_check import residual_metrics
from observables import antichain_width, comparable_fraction, layer_profile
from observables_geo import cover_density, dimension_proxy_views


# ── helpers ──────────────────────────────────────────────────────────────

def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def zscore(arr: np.ndarray) -> np.ndarray:
    m = arr.mean()
    s = arr.std(ddof=0)
    s = s if s > 1e-12 else 1.0
    return (arr - m) / s


def regression_residual(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    xc = x - x.mean()
    yc = y - y.mean()
    d = math.sqrt(float((xc * xc).sum() * (yc * yc).sum()))
    return float((xc * yc).sum() / d) if d > 1e-12 else 0.0


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return pearson(rx, ry)


def permutation_pvalue(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    obs = abs(pearson(x, y))
    count = 1  # continuity correction
    for _ in range(n_perm):
        y_perm = rng.permutation(y)
        if abs(pearson(x, y_perm)) >= obs:
            count += 1
    return float(count / (n_perm + 1))


def partial_corr_and_p(y: np.ndarray, x: np.ndarray, controls: np.ndarray,
                       n_perm: int, seed: int) -> dict:
    """Partial correlation of x with y, controlling for columns in controls."""
    X = np.column_stack([np.ones(len(y)), controls])
    y_res = regression_residual(y, X)
    x_res = regression_residual(x, X)
    r = pearson(x_res, y_res)
    rho = spearman(x_res, y_res)
    p = permutation_pvalue(x_res, y_res, n_perm=n_perm, seed=seed)
    return {"partial_pearson": r, "partial_spearman": rho, "permutation_p": p, "n": len(y)}


# ── Tier 1: All-family frozen_exact analysis ─────────────────────────────

def tier1_all_family(raw_csv_paths: list[str], families: list[str],
                     n_perm: int, seed: int, out_dir: Path) -> pd.DataFrame:
    """Compute HII for every sample in frozen_exact and test partial correlation."""
    print("\n" + "=" * 70)
    print("TIER 1: All-family hierarchy integration analysis")
    print("=" * 70)

    parts = []
    for p in raw_csv_paths:
        df = pd.read_csv(p)
        # frozen_exact has multiple gammas/modes per sample; take unique samples
        df_unique = df.drop_duplicates(subset=["family", "n", "sample_id"]).copy()
        parts.append(df_unique)
    raw = pd.concat(parts, ignore_index=True)
    raw = raw[raw["family"].isin(families)].copy()
    raw = raw.drop_duplicates(subset=["family", "n", "sample_id"]).copy()

    print(f"  Unique samples after dedup: {len(raw)}")
    print(f"  Families: {sorted(raw['family'].unique())}")
    print(f"  N values: {sorted(raw['n'].unique())}")

    metric_rows = []
    total = len(raw)
    for i, row in enumerate(raw.itertuples(index=False)):
        if (i + 1) % 100 == 0:
            print(f"  Computing metrics... {i+1}/{total}")
        fam = str(row.family)
        n = int(row.n)
        sid = int(row.sample_id)
        seed_gen = 1000 * n + sid
        poset = FAMILIES[fam](n=n, seed=seed_gen)
        metrics = residual_metrics(poset)
        aw = antichain_width(poset)
        cf = comparable_fraction(poset)
        dv = dimension_proxy_views(poset)
        metric_rows.append({
            "family": fam,
            "n": n,
            "sample_id": sid,
            "log_H": float(row.log_H_mean),
            "antichain_width": float(aw),
            "comparable_fraction": cf,
            "geo_dim_eff": dv["d_order"],
            **metrics,
        })

    df = pd.DataFrame(metric_rows)

    # Build HII
    for col in ["layer_count", "mean_layer_gap", "long_edge_fraction",
                "adjacent_edge_fraction", "reduction_edge_density"]:
        df[f"z_{col}"] = zscore(df[col].to_numpy())

    df["hierarchy_integration_index"] = (
        df["z_layer_count"]
        + df["z_mean_layer_gap"]
        + df["z_long_edge_fraction"]
        - df["z_adjacent_edge_fraction"]
        - df["z_reduction_edge_density"]
    ) / 5.0

    # Controls
    controls = df[["antichain_width", "comparable_fraction", "geo_dim_eff"]].to_numpy()
    y_logh = df["log_H"].to_numpy()
    x_hii = df["hierarchy_integration_index"].to_numpy()

    # Overall partial correlation
    overall = partial_corr_and_p(y_logh, x_hii, controls, n_perm, seed)

    # Per-family breakdown
    family_rows = []
    for fam, sub in df.groupby("family"):
        if len(sub) < 10:
            continue
        c = sub[["antichain_width", "comparable_fraction", "geo_dim_eff"]].to_numpy()
        y = sub["log_H"].to_numpy()
        x = sub["hierarchy_integration_index"].to_numpy()
        r = partial_corr_and_p(y, x, c, n_perm, seed)
        family_rows.append({"family": fam, **r,
                            "mean_hii": float(x.mean()),
                            "mean_log_H": float(y.mean()),
                            "mean_layer_count": float(sub["layer_count"].mean())})

    # Component-level partial correlations
    component_rows = []
    for comp in ["layer_count", "mean_layer_gap", "long_edge_fraction",
                 "adjacent_edge_fraction", "reduction_edge_density",
                 "cover_density", "layer_signature_redundancy"]:
        x_comp = df[comp].to_numpy()
        r = partial_corr_and_p(y_logh, x_comp, controls, n_perm, seed)
        component_rows.append({"component": comp, **r})

    # Save
    df.to_csv(out_dir / "tier1_all_family_raw.csv", index=False, encoding="utf-8-sig")

    summary = pd.DataFrame([{"scope": "all_families", **overall}])
    summary.to_csv(out_dir / "tier1_overall_summary.csv", index=False, encoding="utf-8-sig")

    fam_df = pd.DataFrame(family_rows)
    fam_df.to_csv(out_dir / "tier1_by_family.csv", index=False, encoding="utf-8-sig")

    comp_df = pd.DataFrame(component_rows)
    comp_df.to_csv(out_dir / "tier1_components.csv", index=False, encoding="utf-8-sig")

    print(f"\n  OVERALL: partial_r(HII, log_H | controls) = {overall['partial_pearson']:.4f}, "
          f"p = {overall['permutation_p']:.4f}")
    print(f"  Per-family results:")
    for _, r in fam_df.iterrows():
        print(f"    {r['family']:30s}  r={r['partial_pearson']:+.4f}  p={r['permutation_p']:.4f}  "
              f"mean_HII={r['mean_hii']:+.3f}  mean_lc={r['mean_layer_count']:.1f}")
    print(f"  Component-level partial correlations:")
    for _, r in comp_df.iterrows():
        print(f"    {r['component']:30s}  r={r['partial_pearson']:+.4f}  p={r['permutation_p']:.4f}")

    return df


# ── Tier 2: Expanded pairwise analysis ───────────────────────────────────

def tier2_pairwise(pairs_csv: str, residual_csv: str,
                   n_perm: int, seed: int, out_dir: Path) -> None:
    """Pairwise delta analysis on expanded matched pairs."""
    print("\n" + "=" * 70)
    print("TIER 2: Expanded pairwise delta analysis")
    print("=" * 70)

    pairs = pd.read_csv(pairs_csv)
    residuals = pd.read_csv(residual_csv)
    print(f"  Pairs: {len(pairs)}, Residual rows: {len(residuals)}")

    METRICS = ["layer_count", "mean_layer_gap", "adjacent_edge_fraction",
               "long_edge_fraction", "reduction_edge_density",
               "cover_density", "layer_signature_redundancy"]

    mlr = residuals[residuals["family"] == "multi_layer_random"].copy()
    lor = residuals[residuals["family"] == "lorentzian_like_2d"].copy()
    mlr = mlr.rename(columns={"seed": "mlr_seed"})
    lor = lor.rename(columns={"seed": "lor_seed"})

    merged = pairs[["n", "mlr_seed", "lor_seed",
                     "log_H_delta_mlr_minus_lor2d",
                     "score_A2_gamma_delta_mlr_minus_lor2d"]].merge(
        mlr[["n", "mlr_seed"] + METRICS], on=["n", "mlr_seed"], how="left"
    ).merge(
        lor[["n", "lor_seed"] + METRICS], on=["n", "lor_seed"], how="left",
        suffixes=("_mlr", "_lor2d"),
    )

    for m in METRICS:
        merged[f"{m}_delta"] = merged[f"{m}_mlr"] - merged[f"{m}_lor2d"]

    merged["hii_delta"] = (
        merged["layer_count_delta"]
        + merged["mean_layer_gap_delta"]
        + merged["long_edge_fraction_delta"]
        + merged["reduction_edge_density_delta"]
        - merged["adjacent_edge_fraction_delta"]
    )

    results = []
    targets = ["log_H_delta_mlr_minus_lor2d", "score_A2_gamma_delta_mlr_minus_lor2d"]
    features = ["hii_delta"] + [f"{m}_delta" for m in METRICS]

    for tgt in targets:
        for feat in features:
            x = merged[feat].to_numpy(dtype=float)
            y = merged[tgt].to_numpy(dtype=float)
            mask = np.isfinite(x) & np.isfinite(y)
            x, y = x[mask], y[mask]
            if len(x) < 5:
                continue
            r_p = pearson(x, y)
            r_s = spearman(x, y)
            p = permutation_pvalue(x, y, n_perm, seed)
            results.append({
                "feature": feat, "target": tgt,
                "n_pairs": len(x),
                "pearson": r_p, "spearman": r_s, "permutation_p": p,
            })

    res_df = pd.DataFrame(results)
    merged.to_csv(out_dir / "tier2_pairwise_raw.csv", index=False, encoding="utf-8-sig")
    res_df.to_csv(out_dir / "tier2_pairwise_summary.csv", index=False, encoding="utf-8-sig")

    print(f"\n  Results (target=log_H_delta):")
    sub = res_df[res_df["target"] == "log_H_delta_mlr_minus_lor2d"]
    for _, r in sub.iterrows():
        print(f"    {r['feature']:40s}  r={r['pearson']:+.4f}  rho={r['spearman']:+.4f}  p={r['permutation_p']:.4f}")


# ── Tier 3: CG-stability linkage ─────────────────────────────────────────

def tier3_cg_linkage(duel_csvs: list[str], n_perm: int, seed: int, out_dir: Path) -> None:
    """Test whether HII connects to coarse-graining stability (family_switch_rate)."""
    print("\n" + "=" * 70)
    print("TIER 3: Hierarchy integration → CG stability linkage")
    print("=" * 70)

    parts = [pd.read_csv(p) for p in duel_csvs]
    raw = pd.concat(parts, ignore_index=True)
    print(f"  Total duel rows: {len(raw)}")
    print(f"  Families: {sorted(raw['family'].unique())}")

    if "cg_family_switch_rate" not in raw.columns:
        print("  WARNING: cg_family_switch_rate not available, skipping tier 3")
        return

    # Compute HII for each sample
    rows = []
    for _, r in raw.iterrows():
        fam = str(r["family"])
        n = int(r["n"])
        s = int(r["seed"])
        poset = FAMILIES[fam](n=n, seed=s)
        metrics = residual_metrics(poset)
        rows.append({
            "family": fam, "n": n, "seed": s,
            "log_H": float(r["log_H"]),
            "cg_switch_rate": float(r["cg_family_switch_rate"]),
            "cg_mean_penalty": float(r.get("cg_mean_penalty", float("nan"))),
            **metrics,
        })

    df = pd.DataFrame(rows)
    for col in ["layer_count", "mean_layer_gap", "long_edge_fraction",
                "adjacent_edge_fraction", "reduction_edge_density"]:
        df[f"z_{col}"] = zscore(df[col].to_numpy())

    df["hii"] = (
        df["z_layer_count"] + df["z_mean_layer_gap"] + df["z_long_edge_fraction"]
        - df["z_adjacent_edge_fraction"] - df["z_reduction_edge_density"]
    ) / 5.0

    results = []
    for target_col in ["cg_switch_rate", "log_H"]:
        for feat_col in ["hii", "layer_count", "mean_layer_gap", "long_edge_fraction"]:
            x = df[feat_col].to_numpy(dtype=float)
            y = df[target_col].to_numpy(dtype=float)
            mask = np.isfinite(x) & np.isfinite(y)
            x, y = x[mask], y[mask]
            if len(x) < 5:
                continue
            r_p = pearson(x, y)
            r_s = spearman(x, y)
            p = permutation_pvalue(x, y, n_perm, seed)
            results.append({
                "feature": feat_col, "target": target_col,
                "n_samples": len(x), "pearson": r_p, "spearman": r_s, "permutation_p": p,
            })

    res_df = pd.DataFrame(results)
    df.to_csv(out_dir / "tier3_cg_linkage_raw.csv", index=False, encoding="utf-8-sig")
    res_df.to_csv(out_dir / "tier3_cg_linkage_summary.csv", index=False, encoding="utf-8-sig")

    print(f"\n  Results:")
    for _, r in res_df.iterrows():
        print(f"    {r['feature']:25s} → {r['target']:20s}  r={r['pearson']:+.4f}  p={r['permutation_p']:.4f}")


# ── main ──────────────────────────────────────────────────────────────────

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Prediction C comprehensive validation")
    p.add_argument("--config", default="config_prediction_c_comprehensive.yaml")
    return p


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)

    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)
    n_perm = int(config["experiment"]["n_permutations"])
    seed = int(config["experiment"]["seed"])

    # Tier 1
    tier1_df = tier1_all_family(
        raw_csv_paths=config["tier1"]["raw_csvs"],
        families=config["tier1"]["families"],
        n_perm=n_perm, seed=seed, out_dir=out_dir,
    )

    # Tier 2
    if "tier2" in config:
        tier2_pairwise(
            pairs_csv=config["tier2"]["pairs_csv"],
            residual_csv=config["tier2"]["residual_csv"],
            n_perm=n_perm, seed=seed, out_dir=out_dir,
        )

    # Tier 3
    if "tier3" in config:
        tier3_cg_linkage(
            duel_csvs=config["tier3"]["duel_csvs"],
            n_perm=n_perm, seed=seed, out_dir=out_dir,
        )

    print("\n" + "=" * 70)
    print(f"All outputs saved to: {out_dir}")
    print("=" * 70)
