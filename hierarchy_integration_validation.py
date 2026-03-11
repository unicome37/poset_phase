from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from experiment import FAMILIES
from matched_residual_freedom_check import residual_metrics


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def zscore(series: pd.Series) -> pd.Series:
    mean = float(series.mean())
    std = float(series.std(ddof=0))
    std = std if std > 1e-12 else 1.0
    return (series - mean) / std


def regression_residual(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = x - x.mean()
    y = y - y.mean()
    denom = math.sqrt(float((x * x).sum()) * float((y * y).sum()))
    if denom <= 1e-12:
        return 0.0
    return float((x * y).sum() / denom)


def permutation_pvalue(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    obs = abs(pearson_corr(x, y))
    count = 1
    for _ in range(n_perm):
        y_perm = rng.permutation(y)
        if abs(pearson_corr(x, y_perm)) >= obs:
            count += 1
    return float(count / (n_perm + 1))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate whether hierarchical integration predicts lower log_H after controlling width/comp/dim_eff.")
    parser.add_argument("--config", default="config_hierarchy_integration_validation.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_parts = [pd.read_csv(path) for path in config["input"]["raw_csvs"]]
    raw_df = pd.concat(raw_parts, ignore_index=True)
    raw_df = raw_df[raw_df["family"].isin(config["experiment"]["families"])].copy()

    metric_rows = []
    for row in raw_df.itertuples(index=False):
        poset = FAMILIES[row.family](n=int(row.n), seed=int(row.seed))
        metrics = residual_metrics(poset)
        metric_rows.append(
            {
                "family": str(row.family),
                "n": int(row.n),
                "seed": int(row.seed),
                "log_H": float(row.log_H),
                "antichain_width": float(row.antichain_width),
                "comparable_fraction": float(row.comparable_fraction),
                "geo_dim_eff": float(row.geo_dim_eff),
                **metrics,
            }
        )

    df = pd.DataFrame(metric_rows)
    df["z_layer_count"] = zscore(df["layer_count"])
    df["z_mean_layer_gap"] = zscore(df["mean_layer_gap"])
    df["z_long_edge_fraction"] = zscore(df["long_edge_fraction"])
    df["z_adjacent_edge_fraction"] = zscore(df["adjacent_edge_fraction"])
    df["z_reduction_edge_density"] = zscore(df["reduction_edge_density"])

    # Higher index = deeper / more globally integrated.
    df["hierarchy_integration_index"] = (
        df["z_layer_count"]
        + df["z_mean_layer_gap"]
        + df["z_long_edge_fraction"]
        - df["z_adjacent_edge_fraction"]
        - df["z_reduction_edge_density"]
    ) / 5.0

    X = np.column_stack(
        [
            np.ones(len(df)),
            df["antichain_width"].to_numpy(dtype=float),
            df["comparable_fraction"].to_numpy(dtype=float),
            df["geo_dim_eff"].to_numpy(dtype=float),
        ]
    )
    y_logh = df["log_H"].to_numpy(dtype=float)
    x_hii = df["hierarchy_integration_index"].to_numpy(dtype=float)

    logh_resid = regression_residual(y_logh, X)
    hii_resid = regression_residual(x_hii, X)
    partial_corr = pearson_corr(hii_resid, logh_resid)
    pvalue = permutation_pvalue(hii_resid, logh_resid, n_perm=int(config["experiment"]["n_permutations"]), seed=int(config["experiment"]["seed"]))

    # Simple regression with hierarchy term added after controls
    X_full = np.column_stack([X, x_hii])
    beta_full, *_ = np.linalg.lstsq(X_full, y_logh, rcond=None)

    summary_df = pd.DataFrame(
        [
            {
                "n_samples": int(len(df)),
                "partial_corr_hii_vs_logH_given_width_comp_dim": float(partial_corr),
                "permutation_pvalue": float(pvalue),
                "beta_antichain_width": float(beta_full[1]),
                "beta_comparable_fraction": float(beta_full[2]),
                "beta_geo_dim_eff": float(beta_full[3]),
                "beta_hierarchy_integration_index": float(beta_full[4]),
            }
        ]
    )

    by_family_df = (
        df.groupby("family")
        .agg(
            mean_log_H=("log_H", "mean"),
            mean_hii=("hierarchy_integration_index", "mean"),
            mean_layer_count=("layer_count", "mean"),
            mean_mean_layer_gap=("mean_layer_gap", "mean"),
            mean_long_edge_fraction=("long_edge_fraction", "mean"),
            mean_adjacent_edge_fraction=("adjacent_edge_fraction", "mean"),
            count=("seed", "count"),
        )
        .reset_index()
    )

    df.to_csv(out_dir / "hierarchy_integration_validation_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "hierarchy_integration_validation_summary.csv", index=False, encoding="utf-8-sig")
    by_family_df.to_csv(out_dir / "hierarchy_integration_validation_by_family.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "hierarchy_integration_validation_summary.csv").as_posix())
    print((out_dir / "hierarchy_integration_validation_by_family.csv").as_posix())
    print()
    print(summary_df.to_string(index=False))
