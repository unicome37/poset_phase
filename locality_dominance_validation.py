from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


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
    parser = argparse.ArgumentParser(description="Validate whether locality-dominance predicts higher log_H after controls.")
    parser.add_argument("--config", default="config_locality_dominance_validation.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(config["input"]["raw_csv"])

    df["z_adjacent_edge_fraction"] = zscore(df["adjacent_edge_fraction"])
    df["z_reduction_edge_density"] = zscore(df["reduction_edge_density"])
    df["z_layer_count"] = zscore(df["layer_count"])
    df["z_mean_layer_gap"] = zscore(df["mean_layer_gap"])
    df["z_long_edge_fraction"] = zscore(df["long_edge_fraction"])

    # Higher = shallower, more local, more adjacent-edge dominated.
    df["locality_dominance_index"] = (
        df["z_adjacent_edge_fraction"]
        + df["z_reduction_edge_density"]
        - df["z_layer_count"]
        - df["z_mean_layer_gap"]
        - df["z_long_edge_fraction"]
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
    x_ldi = df["locality_dominance_index"].to_numpy(dtype=float)

    logh_resid = regression_residual(y_logh, X)
    ldi_resid = regression_residual(x_ldi, X)
    partial_corr = pearson_corr(ldi_resid, logh_resid)
    pvalue = permutation_pvalue(ldi_resid, logh_resid, n_perm=int(config["experiment"]["n_permutations"]), seed=int(config["experiment"]["seed"]))

    X_full = np.column_stack([X, x_ldi])
    beta_full, *_ = np.linalg.lstsq(X_full, y_logh, rcond=None)

    summary_df = pd.DataFrame(
        [
            {
                "n_samples": int(len(df)),
                "partial_corr_ldi_vs_logH_given_width_comp_dim": float(partial_corr),
                "permutation_pvalue": float(pvalue),
                "beta_antichain_width": float(beta_full[1]),
                "beta_comparable_fraction": float(beta_full[2]),
                "beta_geo_dim_eff": float(beta_full[3]),
                "beta_locality_dominance_index": float(beta_full[4]),
            }
        ]
    )

    by_family_df = (
        df.groupby("family")
        .agg(
            mean_log_H=("log_H", "mean"),
            mean_ldi=("locality_dominance_index", "mean"),
            mean_adjacent_edge_fraction=("adjacent_edge_fraction", "mean"),
            mean_reduction_edge_density=("reduction_edge_density", "mean"),
            mean_layer_count=("layer_count", "mean"),
            mean_mean_layer_gap=("mean_layer_gap", "mean"),
            mean_long_edge_fraction=("long_edge_fraction", "mean"),
            count=("seed", "count"),
        )
        .reset_index()
    )

    df.to_csv(out_dir / "locality_dominance_validation_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "locality_dominance_validation_summary.csv", index=False, encoding="utf-8-sig")
    by_family_df.to_csv(out_dir / "locality_dominance_validation_by_family.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "locality_dominance_validation_summary.csv").as_posix())
    print((out_dir / "locality_dominance_validation_by_family.csv").as_posix())
    print()
    print(summary_df.to_string(index=False))
