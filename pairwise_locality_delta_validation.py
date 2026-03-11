from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def zscore(series: pd.Series) -> pd.Series:
    values = series.astype(float)
    std = float(values.std(ddof=0))
    if std == 0.0 or np.isnan(std):
        return pd.Series(np.zeros(len(values)), index=series.index, dtype=float)
    mean = float(values.mean())
    return (values - mean) / std


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    xr = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    yr = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    return pearson_corr(xr, yr)


def permutation_pvalue(x: np.ndarray, y: np.ndarray, n_permutations: int, seed: int) -> float:
    observed = abs(pearson_corr(x, y))
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_permutations):
        shuffled = rng.permutation(y)
        if abs(pearson_corr(x, shuffled)) >= observed:
            hits += 1
    return float((hits + 1) / (n_permutations + 1))


def build_pair_delta_df(raw_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["n", "pair_mlr_seed", "pair_lor_seed"]
    for _, sub in raw_df.groupby(group_cols):
        mlr = sub[sub["family"] == "multi_layer_random"].iloc[0]
        lor = sub[sub["family"] == "lorentzian_like_2d"].iloc[0]
        rows.append(
            {
                "n": int(mlr["n"]),
                "mlr_seed": int(mlr["pair_mlr_seed"]),
                "lor_seed": int(mlr["pair_lor_seed"]),
                "adjacent_edge_fraction_delta": float(mlr["adjacent_edge_fraction"] - lor["adjacent_edge_fraction"]),
                "reduction_edge_density_delta": float(mlr["reduction_edge_density"] - lor["reduction_edge_density"]),
                "layer_count_delta": float(mlr["layer_count"] - lor["layer_count"]),
                "mean_layer_gap_delta": float(mlr["mean_layer_gap"] - lor["mean_layer_gap"]),
                "long_edge_fraction_delta": float(mlr["long_edge_fraction"] - lor["long_edge_fraction"]),
                "cover_density_delta": float(mlr["cover_density"] - lor["cover_density"]),
                "layer_signature_redundancy_delta": float(
                    mlr["layer_signature_redundancy"] - lor["layer_signature_redundancy"]
                ),
            }
        )
    return pd.DataFrame(rows)


def locality_index(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    # Use raw matched-pair deltas here. This reproduces the earlier inline
    # diagnostic where "shallower + more adjacent-layer-dominated" pairs tracked
    # larger log_H deltas most clearly.
    work["locality_dominance_delta_index"] = (
        work["adjacent_edge_fraction_delta"]
        + work["reduction_edge_density_delta"]
        - work["layer_count_delta"]
        - work["mean_layer_gap_delta"]
        - work["long_edge_fraction_delta"]
    )
    work["locality_dominance_delta_index_zscore"] = zscore(work["locality_dominance_delta_index"])
    return work


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate whether pairwise locality-dominance deltas track log_H deltas."
    )
    parser.add_argument("--config", default="config_pairwise_locality_delta_validation.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    pairs_df = pd.read_csv(config["input"]["pairs_csv"])
    raw_df = pd.read_csv(config["input"]["residual_raw_csv"])

    delta_df = build_pair_delta_df(raw_df)
    merged = pairs_df.merge(delta_df, how="inner", on=["n", "mlr_seed", "lor_seed"])
    merged = locality_index(merged)

    n_perm = int(config["experiment"].get("n_permutations", 2000))
    seed = int(config["experiment"].get("seed", 20260311))

    x = merged["locality_dominance_delta_index"].to_numpy(dtype=float)
    xz = merged["locality_dominance_delta_index_zscore"].to_numpy(dtype=float)
    y = merged["log_H_delta_mlr_minus_lor2d"].to_numpy(dtype=float)
    summary_rows = [
        {
            "scope": "all_pairs",
            "n": "all",
            "count": int(len(merged)),
            "pearson_corr": pearson_corr(x, y),
            "pearson_corr_zscore_index": pearson_corr(xz, y),
            "spearman_corr": spearman_corr(x, y),
            "permutation_pvalue": permutation_pvalue(x, y, n_perm, seed),
            "mean_locality_dominance_delta_index": float(np.mean(x)),
            "mean_log_H_delta_mlr_minus_lor2d": float(np.mean(y)),
        }
    ]
    for n, sub in merged.groupby("n", sort=True):
        xi = sub["locality_dominance_delta_index"].to_numpy(dtype=float)
        yi = sub["log_H_delta_mlr_minus_lor2d"].to_numpy(dtype=float)
        summary_rows.append(
            {
                "scope": "by_n",
                "n": int(n),
                "count": int(len(sub)),
                "pearson_corr": pearson_corr(xi, yi),
                "pearson_corr_zscore_index": pearson_corr(
                    sub["locality_dominance_delta_index_zscore"].to_numpy(dtype=float), yi
                ),
                "spearman_corr": spearman_corr(xi, yi),
                "permutation_pvalue": permutation_pvalue(xi, yi, n_perm, seed + int(n)),
                "mean_locality_dominance_delta_index": float(np.mean(xi)),
                "mean_log_H_delta_mlr_minus_lor2d": float(np.mean(yi)),
            }
        )
    summary_df = pd.DataFrame(summary_rows)

    component_summary = (
        merged.groupby("n")
        .agg(
            mean_adjacent_edge_fraction_delta=("adjacent_edge_fraction_delta", "mean"),
            mean_reduction_edge_density_delta=("reduction_edge_density_delta", "mean"),
            mean_layer_count_delta=("layer_count_delta", "mean"),
            mean_mean_layer_gap_delta=("mean_layer_gap_delta", "mean"),
            mean_long_edge_fraction_delta=("long_edge_fraction_delta", "mean"),
            mean_locality_dominance_delta_index=("locality_dominance_delta_index", "mean"),
            mean_log_H_delta_mlr_minus_lor2d=("log_H_delta_mlr_minus_lor2d", "mean"),
            count=("n", "count"),
        )
        .reset_index()
    )

    merged.to_csv(out_dir / "pairwise_locality_delta_validation_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "pairwise_locality_delta_validation_summary.csv", index=False, encoding="utf-8-sig")
    component_summary.to_csv(
        out_dir / "pairwise_locality_delta_validation_components.csv", index=False, encoding="utf-8-sig"
    )

    print((out_dir / "pairwise_locality_delta_validation_summary.csv").as_posix())
    print((out_dir / "pairwise_locality_delta_validation_components.csv").as_posix())
    print()
    print(summary_df.to_string(index=False))
