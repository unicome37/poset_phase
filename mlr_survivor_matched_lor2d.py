from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


MATCH_FEATURES = [
    "antichain_width",
    "comparable_fraction",
    "geo_dim_eff",
    "geo_interval_shape",
]


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def standardized_distance_matrix(left: pd.DataFrame, right: pd.DataFrame, cols: list[str]) -> np.ndarray:
    combined = pd.concat([left[cols], right[cols]], ignore_index=True)
    means = combined.mean()
    stds = combined.std(ddof=0).replace(0.0, 1.0)
    left_z = (left[cols] - means) / stds
    right_z = (right[cols] - means) / stds
    diff = left_z.to_numpy()[:, None, :] - right_z.to_numpy()[None, :, :]
    return np.sqrt((diff * diff).sum(axis=2))


def greedy_match(mlr_df: pd.DataFrame, lor_df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    dist = standardized_distance_matrix(mlr_df, lor_df, cols)
    pairs = []
    used_mlr: set[int] = set()
    used_lor: set[int] = set()

    candidates = []
    for i in range(dist.shape[0]):
        for j in range(dist.shape[1]):
            candidates.append((float(dist[i, j]), i, j))
    candidates.sort(key=lambda x: x[0])

    for d, i, j in candidates:
        if i in used_mlr or j in used_lor:
            continue
        used_mlr.add(i)
        used_lor.add(j)
        pairs.append((i, j, d))
    rows = []
    for i, j, d in pairs:
        mlr = mlr_df.iloc[i]
        lor = lor_df.iloc[j]
        row = {
            "n": int(mlr["n"]),
            "mlr_seed": int(mlr["seed"]),
            "lor_seed": int(lor["seed"]),
            "match_distance": float(d),
        }
        for col in ["log_H", "score_A2_gamma", "geo_total", "cg_family_switch_rate", "cg_mean_penalty", "antichain_width", "comparable_fraction", "geo_dim_eff", "geo_interval_shape"]:
            row[f"{col}_mlr"] = float(mlr[col])
            row[f"{col}_lor2d"] = float(lor[col])
            row[f"{col}_delta_mlr_minus_lor2d"] = float(mlr[col] - lor[col])
        rows.append(row)
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Match MLR survivors to nearest Lor2D samples in the accepted window.")
    parser.add_argument("--config", default="config_mlr_survivor_matched_lor2d.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    parts = []
    for csv_path in config["input"]["raw_csvs"]:
        df = pd.read_csv(csv_path)
        parts.append(df)
    raw_df = pd.concat(parts, ignore_index=True)

    match_rows = []
    for n, sub in raw_df.groupby("n"):
        mlr = sub[sub["family"] == "multi_layer_random"].reset_index(drop=True)
        lor = sub[sub["family"] == "lorentzian_like_2d"].reset_index(drop=True)
        if mlr.empty or lor.empty:
            continue
        matched = greedy_match(mlr, lor, MATCH_FEATURES)
        match_rows.append(matched)

    match_df = pd.concat(match_rows, ignore_index=True)
    summary_df = (
        match_df.groupby("n")
        .agg(
            mean_match_distance=("match_distance", "mean"),
            mean_log_H_delta=("log_H_delta_mlr_minus_lor2d", "mean"),
            mean_score_A2_delta=("score_A2_gamma_delta_mlr_minus_lor2d", "mean"),
            mean_geo_total_delta=("geo_total_delta_mlr_minus_lor2d", "mean"),
            mean_switch_delta=("cg_family_switch_rate_delta_mlr_minus_lor2d", "mean"),
            mean_cg_penalty_delta=("cg_mean_penalty_delta_mlr_minus_lor2d", "mean"),
            mean_width_delta=("antichain_width_delta_mlr_minus_lor2d", "mean"),
            mean_comp_delta=("comparable_fraction_delta_mlr_minus_lor2d", "mean"),
            count=("mlr_seed", "count"),
        )
        .reset_index()
    )

    match_df.to_csv(out_dir / "mlr_survivor_matched_pairs.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "mlr_survivor_matched_summary.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "mlr_survivor_matched_pairs.csv").as_posix())
    print((out_dir / "mlr_survivor_matched_summary.csv").as_posix())
    print()
    print(summary_df.to_string(index=False))
