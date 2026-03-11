from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from experiment import load_config
from normalization import add_normalized_columns


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["n", "gamma", "family", "action_mode", "normalization_method"])
        .agg(
            mean_score_norm=("score_norm", "mean"),
            std_score_norm=("score_norm", "std"),
            count=("score_norm", "count"),
        )
        .reset_index()
        .sort_values(["normalization_method", "action_mode", "n", "gamma", "mean_score_norm"])
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare normalization methods on raw_samples.csv.")
    parser.add_argument("--config", default="config_smallN_exact.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    raw_df = pd.read_csv(out_dir / config["output"]["raw_samples_csv"])

    robust_df = add_normalized_columns(raw_df, method="robust_zscore", group_cols=("n",))
    z_df = add_normalized_columns(raw_df, method="zscore", group_cols=("n",))

    combined = pd.concat([robust_df, z_df], ignore_index=True)
    summary_df = summarize(combined)
    out_path = out_dir / "normalization_sensitivity.csv"
    summary_df.to_csv(out_path, index=False, encoding="utf-8-sig")
    print(summary_df.to_string(index=False))
