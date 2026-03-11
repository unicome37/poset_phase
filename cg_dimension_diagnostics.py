from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import pandas as pd
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DIMENSION_FAMILIES = [
    "lorentzian_like_2d",
    "lorentzian_like_3d",
    "lorentzian_like_4d",
]


def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def summarize_dimension_bias(summary_df: pd.DataFrame, rank_df: pd.DataFrame) -> pd.DataFrame:
    df = rank_df[rank_df["family"].isin(DIMENSION_FAMILIES)].copy()
    grouped = (
        df.groupby("family")
        .agg(
            mean_score_local=("mean_score_local", "mean"),
            mean_penalty_cg=("mean_penalty_cg", "mean"),
            mean_score_augmented=("mean_score_augmented", "mean"),
            mean_self_drift=("mean_self_drift", "mean"),
            family_switch_rate=("family_switch_rate", "mean"),
            mean_rank_shift=("mean_rank_shift", "mean"),
            mean_gc_penalty=("mean_penalty_global_consistency", "mean"),
            best_rank_local=("rank_local", "min"),
            best_rank_augmented=("rank_augmented_cg", "min"),
            mean_rank_local=("rank_local", "mean"),
            mean_rank_augmented=("rank_augmented_cg", "mean"),
        )
        .reset_index()
    )
    grouped["cg_lift"] = grouped["mean_score_augmented"] - grouped["mean_score_local"]
    grouped["local_advantage_over_2d"] = (
        grouped["mean_score_local"]
        - grouped.loc[grouped["family"] == "lorentzian_like_2d", "mean_score_local"].iloc[0]
    )
    grouped["augmented_advantage_over_2d"] = (
        grouped["mean_score_augmented"]
        - grouped.loc[grouped["family"] == "lorentzian_like_2d", "mean_score_augmented"].iloc[0]
    )
    return grouped.sort_values("mean_rank_augmented")


def summarize_by_regime(rank_df: pd.DataFrame) -> pd.DataFrame:
    df = rank_df[rank_df["family"].isin(DIMENSION_FAMILIES)].copy()
    out = (
        df.groupby(["n", "keep_ratio", "family"])
        .agg(
            mean_score_local=("mean_score_local", "mean"),
            mean_penalty_cg=("mean_penalty_cg", "mean"),
            mean_score_augmented=("mean_score_augmented", "mean"),
            mean_self_drift=("mean_self_drift", "mean"),
            family_switch_rate=("family_switch_rate", "mean"),
            mean_rank_local=("rank_local", "mean"),
            mean_rank_augmented=("rank_augmented_cg", "mean"),
        )
        .reset_index()
        .sort_values(["n", "keep_ratio", "mean_rank_augmented", "family"])
    )
    return out


def plot_dimension_bias(summary: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].bar(summary["family"], summary["mean_score_local"])
    axes[0].set_title("Mean Local Score")
    axes[0].tick_params(axis="x", rotation=20)

    axes[1].bar(summary["family"], summary["mean_penalty_cg"])
    axes[1].set_title("Mean CG Penalty")
    axes[1].tick_params(axis="x", rotation=20)

    axes[2].bar(summary["family"], summary["mean_rank_augmented"])
    axes[2].set_title("Mean Rank After CG")
    axes[2].tick_params(axis="x", rotation=20)

    fig.tight_layout()
    fig.savefig(out_dir / "dimension_bias_overview.png", dpi=160)
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Diagnose dimension bias in coarse-grain runs.")
    parser.add_argument("--config", default="config_cg_dimension.yaml", help="Path to YAML config file.")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    summary_path = out_dir / "cg_summary.csv"
    rank_path = out_dir / "cg_rank_summary.csv"

    summary_df = pd.read_csv(summary_path)
    rank_df = pd.read_csv(rank_path)

    dim_summary = summarize_dimension_bias(summary_df, rank_df)
    regime_summary = summarize_by_regime(rank_df)

    dim_summary.to_csv(out_dir / "dimension_bias_summary.csv", index=False, encoding="utf-8-sig")
    regime_summary.to_csv(out_dir / "dimension_bias_by_regime.csv", index=False, encoding="utf-8-sig")

    plot_dimension_bias(dim_summary, out_dir / "plots")

    print("Dimension bias summary:")
    print(dim_summary.to_string(index=False))
    print()
    print("By regime:")
    print(regime_summary.to_string(index=False))


if __name__ == "__main__":
    main()
