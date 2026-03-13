from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


PAIRWISE_COLUMNS = [
    "n",
    "mlr_seed",
    "lor_seed",
    "log_H_delta_mlr_minus_lor2d",
    "score_A2_gamma_delta_mlr_minus_lor2d",
]

RESIDUAL_METRICS = [
    "layer_count",
    "mean_layer_gap",
    "adjacent_edge_fraction",
    "long_edge_fraction",
    "reduction_edge_density",
    "cover_density",
    "layer_signature_redundancy",
]


def load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def permutation_pvalue(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    observed = np.corrcoef(x, y)[0, 1]
    count = 0
    for _ in range(n_perm):
        y_perm = rng.permutation(y)
        corr = np.corrcoef(x, y_perm)[0, 1]
        if abs(corr) >= abs(observed):
            count += 1
    return (count + 1) / (n_perm + 1)


def build_pair_delta_table(pairwise_df: pd.DataFrame, residual_df: pd.DataFrame) -> pd.DataFrame:
    mlr = residual_df.loc[residual_df["family"] == "multi_layer_random"].copy()
    lor = residual_df.loc[residual_df["family"] == "lorentzian_like_2d"].copy()

    mlr = mlr.rename(columns={"seed": "mlr_seed"})
    lor = lor.rename(columns={"seed": "lor_seed"})

    mlr_keep = ["n", "mlr_seed"] + RESIDUAL_METRICS
    lor_keep = ["n", "lor_seed"] + RESIDUAL_METRICS

    merged = pairwise_df[PAIRWISE_COLUMNS].merge(
        mlr[mlr_keep],
        on=["n", "mlr_seed"],
        how="left",
    ).merge(
        lor[lor_keep],
        on=["n", "lor_seed"],
        how="left",
        suffixes=("_mlr", "_lor2d"),
    )

    for metric in RESIDUAL_METRICS:
        merged[f"{metric}_delta_mlr_minus_lor2d"] = (
            merged[f"{metric}_mlr"] - merged[f"{metric}_lor2d"]
        )

    merged["locality_dominance_delta_index"] = (
        merged["adjacent_edge_fraction_delta_mlr_minus_lor2d"]
        + merged["reduction_edge_density_delta_mlr_minus_lor2d"]
        - merged["layer_count_delta_mlr_minus_lor2d"]
        - merged["mean_layer_gap_delta_mlr_minus_lor2d"]
        - merged["long_edge_fraction_delta_mlr_minus_lor2d"]
    )

    merged["hierarchy_integration_delta_index"] = (
        merged["layer_count_delta_mlr_minus_lor2d"]
        + merged["mean_layer_gap_delta_mlr_minus_lor2d"]
        + merged["long_edge_fraction_delta_mlr_minus_lor2d"]
        + merged["reduction_edge_density_delta_mlr_minus_lor2d"]
        - merged["adjacent_edge_fraction_delta_mlr_minus_lor2d"]
    )

    return merged


def summarize_relationship(
    df: pd.DataFrame,
    feature_col: str,
    target_col: str,
    n_perm: int,
    seed: int,
) -> dict:
    x = df[feature_col].to_numpy(dtype=float)
    y = df[target_col].to_numpy(dtype=float)
    pearson = float(np.corrcoef(x, y)[0, 1])
    spearman = float(pd.Series(x).corr(pd.Series(y), method="spearman"))
    p_perm = permutation_pvalue(x, y, n_perm=n_perm, seed=seed)
    return {
        "feature": feature_col,
        "target": target_col,
        "n_pairs": int(len(df)),
        "pearson_corr": pearson,
        "spearman_corr": spearman,
        "permutation_pvalue": float(p_perm),
        "mean_feature": float(np.mean(x)),
        "mean_target": float(np.mean(y)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config_path = Path(args.config)
    config = load_config(config_path)

    pairwise_csv = Path(config["pairwise_csv"])
    residual_csv = Path(config["residual_csv"])
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    n_perm = int(config.get("n_perm", 2000))
    seed = int(config.get("seed", 42))

    pairwise_df = pd.read_csv(pairwise_csv)
    residual_df = pd.read_csv(residual_csv)

    delta_df = build_pair_delta_table(pairwise_df, residual_df)

    summary_rows = []
    component_rows = []
    feature_cols = [
        "locality_dominance_delta_index",
        "hierarchy_integration_delta_index",
        "layer_count_delta_mlr_minus_lor2d",
        "mean_layer_gap_delta_mlr_minus_lor2d",
        "adjacent_edge_fraction_delta_mlr_minus_lor2d",
        "long_edge_fraction_delta_mlr_minus_lor2d",
        "reduction_edge_density_delta_mlr_minus_lor2d",
        "cover_density_delta_mlr_minus_lor2d",
        "layer_signature_redundancy_delta_mlr_minus_lor2d",
    ]
    target_cols = [
        "log_H_delta_mlr_minus_lor2d",
        "score_A2_gamma_delta_mlr_minus_lor2d",
    ]

    for target in target_cols:
        for feature in feature_cols:
            row = summarize_relationship(
                delta_df, feature, target, n_perm=n_perm, seed=seed
            )
            summary_rows.append(row)

    for feature in feature_cols:
        component_rows.append(
            {
                "feature": feature,
                "mean": float(delta_df[feature].mean()),
                "std": float(delta_df[feature].std(ddof=1)),
                "min": float(delta_df[feature].min()),
                "max": float(delta_df[feature].max()),
            }
        )

    delta_df.to_csv(
        output_dir / "prediction_c_pairwise_validation_raw.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(summary_rows).to_csv(
        output_dir / "prediction_c_pairwise_validation_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(component_rows).to_csv(
        output_dir / "prediction_c_pairwise_validation_components.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print(f"Wrote outputs to {output_dir}")


if __name__ == "__main__":
    main()
