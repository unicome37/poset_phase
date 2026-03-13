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

DELTA_METRICS = [
    "log_H",
    "layer_count",
    "mean_layer_gap",
    "adjacent_edge_fraction",
    "long_edge_fraction",
    "reduction_edge_density",
    "cover_density",
    "layer_signature_redundancy",
]


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run large-n matched-pair Prediction C validation from augmented BRONZE cache."
    )
    parser.add_argument(
        "--config",
        default="config_prediction_c_bronze_matched_validation.yaml",
        help="Path to YAML config file.",
    )
    return parser


def standardized_distance_matrix(left: pd.DataFrame, right: pd.DataFrame, cols: list[str]) -> np.ndarray:
    combined = pd.concat([left[cols], right[cols]], ignore_index=True)
    means = combined.mean()
    stds = combined.std(ddof=0).replace(0.0, 1.0)
    left_z = (left[cols] - means) / stds
    right_z = (right[cols] - means) / stds
    diff = left_z.to_numpy()[:, None, :] - right_z.to_numpy()[None, :, :]
    return np.sqrt((diff * diff).sum(axis=2))


def greedy_match(left: pd.DataFrame, right: pd.DataFrame, cols: list[str]) -> list[tuple[int, int, float]]:
    dist = standardized_distance_matrix(left, right, cols)
    candidates: list[tuple[float, int, int]] = []
    for i in range(dist.shape[0]):
        for j in range(dist.shape[1]):
            candidates.append((float(dist[i, j]), i, j))
    candidates.sort(key=lambda x: x[0])

    used_left: set[int] = set()
    used_right: set[int] = set()
    pairs: list[tuple[int, int, float]] = []

    for d, i, j in candidates:
        if i in used_left or j in used_right:
            continue
        used_left.add(i)
        used_right.add(j)
        pairs.append((i, j, d))
    return pairs


def permutation_pvalue(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    observed = np.corrcoef(x, y)[0, 1]
    count = 0
    for _ in range(n_perm):
        corr = np.corrcoef(x, rng.permutation(y))[0, 1]
        if abs(corr) >= abs(observed):
            count += 1
    return (count + 1) / (n_perm + 1)


def summarize_relationship(df: pd.DataFrame, family: str, feature_col: str, target_col: str, n_perm: int, seed: int) -> dict:
    x = df[feature_col].to_numpy(dtype=float)
    y = df[target_col].to_numpy(dtype=float)
    return {
        "family": family,
        "feature": feature_col,
        "target": target_col,
        "n_pairs": int(len(df)),
        "pearson_corr": float(np.corrcoef(x, y)[0, 1]),
        "spearman_corr": float(pd.Series(x).corr(pd.Series(y), method="spearman")),
        "permutation_pvalue": float(permutation_pvalue(x, y, n_perm=n_perm, seed=seed)),
        "mean_feature": float(np.mean(x)),
        "mean_target": float(np.mean(y)),
    }


def add_delta_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for metric in DELTA_METRICS:
        out[f"{metric}_delta_vs_lor2d"] = out[f"{metric}_cmp"] - out[f"{metric}_lor2d"]

    out["locality_dominance_delta_index"] = (
        out["adjacent_edge_fraction_delta_vs_lor2d"]
        + out["reduction_edge_density_delta_vs_lor2d"]
        - out["layer_count_delta_vs_lor2d"]
        - out["mean_layer_gap_delta_vs_lor2d"]
        - out["long_edge_fraction_delta_vs_lor2d"]
    )
    out["hierarchy_integration_delta_index"] = (
        out["layer_count_delta_vs_lor2d"]
        + out["mean_layer_gap_delta_vs_lor2d"]
        + out["long_edge_fraction_delta_vs_lor2d"]
        + out["reduction_edge_density_delta_vs_lor2d"]
        - out["adjacent_edge_fraction_delta_vs_lor2d"]
    )
    return out


def unique_poset_table(df: pd.DataFrame) -> pd.DataFrame:
    geo_dim_col = "geo_dim_eff"
    if geo_dim_col not in df.columns:
        geo_dim_col = "geo_dim_eff_y" if "geo_dim_eff_y" in df.columns else "geo_dim_eff_x"

    geo_interval_shape_col = "geo_interval_shape"
    if geo_interval_shape_col not in df.columns:
        geo_interval_shape_col = (
            "geo_interval_shape_y" if "geo_interval_shape_y" in df.columns else "geo_interval_shape_x"
        )

    keep_cols = [
        "family",
        "n",
        "seed",
        "log_H",
        "antichain_width",
        "comparable_fraction",
        geo_dim_col,
        geo_interval_shape_col,
        "layer_count",
        "mean_layer_gap",
        "adjacent_edge_fraction",
        "long_edge_fraction",
        "reduction_edge_density",
        "cover_density",
        "layer_signature_redundancy",
    ]
    out = df[keep_cols].drop_duplicates(subset=["family", "n", "seed"]).reset_index(drop=True)
    return out.rename(
        columns={
            geo_dim_col: "geo_dim_eff",
            geo_interval_shape_col: "geo_interval_shape",
        }
    )


def build_pair_table(unique_df: pd.DataFrame, families: list[str], reference_family: str) -> pd.DataFrame:
    rows: list[dict] = []

    for family in families:
        cmp_all = unique_df.loc[unique_df["family"] == family].copy()
        ref_all = unique_df.loc[unique_df["family"] == reference_family].copy()
        for n in sorted(set(cmp_all["n"].tolist()) & set(ref_all["n"].tolist())):
            cmp_df = cmp_all.loc[cmp_all["n"] == n].reset_index(drop=True)
            ref_df = ref_all.loc[ref_all["n"] == n].reset_index(drop=True)
            if cmp_df.empty or ref_df.empty:
                continue
            pairs = greedy_match(cmp_df, ref_df, MATCH_FEATURES)
            for i, j, d in pairs:
                cmp_row = cmp_df.iloc[i]
                ref_row = ref_df.iloc[j]
                row = {
                    "family": family,
                    "n": int(n),
                    "cmp_seed": int(cmp_row["seed"]),
                    "lor2d_seed": int(ref_row["seed"]),
                    "match_distance": float(d),
                }
                for metric in DELTA_METRICS:
                    row[f"{metric}_cmp"] = float(cmp_row[metric])
                    row[f"{metric}_lor2d"] = float(ref_row[metric])
                row["antichain_width_cmp"] = float(cmp_row["antichain_width"])
                row["antichain_width_lor2d"] = float(ref_row["antichain_width"])
                row["comparable_fraction_cmp"] = float(cmp_row["comparable_fraction"])
                row["comparable_fraction_lor2d"] = float(ref_row["comparable_fraction"])
                row["geo_dim_eff_cmp"] = float(cmp_row["geo_dim_eff"])
                row["geo_dim_eff_lor2d"] = float(ref_row["geo_dim_eff"])
                row["geo_interval_shape_cmp"] = float(cmp_row["geo_interval_shape"])
                row["geo_interval_shape_lor2d"] = float(ref_row["geo_interval_shape"])
                rows.append(row)

    return add_delta_columns(pd.DataFrame(rows))


def main() -> None:
    args = build_arg_parser().parse_args()
    config_path = Path(args.config)
    config = load_config(config_path)

    input_csv = Path(config["input_csv"])
    if not input_csv.is_absolute():
        input_csv = (config_path.parent / input_csv).resolve()
    output_dir = Path(config["output_dir"])
    if not output_dir.is_absolute():
        output_dir = (config_path.parent / output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_family = str(config.get("reference_family", "lorentzian_like_2d"))
    families = [str(v) for v in config["families"]]
    n_perm = int(config.get("n_perm", 5000))
    seed = int(config.get("seed", 42))

    augmented_df = pd.read_csv(input_csv, low_memory=False)
    unique_df = unique_poset_table(augmented_df)
    pair_df = build_pair_table(unique_df, families=families, reference_family=reference_family)

    summary_rows: list[dict] = []
    feature_cols = [
        "hierarchy_integration_delta_index",
        "locality_dominance_delta_index",
        "layer_count_delta_vs_lor2d",
        "mean_layer_gap_delta_vs_lor2d",
        "adjacent_edge_fraction_delta_vs_lor2d",
        "long_edge_fraction_delta_vs_lor2d",
        "reduction_edge_density_delta_vs_lor2d",
        "cover_density_delta_vs_lor2d",
        "layer_signature_redundancy_delta_vs_lor2d",
    ]
    for family in families:
        sub = pair_df.loc[pair_df["family"] == family].copy()
        for feature in feature_cols:
            summary_rows.append(
                summarize_relationship(
                    sub,
                    family=family,
                    feature_col=feature,
                    target_col="log_H_delta_vs_lor2d",
                    n_perm=n_perm,
                    seed=seed,
                )
            )

    pair_df.to_csv(output_dir / "prediction_c_bronze_matched_pairs.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(summary_rows).to_csv(
        output_dir / "prediction_c_bronze_matched_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    count_df = (
        pair_df.groupby(["family", "n"])
        .size()
        .reset_index(name="n_pairs")
        .sort_values(["family", "n"])
    )
    count_df.to_csv(
        output_dir / "prediction_c_bronze_matched_counts.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print((output_dir / "prediction_c_bronze_matched_summary.csv").as_posix())
    print((output_dir / "prediction_c_bronze_matched_counts.csv").as_posix())
    print()
    print(count_df.to_string(index=False))


if __name__ == "__main__":
    main()
