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


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = x - x.mean()
    y = y - y.mean()
    denom = math.sqrt(float((x * x).sum()) * float((y * y).sum()))
    if denom <= 1e-12:
        return 0.0
    return float((x * y).sum() / denom)


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    xr = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    yr = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    return pearson_corr(xr, yr)


def permutation_pvalue(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    observed = abs(pearson_corr(x, y))
    count = 1
    for _ in range(n_perm):
        shuffled = rng.permutation(y)
        if abs(pearson_corr(x, shuffled)) >= observed:
            count += 1
    return float(count / (n_perm + 1))


def regression_residual(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def build_design_matrix(
    df: pd.DataFrame,
    controls: list[str],
    add_n_fe: bool,
    add_family_fe: bool,
) -> tuple[np.ndarray, list[str]]:
    parts = [np.ones((len(df), 1), dtype=float)]
    names = ["intercept"]

    if controls:
        parts.append(df[controls].to_numpy(dtype=float))
        names.extend(controls)

    if add_n_fe:
        n_dummies = pd.get_dummies(df["n"].astype(str), prefix="n", drop_first=True, dtype=float)
        if not n_dummies.empty:
            parts.append(n_dummies.to_numpy(dtype=float))
            names.extend(n_dummies.columns.tolist())

    if add_family_fe:
        fam_dummies = pd.get_dummies(df["family"], prefix="family", drop_first=True, dtype=float)
        if not fam_dummies.empty:
            parts.append(fam_dummies.to_numpy(dtype=float))
            names.extend(fam_dummies.columns.tolist())

    return np.column_stack(parts), names


def run_partial_regression(
    df: pd.DataFrame,
    feature: str,
    controls: list[str],
    add_n_fe: bool,
    add_family_fe: bool,
    n_perm: int,
    seed: int,
) -> dict:
    X, names = build_design_matrix(df, controls, add_n_fe, add_family_fe)
    y = df["log_H"].to_numpy(dtype=float)
    x = df[feature].to_numpy(dtype=float)

    y_resid = regression_residual(y, X)
    x_resid = regression_residual(x, X)
    partial = pearson_corr(x_resid, y_resid)
    spearman = spearman_corr(x_resid, y_resid)
    pvalue = permutation_pvalue(x_resid, y_resid, n_perm=n_perm, seed=seed)

    X_full = np.column_stack([X, x])
    beta_full, *_ = np.linalg.lstsq(X_full, y, rcond=None)

    return {
        "feature": feature,
        "n_samples": int(len(df)),
        "partial_corr": float(partial),
        "spearman_partial_corr": float(spearman),
        "permutation_pvalue": float(pvalue),
        "beta_feature": float(beta_full[-1]),
        "n_controls": int(len(names) - 1),
    }


def dedupe_augmented_rows(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["family", "n", "seed", "log_H", "antichain_width", "comparable_fraction"]
    for col in ["geo_dim_eff_x", "geo_dim_eff_y", "geo_interval_shape_x", "geo_interval_shape_y"]:
        if col in df.columns:
            cols.append(col)
    out = df[cols].drop_duplicates(["family", "n", "seed"]).copy()

    if "geo_dim_eff" not in out.columns:
        if "geo_dim_eff_y" in out.columns:
            out["geo_dim_eff"] = out["geo_dim_eff_y"]
        elif "geo_dim_eff_x" in out.columns:
            out["geo_dim_eff"] = out["geo_dim_eff_x"]

    if "geo_interval_shape" not in out.columns:
        if "geo_interval_shape_y" in out.columns:
            out["geo_interval_shape"] = out["geo_interval_shape_y"]
        elif "geo_interval_shape_x" in out.columns:
            out["geo_interval_shape"] = out["geo_interval_shape_x"]

    keep = [
        "family",
        "n",
        "seed",
        "log_H",
        "antichain_width",
        "comparable_fraction",
        "geo_dim_eff",
        "geo_interval_shape",
    ]
    return out[keep].copy()


def build_analysis_df(augmented_df: pd.DataFrame, feature_cache_df: pd.DataFrame) -> pd.DataFrame:
    base_df = dedupe_augmented_rows(augmented_df)
    merged = base_df.merge(
        feature_cache_df,
        on=["family", "n", "seed"],
        how="inner",
        validate="one_to_one",
        suffixes=("_base", ""),
    )

    for col in ["antichain_width", "comparable_fraction", "geo_dim_eff", "geo_interval_shape"]:
        base_col = f"{col}_base"
        if col not in merged.columns and base_col in merged.columns:
            merged[col] = merged[base_col]

    merged["z_layer_count"] = zscore(merged["layer_count"])
    merged["z_mean_layer_gap"] = zscore(merged["mean_layer_gap"])
    merged["z_long_edge_fraction"] = zscore(merged["long_edge_fraction"])
    merged["z_adjacent_edge_fraction"] = zscore(merged["adjacent_edge_fraction"])
    merged["z_reduction_edge_density"] = zscore(merged["reduction_edge_density"])
    merged["hierarchy_integration_index"] = (
        merged["z_layer_count"]
        + merged["z_mean_layer_gap"]
        + merged["z_long_edge_fraction"]
        - merged["z_adjacent_edge_fraction"]
        - merged["z_reduction_edge_density"]
    ) / 5.0
    return merged


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run pooled and stratified Prediction C regressions from augmented BRONZE cache."
    )
    parser.add_argument(
        "--config",
        default="config_prediction_c_pooled_regression.yaml",
        help="Path to YAML config.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)

    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    augmented_df = pd.read_csv(config["input"]["augmented_combined_csv"], low_memory=False)
    feature_cache_df = pd.read_csv(config["input"]["feature_cache_csv"], low_memory=False)
    df = build_analysis_df(augmented_df, feature_cache_df)

    controls = list(config["experiment"]["controls"])
    features = list(config["experiment"]["features"])
    n_perm = int(config["experiment"]["n_permutations"])
    seed = int(config["experiment"]["seed"])

    model_specs = [
        ("controls_only", False, False),
        ("controls_plus_n_fe", True, False),
        ("controls_plus_family_fe", False, True),
        ("controls_plus_n_family_fe", True, True),
    ]

    summary_rows = []
    for spec_name, add_n_fe, add_family_fe in model_specs:
        for feature in features:
            result = run_partial_regression(
                df=df,
                feature=feature,
                controls=controls,
                add_n_fe=add_n_fe,
                add_family_fe=add_family_fe,
                n_perm=n_perm,
                seed=seed,
            )
            result["model_spec"] = spec_name
            summary_rows.append(result)

    summary_df = pd.DataFrame(summary_rows)

    family_rows = []
    for family, sub in df.groupby("family", sort=True):
        if len(sub) < max(8, len(controls) + 3):
            continue
        for feature in features:
            result = run_partial_regression(
                df=sub,
                feature=feature,
                controls=controls,
                add_n_fe=True,
                add_family_fe=False,
                n_perm=n_perm,
                seed=seed,
            )
            result["family"] = family
            family_rows.append(result)
    family_df = pd.DataFrame(family_rows)

    n_rows = []
    for n_value, sub in df.groupby("n", sort=True):
        if len(sub) < max(8, len(controls) + 3):
            continue
        for feature in features:
            result = run_partial_regression(
                df=sub,
                feature=feature,
                controls=controls,
                add_n_fe=False,
                add_family_fe=True,
                n_perm=n_perm,
                seed=seed,
            )
            result["n"] = int(n_value)
            n_rows.append(result)
    n_df = pd.DataFrame(n_rows)

    family_means_df = (
        df.groupby("family", sort=True)
        .agg(
            n_samples=("seed", "count"),
            mean_log_H=("log_H", "mean"),
            mean_hierarchy_integration_index=("hierarchy_integration_index", "mean"),
            mean_layer_count=("layer_count", "mean"),
            mean_mean_layer_gap=("mean_layer_gap", "mean"),
            mean_antichain_width=("antichain_width", "mean"),
            mean_comparable_fraction=("comparable_fraction", "mean"),
            mean_geo_dim_eff=("geo_dim_eff", "mean"),
        )
        .reset_index()
    )

    n_means_df = (
        df.groupby("n", sort=True)
        .agg(
            n_samples=("seed", "count"),
            mean_log_H=("log_H", "mean"),
            mean_hierarchy_integration_index=("hierarchy_integration_index", "mean"),
            mean_layer_count=("layer_count", "mean"),
            mean_mean_layer_gap=("mean_layer_gap", "mean"),
        )
        .reset_index()
    )

    df.to_csv(out_dir / "prediction_c_pooled_regression_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "prediction_c_pooled_regression_summary.csv", index=False, encoding="utf-8-sig")
    family_df.to_csv(out_dir / "prediction_c_pooled_regression_by_family.csv", index=False, encoding="utf-8-sig")
    n_df.to_csv(out_dir / "prediction_c_pooled_regression_by_n.csv", index=False, encoding="utf-8-sig")
    family_means_df.to_csv(out_dir / "prediction_c_pooled_regression_family_means.csv", index=False, encoding="utf-8-sig")
    n_means_df.to_csv(out_dir / "prediction_c_pooled_regression_n_means.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "prediction_c_pooled_regression_summary.csv").as_posix())
    print((out_dir / "prediction_c_pooled_regression_by_family.csv").as_posix())
    print((out_dir / "prediction_c_pooled_regression_by_n.csv").as_posix())
    print()
    print(summary_df.to_string(index=False))
