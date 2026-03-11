from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


BASE_WEIGHTS = {
    "geo_width_height": 2.0,
    "geo_dim_proxy_penalty": 8.0,
    "geo_comparability_window": 6.0,
    "geo_cover_density": 3.0,
    "geo_interval_profile": 5.0,
    "geo_interval_shape": 5.0,
    "geo_layer_smoothness": 2.0,
}


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze what remains between A2_full and width_height + dim_consistency variants."
    )
    parser.add_argument("--config", default="config_noncyclic_dim_replacement_gamma_c.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])

    raw_df = pd.read_csv(out_dir / "noncyclic_dim_replacement_raw.csv")
    scale_weight = float(raw_df["scale_matched_weight"].iloc[0])

    unique_cols = [
        "n",
        "family",
        "sample_id",
        "geo_width_height",
        "geo_dim_proxy_penalty",
        "geo_dim_consistency",
        "geo_comparability_window",
        "geo_cover_density",
        "geo_interval_profile",
        "geo_interval_shape",
        "geo_layer_smoothness",
    ]
    base = raw_df[unique_cols].drop_duplicates().copy()

    base["pen_width_height"] = BASE_WEIGHTS["geo_width_height"] * base["geo_width_height"]
    base["pen_dim_proxy"] = BASE_WEIGHTS["geo_dim_proxy_penalty"] * base["geo_dim_proxy_penalty"]
    base["pen_dim_consistency_scaled"] = scale_weight * base["geo_dim_consistency"]
    base["pen_comparability"] = BASE_WEIGHTS["geo_comparability_window"] * base["geo_comparability_window"]
    base["pen_cover"] = BASE_WEIGHTS["geo_cover_density"] * base["geo_cover_density"]
    base["pen_interval_profile"] = BASE_WEIGHTS["geo_interval_profile"] * base["geo_interval_profile"]
    base["pen_interval_shape"] = BASE_WEIGHTS["geo_interval_shape"] * base["geo_interval_shape"]
    base["pen_layer_smoothness"] = BASE_WEIGHTS["geo_layer_smoothness"] * base["geo_layer_smoothness"]

    base["pen_width_plus_consistency"] = base["pen_width_height"] + base["pen_dim_consistency_scaled"]
    base["pen_residual_prior_gap"] = (
        base["pen_dim_proxy"]
        + base["pen_comparability"]
        + base["pen_cover"]
        + base["pen_interval_profile"]
        + base["pen_interval_shape"]
        + base["pen_layer_smoothness"]
        - base["pen_dim_consistency_scaled"]
    )
    base["pen_interval_bundle"] = base["pen_interval_profile"] + base["pen_interval_shape"]
    base["pen_noninterval_residual"] = (
        base["pen_dim_proxy"]
        + base["pen_comparability"]
        + base["pen_cover"]
        + base["pen_layer_smoothness"]
        - base["pen_dim_consistency_scaled"]
    )

    family_summary = (
        base.groupby(["n", "family"])
        .agg(
            mean_pen_width_plus_consistency=("pen_width_plus_consistency", "mean"),
            mean_pen_residual_prior_gap=("pen_residual_prior_gap", "mean"),
            mean_pen_interval_bundle=("pen_interval_bundle", "mean"),
            mean_pen_noninterval_residual=("pen_noninterval_residual", "mean"),
            mean_pen_dim_proxy=("pen_dim_proxy", "mean"),
            mean_pen_dim_consistency_scaled=("pen_dim_consistency_scaled", "mean"),
            mean_pen_comparability=("pen_comparability", "mean"),
            mean_pen_cover=("pen_cover", "mean"),
            mean_pen_layer_smoothness=("pen_layer_smoothness", "mean"),
            count=("sample_id", "count"),
        )
        .reset_index()
    )

    contrast_rows = []
    for n, sub in family_summary.groupby("n", sort=True):
        kr = sub[sub["family"] == "KR_like"].iloc[0]
        lor = sub[sub["family"] == "lorentzian_like_2d"].iloc[0]
        contrast_rows.append(
            {
                "n": int(n),
                "delta_residual_prior_gap_lor_minus_kr": float(
                    lor["mean_pen_residual_prior_gap"] - kr["mean_pen_residual_prior_gap"]
                ),
                "delta_interval_bundle_lor_minus_kr": float(
                    lor["mean_pen_interval_bundle"] - kr["mean_pen_interval_bundle"]
                ),
                "delta_noninterval_residual_lor_minus_kr": float(
                    lor["mean_pen_noninterval_residual"] - kr["mean_pen_noninterval_residual"]
                ),
                "delta_dim_proxy_lor_minus_kr": float(lor["mean_pen_dim_proxy"] - kr["mean_pen_dim_proxy"]),
                "delta_dim_consistency_scaled_lor_minus_kr": float(
                    lor["mean_pen_dim_consistency_scaled"] - kr["mean_pen_dim_consistency_scaled"]
                ),
                "delta_comparability_lor_minus_kr": float(lor["mean_pen_comparability"] - kr["mean_pen_comparability"]),
                "delta_cover_lor_minus_kr": float(lor["mean_pen_cover"] - kr["mean_pen_cover"]),
                "delta_layer_smoothness_lor_minus_kr": float(
                    lor["mean_pen_layer_smoothness"] - kr["mean_pen_layer_smoothness"]
                ),
            }
        )
    contrast_df = pd.DataFrame(contrast_rows)

    family_summary.to_csv(out_dir / "residual_prior_gap_family_summary.csv", index=False, encoding="utf-8-sig")
    contrast_df.to_csv(out_dir / "residual_prior_gap_contrast.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "residual_prior_gap_family_summary.csv").as_posix())
    print((out_dir / "residual_prior_gap_contrast.csv").as_posix())
    print()
    print(contrast_df.to_string(index=False))
