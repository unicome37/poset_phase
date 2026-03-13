from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from experiment import FAMILIES
from matched_residual_freedom_check import residual_metrics
from observables import antichain_width, comparable_fraction
from observables_geo import geometric_components


REQUIRED_COLUMNS = ("family", "n", "seed")


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Augment raw CSVs with Prediction C hierarchy metrics reconstructed from family/n/seed."
    )
    parser.add_argument(
        "--config",
        default="config_augment_prediction_c_features.yaml",
        help="Path to YAML config file.",
    )
    return parser


def normalize_source_path(root: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else root / path


def load_sources(root: Path, source_paths: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    summary_rows: list[dict] = []

    for raw_path in source_paths:
        csv_path = normalize_source_path(root, raw_path)
        df = pd.read_csv(csv_path)
        missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
        if missing:
            raise ValueError(f"{csv_path} missing required columns: {missing}")
        df = df.copy()
        df["source_path"] = csv_path.as_posix()
        frames.append(df)
        summary_rows.append(
            {
                "source_path": csv_path.as_posix(),
                "rows": int(len(df)),
                "unique_posets": int(df.loc[:, ["family", "n", "seed"]].drop_duplicates().shape[0]),
            }
        )

    if not frames:
        return pd.DataFrame(), pd.DataFrame(summary_rows)

    combined = pd.concat(frames, ignore_index=True)
    return combined, pd.DataFrame(summary_rows)


def compute_metric_cache(unique_posets: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []

    for rec in unique_posets.itertuples(index=False):
        family = str(rec.family)
        n = int(rec.n)
        seed = int(rec.seed)
        if family not in FAMILIES:
            raise KeyError(f"Unsupported family: {family}")
        poset = FAMILIES[family](n=n, seed=seed)
        metrics = residual_metrics(poset)
        geo = geometric_components(poset)
        row = {"family": family, "n": n, "seed": seed}
        row["antichain_width"] = float(antichain_width(poset))
        row["comparable_fraction"] = float(comparable_fraction(poset))
        row["geo_dim_eff"] = float(geo["geo_dim_eff"])
        row["geo_interval_shape"] = float(geo["geo_interval_shape"])
        row.update(metrics)
        rows.append(row)

    return pd.DataFrame(rows)


def write_outputs(
    output_dir: Path,
    combined_df: pd.DataFrame,
    source_summary_df: pd.DataFrame,
    cache_df: pd.DataFrame,
    augmented_df: pd.DataFrame,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    source_summary_df.to_csv(
        output_dir / "prediction_c_source_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    cache_df.to_csv(
        output_dir / "prediction_c_feature_cache.csv",
        index=False,
        encoding="utf-8-sig",
    )
    augmented_df.to_csv(
        output_dir / "prediction_c_augmented_combined.csv",
        index=False,
        encoding="utf-8-sig",
    )

    unique_augmented = augmented_df.loc[:, ["source_path", "family", "n", "seed"]].drop_duplicates()
    source_augmented_summary = (
        unique_augmented.groupby("source_path")
        .agg(unique_posets_augmented=("seed", "count"))
        .reset_index()
        .merge(source_summary_df, on="source_path", how="left")
    )
    source_augmented_summary.to_csv(
        output_dir / "prediction_c_augmented_source_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    family_summary = (
        cache_df.groupby("family")
        .agg(
            unique_posets=("seed", "count"),
            n_min=("n", "min"),
            n_max=("n", "max"),
            mean_layer_count=("layer_count", "mean"),
            mean_layer_gap=("mean_layer_gap", "mean"),
        )
        .reset_index()
    )
    family_summary.to_csv(
        output_dir / "prediction_c_augmented_family_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    n_summary = (
        cache_df.groupby("n")
        .agg(
            unique_posets=("seed", "count"),
            mean_layer_count=("layer_count", "mean"),
            mean_layer_gap=("mean_layer_gap", "mean"),
        )
        .reset_index()
        .sort_values("n")
    )
    n_summary.to_csv(
        output_dir / "prediction_c_augmented_n_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print((output_dir / "prediction_c_feature_cache.csv").as_posix())
    print((output_dir / "prediction_c_augmented_combined.csv").as_posix())
    print()
    print(source_augmented_summary.to_string(index=False))
    print()
    print(family_summary.to_string(index=False))


def main() -> None:
    args = build_arg_parser().parse_args()
    config_path = Path(args.config)
    config = load_config(config_path)

    root = config_path.parent
    source_paths = list(config["inputs"]["csvs"])
    output_dir = normalize_source_path(root, str(config["output"]["directory"]))

    combined_df, source_summary_df = load_sources(root, source_paths)
    if combined_df.empty:
        raise ValueError("No input CSVs loaded.")

    unique_posets = combined_df.loc[:, ["family", "n", "seed"]].drop_duplicates()
    cache_df = compute_metric_cache(unique_posets)
    augmented_df = combined_df.merge(cache_df, on=["family", "n", "seed"], how="left")

    write_outputs(output_dir, combined_df, source_summary_df, cache_df, augmented_df)


if __name__ == "__main__":
    main()
