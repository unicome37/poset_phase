from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from experiment import FAMILIES, load_config
from observables import comparable_fraction, neutral_penalty
from observables_geo import (
    cover_density,
    geometric_components,
    interval_empty_fraction,
    mean_interval_size_ratio,
)


def run_diagnostics(
    n_values: tuple[int, ...],
    families: tuple[str, ...],
    samples_per_family: int,
) -> pd.DataFrame:
    rows = []
    for n in n_values:
        for family in families:
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 9000 * n + sample_id
                poset = generator(n=n, seed=seed)
                geo = geometric_components(poset)
                rows.append(
                    {
                        "n": n,
                        "family": family,
                        "sample_id": sample_id,
                        "comparable_fraction": comparable_fraction(poset),
                        "cover_density": cover_density(poset),
                        "interval_empty_fraction": interval_empty_fraction(poset),
                        "mean_interval_size_ratio": mean_interval_size_ratio(poset),
                        "neutral_penalty": neutral_penalty(poset),
                        **geo,
                    }
                )
    return pd.DataFrame(rows)


def summarize_diagnostics(df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = [c for c in df.columns if c not in {"n", "family", "sample_id"}]
    agg = {col: "mean" for col in metric_cols}
    agg["sample_id"] = "count"
    out = (
        df.groupby(["n", "family"])
        .agg(agg)
        .rename(columns={"sample_id": "count"})
        .reset_index()
        .sort_values(["n", "family"])
    )
    return out


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Inspect geometric diagnostic components.")
    parser.add_argument("--config", default="config_medium.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    exp_cfg = config["experiment"]
    out_cfg = config["output"]
    diag_families = tuple(
        config.get(
            "diagnostics",
            {},
        ).get(
            "families",
            ["multi_layer_random", "lorentzian_like_2d", "lorentzian_like_3d", "transitive_percolation"],
        )
    )
    diag_n_values = tuple(config.get("diagnostics", {}).get("n_values", exp_cfg["n_values"]))
    samples_per_family = int(config.get("diagnostics", {}).get("samples_per_family", exp_cfg["samples_per_family"]))

    raw = run_diagnostics(
        n_values=diag_n_values,
        families=diag_families,
        samples_per_family=samples_per_family,
    )
    summary = summarize_diagnostics(raw)

    out_dir = Path(out_cfg["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = out_dir / "diagnostics_raw.csv"
    summary_path = out_dir / "diagnostics_summary.csv"
    raw.to_csv(raw_path, index=False, encoding="utf-8-sig")
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

    print(summary.to_string(index=False))
