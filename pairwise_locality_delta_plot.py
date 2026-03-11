from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot pairwise locality-dominance deltas against log_H deltas.")
    parser.add_argument("--config", default="config_pairwise_locality_delta_validation.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    raw_path = out_dir / "pairwise_locality_delta_validation_raw.csv"
    summary_path = out_dir / "pairwise_locality_delta_validation_summary.csv"

    raw_df = pd.read_csv(raw_path)
    summary_df = pd.read_csv(summary_path)
    overall = summary_df[summary_df["scope"] == "all_pairs"].iloc[0]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    cmap = plt.get_cmap("viridis")
    ns = sorted(raw_df["n"].unique().tolist())
    color_map = {n: cmap(i / max(1, len(ns) - 1)) for i, n in enumerate(ns)}

    for n in ns:
        sub = raw_df[raw_df["n"] == n]
        ax.scatter(
            sub["locality_dominance_delta_index"],
            sub["log_H_delta_mlr_minus_lor2d"],
            s=46,
            alpha=0.85,
            color=color_map[n],
            label=f"N={n}",
        )

    x = raw_df["locality_dominance_delta_index"].to_numpy(dtype=float)
    y = raw_df["log_H_delta_mlr_minus_lor2d"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    x_line = np.linspace(float(x.min()), float(x.max()), 200)
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, color="#d62728", linewidth=2.0, label="Linear fit")

    ax.set_xlabel("Locality-Dominance Delta Index")
    ax.set_ylabel("log_H Delta (MLR - Lor2D)")
    ax.set_title(
        "Matched Pair Locality Dominance vs log_H Delta\n"
        f"Pearson={overall['pearson_corr']:.3f}, Spearman={overall['spearman_corr']:.3f}, p={overall['permutation_pvalue']:.4f}"
    )
    ax.grid(alpha=0.25, linestyle="--")
    ax.legend(frameon=False, ncol=3)

    fig.tight_layout()
    out_path = out_dir / "pairwise_locality_delta_scatter.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    print(out_path.as_posix())
