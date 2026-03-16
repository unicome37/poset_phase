from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot Prediction B gamma_c sensitivity curves.")
    parser.add_argument(
        "--config",
        default="config_prediction_b_weight_sensitivity.yaml",
        help="Path to YAML config.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    report_path = out_dir / "prediction_b_weight_sensitivity_report.csv"
    report_df = pd.read_csv(report_path)
    plot_df = report_df[report_df["status"] == "crossing"].copy()

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), sharey=True)
    colors = {20: "#1f77b4", 28: "#ff7f0e", 36: "#2ca02c", 44: "#d62728"}

    for ax, weight_name, title in [
        (axes[0], "geo_dim_proxy_penalty", "gamma_c vs dim-proxy weight"),
        (axes[1], "geo_width_height", "gamma_c vs width-height weight"),
    ]:
        sub = plot_df[plot_df["weight_name"] == weight_name].copy()
        for n, group in sub.groupby("n", sort=True):
            group = group.sort_values("weight_value")
            ax.plot(
                group["weight_value"],
                group["gamma_c_est"],
                marker="o",
                linewidth=1.8,
                markersize=4.5,
                color=colors.get(int(n), None),
                label=f"N={int(n)}",
            )
        ax.set_title(title)
        ax.set_xlabel("Weight value")
        ax.grid(alpha=0.25, linestyle="--")

    axes[0].set_ylabel("Estimated gamma_c")
    axes[1].legend(frameon=False, loc="best")
    fig.suptitle("Prediction B: gamma_c shifts continuously with backbone weights", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_dir / "prediction_b_weight_sensitivity.png", dpi=170, bbox_inches="tight")

    print((out_dir / "prediction_b_weight_sensitivity.png").as_posix())
