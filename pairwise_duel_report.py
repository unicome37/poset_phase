from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml


METRICS = [
    ("mean_log_H", "log H"),
    ("mean_geo_total", "geo_total"),
    ("mean_cg_penalty", "cg_penalty"),
    ("mean_cg_switch_rate", "cg_switch_rate"),
]


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_comparison_table(rank_df: pd.DataFrame) -> pd.DataFrame:
    pivot = rank_df.pivot(index="n", columns="family")
    rows: list[dict] = []
    for n in sorted(rank_df["n"].unique()):
        row = {"n": int(n)}
        for metric, _label in METRICS + [("mean_score_A2", "score_A2")]:
            lor = float(pivot.loc[n, (metric, "lorentzian_like_2d")])
            mlr = float(pivot.loc[n, (metric, "multi_layer_random")])
            row[f"{metric}_lor2d"] = lor
            row[f"{metric}_mlr"] = mlr
            row[f"{metric}_delta_mlr_minus_lor2d"] = mlr - lor
        rows.append(row)
    return pd.DataFrame(rows)


def plot_duel_metrics(rank_df: pd.DataFrame, out_path: Path) -> None:
    families = ["lorentzian_like_2d", "multi_layer_random"]
    colors = {
        "lorentzian_like_2d": "#1f77b4",
        "multi_layer_random": "#d62728",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    axes = axes.ravel()

    for ax, (metric, label) in zip(axes, METRICS):
        for family in families:
            sub = rank_df[rank_df["family"] == family].sort_values("n")
            ax.plot(sub["n"], sub[metric], marker="o", label=family, color=colors[family])
        ax.set_title(label)
        ax.set_xlabel("N")
        ax.grid(True, alpha=0.3)
    axes[0].legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def write_markdown_summary(comparison_df: pd.DataFrame, out_path: Path) -> None:
    lines = [
        "# Pairwise Duel Summary",
        "",
        "Lor2D vs multi_layer_random under matched width+comparability windows.",
        "",
        "| N | score_A2 delta (MLR-Lor2D) | logH delta | geo_total delta | cg_penalty delta | switch_rate delta |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in comparison_df.itertuples(index=False):
        lines.append(
            "| "
            f"{row.n} | "
            f"{row.mean_score_A2_delta_mlr_minus_lor2d:.3f} | "
            f"{row.mean_log_H_delta_mlr_minus_lor2d:.3f} | "
            f"{row.mean_geo_total_delta_mlr_minus_lor2d:.3f} | "
            f"{row.mean_cg_penalty_delta_mlr_minus_lor2d:.3f} | "
            f"{row.mean_cg_switch_rate_delta_mlr_minus_lor2d:.3f} |"
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "- Negative `score_A2 delta` means `multi_layer_random` scores lower under the current A2.",
            "- Positive `geo_total` / `cg_penalty` / `switch_rate` deltas mean `multi_layer_random` is less geometric and less coarse-grain stable than Lor2D.",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate pairwise duel comparison table and plot.")
    parser.add_argument(
        "--config",
        default="config_pairwise_compressibility_duel.yaml",
        help="Path to YAML config file.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    rank_df = pd.read_csv(out_dir / "pairwise_duel_rankings.csv")

    comparison_df = build_comparison_table(rank_df)
    comparison_df.to_csv(out_dir / "pairwise_duel_comparison.csv", index=False, encoding="utf-8-sig")
    plot_duel_metrics(rank_df, out_dir / "pairwise_duel_metrics.png")
    write_markdown_summary(comparison_df, out_dir / "pairwise_duel_summary.md")

    print((out_dir / "pairwise_duel_comparison.csv").as_posix())
    print((out_dir / "pairwise_duel_metrics.png").as_posix())
    print((out_dir / "pairwise_duel_summary.md").as_posix())
