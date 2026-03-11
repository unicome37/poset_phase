from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from experiment import load_config


def plot_phase_curves(summary_df: pd.DataFrame, out_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for action_mode, mode_df in summary_df.groupby("action_mode"):
        fig, axes = plt.subplots(
            nrows=len(sorted(mode_df["n"].unique())),
            ncols=1,
            figsize=(10, 3.5 * len(sorted(mode_df["n"].unique()))),
            squeeze=False,
        )
        for ax, (n, sub) in zip(axes[:, 0], mode_df.groupby("n")):
            pivot = sub.pivot(index="gamma", columns="family", values="mean_score_norm").sort_index()
            pivot.plot(ax=ax, marker="o")
            ax.set_title(f"{action_mode} normalized phase curves, N={n}")
            ax.set_ylabel("mean_score_norm")
            ax.grid(True, alpha=0.3)
        axes[-1, 0].set_xlabel("gamma")
        fig.tight_layout()
        path = out_dir / f"phase_curves_{action_mode}.png"
        fig.savefig(path, dpi=160, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths


def plot_bootstrap_bars(bootstrap_df: pd.DataFrame, out_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for action_mode, mode_df in bootstrap_df.groupby("action_mode"):
        for n, sub_n in mode_df.groupby("n"):
            gammas = sorted(sub_n["gamma"].unique())
            fig, axes = plt.subplots(
                nrows=len(gammas),
                ncols=1,
                figsize=(10, 3.5 * len(gammas)),
                squeeze=False,
            )
            for ax, gamma in zip(axes[:, 0], gammas):
                sub = sub_n[sub_n["gamma"] == gamma].sort_values("score_norm_mean")
                y = range(len(sub))
                means = sub["score_norm_mean"].to_numpy()
                low = means - sub["score_norm_ci_low"].to_numpy()
                high = sub["score_norm_ci_high"].to_numpy() - means
                ax.errorbar(means, y, xerr=[low, high], fmt="o")
                ax.set_yticks(list(y))
                ax.set_yticklabels(sub["family"])
                ax.set_title(f"{action_mode}, N={n}, gamma={gamma}")
                ax.set_xlabel("score_norm_mean with 95% CI")
                ax.grid(True, axis="x", alpha=0.3)
            fig.tight_layout()
            path = out_dir / f"bootstrap_ci_{action_mode}_N{n}.png"
            fig.savefig(path, dpi=160, bbox_inches="tight")
            plt.close(fig)
            paths.append(path)
    return paths


def plot_calibration(calib_summary_df: pd.DataFrame, out_dir: Path) -> list[Path]:
    paths: list[Path] = []
    fig, ax = plt.subplots(figsize=(10, 5))
    for family, sub in calib_summary_df.groupby("family"):
        sub = sub.sort_values("n")
        ax.plot(sub["n"], sub["mean_abs_error"], marker="o", label=family)
    ax.set_title("SIS calibration: mean absolute error vs N")
    ax.set_xlabel("N")
    ax.set_ylabel("mean_abs_error")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    path = out_dir / "calibration_mean_abs_error.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)
    return paths


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate plots from experiment outputs.")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_cfg = config["output"]
    plot_cfg = config["plots"]

    out_dir = Path(out_cfg["directory"])
    plot_dir = out_dir / plot_cfg["directory"]
    plot_dir.mkdir(parents=True, exist_ok=True)

    summary_df = pd.read_csv(out_dir / out_cfg["summary_csv"])
    bootstrap_df = pd.read_csv(out_dir / out_cfg["bootstrap_csv"])
    calib_summary_df = pd.read_csv(out_dir / out_cfg["calibration_summary_csv"])

    created = []
    created.extend(plot_phase_curves(summary_df, plot_dir))
    created.extend(plot_bootstrap_bars(bootstrap_df, plot_dir))
    created.extend(plot_calibration(calib_summary_df, plot_dir))

    for path in created:
        print(path)
