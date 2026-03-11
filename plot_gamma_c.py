from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from experiment import load_config


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot gamma_c(N) from gamma_c_report.csv.")
    parser.add_argument("--config", default="config_smallN_exact.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    plot_dir = out_dir / config["plots"]["directory"]
    plot_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(out_dir / "gamma_c_report.csv")
    df = df[df["action_mode"] == "A2"].copy()

    fig, ax = plt.subplots(figsize=(8, 5))
    for family_b, sub in df.groupby("family_b"):
        sub = sub.sort_values("n")
        ax.plot(sub["n"], sub["gamma_c_est"], marker="o", label=f"vs {family_b}")

    ax.set_title(r"Estimated $\gamma_c(N)$ for lorentzian_like_2d under $A_2$")
    ax.set_xlabel("N")
    ax.set_ylabel(r"estimated $\gamma_c$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()

    path = plot_dir / "gamma_c_A2.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(path)
