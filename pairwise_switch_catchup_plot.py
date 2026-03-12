from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(r"D:\Kiro")


def main() -> None:
    in_dir = ROOT / "outputs_exploratory" / "pairwise_switch_catchup_analysis"
    out_path = in_dir / "pairwise_switch_catchup_plot.png"

    df = pd.read_csv(in_dir / "pairwise_switch_catchup_summary.csv").sort_values("n")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    ax = axes[0, 0]
    ax.plot(df["n"], df["zeta_cross"], marker="o", linewidth=2, label="zeta_cross")
    ax.set_title("Switch Weight Needed For Reversal")
    ax.set_xlabel("N")
    ax.set_ylabel("zeta_cross")
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(df["n"], df["mean_log_H_delta"], marker="o", linewidth=2, label="log_H_delta")
    ax.plot(df["n"], df["abs_score_gap"], marker="s", linewidth=2, label="|A2 gap|")
    ax.set_title("Growth Of Entropic Advantage")
    ax.set_xlabel("N")
    ax.set_ylabel("delta")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(df["n"], df["mean_switch_delta"], marker="o", linewidth=2, color="#2a9d8f", label="switch_delta")
    ax.set_title("Raw Switch Separation")
    ax.set_xlabel("N")
    ax.set_ylabel("delta switch")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(df["n"], df["zeta_per_logH"], marker="o", linewidth=2, label="zeta_cross / log_H_delta")
    ax.plot(df["n"], df["score_gap_per_switch"], marker="s", linewidth=2, label="|A2 gap| / switch_delta")
    ax.set_title("Catch-Up Ratios")
    ax.set_xlabel("N")
    ax.set_ylabel("ratio")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)

    fig.suptitle("Pairwise Catch-Up Analysis: Lor2D vs Quasi-Lorentzian MLR", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print(out_path.as_posix())


if __name__ == "__main__":
    main()
