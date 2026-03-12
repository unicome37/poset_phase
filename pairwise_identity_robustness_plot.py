from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"D:\Kiro\理论体系\poset_phase\outputs_exploratory\pairwise_blind_identity_scan")


LABELS = {
    "switch_zscore": "switch_zscore",
    "blind_knn3_zscore": "blind_knn3",
    "blind_knn5_zscore": "blind_knn5",
    "blind_self_hit_zscore": "blind_self",
    "blind_self_hit_centered": "blind_self_centered",
}

COLORS = {
    "switch_zscore": "#1d3557",
    "blind_knn3_zscore": "#2a9d8f",
    "blind_knn5_zscore": "#e9c46a",
    "blind_self_hit_zscore": "#e76f51",
    "blind_self_hit_centered": "#8d99ae",
}


def main() -> None:
    cross_df = pd.read_csv(ROOT / "pairwise_blind_identity_crossings.csv")
    cross_df["plot_label"] = cross_df["variant"].map(LABELS).fillna(cross_df["variant"])
    out_path = ROOT / "pairwise_identity_robustness_comparison.png"

    n_values = sorted(cross_df["n"].unique())
    variants = [
        "switch_zscore",
        "blind_knn3_zscore",
        "blind_knn5_zscore",
        "blind_self_hit_zscore",
        "blind_self_hit_centered",
    ]
    x = np.arange(len(n_values))
    width = 0.16

    fig, ax = plt.subplots(figsize=(10, 5.6))
    for i, variant in enumerate(variants):
        sub = cross_df[cross_df["variant"] == variant].set_index("n").reindex(n_values)
        vals = sub["zeta_cross"].to_numpy(dtype=float)
        xpos = x + (i - (len(variants) - 1) / 2) * width

        heights = np.nan_to_num(vals, nan=0.0)
        bars = ax.bar(
            xpos,
            heights,
            width=width,
            color=COLORS[variant],
            alpha=0.88,
            label=LABELS[variant],
            edgecolor="black",
            linewidth=0.4,
        )
        for j, (bar, val) in enumerate(zip(bars, vals)):
            if np.isnan(val):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    0.08,
                    "no cross",
                    ha="center",
                    va="bottom",
                    rotation=90,
                    fontsize=7,
                    color=COLORS[variant],
                )
            else:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.06,
                    f"{val:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    rotation=90,
                )

    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in n_values])
    ax.set_xlabel("N")
    ax.set_ylabel("zeta_cross")
    ax.set_title("Identity Robustness Comparison: switch_zscore vs Blind Alternatives")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(frameon=False, ncol=3, fontsize=8)
    ax.set_ylim(0, max(6.5, np.nanmax(cross_df["zeta_cross"].to_numpy(dtype=float)) + 0.6))

    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(out_path.as_posix())


if __name__ == "__main__":
    main()
