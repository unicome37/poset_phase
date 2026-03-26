"""Generate figures highlighting KR_2layer / KR_4layer control groups.

Fig A: Score ranking heatmap — all 17 families × 5 gammas (A2, N=80)
Fig B: Control group separation — score_norm trajectories across gamma
Fig C: logH vs geometric_penalty scatter — color-coded by family category
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.size": 10,
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.fontsize": 8,
    "legend.framealpha": 0.9,
    "legend.edgecolor": "0.7",
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

OUT_DIR = Path("manuscript_figures")
OUT_DIR.mkdir(exist_ok=True)

# Color palette
BLUE = "#2563eb"
RED = "#dc2626"
GREEN = "#16a34a"
PURPLE = "#7c3aed"
AMBER = "#d97706"
GRAY = "#6b7280"
CYAN = "#0891b2"
PINK = "#db2777"
TEAL = "#0d9488"
ORANGE = "#ea580c"

# Family categories for coloring
LORENTZIAN = ["lorentzian_like_2d", "lorentzian_like_3d",
              "lorentzian_like_4d", "lorentzian_like_5d"]
CONTROLS = ["KR_2layer", "KR_4layer"]
KR_FAMILY = ["KR_like"]
LAYERED = ["absolute_layered", "multi_layer_random",
           "random_layered_k4_uniform", "random_layered_k6_uniform",
           "random_layered_k8_uniform", "random_layered_k6_tapered",
           "random_layered_k6_middle_heavy", "random_layered_k6_longjump"]
OTHER = ["transitive_percolation", "interval_order"]

SHORT_NAMES = {
    "absolute_layered": "AbsLayer",
    "KR_like": "KR (3-layer)",
    "KR_2layer": "KR (2-layer) [ctrl]",
    "KR_4layer": "KR (4-layer) [ctrl]",
    "lorentzian_like_2d": "Lor2D",
    "lorentzian_like_3d": "Lor3D",
    "lorentzian_like_4d": "Lor4D",
    "lorentzian_like_5d": "Lor5D",
    "transitive_percolation": "TransPerc",
    "interval_order": "IntOrder",
    "multi_layer_random": "MLR",
    "random_layered_k4_uniform": "RLk4",
    "random_layered_k6_uniform": "RLk6",
    "random_layered_k8_uniform": "RLk8",
    "random_layered_k6_tapered": "RLk6-tap",
    "random_layered_k6_middle_heavy": "RLk6-mid",
    "random_layered_k6_longjump": "RLk6-lj",
}


def get_category_color(family: str) -> str:
    if family in LORENTZIAN:
        return BLUE
    if family in CONTROLS:
        return RED
    if family in KR_FAMILY:
        return ORANGE
    if family in OTHER:
        return PURPLE
    return GRAY  # layered variants


def load_data() -> pd.DataFrame:
    return pd.read_csv("outputs/raw_samples.csv")


# ── Fig A: Score ranking heatmap ──────────────────────────────────

def fig_a_ranking_heatmap(df: pd.DataFrame):
    """Heatmap: mean score_norm rank per family across gamma (A2, N=80)."""
    sub = df[(df["action_mode"] == "A2") & (df["n"] == 80)]
    agg = sub.groupby(["gamma", "family"])["score_norm"].mean().reset_index()

    gammas = sorted(agg["gamma"].unique())
    families_at_g0 = (
        agg[agg["gamma"] == 0.0]
        .sort_values("score_norm", ascending=False)["family"]
        .tolist()
    )

    # Build rank matrix
    rank_matrix = np.zeros((len(families_at_g0), len(gammas)))
    for j, g in enumerate(gammas):
        sub_g = agg[agg["gamma"] == g].sort_values("score_norm", ascending=False)
        rank_map = {f: r + 1 for r, f in enumerate(sub_g["family"])}
        for i, fam in enumerate(families_at_g0):
            rank_matrix[i, j] = rank_map.get(fam, len(families_at_g0))

    fig, ax = plt.subplots(figsize=(7, 6))
    cmap = plt.cm.RdYlGn_r
    im = ax.imshow(rank_matrix, aspect="auto", cmap=cmap, vmin=1,
                   vmax=len(families_at_g0), interpolation="nearest")

    # Annotate cells
    for i in range(len(families_at_g0)):
        for j in range(len(gammas)):
            r = int(rank_matrix[i, j])
            color = "white" if r > 12 or r < 4 else "black"
            ax.text(j, i, str(r), ha="center", va="center",
                    fontsize=8, color=color, fontweight="bold")

    # Axis labels
    short = [SHORT_NAMES.get(f, f) for f in families_at_g0]
    ax.set_yticks(range(len(families_at_g0)))
    ax.set_yticklabels(short, fontsize=8)
    ax.set_xticks(range(len(gammas)))
    ax.set_xticklabels([f"$\\gamma={g}$" for g in gammas])
    ax.set_xlabel(r"Geometric coupling $\gamma$")
    ax.set_title(r"Family rank by $\langle S_{\rm norm}\rangle$ (A2, $N=80$)")

    # Highlight control rows
    for i, fam in enumerate(families_at_g0):
        if fam in CONTROLS:
            rect = plt.Rectangle((-0.5, i - 0.5), len(gammas), 1,
                                 linewidth=2, edgecolor=RED, facecolor="none",
                                 linestyle="--", zorder=4)
            ax.add_patch(rect)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("Rank (1 = highest score)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fig.savefig(OUT_DIR / "fig_control_ranking_heatmap.png")
    fig.savefig(OUT_DIR / "fig_control_ranking_heatmap.pdf")
    plt.close(fig)
    print("Fig A (ranking heatmap) saved.")


# ── Fig B: Control group separation trajectories ─────────────────

def fig_b_separation_trajectories(df: pd.DataFrame):
    """Score_norm trajectories across gamma for key families (A2, N=60,80)."""
    highlight = [
        "lorentzian_like_2d", "lorentzian_like_4d",
        "KR_like", "KR_2layer", "KR_4layer",
        "transitive_percolation",
    ]
    colors = {
        "lorentzian_like_2d": BLUE,
        "lorentzian_like_4d": CYAN,
        "KR_like": ORANGE,
        "KR_2layer": RED,
        "KR_4layer": PINK,
        "transitive_percolation": GREEN,
    }
    markers = {
        "lorentzian_like_2d": "o",
        "lorentzian_like_4d": "D",
        "KR_like": "^",
        "KR_2layer": "s",
        "KR_4layer": "v",
        "transitive_percolation": "P",
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for ax, n_val in zip(axes, [60, 80]):
        sub = df[(df["action_mode"] == "A2") & (df["n"] == n_val)]
        agg = sub.groupby(["gamma", "family"])["score_norm"].agg(
            ["mean", "std"]).reset_index()

        # Background: other families as thin gray lines
        for fam in agg["family"].unique():
            if fam not in highlight:
                fam_data = agg[agg["family"] == fam].sort_values("gamma")
                ax.plot(fam_data["gamma"], fam_data["mean"],
                        color="#d1d5db", linewidth=0.6, alpha=0.5, zorder=1)

        # Highlighted families
        for fam in highlight:
            fam_data = agg[agg["family"] == fam].sort_values("gamma")
            ax.plot(fam_data["gamma"], fam_data["mean"],
                    f"{markers[fam]}-", color=colors[fam], linewidth=1.5,
                    markersize=6, label=SHORT_NAMES[fam],
                    markeredgecolor="white", markeredgewidth=0.4, zorder=3)

        ax.set_xlabel(r"$\gamma$")
        ax.set_title(f"$N = {n_val}$")
        ax.grid(True, alpha=0.2, linewidth=0.5)
        ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")

    axes[0].set_ylabel(r"Mean $S_{\rm norm}$ (A2)")
    axes[1].legend(loc="best", fontsize=8, ncol=1)

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.suptitle("Control group separation: score trajectories across "
                 r"$\gamma$", fontsize=13)
    fig.savefig(OUT_DIR / "fig_control_separation.png")
    fig.savefig(OUT_DIR / "fig_control_separation.pdf")
    plt.close(fig)
    print("Fig B (separation trajectories) saved.")


# ── Fig C: logH vs geometric penalty scatter ─────────────────────

def fig_c_entropy_vs_geometry(df: pd.DataFrame):
    """Scatter of mean logH vs mean geometric penalty per family (N=80)."""
    sub = df[(df["action_mode"] == "A2") & (df["n"] == 80) & (df["gamma"] == 0.0)]
    agg = (sub.groupby("family")
           .agg(mean_logH=("log_H_mean", "mean"),
                mean_geo=("penalty_geometric", "mean"),
                std_logH=("log_H_mean", "std"),
                std_geo=("penalty_geometric", "std"))
           .reset_index())

    fig, ax = plt.subplots(figsize=(7, 5.5))

    # Manual label offsets for dense clusters
    OFFSETS = {
        "lorentzian_like_2d": (-35, -12),
        "interval_order": (-40, -12),
        "random_layered_k8_uniform": (8, -12),
        "random_layered_k6_uniform": (-38, 6),
        "random_layered_k6_longjump": (8, -12),
        "random_layered_k6_tapered": (8, 8),
        "random_layered_k6_middle_heavy": (8, -12),
        "transitive_percolation": (-30, 8),
        "random_layered_k4_uniform": (-30, -12),
        "multi_layer_random": (-12, 10),
        "KR_4layer": (-50, 10),
        "lorentzian_like_3d": (8, 6),
        "absolute_layered": (-5, 8),
        "KR_like": (8, 6),
        "lorentzian_like_4d": (8, 4),
        "lorentzian_like_5d": (8, 4),
        "KR_2layer": (8, 4),
    }

    for _, row in agg.iterrows():
        fam = row["family"]
        c = get_category_color(fam)
        marker = "s" if fam in CONTROLS else ("^" if fam in KR_FAMILY else "o")
        size = 90 if fam in CONTROLS else 60
        edge = RED if fam in CONTROLS else "white"
        ax.scatter(row["mean_logH"], row["mean_geo"], c=c, s=size,
                   marker=marker, edgecolors=edge, linewidths=1.2, zorder=3)
        oxy = OFFSETS.get(fam, (8, 4))
        ax.annotate(SHORT_NAMES.get(fam, fam),
                    (row["mean_logH"], row["mean_geo"]),
                    textcoords="offset points", xytext=oxy,
                    fontsize=6.5, color=c, alpha=0.85)

    # Legend patches
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE,
               markersize=8, label="Lorentzian-like"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=RED,
               markersize=8, markeredgecolor=RED, label="KR controls (2/4-layer)"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor=ORANGE,
               markersize=8, label="KR (3-layer)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=PURPLE,
               markersize=8, label="TransPerc / IntOrder"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=GRAY,
               markersize=8, label="Layered variants"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=7.5)

    ax.set_xlabel(r"Mean $\log|L(P)|$")
    ax.set_ylabel(r"Mean geometric penalty $\Pi_{\rm geo}$")
    ax.set_title(r"Entropy--geometry plane ($N=80$, A2)")
    ax.grid(True, alpha=0.2, linewidth=0.5)

    fig.savefig(OUT_DIR / "fig_entropy_geometry_scatter.png")
    fig.savefig(OUT_DIR / "fig_entropy_geometry_scatter.pdf")
    plt.close(fig)
    print("Fig C (entropy-geometry scatter) saved.")


# ── Main ─────────────────────────────────────────────────────────

if __name__ == "__main__":
    df = load_data()
    fig_a_ranking_heatmap(df)
    fig_b_separation_trajectories(df)
    fig_c_entropy_vs_geometry(df)
    print("All control-group figures generated.")
