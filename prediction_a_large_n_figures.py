"""
Prediction A — Large-N Closure Figures
=======================================

Generates publication figures for the large-N finite-size scaling results:
  1. Lambda-N heatmap: 4D-winner regions across (N, λ) grid
  2. Xi stability across N (extending beyond N=68)
  3. Link density crossover: 3D vs 4D convergence at large N
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = pathlib.Path("outputs_exploratory/prediction_a_large_n_scaling")
FIG_DIR = OUT_DIR
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ===== Load data =====
raw = pd.read_csv(OUT_DIR / "raw_observables_large_n.csv")
winners = pd.read_csv(OUT_DIR / "winners_large_n.csv")
xi = pd.read_csv(OUT_DIR / "xi_large_n.csv")

agg = raw.groupby(["n", "dim"]).agg(
    link_density=("link_density", "mean"),
    S_link_norm=("S_link_d2_norm", "mean"),
    log_H=("log_H", "mean"),
).reset_index()


# ===== Figure 1: Lambda-N Winner Heatmap =====
def fig_winner_heatmap():
    fig, ax = plt.subplots(figsize=(8, 4))

    n_vals = sorted(winners["n"].unique())
    lam_vals = sorted(winners["lambda"].unique())

    dim_to_num = {"2d": 2, "3d": 3, "4d": 4, "5d": 5}
    dim_colors = {2: "#3498db", 3: "#e67e22", 4: "#e74c3c", 5: "#9b59b6"}
    dim_labels = {2: "2D", 3: "3D", 4: "4D", 5: "5D"}

    for i, lam in enumerate(lam_vals):
        for j, n in enumerate(n_vals):
            row = winners[(winners["lambda"] == lam) & (winners["n"] == n)]
            if row.empty:
                continue
            w = dim_to_num[row.iloc[0]["winner"]]
            rect = plt.Rectangle((j - 0.45, i - 0.45), 0.9, 0.9,
                                  facecolor=dim_colors[w], alpha=0.85,
                                  edgecolor="white", linewidth=1.5)
            ax.add_patch(rect)
            ax.text(j, i, f"{dim_labels[w]}", ha="center", va="center",
                    fontsize=10, fontweight="bold", color="white")

    # Mark the 4D region boundary
    ax.set_xlim(-0.5, len(n_vals) - 0.5)
    ax.set_ylim(-0.5, len(lam_vals) - 0.5)
    ax.set_xticks(range(len(n_vals)))
    ax.set_xticklabels([str(int(n)) for n in n_vals], fontsize=10)
    ax.set_yticks(range(len(lam_vals)))
    ax.set_yticklabels([f"λ={int(l)}" for l in lam_vals], fontsize=10)
    ax.set_xlabel("N (poset size)", fontsize=12)
    ax.set_ylabel("Coupling λ", fontsize=12)
    ax.set_title("4D Selection Window Shifts Right with N", fontsize=13, fontweight="bold")

    # Legend
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=dim_colors[d], alpha=0.85)
               for d in [2, 3, 4, 5]]
    ax.legend(handles, [dim_labels[d] for d in [2, 3, 4, 5]],
              loc="upper left", fontsize=9, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "large_n_winner_heatmap.png", dpi=200)
    fig.savefig(FIG_DIR / "large_n_winner_heatmap.pdf")
    plt.close(fig)
    print("  Saved large_n_winner_heatmap.png/pdf")


# ===== Figure 2: Xi stability across N =====
def fig_xi_stability():
    xi45 = xi[xi["dim_pair"] == "4→5"].sort_values("n")

    fig, ax = plt.subplots(figsize=(7, 4))
    ns = xi45["n"].values
    xis = xi45["Xi"].values

    # Separate small N (≤68) and large N (>68)
    mask_small = ns <= 68
    mask_large = ns > 68

    ax.scatter(ns[mask_small], xis[mask_small], c="#2196F3", s=80, zorder=5,
               label="N ≤ 68 (original)", edgecolors="white", linewidth=1.2)
    ax.scatter(ns[mask_large], xis[mask_large], c="#FF5722", s=100, zorder=5,
               label="N > 68 (new)", marker="D", edgecolors="white", linewidth=1.2)

    # Median band
    med = np.median(xis)
    std_val = np.std(xis)
    ax.axhline(med, color="gray", linestyle="--", alpha=0.5, label=f"Median = {med:.1f}")
    ax.axhspan(med - std_val, med + std_val, alpha=0.1, color="gray")

    # Annotations
    for n, x in zip(ns, xis):
        ax.annotate(f"{x:.1f}", (n, x), textcoords="offset points",
                    xytext=(0, 10), fontsize=8, ha="center", color="#555")

    ax.set_xlabel("N (poset size)", fontsize=12)
    ax.set_ylabel("Ξ₄→₅", fontsize=14)
    ax.set_title("Ξ₄→₅ Stability: N = 20 – 112", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10, framealpha=0.9)
    ax.set_ylim(5, 18)
    ax.set_xlim(10, 120)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "xi_stability_large_n.png", dpi=200)
    fig.savefig(FIG_DIR / "xi_stability_large_n.pdf")
    plt.close(fig)
    print("  Saved xi_stability_large_n.png/pdf")


# ===== Figure 3: Link density crossover =====
def fig_link_density_crossover():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    dim_colors = {"2d": "#3498db", "3d": "#e67e22", "4d": "#e74c3c", "5d": "#9b59b6"}
    dim_labels = {"2d": "2D", "3d": "3D", "4d": "4D", "5d": "5D"}

    # Panel A: Link density vs N for all dims
    for dim in ["2d", "3d", "4d", "5d"]:
        sub = agg[agg["dim"] == dim].sort_values("n")
        ax1.plot(sub["n"], sub["link_density"], "o-", color=dim_colors[dim],
                 label=dim_labels[dim], markersize=6, linewidth=1.5)

    ax1.set_xlabel("N", fontsize=12)
    ax1.set_ylabel("Link density (C₀ / C(N,2))", fontsize=11)
    ax1.set_title("(a) Link Density Profiles", fontsize=12, fontweight="bold")
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Panel B: 3D–4D gap and 4D–5D gap
    ns = sorted(agg["n"].unique())
    gap_34 = []
    gap_45 = []
    for n in ns:
        sub = agg[agg["n"] == n].set_index("dim")
        gap_34.append(sub.loc["3d", "link_density"] - sub.loc["4d", "link_density"])
        gap_45.append(sub.loc["4d", "link_density"] - sub.loc["5d", "link_density"])

    ax2.plot(ns, gap_34, "s-", color="#e67e22", label="Δ(3D – 4D)", markersize=7, linewidth=1.5)
    ax2.plot(ns, gap_45, "D-", color="#9b59b6", label="Δ(4D – 5D)", markersize=7, linewidth=1.5)
    ax2.axhline(0, color="gray", linestyle=":", alpha=0.5)

    # Mark the crossover
    crossover_n = None
    for i, (n, g34) in enumerate(zip(ns, gap_34)):
        if g34 < 0 and crossover_n is None:
            crossover_n = n
    if crossover_n:
        ax2.axvline(crossover_n, color="red", linestyle="--", alpha=0.4, label=f"3D–4D crossover ~ N={crossover_n}")

    ax2.set_xlabel("N", fontsize=12)
    ax2.set_ylabel("Δ(link density)", fontsize=11)
    ax2.set_title("(b) Inter-Dimension Gap Asymmetry", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=9, loc="upper right")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "link_density_crossover.png", dpi=200)
    fig.savefig(FIG_DIR / "link_density_crossover.pdf")
    plt.close(fig)
    print("  Saved link_density_crossover.png/pdf")


# ===== Main =====
if __name__ == "__main__":
    print("Generating large-N closure figures...")
    fig_winner_heatmap()
    fig_xi_stability()
    fig_link_density_crossover()
    print("\nDone — all figures saved to", FIG_DIR)
