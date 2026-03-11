"""Generate publication-quality figures for the manuscript.

Fig 1: gamma_c(N) curve under A2 (confirmatory mainline) with inset
Fig 2: Ablation summary – heatmap matrix
Fig 3: Non-circular replacement comparison
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# ── Publication style ───────────────────────────────────────────────
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

BLUE = "#2563eb"
RED = "#dc2626"
GREEN = "#16a34a"
PURPLE = "#7c3aed"
AMBER = "#d97706"
GRAY = "#6b7280"


# ── Fig 1: gamma_c(N) with inset ──────────────────────────────────

def fig1_gamma_c_curve():
    frozen = pd.read_csv("outputs_confirmatory/frozen_exact/gamma_c_report.csv")
    frozen_a2 = frozen[(frozen["action_mode"] == "A2") & (frozen["family_b"] == "KR_like")]
    medium = pd.read_csv("outputs_confirmatory/medium_exact_scan/gamma_c_report.csv")
    medium_a2 = medium[(medium["action_mode"] == "A2") & (medium["family_b"] == "KR_like")]

    ns = np.concatenate([frozen_a2["n"].values, medium_a2["n"].values])
    gc = np.concatenate([frozen_a2["gamma_c_est"].values, medium_a2["gamma_c_est"].values])

    fig, ax = plt.subplots(figsize=(6.5, 4))

    # Main plot: all N
    ax.plot(ns, gc, "o-", color=BLUE, markersize=5, linewidth=1.2,
            markeredgecolor="white", markeredgewidth=0.5, zorder=3)

    # O(1) band
    ax.fill_between([11, 46], 0, 1.5, alpha=0.06, color=BLUE, zorder=0)
    ax.text(28, 1.35, r"$O(1)$ band", fontsize=8, color=BLUE, alpha=0.6, ha="center")

    # Annotate outlier
    n10_gc = gc[ns == 10][0]
    ax.annotate("finite-size\neffect", xy=(10, n10_gc), xytext=(16, 7.2),
                fontsize=7, color=GRAY,
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=0.7),
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=GRAY, lw=0.5, alpha=0.8))

    ax.set_xlabel(r"Poset size $N$")
    ax.set_ylabel(r"Critical coupling $\gamma_c$")
    ax.set_xticks(ns)
    ax.set_xlim(8, 46)
    ax.set_ylim(0, 10)
    ax.grid(True, alpha=0.2, linewidth=0.5)

    # Inset: N >= 12 detail
    axins = ax.inset_axes([0.42, 0.45, 0.53, 0.48])
    mask = ns >= 12
    axins.plot(ns[mask], gc[mask], "o-", color=BLUE, markersize=5, linewidth=1.2,
               markeredgecolor="white", markeredgewidth=0.5)
    axins.fill_between([11, 46], 0, 1.5, alpha=0.06, color=BLUE)
    axins.set_xlim(11, 46)
    axins.set_ylim(0, 1.5)
    axins.set_xlabel(r"$N$", fontsize=8)
    axins.set_ylabel(r"$\gamma_c$", fontsize=8)
    axins.tick_params(labelsize=7)
    axins.grid(True, alpha=0.2, linewidth=0.4)
    axins.set_title(r"$N \geq 12$ detail", fontsize=8, pad=3)
    ax.indicate_inset_zoom(axins, edgecolor="0.5", linewidth=0.6)

    fig.savefig(OUT_DIR / "fig1_gamma_c_curve.png")
    fig.savefig(OUT_DIR / "fig1_gamma_c_curve.pdf")
    plt.close(fig)
    print("Fig 1 saved.")


# ── Fig 2: Ablation heatmap ────────────────────────────────────────

def fig2_ablation_summary():
    df = pd.read_csv(
        Path("d:/Kiro/outputs_exploratory/geometric_ablation_gamma_c/"
             "geometric_ablation_gamma_c_report.csv")
    )

    variants = [
        ("A2_full",                        "A2 full (baseline)"),
        ("drop_geo_width_height",          r"$-\,g_{\rm wh}$  (width-height)"),
        ("drop_geo_dim_proxy_penalty",     r"$-\,g_{\rm dim}$ (dim proxy)"),
        ("drop_geo_interval_shape",        r"$-\,g_{\rm int}$  (interval shape)"),
        ("drop_geo_interval_profile",      r"$-\,g_{\rm prf}$   (interval profile)"),
        ("drop_geo_comparability_window",  r"$-\,g_{\rm cmp}$ (comparability)"),
        ("drop_geo_cover_density",         r"$-\,g_{\rm cov}$  (cover density)"),
        ("drop_geo_layer_smoothness",      r"$-\,g_{\rm lyr}$   (layer smooth)"),
        ("A1_neutral_only",                "A1 neutral only"),
    ]
    n_values = sorted(df["n"].unique())

    # Build matrix: rows = variants, cols = N
    matrix = np.full((len(variants), len(n_values)), np.nan)
    status_matrix = np.full((len(variants), len(n_values)), "", dtype=object)

    for i, (vname, _) in enumerate(variants):
        sub = df[df["variant"] == vname]
        for j, n in enumerate(n_values):
            row = sub[sub["n"] == n]
            if row.empty or row.iloc[0]["status"] != "crossing":
                matrix[i, j] = 0.0
                status_matrix[i, j] = "no"
            else:
                matrix[i, j] = float(row.iloc[0]["gamma_c_est"])
                status_matrix[i, j] = "yes"

    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Use a diverging colormap; 0 = no crossing gets its own color
    cmap = plt.cm.YlOrRd.copy()
    cmap.set_under("#e5e7eb")  # light gray for no crossing

    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0.01, vmax=2.0,
                   interpolation="nearest")

    # Cell annotations
    for i in range(len(variants)):
        for j in range(len(n_values)):
            if status_matrix[i, j] == "no":
                ax.text(j, i, "—", ha="center", va="center", fontsize=8,
                        color="#9ca3af", fontweight="bold")
            else:
                val = matrix[i, j]
                text_color = "white" if val > 1.2 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=7, color=text_color)

    ax.set_xticks(range(len(n_values)))
    ax.set_xticklabels(n_values)
    ax.set_yticks(range(len(variants)))
    ax.set_yticklabels([lbl for _, lbl in variants], fontsize=8)
    ax.set_xlabel(r"Poset size $N$")
    ax.set_title(r"Ablation: $\gamma_c$ under removal of geometric sub-terms")

    # Highlight critical rows
    for row_idx in [1, 2]:  # width_height, dim_proxy
        rect = plt.Rectangle((-0.5, row_idx - 0.5), len(n_values), 1,
                              linewidth=1.5, edgecolor=RED, facecolor="none",
                              linestyle="--", zorder=4)
        ax.add_patch(rect)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label(r"$\gamma_c$", fontsize=10)
    cbar.ax.tick_params(labelsize=8)

    fig.savefig(OUT_DIR / "fig2_ablation_summary.png")
    fig.savefig(OUT_DIR / "fig2_ablation_summary.pdf")
    plt.close(fig)
    print("Fig 2 saved.")


# ── Fig 3: Non-circular replacement ────────────────────────────────

def fig3_noncircular_replacement():
    df = pd.read_csv(
        Path("d:/Kiro/outputs_exploratory/noncyclic_dim_replacement_gamma_c/"
             "noncyclic_dim_replacement_gamma_c_report.csv")
    )

    configs = [
        ("A2_full",
         r"$A_2$ full (original)", BLUE, "o", 1.5),
        ("replace_dim_with_consistency_scale_matched",
         r"$g_{\rm dim} \to g_{\rm con}$ (matched)", RED, "s", 1.5),
        ("width_height_plus_consistency_only",
         r"$g_{\rm wh} + g_{\rm con}$ only", GREEN, "^", 1.2),
        ("width_height_only",
         r"$g_{\rm wh}$ only", PURPLE, "D", 1.0),
        ("drop_dim_proxy",
         r"$-\,g_{\rm dim}$", AMBER, "v", 1.0),
    ]

    fig, ax = plt.subplots(figsize=(6.5, 4))

    for variant, label, color, marker, lw in configs:
        sub = df[df["variant"] == variant]
        crossing_ns, crossing_gc = [], []
        missing_ns = []
        for _, row in sub.iterrows():
            n = int(row["n"])
            if row["status"] == "crossing" and pd.notna(row["gamma_c_est"]):
                crossing_ns.append(n)
                crossing_gc.append(float(row["gamma_c_est"]))
            else:
                missing_ns.append(n)

        if crossing_gc:
            ax.plot(crossing_ns, crossing_gc, f"{marker}-", color=color,
                    markersize=5, linewidth=lw, label=label, alpha=0.9,
                    markeredgecolor="white", markeredgewidth=0.4)
        for mn in missing_ns:
            ax.scatter(mn, -0.04, marker="x", color=color, s=30,
                       linewidths=1.2, zorder=5, clip_on=False)

    # Horizontal reference
    ax.axhline(y=1.0, color=GRAY, linestyle=":", linewidth=0.6, alpha=0.4)

    ax.set_xlabel(r"Poset size $N$")
    ax.set_ylabel(r"$\gamma_c$")
    ax.legend(loc="upper left", fontsize=7.5, handlelength=2)
    ax.grid(True, alpha=0.2, linewidth=0.5)
    ax.set_ylim(-0.15, 2.0)
    ax.set_xlim(18, 46)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(0.5))

    # Annotation: "x = no crossing"
    ax.text(44.5, -0.08, r"$\times$ = no crossing", fontsize=7, color=GRAY,
            ha="right", va="top")

    fig.savefig(OUT_DIR / "fig3_noncircular_replacement.png")
    fig.savefig(OUT_DIR / "fig3_noncircular_replacement.pdf")
    plt.close(fig)
    print("Fig 3 saved.")


if __name__ == "__main__":
    fig1_gamma_c_curve()
    fig2_ablation_summary()
    fig3_noncircular_replacement()
    print(f"\nAll figures saved to {OUT_DIR.resolve()}")
