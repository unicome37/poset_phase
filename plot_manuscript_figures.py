"""Generate publication-quality figures for the manuscript.

Fig 1: gamma_c(N) curve under A2 (confirmatory mainline) with inset
Fig 2: Ablation summary – heatmap matrix
Fig 3: Non-circular replacement comparison
Fig 4: Exact timing frontier
Fig 5: N=52/56 mixed Lor2D vs KR comparison
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


# ── Fig 4: Exact timing frontier ───────────────────────────────────

def fig4_exact_timing_frontier():
    confirm = pd.read_csv("outputs_confirmatory/exact_timing/exact_timing_summary.csv")
    frontier = pd.read_csv(
        "outputs_exploratory/lor2d_exact_timing_frontier/lor2d_exact_timing_frontier.csv"
    )

    fig, ax = plt.subplots(figsize=(6.5, 4))

    # Lor2D frontier across all currently probed sizes
    ax.semilogy(
        frontier["n"],
        frontier["time_seconds"],
        "o-",
        color=BLUE,
        linewidth=1.6,
        markersize=5,
        markeredgecolor="white",
        markeredgewidth=0.5,
        label="Lor2D exact frontier",
        zorder=3,
    )

    # Confirmatory comparison families
    kr = confirm.dropna(subset=["KR_like"])[["n", "KR_like"]].copy()
    tp = confirm.dropna(subset=["transitive_percolation"])[["n", "transitive_percolation"]].copy()

    ax.semilogy(
        kr["n"],
        kr["KR_like"],
        "s--",
        color=RED,
        linewidth=1.2,
        markersize=4.5,
        markeredgecolor="white",
        markeredgewidth=0.4,
        label="KR exact (confirmatory)",
        alpha=0.95,
    )
    ax.semilogy(
        tp["n"],
        tp["transitive_percolation"],
        "D-.",
        color=GRAY,
        linewidth=1.1,
        markersize=4.2,
        markeredgecolor="white",
        markeredgewidth=0.4,
        label="TP exact (confirmatory)",
        alpha=0.9,
    )

    milestones = frontier[frontier["n"].isin([48, 72, 88, 104])].copy()
    ax.scatter(
        milestones["n"],
        milestones["time_seconds"],
        color=AMBER,
        s=28,
        zorder=4,
        label="Lor2D milestones",
    )
    for _, row in milestones.iterrows():
        ax.annotate(
            f"{int(row['n'])}",
            xy=(row["n"], row["time_seconds"]),
            xytext=(0, 7),
            textcoords="offset points",
            fontsize=7,
            color=GRAY,
            ha="center",
        )

    ax.axvline(44, color=GRAY, linestyle=":", linewidth=0.7, alpha=0.5)
    ax.text(45.2, 70, "confirmatory\nboundary", fontsize=7, color=GRAY, va="center")

    ax.set_xlabel(r"Poset size $N$")
    ax.set_ylabel("Exact runtime (s, log scale)")
    ax.set_title("Asymmetric exact frontier: Lor2D stays tractable beyond KR")
    ax.set_xlim(8, 108)
    ax.grid(True, which="both", alpha=0.2, linewidth=0.5)
    ax.legend(loc="upper left", fontsize=7.5, frameon=True)

    fig.savefig(OUT_DIR / "fig4_exact_timing_frontier.png")
    fig.savefig(OUT_DIR / "fig4_exact_timing_frontier.pdf")
    plt.close(fig)
    print("Fig 4 saved.")


# ── Fig 5: N=52/56 mixed Lor2D vs KR comparison ───────────────────

def fig5_mixed_lor2d_vs_kr():
    n52 = pd.read_csv(
        "outputs_exploratory/prediction_a_n52_mixed/prediction_a_ablation_summary.csv"
    )
    n56 = pd.read_csv(
        "outputs_exploratory/prediction_a_n56_mixed/prediction_a_ablation_summary.csv"
    )
    df = pd.concat([n52, n56], ignore_index=True)
    df = df[
        (df["variant"] == "A2_full")
        & (df["family"].isin(["lorentzian_like_2d", "KR_like"]))
        & (df["n"].isin([52, 56]))
    ].copy()

    pivot_mean = (
        df.pivot_table(index=["n", "gamma"], columns="family", values="mean_score")
        .reset_index()
        .rename_axis(None, axis=1)
    )
    pivot_std = (
        df.pivot_table(index=["n", "gamma"], columns="family", values="std_score")
        .reset_index()
        .rename_axis(None, axis=1)
    )
    merged = pivot_mean.merge(
        pivot_std,
        on=["n", "gamma"],
        suffixes=("_mean", "_std"),
    )
    merged["delta_kr_minus_lor2d"] = (
        merged["KR_like_mean"] - merged["lorentzian_like_2d_mean"]
    )
    merged.to_csv(OUT_DIR / "fig5_mixed_lor2d_vs_kr.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharey=True)

    for ax, n in zip(axes, [52, 56]):
        sub = (
            df[df["n"] == n]
            .sort_values(["family", "gamma"])
            .copy()
        )
        lor = sub[sub["family"] == "lorentzian_like_2d"]
        kr = sub[sub["family"] == "KR_like"]

        ax.plot(
            lor["gamma"], lor["mean_score"],
            "o-", color=BLUE, linewidth=1.6, markersize=4.8,
            markeredgecolor="white", markeredgewidth=0.5,
            label="Lor2D",
            zorder=3,
        )
        ax.fill_between(
            lor["gamma"],
            lor["mean_score"] - lor["std_score"],
            lor["mean_score"] + lor["std_score"],
            color=BLUE, alpha=0.12, linewidth=0,
        )

        ax.plot(
            kr["gamma"], kr["mean_score"],
            "s-", color=RED, linewidth=1.6, markersize=4.6,
            markeredgecolor="white", markeredgewidth=0.5,
            label="KR-like",
            zorder=3,
        )
        ax.fill_between(
            kr["gamma"],
            kr["mean_score"] - kr["std_score"],
            kr["mean_score"] + kr["std_score"],
            color=RED, alpha=0.12, linewidth=0,
        )

        delta_last = float(
            kr.loc[kr["gamma"] == 2.0, "mean_score"].iloc[0]
            - lor.loc[lor["gamma"] == 2.0, "mean_score"].iloc[0]
        )
        ax.axvline(2.0, color=GRAY, linestyle=":", linewidth=0.7, alpha=0.5)
        ax.annotate(
            rf"$\Delta A_{{KR-L2D}}={delta_last:.2f}$",
            xy=(2.0, lor.loc[lor["gamma"] == 2.0, "mean_score"].iloc[0]),
            xytext=(-6, 12 if n == 52 else 16),
            textcoords="offset points",
            fontsize=7,
            color=GRAY,
            ha="right",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="0.8", lw=0.5, alpha=0.95),
        )
        ax.text(
            0.06, 0.08,
            "no crossing up to $\\gamma=2.0$",
            transform=ax.transAxes,
            fontsize=7,
            color=GRAY,
        )
        ax.set_title(rf"$N={n}$")
        ax.set_xlabel(r"$\gamma$")
        ax.grid(True, alpha=0.2, linewidth=0.5)
        ax.set_xlim(-0.03, 2.05)
        ax.set_xticks([0.0, 0.4, 0.8, 1.2, 1.6, 2.0])

    axes[0].set_ylabel(r"Mean action score under $A_2^{\rm full}$")
    axes[0].legend(loc="upper right", fontsize=8, frameon=True)

    fig.suptitle(
        r"Near-wall mixed comparison: Lor2D vs KR under $A_2^{\rm full}$",
        y=1.02,
        fontsize=12,
    )
    fig.savefig(OUT_DIR / "fig5_mixed_lor2d_vs_kr.png")
    fig.savefig(OUT_DIR / "fig5_mixed_lor2d_vs_kr.pdf")
    plt.close(fig)
    print("Fig 5 saved.")


if __name__ == "__main__":
    fig1_gamma_c_curve()
    fig2_ablation_summary()
    fig3_noncircular_replacement()
    fig4_exact_timing_frontier()
    fig5_mixed_lor2d_vs_kr()
    print(f"\nAll figures saved to {OUT_DIR.resolve()}")
