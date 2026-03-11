"""Generate publication-quality figures for the manuscript.

Fig 1: gamma_c(N) curve under A2 (confirmatory mainline)
Fig 2: Ablation summary - which sub-terms are critical
Fig 3: Non-circular replacement comparison
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.size": 11,
    "font.family": "serif",
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

OUT_DIR = Path("manuscript_figures")
OUT_DIR.mkdir(exist_ok=True)

# ── Fig 1: gamma_c(N) confirmatory mainline ────────────────────────

def fig1_gamma_c_curve():
    # Frozen exact (N=10-16)
    frozen = pd.read_csv("outputs_confirmatory/frozen_exact/gamma_c_report.csv")
    frozen_a2 = frozen[(frozen["action_mode"] == "A2") & (frozen["family_b"] == "KR_like")]

    # Medium exact scan (N=20-44)
    medium = pd.read_csv("outputs_confirmatory/medium_exact_scan/gamma_c_report.csv")
    medium_a2 = medium[(medium["action_mode"] == "A2") & (medium["family_b"] == "KR_like")]

    ns_frozen = frozen_a2["n"].values
    gc_frozen = frozen_a2["gamma_c_est"].values
    ns_medium = medium_a2["n"].values
    gc_medium = medium_a2["gamma_c_est"].values

    ns = np.concatenate([ns_frozen, ns_medium])
    gc = np.concatenate([gc_frozen, gc_medium])

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(ns, gc, "o-", color="#2563eb", markersize=7, linewidth=1.5, label=r"$\gamma_c$ (A2, exact)")
    ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("Poset size $N$")
    ax.set_ylabel(r"Critical coupling $\gamma_c$")
    ax.set_title(r"$\gamma_c(N)$ for Lor2D vs KR under $A_2$ (exact entropy)")
    ax.set_xticks(ns)
    ax.set_xlim(8, 46)
    ax.set_ylim(0, 10)

    # Annotate N=10 outlier
    n10_idx = np.where(ns == 10)[0]
    if len(n10_idx) > 0:
        ax.annotate("finite-size\noutlier", xy=(10, gc[n10_idx[0]]),
                     xytext=(14, 7.5), fontsize=8, color="gray",
                     arrowprops=dict(arrowstyle="->", color="gray", lw=0.8))

    # Shade the O(1) band for N>=12
    mask = ns >= 12
    if mask.any():
        ax.fill_between([12, 46], 0, 2.0, alpha=0.08, color="#2563eb")
        ax.text(29, 1.7, "O(1) band", fontsize=9, color="#2563eb", alpha=0.7, ha="center")

    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)
    fig.savefig(OUT_DIR / "fig1_gamma_c_curve.png")
    fig.savefig(OUT_DIR / "fig1_gamma_c_curve.pdf")
    plt.close(fig)
    print(f"Fig 1 saved to {OUT_DIR / 'fig1_gamma_c_curve.png'}")


# ── Fig 2: Ablation summary ────────────────────────────────────────

def fig2_ablation_summary():
    df = pd.read_csv(
        Path("d:/Kiro/outputs_exploratory/geometric_ablation_gamma_c/geometric_ablation_gamma_c_report.csv")
    )

    variants_order = [
        "A2_full",
        "drop_geo_width_height",
        "drop_geo_dim_proxy_penalty",
        "drop_geo_interval_shape",
        "drop_geo_interval_profile",
        "drop_geo_comparability_window",
        "drop_geo_cover_density",
        "drop_geo_layer_smoothness",
        "A1_neutral_only",
    ]
    labels = [
        "A2 full",
        "drop width_height",
        "drop dim_proxy",
        "drop interval_shape",
        "drop interval_profile",
        "drop comparability",
        "drop cover_density",
        "drop layer_smooth",
        "A1 neutral only",
    ]
    n_values = sorted(df["n"].unique())

    fig, ax = plt.subplots(figsize=(10, 5))

    bar_width = 0.08
    x_base = np.arange(len(n_values))
    colors = plt.cm.Set2(np.linspace(0, 1, len(variants_order)))

    for i, (variant, label) in enumerate(zip(variants_order, labels)):
        sub = df[df["variant"] == variant]
        gc_vals = []
        for n in n_values:
            row = sub[sub["n"] == n]
            if row.empty or row.iloc[0]["status"] != "crossing":
                gc_vals.append(0)  # no crossing
            else:
                gc_vals.append(float(row.iloc[0]["gamma_c_est"]))
        
        offset = (i - len(variants_order) / 2) * bar_width
        bars = ax.bar(x_base + offset, gc_vals, bar_width, label=label, color=colors[i], edgecolor="white", linewidth=0.5)

        # Mark missing crossings
        for j, v in enumerate(gc_vals):
            if v == 0:
                ax.scatter(x_base[j] + offset, 0.05, marker="x", color="red", s=20, linewidths=1.2, zorder=5)

    ax.set_xlabel("Poset size $N$")
    ax.set_ylabel(r"$\gamma_c$")
    ax.set_title("Ablation: effect of removing geometric sub-terms on $\\gamma_c$")
    ax.set_xticks(x_base)
    ax.set_xticklabels(n_values)
    ax.legend(loc="upper left", fontsize=7, ncol=2)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_ylim(0, 2.5)

    fig.savefig(OUT_DIR / "fig2_ablation_summary.png")
    fig.savefig(OUT_DIR / "fig2_ablation_summary.pdf")
    plt.close(fig)
    print(f"Fig 2 saved to {OUT_DIR / 'fig2_ablation_summary.png'}")


# ── Fig 3: Non-circular replacement comparison ─────────────────────

def fig3_noncircular_replacement():
    df = pd.read_csv(
        Path("d:/Kiro/outputs_exploratory/noncyclic_dim_replacement_gamma_c/noncyclic_dim_replacement_gamma_c_report.csv")
    )

    variants_to_plot = [
        ("A2_full", "A2 full (original)", "#2563eb", "o"),
        ("replace_dim_with_consistency_scale_matched", "dim→consistency (matched)", "#dc2626", "s"),
        ("width_height_plus_consistency_only", "wh + consistency only", "#16a34a", "^"),
        ("width_height_only", "width_height only", "#9333ea", "D"),
        ("drop_dim_proxy", "drop dim_proxy", "#f59e0b", "v"),
    ]

    fig, ax = plt.subplots(figsize=(7, 4.5))

    for variant, label, color, marker in variants_to_plot:
        sub = df[df["variant"] == variant]
        ns = sub["n"].values
        gc = []
        gc_ns = []
        for _, row in sub.iterrows():
            if row["status"] == "crossing" and pd.notna(row["gamma_c_est"]):
                gc_ns.append(int(row["n"]))
                gc.append(float(row["gamma_c_est"]))

        if gc:
            ax.plot(gc_ns, gc, f"{marker}-", color=color, markersize=6, linewidth=1.2, label=label, alpha=0.85)

        # Mark missing crossings with X
        for _, row in sub.iterrows():
            if row["status"] != "crossing":
                ax.scatter(int(row["n"]), 0.05, marker="x", color=color, s=40, linewidths=1.5, zorder=5)

    ax.set_xlabel("Poset size $N$")
    ax.set_ylabel(r"$\gamma_c$")
    ax.set_title("Non-circular replacement: $\\gamma_c$ under different backbones")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 2.2)
    ax.set_xlim(18, 46)

    fig.savefig(OUT_DIR / "fig3_noncircular_replacement.png")
    fig.savefig(OUT_DIR / "fig3_noncircular_replacement.pdf")
    plt.close(fig)
    print(f"Fig 3 saved to {OUT_DIR / 'fig3_noncircular_replacement.png'}")


if __name__ == "__main__":
    fig1_gamma_c_curve()
    fig2_ablation_summary()
    fig3_noncircular_replacement()
    print(f"\nAll figures saved to {OUT_DIR.resolve()}")
