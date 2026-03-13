#!/usr/bin/env python3
"""Generate manuscript figures for Prediction C paper.

Produces three figures that parallel the preA figure style:
  fig1_simpsons_paradox.png/pdf  — Simpson's Paradox visualisation
  fig2_tier2_delta_correlation.png/pdf — Tier 2 matched-pair Δ-correlation
  fig3_component_decomposition.png/pdf — Cross-tier component |r| bar chart
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

# ── paths ──
ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT  = pathlib.Path(__file__).resolve().parent / "manuscript_figures"
OUT.mkdir(exist_ok=True)

TIER1_RAW = ROOT / "outputs_exploratory" / "prediction_c_comprehensive" / "tier1_all_family_raw.csv"
TIER2_RAW = ROOT / "outputs_exploratory" / "prediction_c_comprehensive" / "tier2_pairwise_raw.csv"
TIER3_RAW = ROOT / "outputs_exploratory" / "prediction_c_comprehensive" / "tier3_cg_linkage_raw.csv"

# ── global style (match preA) ──
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "serif",
})

PALETTE_N = {10: "#1f77b4", 12: "#ff7f0e", 14: "#2ca02c", 16: "#d62728"}
MARKER_N  = {10: "o", 12: "s", 14: "^", 16: "D"}


# ═══════════════════════════════════════════════════════════════════
# Figure 1 — Simpson's Paradox
# ═══════════════════════════════════════════════════════════════════
def fig1_simpsons_paradox():
    df = pd.read_csv(TIER1_RAW)
    fig, ax = plt.subplots(figsize=(7, 5))

    # per-N scatter + regression
    for n_val in sorted(df["n"].unique()):
        sub = df[df["n"] == n_val]
        c = PALETTE_N[n_val]
        m = MARKER_N[n_val]
        ax.scatter(sub["hierarchy_integration_index"], sub["log_H"],
                   c=c, marker=m, s=28, alpha=0.55, edgecolors="none",
                   label=f"N = {n_val}")
        # within-N regression line
        slope, intercept, r, p, se = stats.linregress(
            sub["hierarchy_integration_index"], sub["log_H"])
        xs = np.linspace(sub["hierarchy_integration_index"].min(),
                         sub["hierarchy_integration_index"].max(), 50)
        ax.plot(xs, slope * xs + intercept, color=c, linewidth=1.5, alpha=0.8)

    # aggregate (naïve) regression line — dashed grey
    slope_all, intercept_all, r_all, p_all, se_all = stats.linregress(
        df["hierarchy_integration_index"], df["log_H"])
    xs_all = np.linspace(df["hierarchy_integration_index"].min(),
                         df["hierarchy_integration_index"].max(), 80)
    ax.plot(xs_all, slope_all * xs_all + intercept_all,
            color="grey", linewidth=2, linestyle="--", alpha=0.7,
            label=f"Aggregate (r = +{r_all:.2f})")

    ax.set_xlabel("Hierarchy Integration Index (HII)")
    ax.set_ylabel("log H (combinatorial entropy)")
    ax.set_title("Simpson's Paradox: Aggregate vs Fixed-N Correlations")

    # custom legend
    handles, labels = ax.get_legend_handles_labels()
    # add annotation for within-N direction
    handles.append(Line2D([0], [0], color="black", lw=1.5, alpha=0.8))
    labels.append("Within-N regression (all negative)")
    ax.legend(handles, labels, loc="upper left", framealpha=0.9)

    # annotate sign flip
    ax.annotate("N controls the sign:\naggregate r > 0\nwithin-N  r < 0",
                xy=(0.97, 0.03), xycoords="axes fraction",
                ha="right", va="bottom",
                fontsize=9, fontstyle="italic",
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow",
                          ec="grey", alpha=0.85))

    for fmt in ("png", "pdf"):
        fig.savefig(OUT / f"fig1_simpsons_paradox.{fmt}")
    plt.close(fig)
    print(f"  ✓ fig1_simpsons_paradox saved")


# ═══════════════════════════════════════════════════════════════════
# Figure 2 — Tier 2 matched-pair Δ-correlation
# ═══════════════════════════════════════════════════════════════════
def fig2_tier2_delta():
    df = pd.read_csv(TIER2_RAW)
    fig, ax = plt.subplots(figsize=(7, 5))

    # colour by N
    palette_t2 = {30: "#1f77b4", 40: "#ff7f0e", 44: "#2ca02c",
                  48: "#d62728", 52: "#9467bd", 56: "#8c564b"}

    for n_val in sorted(df["n"].unique()):
        sub = df[df["n"] == n_val]
        c = palette_t2.get(n_val, "grey")
        ax.scatter(sub["hii_delta"], sub["log_H_delta_mlr_minus_lor2d"],
                   c=c, s=50, alpha=0.75, edgecolors="white", linewidth=0.5,
                   label=f"N = {n_val}  (n={len(sub)})", zorder=3)

    # overall regression
    x = df["hii_delta"].values
    y = df["log_H_delta_mlr_minus_lor2d"].values
    slope, intercept, r, p, se = stats.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 80)
    ax.plot(xs, slope * xs + intercept, color="black", lw=2, zorder=2)

    ax.set_xlabel("ΔHII  (MLR − Lor2D)")
    ax.set_ylabel("Δlog H  (MLR − Lor2D)")
    ax.set_title(f"Tier 2: Matched-Pair Δ-Correlation  (r = {r:.3f}, p < 0.001)")

    ax.axhline(0, color="grey", lw=0.5, ls=":")
    ax.axvline(0, color="grey", lw=0.5, ls=":")

    ax.legend(loc="upper left", framealpha=0.9)

    for fmt in ("png", "pdf"):
        fig.savefig(OUT / f"fig2_tier2_delta_correlation.{fmt}")
    plt.close(fig)
    print(f"  ✓ fig2_tier2_delta_correlation saved")


# ═══════════════════════════════════════════════════════════════════
# Figure 3 — Component decomposition cross-tier bar chart
# ═══════════════════════════════════════════════════════════════════
def fig3_component_decomposition():
    """Reproduce Table 12 from the paper as a grouped bar chart."""
    components = ["layer_count", "mean_layer_gap", "long_edge_frac",
                  "adj_edge_frac", "red_edge_dens", "HII (composite)"]
    # values from the paper (Table 12)
    fixed_n14 = [0.787, 0.745, 0.771, 0.771, 0.693, 0.649]
    tier2_delta = [0.816, 0.836, 0.643, 0.643, 0.459, 0.834]
    tier3_cg    = [0.874, 0.847, 0.803, 0.803, 0.571, 0.820]

    x = np.arange(len(components))
    width = 0.25

    fig, ax = plt.subplots(figsize=(9, 5))
    bars1 = ax.bar(x - width, fixed_n14,  width, label="Fixed N = 14",
                   color="#1f77b4", alpha=0.85, edgecolor="white")
    bars2 = ax.bar(x,         tier2_delta, width, label="Tier 2 Δ (46 pairs)",
                   color="#ff7f0e", alpha=0.85, edgecolor="white")
    bars3 = ax.bar(x + width, tier3_cg,    width, label="Tier 3 CG (92 samples)",
                   color="#2ca02c", alpha=0.85, edgecolor="white")

    ax.set_ylabel("|r|  (absolute correlation)")
    ax.set_title("Component-Level |r| Across Three Independent Analyses")
    ax.set_xticks(x)
    ax.set_xticklabels(components, rotation=25, ha="right", fontsize=10)
    ax.set_ylim(0, 1.0)
    ax.axhline(0.8, color="grey", lw=0.8, ls="--", alpha=0.5)
    ax.legend(loc="upper right", framealpha=0.9)

    # value labels
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            h = bar.get_height()
            ax.annotate(f"{h:.2f}",
                        xy=(bar.get_x() + bar.get_width() / 2, h),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7.5)

    for fmt in ("png", "pdf"):
        fig.savefig(OUT / f"fig3_component_decomposition.{fmt}")
    plt.close(fig)
    print(f"  ✓ fig3_component_decomposition saved")


# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("Generating Prediction C manuscript figures …")
    fig1_simpsons_paradox()
    fig2_tier2_delta()
    fig3_component_decomposition()
    print(f"\nAll figures saved to {OUT}/")
