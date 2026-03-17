"""
Prediction A — Ξ Core Figure

Boxplot + scatter of Ξ_{2→3}, Ξ_{3→4}, Ξ_{4→5} across three generators.
One-glance visual: 4→5 boundary is dramatically elevated and tightly clustered.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = Path("outputs_exploratory/prediction_a_xi_parameter")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    xi = pd.read_csv(OUT_DIR / "xi_all_transitions.csv")
    print(f"Loaded {len(xi)} Ξ values")
    print(xi.head())

    # ---- Rename suites for display ----
    suite_labels = {
        "cube_original": "Cube (original)",
        "cube_indep_seed": "Cube (indep seed)",
        "diamond": "Causal diamond",
    }
    xi["Generator"] = xi["suite"].map(suite_labels)

    transitions = ["2→3", "3→4", "4→5"]
    colors = {
        "Cube (original)": "#2176AE",
        "Cube (indep seed)": "#57B8FF",
        "Causal diamond": "#F0803C",
    }
    markers = {
        "Cube (original)": "o",
        "Cube (indep seed)": "s",
        "Causal diamond": "D",
    }

    # ========================================================
    # Figure 1: Strip plot with means + transparent band
    # ========================================================
    fig, ax = plt.subplots(figsize=(7, 5))

    positions = {"2→3": 0, "3→4": 1, "4→5": 2}
    jitter_offsets = {"Cube (original)": -0.15, "Cube (indep seed)": 0.0, "Causal diamond": 0.15}

    for gen_label, gen_group in xi.groupby("Generator"):
        for trans in transitions:
            sub = gen_group[gen_group["dim_pair"] == trans]
            if sub.empty:
                continue
            x_base = positions[trans]
            jitter = jitter_offsets[gen_label]
            x_vals = np.full(len(sub), x_base + jitter) + np.random.default_rng(42).uniform(-0.03, 0.03, len(sub))
            ax.scatter(
                x_vals, sub["Xi"],
                c=colors[gen_label], marker=markers[gen_label],
                s=40, alpha=0.6, edgecolors="white", linewidths=0.5,
                label=gen_label if trans == "2→3" else None,
                zorder=3,
            )
            # Mean marker
            mean_val = sub["Xi"].mean()
            ax.scatter(
                [x_base + jitter], [mean_val],
                c=colors[gen_label], marker=markers[gen_label],
                s=120, alpha=1.0, edgecolors="black", linewidths=1.2,
                zorder=4,
            )

    # Horizontal band for Ξ₄→₅ convergence zone
    xi_45 = xi[xi["dim_pair"] == "4→5"]["Xi"]
    band_lo, band_hi = xi_45.quantile(0.25), xi_45.quantile(0.75)
    ax.axhspan(band_lo, band_hi, xmin=0.6, xmax=1.0, color="#FFD166", alpha=0.25, zorder=1)
    ax.axhline(xi_45.median(), xmin=0.6, xmax=1.0, color="#FFD166", linewidth=2, linestyle="--", alpha=0.7, zorder=2)

    # Annotation
    ax.annotate(
        f"median ≈ {xi_45.median():.1f}\nCV = 13.9%",
        xy=(2.2, xi_45.median()), fontsize=9,
        ha="left", va="center",
        bbox=dict(boxstyle="round,pad=0.3", fc="#FFD166", alpha=0.5),
    )

    # Dashed line separating low-Ξ from high-Ξ regime
    ax.axhline(6.0, color="gray", linestyle=":", linewidth=1, alpha=0.5)
    ax.text(2.35, 6.3, "barrier\nthreshold", fontsize=7, color="gray", ha="left", va="bottom")

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels([r"$\Xi_{2\to3}$", r"$\Xi_{3\to4}$", r"$\Xi_{4\to5}$"], fontsize=12)
    ax.set_ylabel(r"$\Xi$ (link-penalty / entropy ratio)", fontsize=11)
    ax.set_title("Dimensionless Control Parameter Across Generators", fontsize=12, fontweight="bold")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.set_xlim(-0.5, 2.8)
    ax.set_ylim(0, 15)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"xi_core_figure.{ext}", dpi=300)
    print(f"Saved xi_core_figure.png/pdf")
    plt.close()

    # ========================================================
    # Figure 2: Ξ vs N for 4→5 only — convergence across N
    # ========================================================
    fig2, ax2 = plt.subplots(figsize=(7, 4.5))

    xi_45_data = xi[xi["dim_pair"] == "4→5"]
    for gen_label, sub in xi_45_data.groupby("Generator"):
        sub_sorted = sub.sort_values("n")
        ax2.plot(sub_sorted["n"], sub_sorted["Xi"], marker=markers[gen_label],
                 color=colors[gen_label], linewidth=1.5, markersize=7,
                 label=gen_label, alpha=0.85)

    # Band
    overall_median = xi_45_data["Xi"].median()
    ax2.axhline(overall_median, color="#FFD166", linewidth=2, linestyle="--", alpha=0.7,
                label=f"median = {overall_median:.1f}")
    ax2.axhspan(overall_median * 0.85, overall_median * 1.15, color="#FFD166", alpha=0.15)

    ax2.set_xlabel("$N$ (poset size)", fontsize=11)
    ax2.set_ylabel(r"$\Xi_{4\to5}$", fontsize=12)
    ax2.set_title(r"$\Xi_{4\to5}$ Convergence Across Generators and Sizes", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=9, loc="lower right")
    ax2.set_ylim(4, 14)
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig2.savefig(OUT_DIR / f"xi_45_convergence.{ext}", dpi=300)
    print(f"Saved xi_45_convergence.png/pdf")
    plt.close()

    # ========================================================
    # Figure 3: Ξ ratio (4→5)/(3→4) — asymmetric barrier
    # ========================================================
    fig3, ax3 = plt.subplots(figsize=(7, 4.5))

    # Compute ratio for each (suite, n)
    for gen_label in ["Cube (original)", "Cube (indep seed)", "Causal diamond"]:
        sub34 = xi[(xi["Generator"] == gen_label) & (xi["dim_pair"] == "3→4")].set_index("n")
        sub45 = xi[(xi["Generator"] == gen_label) & (xi["dim_pair"] == "4→5")].set_index("n")
        common = sorted(set(sub34.index) & set(sub45.index))
        ratios = [sub45.loc[n, "Xi"] / sub34.loc[n, "Xi"] if sub34.loc[n, "Xi"] > 0.1 else np.nan for n in common]
        ax3.plot(common, ratios, marker=markers[gen_label], color=colors[gen_label],
                 linewidth=1.5, markersize=7, label=gen_label, alpha=0.85)

    ax3.axhline(1.0, color="gray", linestyle=":", linewidth=1, alpha=0.5)
    ax3.text(22, 1.15, "symmetric boundary (no asymmetry)", fontsize=8, color="gray", va="bottom")

    ax3.set_xlabel("$N$ (poset size)", fontsize=11)
    ax3.set_ylabel(r"$\Xi_{4\to5}\, /\, \Xi_{3\to4}$", fontsize=12)
    ax3.set_title("Asymmetric Barrier Ratio: 4→5 is Consistently Harder Than 3→4", fontsize=11, fontweight="bold")
    ax3.legend(fontsize=9, loc="upper right")
    ax3.set_ylim(0, 10)
    ax3.grid(alpha=0.3)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig3.savefig(OUT_DIR / f"xi_barrier_asymmetry.{ext}", dpi=300)
    print(f"Saved xi_barrier_asymmetry.png/pdf")
    plt.close()

    # Summary stats
    print("\n=== Summary Statistics ===")
    for trans in transitions:
        sub = xi[xi["dim_pair"] == trans]
        print(f"\n  {trans}:")
        print(f"    Overall: mean={sub['Xi'].mean():.2f}, median={sub['Xi'].median():.2f}, std={sub['Xi'].std():.2f}")
        for suite in sorted(sub["suite"].unique()):
            ss = sub[sub["suite"] == suite]["Xi"]
            print(f"    {suite}: mean={ss.mean():.2f}, median={ss.median():.2f}, std={ss.std():.2f}")


if __name__ == "__main__":
    main()
