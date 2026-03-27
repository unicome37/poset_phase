from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "manuscript_figures"
UNIFIED_DIR = ROOT / "outputs_unified_functional"
OCC_DIR = ROOT / "outputs_exploratory" / "prediction_a_xi_occupancy_closure"
OUT_DIR.mkdir(exist_ok=True)


plt.rcParams.update({
    "font.size": 10,
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.labelsize": 10.5,
    "axes.titlesize": 10.5,
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


def build_lower_panel(n_value: int = 36) -> pd.DataFrame:
    occ = pd.read_csv(OCC_DIR / "occupancy_variants.csv")
    sub = occ[occ["n"] == n_value].copy()
    sub["dimension"] = sub["d"].astype(int)
    sub["wall_proxy"] = sub["P_occ"]
    return sub[["dimension", "wall_proxy", "nonlink_mc", "P_gam"]].sort_values("dimension")


def build_upper_panel(n_value: int = 36) -> pd.DataFrame:
    pa = pd.read_csv(UNIFIED_DIR / "prediction_A_dimension.csv")
    rf = pd.read_csv(UNIFIED_DIR / "raw_features.csv")

    merged = pa.merge(
        rf[["family", "N", "rep", "sigma_hist"]],
        on=["family", "N", "rep"],
        how="left",
    )
    merged = merged[merged["N"] == n_value].copy()
    merged["wall_reconstructed"] = (
        merged["F_total"]
        - merged["log_H"]
        + 10.0 * merged["sigma_hist"]
        - 0.6 * merged["xi_dim"]
    )

    ordered = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]
    summary = (
        merged[merged["family"].isin(ordered)]
        .groupby("family", as_index=False)
        .agg(
            F_total_mean=("F_total", "mean"),
            xi_dim_mean=("xi_dim", "mean"),
            wall_reconstructed_mean=("wall_reconstructed", "mean"),
            d_eff_mean=("d_eff", "mean"),
        )
    )
    summary["family"] = pd.Categorical(summary["family"], categories=ordered, ordered=True)
    return summary.sort_values("family")


def save_data(lower: pd.DataFrame, upper: pd.DataFrame) -> None:
    lower.to_csv(OUT_DIR / "cqg_fig1_lower_dimensional_proxy.csv", index=False)
    upper.to_csv(OUT_DIR / "cqg_fig1_upper_dimensional_scores.csv", index=False)


def make_figure(lower: pd.DataFrame, upper: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))

    # Lower-dimensional exclusion proxy
    axes[0].plot(lower["dimension"], lower["wall_proxy"], marker="o", markersize=4.5, color="#b91c1c", linewidth=1.7, label=r"$P_{\mathrm{occ}}$")
    axes[0].plot(lower["dimension"], lower["nonlink_mc"], marker="s", markersize=4.0, color="#64748b", linewidth=1.2, linestyle="--", label="nonlink fraction")
    axes[0].set_xlabel("Embedding dimension")
    axes[0].set_ylabel("Occupancy proxy")
    axes[0].set_title("(a) Lower-dimensional exclusion", loc="left", pad=6)
    axes[0].grid(alpha=0.25, linewidth=0.5)
    axes[0].set_xticks(lower["dimension"])
    axes[0].legend(loc="upper right")

    # Upper-dimensional suppression and total score
    x = range(len(upper))
    axes[1].plot(x, upper["F_total_mean"], marker="o", markersize=4.5, color="#1d4ed8", linewidth=1.7, label=r"Total score $\langle \mathcal{F}\rangle$")
    axes[1].plot(x, upper["xi_dim_mean"], marker="^", markersize=4.5, color="#15803d", linewidth=1.4, label=r"$\langle \Xi_d \rangle$")
    axes[1].set_xticks(list(x))
    axes[1].set_xticklabels(upper["family"])
    axes[1].set_xlabel("Candidate class")
    axes[1].set_ylabel("Score contribution")
    axes[1].set_title("(b) Upper-dimensional suppression", loc="left", pad=6)
    axes[1].grid(alpha=0.25, linewidth=0.5)
    axes[1].legend(loc="upper left")

    fig.subplots_adjust(wspace=0.28)
    fig.savefig(OUT_DIR / "cqg_fig1_dimensional_selection.png")
    fig.savefig(OUT_DIR / "cqg_fig1_dimensional_selection.pdf")
    plt.close(fig)


def main() -> None:
    lower = build_lower_panel()
    upper = build_upper_panel()
    save_data(lower, upper)
    make_figure(lower, upper)
    print("Saved:", OUT_DIR / "cqg_fig1_dimensional_selection.png")
    print("Saved:", OUT_DIR / "cqg_fig1_dimensional_selection.pdf")
    print("Saved:", OUT_DIR / "cqg_fig1_lower_dimensional_proxy.csv")
    print("Saved:", OUT_DIR / "cqg_fig1_upper_dimensional_scores.csv")


if __name__ == "__main__":
    main()
