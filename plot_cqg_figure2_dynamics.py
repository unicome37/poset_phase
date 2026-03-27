from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "manuscript_figures"
DATA_DIR = ROOT / "outputs_unified_functional"
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


TRAJECTORIES = [
    {
        "path": DATA_DIR / "metropolis_B_trajectory.csv",
        "label_key": "random_rep0",
        "display": "Unconstrained Metropolis",
        "color": "#dc2626",
    },
    {
        "path": DATA_DIR / "metropolis_M_microcanonical.csv",
        "label_key": "mc_Lor2D",
        "display": "Microcanonical swap",
        "color": "#2563eb",
    },
    {
        "path": DATA_DIR / "metropolis_C_trajectory.csv",
        "label_key": "sediment_rep0",
        "display": "Structured initial condition",
        "color": "#16a34a",
    },
]


def load_curves() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for item in TRAJECTORIES:
        df = pd.read_csv(item["path"])
        sub = df[df["label"] == item["label_key"]].copy()
        if sub.empty:
            raise ValueError(f"Missing label {item['label_key']} in {item['path']}")
        sub["trajectory"] = item["display"]
        sub["color"] = item["color"]
        frames.append(sub[["step", "F", "comp_frac", "d_eff", "trajectory", "color"]])
    merged = pd.concat(frames, ignore_index=True)
    merged.to_csv(OUT_DIR / "cqg_fig2_dynamics_curves.csv", index=False)
    return merged


def make_figure(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))

    for trajectory, sub in df.groupby("trajectory", sort=False):
        color = sub["color"].iloc[0]
        axes[0].plot(sub["step"], sub["comp_frac"], color=color, linewidth=1.7, label=trajectory)
        axes[1].plot(sub["step"], sub["d_eff"], color=color, linewidth=1.7, label=trajectory)

    axes[0].axhspan(0.40, 0.55, color="#fde68a", alpha=0.22, zorder=0)
    axes[0].set_xlabel("Step")
    axes[0].set_ylabel("Comparability fraction")
    axes[0].set_title("(a) Comparability fraction", loc="left", pad=6)
    axes[0].grid(alpha=0.25, linewidth=0.5)

    axes[1].axhspan(3.5, 4.5, color="#bfdbfe", alpha=0.22, zorder=0)
    axes[1].set_xlabel("Step")
    axes[1].set_ylabel(r"Effective dimension $d_{\mathrm{eff}}$")
    axes[1].set_title("(b) Effective dimension", loc="left", pad=6)
    axes[1].grid(alpha=0.25, linewidth=0.5)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02), columnspacing=1.4, handlelength=2.2)
    fig.subplots_adjust(wspace=0.28, top=0.80)
    fig.savefig(OUT_DIR / "cqg_fig2_dynamics.png")
    fig.savefig(OUT_DIR / "cqg_fig2_dynamics.pdf")
    plt.close(fig)


def main() -> None:
    df = load_curves()
    make_figure(df)
    print("Saved:", OUT_DIR / "cqg_fig2_dynamics.png")
    print("Saved:", OUT_DIR / "cqg_fig2_dynamics.pdf")
    print("Saved:", OUT_DIR / "cqg_fig2_dynamics_curves.csv")


if __name__ == "__main__":
    main()
