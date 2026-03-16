from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path("d:/Kiro")
POSSET_ROOT = ROOT / "理论体系" / "poset_phase"
OUT_DIR = POSSET_ROOT / "outputs_exploratory" / "lor2d_exact_timing_frontier"

SOURCE_CSVS = [
    POSSET_ROOT / "outputs_confirmatory" / "exact_timing" / "exact_timing_summary.csv",
    POSSET_ROOT / "outputs_exploratory" / "exact_timing_lor2d_extension" / "exact_timing_summary.csv",
    ROOT / "outputs_exploratory" / "exact_timing_lor2d_extension_n80_n88" / "exact_timing_summary.csv",
    ROOT / "outputs_exploratory" / "exact_timing_lor2d_extension_n96_n104" / "exact_timing_summary.csv",
]


def load_lor2d_frontier() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for csv_path in SOURCE_CSVS:
        df = pd.read_csv(csv_path)
        if "lorentzian_like_2d" in df.columns:
            sub = df[["n", "lorentzian_like_2d"]].copy()
            sub = sub.rename(columns={"lorentzian_like_2d": "time_seconds"})
        else:
            sub = df[["n", "lorentzian_like_2d"]].copy()
            sub = sub.rename(columns={"lorentzian_like_2d": "time_seconds"})
        sub["source_csv"] = str(csv_path)
        frames.append(sub)

    merged = pd.concat(frames, ignore_index=True)
    merged = merged.dropna(subset=["time_seconds"])
    merged = merged.sort_values(["n", "time_seconds"]).drop_duplicates(subset=["n"], keep="last")
    merged["time_class"] = pd.cut(
        merged["time_seconds"],
        bins=[0, 1, 10, 60, float("inf")],
        labels=["sub-second", "seconds", "tens-of-seconds", "minute-scale"],
        include_lowest=True,
    )
    return merged.sort_values("n").reset_index(drop=True)


def plot_frontier(df: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
        }
    )

    fig, ax = plt.subplots(figsize=(6.7, 4.2))
    ax.semilogy(
        df["n"],
        df["time_seconds"],
        color="#2563eb",
        marker="o",
        linewidth=1.8,
        markersize=5,
        markeredgecolor="white",
        markeredgewidth=0.5,
        label="Lor2D exact timing",
    )

    milestones = df[df["n"].isin([48, 72, 88, 104])].copy()
    ax.scatter(
        milestones["n"],
        milestones["time_seconds"],
        color="#dc2626",
        s=28,
        zorder=4,
        label="frontier milestones",
    )
    for _, row in milestones.iterrows():
        ax.annotate(
            f"N={int(row['n'])}\n{row['time_seconds']:.2f}s",
            xy=(row["n"], row["time_seconds"]),
            xytext=(4, 8),
            textcoords="offset points",
            fontsize=7,
            color="#374151",
        )

    ax.axhspan(0.001, 1, color="#93c5fd", alpha=0.12)
    ax.axhspan(1, 10, color="#86efac", alpha=0.10)
    ax.axhspan(10, 60, color="#fde68a", alpha=0.10)
    ax.axhspan(60, 180, color="#fca5a5", alpha=0.09)

    ax.text(103, 0.18, "sub-second", fontsize=7, color="#1d4ed8", ha="right")
    ax.text(103, 3.0, "seconds", fontsize=7, color="#166534", ha="right")
    ax.text(103, 22.0, "tens of seconds", fontsize=7, color="#92400e", ha="right")
    ax.text(103, 92.0, "minute-scale", fontsize=7, color="#991b1b", ha="right")

    ax.set_xlabel(r"Poset size $N$")
    ax.set_ylabel("Exact runtime (s, log scale)")
    ax.set_title("Lor2D single-sided exact timing frontier")
    ax.set_xlim(df["n"].min() - 2, df["n"].max() + 4)
    ax.grid(True, which="both", linestyle="--", alpha=0.25)
    ax.legend(frameon=False, loc="upper left")

    fig.savefig(OUT_DIR / "lor2d_exact_timing_frontier.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_lor2d_frontier()
    df.to_csv(OUT_DIR / "lor2d_exact_timing_frontier.csv", index=False)
    plot_frontier(df)
    print((OUT_DIR / "lor2d_exact_timing_frontier.csv").as_posix())
    print((OUT_DIR / "lor2d_exact_timing_frontier.png").as_posix())


if __name__ == "__main__":
    main()
