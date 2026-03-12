from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


VARIANTS = [
    "A2_full",
    "A2_replace_dim_with_consistency",
    "A2_replace_dim_with_multi_consistency",
]

VARIANT_LABELS = {
    "A2_full": "A2_full",
    "A2_replace_dim_with_consistency": "A2_replace_dim_with_consistency",
    "A2_replace_dim_with_multi_consistency": "A2_replace_dim_with_multi_consistency",
}

VARIANT_COLORS = {
    "A2_full": "#d95f02",
    "A2_replace_dim_with_consistency": "#1b9e77",
    "A2_replace_dim_with_multi_consistency": "#7570b3",
}


def resolve_csv(base: Path, relative: str) -> Path:
    local = base / relative
    if local.exists():
        return local
    root = base.parents[1] / relative
    if root.exists():
        return root
    return local


def load_pairwise_tables(base: Path) -> pd.DataFrame:
    files = [
        resolve_csv(base, "outputs_exploratory/prediction_a_dim_replacement_sp8/prediction_a_ablation_pairwise.csv"),
        resolve_csv(base, "outputs_exploratory/prediction_a_dim_replacement_n44_n48/prediction_a_ablation_pairwise.csv"),
        resolve_csv(base, "outputs_exploratory/prediction_a_n52_mixed/prediction_a_ablation_pairwise.csv"),
        resolve_csv(base, "outputs_exploratory/prediction_a_n56_mixed/prediction_a_ablation_pairwise.csv"),
        resolve_csv(base, "outputs_exploratory/prediction_a_n60_mixed/prediction_a_ablation_pairwise.csv"),
        resolve_csv(base, "outputs_exploratory/prediction_a_n64_mixed/prediction_a_ablation_pairwise.csv"),
        resolve_csv(base, "outputs_exploratory/prediction_a_n68_mixed/prediction_a_ablation_pairwise.csv"),
        resolve_csv(base, "outputs_exploratory/prediction_a_n72_mixed/prediction_a_ablation_pairwise.csv"),
    ]
    frames = [pd.read_csv(path) for path in files]
    df = pd.concat(frames, ignore_index=True)
    return df[df["variant"].isin(VARIANTS)].copy()


def build_margin_summary(df: pd.DataFrame) -> pd.DataFrame:
    # delta_score is left - right; left family is Lor4D in this pipeline.
    # Negative means Lor4D wins. Convert to positive margin-of-victory when Lor4D wins.
    df = df.copy()
    df["lor4d_margin"] = -df["delta_score"]
    summary = (
        df.groupby(["variant", "right_family", "n"], sort=True)
        .agg(
            mean_margin=("lor4d_margin", "mean"),
            min_margin=("lor4d_margin", "min"),
            max_margin=("lor4d_margin", "max"),
            wins=("left_wins", "sum"),
            count=("left_wins", "count"),
        )
        .reset_index()
    )
    return summary


def plot_margin(summary: pd.DataFrame, out_dir: Path) -> None:
    opponents = ["lorentzian_like_2d", "lorentzian_like_3d"]
    labels = {
        "lorentzian_like_2d": "vs Lor2D",
        "lorentzian_like_3d": "vs Lor3D",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=False, constrained_layout=True)

    for ax, opponent in zip(axes, opponents):
        subset = summary[summary["right_family"] == opponent]
        for variant in VARIANTS:
            part = subset[subset["variant"] == variant].sort_values("n")
            color = VARIANT_COLORS[variant]
            ax.fill_between(
                part["n"],
                part["min_margin"],
                part["max_margin"],
                alpha=0.15,
                color=color,
                linewidth=0,
            )
            ax.plot(
                part["n"],
                part["mean_margin"],
                marker="o",
                linewidth=2,
                color=color,
                label=VARIANT_LABELS[variant],
            )
        ax.axhline(0.0, color="#666666", linewidth=1, linestyle="--")
        ax.set_title(labels[opponent])
        ax.set_xlabel("N")
        ax.set_ylabel("Mean ΔA (runner-up minus Lor4D)")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.05))
    png_path = out_dir / "prediction_a_margin_of_victory.png"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    base = Path(__file__).resolve().parent
    out_dir = base / "outputs_exploratory" / "prediction_a_margin_summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_pairwise_tables(base)
    summary = build_margin_summary(df)
    summary.to_csv(out_dir / "prediction_a_margin_summary.csv", index=False, encoding="utf-8-sig")

    # Also extract the hardest end, gamma = 2.0, for quick manuscript use.
    gamma_hard = (
        df[df["gamma"] == 2.0]
        .assign(lor4d_margin=lambda x: -x["delta_score"])
        .groupby(["variant", "right_family", "n"], sort=True)
        .agg(mean_margin=("lor4d_margin", "mean"), wins=("left_wins", "sum"), count=("left_wins", "count"))
        .reset_index()
    )
    gamma_hard.to_csv(out_dir / "prediction_a_margin_gamma2_summary.csv", index=False, encoding="utf-8-sig")

    plot_margin(summary, out_dir)
    print((out_dir / "prediction_a_margin_of_victory.png").as_posix())


if __name__ == "__main__":
    main()
