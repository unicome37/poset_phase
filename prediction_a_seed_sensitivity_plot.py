from __future__ import annotations

import argparse
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

COLORS = {
    "A2_full": "#d95f02",
    "A2_replace_dim_with_consistency": "#1b9e77",
    "A2_replace_dim_with_multi_consistency": "#7570b3",
}


def resolve_dir(base: Path, relative: str | Path) -> Path:
    path = Path(relative)
    if path.is_absolute():
        return path
    script_dir = Path(__file__).resolve().parent
    candidates = [
        script_dir / path,
        script_dir.parent / path,
        base / path,
        base.parent / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return script_dir / path


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot Prediction A seed sensitivity summary.")
    parser.add_argument(
        "--input-dir",
        default="outputs_exploratory/prediction_a_seed_sensitivity_n64",
        help="Directory containing prediction_a_seed_sensitivity_*.csv outputs.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    base = Path(__file__).resolve().parents[2]
    in_dir = resolve_dir(base, args.input_dir)
    winners = pd.read_csv(in_dir / "prediction_a_seed_sensitivity_winners.csv")
    pairwise = pd.read_csv(in_dir / "prediction_a_seed_sensitivity_pairwise.csv")
    n_values = sorted(set(int(x) for x in winners["n"].unique()))
    n_label = ",".join(str(x) for x in n_values)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)

    ax = axes[0]
    pivot = (
        winners.groupby(["seed_label", "variant", "winner_family"]).size().rename("wins").reset_index()
    )
    seed_labels = sorted(pivot["seed_label"].unique())
    x = range(len(seed_labels))
    width = 0.22
    offsets = {
        "A2_full": -width,
        "A2_replace_dim_with_consistency": 0.0,
        "A2_replace_dim_with_multi_consistency": width,
    }
    for variant in VARIANTS:
        part = pivot[(pivot["variant"] == variant) & (pivot["winner_family"] == "lorentzian_like_4d")]
        counts = [int(part.loc[part["seed_label"] == label, "wins"].iloc[0]) if any(part["seed_label"] == label) else 0 for label in seed_labels]
        xpos = [i + offsets[variant] for i in x]
        ax.bar(xpos, counts, width=width, color=COLORS[variant], label=variant)
    ax.set_xticks(list(x))
    ax.set_xticklabels(seed_labels, rotation=15)
    ax.set_ylabel("Lor4D wins across 7 gamma values")
    ax.set_title(f"Winner Robustness at N={n_label}")
    ax.set_ylim(0, 7.5)

    ax = axes[1]
    subset = pairwise[
        (pairwise["right_family"] == "lorentzian_like_3d")
        & (pairwise["variant"].isin(["A2_replace_dim_with_consistency", "A2_replace_dim_with_multi_consistency"]))
    ].copy()
    subset["lor4d_margin"] = -subset["delta_score"]
    for variant in ["A2_replace_dim_with_consistency", "A2_replace_dim_with_multi_consistency"]:
        part = subset[subset["variant"] == variant]
        for seed_label, group in part.groupby("seed_label", sort=True):
            ax.plot(
                group["gamma"],
                group["lor4d_margin"],
                marker="o",
                linewidth=1.8,
                alpha=0.8,
                color=COLORS[variant],
                label=f"{variant} / {seed_label}",
            )
    ax.axhline(0.0, color="#666666", linestyle="--", linewidth=1)
    ax.set_xlabel("gamma")
    ax.set_ylabel("ΔA (Lor3D - Lor4D)")
    ax.set_title("Seed-Sensitive Margin vs Lor3D")

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.08))

    suffix = "_".join(str(x) for x in n_values)
    out_path = in_dir / f"prediction_a_seed_sensitivity_n{suffix}.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(out_path.as_posix())


if __name__ == "__main__":
    main()
