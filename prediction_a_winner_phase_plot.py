from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap


VARIANTS = [
    "A2_full",
    "A2_replace_dim_with_consistency",
    "A2_replace_dim_with_multi_consistency",
]

FAMILY_TO_CODE = {
    "lorentzian_like_2d": 0,
    "lorentzian_like_3d": 1,
    "lorentzian_like_4d": 2,
    "KR_like": 3,
}

CODE_TO_LABEL = {
    0: "Lor2D",
    1: "Lor3D",
    2: "Lor4D",
    3: "KR",
}

COLORS = ["#d95f02", "#7570b3", "#1b9e77", "#666666"]


def resolve_csv(base: Path, relative: str) -> Path:
    local = base / relative
    if local.exists():
        return local
    root = base.parents[1] / relative
    if root.exists():
        return root
    return local


def load_winner_tables(base: Path) -> pd.DataFrame:
    early = pd.read_csv(
        resolve_csv(base, "outputs_exploratory/prediction_a_dim_replacement_sp8/prediction_a_ablation_winners.csv")
    )
    late = pd.read_csv(
        resolve_csv(base, "outputs_exploratory/prediction_a_dim_replacement_n44_n48/prediction_a_ablation_winners.csv")
    )
    n52 = pd.read_csv(
        resolve_csv(base, "outputs_exploratory/prediction_a_n52_mixed/prediction_a_ablation_winners.csv")
    )
    n56 = pd.read_csv(
        resolve_csv(base, "outputs_exploratory/prediction_a_n56_mixed/prediction_a_ablation_winners.csv")
    )
    n60 = pd.read_csv(
        resolve_csv(base, "outputs_exploratory/prediction_a_n60_mixed/prediction_a_ablation_winners.csv")
    )
    n64 = pd.read_csv(
        resolve_csv(base, "outputs_exploratory/prediction_a_n64_mixed/prediction_a_ablation_winners.csv")
    )
    n68 = pd.read_csv(
        resolve_csv(base, "outputs_exploratory/prediction_a_n68_mixed/prediction_a_ablation_winners.csv")
    )
    n72 = pd.read_csv(
        resolve_csv(base, "outputs_exploratory/prediction_a_n72_mixed/prediction_a_ablation_winners.csv")
    )
    keep_cols = ["n", "gamma", "variant", "winner_family"]
    df = pd.concat(
        [early[keep_cols], late[keep_cols], n52[keep_cols], n56[keep_cols], n60[keep_cols], n64[keep_cols], n68[keep_cols], n72[keep_cols]],
        ignore_index=True,
    )
    df = df[df["variant"].isin(VARIANTS)].copy()
    df["winner_code"] = df["winner_family"].map(FAMILY_TO_CODE)
    return df


def build_grid(df: pd.DataFrame, variant: str, n_values: list[int], gammas: list[float]):
    subset = df[df["variant"] == variant]
    pivot = (
        subset.pivot(index="gamma", columns="n", values="winner_code")
        .reindex(index=gammas, columns=n_values)
        .astype(float)
    )
    return pivot


def main() -> None:
    base = Path(__file__).resolve().parent
    out_dir = base / "outputs_exploratory" / "prediction_a_phase_summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_winner_tables(base)
    n_values = sorted(df["n"].unique().tolist())
    gammas = sorted(df["gamma"].unique().tolist())

    cmap = ListedColormap(COLORS)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True, constrained_layout=True)

    for ax, variant in zip(axes, VARIANTS):
        grid = build_grid(df, variant, n_values, gammas)
        im = ax.imshow(
            grid.values,
            cmap=cmap,
            vmin=-0.5,
            vmax=3.5,
            aspect="auto",
            origin="lower",
        )
        ax.set_title(variant)
        ax.set_xticks(range(len(n_values)))
        ax.set_xticklabels(n_values)
        ax.set_yticks(range(len(gammas)))
        ax.set_yticklabels([f"{g:.1f}" for g in gammas])
        ax.set_xlabel("N")

    axes[0].set_ylabel("gamma")

    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", markersize=10, markerfacecolor=COLORS[code], markeredgecolor="none")
        for code in range(4)
    ]
    labels = [CODE_TO_LABEL[code] for code in range(4)]
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.04))

    png_path = out_dir / "prediction_a_winner_phase_comparison.png"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    summary_rows = []
    for variant in VARIANTS:
        counts = (
            df[df["variant"] == variant]["winner_family"]
            .value_counts()
            .to_dict()
        )
        row = {"variant": variant}
        row.update(counts)
        summary_rows.append(row)
    pd.DataFrame(summary_rows).to_csv(out_dir / "prediction_a_winner_phase_counts.csv", index=False, encoding="utf-8-sig")

    print(png_path.as_posix())


if __name__ == "__main__":
    main()
