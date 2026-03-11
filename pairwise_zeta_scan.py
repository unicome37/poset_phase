from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def add_variant_column(raw_df: pd.DataFrame, variant: str) -> pd.DataFrame:
    df = raw_df.copy()
    if variant == "raw_cg":
        df["cg_variant"] = df["cg_mean_penalty"]
    elif variant == "zscore_cg":
        pieces = []
        for n, sub in df.groupby("n"):
            mean = float(sub["cg_mean_penalty"].mean())
            std = float(sub["cg_mean_penalty"].std(ddof=0))
            std = std if std > 1e-12 else 1.0
            sub = sub.copy()
            sub["cg_variant"] = (sub["cg_mean_penalty"] - mean) / std
            pieces.append(sub)
        df = pd.concat(pieces, ignore_index=True)
    else:
        raise ValueError(f"Unsupported variant: {variant}")
    df["variant"] = variant
    return df


def build_scan(raw_df: pd.DataFrame, zetas: list[float], variant: str) -> pd.DataFrame:
    raw_df = add_variant_column(raw_df, variant)
    rows: list[dict] = []
    for zeta in zetas:
        df = raw_df.copy()
        df["score_A2_cg"] = df["score_A2_gamma"] + zeta * df["cg_variant"]
        summary = (
            df.groupby(["n", "family"])
            .agg(
                mean_score_A2_cg=("score_A2_cg", "mean"),
                std_score_A2_cg=("score_A2_cg", "std"),
                mean_score_A2=("score_A2_gamma", "mean"),
                mean_cg_penalty=("cg_variant", "mean"),
                count=("score_A2_gamma", "count"),
            )
            .reset_index()
        )
        summary["zeta"] = zeta
        summary["variant"] = variant
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def crossing_report(scan_df: pd.DataFrame) -> pd.DataFrame:
    pivot = scan_df.pivot_table(
        index=["variant", "n", "zeta"],
        columns="family",
        values="mean_score_A2_cg",
    ).reset_index()
    pivot["delta_mlr_minus_lor2d"] = (
        pivot["multi_layer_random"] - pivot["lorentzian_like_2d"]
    )

    rows: list[dict] = []
    for (variant, n), sub in pivot.groupby(["variant", "n"]):
        sub = sub.sort_values("zeta").reset_index(drop=True)
        zeta_cross = None
        for i in range(1, len(sub)):
            y0 = float(sub.loc[i - 1, "delta_mlr_minus_lor2d"])
            y1 = float(sub.loc[i, "delta_mlr_minus_lor2d"])
            if y0 == 0.0:
                zeta_cross = float(sub.loc[i - 1, "zeta"])
                break
            if y0 < 0.0 <= y1:
                x0 = float(sub.loc[i - 1, "zeta"])
                x1 = float(sub.loc[i, "zeta"])
                frac = (-y0) / max(y1 - y0, 1e-12)
                zeta_cross = x0 + frac * (x1 - x0)
                break
        rows.append(
            {
                "variant": str(variant),
                "n": int(n),
                "zeta_cross": zeta_cross,
                "delta_at_min_zeta": float(sub.iloc[0]["delta_mlr_minus_lor2d"]),
                "delta_at_max_zeta": float(sub.iloc[-1]["delta_mlr_minus_lor2d"]),
            }
        )
    return pd.DataFrame(rows)


def plot_scan(scan_df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(
        nrows=len(sorted(scan_df["n"].unique())),
        ncols=1,
        figsize=(8, 3.5 * len(sorted(scan_df["n"].unique()))),
        squeeze=False,
    )
    for ax, (n, sub) in zip(axes[:, 0], scan_df.groupby("n")):
        for (variant, family), fam_df in sub.groupby(["variant", "family"]):
            fam_df = fam_df.sort_values("zeta")
            ax.plot(
                fam_df["zeta"],
                fam_df["mean_score_A2_cg"],
                marker="o",
                label=f"{family} [{variant}]",
            )
        ax.set_title(f"Controlled-window A2 + zeta*I_cg, N={n}")
        ax.set_xlabel("zeta")
        ax.set_ylabel("mean score")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Zeta scan for Lor2D vs multi-layer under matched compressibility window.")
    parser.add_argument("--config", default="config_pairwise_zeta_scan.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    input_csv = Path(config["input"]["raw_csv"])
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_df = pd.read_csv(input_csv)
    raw_df = raw_df[raw_df["family"].isin(["lorentzian_like_2d", "multi_layer_random"])].copy()
    zetas = [float(v) for v in config["experiment"]["zetas"]]
    variants = [str(v) for v in config["experiment"].get("variants", ["raw_cg"])]

    scan_parts = [build_scan(raw_df, zetas, variant=variant) for variant in variants]
    scan_df = pd.concat(scan_parts, ignore_index=True)
    cross_df = crossing_report(scan_df)

    scan_df.to_csv(out_dir / "pairwise_zeta_scan.csv", index=False, encoding="utf-8-sig")
    cross_df.to_csv(out_dir / "pairwise_zeta_crossings.csv", index=False, encoding="utf-8-sig")
    plot_scan(scan_df, out_dir / "pairwise_zeta_scan.png")

    print((out_dir / "pairwise_zeta_scan.csv").as_posix())
    print((out_dir / "pairwise_zeta_crossings.csv").as_posix())
    print((out_dir / "pairwise_zeta_scan.png").as_posix())
    print()
    print(cross_df.to_string(index=False))
