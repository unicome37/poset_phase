from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


BLUE = "#1f77b4"
RED = "#d62728"
GRAY = "#666666"
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2"]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize the near-wall mixed Lor2D-vs-KR frontier for Prediction B."
    )
    parser.add_argument(
        "--ns",
        nargs="+",
        type=int,
        default=[52, 56, 60, 64, 68, 72],
        help="N values to include.",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs_exploratory/prediction_b_nearwall_mixed_frontier",
        help="Output directory relative to the poset_phase root.",
    )
    return parser


def resolve_summary_csv(poset_root: Path, n: int) -> Path | None:
    local = (
        poset_root
        / "outputs_exploratory"
        / f"prediction_a_n{n}_mixed"
        / "prediction_a_ablation_summary.csv"
    )
    if local.exists():
        return local

    external = (
        poset_root.parent.parent
        / "outputs_exploratory"
        / f"prediction_a_n{n}_mixed"
        / "prediction_a_ablation_summary.csv"
    )
    if external.exists():
        return external

    return None


def load_combined(poset_root: Path, ns: list[int]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for n in ns:
        path = resolve_summary_csv(poset_root, n)
        if path is None:
            continue
        df = pd.read_csv(path)
        df["source_csv"] = path.as_posix()
        frames.append(df)
    if not frames:
        raise FileNotFoundError("No mixed near-wall summary CSVs were found for the requested N values.")
    return pd.concat(frames, ignore_index=True)


def build_frontier_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sub = df[
        (df["variant"] == "A2_full")
        & (df["family"].isin(["lorentzian_like_2d", "KR_like"]))
    ].copy()
    if sub.empty:
        raise ValueError("No A2_full Lor2D/KR rows found in the input summaries.")
    source_map = sub.groupby("n")["source_csv"].first().to_dict()

    pivot = (
        sub.pivot_table(
            index=["n", "gamma"],
            columns="family",
            values=["mean_score", "std_score", "mean_log_H", "mean_penalty", "count"],
            aggfunc="first",
        )
        .sort_index()
    )
    pivot.columns = [f"{metric}_{family}" for metric, family in pivot.columns]
    combined = pivot.reset_index()
    combined["delta_kr_minus_lor2d"] = (
        combined["mean_score_KR_like"] - combined["mean_score_lorentzian_like_2d"]
    )
    combined["abs_gap"] = combined["delta_kr_minus_lor2d"].abs()
    combined["kr_wins"] = combined["delta_kr_minus_lor2d"] < 0
    combined["lor2d_wins"] = combined["delta_kr_minus_lor2d"] > 0

    summaries: list[dict[str, object]] = []
    for n, group in combined.groupby("n", sort=True):
        group = group.sort_values("gamma").reset_index(drop=True)
        min_row = group.loc[group["abs_gap"].idxmin()]
        max_row = group.loc[group["gamma"].idxmax()]
        status = "crossing_observed" if (group["delta_kr_minus_lor2d"] > 0).any() else "no_cross_kr_favored"
        summaries.append(
            {
                "n": int(n),
                "gamma_min": float(group["gamma"].min()),
                "gamma_max": float(group["gamma"].max()),
                "status": status,
                "delta_at_gamma_min": float(group.loc[group["gamma"].idxmin(), "delta_kr_minus_lor2d"]),
                "delta_at_gamma_max": float(max_row["delta_kr_minus_lor2d"]),
                "abs_gap_at_gamma_max": float(abs(max_row["delta_kr_minus_lor2d"])),
                "gamma_at_min_abs_gap": float(min_row["gamma"]),
                "min_abs_gap": float(min_row["abs_gap"]),
                "kr_mean_score_at_gamma_max": float(max_row["mean_score_KR_like"]),
                "lor2d_mean_score_at_gamma_max": float(max_row["mean_score_lorentzian_like_2d"]),
                "kr_mean_penalty": float(max_row["mean_penalty_KR_like"]),
                "lor2d_mean_penalty": float(max_row["mean_penalty_lorentzian_like_2d"]),
                "kr_mean_log_H": float(max_row["mean_log_H_KR_like"]),
                "lor2d_mean_log_H": float(max_row["mean_log_H_lorentzian_like_2d"]),
                "source_csv": str(source_map.get(n, "")),
            }
        )
    summary_df = pd.DataFrame(summaries).sort_values("n").reset_index(drop=True)
    return combined, summary_df


def plot_gap_frontier(combined: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.3))

    for color, (n, group) in zip(PALETTE, combined.groupby("n", sort=True)):
        group = group.sort_values("gamma")
        ax.plot(
            group["gamma"],
            group["delta_kr_minus_lor2d"],
            marker="o",
            linewidth=1.6,
            markersize=4.5,
            color=color,
            label=f"N={int(n)}",
        )
        end = group.iloc[-1]
        ax.annotate(
            f"{int(n)}",
            xy=(end["gamma"], end["delta_kr_minus_lor2d"]),
            xytext=(5, 0),
            textcoords="offset points",
            fontsize=7,
            color=color,
            va="center",
        )

    ax.axhline(0.0, color=GRAY, linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_xlabel(r"Geometric coupling $\gamma$")
    ax.set_ylabel(r"$\Delta A_{\mathrm{KR-L2D}} = A_{\mathrm{KR}} - A_{\mathrm{Lor2D}}$")
    ax.set_title("Prediction B near-wall mixed frontier under $A_2^{full}$")
    ax.grid(True, alpha=0.22, linewidth=0.5)
    ax.legend(loc="upper right", ncol=2, fontsize=7.5, frameon=True)
    fig.tight_layout()

    fig.savefig(out_dir / "prediction_b_nearwall_gap_frontier.png", dpi=180, bbox_inches="tight")
    fig.savefig(out_dir / "prediction_b_nearwall_gap_frontier.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(summary_df: pd.DataFrame, out_dir: Path) -> None:
    lines: list[str] = []
    lines.append("# Prediction B Near-Wall Mixed Frontier")
    lines.append("")
    lines.append("> Scope: exploratory mixed continuation beyond the confirmatory exact window.")
    lines.append("> Families: `lorentzian_like_2d` vs `KR_like`.")
    lines.append("> Action: `A2_full`.")
    lines.append("")

    if (summary_df["status"] == "crossing_observed").any():
        headline = "A crossing is observed in at least one scanned near-wall slice."
    else:
        headline = (
            "No crossing is observed up to the scanned `gamma_max`, but the residual "
            "KR-vs-Lor2D gap stays finite and contracts substantially toward the right boundary."
        )
    lines.append(headline)
    lines.append("")
    lines.append("| N | status | ΔA at γ_max | |ΔA| min | γ at min gap |")
    lines.append("|---|---|---:|---:|---:|")
    for row in summary_df.itertuples(index=False):
        lines.append(
            f"| {row.n} | {row.status} | {row.delta_at_gamma_max:.3f} | {row.min_abs_gap:.3f} | {row.gamma_at_min_abs_gap:.1f} |"
        )
    lines.append("")
    lines.append("Interpretation:")
    lines.append("- Negative `ΔA_KR-L2D` means KR still wins; values near zero indicate near-degeneracy.")
    lines.append("- These mixed results extend the competition frontier, but they do not by themselves establish persistence of a bounded `gamma_c` beyond the exact `N<=44` confirmatory window.")
    lines.append("- The practical value of this report is to separate a stronger claim from a weaker one: the exact bounded-window result remains confirmatory, while the large-`N` mixed frontier shows that the competition remains active rather than collapsing.")
    lines.append("")

    (out_dir / "prediction_b_nearwall_gap_frontier.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    args = build_arg_parser().parse_args()
    poset_root = Path(__file__).resolve().parent
    out_dir = (poset_root / args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = load_combined(poset_root, args.ns)
    combined, summary_df = build_frontier_tables(raw)

    combined.to_csv(out_dir / "prediction_b_nearwall_mixed_combined.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "prediction_b_nearwall_gap_summary.csv", index=False, encoding="utf-8-sig")
    plot_gap_frontier(combined, out_dir)
    write_report(summary_df, out_dir)

    print((out_dir / "prediction_b_nearwall_mixed_combined.csv").as_posix())
    print((out_dir / "prediction_b_nearwall_gap_summary.csv").as_posix())
    print((out_dir / "prediction_b_nearwall_gap_frontier.png").as_posix())
    print((out_dir / "prediction_b_nearwall_gap_frontier.md").as_posix())


if __name__ == "__main__":
    main()
