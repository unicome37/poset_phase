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


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot gamma_c as a function of consistency scale factor for width_height + dim_consistency."
    )
    parser.add_argument(
        "--config",
        default="config_width_height_consistency_scale_scan.yaml",
        help="Path to YAML config.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    report_path = out_dir / "width_height_consistency_scale_scan_report.csv"
    report_df = pd.read_csv(report_path)
    report_df["scale_factor"] = report_df["scale_label"].astype(float)

    # Build a compact summary focused on whether N=44 crossing survives.
    n44 = report_df[report_df["n"] == 44].copy().sort_values("scale_factor")
    threshold_rows = n44[n44["status"] == "crossing"]
    first_scale = float(threshold_rows.iloc[0]["scale_factor"]) if not threshold_rows.empty else None
    threshold_summary = pd.DataFrame(
        [
            {
                "n": 44,
                "first_crossing_scale_factor": first_scale,
                "crossing_scales": ", ".join(f"{x:.3f}" for x in threshold_rows["scale_factor"].tolist()),
            }
        ]
    )
    threshold_summary.to_csv(
        out_dir / "width_height_consistency_scale_threshold_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for n, sub in report_df.groupby("n", sort=True):
        sub = sub.sort_values("scale_factor")
        crossing = sub[sub["status"] == "crossing"]
        failed = sub[sub["status"] != "crossing"]
        if not crossing.empty:
            ax.plot(
                crossing["scale_factor"],
                crossing["gamma_c_est"],
                marker="o",
                linewidth=1.8,
                label=f"N={n}",
            )
        if not failed.empty:
            ax.scatter(
                failed["scale_factor"],
                [0.0] * len(failed),
                marker="x",
                s=42,
                alpha=0.9,
                color=ax.lines[-1].get_color() if ax.lines else None,
            )

    if first_scale is not None:
        ax.axvline(first_scale, color="#d62728", linestyle="--", linewidth=1.5)
        ax.text(first_scale + 0.02, 0.12, f"N=44 recovers at ~{first_scale:.2f}x", color="#d62728")

    ax.set_xlabel("Consistency Weight Scale Factor")
    ax.set_ylabel("Estimated gamma_c")
    ax.set_title("width_height + dim_consistency: gamma_c vs consistency weight")
    ax.grid(alpha=0.25, linestyle="--")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(out_dir / "width_height_consistency_scale_scan.png", dpi=160, bbox_inches="tight")

    print((out_dir / "width_height_consistency_scale_scan.png").as_posix())
    print((out_dir / "width_height_consistency_scale_threshold_summary.csv").as_posix())
