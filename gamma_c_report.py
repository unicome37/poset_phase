from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from experiment import load_config


def estimate_gamma_crossing(
    summary_df: pd.DataFrame,
    family_a: str,
    family_b: str,
    action_mode: str,
) -> pd.DataFrame:
    rows = []
    df = summary_df[summary_df["action_mode"] == action_mode]

    for n in sorted(df["n"].unique()):
        a = (
            df[(df["n"] == n) & (df["family"] == family_a)][["gamma", "mean_score_norm"]]
            .rename(columns={"mean_score_norm": "score_a"})
            .sort_values("gamma")
        )
        b = (
            df[(df["n"] == n) & (df["family"] == family_b)][["gamma", "mean_score_norm"]]
            .rename(columns={"mean_score_norm": "score_b"})
            .sort_values("gamma")
        )
        merged = a.merge(b, on="gamma", how="inner")
        if merged.empty:
            continue

        merged["diff"] = merged["score_a"] - merged["score_b"]
        gamma_c = None
        for i in range(len(merged) - 1):
            d1 = float(merged.iloc[i]["diff"])
            d2 = float(merged.iloc[i + 1]["diff"])
            if d1 == 0.0:
                gamma_c = float(merged.iloc[i]["gamma"])
                break
            if d1 * d2 < 0:
                g1 = float(merged.iloc[i]["gamma"])
                g2 = float(merged.iloc[i + 1]["gamma"])
                gamma_c = g1 + (0.0 - d1) * (g2 - g1) / (d2 - d1)
                break

        rows.append(
            {
                "n": n,
                "action_mode": action_mode,
                "family_a": family_a,
                "family_b": family_b,
                "gamma_c_est": gamma_c,
                "diff_at_min_gamma": float(merged.iloc[0]["diff"]),
                "diff_at_max_gamma": float(merged.iloc[-1]["diff"]),
            }
        )

    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Estimate gamma crossing points from summary.csv.")
    parser.add_argument("--config", default="config_smallN_exact.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    summary_df = pd.read_csv(out_dir / config["output"]["summary_csv"])

    reports = []
    for action_mode in ["A2", "A3"]:
        reports.append(estimate_gamma_crossing(summary_df, "lorentzian_like_2d", "KR_like", action_mode))
        reports.append(estimate_gamma_crossing(summary_df, "lorentzian_like_2d", "multi_layer_random", action_mode))
        reports.append(estimate_gamma_crossing(summary_df, "lorentzian_like_2d", "transitive_percolation", action_mode))

    report_df = pd.concat(reports, ignore_index=True)
    report_path = out_dir / "gamma_c_report.csv"
    report_df.to_csv(report_path, index=False, encoding="utf-8-sig")
    print(report_df.to_string(index=False))
