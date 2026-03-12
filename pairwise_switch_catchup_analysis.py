from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\Kiro")


def load_switch_crossings() -> pd.DataFrame:
    rows = []

    early = pd.read_csv(
        ROOT / "outputs_exploratory" / "pairwise_switch_only_scan" / "pairwise_switch_crossings.csv"
    )
    early = early[early["variant"] == "switch_zscore"].copy()
    rows.append(early[["n", "zeta_cross"]])

    late_44 = pd.read_csv(
        ROOT / "outputs_exploratory" / "frozen_exploratory_submodel" / "scan" / "pairwise_switch_crossings.csv"
    )
    late_44 = late_44[late_44["variant"] == "switch_zscore"].copy()
    rows.append(late_44[["n", "zeta_cross"]])

    late_48 = pd.read_csv(
        ROOT / "outputs_exploratory" / "frozen_exploratory_submodel_48" / "scan" / "pairwise_switch_crossings.csv"
    )
    late_48 = late_48[late_48["variant"] == "switch_zscore"].copy()
    rows.append(late_48[["n", "zeta_cross"]])

    return pd.concat(rows, ignore_index=True).sort_values("n").reset_index(drop=True)


def main() -> None:
    out_dir = ROOT / "outputs_exploratory" / "pairwise_switch_catchup_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    crossings = load_switch_crossings()
    matched = pd.read_csv(
        ROOT / "outputs_exploratory" / "mlr_survivor_matched_lor2d" / "mlr_survivor_matched_summary.csv"
    )

    merged = crossings.merge(matched, on="n", how="inner").sort_values("n").reset_index(drop=True)
    merged["abs_score_gap"] = -merged["mean_score_A2_delta"]
    merged["zeta_per_logH"] = merged["zeta_cross"] / merged["mean_log_H_delta"]
    merged["logH_per_switch"] = merged["mean_log_H_delta"] / merged["mean_switch_delta"]
    merged["score_gap_per_switch"] = merged["abs_score_gap"] / merged["mean_switch_delta"]
    merged["switch_per_logH"] = merged["mean_switch_delta"] / merged["mean_log_H_delta"]

    summary = merged[
        [
            "n",
            "zeta_cross",
            "mean_log_H_delta",
            "abs_score_gap",
            "mean_switch_delta",
            "mean_geo_total_delta",
            "zeta_per_logH",
            "logH_per_switch",
            "score_gap_per_switch",
            "switch_per_logH",
            "count",
        ]
    ].copy()

    # Simple trend diagnostics.
    trend = []
    for col in ["zeta_cross", "mean_log_H_delta", "abs_score_gap", "mean_switch_delta", "zeta_per_logH"]:
        s = summary[["n", col]].dropna()
        if len(s) >= 2:
            slope = (float(s.iloc[-1][col]) - float(s.iloc[0][col])) / (float(s.iloc[-1]["n"]) - float(s.iloc[0]["n"]))
        else:
            slope = float("nan")
        trend.append(
            {
                "metric": col,
                "start_n": int(s.iloc[0]["n"]) if len(s) else None,
                "end_n": int(s.iloc[-1]["n"]) if len(s) else None,
                "start_value": float(s.iloc[0][col]) if len(s) else None,
                "end_value": float(s.iloc[-1][col]) if len(s) else None,
                "secant_slope_per_N": slope,
            }
        )
    trend_df = pd.DataFrame(trend)

    summary.to_csv(out_dir / "pairwise_switch_catchup_summary.csv", index=False, encoding="utf-8-sig")
    trend_df.to_csv(out_dir / "pairwise_switch_catchup_trend.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "pairwise_switch_catchup_summary.csv").as_posix())
    print((out_dir / "pairwise_switch_catchup_trend.csv").as_posix())
    print()
    print(summary.to_string(index=False))
    print()
    print(trend_df.to_string(index=False))


if __name__ == "__main__":
    main()
