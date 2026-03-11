from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def zscore_by_n(df: pd.DataFrame, value_col: str) -> pd.Series:
    out = []
    for _n, sub in df.groupby("n", sort=False):
        mean = float(sub[value_col].mean())
        std = float(sub[value_col].std(ddof=0))
        std = std if std > 1e-12 else 1.0
        out.extend(((sub[value_col] - mean) / std).tolist())
    return pd.Series(out, index=df.index)


def crossing_for_variant(df: pd.DataFrame, value_col: str, zetas: list[float], variant: str) -> pd.DataFrame:
    rows = []
    for n, sub_n in df.groupby("n"):
        scan_rows = []
        for zeta in zetas:
            fam_means = (
                sub_n.assign(score_eval=sub_n["score_A2_gamma"] + zeta * sub_n[value_col])
                .groupby("family")["score_eval"]
                .mean()
            )
            delta = float(fam_means["multi_layer_random"] - fam_means["lorentzian_like_2d"])
            scan_rows.append({"n": int(n), "zeta": float(zeta), "delta": delta})
        scan_df = pd.DataFrame(scan_rows).sort_values("zeta").reset_index(drop=True)
        cross = None
        for i in range(1, len(scan_df)):
            y0 = float(scan_df.loc[i - 1, "delta"])
            y1 = float(scan_df.loc[i, "delta"])
            if y0 < 0.0 <= y1:
                x0 = float(scan_df.loc[i - 1, "zeta"])
                x1 = float(scan_df.loc[i, "zeta"])
                frac = (-y0) / max(y1 - y0, 1e-12)
                cross = x0 + frac * (x1 - x0)
                break
        rows.append(
            {
                "variant": variant,
                "n": int(n),
                "zeta_cross": cross,
                "delta_at_min_zeta": float(scan_df.iloc[0]["delta"]),
                "delta_at_max_zeta": float(scan_df.iloc[-1]["delta"]),
            }
        )
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Robustness check for switch_zscore under blinder family-identity variants.")
    parser.add_argument("--config", default="config_pairwise_switch_robustness_check.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    input_csv = Path(config["input"]["raw_csv"])
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv)
    df = df[df["family"].isin(["lorentzian_like_2d", "multi_layer_random"])].copy()

    # Baseline: current switch_zscore.
    df["switch_zscore"] = zscore_by_n(df, "cg_family_switch_rate")

    # Blind variant 1: thresholded family retention (only asks whether switch rate is above 0.5).
    df["switch_threshold"] = (df["cg_family_switch_rate"] > 0.5).astype(float)
    df["switch_threshold_zscore"] = zscore_by_n(df, "switch_threshold")

    # Blind variant 2: use coarse-grain penalty rank surrogate without explicit centroid distance.
    # Here we retain only switch/no-switch information and discard original weighting.
    df["switch_centered"] = df.groupby("n")["cg_family_switch_rate"].transform(lambda s: s - float(s.mean()))

    zetas = [float(v) for v in config["experiment"]["zetas"]]

    reports = [
        crossing_for_variant(df, "switch_zscore", zetas, "switch_zscore"),
        crossing_for_variant(df, "switch_threshold_zscore", zetas, "switch_threshold_zscore"),
        crossing_for_variant(df, "switch_centered", zetas, "switch_centered"),
    ]
    report_df = pd.concat(reports, ignore_index=True)

    summary = (
        df.groupby(["n", "family"])
        .agg(
            mean_switch=("cg_family_switch_rate", "mean"),
            mean_switch_z=("switch_zscore", "mean"),
            mean_switch_threshold=("switch_threshold", "mean"),
            mean_switch_centered=("switch_centered", "mean"),
            mean_score_A2=("score_A2_gamma", "mean"),
            count=("score_A2_gamma", "count"),
        )
        .reset_index()
    )

    summary.to_csv(out_dir / "pairwise_switch_robustness_summary.csv", index=False, encoding="utf-8-sig")
    report_df.to_csv(out_dir / "pairwise_switch_robustness_crossings.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "pairwise_switch_robustness_summary.csv").as_posix())
    print((out_dir / "pairwise_switch_robustness_crossings.csv").as_posix())
    print()
    print(report_df.to_string(index=False))
