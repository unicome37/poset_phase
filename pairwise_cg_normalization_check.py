from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def zeta_crossing(df: pd.DataFrame, score_col: str, zetas: list[float]) -> pd.DataFrame:
    rows = []
    for n, sub_n in df.groupby("n"):
        pivot = sub_n.pivot(index="family", columns="seed", values=[score_col, "cg_mean_penalty"])
        # build mean curve directly from sample table
        scan_rows = []
        for zeta in zetas:
            fam_means = (
                sub_n.assign(score_eval=sub_n[score_col] + zeta * sub_n["cg_variant"])
                .groupby("family")["score_eval"]
                .mean()
            )
            scan_rows.append(
                {
                    "n": int(n),
                    "zeta": float(zeta),
                    "lor2d": float(fam_means["lorentzian_like_2d"]),
                    "mlr": float(fam_means["multi_layer_random"]),
                    "delta_mlr_minus_lor2d": float(fam_means["multi_layer_random"] - fam_means["lorentzian_like_2d"]),
                }
            )
        scan_df = pd.DataFrame(scan_rows)
        cross = None
        for i in range(1, len(scan_df)):
            y0 = float(scan_df.loc[i - 1, "delta_mlr_minus_lor2d"])
            y1 = float(scan_df.loc[i, "delta_mlr_minus_lor2d"])
            if y0 < 0.0 <= y1:
                x0 = float(scan_df.loc[i - 1, "zeta"])
                x1 = float(scan_df.loc[i, "zeta"])
                frac = (-y0) / max(y1 - y0, 1e-12)
                cross = x0 + frac * (x1 - x0)
                break
        rows.append(
            {
                "n": int(n),
                "variant": str(sub_n["variant"].iloc[0]),
                "delta_at_min_zeta": float(scan_df.iloc[0]["delta_mlr_minus_lor2d"]),
                "delta_at_max_zeta": float(scan_df.iloc[-1]["delta_mlr_minus_lor2d"]),
                "zeta_cross": cross,
            }
        )
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check CG penalty decomposition and normalization.")
    parser.add_argument("--config", default="config_pairwise_cg_normalization_check.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    input_csv = Path(config["input"]["raw_csv"])
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(input_csv)
    raw = raw[raw["family"].isin(["lorentzian_like_2d", "multi_layer_random"])].copy()
    raw["cg_drift_component"] = 10.0 * raw["cg_mean_self_drift"]
    raw["cg_switch_component"] = 1.5 * raw["cg_family_switch_rate"]

    # Build three variants: raw cg, drift only, z-scored cg within n.
    variant_rows = []
    for variant, source_col in [
        ("raw_cg", "cg_mean_penalty"),
        ("drift_only", "cg_drift_component"),
        ("switch_only", "cg_switch_component"),
    ]:
        df = raw.copy()
        df["variant"] = variant
        df["cg_variant"] = df[source_col]
        variant_rows.append(df)

    z_rows = []
    for n, sub in raw.groupby("n"):
        mean = float(sub["cg_mean_penalty"].mean())
        std = float(sub["cg_mean_penalty"].std(ddof=0))
        std = std if std > 1e-12 else 1.0
        sub = sub.copy()
        sub["variant"] = "zscore_cg"
        sub["cg_variant"] = (sub["cg_mean_penalty"] - mean) / std
        z_rows.append(sub)
    variant_rows.append(pd.concat(z_rows, ignore_index=True))

    all_df = pd.concat(variant_rows, ignore_index=True)
    zetas = [float(v) for v in config["experiment"]["zetas"]]

    reports = []
    for variant, sub in all_df.groupby("variant"):
        reports.append(zeta_crossing(sub, "score_A2_gamma", zetas))
    report_df = pd.concat(reports, ignore_index=True)

    summary = (
        all_df.groupby(["variant", "n", "family"])
        .agg(
            mean_cg_variant=("cg_variant", "mean"),
            mean_drift_component=("cg_drift_component", "mean"),
            mean_switch_component=("cg_switch_component", "mean"),
            mean_score_A2=("score_A2_gamma", "mean"),
            count=("score_A2_gamma", "count"),
        )
        .reset_index()
    )

    summary.to_csv(out_dir / "pairwise_cg_component_summary.csv", index=False, encoding="utf-8-sig")
    report_df.to_csv(out_dir / "pairwise_cg_normalization_crossings.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "pairwise_cg_component_summary.csv").as_posix())
    print((out_dir / "pairwise_cg_normalization_crossings.csv").as_posix())
    print()
    print(report_df.to_string(index=False))
