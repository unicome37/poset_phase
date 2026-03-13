from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_switch_std_map(config: dict) -> dict[int, float]:
    rows: list[dict] = []
    for item in config["input"]["switch_duel_raw_files"]:
        n = int(item["n"])
        csv_path = Path(item["path"])
        df = pd.read_csv(csv_path)
        df = df[df["n"] == n].copy()
        df = df[df["family"].isin(["lorentzian_like_2d", "multi_layer_random"])].copy()
        std = float(df["cg_family_switch_rate"].std(ddof=0))
        if std <= 1e-12:
            std = 1.0
        rows.append({"n": n, "switch_std": std, "raw_count": len(df)})
    return {int(row["n"]): float(row["switch_std"]) for row in rows}


def merge_inputs(config: dict) -> pd.DataFrame:
    depth_df = pd.read_csv(Path(config["input"]["matched_pair_raw_csv"]))
    switch_df = pd.read_csv(Path(config["input"]["matched_pair_switch_csv"]))
    keep_cols = [
        "n",
        "mlr_seed",
        "lor_seed",
        "cg_family_switch_rate_delta_mlr_minus_lor2d",
    ]
    merged = depth_df.merge(
        switch_df[keep_cols],
        on=["n", "mlr_seed", "lor_seed"],
        how="inner",
        validate="one_to_one",
    )
    return merged


def add_scaled_features(df: pd.DataFrame, feature_specs: list[dict], switch_std_map: dict[int, float]) -> pd.DataFrame:
    out = df.copy()
    out["switch_std_ref"] = out["n"].map(switch_std_map)
    out["switch_delta_z_equiv"] = (
        out["cg_family_switch_rate_delta_mlr_minus_lor2d"] / out["switch_std_ref"]
    )

    for spec in feature_specs:
        target = spec["name"]
        if "source_col" in spec:
            source = spec["source_col"]
            sign = float(spec.get("sign", 1.0))
            base = sign * out[source]
        else:
            base = pd.Series(0.0, index=out.index, dtype=float)
            for component in spec["components"]:
                comp_name = component["feature"]
                weight = float(component.get("weight", 1.0))
                base = base + weight * out[comp_name]
        std = float(base.std(ddof=0))
        std = std if std > 1e-12 else 1.0
        out[target] = base / std

    return out


def summarize_crossings(df: pd.DataFrame, feature_specs: list[dict], etas: list[float]) -> pd.DataFrame:
    baseline = (
        df.groupby("n")
        .agg(
            pair_count=("n", "size"),
            mean_score_delta=("score_A2_gamma_delta_mlr_minus_lor2d", "mean"),
            mean_switch_delta_raw=("cg_family_switch_rate_delta_mlr_minus_lor2d", "mean"),
            mean_switch_delta_z=("switch_delta_z_equiv", "mean"),
        )
        .reset_index()
    )
    baseline["baseline_zeta_cross"] = -baseline["mean_score_delta"] / baseline["mean_switch_delta_z"]

    rows: list[dict] = []
    for spec in feature_specs:
        feature = spec["name"]
        for eta in etas:
            tmp = (
                df.assign(
                    adjusted_score_delta=(
                        df["score_A2_gamma_delta_mlr_minus_lor2d"] + eta * df[feature]
                    )
                )
                .groupby("n")
                .agg(
                    pair_count=("n", "size"),
                    mean_adjusted_score_delta=("adjusted_score_delta", "mean"),
                    mean_depth_term=(feature, "mean"),
                    mean_switch_delta_z=("switch_delta_z_equiv", "mean"),
                )
                .reset_index()
            )
            tmp["feature"] = feature
            tmp["eta_depth"] = eta
            tmp["zeta_cross"] = -tmp["mean_adjusted_score_delta"] / tmp["mean_switch_delta_z"]
            rows.extend(tmp.to_dict("records"))

    summary = pd.DataFrame(rows).merge(
        baseline[["n", "baseline_zeta_cross", "mean_score_delta", "mean_switch_delta_raw", "mean_switch_delta_z"]],
        on="n",
        how="left",
    )
    summary["zeta_reduction_abs"] = summary["baseline_zeta_cross"] - summary["zeta_cross"]
    summary["zeta_reduction_pct"] = summary["zeta_reduction_abs"] / summary["baseline_zeta_cross"]
    return summary.sort_values(["feature", "eta_depth", "n"]).reset_index(drop=True)


def build_best_eta_table(summary: pd.DataFrame) -> pd.DataFrame:
    best = (
        summary.sort_values(["feature", "n", "zeta_cross", "eta_depth"])
        .groupby(["feature", "n"], as_index=False)
        .first()
    )
    return best[
        [
            "feature",
            "n",
            "eta_depth",
            "baseline_zeta_cross",
            "zeta_cross",
            "zeta_reduction_abs",
            "zeta_reduction_pct",
            "mean_depth_term",
            "pair_count",
        ]
    ].copy()


def build_feature_overview(best_df: pd.DataFrame) -> pd.DataFrame:
    overview = (
        best_df.groupby("feature")
        .agg(
            mean_baseline_zeta_cross=("baseline_zeta_cross", "mean"),
            mean_best_zeta_cross=("zeta_cross", "mean"),
            mean_reduction_abs=("zeta_reduction_abs", "mean"),
            mean_reduction_pct=("zeta_reduction_pct", "mean"),
            min_best_zeta_cross=("zeta_cross", "min"),
            max_best_zeta_cross=("zeta_cross", "max"),
        )
        .reset_index()
    )
    return overview.sort_values("mean_reduction_pct", ascending=False).reset_index(drop=True)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Test whether Prediction C depth terms reduce the switch_zscore crossing weight."
    )
    parser.add_argument(
        "--config",
        default="config_prediction_c_switch_enhancement_scan.yaml",
        help="Path to YAML config.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)

    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    feature_specs = list(config["experiment"]["features"])
    etas = [float(v) for v in config["experiment"]["eta_depth_grid"]]

    df = merge_inputs(config)
    switch_std_map = load_switch_std_map(config)
    df = add_scaled_features(df, feature_specs, switch_std_map)

    summary = summarize_crossings(df, feature_specs, etas)
    best = build_best_eta_table(summary)
    overview = build_feature_overview(best)

    summary.to_csv(
        out_dir / "prediction_c_switch_enhancement_scan.csv",
        index=False,
        encoding="utf-8-sig",
    )
    best.to_csv(
        out_dir / "prediction_c_switch_enhancement_best.csv",
        index=False,
        encoding="utf-8-sig",
    )
    overview.to_csv(
        out_dir / "prediction_c_switch_enhancement_overview.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print((out_dir / "prediction_c_switch_enhancement_scan.csv").as_posix())
    print((out_dir / "prediction_c_switch_enhancement_best.csv").as_posix())
    print((out_dir / "prediction_c_switch_enhancement_overview.csv").as_posix())
    print()
    print(best.to_string(index=False))
    print()
    print(overview.to_string(index=False))
