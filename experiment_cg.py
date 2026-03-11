from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from action import action_value, get_action_penalty
from coarse_grain import coarse_grain_delete_nodes
from experiment import FAMILIES, load_config
from runtime_utils import estimate_entropy
from stability import (
    coarse_grain_penalty,
    consistency_reference_windows,
    family_centroids,
    global_consistency_penalty_with_reference,
    nearest_family,
    self_drift,
    signature_dict,
)


def build_original_samples(
    *,
    n_values: tuple[int, ...],
    families: tuple[str, ...],
    samples_per_family: int,
    sis_runs: int,
    exact_threshold: int,
    gammas: tuple[float, ...],
    beta: float,
    action_modes: tuple[str, ...],
) -> tuple[pd.DataFrame, dict[tuple[int, str, int], object]]:
    rows: list[dict] = []
    posets: dict[tuple[int, str, int], object] = {}

    for n in n_values:
        for family in families:
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 1000 * n + sample_id
                poset = generator(n=n, seed=seed)
                posets[(n, family, sample_id)] = poset
                log_h, entropy_method = estimate_entropy(
                    poset,
                    sis_runs=sis_runs,
                    seed=seed,
                    exact_threshold=exact_threshold,
                )
                sig = signature_dict(poset)
                for gamma in gammas:
                    for action_mode in action_modes:
                        penalty_local = get_action_penalty(poset, action_mode)
                        score_local = action_value(
                            log_extensions=log_h,
                            penalty=penalty_local,
                            beta=beta,
                            gamma=gamma,
                        )
                        row = {
                            "stage": "original",
                            "n": n,
                            "n_cg": n,
                            "family": family,
                            "sample_id": sample_id,
                            "keep_ratio": 1.0,
                            "cg_repeat": 0,
                            "gamma": gamma,
                            "beta": beta,
                            "action_mode": action_mode,
                            "entropy_method": entropy_method,
                            "log_H": log_h,
                            "penalty_local": penalty_local,
                            "score_local": score_local,
                            "penalty_global_consistency": 0.0,
                            "score_augmented": score_local,
                            "score_augmented_full": score_local,
                        }
                        row.update(sig)
                        rows.append(row)

    return pd.DataFrame(rows), posets


def add_family_ranks(df: pd.DataFrame, score_col: str, rank_col: str) -> pd.DataFrame:
    family_means = (
        df.groupby(["n", "gamma", "action_mode", "family"])[score_col]
        .mean()
        .reset_index()
    )
    family_means[rank_col] = (
        family_means.groupby(["n", "gamma", "action_mode"])[score_col]
        .rank(method="dense", ascending=True)
        .astype(int)
    )
    return df.merge(
        family_means[["n", "gamma", "action_mode", "family", rank_col]],
        on=["n", "gamma", "action_mode", "family"],
        how="left",
    )


def build_coarse_grained_samples(
    *,
    original_df: pd.DataFrame,
    posets: dict[tuple[int, str, int], object],
    keep_ratios: tuple[float, ...],
    cg_repeats: int,
    sis_runs: int,
    exact_threshold: int,
    zeta: float,
    eta: float,
    gc_reference_windows: dict[int, dict[str, tuple[float, float]]] | None,
) -> pd.DataFrame:
    base_sig_df = (
        original_df[
            ["n", "family", "sample_id", "gamma", "action_mode", "base_family_rank"]
            + [c for c in original_df.columns if c.startswith("sig_")]
        ]
        .drop_duplicates(subset=["n", "family", "sample_id"])
    )
    sig_lookup = {
        (int(row.n), str(row.family), int(row.sample_id)): {col: float(getattr(row, col)) for col in base_sig_df.columns if col.startswith("sig_")}
        for row in base_sig_df.itertuples(index=False)
    }
    centroids = family_centroids(base_sig_df[["n", "family"] + [c for c in base_sig_df.columns if c.startswith("sig_")]])

    rows: list[dict] = []

    unique_original = original_df[["n", "family", "sample_id", "gamma", "beta", "action_mode", "base_family_rank"]].drop_duplicates()
    for row in unique_original.itertuples(index=False):
        key = (int(row.n), str(row.family), int(row.sample_id))
        poset = posets[key]
        before_sig = sig_lookup[key]
        for keep_ratio in keep_ratios:
            for cg_repeat in range(cg_repeats):
                seed = 500_000 + row.n * 1000 + row.sample_id * 100 + int(round(keep_ratio * 100)) * 10 + cg_repeat
                cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=keep_ratio, seed=seed)
                log_h, entropy_method = estimate_entropy(
                    cg_poset,
                    sis_runs=sis_runs,
                    seed=seed,
                    exact_threshold=exact_threshold,
                )
                sig_after = signature_dict(cg_poset)
                drift_value = self_drift(before_sig, sig_after)
                nearest, nearest_dist = nearest_family(sig_after, centroids, n_value=int(row.n))
                family_switch = 0.0 if nearest == row.family else 1.0
                rank_after = int(
                    original_df[
                        (original_df["n"] == row.n)
                        & (original_df["gamma"] == row.gamma)
                        & (original_df["action_mode"] == row.action_mode)
                        & (original_df["family"] == nearest)
                    ]["base_family_rank"].iloc[0]
                )
                rank_shift = abs(int(row.base_family_rank) - rank_after)
                penalty_cg = coarse_grain_penalty(
                    drift_value=drift_value,
                    family_switch_penalty=family_switch,
                    rank_shift_penalty=rank_shift,
                )
                penalty_global, gc_components = global_consistency_penalty_with_reference(
                    cg_poset,
                    reference_windows=gc_reference_windows,
                )
                penalty_local = get_action_penalty(cg_poset, str(row.action_mode))
                score_local = action_value(
                    log_extensions=log_h,
                    penalty=penalty_local,
                    beta=float(row.beta),
                    gamma=float(row.gamma),
                )
                score_augmented = score_local + zeta * penalty_cg
                score_augmented_full = score_augmented + eta * penalty_global
                out_row = {
                    "stage": "coarse_grained",
                    "n": int(row.n),
                    "n_cg": int(cg_poset.n),
                    "family": str(row.family),
                    "sample_id": int(row.sample_id),
                    "keep_ratio": float(keep_ratio),
                    "cg_repeat": int(cg_repeat),
                    "gamma": float(row.gamma),
                    "beta": float(row.beta),
                    "action_mode": str(row.action_mode),
                    "entropy_method": entropy_method,
                    "log_H": log_h,
                    "penalty_local": penalty_local,
                    "score_local": score_local,
                    "penalty_cg": penalty_cg,
                    "penalty_global_consistency": penalty_global,
                    "score_augmented": score_augmented,
                    "score_augmented_full": score_augmented_full,
                    "self_drift": drift_value,
                    "nearest_family": nearest,
                    "nearest_family_dist": nearest_dist,
                    "family_switch_penalty": family_switch,
                    "rank_before": int(row.base_family_rank),
                    "rank_after": int(rank_after),
                    "rank_shift_penalty": float(rank_shift),
                }
                out_row.update(sig_after)
                out_row.update(gc_components)
                rows.append(out_row)

    return pd.DataFrame(rows)


def build_reference_windows_from_reference_family(
    *,
    posets: dict[tuple[int, str, int], object],
    reference_family: str,
    keep_ratios: tuple[float, ...],
    cg_repeats: int,
    max_pairs: int,
    sigma_scale: float,
) -> dict[int, dict[str, tuple[float, float]]]:
    reference_posets_by_ncg: dict[int, list[object]] = {}
    for (n_value, family, sample_id), poset in posets.items():
        if family != reference_family:
            continue
        reference_posets_by_ncg.setdefault(int(poset.n), []).append(poset)
        for keep_ratio in keep_ratios:
            for cg_repeat in range(cg_repeats):
                seed = 900_000 + n_value * 1000 + sample_id * 100 + int(round(keep_ratio * 100)) * 10 + cg_repeat
                cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=keep_ratio, seed=seed)
                reference_posets_by_ncg.setdefault(int(cg_poset.n), []).append(cg_poset)
    return consistency_reference_windows(
        reference_posets_by_ncg,
        max_pairs=max_pairs,
        sigma_scale=sigma_scale,
    )


def summarize_cg(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["n", "n_cg", "keep_ratio", "gamma", "action_mode", "family"])
        .agg(
            mean_score_local=("score_local", "mean"),
            mean_score_augmented=("score_augmented", "mean"),
            mean_score_augmented_full=("score_augmented_full", "mean"),
            mean_penalty_cg=("penalty_cg", "mean"),
            mean_penalty_global_consistency=("penalty_global_consistency", "mean"),
            mean_self_drift=("self_drift", "mean"),
            family_switch_rate=("family_switch_penalty", "mean"),
            mean_rank_shift=("rank_shift_penalty", "mean"),
            mean_gc_var_comp=("gc_var_comp", "mean"),
            mean_gc_var_d_eff=("gc_var_d_eff", "mean"),
            count=("score_augmented", "count"),
        )
        .reset_index()
        .sort_values(["n", "keep_ratio", "action_mode", "gamma", "mean_score_augmented_full"])
    )


def summarize_rank_shift(summary_df: pd.DataFrame) -> pd.DataFrame:
    base = summary_df.copy()
    base["rank_local"] = (
        base.groupby(["n", "n_cg", "keep_ratio", "gamma", "action_mode"])["mean_score_local"]
        .rank(method="dense", ascending=True)
        .astype(int)
    )
    base["rank_augmented_cg"] = (
        base.groupby(["n", "n_cg", "keep_ratio", "gamma", "action_mode"])["mean_score_augmented"]
        .rank(method="dense", ascending=True)
        .astype(int)
    )
    base["rank_augmented_full"] = (
        base.groupby(["n", "n_cg", "keep_ratio", "gamma", "action_mode"])["mean_score_augmented_full"]
        .rank(method="dense", ascending=True)
        .astype(int)
    )
    base["rank_delta_due_to_cg_penalty"] = base["rank_augmented_cg"] - base["rank_local"]
    base["rank_delta_due_to_full_penalty"] = base["rank_augmented_full"] - base["rank_local"]
    return base


def write_outputs(original_df: pd.DataFrame, cg_df: pd.DataFrame, summary_df: pd.DataFrame, rank_df: pd.DataFrame, config: dict) -> None:
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)
    original_df.to_csv(out_dir / "cg_original_samples.csv", index=False, encoding="utf-8-sig")
    cg_df.to_csv(out_dir / "cg_raw_samples.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "cg_summary.csv", index=False, encoding="utf-8-sig")
    rank_df.to_csv(out_dir / "cg_rank_summary.csv", index=False, encoding="utf-8-sig")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run coarse-graining stability pilot.")
    parser.add_argument("--config", default="config_cg_smoke.yaml", help="Path to YAML config file.")
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    config = load_config(args.config)
    exp_cfg = config["experiment"]
    cg_cfg = config["coarse_grain"]

    families = tuple(exp_cfg.get("families", list(FAMILIES.keys())))
    original_df, posets = build_original_samples(
        n_values=tuple(exp_cfg["n_values"]),
        families=families,
        samples_per_family=int(exp_cfg["samples_per_family"]),
        sis_runs=int(exp_cfg["sis_runs"]),
        exact_threshold=int(exp_cfg.get("exact_threshold", 0)),
        gammas=tuple(exp_cfg["gammas"]),
        beta=float(exp_cfg["beta"]),
        action_modes=tuple(exp_cfg.get("action_modes", ["A2"])),
    )
    original_df = add_family_ranks(original_df, score_col="score_local", rank_col="base_family_rank")
    reference_family = str(cg_cfg.get("reference_family", "lorentzian_like_2d"))
    gc_reference_windows = build_reference_windows_from_reference_family(
        posets=posets,
        reference_family=reference_family,
        keep_ratios=tuple(cg_cfg["keep_ratios"]),
        cg_repeats=int(cg_cfg["repeats"]),
        max_pairs=int(cg_cfg.get("gc_reference_max_pairs", 32)),
        sigma_scale=float(cg_cfg.get("gc_reference_sigma_scale", 1.0)),
    )

    cg_df = build_coarse_grained_samples(
        original_df=original_df,
        posets=posets,
        keep_ratios=tuple(cg_cfg["keep_ratios"]),
        cg_repeats=int(cg_cfg["repeats"]),
        sis_runs=int(exp_cfg["sis_runs"]),
        exact_threshold=int(exp_cfg.get("exact_threshold", 0)),
        zeta=float(cg_cfg.get("zeta", 1.0)),
        eta=float(cg_cfg.get("eta", 1.0)),
        gc_reference_windows=gc_reference_windows,
    )
    summary_df = summarize_cg(cg_df)
    rank_df = summarize_rank_shift(summary_df)
    write_outputs(original_df, cg_df, summary_df, rank_df, config)

    print("Coarse-grain summary:")
    print(summary_df.head(40).to_string(index=False))
    print("\nRank summary:")
    print(rank_df.head(40).to_string(index=False))


if __name__ == "__main__":
    main()
