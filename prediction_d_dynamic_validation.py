"""Prediction D dynamic-process validation: coarse-grained identity retention as a global filter.

Implements the three-layer evaluation described in:
  D:\\Kiro\\理论体系\\结构存在论_推论D验证方案.md

Layers:
  A) Identity retention: cg_family_switch_rate (lower is better)
  B) Rank stability: rank correlations / winner persistence under CG penalty
  C) Independent gain: does CG stability retain explanatory power after controlling
     local geometric terms? (partial correlation + permutation test)

This script is designed to run on existing outputs already present in this repo:
  - outputs_confirmatory/frozen_cg/cg_rank_summary.csv
  - outputs_exploratory/**/duel/pairwise_duel_raw.csv
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def zscore(arr: np.ndarray) -> np.ndarray:
    mean = float(arr.mean())
    std = float(arr.std(ddof=0))
    std = std if std > 1e-12 else 1.0
    return (arr - mean) / std


def regression_residual(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    xc = x - x.mean()
    yc = y - y.mean()
    denom = math.sqrt(float((xc * xc).sum() * (yc * yc).sum()))
    return float((xc * yc).sum() / denom) if denom > 1e-12 else 0.0


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return pearson(rx, ry)


def permutation_pvalue_abs_pearson(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    obs = abs(pearson(x, y))
    count = 1  # continuity correction
    for _ in range(int(n_perm)):
        y_perm = rng.permutation(y)
        if abs(pearson(x, y_perm)) >= obs:
            count += 1
    return float(count / (n_perm + 1))


def partial_corr_report(
    *,
    y: np.ndarray,
    x: np.ndarray,
    controls: np.ndarray,
    n_perm: int,
    seed: int,
) -> dict[str, float]:
    X = np.column_stack([np.ones(len(y)), controls])
    y_res = regression_residual(y, X)
    x_res = regression_residual(x, X)
    return {
        "partial_pearson": pearson(x_res, y_res),
        "partial_spearman": spearman(x_res, y_res),
        "permutation_p_abs_pearson": permutation_pvalue_abs_pearson(x_res, y_res, n_perm=n_perm, seed=seed),
        "n": float(len(y)),
    }


def discover_pairwise_duel_csvs(root: Path) -> list[Path]:
    return sorted(root.rglob("pairwise_duel_raw.csv"))


def load_pairwise_duel(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["dataset"] = csv_path.parent.as_posix()
    return df


def summarize_identity_retention_pairwise(df: pd.DataFrame) -> pd.DataFrame:
    out = (
        df.groupby(["dataset", "n", "family"])
        .agg(
            mean_switch_rate=("cg_family_switch_rate", "mean"),
            mean_cg_penalty=("cg_mean_penalty", "mean"),
            mean_geo_total=("geo_total", "mean"),
            mean_log_H=("log_H", "mean"),
            count=("seed", "count"),
        )
        .reset_index()
    )
    out["identity_retention"] = 1.0 - out["mean_switch_rate"]
    return out.sort_values(["n", "dataset", "identity_retention"], ascending=[True, True, False])


def zeta_crossings_pairwise(df: pd.DataFrame, zetas: list[float]) -> pd.DataFrame:
    rows: list[dict[str, float | str | int | None]] = []
    for (dataset, n), sub in df.groupby(["dataset", "n"]):
        families = sorted(sub["family"].unique().tolist())
        if len(families) != 2:
            # Some directories may contain incomplete or non-pairwise exports.
            continue
        if "lorentzian_like_2d" in families:
            fam_a = "lorentzian_like_2d"
            fam_b = families[0] if families[1] == fam_a else families[1]
        else:
            fam_a, fam_b = families[0], families[1]

        scan = []
        for zeta in zetas:
            score_eval = sub["score_A2_gamma"] + float(zeta) * sub["cg_mean_penalty"]
            means = score_eval.groupby(sub["family"]).mean()
            if fam_a not in means.index or fam_b not in means.index:
                continue
            delta = float(means[fam_b] - means[fam_a])
            scan.append({"zeta": float(zeta), "delta": delta})
        if not scan:
            continue
        scan_df = pd.DataFrame(scan).sort_values("zeta").reset_index(drop=True)
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
                "dataset": str(dataset),
                "n": int(n),
                "family_a": fam_a,
                "family_b": fam_b,
                "zeta_cross": cross,
                "delta_at_min_zeta": float(scan_df.iloc[0]["delta"]),
                "delta_at_max_zeta": float(scan_df.iloc[-1]["delta"]),
            }
        )
    return pd.DataFrame(rows).sort_values(["n", "dataset"]).reset_index(drop=True)


def load_cg_rank_summary(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["source"] = csv_path.parent.as_posix()
    return df


def discover_cg_raw_csvs(root: Path) -> list[Path]:
    # experiment_cg.py writes these in multiple output trees.
    candidates = []
    for p in root.rglob("cg_raw_samples.csv"):
        candidates.append(p)
    return sorted(candidates)


def load_cg_raw(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["source"] = csv_path.parent.as_posix()
    return df


def icg_from_penalty(series: pd.Series) -> pd.Series:
    # Higher I_cg = more stable = lower penalty_cg.
    return pd.Series(-zscore(series.to_numpy(float)), index=series.index)


def summarize_icg_redteam(cg_raw: pd.DataFrame) -> pd.DataFrame:
    cg_only = cg_raw[cg_raw["stage"] == "coarse_grained"].copy()
    group_cols = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode"]
    cg_only["I_cg"] = cg_only.groupby(group_cols)["penalty_cg"].transform(icg_from_penalty)

    summary = (
        cg_only.groupby(group_cols + ["family"])
        .agg(
            mean_I_cg=("I_cg", "mean"),
            mean_switch_rate=("family_switch_penalty", "mean"),
            mean_drift=("self_drift", "mean"),
            mean_rank_shift=("rank_shift_penalty", "mean"),
            mean_penalty_cg=("penalty_cg", "mean"),
            mean_score_local=("score_local", "mean"),
            mean_penalty_local=("penalty_local", "mean"),
            mean_log_H=("log_H", "mean"),
            count=("I_cg", "count"),
        )
        .reset_index()
    )
    summary["identity_retention"] = 1.0 - summary["mean_switch_rate"]
    return summary.sort_values(group_cols + ["mean_I_cg"], ascending=[True, True, True, True, True, True, False])


def independent_gain_cg_raw(
    cg_raw: pd.DataFrame,
    *,
    n_perm: int,
    seed: int,
) -> pd.DataFrame:
    cg_only = cg_raw[cg_raw["stage"] == "coarse_grained"].copy()
    group_cols = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode"]
    cg_only["I_cg"] = cg_only.groupby(group_cols)["penalty_cg"].transform(icg_from_penalty)

    # Controls: local action + a minimal geometric signature set.
    geo_cols = [c for c in ["penalty_local", "sig_comp", "sig_d_eff", "sig_height_ratio", "sig_width_ratio", "sig_degree_var"] if c in cg_only.columns]
    rows: list[dict[str, float | str | int]] = []
    for keys, sub in cg_only.groupby(group_cols, sort=False):
        sub = sub.dropna(subset=["score_local", "I_cg"] + geo_cols).copy()
        if len(sub) < 12:
            continue
        rep = partial_corr_report(
            y=sub["score_local"].to_numpy(float),
            x=sub["I_cg"].to_numpy(float),
            controls=sub[geo_cols].to_numpy(float),
            n_perm=n_perm,
            seed=seed + int(keys[1]) * 31 + int(round(float(keys[3]) * 1000)),
        )
        rep.update(
            {
                "source": str(keys[0]),
                "n": int(keys[1]),
                "n_cg": int(keys[2]),
                "keep_ratio": float(keys[3]),
                "gamma": float(keys[4]),
                "action_mode": str(keys[5]),
                "target_col": "score_local",
                "x_col": "I_cg",
                "n_controls": float(len(geo_cols)),
            }
        )
        rows.append(rep)
    return pd.DataFrame(rows).sort_values(["n", "keep_ratio", "gamma", "n_cg", "source"]).reset_index(drop=True)


def cg_rank_stability(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | str | int]] = []
    group_cols = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode"]
    for keys, sub in df.groupby(group_cols, sort=False):
        sub = sub.copy()
        r_local_vs_cg = spearman(sub["rank_local"].to_numpy(float), sub["rank_augmented_cg"].to_numpy(float))
        r_local_vs_full = spearman(sub["rank_local"].to_numpy(float), sub["rank_augmented_full"].to_numpy(float))
        winner_local = str(sub.sort_values("rank_local").iloc[0]["family"])
        winner_cg = str(sub.sort_values("rank_augmented_cg").iloc[0]["family"])
        winner_full = str(sub.sort_values("rank_augmented_full").iloc[0]["family"])
        rows.append(
            {
                "source": str(keys[0]),
                "n": int(keys[1]),
                "n_cg": int(keys[2]),
                "keep_ratio": float(keys[3]),
                "gamma": float(keys[4]),
                "action_mode": str(keys[5]),
                "spearman_rank_local_vs_aug_cg": float(r_local_vs_cg),
                "spearman_rank_local_vs_aug_full": float(r_local_vs_full),
                "winner_local": winner_local,
                "winner_aug_cg": winner_cg,
                "winner_aug_full": winner_full,
                "winner_switch_due_to_cg": float(winner_local != winner_cg),
                "winner_switch_due_to_full": float(winner_local != winner_full),
                "n_families": float(len(sub)),
            }
        )
    return pd.DataFrame(rows).sort_values(["n", "gamma", "keep_ratio", "n_cg"])


def _family_rank_frame(
    cg_only: pd.DataFrame,
    *,
    zeta: float,
    group_cols: list[str],
) -> pd.DataFrame:
    sub = cg_only.copy()
    sub["score_eval"] = sub["score_local"].astype(float) + float(zeta) * sub["penalty_cg"].astype(float)
    fam_means = (
        sub.groupby(group_cols + ["family"])
        .agg(
            mean_score_local=("score_local", "mean"),
            mean_penalty_cg=("penalty_cg", "mean"),
            mean_score_eval=("score_eval", "mean"),
            mean_I_cg=("I_cg", "mean"),
            mean_penalty_local=("penalty_local", "mean"),
            mean_log_H=("log_H", "mean"),
            mean_sig_comp=("sig_comp", "mean"),
            mean_sig_d_eff=("sig_d_eff", "mean"),
            mean_sig_height_ratio=("sig_height_ratio", "mean"),
            mean_sig_width_ratio=("sig_width_ratio", "mean"),
            mean_sig_degree_var=("sig_degree_var", "mean"),
            mean_retain_identity=("retain_identity", "mean"),
            count=("score_eval", "count"),
        )
        .reset_index()
    )
    fam_means["zeta"] = float(zeta)
    fam_means["rank_eval"] = (
        fam_means.groupby(group_cols)["mean_score_eval"]
        .rank(method="dense", ascending=True)
        .astype(int)
    )
    fam_means["rank_local_eval"] = (
        fam_means.groupby(group_cols)["mean_score_local"]
        .rank(method="dense", ascending=True)
        .astype(int)
    )
    fam_means["is_winner_eval"] = (fam_means["rank_eval"] == 1).astype(float)
    fam_means["is_winner_local"] = (fam_means["rank_local_eval"] == 1).astype(float)
    return fam_means


def zeta_scan_rankings(cg_raw: pd.DataFrame, zetas: list[float]) -> pd.DataFrame:
    cg_only = cg_raw[cg_raw["stage"] == "coarse_grained"].copy()
    group_cols = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode"]
    cg_only["I_cg"] = cg_only.groupby(group_cols)["penalty_cg"].transform(icg_from_penalty)
    cg_only["retain_identity"] = 1.0 - cg_only["family_switch_penalty"].astype(float)

    frames = []
    for zeta in zetas:
        frames.append(_family_rank_frame(cg_only, zeta=zeta, group_cols=group_cols))
    return pd.concat(frames, ignore_index=True)


def compute_penalty_variants(cg_only: pd.DataFrame) -> pd.DataFrame:
    """Return cg_only with additional penalty/I_cg variants.

    Purpose: red-team the possibility that D is an artifact of the centroid-based
    family-switch classifier. Variants excluding family_switch provide a check.
    """
    out = cg_only.copy()
    # Default weights in stability.coarse_grain_penalty.
    out["penalty_full"] = out["penalty_cg"].astype(float)
    out["penalty_drift"] = 10.0 * out["self_drift"].astype(float)
    out["penalty_switch"] = 1.5 * out["family_switch_penalty"].astype(float)
    out["penalty_rank"] = 0.5 * out["rank_shift_penalty"].astype(float)
    out["penalty_no_switch"] = out["penalty_drift"] + out["penalty_rank"]
    return out


def zeta_scan_rankings_variants(cg_raw: pd.DataFrame, zetas: list[float]) -> pd.DataFrame:
    cg_only = cg_raw[cg_raw["stage"] == "coarse_grained"].copy()
    group_cols = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode"]
    cg_only["retain_identity"] = 1.0 - cg_only["family_switch_penalty"].astype(float)
    cg_only = compute_penalty_variants(cg_only)

    variants = [
        ("full", "penalty_full"),
        ("no_switch", "penalty_no_switch"),
        ("drift", "penalty_drift"),
        ("rank", "penalty_rank"),
        ("switch", "penalty_switch"),
    ]

    frames = []
    for variant_name, penalty_col in variants:
        sub = cg_only.copy()
        sub["penalty_variant"] = sub[penalty_col].astype(float)
        sub["I_cg_variant"] = sub.groupby(group_cols)["penalty_variant"].transform(icg_from_penalty)
        for zeta in zetas:
            sub2 = sub.copy()
            sub2["I_cg"] = sub2["I_cg_variant"]
            sub2["penalty_cg"] = sub2["penalty_variant"]
            fam = _family_rank_frame(sub2, zeta=zeta, group_cols=group_cols)
            fam["icg_variant"] = str(variant_name)
            frames.append(fam)
    return pd.concat(frames, ignore_index=True)


def block_permutation_delta_rank_test(
    zeta_rank_df: pd.DataFrame,
    *,
    n_perm: int,
    seed: int,
) -> pd.DataFrame:
    """Block permutation test for 'I_cg predicts rank improvement due to CG penalty'.

    Target: delta_rank = rank_eval - rank_local_eval (negative = improved rank).
    Statistic: mean Spearman across blocks between mean_I_cg and (-delta_rank).
    Permutation: within each block, permute mean_I_cg across families.
    """
    rng = np.random.default_rng(seed)
    group_cols = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode", "zeta", "icg_variant"]

    base = zeta_rank_df.copy()
    base["delta_rank"] = base["rank_eval"].astype(float) - base["rank_local_eval"].astype(float)
    base["improve_rank"] = -base["delta_rank"]

    blocks = []
    for keys, sub in base.groupby(group_cols, sort=False):
        if sub["family"].nunique() < 4:
            continue
        x = sub["mean_I_cg"].to_numpy(float)
        y = sub["improve_rank"].to_numpy(float)
        blocks.append((keys, x, y))

    if not blocks:
        return pd.DataFrame()

    def stat(block_list: list[tuple]) -> float:
        rs = []
        for _keys, x_, y_ in block_list:
            rs.append(spearman(x_, y_))
        return float(np.mean(rs)) if rs else 0.0

    obs = stat(blocks)
    count = 1
    for _ in range(int(n_perm)):
        perm_blocks = []
        for keys, x, y in blocks:
            perm_blocks.append((keys, rng.permutation(x), y))
        if abs(stat(perm_blocks)) >= abs(obs):
            count += 1
    p = float(count / (int(n_perm) + 1))

    # Also report per-variant obs (no permutation, descriptive).
    variant_rows = []
    for variant, sub in base.groupby("icg_variant", sort=False):
        sub_blocks = []
        for keys, s2 in sub.groupby(["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode", "zeta"], sort=False):
            if s2["family"].nunique() < 4:
                continue
            sub_blocks.append(
                (
                    keys,
                    s2["mean_I_cg"].to_numpy(float),
                    s2["improve_rank"].to_numpy(float),
                )
            )
        variant_rows.append(
            {
                "icg_variant": str(variant),
                "n_blocks": float(len(sub_blocks)),
                "mean_spearman_obs": float(stat(sub_blocks)) if sub_blocks else float("nan"),
            }
        )
    out = pd.DataFrame(variant_rows).sort_values("icg_variant").reset_index(drop=True)
    out.insert(0, "blockperm_p_abs_mean_spearman", p)
    out.insert(0, "blockperm_obs_abs_mean_spearman", float(abs(obs)))
    out.insert(0, "n_perm", float(n_perm))
    out.insert(0, "n_blocks_total", float(len(blocks)))
    return out


def stratified_blockperm_delta_rank_test(
    zeta_rank_df: pd.DataFrame,
    *,
    strata_cols: tuple[str, ...] = ("icg_variant", "n", "keep_ratio", "gamma"),
    n_perm: int,
    seed: int,
) -> pd.DataFrame:
    """Stratified block permutation test on delta-rank improvement.

    Target: delta_rank = rank_eval - rank_local_eval (negative = improved rank).
    Statistic per stratum: mean Spearman(mean_I_cg, -delta_rank) across blocks.
    Permutation: within each block, permute mean_I_cg across families.

    Implementation: Spearman is computed via Pearson correlation of within-block ranks.
    Pre-centered rank vectors make each permutation a dot product.
    """
    rng = np.random.default_rng(seed)

    base = zeta_rank_df.copy()
    base["delta_rank"] = base["rank_eval"].astype(float) - base["rank_local_eval"].astype(float)
    base["improve_rank"] = -base["delta_rank"]

    block_cols = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode", "zeta", "icg_variant"]
    needed = set(block_cols + list(strata_cols) + ["family", "mean_I_cg", "improve_rank"])
    missing = [c for c in needed if c not in base.columns]
    if missing:
        raise ValueError(f"Missing required columns for stratified test: {missing}")

    blocks = []
    for _keys, sub in base.groupby(block_cols, sort=False):
        if sub["family"].nunique() < 4:
            continue
        x = sub["mean_I_cg"].to_numpy(float)
        y = sub["improve_rank"].to_numpy(float)

        rx = np.argsort(np.argsort(x)).astype(float)
        ry = np.argsort(np.argsort(y)).astype(float)
        rx_c = rx - float(rx.mean())
        ry_c = ry - float(ry.mean())
        denom = math.sqrt(float((rx_c * rx_c).sum() * (ry_c * ry_c).sum()))
        denom = denom if denom > 1e-12 else 1.0

        strata_key = tuple(sub.iloc[0][c] for c in strata_cols)
        blocks.append((strata_key, rx_c, ry_c, denom))

    if not blocks:
        return pd.DataFrame()

    strata = sorted({b[0] for b in blocks})
    idx = {k: i for i, k in enumerate(strata)}
    n_blocks = np.zeros(len(strata), dtype=int)
    obs_sum = np.zeros(len(strata), dtype=float)

    for sk, rx_c, ry_c, denom in blocks:
        i = idx[sk]
        n_blocks[i] += 1
        obs_sum[i] += float(np.dot(rx_c, ry_c) / denom)

    obs_mean = obs_sum / np.maximum(n_blocks, 1)
    count = np.ones(len(strata), dtype=int)

    for _ in range(int(n_perm)):
        perm_sum = np.zeros(len(strata), dtype=float)
        for sk, rx_c, ry_c, denom in blocks:
            i = idx[sk]
            perm_idx = rng.permutation(len(rx_c))
            perm_sum[i] += float(np.dot(rx_c[perm_idx], ry_c) / denom)
        perm_mean = perm_sum / np.maximum(n_blocks, 1)
        count += (np.abs(perm_mean) >= np.abs(obs_mean)).astype(int)

    p = count / float(int(n_perm) + 1)

    rows = []
    for sk in strata:
        i = idx[sk]
        row: dict[str, float | str] = {
            "n_perm": float(n_perm),
            "n_blocks": float(n_blocks[i]),
            "obs_mean_spearman": float(obs_mean[i]),
            "p_abs_mean_spearman": float(p[i]),
        }
        for j, col in enumerate(strata_cols):
            val = sk[j]
            if col in {"n", "n_cg"}:
                row[col] = float(val)
            elif col in {"keep_ratio", "gamma"}:
                row[col] = float(val)
            else:
                row[col] = str(val)
        rows.append(row)

    out = pd.DataFrame(rows)
    sort_cols = [c for c in ["icg_variant", "n", "keep_ratio", "gamma"] if c in out.columns]
    if sort_cols:
        out = out.sort_values(sort_cols).reset_index(drop=True)
    return out


def zeta_scan_stability(zeta_rank_df: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode"]
    rows: list[dict[str, float | str | int]] = []
    for keys, sub in zeta_rank_df.groupby(group_cols, sort=False):
        winners = (
            sub[sub["is_winner_eval"] == 1.0][["zeta", "family"]]
            .sort_values("zeta")
            .drop_duplicates(subset=["zeta"])
        )
        winner_list = winners["family"].tolist()
        winner_switches = sum(1 for i in range(1, len(winner_list)) if winner_list[i] != winner_list[i - 1])

        zetas = sorted(sub["zeta"].unique().tolist())
        z0 = zetas[0]
        r0 = sub[sub["zeta"] == z0].set_index("family")["rank_eval"].to_dict()
        spears = []
        for z in zetas[1:]:
            rz = sub[sub["zeta"] == z].set_index("family")["rank_eval"].to_dict()
            fams = sorted(set(r0.keys()) & set(rz.keys()))
            if len(fams) < 3:
                continue
            spears.append(
                spearman(
                    np.asarray([r0[f] for f in fams], dtype=float),
                    np.asarray([rz[f] for f in fams], dtype=float),
                )
            )
        spears = [float(v) for v in spears]
        winner_persistence = float(winners["family"].value_counts().max() / max(len(zetas), 1))

        lor2d_rows = sub[(sub["family"] == "lorentzian_like_2d") & (sub["is_winner_eval"] == 1.0)]
        zeta_lor2d_win = float(lor2d_rows["zeta"].min()) if len(lor2d_rows) else float("nan")
        rows.append(
            {
                "source": str(keys[0]),
                "n": int(keys[1]),
                "n_cg": int(keys[2]),
                "keep_ratio": float(keys[3]),
                "gamma": float(keys[4]),
                "action_mode": str(keys[5]),
                "n_families": float(sub["family"].nunique()),
                "n_zetas": float(len(zetas)),
                "winner_switches": float(winner_switches),
                "winner_persistence_frac": float(winner_persistence),
                "spearman_rank_continuity_mean": float(np.mean(spears)) if spears else float("nan"),
                "spearman_rank_continuity_min": float(np.min(spears)) if spears else float("nan"),
                "zeta_lor2d_first_win": float(zeta_lor2d_win),
            }
        )
    return pd.DataFrame(rows).sort_values(["n", "keep_ratio", "gamma", "n_cg", "source"]).reset_index(drop=True)


def independent_gain_winner_rank(zeta_rank_df: pd.DataFrame, *, n_perm: int, seed: int) -> pd.DataFrame:
    controls = [
        "mean_penalty_local",
        "mean_sig_comp",
        "mean_sig_d_eff",
        "mean_sig_height_ratio",
        "mean_sig_width_ratio",
        "mean_sig_degree_var",
        "mean_score_local",
        "mean_log_H",
        "keep_ratio",
        "gamma",
        "n",
        "n_cg",
    ]
    base = zeta_rank_df.copy()
    x_col = "mean_I_cg"
    base = base.dropna(subset=[x_col, "rank_eval", "is_winner_eval"] + [c for c in controls if c in base.columns])
    X_cols = [c for c in controls if c in base.columns]
    if len(base) < 30 or len(X_cols) < 3:
        return pd.DataFrame()

    out = []
    for target_col in ["is_winner_eval", "rank_eval"]:
        rep = partial_corr_report(
            y=base[target_col].to_numpy(float),
            x=base[x_col].to_numpy(float),
            controls=base[X_cols].to_numpy(float),
            n_perm=n_perm,
            seed=seed + (0 if target_col == "is_winner_eval" else 1),
        )
        rep.update(
            {
                "target_col": target_col,
                "x_col": x_col,
                "n_controls": float(len(X_cols)),
                "n_rows": float(len(base)),
            }
        )
        out.append(rep)
    return pd.DataFrame(out)


def independent_gain_pairwise(
    df: pd.DataFrame,
    *,
    target_col: str,
    cg_col: str,
    geo_cols: list[str],
    n_perm: int,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict[str, float | str | int]] = []
    for (dataset, n), sub in df.groupby(["dataset", "n"]):
        sub = sub.dropna(subset=[target_col, cg_col] + geo_cols).copy()
        if len(sub) < 8:
            continue
        y = sub[target_col].to_numpy(float)
        x = sub[cg_col].to_numpy(float)
        controls = sub[geo_cols].to_numpy(float)
        rep = partial_corr_report(y=y, x=x, controls=controls, n_perm=n_perm, seed=seed + int(n))
        rep.update(
            {
                "dataset": str(dataset),
                "n": int(n),
                "target_col": target_col,
                "cg_col": cg_col,
                "n_controls": float(len(geo_cols)),
            }
        )
        rows.append(rep)
    return pd.DataFrame(rows).sort_values(["n", "dataset"]).reset_index(drop=True)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prediction D dynamic validation: coarse-graining identity retention.")
    parser.add_argument("--out", default="outputs_exploratory/prediction_d_dynamic", help="Output directory.")
    parser.add_argument("--n-perm", type=int, default=5000, help="Permutation count for p-values.")
    parser.add_argument(
        "--n-perm-blockperm",
        type=int,
        default=0,
        help="Permutation count for pooled block-permutation tests (0 = use --n-perm).",
    )
    parser.add_argument(
        "--n-perm-stratified",
        type=int,
        default=0,
        help="Permutation count for stratified block-permutation tests (0 = use --n-perm).",
    )
    parser.add_argument("--seed", type=int, default=7, help="RNG seed base.")
    parser.add_argument(
        "--only-blockperm",
        action="store_true",
        help="If set, only run the block-permutation red-team tests (fast path).",
    )
    parser.add_argument(
        "--cg-source-contains",
        default="",
        help="If set, only include cg_*.csv whose path contains this substring.",
    )
    parser.add_argument(
        "--zeta-scan",
        default="0,0.5,1,1.5,2,2.5,3,3.5,4,5,6",
        help="Comma-separated zeta values for winner/rank scan.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    zeta_scan = [float(x.strip()) for x in str(args.zeta_scan).split(",") if x.strip()]

    n_perm_default = int(args.n_perm)
    n_perm_blockperm = int(args.n_perm_blockperm) if int(args.n_perm_blockperm) > 0 else n_perm_default
    n_perm_stratified = int(args.n_perm_stratified) if int(args.n_perm_stratified) > 0 else n_perm_default

    # Pairwise-duel (Lor2D vs MLR) stability summaries.
    duel_csvs = discover_pairwise_duel_csvs(Path("outputs_exploratory"))
    duel_csvs += [Path("outputs_exploratory/pairwise_compressibility_duel_expanded/pairwise_duel_raw.csv")]
    duel_csvs = [p for p in duel_csvs if p.exists()]

    duel_frames = [load_pairwise_duel(p) for p in duel_csvs]
    duel_df = pd.concat(duel_frames, ignore_index=True) if duel_frames else pd.DataFrame()

    if not duel_df.empty:
        identity = summarize_identity_retention_pairwise(duel_df)
        identity.to_csv(out_dir / "pairwise_identity_retention.csv", index=False, encoding="utf-8-sig")

        zetas = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        crossings = zeta_crossings_pairwise(duel_df, zetas=zetas)
        crossings.to_csv(out_dir / "pairwise_zeta_crossings.csv", index=False, encoding="utf-8-sig")

        duel_df = duel_df.copy()
        duel_df["switch_zscore"] = duel_df.groupby("n")["cg_family_switch_rate"].transform(
            lambda s: pd.Series(zscore(s.to_numpy(float)), index=s.index)
        )
        # Independent gain: does CG stability still relate to entropy/score after controlling geo terms?
        geo_cols = [
            "antichain_width",
            "comparable_fraction",
            "geo_dim_eff",
            "geo_cover_density",
            "geo_layer_smoothness",
            "geo_total",
        ]
        geo_cols = [c for c in geo_cols if c in duel_df.columns]
        ind_score = independent_gain_pairwise(
            duel_df,
            target_col="score_A2_gamma",
            cg_col="switch_zscore",
            geo_cols=geo_cols,
            n_perm=int(args.n_perm),
            seed=int(args.seed),
        )
        ind_score.to_csv(out_dir / "pairwise_independent_gain_score.csv", index=False, encoding="utf-8-sig")

        ind_entropy = independent_gain_pairwise(
            duel_df,
            target_col="log_H",
            cg_col="switch_zscore",
            geo_cols=geo_cols,
            n_perm=int(args.n_perm),
            seed=int(args.seed) + 1000,
        )
        ind_entropy.to_csv(out_dir / "pairwise_independent_gain_entropy.csv", index=False, encoding="utf-8-sig")

    # CG experiment outputs (multi-family): rank stability + red-team I_cg.
    cg_rank_paths = []
    for candidate in [
        Path("outputs_confirmatory/frozen_cg/cg_rank_summary.csv"),
        Path("outputs_frozen_cg/cg_rank_summary.csv"),
        Path("outputs_cg_dimension/cg_rank_summary.csv"),
        Path("outputs_exploratory/cg_dimension/cg_rank_summary.csv"),
        Path("outputs_exploratory/cg_pilot/cg_rank_summary.csv"),
        Path("outputs_exploratory/cg_smoke/cg_rank_summary.csv"),
    ]:
        if candidate.exists():
            cg_rank_paths.append(candidate)

    cg_frames = [load_cg_rank_summary(p) for p in cg_rank_paths]
    cg_df = pd.concat(cg_frames, ignore_index=True) if cg_frames else pd.DataFrame()

    if not cg_df.empty:
        stability = cg_rank_stability(cg_df)
        stability.to_csv(out_dir / "cg_rank_stability.csv", index=False, encoding="utf-8-sig")

        cg_df.to_csv(out_dir / "cg_rank_summary_merged.csv", index=False, encoding="utf-8-sig")

    cg_raw_paths = discover_cg_raw_csvs(Path("."))
    if args.cg_source_contains:
        cg_raw_paths = [p for p in cg_raw_paths if args.cg_source_contains in p.as_posix()]
    cg_raw_frames = [load_cg_raw(p) for p in cg_raw_paths]
    cg_raw_df = pd.concat(cg_raw_frames, ignore_index=True) if cg_raw_frames else pd.DataFrame()

    if not cg_raw_df.empty:
        if bool(args.only_blockperm):
            # Fast path: just produce the dynamic-process red-team (variants + block permutation).
            zeta_rank_var = zeta_scan_rankings_variants(cg_raw_df, zetas=zeta_scan)
            zeta_rank_var.to_csv(out_dir / "cg_zeta_scan_rankings_variants.csv", index=False, encoding="utf-8-sig")
            zeta_stab_var = zeta_scan_stability(zeta_rank_var)
            zeta_stab_var.to_csv(out_dir / "cg_zeta_scan_stability_variants.csv", index=False, encoding="utf-8-sig")

            blockperm = block_permutation_delta_rank_test(
                zeta_rank_var,
                n_perm=n_perm_blockperm,
                seed=int(args.seed) + 9000,
            )
            if not blockperm.empty:
                blockperm.to_csv(out_dir / "cg_blockperm_delta_rank.csv", index=False, encoding="utf-8-sig")

            strat = stratified_blockperm_delta_rank_test(
                zeta_rank_var,
                strata_cols=("icg_variant", "n", "keep_ratio", "gamma"),
                n_perm=n_perm_stratified,
                seed=int(args.seed) + 9100,
            )
            if not strat.empty:
                strat.to_csv(out_dir / "cg_blockperm_delta_rank_stratified.csv", index=False, encoding="utf-8-sig")

            lines = ["# Prediction D Dynamic Validation (Auto Report)", ""]
            lines.append(f"- Generated at: `{out_dir.as_posix()}`")
            lines.append(f"- CG raw rows: `{len(cg_raw_df)}` from `{len(cg_raw_paths)}` CSVs")
            if args.cg_source_contains:
                lines.append(f"- Filter: `cg_source_contains={args.cg_source_contains}`")
            lines.append(f"- n_perm_blockperm: `{n_perm_blockperm}`")
            lines.append(f"- n_perm_stratified: `{n_perm_stratified}`")
            lines.append("")
            lines.append("Outputs:")
            lines.append(f"- `cg_zeta_scan_rankings_variants.csv`")
            lines.append(f"- `cg_zeta_scan_stability_variants.csv`")
            lines.append(f"- `cg_blockperm_delta_rank.csv`")
            lines.append(f"- `cg_blockperm_delta_rank_stratified.csv`")
            (out_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
            return

        redteam = summarize_icg_redteam(cg_raw_df)
        redteam.to_csv(out_dir / "cg_redteam_icg_summary.csv", index=False, encoding="utf-8-sig")
        ind_gain = independent_gain_cg_raw(cg_raw_df, n_perm=n_perm_default, seed=int(args.seed) + 2000)
        ind_gain.to_csv(out_dir / "cg_independent_gain_score_local.csv", index=False, encoding="utf-8-sig")

        # Retention as target (linear probability proxy): does I_cg keep explanatory power after controls?
        cg_only = cg_raw_df[cg_raw_df["stage"] == "coarse_grained"].copy()
        group_cols = ["source", "n", "n_cg", "keep_ratio", "gamma", "action_mode"]
        cg_only["I_cg"] = cg_only.groupby(group_cols)["penalty_cg"].transform(icg_from_penalty)
        cg_only["retain_identity"] = 1.0 - cg_only["family_switch_penalty"].astype(float)
        geo_cols = [c for c in ["penalty_local", "sig_comp", "sig_d_eff", "sig_height_ratio", "sig_width_ratio", "sig_degree_var"] if c in cg_only.columns]
        rows = []
        for keys, sub in cg_only.groupby(group_cols, sort=False):
            sub = sub.dropna(subset=["retain_identity", "I_cg"] + geo_cols).copy()
            if len(sub) < 12:
                continue
            rep = partial_corr_report(
                y=sub["retain_identity"].to_numpy(float),
                x=sub["I_cg"].to_numpy(float),
                controls=sub[geo_cols].to_numpy(float),
                n_perm=n_perm_default,
                seed=int(args.seed) + 5000 + int(keys[1]) * 37,
            )
            rep.update(
                {
                    "source": str(keys[0]),
                    "n": int(keys[1]),
                    "n_cg": int(keys[2]),
                    "keep_ratio": float(keys[3]),
                    "gamma": float(keys[4]),
                    "action_mode": str(keys[5]),
                    "target_col": "retain_identity",
                    "x_col": "I_cg",
                    "n_controls": float(len(geo_cols)),
                }
            )
            rows.append(rep)
        pd.DataFrame(rows).sort_values(["n", "keep_ratio", "gamma", "n_cg", "source"]).to_csv(
            out_dir / "cg_independent_gain_retention.csv",
            index=False,
            encoding="utf-8-sig",
        )

        zeta_rank_df = zeta_scan_rankings(cg_raw_df, zetas=zeta_scan)
        zeta_rank_df.to_csv(out_dir / "cg_zeta_scan_rankings.csv", index=False, encoding="utf-8-sig")
        zeta_stab = zeta_scan_stability(zeta_rank_df)
        zeta_stab.to_csv(out_dir / "cg_zeta_scan_stability.csv", index=False, encoding="utf-8-sig")
        ind_wr = independent_gain_winner_rank(zeta_rank_df, n_perm=n_perm_default, seed=int(args.seed) + 7000)
        if not ind_wr.empty:
            ind_wr.to_csv(out_dir / "cg_independent_gain_winner_rank.csv", index=False, encoding="utf-8-sig")

        # Red-team 2.0: variants excluding family-switch, and block permutation on delta-rank.
        zeta_rank_var = zeta_scan_rankings_variants(cg_raw_df, zetas=zeta_scan)
        zeta_rank_var.to_csv(out_dir / "cg_zeta_scan_rankings_variants.csv", index=False, encoding="utf-8-sig")
        zeta_stab_var = zeta_scan_stability(zeta_rank_var)
        zeta_stab_var.to_csv(out_dir / "cg_zeta_scan_stability_variants.csv", index=False, encoding="utf-8-sig")
        blockperm = block_permutation_delta_rank_test(
            zeta_rank_var,
            n_perm=n_perm_blockperm,
            seed=int(args.seed) + 9000,
        )
        if not blockperm.empty:
            blockperm.to_csv(out_dir / "cg_blockperm_delta_rank.csv", index=False, encoding="utf-8-sig")

        strat = stratified_blockperm_delta_rank_test(
            zeta_rank_var,
            strata_cols=("icg_variant", "n", "keep_ratio", "gamma"),
            n_perm=n_perm_stratified,
            seed=int(args.seed) + 9100,
        )
        if not strat.empty:
            strat.to_csv(out_dir / "cg_blockperm_delta_rank_stratified.csv", index=False, encoding="utf-8-sig")

    # Short Markdown report for quick reading.
    lines = ["# Prediction D Dynamic Validation (Auto Report)", ""]
    lines.append(f"- Generated at: `{out_dir.as_posix()}`")
    if not duel_df.empty:
        lines.append(f"- Pairwise duel rows: `{len(duel_df)}` from `{len(duel_csvs)}` CSVs")
    if not cg_df.empty:
        lines.append(f"- CG rank rows: `{len(cg_df)}` from `{len(cg_rank_paths)}` CSVs")
    if not cg_raw_df.empty:
        lines.append(f"- CG raw rows: `{len(cg_raw_df)}` from `{len(cg_raw_paths)}` CSVs")
    lines.append("")
    lines.append("Outputs:")
    if not duel_df.empty:
        lines.append(f"- `pairwise_identity_retention.csv`")
        lines.append(f"- `pairwise_zeta_crossings.csv`")
        lines.append(f"- `pairwise_independent_gain_score.csv`")
        lines.append(f"- `pairwise_independent_gain_entropy.csv`")
    if not cg_df.empty:
        lines.append(f"- `cg_rank_stability.csv`")
        lines.append(f"- `cg_rank_summary_merged.csv`")
    if not cg_raw_df.empty:
        lines.append(f"- `cg_redteam_icg_summary.csv`")
        lines.append(f"- `cg_independent_gain_score_local.csv`")
        lines.append(f"- `cg_independent_gain_retention.csv`")
        lines.append(f"- `cg_zeta_scan_rankings.csv`")
        lines.append(f"- `cg_zeta_scan_stability.csv`")
        lines.append(f"- `cg_independent_gain_winner_rank.csv`")
        lines.append(f"- `cg_zeta_scan_rankings_variants.csv`")
        lines.append(f"- `cg_zeta_scan_stability_variants.csv`")
        lines.append(f"- `cg_blockperm_delta_rank.csv`")
        lines.append(f"- `cg_blockperm_delta_rank_stratified.csv`")
    (out_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
