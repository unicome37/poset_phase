from __future__ import annotations

import numpy as np
import pandas as pd


def _safe_center_scale(values: pd.Series, method: str) -> tuple[pd.Series, float, float]:
    arr = values.astype(float)
    if method == "zscore":
        center = float(arr.mean())
        scale = arr.std(ddof=0)
        if scale == 0 or np.isnan(scale):
            return arr * 0.0, center, 0.0
        return (arr - center) / scale, center, float(scale)

    if method == "robust_zscore":
        median = arr.median()
        mad = (arr - median).abs().median()
        if mad == 0 or np.isnan(mad):
            return arr * 0.0, float(median), 0.0
        scale = 1.4826 * mad
        return (arr - median) / scale, float(median), float(scale)

    raise ValueError(f"Unsupported normalization method: {method}")


def add_normalized_columns(
    df: pd.DataFrame,
    method: str = "robust_zscore",
    group_cols: tuple[str, ...] = ("n",),
) -> pd.DataFrame:
    """Normalize H and penalty within each N-slice before recomputing score.

    The normalization is intentionally performed on the two components
    separately, rather than directly on the total score.
    """
    out = df.copy()
    log_h_norm_parts = []
    penalty_norm_parts = []
    h_scales = []
    penalty_scales = []

    for _, group in out.groupby(list(group_cols), sort=False):
        h_norm, _, h_scale = _safe_center_scale(group["log_H_mean"], method)
        p_norm, _, p_scale = _safe_center_scale(group["penalty_effective"], method)
        log_h_norm_parts.append(h_norm)
        penalty_norm_parts.append(p_norm)
        h_scales.append(pd.Series(h_scale, index=group.index))
        penalty_scales.append(pd.Series(p_scale, index=group.index))

    out["log_H_norm"] = pd.concat(log_h_norm_parts).sort_index()
    out["penalty_norm"] = pd.concat(penalty_norm_parts).sort_index()
    out["log_H_norm_scale"] = pd.concat(h_scales).sort_index()
    out["penalty_norm_scale"] = pd.concat(penalty_scales).sort_index()
    out["score_norm"] = -out["beta"] * out["log_H_norm"] + out["gamma"] * out["penalty_norm"]
    out["score_norm_std_est"] = np.where(
        out["log_H_norm_scale"] > 0,
        (out["beta"].abs() * out["log_H_std"]) / out["log_H_norm_scale"],
        0.0,
    )
    out["normalization_method"] = method
    out["normalization_group"] = ",".join(group_cols)
    return out


def add_size_scaled_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add simple explicit size normalizations as sensitivity-analysis columns."""
    out = df.copy()
    n = out["n"].astype(float)
    out["score_per_node"] = out["score"] / n
    out["score_per_nlogn"] = out["score"] / (n * np.log(np.maximum(n, 2.0)))
    return out
