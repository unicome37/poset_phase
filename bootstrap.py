from __future__ import annotations

import numpy as np
import pandas as pd


def _bootstrap_mean(values: np.ndarray, n_bootstrap: int, rng: np.random.Generator) -> tuple[float, float]:
    draws = []
    n = len(values)
    for _ in range(n_bootstrap):
        sample = rng.choice(values, size=n, replace=True)
        draws.append(float(sample.mean()))
    arr = np.asarray(draws, dtype=float)
    return float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))


def _hierarchical_bootstrap_mean(
    means: np.ndarray,
    stds: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    draws = []
    n = len(means)
    for _ in range(n_bootstrap):
        idx = rng.choice(np.arange(n), size=n, replace=True)
        sampled_means = means[idx]
        sampled_stds = stds[idx]
        inner = rng.normal(loc=sampled_means, scale=sampled_stds)
        draws.append(float(inner.mean()))
    arr = np.asarray(draws, dtype=float)
    return float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))


def bootstrap_group_summary(
    df: pd.DataFrame,
    value_col: str = "score_norm",
    group_cols: tuple[str, ...] = ("n", "gamma", "family"),
    n_bootstrap: int = 1000,
    seed: int = 12345,
    value_std_col: str | None = None,
) -> pd.DataFrame:
    """Bootstrap using poset samples as the statistical unit.

    If `value_std_col` is provided, perform a two-level bootstrap:
    outer resampling over posets and inner Gaussian perturbation using the
    per-poset measurement error estimate.
    """
    rng = np.random.default_rng(seed)
    rows = []

    for key, group in df.groupby(list(group_cols)):
        values = group[value_col].to_numpy(dtype=float)
        mean = float(values.mean())
        std = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        if value_std_col is None:
            ci_low, ci_high = _bootstrap_mean(values, n_bootstrap=n_bootstrap, rng=rng)
            mode = "outer_only"
        else:
            stds = group[value_std_col].fillna(0.0).to_numpy(dtype=float)
            stds = np.clip(stds, a_min=0.0, a_max=None)
            ci_low, ci_high = _hierarchical_bootstrap_mean(
                values,
                stds,
                n_bootstrap=n_bootstrap,
                rng=rng,
            )
            mode = "outer_plus_measurement"

        if not isinstance(key, tuple):
            key = (key,)

        row = dict(zip(group_cols, key))
        row.update(
            {
                f"{value_col}_mean": mean,
                f"{value_col}_std": std,
                f"{value_col}_ci_low": ci_low,
                f"{value_col}_ci_high": ci_high,
                "bootstrap_mode": mode,
                "count": int(len(values)),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows).sort_values(list(group_cols))
