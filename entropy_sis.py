from __future__ import annotations

import math

import numpy as np

from generators import Poset


def estimate_log_linear_extensions_sis(
    poset: Poset,
    n_runs: int = 256,
    seed: int | None = None,
) -> tuple[float, float]:
    """SIS 估计 log(number of linear extensions).

    返回:
    - log_mean_estimate
    - std_of_log_weights
    """
    rng = np.random.default_rng(seed)
    n = poset.n
    preds_template = [set(np.where(poset.closure[:, j])[0]) for j in range(n)]
    estimates = []

    for _ in range(n_runs):
        preds = [s.copy() for s in preds_template]
        remaining = set(range(n))
        available = {v for v in range(n) if not preds[v]}
        log_w = 0.0

        while remaining:
            k = len(available)
            if k == 0:
                log_w = float("-inf")
                break

            log_w += math.log(k)
            chosen = rng.choice(list(available))
            remaining.remove(chosen)
            available.remove(chosen)

            for v in list(remaining):
                if chosen in preds[v]:
                    preds[v].remove(chosen)
                    if not preds[v]:
                        available.add(v)

        estimates.append(log_w)

    arr = np.asarray(estimates, dtype=float)
    m = float(np.max(arr))
    log_mean = m + math.log(np.mean(np.exp(arr - m)))
    log_std = float(arr.std())
    return log_mean, log_std
