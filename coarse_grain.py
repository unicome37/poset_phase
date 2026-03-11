from __future__ import annotations

import numpy as np

from generators import Poset


def coarse_grain_delete_nodes(
    poset: Poset,
    keep_ratio: float = 0.75,
    seed: int | None = None,
) -> Poset:
    """Random node deletion coarse graining via induced sub-poset.

    The closure restricted to an induced subset remains a valid transitive
    closure for the sub-poset, so no extra reconstruction is required.
    """
    n = poset.n
    if n <= 1:
        return Poset(poset.closure.copy())

    keep_ratio = float(np.clip(keep_ratio, 0.1, 1.0))
    keep_n = max(2, int(round(n * keep_ratio)))
    keep_n = min(keep_n, n)

    rng = np.random.default_rng(seed)
    keep_idx = np.sort(rng.choice(n, size=keep_n, replace=False))
    closure = poset.closure[np.ix_(keep_idx, keep_idx)].copy()
    np.fill_diagonal(closure, False)
    return Poset(closure)
