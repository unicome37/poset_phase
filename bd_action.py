from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from generators import Poset


@dataclass(frozen=True)
class IntervalCounts:
    """Causal interval size histogram for a poset.

    We define C_k as the number of related pairs (i ≺ j) with exactly k elements
    strictly between i and j. In particular:
      C_0 = number of links.
    """

    counts: dict[int, int]
    total_relations: int

    def get(self, k: int) -> int:
        return int(self.counts.get(k, 0))


def interval_size_matrix(poset: Poset) -> np.ndarray:
    """Return matrix M where M[i,j] = #{k : i≺k≺j}, for strict order i≺j.

    Computed as (C @ C) where C is the strict transitive closure (0/1).
    """
    c = poset.closure.astype(np.int32)
    return c @ c


def count_intervals_fast(poset: Poset, *, k_max: int | None = None) -> IntervalCounts:
    """Count causal intervals C_k using vectorized matrix operations.

    Args:
      k_max: if provided, only retain counts for k <= k_max (tail is dropped).
             This is useful when only low-k terms are needed (e.g. BDG d=4 uses k<=3).
    """
    if poset.n <= 1:
        return IntervalCounts(counts={}, total_relations=0)

    sizes = interval_size_matrix(poset)
    ks = sizes[poset.closure].astype(np.int64, copy=False)
    total = int(ks.size)
    if total == 0:
        return IntervalCounts(counts={}, total_relations=0)

    bc = np.bincount(ks)
    if k_max is None:
        nonzero = np.nonzero(bc)[0]
        out = {int(k): int(bc[k]) for k in nonzero}
        return IntervalCounts(counts=out, total_relations=total)

    k_max = int(k_max)
    upper = min(k_max, int(bc.size) - 1)
    out = {k: int(bc[k]) for k in range(upper + 1) if bc[k] != 0}
    return IntervalCounts(counts=out, total_relations=total)


def bdg_action_d2_link(counts: IntervalCounts, n: int, *, normalized: bool = True) -> float:
    """Minimal d=2 BD/BDG link action: S = N - 2*C_0."""
    s = float(n) - 2.0 * float(counts.get(0))
    return s / float(n) if normalized and n > 0 else s


def bdg_action_d2_corrected(counts: IntervalCounts, n: int, *, normalized: bool = True) -> float:
    """d=2 with first interval correction: S = N - 2*C_0 + 2*C_1."""
    s = float(n) - 2.0 * float(counts.get(0)) + 2.0 * float(counts.get(1))
    return s / float(n) if normalized and n > 0 else s


def bdg_action_d4_standard(counts: IntervalCounts, n: int, *, normalized: bool = True) -> float:
    """Standard literature d=4 BDG combination:

      S^(4) = N - C_0 + 9*C_1 - 16*C_2 + 8*C_3
    """
    s = (
        float(n)
        - 1.0 * float(counts.get(0))
        + 9.0 * float(counts.get(1))
        - 16.0 * float(counts.get(2))
        + 8.0 * float(counts.get(3))
    )
    return s / float(n) if normalized and n > 0 else s


def bd_action_d4_truncated(counts: IntervalCounts, n: int, *, normalized: bool = True) -> float:
    """BD scalar-curvature action proxy used in early bridge experiments.

    Matches the truncated coefficient form used in `prediction_a_bd_bridge.py`:

      S = N - (2/√6)C0 + (8/(3√6))C1 - (4/√6)C2
    """
    sqrt6 = math.sqrt(6.0)
    s = (
        float(n)
        - (2.0 / sqrt6) * float(counts.get(0))
        + (8.0 / (3.0 * sqrt6)) * float(counts.get(1))
        - (4.0 / sqrt6) * float(counts.get(2))
    )
    return s / float(n) if normalized and n > 0 else s


def bd_ratio_metric(counts: IntervalCounts) -> float:
    """A simple BD-inspired 'interval richness' metric.

    R = (fraction of non-links) * (1 + mean interval size),
    where interval size is k for pairs with k intervening elements.
    """
    if counts.total_relations <= 0:
        return 0.0
    total = float(counts.total_relations)
    c0 = float(counts.get(0))
    non_link_frac = 1.0 - c0 / total
    mean_k = 0.0
    for k, v in counts.counts.items():
        mean_k += float(k) * float(v)
    mean_k = mean_k / total
    return non_link_frac * (1.0 + mean_k)

