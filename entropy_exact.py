from __future__ import annotations

from functools import lru_cache
import math

from generators import Poset
from entropy_exact_c import ExactCBackendUnavailable, count_linear_extensions_exact_c


def count_linear_extensions_exact_python(poset: Poset) -> int:
    n = poset.n
    preds = []
    for j in range(n):
        mask = 0
        for i in range(n):
            if poset.closure[i, j]:
                mask |= 1 << i
        preds.append(mask)

    @lru_cache(maxsize=None)
    def dp(used: int) -> int:
        if used == (1 << n) - 1:
            return 1

        total = 0
        for v in range(n):
            if (used >> v) & 1:
                continue
            if preds[v] & ~used:
                continue
            total += dp(used | (1 << v))
        return total

    return dp(0)


def count_linear_extensions_exact(poset: Poset, prefer_c: bool = True) -> int:
    if prefer_c:
        try:
            return count_linear_extensions_exact_c(poset)
        except ExactCBackendUnavailable:
            pass
    return count_linear_extensions_exact_python(poset)


def log_linear_extensions_exact(poset: Poset) -> float:
    return math.log(count_linear_extensions_exact(poset))
