from __future__ import annotations

import numpy as np

from generators import Poset


def comparable_fraction(poset: Poset) -> float:
    c = poset.closure
    n = poset.n
    total = n * (n - 1) / 2
    comp = np.triu(c | c.T, 1).sum()
    return float(comp / total) if total else 0.0


def degree_stats(poset: Poset) -> dict[str, float]:
    c = poset.closure.astype(np.int32)
    out_deg = c.sum(axis=1)
    in_deg = c.sum(axis=0)
    return {
        "in_mean": float(in_deg.mean()),
        "out_mean": float(out_deg.mean()),
        "in_var": float(in_deg.var()),
        "out_var": float(out_deg.var()),
    }


def layer_profile(poset: Poset) -> np.ndarray:
    """按最小元分层，返回每层元素数量。"""
    c = poset.closure
    indeg = c.sum(axis=0).astype(int)
    remaining = set(range(poset.n))
    layers = []

    while remaining:
        mins = [i for i in remaining if indeg[i] == 0]
        if not mins:
            raise ValueError("Input relation is not a poset: cycle detected.")
        layers.append(len(mins))
        for u in mins:
            remaining.remove(u)
            for v in np.where(c[u])[0]:
                indeg[v] -= 1

    return np.asarray(layers, dtype=float)


def normalized_degree_variance(poset: Poset) -> float:
    stats = degree_stats(poset)
    n = max(poset.n, 1)
    return (stats["in_var"] + stats["out_var"]) / (2.0 * n * n)


def layer_imbalance(poset: Poset) -> float:
    profile = layer_profile(poset)
    profile = profile / profile.sum()
    return float(((profile - profile.mean()) ** 2).mean())


def antichain_width(poset: Poset) -> int:
    """Maximum antichain size via Dilworth reduction to bipartite matching.

    For a finite poset, width = n - maximum matching size in the bipartite graph
    with left/right copies of vertices and edges i_L -> j_R whenever i < j.
    """
    n = poset.n
    closure = poset.closure
    adj = [list(np.where(closure[i])[0]) for i in range(n)]
    match_r = [-1] * n

    def try_augment(u: int, seen: list[bool]) -> bool:
        for v in adj[u]:
            if seen[v]:
                continue
            seen[v] = True
            if match_r[v] == -1 or try_augment(match_r[v], seen):
                match_r[v] = u
                return True
        return False

    match_size = 0
    for u in range(n):
        seen = [False] * n
        if try_augment(u, seen):
            match_size += 1
    return int(n - match_size)


def antichain_width_ratio(poset: Poset) -> float:
    return float(antichain_width(poset) / max(poset.n, 1))


def extreme_comparability_penalty(poset: Poset, lower: float = 0.1, upper: float = 0.9) -> float:
    """只惩罚极端尾部，避免围绕某个人工中心值强行塑形。"""
    comp = comparable_fraction(poset)
    if comp < lower:
        return float((lower - comp) ** 2)
    if comp > upper:
        return float((comp - upper) ** 2)
    return 0.0


def neutral_penalty(poset: Poset) -> float:
    return (
        1.0 * normalized_degree_variance(poset)
        + 5.0 * layer_imbalance(poset)
        + 2.0 * extreme_comparability_penalty(poset)
    )
