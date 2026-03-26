from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Poset:
    """用传递闭包表示有限偏序。"""

    closure: np.ndarray  # shape (n, n), dtype=bool

    @property
    def n(self) -> int:
        return self.closure.shape[0]


def transitive_closure(adj: np.ndarray) -> np.ndarray:
    """Warshall-style 布尔传递闭包。"""
    c = adj.copy()
    n = c.shape[0]
    for k in range(n):
        c |= c[:, [k]] & c[[k], :]
    np.fill_diagonal(c, False)
    return c


def _random_layer_split(n: int, n_layers: int, rng: np.random.Generator) -> list[np.ndarray]:
    cuts = np.sort(rng.choice(np.arange(1, n), size=n_layers - 1, replace=False))
    return list(np.split(np.arange(n), cuts))


def _split_with_sizes(indices: np.ndarray, sizes: list[int]) -> list[np.ndarray]:
    cuts = np.cumsum(sizes[:-1])
    return list(np.split(indices, cuts))


def _layer_sizes_by_mode(
    n: int,
    n_layers: int,
    mode: str,
    rng: np.random.Generator,
) -> list[int]:
    if n_layers > n:
        raise ValueError(f"n_layers={n_layers} cannot exceed n={n}")

    if mode == "uniform":
        base = [1] * n_layers
        for _ in range(n - n_layers):
            base[int(rng.integers(0, n_layers))] += 1
        return base

    if mode == "tapered":
        weights = np.linspace(n_layers, 1.0, n_layers)
    elif mode == "middle_heavy":
        center = (n_layers - 1) / 2.0
        idx = np.arange(n_layers)
        weights = np.maximum(1.0, n_layers - np.abs(idx - center) * 2.0)
    else:
        raise ValueError(f"Unsupported imbalance mode: {mode}")

    sizes = np.ones(n_layers, dtype=int)
    remaining = n - n_layers
    if remaining > 0:
        probs = weights / weights.sum()
        sizes += rng.multinomial(remaining, probs)
    return sizes.tolist()


def generate_absolute_layered(
    n: int,
    n_layers: int = 4,
    seed: int | None = None,
) -> Poset:
    rng = np.random.default_rng(seed)
    layers = _random_layer_split(n, n_layers, rng)
    adj = np.zeros((n, n), dtype=bool)
    for i in range(len(layers)):
        for j in range(i + 1, len(layers)):
            adj[np.ix_(layers[i], layers[j])] = True
    return Poset(transitive_closure(adj))


def generate_kr_like(n: int, p: float = 0.5, seed: int | None = None) -> Poset:
    rng = np.random.default_rng(seed)
    sizes = [n // 4, n // 2, n - (n // 4 + n // 2)]
    idx = np.arange(n)
    l1, l2, l3 = np.split(idx, [sizes[0], sizes[0] + sizes[1]])

    adj = np.zeros((n, n), dtype=bool)
    adj[np.ix_(l1, l2)] = rng.random((len(l1), len(l2))) < p
    adj[np.ix_(l2, l3)] = rng.random((len(l2), len(l3))) < p
    return Poset(transitive_closure(adj))


def generate_transitive_percolation(
    n: int,
    p: float = 0.08,
    seed: int | None = None,
) -> Poset:
    rng = np.random.default_rng(seed)
    adj = np.triu(rng.random((n, n)) < p, k=1)
    return Poset(transitive_closure(adj))


def generate_lorentzian_like_2d(n: int, seed: int | None = None) -> Poset:
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random(n)

    dt = t[None, :] - t[:, None]
    dx = np.abs(x[None, :] - x[:, None])
    adj = (dt > 0.0) & (dt >= dx)
    return Poset(transitive_closure(adj))


def generate_lorentzian_like_3d(n: int, seed: int | None = None) -> Poset:
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random(n)
    y = rng.random(n)

    dt = t[None, :] - t[:, None]
    dx2 = (x[None, :] - x[:, None]) ** 2
    dy2 = (y[None, :] - y[:, None]) ** 2
    adj = (dt > 0.0) & ((dt * dt) >= (dx2 + dy2))
    return Poset(transitive_closure(adj))


def generate_lorentzian_like_4d(n: int, seed: int | None = None) -> Poset:
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random(n)
    y = rng.random(n)
    z = rng.random(n)

    dt = t[None, :] - t[:, None]
    dx2 = (x[None, :] - x[:, None]) ** 2
    dy2 = (y[None, :] - y[:, None]) ** 2
    dz2 = (z[None, :] - z[:, None]) ** 2
    adj = (dt > 0.0) & ((dt * dt) >= (dx2 + dy2 + dz2))
    return Poset(transitive_closure(adj))


def generate_lorentzian_like_5d(n: int, seed: int | None = None) -> Poset:
    """5D Lorentzian-like poset: 1 time + 4 spatial dimensions."""
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random(n)
    y = rng.random(n)
    z = rng.random(n)
    w = rng.random(n)

    dt = t[None, :] - t[:, None]
    dx2 = (x[None, :] - x[:, None]) ** 2
    dy2 = (y[None, :] - y[:, None]) ** 2
    dz2 = (z[None, :] - z[:, None]) ** 2
    dw2 = (w[None, :] - w[:, None]) ** 2
    adj = (dt > 0.0) & ((dt * dt) >= (dx2 + dy2 + dz2 + dw2))
    return Poset(transitive_closure(adj))


def generate_lorentzian_like_6d(n: int, seed: int | None = None) -> Poset:
    """6D Lorentzian-like poset: 1 time + 5 spatial dimensions."""
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random(n)
    y = rng.random(n)
    z = rng.random(n)
    w = rng.random(n)
    v = rng.random(n)

    dt = t[None, :] - t[:, None]
    dx2 = (x[None, :] - x[:, None]) ** 2
    dy2 = (y[None, :] - y[:, None]) ** 2
    dz2 = (z[None, :] - z[:, None]) ** 2
    dw2 = (w[None, :] - w[:, None]) ** 2
    dv2 = (v[None, :] - v[:, None]) ** 2
    adj = (dt > 0.0) & ((dt * dt) >= (dx2 + dy2 + dz2 + dw2 + dv2))
    return Poset(transitive_closure(adj))


def generate_lorentzian_like_7d(n: int, seed: int | None = None) -> Poset:
    """7D Lorentzian-like poset: 1 time + 6 spatial dimensions."""
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random(n)
    y = rng.random(n)
    z = rng.random(n)
    w = rng.random(n)
    v = rng.random(n)
    u = rng.random(n)

    dt = t[None, :] - t[:, None]
    dx2 = (x[None, :] - x[:, None]) ** 2
    dy2 = (y[None, :] - y[:, None]) ** 2
    dz2 = (z[None, :] - z[:, None]) ** 2
    dw2 = (w[None, :] - w[:, None]) ** 2
    dv2 = (v[None, :] - v[:, None]) ** 2
    du2 = (u[None, :] - u[:, None]) ** 2
    adj = (dt > 0.0) & ((dt * dt) >= (dx2 + dy2 + dz2 + dw2 + dv2 + du2))
    return Poset(transitive_closure(adj))


def generate_interval_order(
    n: int,
    mean_length: float = 0.3,
    seed: int | None = None,
) -> Poset:
    rng = np.random.default_rng(seed)
    starts = rng.random(n)
    lengths = rng.exponential(mean_length, size=n)
    ends = starts + lengths

    adj = ends[:, None] < starts[None, :]
    return Poset(transitive_closure(adj))


def generate_multi_layer_random(
    n: int,
    n_layers: int = 5,
    p: float = 0.3,
    seed: int | None = None,
) -> Poset:
    rng = np.random.default_rng(seed)
    layers = _random_layer_split(n, n_layers, rng)

    adj = np.zeros((n, n), dtype=bool)
    for i in range(len(layers)):
        for j in range(i + 1, len(layers)):
            mask = rng.random((len(layers[i]), len(layers[j]))) < p
            adj[np.ix_(layers[i], layers[j])] = mask
    return Poset(transitive_closure(adj))


def generate_random_layered(
    n: int,
    n_layers: int = 6,
    imbalance_mode: str = "uniform",
    adjacent_p: float = 0.35,
    skip_p: float = 0.08,
    seed: int | None = None,
) -> Poset:
    rng = np.random.default_rng(seed)
    permuted = rng.permutation(n)
    sizes = _layer_sizes_by_mode(n, n_layers, imbalance_mode, rng)
    layers = _split_with_sizes(permuted, sizes)

    adj = np.zeros((n, n), dtype=bool)
    for i in range(len(layers)):
        for j in range(i + 1, len(layers)):
            p = adjacent_p if j == i + 1 else skip_p
            mask = rng.random((len(layers[i]), len(layers[j]))) < p
            adj[np.ix_(layers[i], layers[j])] = mask
    return Poset(transitive_closure(adj))


def generate_random_layered_k4_uniform(n: int, seed: int | None = None) -> Poset:
    return generate_random_layered(n=n, n_layers=4, imbalance_mode="uniform", adjacent_p=0.38, skip_p=0.06, seed=seed)


def generate_random_layered_k6_uniform(n: int, seed: int | None = None) -> Poset:
    return generate_random_layered(n=n, n_layers=6, imbalance_mode="uniform", adjacent_p=0.34, skip_p=0.08, seed=seed)


def generate_random_layered_k8_uniform(n: int, seed: int | None = None) -> Poset:
    return generate_random_layered(n=n, n_layers=8, imbalance_mode="uniform", adjacent_p=0.30, skip_p=0.10, seed=seed)


def generate_random_layered_k6_tapered(n: int, seed: int | None = None) -> Poset:
    return generate_random_layered(n=n, n_layers=6, imbalance_mode="tapered", adjacent_p=0.34, skip_p=0.08, seed=seed)


def generate_random_layered_k6_middle_heavy(n: int, seed: int | None = None) -> Poset:
    return generate_random_layered(n=n, n_layers=6, imbalance_mode="middle_heavy", adjacent_p=0.34, skip_p=0.08, seed=seed)


def generate_random_layered_k6_longjump(n: int, seed: int | None = None) -> Poset:
    return generate_random_layered(n=n, n_layers=6, imbalance_mode="uniform", adjacent_p=0.32, skip_p=0.16, seed=seed)


def generate_kr_2layer(n: int, p: float = 0.5, seed: int | None = None) -> Poset:
    """2-layer bipartite poset (KR control): bottom ~n/4, top ~3n/4."""
    rng = np.random.default_rng(seed)
    n_bot = max(1, n // 4)
    n_top = n - n_bot
    idx = rng.permutation(n)
    bot, top = idx[:n_bot], idx[n_bot:]
    adj = np.zeros((n, n), dtype=bool)
    adj[np.ix_(bot, top)] = rng.random((n_bot, n_top)) < p
    return Poset(transitive_closure(adj))


def generate_kr_4layer(n: int, p: float = 0.5, seed: int | None = None) -> Poset:
    """4-layer KR-extended poset: ratio ~1:3:3:1, adjacent-layer edges only."""
    rng = np.random.default_rng(seed)
    s1 = max(1, n // 8)
    s4 = max(1, n // 8)
    mid = n - s1 - s4
    s2 = mid // 2
    s3 = mid - s2
    sizes = [s1, s2, s3, s4]
    idx = rng.permutation(n)
    layers = _split_with_sizes(idx, sizes)
    adj = np.zeros((n, n), dtype=bool)
    for i in range(len(layers) - 1):
        mask = rng.random((len(layers[i]), len(layers[i + 1]))) < p
        adj[np.ix_(layers[i], layers[i + 1])] = mask
    return Poset(transitive_closure(adj))
