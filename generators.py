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
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p:
                adj[i, j] = True
    return Poset(transitive_closure(adj))


def generate_lorentzian_like_2d(n: int, seed: int | None = None) -> Poset:
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random(n)

    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(n):
            dt = t[j] - t[i]
            if dt > 0 and dt >= abs(x[j] - x[i]):
                adj[i, j] = True
    return Poset(transitive_closure(adj))


def generate_lorentzian_like_3d(n: int, seed: int | None = None) -> Poset:
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random(n)
    y = rng.random(n)

    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(n):
            dt = t[j] - t[i]
            if dt > 0 and dt * dt >= (x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2:
                adj[i, j] = True
    return Poset(transitive_closure(adj))


def generate_lorentzian_like_4d(n: int, seed: int | None = None) -> Poset:
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random(n)
    y = rng.random(n)
    z = rng.random(n)

    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(n):
            dt = t[j] - t[i]
            if dt > 0 and dt * dt >= (
                (x[j] - x[i]) ** 2
                + (y[j] - y[i]) ** 2
                + (z[j] - z[i]) ** 2
            ):
                adj[i, j] = True
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

    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(n):
            if ends[i] < starts[j]:
                adj[i, j] = True
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
