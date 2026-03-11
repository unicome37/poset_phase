from __future__ import annotations

import numpy as np

from generators import Poset
from observables import comparable_fraction, layer_profile


R2_REF = 0.5009


def height_ratio(poset: Poset) -> float:
    """Number of layers normalized by n."""
    layers = layer_profile(poset)
    return float(len(layers) / max(poset.n, 1))


def width_ratio(poset: Poset) -> float:
    """Max layer occupancy normalized by n."""
    layers = layer_profile(poset)
    return float(layers.max() / max(poset.n, 1))


def width_height_balance_penalty(poset: Poset, target_ratio: float = 1.0) -> float:
    """Reward width/height balance instead of extreme flattening or needle-like towers."""
    layers = layer_profile(poset)
    width = float(layers.max())
    height = float(len(layers))
    ratio = width / max(height, 1.0)
    return float((np.log(max(ratio, 1e-8)) - np.log(target_ratio)) ** 2)


def window_penalty(value: float, lower: float, upper: float) -> float:
    """Zero inside a target window, quadratic penalty outside."""
    if value < lower:
        return float((lower - value) ** 2)
    if value > upper:
        return float((value - upper) ** 2)
    return 0.0


def comparability_window_penalty(
    poset: Poset,
    lower: float = 0.10,
    upper: float = 0.38,
) -> float:
    """Geometric proxy via order-fraction window.

    The aim is not to force a single target value, but to suppress
    structures that are too chain-like or too antichain-like.
    """
    comp = comparable_fraction(poset)
    return window_penalty(comp, lower=lower, upper=upper)


def estimate_dimension_proxy_from_order_fraction(r: float) -> float:
    """Small-N proxy inspired by the review note's 2D Minkowski reference design."""
    r = max(1e-6, min(1.0 - 1e-6, r))
    delta = R2_REF - r
    d_eff = 2.0 + 6.0 * delta
    return float(np.clip(d_eff, 0.5, 8.0))


def dimension_penalty_from_order_fraction(poset: Poset, d_target: float = 2.0) -> tuple[float, float]:
    r = comparable_fraction(poset)
    d_eff = estimate_dimension_proxy_from_order_fraction(r)
    return float((d_eff - d_target) ** 2), float(d_eff)


def _comparable_fraction_from_closure(closure: np.ndarray) -> float:
    m = int(closure.shape[0])
    total = m * (m - 1) / 2
    if total <= 0:
        return 0.0
    return float(closure.sum() / total)


def dimension_consistency_penalty(
    poset: Poset,
    max_pairs: int = 64,
    min_interval_size: int = 4,
) -> tuple[float, float, float, float, int]:
    """Non-targeted dimensional self-consistency penalty.

    Instead of pulling the structure toward a fixed target dimension, this asks
    whether locally sampled intervals agree with each other and with the global
    order-fraction dimension proxy.
    """
    global_r = comparable_fraction(poset)
    global_d_eff = estimate_dimension_proxy_from_order_fraction(global_r)
    intervals = sample_intervals(poset, max_pairs=max_pairs)
    local_d_effs = []

    for interior in intervals:
        if len(interior) < min_interval_size:
            continue
        sub = poset.closure[np.ix_(interior, interior)]
        local_r = _comparable_fraction_from_closure(sub)
        local_d_effs.append(estimate_dimension_proxy_from_order_fraction(local_r))

    if not local_d_effs:
        return 0.0, float(global_d_eff), float(global_d_eff), 0.0, 0

    local = np.asarray(local_d_effs, dtype=float)
    mean_local = float(local.mean())
    var_local = float(local.var(ddof=0))
    mismatch = float((mean_local - global_d_eff) ** 2)
    total = var_local + mismatch
    return float(total), float(global_d_eff), mean_local, var_local, int(len(local))


def cover_density(poset: Poset) -> float:
    """Density of the transitive reduction (Hasse diagram) edges."""
    c = poset.closure
    n = poset.n
    cover_count = 0

    for i in range(n):
        succ = np.where(c[i])[0]
        for j in succ:
            is_cover = True
            mid = np.where(c[i] & c[:, j])[0]
            if len(mid) > 0:
                is_cover = False
            if is_cover:
                cover_count += 1

    total = n * (n - 1) / 2
    return float(cover_count / total) if total else 0.0


def cover_density_penalty(
    poset: Poset,
    lower: float = 0.03,
    upper: float = 0.20,
) -> float:
    """Suppress both over-connected and overly sparse Hasse diagrams."""
    return window_penalty(cover_density(poset), lower=lower, upper=upper)


def sample_intervals(
    poset: Poset,
    max_pairs: int = 64,
    seed: int = 0,
) -> list[np.ndarray]:
    """Sample intervals I(x, y) = {z | x < z < y} from comparable pairs."""
    c = poset.closure
    comparable_pairs = np.argwhere(c)
    if len(comparable_pairs) == 0:
        return []

    rng = np.random.default_rng(seed + poset.n)
    if len(comparable_pairs) > max_pairs:
        idx = rng.choice(len(comparable_pairs), size=max_pairs, replace=False)
        comparable_pairs = comparable_pairs[idx]

    intervals = []
    for x, y in comparable_pairs:
        interior = np.where(c[x] & c[:, y])[0]
        interior = interior[(interior != x) & (interior != y)]
        intervals.append(interior)
    return intervals


def _layer_profile_from_closure(closure: np.ndarray) -> np.ndarray:
    indeg = closure.sum(axis=0).astype(int)
    remaining = set(range(closure.shape[0]))
    layers = []

    while remaining:
        mins = [i for i in remaining if indeg[i] == 0]
        if not mins:
            break
        layers.append(len(mins))
        for u in mins:
            remaining.remove(u)
            for v in np.where(closure[u])[0]:
                indeg[v] -= 1

    return np.asarray(layers, dtype=float)


def interval_empty_fraction(poset: Poset, max_pairs: int = 64) -> float:
    intervals = sample_intervals(poset, max_pairs=max_pairs)
    if not intervals:
        return 1.0
    empty = sum(1 for interior in intervals if len(interior) == 0)
    return float(empty / len(intervals))


def mean_interval_size_ratio(poset: Poset, max_pairs: int = 64) -> float:
    intervals = sample_intervals(poset, max_pairs=max_pairs)
    if not intervals:
        return 0.0
    sizes = np.array([len(interior) for interior in intervals], dtype=float)
    return float(sizes.mean() / max(poset.n, 1))


def interval_profile_penalty(
    poset: Poset,
    empty_upper: float = 0.55,
    mean_size_lower: float = 0.04,
    mean_size_upper: float = 0.28,
) -> float:
    """Penalize interval statistics inconsistent with locally rich causal structure."""
    empty_frac = interval_empty_fraction(poset)
    mean_size_ratio = mean_interval_size_ratio(poset)
    return (
        4.0 * window_penalty(mean_size_ratio, lower=mean_size_lower, upper=mean_size_upper)
        + 5.0 * window_penalty(1.0 - empty_frac, lower=1.0 - empty_upper, upper=1.0)
    )


def interval_shape_penalty(
    poset: Poset,
    max_pairs: int = 64,
    min_interval_size: int = 4,
    height_ratio_lower: float = 0.30,
    peak_ratio_upper: float = 0.72,
) -> float:
    """Penalize blocky interval interiors that look like artificial layer stacks.

    For each sampled interval, build the induced sub-poset on its interior and
    inspect its own layer profile. Lorentzian-like intervals should typically
    have multiple interior layers and avoid a single dominant layer.
    """
    c = poset.closure
    intervals = sample_intervals(poset, max_pairs=max_pairs)
    if not intervals:
        return 0.0

    penalties = []
    for interior in intervals:
        if len(interior) < min_interval_size:
            continue
        sub = c[np.ix_(interior, interior)]
        layers = _layer_profile_from_closure(sub)
        if len(layers) == 0:
            continue
        size = float(len(interior))
        height_ratio = len(layers) / size
        peak_ratio = float(layers.max() / size)
        penalties.append(
            3.0 * window_penalty(height_ratio, lower=height_ratio_lower, upper=1.0)
            + 4.0 * window_penalty(peak_ratio, lower=0.0, upper=peak_ratio_upper)
        )

    if not penalties:
        return 0.0
    return float(np.mean(penalties))


def local_layer_smoothness_penalty(poset: Poset) -> float:
    """Penalize violent changes in adjacent layer sizes."""
    layers = layer_profile(poset).astype(float)
    if len(layers) <= 1:
        return 0.0
    layers = layers / layers.sum()
    diffs = np.diff(layers)
    return float((diffs * diffs).mean())


def geometric_penalty(poset: Poset) -> float:
    dim_penalty, _ = dimension_penalty_from_order_fraction(poset)
    return (
        2.0 * width_height_balance_penalty(poset)
        + 8.0 * dim_penalty
        + 6.0 * comparability_window_penalty(poset)
        + 3.0 * cover_density_penalty(poset)
        + 5.0 * interval_profile_penalty(poset)
        + 5.0 * interval_shape_penalty(poset)
        + 2.0 * local_layer_smoothness_penalty(poset)
    )


def geometric_components(poset: Poset) -> dict[str, float]:
    width_height = width_height_balance_penalty(poset)
    dim_penalty, d_eff = dimension_penalty_from_order_fraction(poset)
    dim_consistency, global_d_eff, local_mean_d_eff, local_var_d_eff, local_dim_count = dimension_consistency_penalty(poset)
    comparability = comparability_window_penalty(poset)
    cover = cover_density_penalty(poset)
    interval_profile = interval_profile_penalty(poset)
    interval_shape = interval_shape_penalty(poset)
    layer_smoothness = local_layer_smoothness_penalty(poset)
    return {
        "geo_width_height": width_height,
        "geo_dim_proxy_penalty": dim_penalty,
        "geo_dim_eff": d_eff,
        "geo_dim_consistency": dim_consistency,
        "geo_dim_consistency_global_d_eff": global_d_eff,
        "geo_dim_consistency_local_mean_d_eff": local_mean_d_eff,
        "geo_dim_consistency_local_var_d_eff": local_var_d_eff,
        "geo_dim_consistency_local_count": float(local_dim_count),
        "geo_comparability_window": comparability,
        "geo_cover_density": cover,
        "geo_interval_profile": interval_profile,
        "geo_interval_shape": interval_shape,
        "geo_layer_smoothness": layer_smoothness,
        "geo_total": (
            2.0 * width_height
            + 8.0 * dim_penalty
            + 6.0 * comparability
            + 3.0 * cover
            + 5.0 * interval_profile
            + 5.0 * interval_shape
            + 2.0 * layer_smoothness
        ),
    }
