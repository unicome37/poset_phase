from __future__ import annotations

from collections.abc import Mapping

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


def estimate_dimension_proxy_from_chain_depth(n: int, height_layers: float) -> float:
    """Rough dimension proxy from chain-depth scaling."""
    n = max(int(n), 2)
    h = max(float(height_layers), 1.000001)
    d_eff = np.log(n) / np.log(h)
    return float(np.clip(d_eff, 0.5, 8.0))


def estimate_dimension_proxy_from_width_scaling(n: int, width_layers: float) -> float:
    """Rough dimension proxy from width scaling."""
    n = max(int(n), 2)
    w = min(max(float(width_layers), 1.000001), float(n) - 1e-6)
    alpha = np.log(w) / np.log(n)
    alpha = float(np.clip(alpha, 1e-6, 1.0 - 1e-6))
    d_eff = 1.0 / (1.0 - alpha)
    return float(np.clip(d_eff, 0.5, 8.0))


def dimension_penalty_from_order_fraction(poset: Poset, d_target: float = 2.0) -> tuple[float, float]:
    r = comparable_fraction(poset)
    d_eff = estimate_dimension_proxy_from_order_fraction(r)
    return float((d_eff - d_target) ** 2), float(d_eff)


def dimension_proxy_views(poset: Poset) -> dict[str, float]:
    layers = layer_profile(poset)
    height_layers = float(len(layers)) if len(layers) else 1.0
    width_layers = float(layers.max()) if len(layers) else 1.0
    comp = comparable_fraction(poset)
    return {
        "d_order": estimate_dimension_proxy_from_order_fraction(comp),
        "d_chain": estimate_dimension_proxy_from_chain_depth(poset.n, height_layers),
        "d_width": estimate_dimension_proxy_from_width_scaling(poset.n, width_layers),
    }


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


def multi_estimator_dimension_consistency_penalty(
    poset: Poset,
    max_pairs: int = 64,
    min_interval_size: int = 4,
) -> tuple[float, float, float, float, int]:
    """Strengthen dim_consistency by requiring cross-estimator agreement."""
    views = dimension_proxy_views(poset)
    global_values = np.asarray(list(views.values()), dtype=float)
    cross_view_var = float(global_values.var(ddof=0))

    base_total, global_d_eff, mean_local, var_local, local_count = dimension_consistency_penalty(
        poset, max_pairs=max_pairs, min_interval_size=min_interval_size
    )
    total = float(base_total + cross_view_var)
    return total, float(global_d_eff), float(mean_local), float(var_local + cross_view_var), int(local_count)


def cover_density(poset: Poset) -> float:
    """Density of the transitive reduction (Hasse diagram) edges."""
    c = poset.closure.astype(np.uint8, copy=False)
    n = poset.n
    has_intermediate = (c @ c).astype(bool, copy=False)
    cover_mask = poset.closure & ~has_intermediate
    np.fill_diagonal(cover_mask, False)
    total = n * (n - 1) / 2
    return float(cover_mask.sum() / total) if total else 0.0


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


DEFAULT_GEOMETRIC_WEIGHTS: dict[str, float] = {
    "geo_width_height": 2.0,
    "geo_dim_proxy_penalty": 8.0,
    "geo_comparability_window": 6.0,
    "geo_cover_density": 3.0,
    "geo_interval_profile": 5.0,
    "geo_interval_shape": 5.0,
    "geo_layer_smoothness": 2.0,
}


def weighted_geometric_total(
    components: Mapping[str, float],
    weights: Mapping[str, float] | None = None,
) -> float:
    active_weights = DEFAULT_GEOMETRIC_WEIGHTS if weights is None else weights
    return float(sum(active_weights.get(key, 0.0) * components[key] for key in DEFAULT_GEOMETRIC_WEIGHTS))


def geometric_penalty(poset: Poset, weights: Mapping[str, float] | None = None) -> float:
    components = geometric_components(poset)
    return weighted_geometric_total(components, weights=weights)


def geometric_penalty_from_components(components: Mapping[str, float], weights: Mapping[str, float] | None = None) -> float:
    return weighted_geometric_total(components, weights=weights)


def geometric_components(poset: Poset) -> dict[str, float]:
    dim_penalty, _ = dimension_penalty_from_order_fraction(poset)
    width_height = width_height_balance_penalty(poset)
    dim_penalty, d_eff = dimension_penalty_from_order_fraction(poset)
    dim_consistency, global_d_eff, local_mean_d_eff, local_var_d_eff, local_dim_count = dimension_consistency_penalty(poset)
    dim_multi_consistency, _, _, multi_var_d_eff, _ = multi_estimator_dimension_consistency_penalty(poset)
    views = dimension_proxy_views(poset)
    comparability = comparability_window_penalty(poset)
    cover = cover_density_penalty(poset)
    interval_profile = interval_profile_penalty(poset)
    interval_shape = interval_shape_penalty(poset)
    layer_smoothness = local_layer_smoothness_penalty(poset)
    components = {
        "geo_width_height": width_height,
        "geo_dim_proxy_penalty": dim_penalty,
        "geo_dim_eff": d_eff,
        "geo_dim_consistency": dim_consistency,
        "geo_dim_consistency_global_d_eff": global_d_eff,
        "geo_dim_consistency_local_mean_d_eff": local_mean_d_eff,
        "geo_dim_consistency_local_var_d_eff": local_var_d_eff,
        "geo_dim_consistency_local_count": float(local_dim_count),
        "geo_dim_multi_consistency": dim_multi_consistency,
        "geo_dim_multi_consistency_var": multi_var_d_eff,
        "geo_d_order": views["d_order"],
        "geo_d_chain": views["d_chain"],
        "geo_d_width": views["d_width"],
        "geo_comparability_window": comparability,
        "geo_cover_density": cover,
        "geo_interval_profile": interval_profile,
        "geo_interval_shape": interval_shape,
        "geo_layer_smoothness": layer_smoothness,
    }
    components["geo_total"] = weighted_geometric_total(components)
    return components
