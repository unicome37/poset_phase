from __future__ import annotations

import numpy as np
import pandas as pd

from generators import Poset
from observables import comparable_fraction, normalized_degree_variance
from observables_geo import (
    dimension_penalty_from_order_fraction,
    height_ratio,
    interval_empty_fraction,
    sample_intervals,
    mean_interval_size_ratio,
    width_ratio,
)


SIGNATURE_COLUMNS = [
    "sig_comp",
    "sig_d_eff",
    "sig_height_ratio",
    "sig_width_ratio",
    "sig_degree_var",
    "sig_interval_empty",
    "sig_interval_size_ratio",
]

GC_COMPONENT_COLUMNS = [
    "gc_var_comp",
    "gc_var_d_eff",
    "gc_var_height_ratio",
    "gc_var_width_ratio",
    "gc_var_size_ratio",
]


def signature_dict(poset: Poset) -> dict[str, float]:
    _, d_eff = dimension_penalty_from_order_fraction(poset)
    return {
        "sig_comp": comparable_fraction(poset),
        "sig_d_eff": d_eff,
        "sig_height_ratio": height_ratio(poset),
        "sig_width_ratio": width_ratio(poset),
        "sig_degree_var": normalized_degree_variance(poset),
        "sig_interval_empty": interval_empty_fraction(poset),
        "sig_interval_size_ratio": mean_interval_size_ratio(poset),
    }


def signature_vector_from_dict(signature: dict[str, float]) -> np.ndarray:
    return np.asarray([float(signature[col]) for col in SIGNATURE_COLUMNS], dtype=float)


def self_drift(before: dict[str, float], after: dict[str, float]) -> float:
    delta = signature_vector_from_dict(after) - signature_vector_from_dict(before)
    return float(np.mean(delta * delta))


def family_centroids(df: pd.DataFrame, group_col: str = "n") -> dict[tuple[int, str], np.ndarray]:
    centroids: dict[tuple[int, str], np.ndarray] = {}
    grouped = df.groupby([group_col, "family"])[SIGNATURE_COLUMNS].mean()
    for (n_value, family), row in grouped.iterrows():
        centroids[(int(n_value), str(family))] = row.to_numpy(dtype=float)
    return centroids


def nearest_family(
    signature: dict[str, float],
    centroids: dict[tuple[int, str], np.ndarray],
    n_value: int,
) -> tuple[str, float]:
    vec = signature_vector_from_dict(signature)
    candidates = [(family, centroid) for (n_ref, family), centroid in centroids.items() if n_ref == n_value]
    if not candidates:
        raise ValueError(f"No family centroids found for n={n_value}")

    best_family = ""
    best_dist = float("inf")
    for family, centroid in candidates:
        dist = float(np.linalg.norm(vec - centroid))
        if dist < best_dist:
            best_family = family
            best_dist = dist
    return best_family, best_dist


def coarse_grain_penalty(
    drift_value: float,
    family_switch_penalty: float,
    rank_shift_penalty: float,
    w_drift: float = 10.0,
    w_family: float = 1.5,
    w_rank: float = 0.5,
) -> float:
    return (
        w_drift * float(drift_value)
        + w_family * float(family_switch_penalty)
        + w_rank * float(rank_shift_penalty)
    )


def _subposet_signature_from_indices(poset: Poset, indices: np.ndarray) -> dict[str, float] | None:
    if len(indices) < 3:
        return None
    closure = poset.closure[np.ix_(indices, indices)].copy()
    np.fill_diagonal(closure, False)
    sub = Poset(closure)
    _, d_eff = dimension_penalty_from_order_fraction(sub)
    return {
        "comp": comparable_fraction(sub),
        "d_eff": d_eff,
        "height_ratio": height_ratio(sub),
        "width_ratio": width_ratio(sub),
        "size_ratio": float(len(indices) / max(poset.n, 1)),
    }


def local_interval_signatures(
    poset: Poset,
    max_pairs: int = 32,
) -> list[dict[str, float]]:
    signatures: list[dict[str, float]] = []
    for interior in sample_intervals(poset, max_pairs=max_pairs):
        sig = _subposet_signature_from_indices(poset, interior)
        if sig is not None:
            signatures.append(sig)
    return signatures


def global_consistency_penalty(
    poset: Poset,
    max_pairs: int = 32,
    w_comp: float = 4.0,
    w_dim: float = 5.0,
    w_height: float = 2.0,
    w_width: float = 2.0,
    w_size: float = 1.0,
) -> tuple[float, dict[str, float]]:
    sigs = local_interval_signatures(poset, max_pairs=max_pairs)
    if len(sigs) < 2:
        components = {
            "gc_var_comp": 0.0,
            "gc_var_d_eff": 0.0,
            "gc_var_height_ratio": 0.0,
            "gc_var_width_ratio": 0.0,
            "gc_var_size_ratio": 0.0,
            "gc_count": float(len(sigs)),
            "gc_total": 0.0,
        }
        return 0.0, components

    arr = {key: np.asarray([sig[key] for sig in sigs], dtype=float) for key in sigs[0]}
    var_comp = float(np.var(arr["comp"]))
    var_dim = float(np.var(arr["d_eff"]))
    var_height = float(np.var(arr["height_ratio"]))
    var_width = float(np.var(arr["width_ratio"]))
    var_size = float(np.var(arr["size_ratio"]))
    total = (
        w_comp * var_comp
        + w_dim * var_dim
        + w_height * var_height
        + w_width * var_width
        + w_size * var_size
    )
    components = {
        "gc_var_comp": var_comp,
        "gc_var_d_eff": var_dim,
        "gc_var_height_ratio": var_height,
        "gc_var_width_ratio": var_width,
        "gc_var_size_ratio": var_size,
        "gc_count": float(len(sigs)),
        "gc_total": float(total),
    }
    return float(total), components


def consistency_reference_windows(
    posets_by_n: dict[int, list[Poset]],
    max_pairs: int = 32,
    sigma_scale: float = 1.0,
) -> dict[int, dict[str, tuple[float, float]]]:
    """Build per-n reference windows from a reference family, typically lorentzian_like_2d."""
    windows: dict[int, dict[str, tuple[float, float]]] = {}
    for n_value, posets in posets_by_n.items():
        component_rows = []
        for poset in posets:
            _, comps = global_consistency_penalty(poset, max_pairs=max_pairs)
            component_rows.append(comps)
        if not component_rows:
            continue
        frame = pd.DataFrame(component_rows)
        per_n: dict[str, tuple[float, float]] = {}
        for col in GC_COMPONENT_COLUMNS:
            mean = float(frame[col].mean())
            std = float(frame[col].std(ddof=0))
            lower = max(0.0, mean - sigma_scale * std)
            upper = mean + sigma_scale * std
            per_n[col] = (lower, upper)
        windows[int(n_value)] = per_n
    return windows


def _window_penalty(value: float, lower: float, upper: float) -> float:
    if value < lower:
        return float((lower - value) ** 2)
    if value > upper:
        return float((value - upper) ** 2)
    return 0.0


def global_consistency_penalty_with_reference(
    poset: Poset,
    reference_windows: dict[int, dict[str, tuple[float, float]]] | None,
    max_pairs: int = 32,
    weights: dict[str, float] | None = None,
) -> tuple[float, dict[str, float]]:
    total_raw, components = global_consistency_penalty(poset, max_pairs=max_pairs)
    if not reference_windows or poset.n not in reference_windows:
        return total_raw, components

    weights = weights or {
        "gc_var_comp": 1.5,
        "gc_var_d_eff": 0.5,
        "gc_var_height_ratio": 1.0,
        "gc_var_width_ratio": 1.0,
        "gc_var_size_ratio": 0.5,
    }
    windows = reference_windows[poset.n]
    penalties = {}
    total = 0.0
    for col in GC_COMPONENT_COLUMNS:
        if col not in windows or col not in components:
            penalties[f"{col}_ref_penalty"] = 0.0
            continue
        lower, upper = windows[col]
        penalty = weights.get(col, 1.0) * _window_penalty(float(components[col]), lower, upper)
        penalties[f"{col}_ref_penalty"] = float(penalty)
        total += penalty
    components = dict(components)
    components.update(penalties)
    components["gc_total_ref"] = float(total)
    return float(total), components
