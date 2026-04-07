"""Pure information-theoretic observables for poset phase experiment.

All penalties are derived from Shannon / von Neumann entropy measures
applied to the poset's combinatorial structure.  **No geometric targets**
(dimension, comparability windows, Minkowski reference values ...) appear
anywhere in this module.

Design principle
----------------
Each penalty detects *structural degeneracy* (low information content)
rather than *deviation from a geometric ideal*.

* Chains, antichains, hub-dominated graphs  → high penalty
* Manifold-like, well-organised graphs      → low penalty

This allows testing whether the action  A = −β·logH + γ·I_info
still selects Lorentzian-like posets as the dominant phase even when
the penalty encodes zero geometric knowledge — thereby defeating
the "circular argument" criticism.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from generators import Poset
from observables import layer_profile


# =========================================================================
# Internal helpers (self-contained — no imports from observables_geo)
# =========================================================================

def _cover_relation(poset: Poset) -> np.ndarray:
    """Boolean matrix of the Hasse diagram (transitive reduction)."""
    c = poset.closure.astype(np.uint8, copy=False)
    has_intermediate = (c @ c).astype(bool, copy=False)
    cover = poset.closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return cover


def _cover_density(poset: Poset) -> float:
    """Fraction of directed cover edges over n(n−1)/2."""
    cover = _cover_relation(poset)
    total = poset.n * (poset.n - 1) / 2
    return float(cover.sum() / total) if total > 0 else 0.0


def _sample_intervals(
    poset: Poset,
    max_pairs: int = 64,
    seed: int = 0,
) -> list[np.ndarray]:
    """Sample causal intervals I(x,y) = {z : x < z < y}."""
    c = poset.closure
    pairs = np.argwhere(c)
    if len(pairs) == 0:
        return []
    rng = np.random.default_rng(seed + poset.n)
    if len(pairs) > max_pairs:
        idx = rng.choice(len(pairs), size=max_pairs, replace=False)
        pairs = pairs[idx]
    intervals: list[np.ndarray] = []
    for x, y in pairs:
        interior = np.where(c[x] & c[:, y])[0]
        interior = interior[(interior != x) & (interior != y)]
        intervals.append(interior)
    return intervals


def _shannon_entropy(probs: np.ndarray) -> float:
    """H = −Σ pᵢ ln pᵢ (natural log, zero-safe)."""
    p = probs[probs > 0]
    return float(-np.sum(p * np.log(p)))


def _binary_entropy(p: float) -> float:
    """H_b(p) = −p ln p − (1−p) ln(1−p)."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return float(-p * np.log(p) - (1 - p) * np.log(1 - p))


# =========================================================================
# Raw information-theoretic measures  (not penalties, for diagnostics)
# =========================================================================

def von_neumann_graph_entropy(poset: Poset) -> float:
    r"""Von Neumann entropy of the normalised Laplacian of the Hasse diagram.

    .. math::
        S_{vN} = -\sum_i \frac{\lambda_i}{\mathrm{tr}\,L}
                       \ln \frac{\lambda_i}{\mathrm{tr}\,L}

    Computed from the **undirected** Hasse graph  L = D − A.
    """
    cover = _cover_relation(poset)
    adj = (cover | cover.T).astype(float)
    deg = adj.sum(axis=1)
    laplacian = np.diag(deg) - adj

    eigs = np.linalg.eigvalsh(laplacian)
    eigs = eigs[eigs > 1e-12]  # discard zero eigenvalues
    if len(eigs) == 0:
        return 0.0
    total = eigs.sum()
    probs = eigs / total
    return _shannon_entropy(probs)


def degree_gini(poset: Poset) -> float:
    """Gini coefficient of the total (in+out) degree sequence from closure.

    Self-loops (diagonal) are excluded to avoid inflating all degrees
    by a constant.

    * G = 0 → perfectly uniform degrees
    * G → 1 → extreme hub concentration (star-like)
    """
    c = poset.closure.copy()
    np.fill_diagonal(c, False)
    c = c.astype(float)
    total_deg = c.sum(axis=1) + c.sum(axis=0)
    n = len(total_deg)
    if n <= 1:
        return 0.0
    sorted_deg = np.sort(total_deg)
    denom = float(n * sorted_deg.sum())
    if denom <= 0:
        return 0.0
    numerator = 2.0 * np.sum(np.arange(1, n + 1) * sorted_deg)
    gini = numerator / denom - (n + 1.0) / n
    return float(np.clip(gini, 0.0, 1.0))


def layer_size_entropy_normalised(poset: Poset) -> float:
    """Normalised Shannon entropy of layer sizes:  H / ln K  ∈ [0, 1].

    * 1.0  →  all layers have equal size
    * → 0  →  extreme concentration in one or two layers
    """
    profile = layer_profile(poset).astype(float)
    K = len(profile)
    if K <= 1:
        return 0.0
    probs = profile / profile.sum()
    H = _shannon_entropy(probs)
    H_max = np.log(K)
    return float(np.clip(H / H_max, 0.0, 1.0)) if H_max > 0 else 0.0


def cover_binary_entropy_normalised(poset: Poset) -> float:
    """Normalised binary entropy of cover density:  H_b(p) / ln 2  ∈ [0, 1].

    * 1.0  →  cover density exactly 0.5 (maximally informative edge set)
    * → 0  →  extremely sparse or dense Hasse diagram
    """
    p = np.clip(_cover_density(poset), 0.0, 1.0)
    H = _binary_entropy(p)
    return float(H / np.log(2))


def interval_size_entropy_normalised(poset: Poset, max_pairs: int = 64) -> float:
    """Normalised Shannon entropy of the interval-size distribution.

    * High  →  rich variety of interval sizes (manifold-like causal structure)
    * Low   →  monotonous / degenerate intervals
    """
    intervals = _sample_intervals(poset, max_pairs=max_pairs)
    if not intervals:
        return 0.0
    sizes = np.array([len(interior) for interior in intervals])
    max_size = int(sizes.max())
    if max_size <= 0:
        return 0.0
    bins = np.bincount(sizes, minlength=max_size + 1).astype(float)
    probs = bins / bins.sum()
    H = _shannon_entropy(probs)
    H_max = np.log(max_size + 1)
    return float(np.clip(H / H_max, 0.0, 1.0)) if H_max > 0 else 0.0


# =========================================================================
# Penalty functions  (higher = structurally worse)
# =========================================================================

def spectral_entropy_deficit_penalty(poset: Poset) -> float:
    r"""$(1 - S_{vN} / \ln n)^2$.  Penalises spectrally degenerate structures."""
    n = poset.n
    if n <= 2:
        return 0.0
    S = von_neumann_graph_entropy(poset)
    S_max = np.log(n)
    deficit = 1.0 - S / S_max if S_max > 0 else 1.0
    return float(deficit ** 2)


def degree_heterogeneity_penalty(poset: Poset) -> float:
    r"""$G^2_{\mathrm{degree}}$.  Penalises hub-dominated / star-like structures."""
    return float(degree_gini(poset) ** 2)


def layer_concentration_penalty(poset: Poset) -> float:
    r"""$(1 - H_{\mathrm{layer}} / \ln K)^2$.  Penalises extreme layer imbalance."""
    return float((1.0 - layer_size_entropy_normalised(poset)) ** 2)


def edge_density_extremity_penalty(poset: Poset) -> float:
    r"""$(1 - H_b(p) / \ln 2)^2$.  Penalises extreme cover density."""
    return float((1.0 - cover_binary_entropy_normalised(poset)) ** 2)


def interval_diversity_deficit_penalty(poset: Poset, max_pairs: int = 64) -> float:
    r"""$(1 - H_{\mathrm{interval}} / \ln m)^2$.  Penalises monotonous intervals."""
    return float((1.0 - interval_size_entropy_normalised(poset, max_pairs)) ** 2)


# =========================================================================
# Composite system  (mirrors geometric_components / geometric_penalty)
# =========================================================================

DEFAULT_INFO_WEIGHTS: dict[str, float] = {
    "info_spectral_entropy_deficit": 5.0,
    "info_degree_heterogeneity": 5.0,
    "info_layer_concentration": 5.0,
    "info_edge_density_extremity": 5.0,
    "info_interval_diversity_deficit": 5.0,
}


def weighted_info_total(
    components: Mapping[str, float],
    weights: Mapping[str, float] | None = None,
) -> float:
    w = DEFAULT_INFO_WEIGHTS if weights is None else weights
    return float(sum(w.get(k, 0.0) * components[k] for k in DEFAULT_INFO_WEIGHTS))


def info_components(poset: Poset) -> dict[str, float]:
    """Compute all information-theoretic observables and their penalties."""
    svn = von_neumann_graph_entropy(poset)
    gini = degree_gini(poset)
    layer_ent = layer_size_entropy_normalised(poset)
    cover_ent = cover_binary_entropy_normalised(poset)
    interval_ent = interval_size_entropy_normalised(poset)

    p_spectral = spectral_entropy_deficit_penalty(poset)
    p_degree = degree_heterogeneity_penalty(poset)
    p_layer = layer_concentration_penalty(poset)
    p_edge = edge_density_extremity_penalty(poset)
    p_interval = interval_diversity_deficit_penalty(poset)

    components: dict[str, float] = {
        # Raw measures (diagnostics; not included in weighted total)
        "info_von_neumann_entropy": svn,
        "info_degree_gini": gini,
        "info_layer_entropy_norm": layer_ent,
        "info_cover_binary_entropy_norm": cover_ent,
        "info_interval_size_entropy_norm": interval_ent,
        # Penalty terms (these contribute to info_total)
        "info_spectral_entropy_deficit": p_spectral,
        "info_degree_heterogeneity": p_degree,
        "info_layer_concentration": p_layer,
        "info_edge_density_extremity": p_edge,
        "info_interval_diversity_deficit": p_interval,
    }
    components["info_total"] = weighted_info_total(components)
    return components


def info_penalty(poset: Poset, weights: Mapping[str, float] | None = None) -> float:
    """Aggregate info-theoretic penalty.  Drop-in replacement for geometric_penalty."""
    components = info_components(poset)
    return weighted_info_total(components, weights=weights)


def info_penalty_from_components(
    components: Mapping[str, float],
    weights: Mapping[str, float] | None = None,
) -> float:
    return weighted_info_total(components, weights=weights)
