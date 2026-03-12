from __future__ import annotations

from entropy_exact import log_linear_extensions_exact
from entropy_sis import estimate_log_linear_extensions_sis


def estimate_entropy(poset, sis_runs: int, seed: int, exact_threshold: int) -> tuple[float, str]:
    """Shared entropy entrypoint for exact/SIS fallback."""
    if poset.n <= exact_threshold:
        return log_linear_extensions_exact(poset), "exact"
    mean, _ = estimate_log_linear_extensions_sis(poset, n_runs=sis_runs, seed=seed)
    return mean, "sis"


def estimate_entropy_by_family(
    poset,
    family: str,
    sis_runs: int,
    seed: int,
    default_exact_threshold: int,
    family_exact_thresholds: dict[str, int] | None = None,
) -> tuple[float, str]:
    """Family-aware entropy entrypoint.

    The current exact wall is family-dependent: Lorentzian-like 2D remains cheap
    at larger N, while KR-like quickly becomes the bottleneck. This helper keeps
    one shared policy surface while allowing per-family exact thresholds.
    """
    thresholds = family_exact_thresholds or {}
    exact_threshold = int(thresholds.get(family, default_exact_threshold))
    return estimate_entropy(poset, sis_runs=sis_runs, seed=seed, exact_threshold=exact_threshold)
