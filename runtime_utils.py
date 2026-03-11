from __future__ import annotations

from entropy_exact import log_linear_extensions_exact
from entropy_sis import estimate_log_linear_extensions_sis


def estimate_entropy(poset, sis_runs: int, seed: int, exact_threshold: int) -> tuple[float, str]:
    """Shared entropy entrypoint for exact/SIS fallback."""
    if poset.n <= exact_threshold:
        return log_linear_extensions_exact(poset), "exact"
    mean, _ = estimate_log_linear_extensions_sis(poset, n_runs=sis_runs, seed=seed)
    return mean, "sis"
