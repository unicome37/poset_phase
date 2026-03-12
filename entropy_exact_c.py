from __future__ import annotations

"""Optional C-backed exact entropy interface.

This module is intentionally conservative. The current project can benefit from
moving the downset DP core out of CPython, but an actual exact fast path must
still preserve arbitrary-precision counts. Since the local environment may not
have a compiler available and the bigint-preserving C backend is not yet
implemented, this module currently exposes a stable integration surface:

- if a future compiled backend exists, it can be loaded here;
- if not, callers get a clean "backend unavailable" signal and can fall back
  to the Python exact implementation.

The goal is to freeze the interface now, so future C/Cython work does not
require another round of project-wide entrypoint changes.
"""

from ctypes import CDLL
from pathlib import Path

from generators import Poset


class ExactCBackendUnavailable(RuntimeError):
    pass


def _candidate_library_paths() -> list[Path]:
    here = Path(__file__).resolve().parent
    return [
        here / "entropy_dp.dll",
        here / "build" / "entropy_dp.dll",
    ]


def load_backend() -> CDLL:
    for path in _candidate_library_paths():
        if path.exists():
            return CDLL(str(path))
    raise ExactCBackendUnavailable(
        "No compiled entropy_dp backend found. Expected one of: "
        + ", ".join(str(p) for p in _candidate_library_paths())
    )


def count_linear_extensions_exact_c(poset: Poset) -> int:
    """Placeholder for a future bigint-preserving C implementation.

    We intentionally do not provide an unsafe uint64-only exact path here,
    because that would silently corrupt large counts. Until the compiled backend
    supports arbitrary-precision exact counts, this function remains a guarded
    stub.
    """
    _ = poset
    load_backend()
    raise ExactCBackendUnavailable(
        "Compiled backend detected, but exact bigint-preserving entrypoints are not yet wired."
    )


def c_backend_status() -> dict[str, str]:
    try:
        backend = load_backend()
        return {"available": "yes", "path": str(backend._name)}
    except ExactCBackendUnavailable as exc:
        return {"available": "no", "reason": str(exc)}
