# Exact Backend Notes

## Current status

The project now has a stable plug-in surface for a future C-backed exact
linear-extension counter:

- [entropy_exact.py](/d:/Kiro/理论体系/poset_phase/entropy_exact.py) now treats the C backend as an optional fast path and cleanly falls back to Python exact DP.
- [entropy_exact_c.py](/d:/Kiro/理论体系/poset_phase/entropy_exact_c.py) defines the backend-loading interface.
- [entropy_dp.c](/d:/Kiro/理论体系/poset_phase/entropy_dp.c) is a placeholder compilation unit.

## Why there is no real C fast path yet

The bottleneck is no longer conceptual. It is representational.

The current Python implementation returns **exact arbitrary-precision counts**.
A naive `ctypes` C implementation using only `uint64_t` would be wrong for the
relevant regime, because the number of linear extensions quickly exceeds 64-bit
range.

So a real backend must provide both:

1. downset-DP acceleration outside CPython;
2. bigint-preserving accumulation and memoization.

Until that exists, the project keeps the Python exact path as the only exact
backend.

## Most realistic next implementation target

The smallest credible exact backend is:

- a C or Cython memoized downset DP;
- a compact bigint representation for memo values;
- the same public interface as `count_linear_extensions_exact()`.

At that point `entropy_exact_c.py` can expose a real fast path without changing
any higher-level experiment scripts.
