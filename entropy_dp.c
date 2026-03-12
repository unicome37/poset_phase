#include <stdint.h>

/*
 * Skeleton for a future C-backed exact DP kernel.
 *
 * Important constraint:
 * The project requires exact linear-extension counts, not merely floating-point
 * log-counts. Since exact counts quickly exceed uint64_t, a real backend must
 * provide arbitrary-precision accumulation (or an equivalent exact interface).
 *
 * For that reason this file is currently only a compilation placeholder. It
 * establishes the expected backend artifact name (`entropy_dp.dll`) without
 * pretending that a uint64-only implementation would be correct.
 *
 * Once a bigint-preserving memoized DP is implemented here, it can be exposed
 * through entropy_exact_c.py without further changes to the project entrypoint.
 */

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

EXPORT int entropy_dp_backend_placeholder(void) {
    return 1;
}
