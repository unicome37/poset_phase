"""N=64 large-scale cluster Metropolis verification.

Tests whether cluster moves can stabilize the Lorentzian window at N=64,
where single-swap Metropolis previously failed due to mixing time bottleneck.

Configurations tested:
  1. single-300:   single swap only, 300 steps (baseline, expected to fail)
  2. single-1500:  single swap only, 1500 steps (more steps, still local)
  3. clust30-1500: 30% cluster, 1500 steps
  4. clust50-1500: 50% cluster, 1500 steps (recommended from §5.7.3)
  5. clust50-3000: 50% cluster, 3000 steps (extended)

Each config runs 10 repetitions. Reports:
  - Lorentzian window hit rate (comp_frac ∈ [0.30, 0.55])
  - Mean/std of comp_frac at final step
  - Acceptance rates (single vs cluster)
  - F drop from initial
"""
from __future__ import annotations

import sys
import time
from dataclasses import dataclass

import numpy as np

# Ensure parent dir is in path
sys.path.insert(0, ".")

from generators import Poset, generate_lorentzian_like_2d
from unified_functional import FunctionalWeights
from unified_functional_metropolis import (
    compute_F_fast,
    compute_observables_snapshot,
    count_comparable_pairs,
    microcanonical_metropolis_sample,
    CALIBRATED_WEIGHTS,
)
from _cluster_move import cluster_metropolis_sample


# Lorentzian window
WIN_LO, WIN_HI = 0.30, 0.55

N = 64
N_REPS = 5
T = 2.0
SEED_BASE = 12345


@dataclass
class ConfigResult:
    name: str
    comp_fracs: list[float]
    f_finals: list[float]
    f_drops: list[float]
    accept_rates: list[float]
    cluster_rates: list[float]
    single_rates: list[float]
    n_layers_list: list[int]
    d_effs: list[float]
    elapsed_secs: list[float]

    @property
    def hit_rate(self) -> float:
        return np.mean([WIN_LO <= c <= WIN_HI for c in self.comp_fracs])

    @property
    def comp_mean(self) -> float:
        return float(np.mean(self.comp_fracs))

    @property
    def comp_std(self) -> float:
        return float(np.std(self.comp_fracs))


def run_single_config(
    name: str,
    n_steps: int,
    cluster_fraction: float,
) -> ConfigResult:
    """Run one configuration across N_REPS repetitions."""
    result = ConfigResult(
        name=name,
        comp_fracs=[], f_finals=[], f_drops=[], accept_rates=[],
        cluster_rates=[], single_rates=[], n_layers_list=[], d_effs=[],
        elapsed_secs=[],
    )

    for rep in range(N_REPS):
        seed = SEED_BASE + rep * 1000
        init = generate_lorentzian_like_2d(N, seed=seed)
        init_F = compute_F_fast(init, CALIBRATED_WEIGHTS)

        t0 = time.time()

        if cluster_fraction <= 0:
            # Pure single-swap
            res = microcanonical_metropolis_sample(
                init, CALIBRATED_WEIGHTS, n_steps=n_steps,
                temperature=T, seed=seed,
                label=f"{name}_r{rep}",
            )
            cluster_rate = 0.0
            single_rate = res.acceptance_rate
        else:
            res = cluster_metropolis_sample(
                init, CALIBRATED_WEIGHTS, n_steps=n_steps,
                temperature=T, seed=seed,
                label=f"{name}_r{rep}",
                cluster_fraction=cluster_fraction,
                multi_k=3,
            )
            cluster_rate = res.move_stats["cluster_rate"]
            single_rate = res.move_stats["single_rate"]

        elapsed = time.time() - t0
        final = res.trajectory[-1]

        result.comp_fracs.append(final["comp_frac"])
        result.f_finals.append(final["F"])
        result.f_drops.append(final["F"] - init_F)
        result.accept_rates.append(res.acceptance_rate)
        result.cluster_rates.append(cluster_rate)
        result.single_rates.append(single_rate)
        result.n_layers_list.append(final["n_layers"])
        result.d_effs.append(final["d_eff"])
        result.elapsed_secs.append(elapsed)

        # Progress indicator
        in_win = WIN_LO <= final["comp_frac"] <= WIN_HI
        print(f"  [{name}] rep {rep+1}/{N_REPS}: "
              f"comp={final['comp_frac']:.3f} "
              f"{'✓' if in_win else '✗'} "
              f"F={final['F']:.1f} "
              f"ΔF={final['F']-init_F:+.1f} "
              f"({elapsed:.1f}s)")

    return result


def print_summary(results: list[ConfigResult]) -> None:
    """Print formatted summary table."""
    print("\n" + "=" * 90)
    print(f"N=64 CLUSTER METROPOLIS VERIFICATION — {N_REPS} reps each, T={T}")
    print("=" * 90)

    header = (f"  {'Config':<18s} {'Hit%':>6s} {'comp_μ':>8s} {'comp_σ':>8s} "
              f"{'ΔF_μ':>8s} {'acc%':>6s} {'clust%':>7s} {'d_eff':>7s} "
              f"{'layers':>7s} {'sec':>6s}")
    print(header)
    print("  " + "─" * 86)

    for r in results:
        print(f"  {r.name:<18s} "
              f"{r.hit_rate:6.0%} "
              f"{r.comp_mean:8.3f} "
              f"{r.comp_std:8.3f} "
              f"{np.mean(r.f_drops):+8.1f} "
              f"{np.mean(r.accept_rates):6.1%} "
              f"{np.mean(r.cluster_rates):7.1%} "
              f"{np.mean(r.d_effs):7.2f} "
              f"{np.mean(r.n_layers_list):7.1f} "
              f"{np.mean(r.elapsed_secs):6.1f}")

    # N=36 reference from §5.7
    print("\n  (§5.7 reference: N=36 clust50-600 → 100% hit, comp=0.524±0.012)")

    # Verdict
    best = max(results, key=lambda r: r.hit_rate)
    print(f"\n  Best config: {best.name} — hit rate {best.hit_rate:.0%}, "
          f"comp = {best.comp_mean:.3f} ± {best.comp_std:.3f}")

    if best.hit_rate >= 0.8:
        print("  ✓ VERDICT: Cluster moves successfully stabilize Lorentzian window at N=64")
    elif best.hit_rate >= 0.5:
        print("  △ VERDICT: Partial success — cluster moves help but more steps/tuning needed")
    else:
        print("  ✗ VERDICT: N=64 remains challenging — may need new move types or annealing")


def main():
    print(f"N=64 Cluster Metropolis Verification")
    print(f"  N={N}, T={T}, reps={N_REPS}, seed_base={SEED_BASE}")
    print(f"  Lorentzian window: comp_frac ∈ [{WIN_LO}, {WIN_HI}]")
    print()

    configs = [
        ("single-300",    300,  0.0),
        ("clust50-600",   600,  0.5),
        ("clust50-1500", 1500,  0.5),
    ]

    all_results = []
    for name, steps, cf in configs:
        print(f"\n--- {name} (steps={steps}, cluster_frac={cf}) ---")
        result = run_single_config(name, steps, cf)
        all_results.append(result)

    print_summary(all_results)


if __name__ == "__main__":
    main()
