"""N=64 Parallel Tempering Metropolis for Lorentzian window stabilization.

The N=64 bottleneck (§5.7.4): single-swap Metropolis gets trapped in local
minima with comp_frac ~0.57-0.64, above the Lorentzian window (0.30-0.55).
Cluster moves improved to 40% hit rate but are insufficient.

Parallel Tempering (Replica Exchange):
  - Run K replicas at temperatures T_0 < T_1 < ... < T_{K-1}
  - Each replica does standard Metropolis (swap + cluster moves)
  - Periodically attempt to swap configurations between adjacent replicas
  - Swap accepted with probability min(1, exp((β_i - β_j)(E_i - E_j)))
  - High-T replicas explore freely, low-T replicas exploit good basins
  - Configurations diffuse through temperature space, escaping local traps

Additionally implements simulated annealing as a simpler baseline.

Usage:
    python _n64_parallel_tempering.py
"""
from __future__ import annotations

import math
import sys
import time
from dataclasses import dataclass, field

import numpy as np

sys.path.insert(0, ".")

from generators import Poset, generate_lorentzian_like_2d, transitive_closure
from unified_functional import FunctionalWeights
from unified_functional_metropolis import (
    compute_F_fast,
    compute_observables_snapshot,
    count_comparable_pairs,
    propose_swap_move,
    hasse_diagram,
    CALIBRATED_WEIGHTS,
)
from _cluster_move import (
    propose_cluster_move,
    propose_multi_swap,
)


# ── Lorentzian window ──────────────────────────────────────────────────
WIN_LO, WIN_HI = 0.30, 0.55


# ── F7 computation (with sigmoid wall) ─────────────────────────────────

def compute_F7(poset: Poset, weights: FunctionalWeights) -> float:
    """Compute F7 with sigmoid wall (the definitive main model).

    F7 = logH + γΠ_geo - λΣ_hist + ηΞ_d + α(N)·σ((R - Rc)/w)

    Uses calibrated parameters from §5.10.7:
      α₀=16, q=-0.5, Rc=0.25, w=0.015, N₀=20
    """
    from unified_functional import (
        compute_log_H, compute_pi_geo, compute_sigma_hist, compute_xi_dim,
    )
    from observables import comparable_fraction

    log_H = compute_log_H(poset, n_runs=32)  # reduced for speed at N=64
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim, _ = compute_xi_dim(poset)

    # Sigmoid wall parameters (§5.10.7 main model)
    alpha_0 = 16.0
    q = -0.5
    N0 = 20.0
    Rc = 0.25
    w_sig = 0.015
    N = poset.n

    # Interval richness proxy: use comparable_fraction as R proxy
    # (R ≈ fraction of pairs that are comparable minus link fraction)
    comp = comparable_fraction(poset)
    # For the sigmoid wall, use comp_frac as R (simplification consistent
    # with the prediction experiments which use this proxy)
    R = comp

    alpha_N = alpha_0 * (N0 / max(N, 1)) ** abs(q)
    sigmoid_val = 1.0 / (1.0 + math.exp(-(R - Rc) / w_sig))
    wall = alpha_N * sigmoid_val

    F_base = (
        weights.beta * log_H
        + weights.gamma * pi_geo
        - weights.lam * sigma_hist
        + weights.eta * xi_dim
    )
    return F_base + wall


# ── Parallel Tempering ─────────────────────────────────────────────────

@dataclass
class PTReplica:
    """State of a single replica in parallel tempering."""
    poset: Poset
    F: float
    temperature: float
    accepted: int = 0
    proposed: int = 0
    swaps_accepted: int = 0
    swaps_proposed: int = 0


@dataclass
class PTResult:
    """Result of a parallel tempering run."""
    trajectories: list[list[dict]]  # per-replica trajectory
    final_posets: list[Poset]
    final_Fs: list[float]
    temperatures: list[float]
    swap_acceptance_rates: list[float]
    move_acceptance_rates: list[float]
    best_replica_idx: int
    elapsed_seconds: float


def parallel_tempering(
    initial_poset: Poset,
    weights: FunctionalWeights,
    n_steps: int = 3000,
    temperatures: list[float] | None = None,
    cluster_fraction: float = 0.4,
    swap_interval: int = 10,
    seed: int = 42,
    log_every: int = 50,
    use_F7: bool = False,
    label: str = "",
) -> PTResult:
    """Run parallel tempering with replica exchange.

    Parameters
    ----------
    initial_poset : starting configuration (cloned to all replicas)
    weights : functional weights
    n_steps : total Metropolis steps per replica
    temperatures : list of temperatures (default: geometric ladder)
    cluster_fraction : fraction of steps using cluster moves
    swap_interval : attempt replica swaps every this many steps
    seed : random seed
    log_every : log observables every this many steps
    use_F7 : if True, use F7 with sigmoid wall; else use F_fast (old F5)
    label : run label for logging
    """
    rng = np.random.default_rng(seed)

    if temperatures is None:
        # Geometric temperature ladder: 5 replicas from T=1.5 to T=12
        temperatures = [1.5, 2.5, 4.0, 7.0, 12.0]

    K = len(temperatures)
    compute_F = compute_F7 if use_F7 else lambda p, w: compute_F_fast(p, w)

    # Initialize replicas (all start from same configuration)
    replicas = []
    for k in range(K):
        p = Poset(initial_poset.closure.copy())
        F = compute_F(p, weights)
        replicas.append(PTReplica(poset=p, F=F, temperature=temperatures[k]))

    trajectories = [[] for _ in range(K)]

    def log_snapshot(k: int, step: int):
        r = replicas[k]
        obs = compute_observables_snapshot(r.poset)
        n_comp = count_comparable_pairs(r.poset.closure)
        trajectories[k].append({
            "step": step, "F": r.F, "T": r.temperature,
            "accepted": r.accepted, "n_comparable": n_comp,
            "replica": k, "label": label, **obs,
        })

    # Log initial state
    for k in range(K):
        log_snapshot(k, 0)

    swap_accept_counts = [0] * (K - 1)
    swap_attempt_counts = [0] * (K - 1)

    t0 = time.time()

    for step in range(1, n_steps + 1):
        # ── Step 1: Metropolis move for each replica ──
        for k in range(K):
            r = replicas[k]
            r.proposed += 1

            # Choose move type
            if rng.random() < cluster_fraction:
                new_closure = propose_cluster_move(r.poset.closure, rng)
            else:
                new_closure = propose_swap_move(r.poset.closure, rng)

            if new_closure is None:
                continue

            new_poset = Poset(new_closure)
            new_F = compute_F(new_poset, weights)

            delta_F = new_F - r.F
            beta_k = 1.0 / max(r.temperature, 1e-10)
            if delta_F < 0 or rng.random() < math.exp(-beta_k * delta_F):
                r.poset = new_poset
                r.F = new_F
                r.accepted += 1

        # ── Step 2: Replica exchange (every swap_interval steps) ──
        if step % swap_interval == 0:
            # Attempt swaps between adjacent replicas (random direction)
            if rng.random() < 0.5:
                pairs = range(0, K - 1)
            else:
                pairs = range(K - 2, -1, -1)

            for k in pairs:
                swap_attempt_counts[k] += 1
                r_lo = replicas[k]
                r_hi = replicas[k + 1]

                beta_lo = 1.0 / max(r_lo.temperature, 1e-10)
                beta_hi = 1.0 / max(r_hi.temperature, 1e-10)
                delta_beta = beta_lo - beta_hi
                delta_E = r_lo.F - r_hi.F

                # Metropolis criterion for swap
                log_accept = delta_beta * delta_E
                if log_accept > 0 or rng.random() < math.exp(log_accept):
                    # Swap configurations (not temperatures!)
                    r_lo.poset, r_hi.poset = r_hi.poset, r_lo.poset
                    r_lo.F, r_hi.F = r_hi.F, r_lo.F
                    swap_accept_counts[k] += 1

        # ── Step 3: Log ──
        if step % log_every == 0:
            for k in range(K):
                log_snapshot(k, step)

    # Final snapshot
    for k in range(K):
        log_snapshot(k, n_steps)

    elapsed = time.time() - t0

    # Find best replica (lowest T that ended in window)
    best_idx = 0
    best_F = float("inf")
    for k in range(K):
        obs = compute_observables_snapshot(replicas[k].poset)
        if WIN_LO <= obs["comp_frac"] <= WIN_HI:
            if replicas[k].F < best_F:
                best_F = replicas[k].F
                best_idx = k

    swap_rates = [
        swap_accept_counts[k] / max(swap_attempt_counts[k], 1)
        for k in range(K - 1)
    ]
    move_rates = [
        replicas[k].accepted / max(replicas[k].proposed, 1)
        for k in range(K)
    ]

    return PTResult(
        trajectories=trajectories,
        final_posets=[r.poset for r in replicas],
        final_Fs=[r.F for r in replicas],
        temperatures=temperatures,
        swap_acceptance_rates=swap_rates,
        move_acceptance_rates=move_rates,
        best_replica_idx=best_idx,
        elapsed_seconds=elapsed,
    )


# ── Simulated Annealing baseline ───────────────────────────────────────

@dataclass
class AnnealResult:
    trajectory: list[dict]
    final_poset: Poset
    final_F: float
    acceptance_rate: float
    elapsed_seconds: float


def simulated_annealing(
    initial_poset: Poset,
    weights: FunctionalWeights,
    n_steps: int = 5000,
    T_start: float = 10.0,
    T_end: float = 1.0,
    schedule: str = "geometric",
    cluster_fraction: float = 0.3,
    seed: int = 42,
    log_every: int = 50,
    use_F7: bool = False,
    label: str = "",
) -> AnnealResult:
    """Simulated annealing: cool from T_start to T_end over n_steps."""
    rng = np.random.default_rng(seed)
    compute_F = compute_F7 if use_F7 else lambda p, w: compute_F_fast(p, w)

    current = Poset(initial_poset.closure.copy())
    current_F = compute_F(current, weights)

    trajectory = []
    accepted = 0
    t0 = time.time()

    for step in range(n_steps):
        # Temperature schedule
        frac = step / max(n_steps - 1, 1)
        if schedule == "geometric":
            T = T_start * (T_end / T_start) ** frac
        elif schedule == "linear":
            T = T_start + (T_end - T_start) * frac
        else:
            T = T_start * (T_end / T_start) ** frac

        # Propose move
        if rng.random() < cluster_fraction:
            new_closure = propose_cluster_move(current.closure, rng)
        else:
            new_closure = propose_swap_move(current.closure, rng)

        if new_closure is not None:
            new_poset = Poset(new_closure)
            new_F = compute_F(new_poset, weights)
            delta_F = new_F - current_F
            beta = 1.0 / max(T, 1e-10)
            if delta_F < 0 or rng.random() < math.exp(-beta * delta_F):
                current = new_poset
                current_F = new_F
                accepted += 1

        if step % log_every == 0:
            obs = compute_observables_snapshot(current)
            trajectory.append({
                "step": step, "F": current_F, "T": T,
                "accepted": accepted, "label": label, **obs,
            })

    # Final
    obs = compute_observables_snapshot(current)
    trajectory.append({
        "step": n_steps, "F": current_F, "T": T_end,
        "accepted": accepted, "label": label, **obs,
    })

    elapsed = time.time() - t0
    return AnnealResult(
        trajectory=trajectory,
        final_poset=current,
        final_F=current_F,
        acceptance_rate=accepted / max(n_steps, 1),
        elapsed_seconds=elapsed,
    )


# ── Main experiment ────────────────────────────────────────────────────

def main():
    N = 64
    N_REPS = 5
    SEED_BASE = 77777

    print("=" * 80)
    print(f"N={N} PARALLEL TEMPERING & SIMULATED ANNEALING EXPERIMENT")
    print("=" * 80)
    print(f"  Reps: {N_REPS}")
    print(f"  Lorentzian window: comp_frac ∈ [{WIN_LO}, {WIN_HI}]")
    print()

    configs = [
        # (name, method, kwargs)
        ("SA-5k", "anneal", dict(
            n_steps=5000, T_start=10.0, T_end=1.5,
            cluster_fraction=0.3,
        )),
        ("PT-3k-5T", "pt", dict(
            n_steps=3000, temperatures=[1.5, 3.0, 5.0, 8.0, 14.0],
            cluster_fraction=0.4, swap_interval=10,
        )),
        ("PT-5k-7T", "pt", dict(
            n_steps=5000, temperatures=[1.5, 2.5, 4.0, 6.0, 9.0, 13.0, 20.0],
            cluster_fraction=0.4, swap_interval=8,
        )),
    ]

    all_results = {}

    for config_name, method, kwargs in configs:
        print(f"\n{'─' * 70}")
        print(f"Config: {config_name} ({method})")
        print(f"{'─' * 70}")

        hits = 0
        comp_fracs = []
        f_finals = []
        times = []

        for rep in range(N_REPS):
            seed = SEED_BASE + rep * 1000
            init = generate_lorentzian_like_2d(N, seed=seed)

            if method == "anneal":
                res = simulated_annealing(
                    init, CALIBRATED_WEIGHTS,
                    seed=seed, label=f"{config_name}_r{rep}",
                    **kwargs,
                )
                final_obs = res.trajectory[-1]
                elapsed = res.elapsed_seconds

            elif method == "pt":
                res = parallel_tempering(
                    init, CALIBRATED_WEIGHTS,
                    seed=seed, label=f"{config_name}_r{rep}",
                    **kwargs,
                )
                # Report the lowest-temperature replica
                final_obs = res.trajectories[0][-1]
                elapsed = res.elapsed_seconds

            comp = final_obs["comp_frac"]
            in_win = WIN_LO <= comp <= WIN_HI
            if in_win:
                hits += 1

            comp_fracs.append(comp)
            f_finals.append(final_obs["F"])
            times.append(elapsed)

            marker = "✓" if in_win else "✗"
            print(f"  rep {rep}: comp={comp:.3f} {marker}  "
                  f"F={final_obs['F']:.1f}  "
                  f"d_eff={final_obs['d_eff']:.2f}  "
                  f"layers={final_obs['n_layers']}  "
                  f"({elapsed:.1f}s)")

            # For PT, also report swap acceptance rates
            if method == "pt":
                sr = res.swap_acceptance_rates
                mr = res.move_acceptance_rates
                print(f"         swap_rates: {[f'{r:.2f}' for r in sr]}")
                print(f"         move_rates: {[f'{r:.2f}' for r in mr]}")

        hit_rate = hits / N_REPS
        all_results[config_name] = {
            "hit_rate": hit_rate,
            "comp_mean": np.mean(comp_fracs),
            "comp_std": np.std(comp_fracs),
            "F_mean": np.mean(f_finals),
            "time_mean": np.mean(times),
        }

        print(f"\n  → {config_name}: hit_rate={hit_rate:.0%}, "
              f"comp={np.mean(comp_fracs):.3f}±{np.std(comp_fracs):.3f}, "
              f"F={np.mean(f_finals):.1f}, "
              f"time={np.mean(times):.1f}s")

    # ── Summary ──
    print(f"\n{'=' * 80}")
    print("SUMMARY")
    print(f"{'=' * 80}")

    # Reference: N=64 cluster-only from §5.7.4
    print(f"\n  {'Config':<16s} {'Hit%':>6s} {'comp_μ':>8s} {'comp_σ':>8s} "
          f"{'F_μ':>8s} {'time_s':>8s}")
    print(f"  {'─'*16} {'─'*6} {'─'*8} {'─'*8} {'─'*8} {'─'*8}")
    print(f"  {'§5.7.4 clust50':<16s} {'40%':>6s} {'0.570':>8s} {'~0.04':>8s} "
          f"{'—':>8s} {'—':>8s}")

    for name, r in all_results.items():
        print(f"  {name:<16s} {r['hit_rate']:6.0%} "
              f"{r['comp_mean']:8.3f} {r['comp_std']:8.3f} "
              f"{r['F_mean']:8.1f} {r['time_mean']:8.1f}")

    best_name = max(all_results, key=lambda k: all_results[k]["hit_rate"])
    best = all_results[best_name]
    print(f"\n  Best: {best_name} — {best['hit_rate']:.0%} hit rate")

    if best["hit_rate"] >= 0.8:
        print("  ✓ VERDICT: N=64 Lorentzian window successfully stabilized!")
    elif best["hit_rate"] >= 0.5:
        print("  △ VERDICT: Substantial improvement over §5.7.4 baseline (40%)")
    else:
        print("  ✗ VERDICT: Still below 50% — need further optimization")


if __name__ == "__main__":
    main()
