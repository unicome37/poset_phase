"""Cluster moves for microcanonical Metropolis on poset space.

Implements three types of large-scale moves to improve mixing time
for N ≥ 36 posets, where single-swap moves are too local.

Move types:
  1. layer_block_swap:  Swap positions of two node groups in the layer order
  2. layer_resplit:     Merge two adjacent layers, re-partition, rewire
  3. multi_swap:        Perform k independent single swaps as one move

All moves preserve the poset property (transitivity, acyclicity).
The microcanonical constraint (approximately fixed comparable pair count)
is enforced by the acceptance criterion, not by the move itself.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from generators import Poset, transitive_closure
from unified_functional_metropolis import (
    hasse_diagram,
    propose_swap_move,
    count_comparable_pairs,
    compute_F_fast,
    compute_observables_snapshot,
    CALIBRATED_WEIGHTS,
)
from unified_functional import FunctionalWeights


# ---------------------------------------------------------------------------
# Helper: compute layer assignment for each node
# ---------------------------------------------------------------------------

def compute_layer_assignment(closure: np.ndarray) -> list[list[int]]:
    """Assign each node to its layer (topological sort by longest incoming path)."""
    n = closure.shape[0]
    indeg = closure.sum(axis=0).astype(int)
    remaining = set(range(n))
    layers = []
    
    while remaining:
        mins = [i for i in remaining if indeg[i] == 0]
        if not mins:
            break  # cycle detected, shouldn't happen for valid poset
        layers.append(mins)
        for u in mins:
            remaining.remove(u)
            for v in np.where(closure[u])[0]:
                indeg[v] -= 1
    
    return layers


# ---------------------------------------------------------------------------
# Move 1: Layer block swap
# ---------------------------------------------------------------------------

def propose_layer_block_swap(
    closure: np.ndarray, rng: np.random.Generator
) -> np.ndarray | None:
    """Swap two adjacent layer blocks.
    
    Pick two adjacent layers L_i, L_{i+1}. For each pair (a in L_i, b in L_{i+1}):
    - If a < b in original, with probability p reverse to b < a (or make incomparable)
    - Recompute transitive closure
    
    This creates a "macroscopic" rearrangement of O(|L_i|·|L_{i+1}|) relations.
    """
    layers = compute_layer_assignment(closure)
    if len(layers) < 2:
        return None
    
    # Pick a random adjacent pair
    idx = rng.integers(0, len(layers) - 1)
    L_top = layers[idx]
    L_bot = layers[idx + 1]
    
    if len(L_top) == 0 or len(L_bot) == 0:
        return None
    
    # Extract Hasse diagram
    cover = hasse_diagram(closure)
    new_cover = cover.copy()
    
    # For each covering relation between L_top and L_bot, flip with prob 0.5
    flipped = 0
    for a in L_top:
        for b in L_bot:
            if cover[a, b]:
                if rng.random() < 0.5:
                    new_cover[a, b] = False
                    new_cover[b, a] = True
                    flipped += 1
    
    if flipped == 0:
        return None  # No change made
    
    # Recompute transitive closure
    new_closure = transitive_closure(new_cover)
    
    # Verify acyclicity
    if np.any(new_closure & new_closure.T):
        return None
    
    return new_closure


# ---------------------------------------------------------------------------
# Move 2: Layer merge-resplit
# ---------------------------------------------------------------------------

def propose_layer_resplit(
    closure: np.ndarray, rng: np.random.Generator
) -> np.ndarray | None:
    """Merge two adjacent layers and re-split into new layers.
    
    Pick adjacent layers L_i, L_{i+1}. Pool their nodes. Remove all internal
    relations among pooled nodes. Randomly re-partition the pool into two new
    sub-layers. Add random inter-layer relations. Recompute closure.
    """
    layers = compute_layer_assignment(closure)
    if len(layers) < 2:
        return None
    
    idx = rng.integers(0, len(layers) - 1)
    pool = layers[idx] + layers[idx + 1]
    pool_set = set(pool)
    n = closure.shape[0]
    
    if len(pool) < 2:
        return None
    
    # Start from Hasse diagram
    cover = hasse_diagram(closure)
    new_cover = cover.copy()
    
    # Remove all internal relations among pooled nodes
    for a in pool:
        for b in pool:
            new_cover[a, b] = False
    
    # Randomly split pool into two new sub-layers
    rng.shuffle(pool)
    split_point = rng.integers(1, len(pool))
    new_top = pool[:split_point]
    new_bot = pool[split_point:]
    
    # Add random covering relations from new_top to new_bot
    # Use a density similar to the original inter-layer density
    orig_density = sum(1 for a in layers[idx] for b in layers[idx+1] if cover[a, b])
    orig_total = max(len(layers[idx]) * len(layers[idx+1]), 1)
    target_density = orig_density / orig_total
    target_density = max(target_density, 0.2)  # at least 20% to avoid disconnection
    
    for a in new_top:
        for b in new_bot:
            if rng.random() < target_density:
                new_cover[a, b] = True
    
    # Also reconnect external edges: nodes above pool → new_top, new_bot → nodes below pool
    # (External edges are already preserved in new_cover since we only cleared internal ones)
    
    new_closure = transitive_closure(new_cover)
    
    if np.any(new_closure & new_closure.T):
        return None
    
    return new_closure


# ---------------------------------------------------------------------------
# Move 3: Multi-swap (k independent single swaps)
# ---------------------------------------------------------------------------

def propose_multi_swap(
    closure: np.ndarray, rng: np.random.Generator, k: int = 3
) -> np.ndarray | None:
    """Perform k sequential single swap moves as one compound move.
    
    This is the simplest cluster move: just chain multiple single swaps.
    The compound move is accepted or rejected as a whole.
    """
    current = closure.copy()
    any_changed = False
    
    for _ in range(k):
        result = propose_swap_move(current, rng)
        if result is not None:
            current = result
            any_changed = True
    
    if not any_changed:
        return None
    
    return current


# ---------------------------------------------------------------------------
# Cluster move dispatcher
# ---------------------------------------------------------------------------

def propose_cluster_move(
    closure: np.ndarray, 
    rng: np.random.Generator,
    move_weights: tuple[float, float, float] = (0.3, 0.3, 0.4),
) -> np.ndarray | None:
    """Propose a cluster move from the repertoire.
    
    Weights control the probability of each move type:
      [0] layer_block_swap
      [1] layer_resplit
      [2] multi_swap (k=3)
    """
    w = np.array(move_weights) / sum(move_weights)
    choice = rng.choice(3, p=w)
    
    if choice == 0:
        return propose_layer_block_swap(closure, rng)
    elif choice == 1:
        return propose_layer_resplit(closure, rng)
    else:
        return propose_multi_swap(closure, rng, k=3)


# ---------------------------------------------------------------------------
# Cluster Metropolis sampler
# ---------------------------------------------------------------------------

@dataclass
class ClusterMetropolisResult:
    trajectory: list[dict]
    final_poset: Poset
    acceptance_rate: float
    n_steps: int
    move_stats: dict


def cluster_metropolis_sample(
    initial_poset: Poset,
    weights: FunctionalWeights,
    n_steps: int = 500,
    temperature: float = 1.0,
    seed: int = 42,
    log_every: int = 10,
    label: str = "",
    cluster_fraction: float = 0.3,
    multi_k: int = 3,
) -> ClusterMetropolisResult:
    """Metropolis with mixed single-swap + cluster moves.
    
    At each step:
    - With probability cluster_fraction: propose a cluster move
    - Otherwise: propose a single swap move (original microcanonical)
    
    This combines the precision of single swaps with the long-range
    exploration of cluster moves.
    """
    rng = np.random.default_rng(seed)
    current = initial_poset
    current_F = compute_F_fast(current, weights)
    
    trajectory = []
    accepted = 0
    cluster_proposed = 0
    cluster_accepted = 0
    single_proposed = 0
    single_accepted = 0
    
    for step in range(n_steps):
        use_cluster = rng.random() < cluster_fraction
        
        if use_cluster:
            new_closure = propose_cluster_move(current.closure, rng)
            cluster_proposed += 1
        else:
            new_closure = propose_swap_move(current.closure, rng)
            single_proposed += 1
        
        if new_closure is None:
            if step % log_every == 0:
                obs = compute_observables_snapshot(current)
                n_comp = count_comparable_pairs(current.closure)
                trajectory.append({
                    "step": step, "F": current_F, "accepted": accepted,
                    "n_comparable": n_comp, "label": label, **obs,
                })
            continue
        
        new_poset = Poset(new_closure)
        new_F = compute_F_fast(new_poset, weights)
        
        delta_F = new_F - current_F
        if delta_F < 0 or rng.random() < math.exp(-delta_F / max(temperature, 1e-10)):
            current = new_poset
            current_F = new_F
            accepted += 1
            if use_cluster:
                cluster_accepted += 1
            else:
                single_accepted += 1
        
        if step % log_every == 0:
            obs = compute_observables_snapshot(current)
            n_comp = count_comparable_pairs(current.closure)
            trajectory.append({
                "step": step, "F": current_F, "accepted": accepted,
                "n_comparable": n_comp, "label": label, **obs,
            })
    
    # Final snapshot
    obs = compute_observables_snapshot(current)
    n_comp = count_comparable_pairs(current.closure)
    trajectory.append({
        "step": n_steps, "F": current_F, "accepted": accepted,
        "n_comparable": n_comp, "label": label, **obs,
    })
    
    move_stats = {
        "cluster_proposed": cluster_proposed,
        "cluster_accepted": cluster_accepted,
        "cluster_rate": cluster_accepted / max(cluster_proposed, 1),
        "single_proposed": single_proposed,
        "single_accepted": single_accepted,
        "single_rate": single_accepted / max(single_proposed, 1),
    }
    
    return ClusterMetropolisResult(
        trajectory=trajectory,
        final_poset=current,
        acceptance_rate=accepted / max(n_steps, 1),
        n_steps=n_steps,
        move_stats=move_stats,
    )


# ---------------------------------------------------------------------------
# Benchmark: compare single-swap vs cluster for various N
# ---------------------------------------------------------------------------

def run_benchmark(
    n_values: list[int] = [16, 20, 28, 36],
    steps: int = 300,
    T: float = 2.0,
    seed: int = 42,
) -> None:
    """Compare single-swap microcanonical vs cluster Metropolis."""
    from generators import generate_lorentzian_like_2d
    from unified_functional_metropolis import microcanonical_metropolis_sample
    
    print("="*75)
    print("CLUSTER MOVE BENCHMARK: Single-Swap vs Cluster Metropolis")
    print("="*75)
    
    print(f"\n  Settings: steps={steps}, T={T}")
    print(f"  Cluster: 30% cluster moves (layer_swap/resplit/multi_swap), 70% single swap")
    
    header = (f"  {'N':>3s} {'Method':>12s} {'accept':>8s} {'comp_f':>8s} "
              f"{'d_eff':>8s} {'layers':>8s} {'F_final':>8s} {'in_win':>8s} "
              f"{'F_drop':>8s}")
    
    print(f"\n{header}")
    print(f"  {'─'*3} {'─'*12} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*8}")
    
    results = []
    
    for n in n_values:
        for rep in range(3):
            s = seed + rep * 1000
            init = generate_lorentzian_like_2d(n, seed=s)
            init_F = compute_F_fast(init, CALIBRATED_WEIGHTS)
            
            # Single-swap (original microcanonical)
            res_single = microcanonical_metropolis_sample(
                init, CALIBRATED_WEIGHTS, n_steps=steps,
                temperature=T, seed=s, label=f"single_N{n}_r{rep}",
            )
            final_s = res_single.trajectory[-1]
            in_win_s = 0.30 <= final_s['comp_frac'] <= 0.55
            
            # Cluster
            res_cluster = cluster_metropolis_sample(
                init, CALIBRATED_WEIGHTS, n_steps=steps,
                temperature=T, seed=s, label=f"cluster_N{n}_r{rep}",
                cluster_fraction=0.3, multi_k=3,
            )
            final_c = res_cluster.trajectory[-1]
            in_win_c = 0.30 <= final_c['comp_frac'] <= 0.55
            
            if rep == 0:  # only print first rep for readability
                print(f"  {n:3d} {'single':>12s} {res_single.acceptance_rate:8.3f} "
                      f"{final_s['comp_frac']:8.3f} {final_s['d_eff']:8.2f} "
                      f"{final_s['n_layers']:8d} {final_s['F']:8.2f} "
                      f"{'✓' if in_win_s else '✗':>8s} {final_s['F']-init_F:+8.2f}")
                print(f"  {'':>3s} {'cluster':>12s} {res_cluster.acceptance_rate:8.3f} "
                      f"{final_c['comp_frac']:8.3f} {final_c['d_eff']:8.2f} "
                      f"{final_c['n_layers']:8d} {final_c['F']:8.2f} "
                      f"{'✓' if in_win_c else '✗':>8s} {final_c['F']-init_F:+8.2f}")
                ms = res_cluster.move_stats
                print(f"  {'':>3s} {'':>12s} cluster_accept={ms['cluster_rate']:.3f}, "
                      f"single_accept={ms['single_rate']:.3f}")
            
            results.append({
                "N": n, "rep": rep, "method": "single",
                "accept": res_single.acceptance_rate,
                "comp_frac": final_s['comp_frac'],
                "F_final": final_s['F'],
                "F_drop": final_s['F'] - init_F,
                "in_window": in_win_s,
                "n_layers": final_s['n_layers'],
            })
            results.append({
                "N": n, "rep": rep, "method": "cluster",
                "accept": res_cluster.acceptance_rate,
                "comp_frac": final_c['comp_frac'],
                "F_final": final_c['F'],
                "F_drop": final_c['F'] - init_F,
                "in_window": in_win_c,
                "n_layers": final_c['n_layers'],
            })
    
    # Summary
    print(f"\n{'='*75}")
    print("SUMMARY")
    print(f"{'='*75}")
    
    for n in n_values:
        for method in ["single", "cluster"]:
            rows = [r for r in results if r["N"] == n and r["method"] == method]
            mean_accept = np.mean([r["accept"] for r in rows])
            mean_F_drop = np.mean([r["F_drop"] for r in rows])
            win_rate = np.mean([r["in_window"] for r in rows])
            mean_comp = np.mean([r["comp_frac"] for r in rows])
            
            print(f"  N={n:2d} {method:>8s}: accept={mean_accept:.3f}, "
                  f"ΔF={mean_F_drop:+.2f}, window_hit={win_rate:.0%}, "
                  f"comp={mean_comp:.3f}")


if __name__ == "__main__":
    run_benchmark(n_values=[16, 20, 28, 36], steps=300, T=2.0)
