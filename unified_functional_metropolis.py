"""Unified Structural Functional — Metropolis Dynamics

Implements a Metropolis-Hastings sampler on the space of finite posets,
using the calibrated unified functional F[X] as the energy function.

Local moves:
  - Add a covering relation (edge in Hasse diagram)
  - Remove a covering relation
  - Swap: remove one edge, add another

The sampler evolves a poset toward low-F structures. Four experiments
test whether the equilibrium distribution reproduces the four inferences:

  Exp-B: Starting from random poset, does the system evolve toward
         Lorentzian-like observables (not KR-like)?
  Exp-C: Does layer depth increase during evolution (sedimentation)?
  Exp-A: Does effective dimension converge near d=4 from various starts?
  Exp-D: Are equilibrium structures more CG-stable than initial ones?

Usage:
    python unified_functional_metropolis.py [--n N] [--steps STEPS] [--T TEMP]
"""
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from generators import Poset, transitive_closure
from unified_functional import (
    FunctionalWeights,
    FunctionalComponents,
    compute_log_H,
    compute_pi_geo,
    compute_sigma_hist,
    compute_xi_dim,
    compute_pi_cg,
)
from observables import comparable_fraction, layer_profile, antichain_width
from observables_geo import height_ratio, width_ratio, estimate_dimension_proxy_from_order_fraction


# ---------------------------------------------------------------------------
# Calibrated weights from grid search
# ---------------------------------------------------------------------------

CALIBRATED_WEIGHTS = FunctionalWeights(
    beta=2.0, gamma=0.5, lam=1.5, eta=0.1, kappa=0.05
)


# ---------------------------------------------------------------------------
# Hasse diagram (covering relations) extraction
# ---------------------------------------------------------------------------

def hasse_diagram(closure: np.ndarray) -> np.ndarray:
    """Extract the Hasse diagram (transitive reduction) from a closure.
    
    cover[i,j] = True iff i < j and there is no k with i < k < j.
    """
    n = closure.shape[0]
    c = closure.astype(np.uint8)
    has_intermediate = (c @ c).astype(bool)
    cover = closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return cover


# ---------------------------------------------------------------------------
# Local move operators
# ---------------------------------------------------------------------------

def propose_add_edge(closure: np.ndarray, rng: np.random.Generator) -> np.ndarray | None:
    """Propose adding a new covering relation between two incomparable elements."""
    n = closure.shape[0]
    # Find incomparable pairs (neither i<j nor j<i)
    comparable = closure | closure.T
    np.fill_diagonal(comparable, True)
    incomparable = ~comparable
    candidates = np.argwhere(np.triu(incomparable, k=1))
    
    if len(candidates) == 0:
        return None
    
    idx = rng.integers(len(candidates))
    i, j = candidates[idx]
    
    # Add edge i → j (or j → i, chosen randomly)
    new_closure = closure.copy()
    if rng.random() < 0.5:
        new_closure[i, j] = True
    else:
        new_closure[j, i] = True
    
    # Recompute transitive closure
    new_closure = transitive_closure(new_closure)
    
    # Verify it's still a valid partial order (no cycles)
    if np.any(new_closure & new_closure.T):
        return None  # Would create a cycle
    
    return new_closure


def propose_remove_edge(closure: np.ndarray, rng: np.random.Generator) -> np.ndarray | None:
    """Propose removing a covering relation."""
    cover = hasse_diagram(closure)
    edges = np.argwhere(cover)
    
    if len(edges) == 0:
        return None
    
    idx = rng.integers(len(edges))
    i, j = edges[idx]
    
    # Remove this covering relation and recompute closure from remaining covers
    new_cover = cover.copy()
    new_cover[i, j] = False
    new_closure = transitive_closure(new_cover)
    
    return new_closure


def propose_move(closure: np.ndarray, rng: np.random.Generator, remove_bias: float = 0.5) -> np.ndarray | None:
    """Propose a local move: add or remove an edge.
    
    remove_bias controls the probability of proposing a removal vs addition.
    Default 0.5 = equal probability. Set higher to counteract the tendency
    of sparse posets to over-order (since there are many more incomparable
    pairs than covering edges at low density).
    """
    if rng.random() < remove_bias:
        return propose_remove_edge(closure, rng)
    else:
        return propose_add_edge(closure, rng)


# ---------------------------------------------------------------------------
# Microcanonical move: swap (remove one edge + add one) to fix relation count
# ---------------------------------------------------------------------------

def count_comparable_pairs(closure: np.ndarray) -> int:
    """Count the number of comparable pairs (i<j or j<i) in the poset."""
    return int(np.sum(closure | closure.T)) // 2


def propose_swap_move(closure: np.ndarray, rng: np.random.Generator) -> np.ndarray | None:
    """Propose a microcanonical swap: remove one Hasse edge, add one new relation.
    
    This keeps the total number of comparable pairs approximately fixed,
    preventing the over-ordering that occurs in canonical (add/remove) moves.
    The swap is: remove a covering relation, then add a relation between
    a randomly chosen incomparable pair.
    """
    cover = hasse_diagram(closure)
    edges = np.argwhere(cover)
    
    if len(edges) == 0:
        return None
    
    # Step 1: Remove a random covering relation
    rm_idx = rng.integers(len(edges))
    ri, rj = edges[rm_idx]
    new_cover = cover.copy()
    new_cover[ri, rj] = False
    intermediate_closure = transitive_closure(new_cover)
    
    # Step 2: Find incomparable pairs in the intermediate state
    comparable = intermediate_closure | intermediate_closure.T
    np.fill_diagonal(comparable, True)
    incomparable = ~comparable
    candidates = np.argwhere(np.triu(incomparable, k=1))
    
    if len(candidates) == 0:
        return None
    
    # Step 3: Add a random new relation
    add_idx = rng.integers(len(candidates))
    ai, aj = candidates[add_idx]
    
    if rng.random() < 0.5:
        intermediate_closure[ai, aj] = True
    else:
        intermediate_closure[aj, ai] = True
    
    new_closure = transitive_closure(intermediate_closure)
    
    # Verify no cycles
    if np.any(new_closure & new_closure.T):
        return None
    
    return new_closure


def microcanonical_metropolis_sample(
    initial_poset: Poset,
    weights: FunctionalWeights,
    n_steps: int = 500,
    temperature: float = 1.0,
    seed: int = 42,
    log_every: int = 10,
    label: str = "",
) -> "MetropolisResult":
    """Run microcanonical Metropolis: swap moves only (fixed relation count).
    
    At each step:
    1. Propose a swap move (remove one edge + add one relation)
    2. Compute ΔF = F[X'] - F[X]
    3. Accept with probability min(1, exp(-ΔF/T))
    
    The microcanonical constraint prevents over-ordering by keeping
    the comparable pair count approximately fixed.
    """
    rng = np.random.default_rng(seed)
    current = initial_poset
    current_F = compute_F_fast(current, weights)
    
    trajectory = []
    accepted = 0
    
    for step in range(n_steps):
        new_closure = propose_swap_move(current.closure, rng)
        
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
        
        if step % log_every == 0:
            obs = compute_observables_snapshot(current)
            n_comp = count_comparable_pairs(current.closure)
            trajectory.append({
                "step": step, "F": current_F, "accepted": accepted,
                "n_comparable": n_comp, "label": label, **obs,
            })
    
    obs = compute_observables_snapshot(current)
    n_comp = count_comparable_pairs(current.closure)
    trajectory.append({
        "step": n_steps, "F": current_F, "accepted": accepted,
        "n_comparable": n_comp, "label": label, **obs,
    })
    
    return MetropolisResult(
        trajectory=trajectory,
        final_poset=current,
        acceptance_rate=accepted / max(n_steps, 1),
        n_steps=n_steps,
    )


# ---------------------------------------------------------------------------
# Fast F[X] computation (skip expensive pi_cg for inner loop)
# ---------------------------------------------------------------------------

def compute_F_fast(poset: Poset, weights: FunctionalWeights) -> float:
    """Fast F[X] computation skipping the expensive pi_cg term.
    
    For Metropolis inner loop, we approximate F by dropping κ·Π_cg
    since κ=0.05 is the smallest weight and Π_cg requires multiple
    coarse-graining samples.
    """
    log_H = compute_log_H(poset, n_runs=64)  # fewer SIS runs for speed
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim, _ = compute_xi_dim(poset)
    
    return (
        weights.beta * log_H
        + weights.gamma * pi_geo
        - weights.lam * sigma_hist
        + weights.eta * xi_dim
    )


def compute_observables_snapshot(poset: Poset) -> dict:
    """Compute a snapshot of key observables for logging."""
    comp = comparable_fraction(poset)
    profile = layer_profile(poset)
    n_layers = len(profile)
    d_eff = estimate_dimension_proxy_from_order_fraction(comp)
    aw = antichain_width(poset)
    
    return {
        "comp_frac": comp,
        "n_layers": n_layers,
        "d_eff": d_eff,
        "aw": aw,
        "aw_ratio": aw / max(poset.n, 1),
        "height_ratio": height_ratio(poset),
        "width_ratio": width_ratio(poset),
    }


# ---------------------------------------------------------------------------
# Metropolis sampler
# ---------------------------------------------------------------------------

@dataclass
class MetropolisResult:
    """Result of a Metropolis run."""
    trajectory: list[dict]
    final_poset: Poset
    acceptance_rate: float
    n_steps: int


def metropolis_sample(
    initial_poset: Poset,
    weights: FunctionalWeights,
    n_steps: int = 500,
    temperature: float = 1.0,
    seed: int = 42,
    log_every: int = 10,
    label: str = "",
    remove_bias: float = 0.5,
) -> MetropolisResult:
    """Run Metropolis-Hastings sampling on poset space.
    
    At each step:
    1. Propose a local move (add/remove covering relation)
    2. Compute ΔF = F[X'] - F[X]
    3. Accept with probability min(1, exp(-ΔF/T))
    """
    rng = np.random.default_rng(seed)
    current = initial_poset
    current_F = compute_F_fast(current, weights)
    
    trajectory = []
    accepted = 0
    
    for step in range(n_steps):
        # Propose move
        new_closure = propose_move(current.closure, rng, remove_bias=remove_bias)
        
        if new_closure is None:
            # Invalid move, skip
            if step % log_every == 0:
                obs = compute_observables_snapshot(current)
                trajectory.append({
                    "step": step, "F": current_F, "accepted": accepted,
                    "label": label, **obs,
                })
            continue
        
        new_poset = Poset(new_closure)
        new_F = compute_F_fast(new_poset, weights)
        
        # Metropolis acceptance
        delta_F = new_F - current_F
        if delta_F < 0 or rng.random() < math.exp(-delta_F / max(temperature, 1e-10)):
            current = new_poset
            current_F = new_F
            accepted += 1
        
        # Log
        if step % log_every == 0:
            obs = compute_observables_snapshot(current)
            trajectory.append({
                "step": step, "F": current_F, "accepted": accepted,
                "label": label, **obs,
            })
    
    # Final snapshot
    obs = compute_observables_snapshot(current)
    trajectory.append({
        "step": n_steps, "F": current_F, "accepted": accepted,
        "label": label, **obs,
    })
    
    return MetropolisResult(
        trajectory=trajectory,
        final_poset=current,
        acceptance_rate=accepted / max(n_steps, 1),
        n_steps=n_steps,
    )


# ---------------------------------------------------------------------------
# Experiment B: Phase competition
# ---------------------------------------------------------------------------

def experiment_B(n: int, steps: int, T: float, seed: int) -> list[dict]:
    """Start from random transitive percolation (structureless),
    evolve under F, check if result resembles Lorentzian-like."""
    from generators import generate_transitive_percolation, generate_lorentzian_like_2d, generate_kr_like
    
    print(f"\n{'='*50}")
    print(f"Experiment B: Phase Competition (N={n}, T={T})")
    print(f"{'='*50}")
    
    all_traj = []
    
    for rep in range(3):
        s = seed + rep * 1000
        
        # Start from random percolation
        init = generate_transitive_percolation(n, p=0.12, seed=s)
        result = metropolis_sample(
            init, CALIBRATED_WEIGHTS, n_steps=steps,
            temperature=T, seed=s, label=f"random_rep{rep}",
        )
        all_traj.extend(result.trajectory)
        
        # Reference: Lorentzian-like observables
        lor_ref = generate_lorentzian_like_2d(n, seed=s)
        lor_obs = compute_observables_snapshot(lor_ref)
        lor_F = compute_F_fast(lor_ref, CALIBRATED_WEIGHTS)
        
        kr_ref = generate_kr_like(n, seed=s)
        kr_obs = compute_observables_snapshot(kr_ref)
        kr_F = compute_F_fast(kr_ref, CALIBRATED_WEIGHTS)
        
        final_obs = result.trajectory[-1]
        
        print(f"\n  Rep {rep}: acceptance={result.acceptance_rate:.2f}")
        print(f"    Initial  → F={result.trajectory[0]['F']:.2f}, comp={result.trajectory[0]['comp_frac']:.3f}, d_eff={result.trajectory[0]['d_eff']:.2f}")
        print(f"    Final    → F={final_obs['F']:.2f}, comp={final_obs['comp_frac']:.3f}, d_eff={final_obs['d_eff']:.2f}")
        print(f"    Lor2D ref→ F={lor_F:.2f}, comp={lor_obs['comp_frac']:.3f}, d_eff={lor_obs['d_eff']:.2f}")
        print(f"    KR ref   → F={kr_F:.2f}, comp={kr_obs['comp_frac']:.3f}, d_eff={kr_obs['d_eff']:.2f}")
        
        # Check: is final closer to Lor or KR?
        lor_dist = abs(final_obs['comp_frac'] - lor_obs['comp_frac'])
        kr_dist = abs(final_obs['comp_frac'] - kr_obs['comp_frac'])
        closer = "Lor2D" if lor_dist < kr_dist else "KR"
        print(f"    Final comp_frac closer to: {closer} (dist_Lor={lor_dist:.3f}, dist_KR={kr_dist:.3f})")
    
    return all_traj


# ---------------------------------------------------------------------------
# Experiment C: Sedimentation
# ---------------------------------------------------------------------------

def experiment_C(n: int, steps: int, T: float, seed: int) -> list[dict]:
    """Track layer depth growth during Metropolis evolution.
    Prediction C says: evolution under F should increase layer depth."""
    from generators import generate_transitive_percolation
    
    print(f"\n{'='*50}")
    print(f"Experiment C: Sedimentation (N={n}, T={T})")
    print(f"{'='*50}")
    
    all_traj = []
    
    for rep in range(3):
        s = seed + rep * 1000
        init = generate_transitive_percolation(n, p=0.10, seed=s)
        result = metropolis_sample(
            init, CALIBRATED_WEIGHTS, n_steps=steps,
            temperature=T, seed=s, label=f"sediment_rep{rep}",
        )
        all_traj.extend(result.trajectory)
        
        init_layers = result.trajectory[0]['n_layers']
        final_layers = result.trajectory[-1]['n_layers']
        init_F = result.trajectory[0]['F']
        final_F = result.trajectory[-1]['F']
        
        print(f"  Rep {rep}: layers {init_layers} → {final_layers} "
              f"(Δ={final_layers - init_layers:+d}), "
              f"F {init_F:.2f} → {final_F:.2f} "
              f"(Δ={final_F - init_F:+.2f}), "
              f"accept={result.acceptance_rate:.2f}")
    
    return all_traj


# ---------------------------------------------------------------------------
# Experiment A: Dimensional selection
# ---------------------------------------------------------------------------

def experiment_A(n: int, steps: int, T: float, seed: int) -> list[dict]:
    """Initialize Lorentzian-like at various dimensions, evolve under F.
    Check if d_eff converges toward a common basin."""
    from generators import (
        generate_lorentzian_like_2d,
        generate_lorentzian_like_3d,
        generate_lorentzian_like_4d,
        generate_lorentzian_like_5d,
    )
    
    print(f"\n{'='*50}")
    print(f"Experiment A: Dimensional Selection (N={n}, T={T})")
    print(f"{'='*50}")
    
    gens = {
        "Lor2D": generate_lorentzian_like_2d,
        "Lor3D": generate_lorentzian_like_3d,
        "Lor4D": generate_lorentzian_like_4d,
        "Lor5D": generate_lorentzian_like_5d,
    }
    
    all_traj = []
    
    for dim_name, gen in gens.items():
        s = seed + hash(dim_name) % 10000
        init = gen(n, seed=s)
        result = metropolis_sample(
            init, CALIBRATED_WEIGHTS, n_steps=steps,
            temperature=T, seed=s, label=f"dim_{dim_name}",
        )
        all_traj.extend(result.trajectory)
        
        init_d = result.trajectory[0]['d_eff']
        final_d = result.trajectory[-1]['d_eff']
        init_F = result.trajectory[0]['F']
        final_F = result.trajectory[-1]['F']
        
        print(f"  {dim_name}: d_eff {init_d:.2f} → {final_d:.2f} "
              f"(Δ={final_d - init_d:+.2f}), "
              f"F {init_F:.2f} → {final_F:.2f}, "
              f"accept={result.acceptance_rate:.2f}")
    
    return all_traj


# ---------------------------------------------------------------------------
# Experiment D: CG stability
# ---------------------------------------------------------------------------

def experiment_D(n: int, steps: int, T: float, seed: int) -> list[dict]:
    """Evolve structures, then compare CG drift before vs after evolution.
    Prediction D: equilibrium structures should be more CG-stable."""
    from generators import generate_transitive_percolation
    from coarse_grain import coarse_grain_delete_nodes
    
    print(f"\n{'='*50}")
    print(f"Experiment D: CG Stability (N={n}, T={T})")
    print(f"{'='*50}")
    
    results = []
    
    for rep in range(3):
        s = seed + rep * 1000
        init = generate_transitive_percolation(n, p=0.12, seed=s)
        
        # Measure initial CG drift
        init_pi_cg = compute_pi_cg(init, keep_ratio=0.7, n_cg_samples=5)
        
        # Evolve
        result = metropolis_sample(
            init, CALIBRATED_WEIGHTS, n_steps=steps,
            temperature=T, seed=s, label=f"cg_rep{rep}",
        )
        
        # Measure final CG drift
        final_pi_cg = compute_pi_cg(result.final_poset, keep_ratio=0.7, n_cg_samples=5)
        
        row = {
            "rep": rep,
            "pi_cg_init": init_pi_cg,
            "pi_cg_final": final_pi_cg,
            "delta_pi_cg": final_pi_cg - init_pi_cg,
            "F_init": result.trajectory[0]['F'],
            "F_final": result.trajectory[-1]['F'],
            "accept": result.acceptance_rate,
        }
        results.append(row)
        
        print(f"  Rep {rep}: Π_cg {init_pi_cg:.4f} → {final_pi_cg:.4f} "
              f"(Δ={final_pi_cg - init_pi_cg:+.4f}), "
              f"F {row['F_init']:.2f} → {row['F_final']:.2f}, "
              f"accept={result.acceptance_rate:.2f}")
    
    # Summary
    mean_init = np.mean([r['pi_cg_init'] for r in results])
    mean_final = np.mean([r['pi_cg_final'] for r in results])
    print(f"\n  Mean Π_cg: {mean_init:.4f} → {mean_final:.4f} "
          f"({'↓ MORE STABLE' if mean_final < mean_init else '↑ less stable'})")
    
    return results


# ---------------------------------------------------------------------------
# Experiment M: Microcanonical phase competition
# ---------------------------------------------------------------------------

def experiment_microcanonical_B(n: int, steps: int, T: float, seed: int, outdir: Path) -> None:
    """Microcanonical version of Exp-B: start from various initial densities,
    evolve under swap moves (fixed relation count), check if observables
    converge to Lorentzian-like window.
    
    We initialize Lorentzian-like posets at different dimensions (which have
    different comp_frac values) and evolve them under swap moves. Since the
    comparable pair count is approximately conserved, the system cannot
    over-order — it must redistribute relations to minimize F.
    """
    from generators import (
        generate_lorentzian_like_2d,
        generate_lorentzian_like_3d,
        generate_lorentzian_like_4d,
        generate_transitive_percolation,
        generate_kr_like,
    )
    
    print(f"\n{'='*60}")
    print(f"Experiment M: Microcanonical Phase Competition (N={n}, T={T})")
    print(f"{'='*60}")
    
    # Different initial conditions spanning the comp_frac range
    inits = {}
    for s_offset, (name, gen, kwargs) in enumerate([
        ("Lor2D", generate_lorentzian_like_2d, {}),
        ("Lor3D", generate_lorentzian_like_3d, {}),
        ("Lor4D", generate_lorentzian_like_4d, {}),
        ("KR", generate_kr_like, {}),
        ("TP_sparse", generate_transitive_percolation, {"p": 0.08}),
        ("TP_dense", generate_transitive_percolation, {"p": 0.20}),
    ]):
        s = seed + s_offset * 1000
        p = Poset(gen(n, seed=s).closure) if not kwargs else Poset(gen(n, seed=s, **kwargs).closure)
        inits[name] = (p, s)
    
    all_traj = []
    
    # Reference observables
    lor2d_ref = generate_lorentzian_like_2d(n, seed=seed)
    lor_obs = compute_observables_snapshot(lor2d_ref)
    lor_F = compute_F_fast(lor2d_ref, CALIBRATED_WEIGHTS)
    print(f"\n  Lor2D reference: comp={lor_obs['comp_frac']:.3f}, d_eff={lor_obs['d_eff']:.2f}, F={lor_F:.2f}")
    
    for name, (init_poset, s) in inits.items():
        init_obs = compute_observables_snapshot(init_poset)
        init_comp_pairs = count_comparable_pairs(init_poset.closure)
        
        result = microcanonical_metropolis_sample(
            init_poset, CALIBRATED_WEIGHTS, n_steps=steps,
            temperature=T, seed=s, label=f"mc_{name}",
        )
        all_traj.extend(result.trajectory)
        
        final = result.trajectory[-1]
        final_comp_pairs = final.get('n_comparable', '?')
        
        # Check proximity to Lor2D window
        lor_dist = abs(final['comp_frac'] - lor_obs['comp_frac'])
        in_window = 0.30 <= final['comp_frac'] <= 0.55
        
        print(f"\n  {name}:")
        print(f"    comp_pairs: {init_comp_pairs} → {final_comp_pairs}")
        print(f"    comp_frac:  {init_obs['comp_frac']:.3f} → {final['comp_frac']:.3f} "
              f"(Lor2D ref={lor_obs['comp_frac']:.3f}, dist={lor_dist:.3f})")
        print(f"    d_eff:      {init_obs['d_eff']:.2f} → {final['d_eff']:.2f}")
        print(f"    F:          {init_obs.get('F', '?')} → {final['F']:.2f}")
        print(f"    n_layers:   {init_obs['n_layers']} → {final['n_layers']}")
        print(f"    accept:     {result.acceptance_rate:.2f}")
        print(f"    IN LORENTZIAN WINDOW (0.30-0.55): {'✓ YES' if in_window else '✗ NO'}")
    
    save_csv(all_traj, str(outdir / "metropolis_M_microcanonical.csv"))
    
    # Temperature scan for microcanonical
    print(f"\n  --- Temperature scan (microcanonical, Lor2D start) ---")
    print(f"  {'T':>6s} {'comp':>8s} {'d_eff':>8s} {'F':>8s} {'layers':>8s} {'accept':>8s} {'comp_pairs':>10s}")
    
    for T_scan in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]:
        init = generate_lorentzian_like_2d(n, seed=seed)
        result = microcanonical_metropolis_sample(
            init, CALIBRATED_WEIGHTS, n_steps=steps,
            temperature=T_scan, seed=seed, log_every=steps,
        )
        final = result.trajectory[-1]
        print(f"  {T_scan:6.1f} {final['comp_frac']:8.3f} {final['d_eff']:8.2f} "
              f"{final['F']:8.2f} {final['n_layers']:8d} {result.acceptance_rate:8.2f} "
              f"{final.get('n_comparable', '?'):>10}")


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def save_csv(rows: list[dict], path: str) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  → saved {path} ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Unified Functional — Metropolis Dynamics")
    parser.add_argument("--n", type=int, default=16, help="Poset size")
    parser.add_argument("--steps", type=int, default=300, help="Metropolis steps per run")
    parser.add_argument("--T", type=float, default=2.0, help="Temperature")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--outdir", type=str, default="outputs_unified_functional",
                        help="Output directory")
    parser.add_argument("--experiments", type=str, default="BCAD",
                        help="Which experiments to run (e.g. 'B', 'BC', 'BCAD')")
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    print(f"Metropolis Dynamics — N={args.n}, steps={args.steps}, T={args.T}")
    print(f"Calibrated weights: β=2.0, γ=0.5, λ=1.5, η=0.1, κ=0.05")
    
    if "B" in args.experiments:
        traj_B = experiment_B(args.n, args.steps, args.T, args.seed)
        save_csv(traj_B, str(outdir / "metropolis_B_trajectory.csv"))
    
    if "C" in args.experiments:
        traj_C = experiment_C(args.n, args.steps, args.T, args.seed + 5000)
        save_csv(traj_C, str(outdir / "metropolis_C_trajectory.csv"))
    
    if "A" in args.experiments:
        traj_A = experiment_A(args.n, args.steps, args.T, args.seed + 10000)
        save_csv(traj_A, str(outdir / "metropolis_A_trajectory.csv"))
    
    if "D" in args.experiments:
        results_D = experiment_D(args.n, args.steps, args.T, args.seed + 15000)
        save_csv(results_D, str(outdir / "metropolis_D_cg_stability.csv"))
    
    # --- Microcanonical experiments ---
    if "M" in args.experiments:
        experiment_microcanonical_B(args.n, args.steps, args.T, args.seed, outdir)
    
    print(f"\n{'='*50}")
    print(f"Metropolis dynamics complete. Results in: {outdir}")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()
