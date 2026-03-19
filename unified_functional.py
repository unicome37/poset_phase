"""Unified Structural Functional — Toy Model

Computes the unified structural functional F[X] for finite posets,
connecting all four inferences (B/C/A/D) through a single cost functional.

F[X] = β·log H(X) + γ·Π_geo(X) − λ·Σ_hist(X) + η·Ξ_d(X) + κ·Π_cg(X)

Usage:
    python unified_functional.py [--n N] [--reps REPS] [--seed SEED]
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np

from generators import (
    Poset,
    generate_kr_like,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)
from observables import (
    comparable_fraction,
    layer_profile,
    neutral_penalty,
    antichain_width,
)
from observables_geo import (
    geometric_penalty,
    dimension_consistency_penalty,
    estimate_dimension_proxy_from_order_fraction,
    height_ratio,
    width_ratio,
)
from coarse_grain import coarse_grain_delete_nodes
from entropy_sis import estimate_log_linear_extensions_sis


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class FunctionalWeights:
    """Weights (β, γ, λ, η, κ) for the unified structural functional."""
    beta: float = 1.0    # combinatorial entropy weight
    gamma: float = 1.0   # geometric constraint weight
    lam: float = 0.5     # historical sedimentation benefit weight
    eta: float = 0.3     # dimensional barrier weight
    kappa: float = 0.2   # coarse-graining drift penalty weight


@dataclass
class FunctionalComponents:
    """Individual components of the unified functional."""
    log_H: float = 0.0
    pi_geo: float = 0.0
    sigma_hist: float = 0.0
    xi_dim: float = 0.0
    pi_cg: float = 0.0
    # Derived
    F_total: float = 0.0
    # Diagnostics
    d_eff: float = 0.0
    n_layers: int = 0
    mean_layer_gap: float = 0.0
    comp_frac: float = 0.0
    aw: int = 0


# ---------------------------------------------------------------------------
# Component computations
# ---------------------------------------------------------------------------

def compute_log_H(poset: Poset, n_runs: int = 256) -> float:
    """Combinatorial entropy via SIS estimator (fast, works for all N)."""
    log_mean, _ = estimate_log_linear_extensions_sis(poset, n_runs=n_runs)
    return log_mean


def compute_pi_geo(poset: Poset) -> float:
    """Geometric constraint penalty — neutral + geometric."""
    return neutral_penalty(poset) + geometric_penalty(poset)


def compute_sigma_hist(poset: Poset) -> float:
    """Historical sedimentation benefit — based on layer depth.

    Σ_hist = α₁·ℓ(X) + α₂·ḡ(X)
    where ℓ = layer count, ḡ = mean layer gap (normalized).
    """
    profile = layer_profile(poset)
    n_layers = len(profile)
    n = poset.n

    # Normalized layer count: fraction of N that are layers
    norm_layers = n_layers / max(n, 1)

    # Mean layer gap: average elements per layer (inverse = gap thinness)
    mean_size = profile.mean() if len(profile) > 0 else n
    # Smaller mean size = deeper layering = more sedimentation
    mean_gap = 1.0 / max(mean_size, 1.0)

    alpha1 = 1.0
    alpha2 = 0.5
    return alpha1 * norm_layers + alpha2 * mean_gap


def compute_xi_dim(poset: Poset) -> tuple[float, float]:
    """Dimensional barrier — cost of being far from stable dimension.

    Uses non-target-anchored consistency penalty: measures how self-consistent
    the structure's dimension estimate is, without pulling toward a fixed target.
    """
    total_pen, d_eff, _, _, _ = dimension_consistency_penalty(poset)
    return total_pen, d_eff


def compute_pi_cg(
    poset: Poset,
    keep_ratio: float = 0.7,
    n_cg_samples: int = 5,
    n_sis_samples: int = 200,
) -> float:
    """Coarse-graining drift penalty.

    Measures how much the observable signature changes after random node deletion.
    """
    if poset.n < 6:
        return 0.0

    # Compute original observables
    orig_comp = comparable_fraction(poset)
    orig_hr = height_ratio(poset)
    orig_wr = width_ratio(poset)

    drifts = []
    for i in range(n_cg_samples):
        cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=keep_ratio, seed=42 + i)
        if cg_poset.n < 3:
            continue
        cg_comp = comparable_fraction(cg_poset)
        cg_hr = height_ratio(cg_poset)
        cg_wr = width_ratio(cg_poset)

        drift = math.sqrt(
            (cg_comp - orig_comp) ** 2
            + (cg_hr - orig_hr) ** 2
            + (cg_wr - orig_wr) ** 2
        )
        drifts.append(drift)

    return float(np.mean(drifts)) if drifts else 0.0


# ---------------------------------------------------------------------------
# Unified functional
# ---------------------------------------------------------------------------

def compute_unified_functional(
    poset: Poset,
    weights: FunctionalWeights | None = None,
) -> FunctionalComponents:
    """Compute the full unified structural functional F[X].

    F[X] = β·log_H + γ·Π_geo − λ·Σ_hist + η·Ξ_d + κ·Π_cg
    """
    if weights is None:
        weights = FunctionalWeights()

    log_H = compute_log_H(poset)
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim, d_eff = compute_xi_dim(poset)
    pi_cg = compute_pi_cg(poset)

    F_total = (
        weights.beta * log_H
        + weights.gamma * pi_geo
        - weights.lam * sigma_hist
        + weights.eta * xi_dim
        + weights.kappa * pi_cg
    )

    profile = layer_profile(poset)
    n_layers = len(profile)
    mean_size = float(profile.mean()) if len(profile) > 0 else poset.n
    mean_gap = 1.0 / max(mean_size, 1.0)

    return FunctionalComponents(
        log_H=log_H,
        pi_geo=pi_geo,
        sigma_hist=sigma_hist,
        xi_dim=xi_dim,
        pi_cg=pi_cg,
        F_total=F_total,
        d_eff=d_eff,
        n_layers=n_layers,
        mean_layer_gap=mean_gap,
        comp_frac=comparable_fraction(poset),
        aw=antichain_width(poset),
    )


# ---------------------------------------------------------------------------
# Projection experiments: reproduce each inference
# ---------------------------------------------------------------------------

def experiment_prediction_B(
    n: int,
    reps: int,
    weights: FunctionalWeights,
    seed: int = 0,
) -> list[dict]:
    """Prediction B: Lorentzian-like vs KR-like phase competition.

    Shows F[X_Lor] < F[X_KR] when geometric constraint is active.
    """
    results = []
    for rep in range(reps):
        s = seed + rep * 1000
        lor = generate_lorentzian_like_2d(n, seed=s)
        kr = generate_kr_like(n, seed=s)

        fc_lor = compute_unified_functional(lor, weights)
        fc_kr = compute_unified_functional(kr, weights)

        results.append({
            "rep": rep,
            "family": "Lor2D",
            "N": n,
            "F_total": fc_lor.F_total,
            "log_H": fc_lor.log_H,
            "pi_geo": fc_lor.pi_geo,
            "sigma_hist": fc_lor.sigma_hist,
            "xi_dim": fc_lor.xi_dim,
            "pi_cg": fc_lor.pi_cg,
            "d_eff": fc_lor.d_eff,
            "comp_frac": fc_lor.comp_frac,
        })
        results.append({
            "rep": rep,
            "family": "KR_like",
            "N": n,
            "F_total": fc_kr.F_total,
            "log_H": fc_kr.log_H,
            "pi_geo": fc_kr.pi_geo,
            "sigma_hist": fc_kr.sigma_hist,
            "xi_dim": fc_kr.xi_dim,
            "pi_cg": fc_kr.pi_cg,
            "d_eff": fc_kr.d_eff,
            "comp_frac": fc_kr.comp_frac,
        })
    return results


def experiment_prediction_C(
    n: int,
    reps: int,
    weights: FunctionalWeights,
    seed: int = 0,
) -> list[dict]:
    """Prediction C: deeper hierarchy → lower effective cost.

    Compares Lorentzian-like posets and checks correlation between
    layer depth (Σ_hist) and log_H.
    """
    results = []
    for rep in range(reps):
        s = seed + rep * 1000
        poset = generate_lorentzian_like_2d(n, seed=s)
        fc = compute_unified_functional(poset, weights)
        results.append({
            "rep": rep,
            "family": "Lor2D",
            "N": n,
            "n_layers": fc.n_layers,
            "sigma_hist": fc.sigma_hist,
            "log_H": fc.log_H,
            "F_total": fc.F_total,
            "mean_layer_gap": fc.mean_layer_gap,
        })
    return results


def experiment_prediction_A(
    n: int,
    reps: int,
    weights: FunctionalWeights,
    seed: int = 0,
) -> list[dict]:
    """Prediction A: dimensional selection — 4D as local minimum.

    Compares F[X] across Lorentzian-like posets of different dimensions.
    """
    generators = {
        "Lor2D": generate_lorentzian_like_2d,
        "Lor3D": generate_lorentzian_like_3d,
        "Lor4D": generate_lorentzian_like_4d,
        "Lor5D": generate_lorentzian_like_5d,
    }
    results = []
    for rep in range(reps):
        s = seed + rep * 1000
        for name, gen in generators.items():
            poset = gen(n, seed=s)
            fc = compute_unified_functional(poset, weights)
            results.append({
                "rep": rep,
                "family": name,
                "N": n,
                "d_eff": fc.d_eff,
                "F_total": fc.F_total,
                "log_H": fc.log_H,
                "pi_geo": fc.pi_geo,
                "xi_dim": fc.xi_dim,
                "comp_frac": fc.comp_frac,
            })
    return results


def experiment_prediction_D(
    n: int,
    reps: int,
    weights: FunctionalWeights,
    seed: int = 0,
) -> list[dict]:
    """Prediction D: coarse-graining stability.

    For each structure, measures F before and after coarse-graining.
    Stable structures should have |ΔF| small.
    """
    generators = {
        "Lor2D": generate_lorentzian_like_2d,
        "KR_like": generate_kr_like,
        "Lor4D": generate_lorentzian_like_4d,
    }
    results = []
    for rep in range(reps):
        s = seed + rep * 1000
        for name, gen in generators.items():
            poset = gen(n, seed=s)
            fc_orig = compute_unified_functional(poset, weights)

            cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=0.7, seed=s + 99)
            fc_cg = compute_unified_functional(cg_poset, weights)

            delta_F = fc_cg.F_total - fc_orig.F_total
            results.append({
                "rep": rep,
                "family": name,
                "N": n,
                "F_orig": fc_orig.F_total,
                "F_cg": fc_cg.F_total,
                "delta_F": delta_F,
                "abs_delta_F": abs(delta_F),
                "pi_cg": fc_orig.pi_cg,
                "d_eff_orig": fc_orig.d_eff,
                "d_eff_cg": fc_cg.d_eff,
            })
    return results


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


def print_summary_B(rows: list[dict]) -> None:
    lor = [r for r in rows if r["family"] == "Lor2D"]
    kr = [r for r in rows if r["family"] == "KR_like"]
    if not lor or not kr:
        return

    mean_F_lor = np.mean([r["F_total"] for r in lor])
    mean_F_kr = np.mean([r["F_total"] for r in kr])
    lor_wins = sum(1 for l, k in zip(lor, kr) if l["F_total"] < k["F_total"])

    print(f"\n=== Prediction B: Phase Competition (N={lor[0]['N']}) ===")
    print(f"  mean F[Lor2D] = {mean_F_lor:.4f}")
    print(f"  mean F[KR]    = {mean_F_kr:.4f}")
    print(f"  Lor2D wins    = {lor_wins}/{len(lor)} ({100*lor_wins/len(lor):.0f}%)")
    if mean_F_lor < mean_F_kr:
        print(f"  → B CONFIRMED: Lor2D has lower structural cost (ΔF = {mean_F_kr - mean_F_lor:.4f})")
    else:
        print(f"  → B not confirmed at these weights (ΔF = {mean_F_kr - mean_F_lor:.4f})")


def print_summary_C(rows: list[dict]) -> None:
    if not rows:
        return
    layers = np.array([r["n_layers"] for r in rows])
    log_Hs = np.array([r["log_H"] for r in rows])
    if len(layers) < 3 or layers.std() < 1e-8:
        print("\n=== Prediction C: too little variance in layers ===")
        return

    r = np.corrcoef(layers, log_Hs)[0, 1]
    print(f"\n=== Prediction C: Historical Sedimentation (N={rows[0]['N']}) ===")
    print(f"  corr(n_layers, log_H) = {r:.4f}")
    print(f"  layer range: {int(layers.min())}–{int(layers.max())}")
    if r < -0.1:
        print(f"  → C direction CONFIRMED: deeper layers correlate with lower entropy")
    else:
        print(f"  → C direction not confirmed (r = {r:.4f})")


def print_summary_A(rows: list[dict]) -> None:
    if not rows:
        return
    families = sorted(set(r["family"] for r in rows))
    print(f"\n=== Prediction A: Dimensional Selection (N={rows[0]['N']}) ===")

    mean_Fs = {}
    for fam in families:
        fam_rows = [r for r in rows if r["family"] == fam]
        mean_F = np.mean([r["F_total"] for r in fam_rows])
        mean_d = np.mean([r["d_eff"] for r in fam_rows])
        mean_Fs[fam] = mean_F
        print(f"  {fam:8s}: mean F = {mean_F:.4f}, mean d_eff = {mean_d:.2f}")

    # Prediction A's real content: asymmetric barrier Ξ_{d→d+1}
    # The marginal cost of going one dimension higher should INCREASE at d=4
    dim_order = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]
    present = [d for d in dim_order if d in mean_Fs]
    if len(present) >= 3:
        print(f"\n  Marginal cost ΔF per dimension step:")
        deltas = []
        for i in range(len(present) - 1):
            delta = mean_Fs[present[i + 1]] - mean_Fs[present[i]]
            deltas.append((present[i], present[i + 1], delta))
            print(f"    {present[i]} → {present[i+1]}: ΔF = {delta:+.4f}")

        # Check asymmetric barrier: Ξ_{4→5} > Ξ_{3→4}
        if len(deltas) >= 3:
            d3to4 = deltas[1][2]  # Lor3D → Lor4D
            d4to5 = deltas[2][2]  # Lor4D → Lor5D
            if d4to5 > d3to4:
                print(f"  → A barrier asymmetry: Ξ(4→5) = {d4to5:.4f} > Ξ(3→4) = {d3to4:.4f}")
                print(f"    Barrier INCREASES at d=4 — consistent with 4D locking")
            else:
                print(f"  → A barrier asymmetry not found: Ξ(4→5) = {d4to5:.4f} ≤ Ξ(3→4) = {d3to4:.4f}")

        # Also check if marginal cost is decelerating toward 4D
        if len(deltas) >= 2:
            d2to3 = deltas[0][2]
            d3to4 = deltas[1][2]
            if d3to4 < d2to3:
                print(f"    Marginal cost DECELERATING into 4D (ΔF shrinks): {d2to3:.4f} → {d3to4:.4f}")
            else:
                print(f"    Marginal cost accelerating into 4D: {d2to3:.4f} → {d3to4:.4f}")


def print_summary_D(rows: list[dict]) -> None:
    if not rows:
        return
    families = sorted(set(r["family"] for r in rows))
    print(f"\n=== Prediction D: Coarse-Graining Stability (N={rows[0]['N']}) ===")
    for fam in families:
        fam_rows = [r for r in rows if r["family"] == fam]
        mean_abs_dF = np.mean([r["abs_delta_F"] for r in fam_rows])
        mean_pi_cg = np.mean([r["pi_cg"] for r in fam_rows])
        print(f"  {fam:8s}: mean |ΔF| = {mean_abs_dF:.4f}, mean Π_cg = {mean_pi_cg:.4f}")

    # Check if Lor has lower drift
    lor_drift = np.mean([r["abs_delta_F"] for r in rows if r["family"] == "Lor2D"])
    kr_drift = np.mean([r["abs_delta_F"] for r in rows if r["family"] == "KR_like"])
    if lor_drift < kr_drift:
        print(f"  → D direction CONFIRMED: Lor2D more stable under coarse-graining")
    else:
        print(f"  → D direction not confirmed")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Unified Structural Functional — Toy Model")
    parser.add_argument("--n", type=int, default=20, help="Number of poset elements")
    parser.add_argument("--reps", type=int, default=10, help="Repetitions per experiment")
    parser.add_argument("--seed", type=int, default=42, help="Base random seed")
    parser.add_argument("--outdir", type=str, default="outputs_unified_functional",
                        help="Output directory for CSV results")
    # Weight overrides
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--lam", type=float, default=0.5)
    parser.add_argument("--eta", type=float, default=0.3)
    parser.add_argument("--kappa", type=float, default=0.2)
    args = parser.parse_args()

    weights = FunctionalWeights(
        beta=args.beta,
        gamma=args.gamma,
        lam=args.lam,
        eta=args.eta,
        kappa=args.kappa,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Unified Structural Functional — N={args.n}, reps={args.reps}")
    print(f"Weights: β={weights.beta}, γ={weights.gamma}, λ={weights.lam}, η={weights.eta}, κ={weights.kappa}")
    print()

    # --- Prediction B ---
    print("Running Prediction B (phase competition)...")
    rows_B = experiment_prediction_B(args.n, args.reps, weights, seed=args.seed)
    save_csv(rows_B, str(outdir / "prediction_B_phase_competition.csv"))
    print_summary_B(rows_B)

    # --- Prediction C ---
    print("\nRunning Prediction C (historical sedimentation)...")
    rows_C = experiment_prediction_C(args.n, args.reps * 3, weights, seed=args.seed + 10000)
    save_csv(rows_C, str(outdir / "prediction_C_sedimentation.csv"))
    print_summary_C(rows_C)

    # --- Prediction A ---
    print("\nRunning Prediction A (dimensional selection)...")
    rows_A = experiment_prediction_A(args.n, args.reps, weights, seed=args.seed + 20000)
    save_csv(rows_A, str(outdir / "prediction_A_dimension.csv"))
    print_summary_A(rows_A)

    # --- Prediction D ---
    print("\nRunning Prediction D (coarse-graining stability)...")
    rows_D = experiment_prediction_D(args.n, args.reps, weights, seed=args.seed + 30000)
    save_csv(rows_D, str(outdir / "prediction_D_cg_stability.csv"))
    print_summary_D(rows_D)

    print("\n" + "=" * 60)
    print("All four projections computed. Results saved to:", outdir)
    print("=" * 60)


if __name__ == "__main__":
    main()
