"""Unified Structural Functional — Weight Fitting

Fits optimal weights (β, γ, λ, η, κ) by generating a multi-N, multi-family
dataset and optimizing a composite loss that encodes all four inference criteria:

  L_B: Lor2D should beat KR_like (phase competition)
  L_C: corr(Σ_hist, log_H) should be negative within fixed N (sedimentation)
  L_A: marginal ΔF should show decelerating or saturating pattern near d=4
  L_D: Lor should have lower |ΔF_cg| than KR (coarse-graining stability)

Usage:
    python unified_functional_weight_fit.py [--n_list 16 20 28 36] [--reps 8]
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from unified_functional import (
    FunctionalWeights,
    compute_log_H,
    compute_pi_geo,
    compute_sigma_hist,
    compute_xi_dim,
    compute_pi_cg,
)
from generators import (
    Poset,
    generate_kr_like,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)
from observables import comparable_fraction, layer_profile
from observables_geo import height_ratio, width_ratio
from coarse_grain import coarse_grain_delete_nodes


# ---------------------------------------------------------------------------
# Raw feature extraction (weight-independent)
# ---------------------------------------------------------------------------

@dataclass
class RawFeatures:
    """Weight-independent features for a single poset."""
    family: str
    N: int
    rep: int
    seed: int
    log_H: float
    pi_geo: float
    sigma_hist: float
    xi_dim: float
    pi_cg: float
    d_eff: float
    n_layers: int
    comp_frac: float
    # Coarse-grained versions
    log_H_cg: float = 0.0
    pi_geo_cg: float = 0.0
    sigma_hist_cg: float = 0.0
    xi_dim_cg: float = 0.0
    pi_cg_cg: float = 0.0
    d_eff_cg: float = 0.0


def extract_features(
    family: str,
    gen_func,
    n: int,
    rep: int,
    seed: int,
    keep_ratio: float = 0.7,
) -> RawFeatures:
    """Extract all weight-independent features from a poset."""
    poset = gen_func(n, seed=seed)

    log_H = compute_log_H(poset)
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim, d_eff = compute_xi_dim(poset)
    pi_cg = compute_pi_cg(poset, keep_ratio=keep_ratio, n_cg_samples=3)

    profile = layer_profile(poset)
    n_layers = len(profile)
    comp_frac = comparable_fraction(poset)

    # Coarse-grained features
    cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=keep_ratio, seed=seed + 99)
    log_H_cg = compute_log_H(cg_poset)
    pi_geo_cg = compute_pi_geo(cg_poset)
    sigma_hist_cg = compute_sigma_hist(cg_poset)
    xi_dim_cg, d_eff_cg = compute_xi_dim(cg_poset)
    pi_cg_cg = compute_pi_cg(cg_poset, keep_ratio=keep_ratio, n_cg_samples=3)

    return RawFeatures(
        family=family, N=n, rep=rep, seed=seed,
        log_H=log_H, pi_geo=pi_geo, sigma_hist=sigma_hist,
        xi_dim=xi_dim, pi_cg=pi_cg, d_eff=d_eff,
        n_layers=n_layers, comp_frac=comp_frac,
        log_H_cg=log_H_cg, pi_geo_cg=pi_geo_cg,
        sigma_hist_cg=sigma_hist_cg, xi_dim_cg=xi_dim_cg,
        pi_cg_cg=pi_cg_cg, d_eff_cg=d_eff_cg,
    )


def F_from_features(f: RawFeatures, w: FunctionalWeights) -> float:
    """Compute F[X] from pre-extracted features and weights."""
    return (
        w.beta * f.log_H
        + w.gamma * f.pi_geo
        - w.lam * f.sigma_hist
        + w.eta * f.xi_dim
        + w.kappa * f.pi_cg
    )


def F_cg_from_features(f: RawFeatures, w: FunctionalWeights) -> float:
    """Compute F[R(X)] from pre-extracted coarse-grained features."""
    return (
        w.beta * f.log_H_cg
        + w.gamma * f.pi_geo_cg
        - w.lam * f.sigma_hist_cg
        + w.eta * f.xi_dim_cg
        + w.kappa * f.pi_cg_cg
    )


# ---------------------------------------------------------------------------
# Loss functions for each inference
# ---------------------------------------------------------------------------

def loss_B(features: list[RawFeatures], w: FunctionalWeights) -> float:
    """Prediction B: Lor2D should have lower F than KR_like.

    Loss = mean over (N, rep) pairs of max(0, F[Lor2D] - F[KR] + margin).
    We want F[Lor2D] < F[KR], so penalize when Lor2D >= KR.
    """
    margin = 0.5
    total_loss = 0.0
    count = 0
    for n_val in set(f.N for f in features):
        lor_feats = [f for f in features if f.family == "Lor2D" and f.N == n_val]
        kr_feats = [f for f in features if f.family == "KR_like" and f.N == n_val]
        for lf in lor_feats:
            F_lor = F_from_features(lf, w)
            for kf in kr_feats:
                if kf.rep == lf.rep:
                    F_kr = F_from_features(kf, w)
                    violation = F_lor - F_kr + margin
                    total_loss += max(0.0, violation) ** 2
                    count += 1
    return total_loss / max(count, 1)


def loss_C(features: list[RawFeatures], w: FunctionalWeights) -> float:
    """Prediction C: within fixed N, deeper layers → lower F.

    Loss = -corr(sigma_hist, F) across Lor2D samples at each N.
    We want negative correlation (more sedimentation → lower F),
    so loss is high when correlation is positive.
    """
    total_loss = 0.0
    count = 0
    for n_val in set(f.N for f in features):
        lor_feats = [f for f in features if f.family == "Lor2D" and f.N == n_val]
        if len(lor_feats) < 4:
            continue
        sigs = np.array([f.sigma_hist for f in lor_feats])
        Fs = np.array([F_from_features(f, w) for f in lor_feats])
        if sigs.std() < 1e-8 or Fs.std() < 1e-8:
            continue
        r = np.corrcoef(sigs, Fs)[0, 1]
        # We want r < 0 (more sedimentation → lower F)
        # Penalize r > -0.3
        target = -0.3
        if r > target:
            total_loss += (r - target) ** 2
        count += 1
    return total_loss / max(count, 1)


def loss_A(features: list[RawFeatures], w: FunctionalWeights) -> float:
    """Prediction A: marginal ΔF should decelerate approaching d=4.

    The marginal cost of going from d to d+1 should decrease or stabilize
    near d=4, creating a "diminishing returns" pattern that locks 4D.
    
    Specifically: we want ΔF(3→4) < ΔF(2→3) (deceleration into 4D)
    This is the signature of an approaching local minimum in the marginal cost.
    """
    dim_families = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]
    total_loss = 0.0
    count = 0

    for n_val in set(f.N for f in features):
        mean_Fs = {}
        for fam in dim_families:
            fam_feats = [f for f in features if f.family == fam and f.N == n_val]
            if fam_feats:
                mean_Fs[fam] = np.mean([F_from_features(f, w) for f in fam_feats])

        if len(mean_Fs) >= 3:
            present = [d for d in dim_families if d in mean_Fs]
            deltas = []
            for i in range(len(present) - 1):
                deltas.append(mean_Fs[present[i + 1]] - mean_Fs[present[i]])

            # Want deceleration: delta[1] < delta[0] (ΔF shrinks approaching 4D)
            if len(deltas) >= 2:
                accel = deltas[1] - deltas[0]  # negative = deceleration
                if accel > 0:  # accelerating = bad
                    total_loss += accel ** 2
                count += 1

            # Also want all deltas positive (higher dim = higher cost)
            for d in deltas:
                if d < 0:
                    total_loss += d ** 2
                    count += 1

    return total_loss / max(count, 1)


def loss_D(features: list[RawFeatures], w: FunctionalWeights) -> float:
    """Prediction D: Lor2D should be more stable under coarse-graining than KR.

    Loss penalizes when |ΔF_cg| for Lor2D >= |ΔF_cg| for KR.
    """
    total_loss = 0.0
    count = 0

    for n_val in set(f.N for f in features):
        lor_feats = [f for f in features if f.family == "Lor2D" and f.N == n_val]
        kr_feats = [f for f in features if f.family == "KR_like" and f.N == n_val]

        if lor_feats and kr_feats:
            lor_drifts = []
            for lf in lor_feats:
                F_orig = F_from_features(lf, w)
                F_cg = F_cg_from_features(lf, w)
                lor_drifts.append(abs(F_cg - F_orig))

            kr_drifts = []
            for kf in kr_feats:
                F_orig = F_from_features(kf, w)
                F_cg = F_cg_from_features(kf, w)
                kr_drifts.append(abs(F_cg - F_orig))

            mean_lor = np.mean(lor_drifts)
            mean_kr = np.mean(kr_drifts)

            # Want mean_lor < mean_kr
            if mean_lor >= mean_kr:
                total_loss += (mean_lor - mean_kr + 0.1) ** 2
            count += 1

    return total_loss / max(count, 1)


def composite_loss(
    features: list[RawFeatures],
    w: FunctionalWeights,
    w_B: float = 1.0,
    w_C: float = 1.0,
    w_A: float = 0.5,
    w_D: float = 0.5,
) -> tuple[float, dict[str, float]]:
    """Composite loss combining all four inference criteria."""
    lB = loss_B(features, w)
    lC = loss_C(features, w)
    lA = loss_A(features, w)
    lD = loss_D(features, w)

    total = w_B * lB + w_C * lC + w_A * lA + w_D * lD
    return total, {"L_B": lB, "L_C": lC, "L_A": lA, "L_D": lD, "L_total": total}


# ---------------------------------------------------------------------------
# Grid search
# ---------------------------------------------------------------------------

def grid_search(
    features: list[RawFeatures],
    beta_range: list[float],
    gamma_range: list[float],
    lam_range: list[float],
    eta_range: list[float],
    kappa_range: list[float],
) -> tuple[FunctionalWeights, float, dict, list[dict]]:
    """Exhaustive grid search over weight combinations."""
    best_w = None
    best_loss = float("inf")
    best_components = {}
    all_results = []

    total_combos = (
        len(beta_range) * len(gamma_range) * len(lam_range)
        * len(eta_range) * len(kappa_range)
    )
    print(f"Grid search: {total_combos} combinations")

    for i, (b, g, l, e, k) in enumerate(
        itertools.product(beta_range, gamma_range, lam_range, eta_range, kappa_range)
    ):
        w = FunctionalWeights(beta=b, gamma=g, lam=l, eta=e, kappa=k)
        loss, components = composite_loss(features, w)

        row = {"beta": b, "gamma": g, "lam": l, "eta": e, "kappa": k}
        row.update(components)
        all_results.append(row)

        if loss < best_loss:
            best_loss = loss
            best_w = w
            best_components = components

        if (i + 1) % 500 == 0:
            print(f"  [{i+1}/{total_combos}] best loss so far: {best_loss:.6f}")

    return best_w, best_loss, best_components, all_results


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_weights(features: list[RawFeatures], w: FunctionalWeights) -> None:
    """Run all four projection checks with the given weights."""
    print(f"\n{'='*60}")
    print(f"Validating weights: β={w.beta}, γ={w.gamma}, λ={w.lam}, η={w.eta}, κ={w.kappa}")
    print(f"{'='*60}")

    # --- B ---
    print(f"\n--- Prediction B: Phase Competition ---")
    for n_val in sorted(set(f.N for f in features)):
        lor = [f for f in features if f.family == "Lor2D" and f.N == n_val]
        kr = [f for f in features if f.family == "KR_like" and f.N == n_val]
        if not lor or not kr:
            continue
        wins = 0
        for lf in lor:
            F_l = F_from_features(lf, w)
            matched_kr = [kf for kf in kr if kf.rep == lf.rep]
            if matched_kr:
                F_k = F_from_features(matched_kr[0], w)
                if F_l < F_k:
                    wins += 1
        mean_F_lor = np.mean([F_from_features(f, w) for f in lor])
        mean_F_kr = np.mean([F_from_features(f, w) for f in kr])
        print(f"  N={n_val:3d}: F[Lor]={mean_F_lor:8.2f}, F[KR]={mean_F_kr:8.2f}, "
              f"Lor wins={wins}/{len(lor)}, ΔF={mean_F_kr-mean_F_lor:+.2f}")

    # --- C ---
    print(f"\n--- Prediction C: Sedimentation ---")
    for n_val in sorted(set(f.N for f in features)):
        lor = [f for f in features if f.family == "Lor2D" and f.N == n_val]
        if len(lor) < 4:
            continue
        sigs = np.array([f.sigma_hist for f in lor])
        Fs = np.array([F_from_features(f, w) for f in lor])
        if sigs.std() < 1e-8 or Fs.std() < 1e-8:
            continue
        r = np.corrcoef(sigs, Fs)[0, 1]
        r_logH = np.corrcoef(
            np.array([f.n_layers for f in lor]),
            np.array([f.log_H for f in lor])
        )[0, 1]
        print(f"  N={n_val:3d}: corr(Σ_hist, F)={r:+.4f}, corr(layers, logH)={r_logH:+.4f}")

    # --- A ---
    print(f"\n--- Prediction A: Dimensional Selection ---")
    dim_fams = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]
    for n_val in sorted(set(f.N for f in features)):
        line = f"  N={n_val:3d}: "
        mean_Fs = {}
        for fam in dim_fams:
            ff = [f for f in features if f.family == fam and f.N == n_val]
            if ff:
                mF = np.mean([F_from_features(f, w) for f in ff])
                mean_Fs[fam] = mF
                line += f"{fam}={mF:.1f} "
        print(line)
        present = [d for d in dim_fams if d in mean_Fs]
        if len(present) >= 2:
            deltas_str = "    ΔF: "
            for i in range(len(present) - 1):
                d = mean_Fs[present[i + 1]] - mean_Fs[present[i]]
                deltas_str += f"{present[i][3:]}→{present[i+1][3:]}={d:+.1f} "
            print(deltas_str)

    # --- D ---
    print(f"\n--- Prediction D: Coarse-Graining Stability ---")
    for n_val in sorted(set(f.N for f in features)):
        for fam in ["Lor2D", "KR_like", "Lor4D"]:
            ff = [f for f in features if f.family == fam and f.N == n_val]
            if not ff:
                continue
            drifts = [abs(F_cg_from_features(f, w) - F_from_features(f, w)) for f in ff]
            print(f"  N={n_val:3d} {fam:8s}: mean |ΔF_cg|={np.mean(drifts):.4f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Unified Functional — Weight Fitting")
    parser.add_argument("--n_list", type=int, nargs="+", default=[16, 20, 28],
                        help="List of N values to sample")
    parser.add_argument("--reps", type=int, default=8, help="Repetitions per (family, N)")
    parser.add_argument("--seed", type=int, default=42, help="Base random seed")
    parser.add_argument("--outdir", type=str, default="outputs_unified_functional",
                        help="Output directory")
    parser.add_argument("--grid_resolution", type=int, default=4,
                        help="Number of grid points per weight dimension")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Step 1: Generate dataset ---
    generators = {
        "Lor2D": generate_lorentzian_like_2d,
        "Lor3D": generate_lorentzian_like_3d,
        "Lor4D": generate_lorentzian_like_4d,
        "Lor5D": generate_lorentzian_like_5d,
        "KR_like": generate_kr_like,
    }

    print(f"Step 1: Generating dataset")
    print(f"  N values: {args.n_list}")
    print(f"  Families: {list(generators.keys())}")
    print(f"  Reps per cell: {args.reps}")

    all_features: list[RawFeatures] = []
    total_cells = len(args.n_list) * len(generators) * args.reps
    done = 0

    for n_val in args.n_list:
        for fam_name, gen_func in generators.items():
            for rep in range(args.reps):
                s = args.seed + rep * 1000 + n_val * 100
                feat = extract_features(fam_name, gen_func, n_val, rep, s)
                all_features.append(feat)
                done += 1
                if done % 20 == 0:
                    print(f"  [{done}/{total_cells}] extracted")

    print(f"  Total: {len(all_features)} feature vectors")

    # Save raw features
    feat_rows = []
    for f in all_features:
        feat_rows.append({
            "family": f.family, "N": f.N, "rep": f.rep,
            "log_H": f.log_H, "pi_geo": f.pi_geo, "sigma_hist": f.sigma_hist,
            "xi_dim": f.xi_dim, "pi_cg": f.pi_cg, "d_eff": f.d_eff,
            "n_layers": f.n_layers, "comp_frac": f.comp_frac,
            "log_H_cg": f.log_H_cg, "pi_geo_cg": f.pi_geo_cg,
            "sigma_hist_cg": f.sigma_hist_cg, "xi_dim_cg": f.xi_dim_cg,
            "pi_cg_cg": f.pi_cg_cg, "d_eff_cg": f.d_eff_cg,
        })
    feat_path = str(outdir / "raw_features.csv")
    with open(feat_path, "w", newline="", encoding="utf-8") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=list(feat_rows[0].keys()))
        writer.writeheader()
        writer.writerows(feat_rows)
    print(f"  → saved {feat_path}")

    # --- Step 2: Grid search ---
    print(f"\nStep 2: Grid search over weights")
    res = args.grid_resolution
    beta_range = np.linspace(0.5, 2.0, res).tolist()
    gamma_range = np.linspace(0.5, 3.0, res).tolist()
    lam_range = np.linspace(0.1, 1.5, res).tolist()
    eta_range = np.linspace(0.1, 1.0, res).tolist()
    kappa_range = np.linspace(0.05, 0.5, res).tolist()

    best_w, best_loss, best_comp, all_grid = grid_search(
        all_features, beta_range, gamma_range, lam_range, eta_range, kappa_range,
    )

    print(f"\n  BEST weights: β={best_w.beta:.3f}, γ={best_w.gamma:.3f}, "
          f"λ={best_w.lam:.3f}, η={best_w.eta:.3f}, κ={best_w.kappa:.3f}")
    print(f"  BEST loss: {best_loss:.6f}")
    print(f"  Components: {best_comp}")

    # Save grid results
    grid_path = str(outdir / "weight_grid_search.csv")
    with open(grid_path, "w", newline="", encoding="utf-8") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=list(all_grid[0].keys()))
        writer.writeheader()
        writer.writerows(all_grid)
    print(f"  → saved {grid_path}")

    # --- Step 3: Validate ---
    print(f"\nStep 3: Validate best weights")
    validate_weights(all_features, best_w)

    # Also validate with default weights for comparison
    print(f"\n--- Comparison: Default weights ---")
    default_w = FunctionalWeights()
    _, default_comp = composite_loss(all_features, default_w)
    print(f"  Default loss: {default_comp['L_total']:.6f}")
    print(f"  Default components: {default_comp}")
    _, best_comp2 = composite_loss(all_features, best_w)
    print(f"  Optimized loss: {best_comp2['L_total']:.6f}")
    improvement = (default_comp['L_total'] - best_comp2['L_total']) / max(default_comp['L_total'], 1e-8) * 100
    print(f"  Improvement: {improvement:.1f}%")

    # Save summary
    summary_path = str(outdir / "weight_fit_summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=[
            "type", "beta", "gamma", "lam", "eta", "kappa",
            "L_B", "L_C", "L_A", "L_D", "L_total",
        ])
        writer.writeheader()
        writer.writerow({
            "type": "default",
            "beta": 1.0, "gamma": 1.0, "lam": 0.5, "eta": 0.3, "kappa": 0.2,
            **default_comp,
        })
        writer.writerow({
            "type": "optimized",
            "beta": best_w.beta, "gamma": best_w.gamma,
            "lam": best_w.lam, "eta": best_w.eta, "kappa": best_w.kappa,
            **best_comp2,
        })
    print(f"  → saved {summary_path}")

    print(f"\n{'='*60}")
    print(f"Weight fitting complete.")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
