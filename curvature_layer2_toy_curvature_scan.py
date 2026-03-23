"""Conjecture E - Layer 2 toy-curvature scan.

This script introduces a minimal curved baseline for the second layer of
Conjecture E.

We use a conformally flat toy metric on the unit causal diamond:

    ds^2 = Omega(t)^2 (-dt^2 + dx^2 + ...)

with Omega(t) = exp(kappa * (t - 1/2)^2).

This is not claimed to be a full de Sitter or Schwarzschild construction.
It is a controlled, positive conformal deformation that:
  - preserves causal order,
  - changes the sampling density,
  - induces systematic deviations in interval statistics.

The purpose is to verify that the layer-2 pipeline can detect a tunable
"curvature-like" signal before we attempt a more realistic geometry.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from curvature_layer2_baseline import (
    causal_pair_data,
    fit_curvature_proxy,
    fit_flat_volume_law,
    poset_from_points,
)


def unit_ball_volume(m: int) -> float:
    return math.pi ** (m / 2.0) / math.gamma(m / 2.0 + 1.0)


def conformal_weight(t: np.ndarray, kappa: float) -> np.ndarray:
    """Positive conformal factor Omega(t)."""
    return np.exp(kappa * (t - 0.5) ** 2)


def sprinkle_toy_curved_diamond(
    n: int,
    d_spatial: int,
    kappa: float,
    seed: int | None = None,
) -> np.ndarray:
    """Rejection sample points with a conformal time-dependent density.

    The causal region is still the unit diamond in coordinate space. The
    measure is distorted by Omega(t)^d_spatial, so the point cloud is denser
    in time regions with larger conformal weight.
    """
    rng = np.random.default_rng(seed)
    points = np.empty((0, 1 + d_spatial), dtype=float)
    w_max = float(np.exp(kappa * 0.25))
    while points.shape[0] < n:
        batch_size = max(500, 20 * n)
        t = rng.random(batch_size)
        x = rng.random((batch_size, d_spatial)) - 0.5
        r2 = np.sum(x * x, axis=1)
        r_max = np.minimum(t, 1.0 - t)
        causal_region = r2 <= (r_max * r_max)
        if not np.any(causal_region):
            continue
        w = conformal_weight(t[causal_region], kappa) ** d_spatial
        accept = rng.random(w.size) < np.clip(w / w_max, 0.0, 1.0)
        accepted = np.column_stack([t[causal_region][accept], x[causal_region][accept]])
        points = np.vstack([points, accepted])
    return points[:n]


def interval_fraction(poset_n: int, k: np.ndarray) -> np.ndarray:
    return k / max(poset_n - 2, 1)


@dataclass(frozen=True)
class ToyRow:
    kappa: float
    d: int
    N: int
    reps: int
    tau_quantile: float
    n_causal_pairs: int
    flat_prefactor_fit: float
    flat_rmse: float
    flat_r2: float
    curvature_intercept: float
    curvature_r_hat: float
    mean_tau: float
    median_tau: float


def main() -> int:
    ap = argparse.ArgumentParser(description="Conjecture E layer-2 toy curvature scan")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 4])
    ap.add_argument("--n", type=int, default=192)
    ap.add_argument("--reps", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tau-quantile", type=float, default=0.75)
    ap.add_argument("--kappas", nargs="*", type=float, default=[0.0, 0.5, 1.0, 2.0])
    ap.add_argument("--out", default="outputs_unified_functional/curvature_layer2_toy_curvature_scan.csv")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[ToyRow] = []
    for d in args.dims:
        if d < 2:
            raise ValueError("d must be >= 2")
        for kappa in args.kappas:
            tau_parts: list[np.ndarray] = []
            k_parts: list[np.ndarray] = []
            for rep in range(args.reps):
                points = sprinkle_toy_curved_diamond(
                    args.n,
                    d - 1,
                    kappa=kappa,
                    seed=args.seed + d * 1000 + rep + int(kappa * 100),
                )
                poset = poset_from_points(points)
                tau, k = causal_pair_data(points, poset)
                if tau.size == 0:
                    continue
                tau_parts.append(tau)
                k_parts.append(k)

            if not tau_parts:
                raise RuntimeError(f"No causal pairs for d={d}, N={args.n}, kappa={kappa}")

            tau = np.concatenate(tau_parts)
            k = np.concatenate(k_parts)
            y = interval_fraction(args.n, k)
            x = np.power(tau, d)
            ratio = y / np.clip(x, 1e-12, None)
            fit_mask = tau <= np.quantile(tau, args.tau_quantile)
            if fit_mask.sum() < 10:
                fit_mask = np.ones_like(tau, dtype=bool)
            flat_prefactor, flat_rmse, flat_r2 = fit_flat_volume_law(tau[fit_mask], y[fit_mask], d)
            curvature_intercept, curvature_r_hat = fit_curvature_proxy(tau, ratio, d)

            rows.append(
                ToyRow(
                    kappa=kappa,
                    d=d,
                    N=args.n,
                    reps=args.reps,
                    tau_quantile=args.tau_quantile,
                    n_causal_pairs=int(tau.size),
                    flat_prefactor_fit=flat_prefactor,
                    flat_rmse=flat_rmse,
                    flat_r2=flat_r2,
                    curvature_intercept=curvature_intercept,
                    curvature_r_hat=curvature_r_hat,
                    mean_tau=float(np.mean(tau)),
                    median_tau=float(np.median(tau)),
                )
            )

    with out_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(ToyRow.__dataclass_fields__.keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row.__dict__)

    print(f"Saved: {out_path}")
    print("\nToy curvature summary:")
    for d in args.dims:
        print(f"\n  d={d}:")
        for row in [r for r in rows if r.d == d]:
            print(
                f"    kappa={row.kappa:>4.1f} pairs={row.n_causal_pairs:6d} "
                f"pref={row.flat_prefactor_fit:.4f} r2={row.flat_r2:.3f} "
                f"R_hat={row.curvature_r_hat:+.3f}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

