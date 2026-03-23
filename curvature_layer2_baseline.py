"""Conjecture E - Layer 2 baseline: flat-space interval-volume / curvature proxy.

This script builds the smallest coordinate-aware experiment needed for the
second layer of Conjecture E:

1. Sprinkle points uniformly in a causal diamond.
2. Build the induced poset from the Minkowski causal order.
3. For each causal pair, compute:
   - proper time tau
   - interval size k = |I(x, y)|
4. Fit the flat-space Alexandrov-volume law

      E[k/(N-2)] ~= kappa_d * tau^d

   and extract a first curvature proxy from the residual term

      k/(N-2) / tau^d ~= 1 - c_d * R_hat * tau^2.

The intent is not to prove a curved-spacetime theorem. The goal is to
establish a clean, coordinate-aware pipeline that can later be reused on
curved sprinklings once a curved generator is added.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from generators import Poset, transitive_closure


def unit_ball_volume(m: int) -> float:
    return math.pi ** (m / 2.0) / math.gamma(m / 2.0 + 1.0)


def alexandrov_volume_constant(d: int) -> float:
    """Flat d-dimensional Alexandrov interval volume prefactor."""
    return unit_ball_volume(d - 1) / (2 ** (d - 1) * d)


def sprinkle_causal_diamond(n: int, d_spatial: int, seed: int | None = None) -> np.ndarray:
    """Sample points uniformly in a unit causal diamond.

    The diamond is J^+(0) ∩ J^-(1) in (1 + d_spatial)-dimensional Minkowski
    space, represented in coordinates [t, x1, ..., xd].
    """
    rng = np.random.default_rng(seed)
    points = np.empty((0, 1 + d_spatial), dtype=float)
    while points.shape[0] < n:
        batch_size = max(500, 20 * n)
        t = rng.random(batch_size)
        x = rng.random((batch_size, d_spatial)) - 0.5
        r2 = np.sum(x * x, axis=1)
        r_max = np.minimum(t, 1.0 - t)
        accept = r2 <= (r_max * r_max)
        accepted = np.column_stack([t[accept], x[accept]])
        points = np.vstack([points, accepted])
    return points[:n]


def poset_from_points(points: np.ndarray) -> Poset:
    t = points[:, 0]
    spatial = points[:, 1:]
    dt = t[None, :] - t[:, None]
    if spatial.shape[1] == 0:
        spatial_d2 = 0.0
    else:
        spatial_d2 = np.sum((spatial[None, :, :] - spatial[:, None, :]) ** 2, axis=2)
    adj = (dt > 0.0) & (dt * dt >= spatial_d2)
    return Poset(transitive_closure(adj))


def interval_sizes(poset: Poset) -> np.ndarray:
    c = poset.closure.astype(np.int32)
    return c @ c


def causal_pair_data(points: np.ndarray, poset: Poset) -> tuple[np.ndarray, np.ndarray]:
    """Return tau and interval size k for all causal pairs."""
    t = points[:, 0]
    spatial = points[:, 1:]
    dt = t[None, :] - t[:, None]
    if spatial.shape[1] == 0:
        spatial_d2 = 0.0
    else:
        spatial_d2 = np.sum((spatial[None, :, :] - spatial[:, None, :]) ** 2, axis=2)
    causal = poset.closure
    tau2 = np.clip(dt * dt - spatial_d2, 0.0, None)
    tau = np.sqrt(tau2[causal])
    k = interval_sizes(poset)[causal].astype(float)
    return tau, k


def fit_flat_volume_law(tau: np.ndarray, y: np.ndarray, d: int) -> tuple[float, float, float]:
    """Fit y ~= a * tau^d and return (a, rmse, r2)."""
    x = np.power(tau, d)
    a = float(np.dot(x, y) / max(np.dot(x, x), 1e-12))
    pred = a * x
    resid = y - pred
    rmse = float(np.sqrt(np.mean(resid * resid)))
    denom = float(np.sum((y - float(np.mean(y))) ** 2))
    r2 = 1.0 - float(np.sum(resid * resid)) / denom if denom > 0 else float("nan")
    return a, rmse, r2


def fit_curvature_proxy(tau: np.ndarray, ratio: np.ndarray, d: int) -> tuple[float, float]:
    """Fit ratio ~= 1 - c * tau^2 and map c to an R_hat proxy.

    For the flat-space expectation, ratio should be close to 1 and the fitted
    slope should be near zero.
    """
    x = tau * tau
    y = ratio
    A = np.c_[np.ones_like(x), x]
    intercept, slope = np.linalg.lstsq(A, y, rcond=None)[0]
    coeff = 6.0 * (d + 1) * (d + 2)
    r_hat = -slope * coeff
    return float(intercept), float(r_hat)


@dataclass(frozen=True)
class Summary:
    d: int
    N: int
    seed: int
    n_causal_pairs: int
    mean_tau: float
    mean_interval_size: float
    flat_prefactor_fit: float
    flat_rmse: float
    flat_r2: float
    curvature_intercept: float
    curvature_r_hat: float
    kappa_d: float


def main() -> int:
    ap = argparse.ArgumentParser(description="Conjecture E layer-2 flat-space baseline")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 4], help="spacetime dimensions to test")
    ap.add_argument("--n", type=int, default=192)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tau-quantile", type=float, default=0.75, help="fit window on the lower-tau quantile")
    ap.add_argument("--out", default="outputs_unified_functional/curvature_layer2_baseline.csv")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[Summary] = []
    for d in args.dims:
        if d < 2:
            raise ValueError("d must be >= 2")
        tau_parts: list[np.ndarray] = []
        k_parts: list[np.ndarray] = []
        for rep in range(args.reps):
            points = sprinkle_causal_diamond(args.n, d - 1, seed=args.seed + d * 1000 + rep)
            poset = poset_from_points(points)
            tau, k = causal_pair_data(points, poset)
            if tau.size == 0:
                continue
            tau_parts.append(tau)
            k_parts.append(k)
        if not tau_parts:
            raise RuntimeError(f"No causal pairs for d={d}, N={args.n}")
        tau = np.concatenate(tau_parts)
        k = np.concatenate(k_parts)

        y = k / max(args.n - 2, 1)
        x = np.power(tau, d)

        # Flat-space volume law fit in normalized units.
        flat_prefactor, flat_rmse, flat_r2 = fit_flat_volume_law(tau, y, d)

        # Curvature proxy: normalize by the flat diamond law, then fit residual slope.
        ratio = y / np.clip(x, 1e-12, None)
        fit_mask = tau <= np.quantile(tau, args.tau_quantile)
        if fit_mask.sum() < 10:
            fit_mask = np.ones_like(tau, dtype=bool)
        flat_prefactor, flat_rmse, flat_r2 = fit_flat_volume_law(tau[fit_mask], y[fit_mask], d)
        curvature_intercept, curvature_r_hat = fit_curvature_proxy(tau, ratio, d)

        rows.append(
            Summary(
                d=d,
                N=args.n,
                seed=args.seed + d * 1000,
                n_causal_pairs=int(tau.size),
                mean_tau=float(np.mean(tau)),
                mean_interval_size=float(np.mean(k)),
                flat_prefactor_fit=flat_prefactor,
                flat_rmse=flat_rmse,
                flat_r2=flat_r2,
                curvature_intercept=curvature_intercept,
                curvature_r_hat=curvature_r_hat,
                kappa_d=alexandrov_volume_constant(d),
            )
        )

    with out_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(Summary.__dataclass_fields__.keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row.__dict__)

    print(f"Saved: {out_path}")
    print("\nFlat-space baseline:")
    for row in rows:
        print(
            f"  d={row.d}: causal_pairs={row.n_causal_pairs}, "
            f"flat_prefactor_fit={row.flat_prefactor_fit:.6f}, "
            f"flat_r2={row.flat_r2:.4f}, curvature_R_hat={row.curvature_r_hat:+.4f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
