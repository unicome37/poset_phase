"""Conjecture E - Layer 2 scaling scan.

This script extends `curvature_layer2_baseline.py` into a finite-size scaling
diagnostic:

  - sweep several N values,
  - average over multiple seeds,
  - fit the flat Alexandrov-volume law on the lower-tau window,
  - track how the fitted prefactor and curvature proxy evolve with N.

The goal is to establish a stable baseline before attempting a genuinely
curved generator.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from curvature_layer2_baseline import (
    alexandrov_volume_constant,
    causal_pair_data,
    fit_curvature_proxy,
    fit_flat_volume_law,
    poset_from_points,
    sprinkle_causal_diamond,
)


@dataclass(frozen=True)
class ScaleRow:
    d: int
    N: int
    reps: int
    tau_quantile: float
    n_causal_pairs: int
    mean_tau: float
    median_tau: float
    flat_prefactor_fit: float
    flat_rmse: float
    flat_r2: float
    curvature_intercept: float
    curvature_r_hat: float
    kappa_d: float


def parse_int_list(values: list[str]) -> list[int]:
    return [int(v) for v in values]


def main() -> int:
    ap = argparse.ArgumentParser(description="Conjecture E layer-2 scaling scan")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 4])
    ap.add_argument("--n-values", nargs="*", type=int, default=[64, 96, 128, 160, 192])
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tau-quantile", type=float, default=0.75)
    ap.add_argument("--out", default="outputs_unified_functional/curvature_layer2_scaling_scan.csv")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[ScaleRow] = []
    for d in args.dims:
        if d < 2:
            raise ValueError("d must be >= 2")
        for n in args.n_values:
            tau_parts: list[np.ndarray] = []
            k_parts: list[np.ndarray] = []
            seed_base = args.seed + d * 1000 + n * 10
            for rep in range(args.reps):
                points = sprinkle_causal_diamond(n, d - 1, seed=seed_base + rep)
                poset = poset_from_points(points)
                tau, k = causal_pair_data(points, poset)
                if tau.size == 0:
                    continue
                tau_parts.append(tau)
                k_parts.append(k)

            if not tau_parts:
                raise RuntimeError(f"No causal pairs for d={d}, N={n}")

            tau = np.concatenate(tau_parts)
            k = np.concatenate(k_parts)
            y = k / max(n - 2, 1)
            x = np.power(tau, d)
            ratio = y / np.clip(x, 1e-12, None)
            fit_mask = tau <= np.quantile(tau, args.tau_quantile)
            if fit_mask.sum() < 10:
                fit_mask = np.ones_like(tau, dtype=bool)

            flat_prefactor, flat_rmse, flat_r2 = fit_flat_volume_law(tau[fit_mask], y[fit_mask], d)
            curvature_intercept, curvature_r_hat = fit_curvature_proxy(tau, ratio, d)

            rows.append(
                ScaleRow(
                    d=d,
                    N=n,
                    reps=args.reps,
                    tau_quantile=args.tau_quantile,
                    n_causal_pairs=int(tau.size),
                    mean_tau=float(np.mean(tau)),
                    median_tau=float(np.median(tau)),
                    flat_prefactor_fit=flat_prefactor,
                    flat_rmse=flat_rmse,
                    flat_r2=flat_r2,
                    curvature_intercept=curvature_intercept,
                    curvature_r_hat=curvature_r_hat,
                    kappa_d=alexandrov_volume_constant(d),
                )
            )

    with out_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(ScaleRow.__dataclass_fields__.keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row.__dict__)

    print(f"Saved: {out_path}")
    print("\nLayer-2 scaling summary:")
    for d in args.dims:
        print(f"\n  d={d}:")
        for row in [r for r in rows if r.d == d]:
            print(
                f"    N={row.N:3d} pairs={row.n_causal_pairs:6d} "
                f"pref={row.flat_prefactor_fit:.4f} r2={row.flat_r2:.3f} "
                f"R_hat={row.curvature_r_hat:+.3f}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

