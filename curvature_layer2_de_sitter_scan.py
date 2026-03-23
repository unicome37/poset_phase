"""Conjecture E - Layer 2 de Sitter-like scan.

This script is a more physical follow-up to the toy conformal scan.

We sample points in a flat-slicing FRW / de Sitter-like geometry with
scale factor

    a(t) = exp(H * t)

and causal relation

    x2 is in the causal future of x1 iff
      t2 > t1 and ||Δx|| <= ∫_{t1}^{t2} dt / a(t).

For H = 0 this reduces to the Minkowski causal diamond used in the
flat baseline. For H > 0 it produces a curvature-controlled causal
deformation while keeping the causal order explicit.

The output is still a finite-size diagnostic, not a proof of curvature
recovery. The aim is to see whether the layer-2 pipeline responds in a
systematic way to a more realistic curved generator.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from curvature_layer2_baseline import fit_curvature_proxy, fit_flat_volume_law
from generators import Poset, transitive_closure


def unit_ball_volume(m: int) -> float:
    return math.pi ** (m / 2.0) / math.gamma(m / 2.0 + 1.0)


def de_sitter_scale_factor(t: np.ndarray, hubble: float) -> np.ndarray:
    return np.exp(hubble * (t - 0.5))


def sprinkle_de_sitter_like_diamond(
    n: int,
    d_spatial: int,
    hubble: float,
    seed: int | None = None,
) -> np.ndarray:
    """Sample points in a finite FRW/de Sitter-like box with volume weighting."""
    rng = np.random.default_rng(seed)
    points = np.empty((0, 1 + d_spatial), dtype=float)
    a_max = float(np.exp(abs(hubble) * 0.5))
    while points.shape[0] < n:
        batch_size = max(500, 20 * n)
        t = rng.random(batch_size)
        x = rng.random((batch_size, d_spatial)) - 0.5
        a = de_sitter_scale_factor(t, hubble)
        if hubble >= 0:
            accept_prob = np.clip((a ** d_spatial) / (a_max ** d_spatial), 0.0, 1.0)
        else:
            accept_prob = np.clip((a ** d_spatial), 0.0, 1.0)
        accept = rng.random(batch_size) < accept_prob
        accepted = np.column_stack([t[accept], x[accept]])
        points = np.vstack([points, accepted])
    return points[:n]


def poset_from_de_sitter_points(points: np.ndarray, hubble: float) -> Poset:
    t = points[:, 0]
    spatial = points[:, 1:]
    dt = t[None, :] - t[:, None]

    if hubble == 0.0:
        horizon = np.clip(dt, 0.0, None)
    else:
        ti = t[:, None]
        tj = t[None, :]
        horizon = (np.exp(-hubble * ti) - np.exp(-hubble * tj)) / hubble
        horizon = np.clip(horizon, 0.0, None)

    spatial_d2 = np.sum((spatial[None, :, :] - spatial[:, None, :]) ** 2, axis=2)
    adj = (dt > 0.0) & (spatial_d2 <= horizon * horizon)
    return Poset(transitive_closure(adj))


def de_sitter_interval_proxy(points: np.ndarray, hubble: float, poset: Poset) -> np.ndarray:
    """Return a causal-distance proxy that reduces to dt at H=0."""
    t = points[:, 0]
    dt = t[None, :] - t[:, None]
    if hubble == 0.0:
        return np.clip(dt, 0.0, None)[poset.closure]
    ti = t[:, None]
    tj = t[None, :]
    chi = (np.exp(-hubble * ti) - np.exp(-hubble * tj)) / hubble
    return np.clip(chi, 0.0, None)[poset.closure]


@dataclass(frozen=True)
class DeSitterRow:
    d: int
    N: int
    reps: int
    hubble: float
    tau_quantile: float
    n_causal_pairs: int
    mean_tau: float
    median_tau: float
    flat_prefactor_fit: float
    flat_rmse: float
    flat_r2: float
    curvature_intercept: float
    curvature_r_hat: float


def main() -> int:
    ap = argparse.ArgumentParser(description="Conjecture E layer-2 de Sitter-like scan")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 4])
    ap.add_argument("--n", type=int, default=320)
    ap.add_argument("--reps", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tau-quantile", type=float, default=0.8)
    ap.add_argument("--hubbles", nargs="*", type=float, default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--out", default="outputs_unified_functional/curvature_layer2_de_sitter_scan.csv")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[DeSitterRow] = []
    for d in args.dims:
        if d < 2:
            raise ValueError("d must be >= 2")
        for hubble in args.hubbles:
            tau_parts: list[np.ndarray] = []
            k_parts: list[np.ndarray] = []
            for rep in range(args.reps):
                points = sprinkle_de_sitter_like_diamond(
                    args.n,
                    d - 1,
                    hubble=hubble,
                    seed=args.seed + d * 1000 + rep + int(hubble * 100),
                )
                poset = poset_from_de_sitter_points(points, hubble)
                tau = de_sitter_interval_proxy(points, hubble, poset)
                k = (poset.closure.astype(np.int32) @ poset.closure.astype(np.int32))[poset.closure].astype(float)
                if tau.size == 0:
                    continue
                tau_parts.append(tau)
                k_parts.append(k)

            if not tau_parts:
                raise RuntimeError(f"No causal pairs for d={d}, N={args.n}, H={hubble}")

            tau = np.concatenate(tau_parts)
            k = np.concatenate(k_parts)
            y = k / max(args.n - 2, 1)
            x = np.power(tau, d)
            ratio = y / np.clip(x, 1e-12, None)
            fit_mask = tau <= np.quantile(tau, args.tau_quantile)
            if fit_mask.sum() < 10:
                fit_mask = np.ones_like(tau, dtype=bool)
            flat_prefactor, flat_rmse, flat_r2 = fit_flat_volume_law(tau[fit_mask], y[fit_mask], d)
            curvature_intercept, curvature_r_hat = fit_curvature_proxy(tau, ratio, d)

            rows.append(
                DeSitterRow(
                    d=d,
                    N=args.n,
                    reps=args.reps,
                    hubble=hubble,
                    tau_quantile=args.tau_quantile,
                    n_causal_pairs=int(tau.size),
                    mean_tau=float(np.mean(tau)),
                    median_tau=float(np.median(tau)),
                    flat_prefactor_fit=flat_prefactor,
                    flat_rmse=flat_rmse,
                    flat_r2=flat_r2,
                    curvature_intercept=curvature_intercept,
                    curvature_r_hat=curvature_r_hat,
                )
            )

    with out_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(DeSitterRow.__dataclass_fields__.keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row.__dict__)

    print(f"Saved: {out_path}")
    print("\nde Sitter-like summary:")
    for d in args.dims:
        print(f"\n  d={d}:")
        for row in [r for r in rows if r.d == d]:
            print(
                f"    H={row.hubble:>4.2f} pairs={row.n_causal_pairs:6d} "
                f"pref={row.flat_prefactor_fit:.4f} r2={row.flat_r2:.3f} "
                f"R_hat={row.curvature_r_hat:+.3f}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

