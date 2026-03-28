"""FLRW and Schwarzschild background sprinklers.

FLRW_matter: a(t) = (1 + kappa*t)^(2/3), decelerating matter-dominated.
  kappa = 0 → flat Minkowski (Lor4D).
  kappa > 0 → expanding universe with Ricci curvature from matter.
  Comoving horizon: chi = (3/kappa)*[(1+kappa*t2)^(1/3) - (1+kappa*t1)^(1/3)].

Schwarzschild_weak: linearised Schwarzschild (weak-field approximation).
  ds^2 ≈ -(1+2Φ)dt^2 + (1-2Φ)|dx|^2
  Φ(r) = -phi0 * 0.5 / sqrt(r^2 + eps^2),  softened point-mass at spatial centre.
  phi0 = 0 → flat Minkowski.
  phi0 > 0 → static gravitational field, light-cone narrowing near centre.
"""
from __future__ import annotations

import numpy as np

from generators import Poset, transitive_closure


# ── FLRW matter-dominated ──────────────────────────────────────────


def sprinkle_flrw_matter_diamond(
    n: int,
    d_spatial: int = 3,
    kappa: float = 0.0,
    seed: int | None = None,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    points = np.empty((0, 1 + d_spatial), dtype=float)
    a_max = (1.0 + kappa) ** (2.0 / 3.0)  # max at t=1

    while points.shape[0] < n:
        batch = max(500, 20 * n)
        t = rng.random(batch)
        x = rng.random((batch, d_spatial)) - 0.5
        a = (1.0 + kappa * t) ** (2.0 / 3.0)
        accept_prob = np.clip((a / a_max) ** d_spatial, 0.0, 1.0)
        accept = rng.random(batch) < accept_prob
        accepted = np.column_stack([t[accept], x[accept]])
        points = np.vstack([points, accepted])

    return points[:n]


def poset_from_flrw_matter_points(pts: np.ndarray, kappa: float) -> Poset:
    t = pts[:, 0]
    sp = pts[:, 1:]
    dt = t[None, :] - t[:, None]

    if kappa == 0.0:
        horizon = np.clip(dt, 0.0, None)
    else:
        ti, tj = t[:, None], t[None, :]
        horizon = (3.0 / kappa) * (
            (1.0 + kappa * tj) ** (1.0 / 3.0)
            - (1.0 + kappa * ti) ** (1.0 / 3.0)
        )
        horizon = np.clip(horizon, 0.0, None)

    sp_d2 = np.sum((sp[None, :, :] - sp[:, None, :]) ** 2, axis=2)
    adj = (dt > 0.0) & (sp_d2 <= horizon * horizon)
    return Poset(transitive_closure(adj))


# ── Schwarzschild weak-field ───────────────────────────────────────

_SOFTENING2 = 0.01 ** 2  # avoid 1/r singularity


def sprinkle_schwarzschild_diamond(
    n: int,
    d_spatial: int = 3,
    phi0: float = 0.0,
    seed: int | None = None,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    t = rng.random(n)
    x = rng.random((n, d_spatial)) - 0.5

    if phi0 == 0.0:
        return np.column_stack([t, x])

    # Volume weight: sqrt(-g) ≈ (1-2Φ)^{d/2} (1+2Φ)^{1/2}
    # Clamp Φ to weak-field regime to avoid sqrt of negative
    r2 = np.sum(x ** 2, axis=1)
    phi = -phi0 * 0.5 / np.sqrt(r2 + _SOFTENING2)
    phi = np.clip(phi, -0.4, 0.0)  # keep 1+2Φ ≥ 0.2
    w = (1.0 - 2.0 * phi) ** (d_spatial / 2.0) * (1.0 + 2.0 * phi) ** 0.5
    w_max = w.max()

    # Rejection on volume weight
    mask = rng.random(n) < (w / w_max)
    accepted = np.column_stack([t[mask], x[mask]])

    # Top-up if we lost too many
    while accepted.shape[0] < n:
        extra = max(500, 5 * n)
        t2 = rng.random(extra)
        x2 = rng.random((extra, d_spatial)) - 0.5
        r2_ = np.sum(x2 ** 2, axis=1)
        phi_ = -phi0 * 0.5 / np.sqrt(r2_ + _SOFTENING2)
        phi_ = np.clip(phi_, -0.4, 0.0)
        w_ = (1.0 - 2.0 * phi_) ** (d_spatial / 2.0) * (1.0 + 2.0 * phi_) ** 0.5
        m2 = rng.random(extra) < (w_ / w_max)
        accepted = np.vstack([accepted, np.column_stack([t2[m2], x2[m2]])])

    return accepted[:n]


def poset_from_schwarzschild_points(pts: np.ndarray, phi0: float) -> Poset:
    t = pts[:, 0]
    sp = pts[:, 1:]
    dt = t[None, :] - t[:, None]

    if phi0 == 0.0:
        sp_d2 = np.sum((sp[None, :, :] - sp[:, None, :]) ** 2, axis=2)
        adj = (dt > 0.0) & (sp_d2 <= dt * dt)
        return Poset(transitive_closure(adj))

    # Midpoint potential (clamped to weak-field regime)
    mid_sp = 0.5 * (sp[None, :, :] + sp[:, None, :])
    r_mid2 = np.sum(mid_sp ** 2, axis=2)
    phi = -phi0 * 0.5 / np.sqrt(r_mid2 + _SOFTENING2)
    phi = np.clip(phi, -0.4, 0.0)

    sp_d2 = np.sum((sp[None, :, :] - sp[:, None, :]) ** 2, axis=2)
    # (1+2Φ)(Δt)^2 > (1-2Φ)|Δx|^2,  Δt > 0
    lhs = (1.0 + 2.0 * phi) * dt * dt
    rhs = (1.0 - 2.0 * phi) * sp_d2
    adj = (dt > 0.0) & (lhs > rhs)
    return Poset(transitive_closure(adj))
