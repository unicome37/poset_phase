#!/usr/bin/env python3
"""Prediction C — Edge Density Phase Diagram.

How does the slope d(log_H)/dk depend on edge density p?

We sweep p from 0.05 to 1.0, generating k-layered posets at each p,
and measure log_H vs k. This maps out the complete phase diagram:

- At p → 0: layers barely connected, log_H ≈ log(N!), slope ≈ 0
- At p = 1: complete layered, slope is maximally negative (theorem)
- In between: smooth transition

The phase diagram reveals:
1. Whether the Prediction C effect has a critical p below which it vanishes
2. How the sparse-regime slope compares to the analytical bound
3. Whether there's a non-monotone region (sign reversal)
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from generators import Poset, transitive_closure
from entropy_exact import count_linear_extensions_exact
from matched_residual_freedom_check import layer_index_by_minima

OUT_DIR = Path("outputs_exploratory/prediction_c_edge_density_phase")
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_VALUES = [14, 16]
K_VALUES = [2, 3, 4, 5, 6]
P_VALUES = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
REPS_PER_CELL = 40
SEED_BASE = 20260318


def generate_layered_at_p(n: int, k: int, p: float, seed: int) -> Poset:
    """Generate k-layered poset with uniform edge probability p between all layer pairs."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    m = n // k
    remainder = n % k
    sizes = [m + 1] * remainder + [m] * (k - remainder)
    layers = []
    offset = 0
    for s in sizes:
        layers.append(perm[offset:offset + s])
        offset += s

    adj = np.zeros((n, n), dtype=bool)
    for i in range(len(layers)):
        for j in range(i + 1, len(layers)):
            if p >= 1.0:
                adj[np.ix_(layers[i], layers[j])] = True
            else:
                mask = rng.random((len(layers[i]), len(layers[j]))) < p
                adj[np.ix_(layers[i], layers[j])] = mask
    return Poset(transitive_closure(adj))


def log_H(poset: Poset) -> float:
    return math.log(count_linear_extensions_exact(poset, prefer_c=False))


def actual_layer_count(poset: Poset) -> int:
    idx = layer_index_by_minima(poset.closure)
    return int(idx.max() + 1) if len(idx) > 0 else 0


def log_H_complete(N: int, k: int) -> float:
    if N % k != 0:
        m = N // k
        r = N % k
        sizes = [m + 1] * r + [m] * (k - r)
        return sum(math.lgamma(s + 1) for s in sizes)
    m = N // k
    return k * math.lgamma(m + 1)


def run():
    rows = []
    t0 = time.time()

    for N in N_VALUES:
        for p in P_VALUES:
            for k in K_VALUES:
                if k > N:
                    continue
                for rep in range(REPS_PER_CELL):
                    seed = SEED_BASE + N * 100000 + int(p * 100) * 10000 + k * 1000 + rep
                    poset = generate_layered_at_p(N, k, p, seed)
                    lh = log_H(poset)
                    ak = actual_layer_count(poset)
                    rows.append({
                        "N": N, "k_nominal": k, "p": p,
                        "k_actual": ak,
                        "log_H": lh,
                        "seed": seed,
                    })
            elapsed = time.time() - t0
            print(f"  N={N}, p={p:.2f}: done ({elapsed:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "edge_density_phase_raw.csv", index=False)
    print(f"\nTotal samples: {len(df)}")

    # ── Analysis ────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("PHASE DIAGRAM: d(log_H)/dk vs edge density p")
    print("=" * 70)

    # 1. Compute slope (d log_H / d k) at each (N, p) by regression
    slope_rows = []
    for N in N_VALUES:
        log_nfact = math.lgamma(N + 1)
        for p in P_VALUES:
            sub = df[(df["N"] == N) & (df["p"] == p)]
            if len(sub) < 10:
                continue
            # Regression on actual k
            slope_a, _, r_a, p_a, se_a = stats.linregress(sub["k_actual"], sub["log_H"])
            # Regression on nominal k
            slope_n, _, r_n, p_n, se_n = stats.linregress(sub["k_nominal"], sub["log_H"])
            slope_rows.append({
                "N": N, "p": p,
                "slope_actual_k": slope_a, "r_actual_k": r_a, "p_actual_k": p_a,
                "slope_nominal_k": slope_n, "r_nominal_k": r_n, "p_nominal_k": p_n,
            })

    slope_df = pd.DataFrame(slope_rows)
    slope_df.to_csv(OUT_DIR / "edge_density_slope.csv", index=False)

    print("\n--- Slope d(log_H)/dk_nominal vs p ---")
    for N in N_VALUES:
        print(f"\n  N={N}:")
        sub = slope_df[slope_df["N"] == N]
        for _, row in sub.iterrows():
            sign = "✓ negative" if row.slope_nominal_k < 0 else "✗ positive"
            print(f"    p={row.p:.2f}: slope={row.slope_nominal_k:.4f}, r={row.r_nominal_k:.3f} ({sign})")

    # 2. Analytical comparison at p=1
    print("\n--- Analytical vs Empirical at p=1.0 ---")
    for N in N_VALUES:
        print(f"\n  N={N}:")
        for k in K_VALUES:
            if k > N:
                continue
            expected = log_H_complete(N, k)
            sub = df[(df["N"] == N) & (df["p"] == 1.0) & (df["k_nominal"] == k)]
            if len(sub) > 0:
                empirical = sub["log_H"].mean()
                err = abs(empirical - expected)
                print(f"    k={k}: analytical={expected:.4f}, empirical={empirical:.4f}, err={err:.4f}")

    # 3. Plot phase diagram
    fig, axes = plt.subplots(1, len(N_VALUES), figsize=(6 * len(N_VALUES), 4.5))
    if len(N_VALUES) == 1:
        axes = [axes]

    for ax, N in zip(axes, N_VALUES):
        sub = slope_df[slope_df["N"] == N]
        ax.plot(sub["p"], sub["slope_nominal_k"], "o-", color="C0", label="slope d(log_H)/dk_nom")
        ax.axhline(0, color="gray", linestyle=":", alpha=0.5)
        ax.set_xlabel("Edge probability p")
        ax.set_ylabel("Slope d(log_H)/dk")
        ax.set_title(f"N = {N}")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.suptitle("Prediction C Phase Diagram: Slope vs Edge Density", fontsize=12)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "phase_diagram_slope_vs_p.png", dpi=150)
    print(f"\n  Saved: phase_diagram_slope_vs_p.png")

    # 4. Heatmap: mean log_H as function of (k, p)
    fig2, axes2 = plt.subplots(1, len(N_VALUES), figsize=(6 * len(N_VALUES), 4.5))
    if len(N_VALUES) == 1:
        axes2 = [axes2]

    for ax, N in zip(axes2, N_VALUES):
        pivot = df[df["N"] == N].groupby(["p", "k_nominal"])["log_H"].mean().reset_index()
        pivot_table = pivot.pivot(index="p", columns="k_nominal", values="log_H")
        im = ax.imshow(pivot_table.values, aspect="auto", origin="lower",
                       extent=[K_VALUES[0]-0.5, K_VALUES[-1]+0.5,
                               P_VALUES[0]-0.025, P_VALUES[-1]+0.025],
                       cmap="viridis")
        ax.set_xlabel("k (nominal layers)")
        ax.set_ylabel("Edge probability p")
        ax.set_title(f"N = {N}: mean log_H(k, p)")
        plt.colorbar(im, ax=ax, label="log H")

    plt.suptitle("Prediction C: Entropy Landscape", fontsize=12)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "entropy_landscape_heatmap.png", dpi=150)
    print(f"  Saved: entropy_landscape_heatmap.png")

    # 5. Key finding
    print("\n" + "=" * 70)
    # Check if slope is negative at all p values
    all_negative = slope_df["slope_nominal_k"].max() < 0
    if all_negative:
        print("  ✓ Slope d(log_H)/dk < 0 at ALL edge densities tested")
        print("    → Prediction C holds universally from p=0.05 to p=1.0")
        print("    → No phase transition or sign reversal detected")
    else:
        positive_p = slope_df[slope_df["slope_nominal_k"] >= 0]["p"].values
        print(f"  ⚠ Sign reversal detected at p = {positive_p}")

    # Magnitude trend
    print("\n  Magnitude trend:")
    for N in N_VALUES:
        sub = slope_df[slope_df["N"] == N]
        print(f"    N={N}: slope ranges from {sub['slope_nominal_k'].min():.3f} (p={sub.loc[sub['slope_nominal_k'].idxmin(), 'p']:.2f}) "
              f"to {sub['slope_nominal_k'].max():.3f} (p={sub.loc[sub['slope_nominal_k'].idxmax(), 'p']:.2f})")
    print("=" * 70)


if __name__ == "__main__":
    run()
