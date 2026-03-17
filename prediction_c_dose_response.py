#!/usr/bin/env python3
"""Prediction C — Layered k-family dose-response experiment.

Natural experiment: generate posets with controlled layer counts k = 2..8 at
fixed N, using exact entropy counting. Test whether increasing k (more
hierarchical layers) systematically lowers log_H.

This is the strongest possible dose-response design because:
- k is the CONTROLLED treatment variable
- N is held constant
- exact log_H eliminates estimation error
- many replications per (N, k) cell
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

from generators import generate_random_layered
from entropy_exact import count_linear_extensions_exact
from matched_residual_freedom_check import layer_index_by_minima

# ── Output dir ──────────────────────────────────────────────
OUT_DIR = Path("outputs_exploratory/prediction_c_dose_response")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Config ──────────────────────────────────────────────────
N_VALUES = [14, 16, 18]                 # exact counting feasible
K_VALUES = [2, 3, 4, 5, 6, 7, 8]       # layer counts (treatment levels)
REPS_PER_CELL = 60                      # replications per (N, k)
ADJACENT_P = 0.35
SKIP_P = 0.08
SEED_BASE = 20260317


def log_H(poset) -> float:
    return math.log(count_linear_extensions_exact(poset, prefer_c=False))


def actual_layer_count(poset) -> int:
    idx = layer_index_by_minima(poset.closure)
    return int(idx.max() + 1) if len(idx) > 0 else 0


def run():
    rng = np.random.default_rng(SEED_BASE)
    rows = []
    t0 = time.time()

    for N in N_VALUES:
        print(f"\n=== N = {N} ===")
        for k in K_VALUES:
            if k > N:
                continue
            for rep in range(REPS_PER_CELL):
                seed = SEED_BASE + N * 10000 + k * 1000 + rep
                poset = generate_random_layered(
                    n=N, n_layers=k,
                    imbalance_mode="uniform",
                    adjacent_p=ADJACENT_P,
                    skip_p=SKIP_P,
                    seed=seed,
                )
                lh = log_H(poset)
                actual_k = actual_layer_count(poset)
                rows.append({
                    "N": N,
                    "k_nominal": k,
                    "k_actual": actual_k,
                    "log_H": lh,
                    "log_H_per_element": lh / N,
                    "seed": seed,
                })
            elapsed = time.time() - t0
            print(f"  N={N}, k={k}: {REPS_PER_CELL} reps done ({elapsed:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "dose_response_raw.csv", index=False)
    print(f"\nTotal rows: {len(df)}")

    # ── Analysis ────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("ANALYSIS: k-family Dose-Response")
    print("=" * 60)

    # 1. Per-N summary
    summary = (
        df.groupby(["N", "k_nominal"])
        .agg(
            mean_log_H=("log_H", "mean"),
            std_log_H=("log_H", "std"),
            mean_k_actual=("k_actual", "mean"),
            mean_log_H_per_elem=("log_H_per_element", "mean"),
            n_reps=("seed", "count"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_DIR / "dose_response_summary.csv", index=False)
    print("\n", summary.to_string(index=False))

    # 2. Per-N linear regression: log_H ~ k_nominal
    print("\n--- Per-N Regression: log_H ~ k_nominal ---")
    reg_rows = []
    for N in N_VALUES:
        sub = df[df["N"] == N]
        slope, intercept, r, p_val, se = stats.linregress(
            sub["k_nominal"], sub["log_H"]
        )
        print(f"  N={N}: slope={slope:.4f}, r={r:.3f}, p={p_val:.2e}")
        reg_rows.append({
            "N": N,
            "slope": slope,
            "intercept": intercept,
            "r": r,
            "p_value": p_val,
            "se_slope": se,
            "n_obs": len(sub),
        })

    reg_df = pd.DataFrame(reg_rows)
    reg_df.to_csv(OUT_DIR / "dose_response_regression.csv", index=False)

    # 3. Jonckheere-Terpstra-like test: monotone decreasing trend?
    print("\n--- Monotone Trend Test ---")
    for N in N_VALUES:
        sub = df[df["N"] == N]
        groups = [sub[sub["k_nominal"] == k]["log_H"].values for k in K_VALUES if k <= N]
        # Page's L test approximation: count pairwise concordances
        concordant = 0
        total_pairs = 0
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                for xi in groups[i]:
                    for xj in groups[j]:
                        total_pairs += 1
                        if xi > xj:
                            concordant += 1
                        elif xi == xj:
                            concordant += 0.5
        prop = concordant / total_pairs if total_pairs > 0 else 0
        # Use Spearman on cell means as summary
        means = [g.mean() for g in groups]
        ks = list(range(len(means)))
        rho, p_sp = stats.spearmanr(ks, means)
        print(f"  N={N}: Spearman ρ(k, mean_log_H) = {rho:.3f}, p = {p_sp:.4f}")
        print(f"         Pairwise concordance (k↑ → log_H↓): {prop:.3f}")

    # 4. Plot: dose-response curves
    fig, axes = plt.subplots(1, len(N_VALUES), figsize=(5 * len(N_VALUES), 4), sharey=False)
    if len(N_VALUES) == 1:
        axes = [axes]
    for ax, N in zip(axes, N_VALUES):
        sub = summary[summary["N"] == N]
        ax.errorbar(
            sub["k_nominal"], sub["mean_log_H"],
            yerr=sub["std_log_H"],
            fmt="o-", capsize=4, color="C0"
        )
        ax.set_xlabel("k (number of layers)")
        ax.set_ylabel("log H (exact)")
        ax.set_title(f"N = {N}")
        ax.grid(True, alpha=0.3)
    plt.suptitle("Prediction C Dose-Response: k → log H", fontsize=13)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "dose_response_curves.png", dpi=150)
    print(f"\n  Saved: dose_response_curves.png")

    # 5. Normalized plot: log_H / log(N!) to show relative constraint
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    for N in N_VALUES:
        sub = summary[summary["N"] == N]
        log_n_fact = sum(math.log(i) for i in range(1, N + 1))
        ax2.plot(
            sub["k_nominal"],
            sub["mean_log_H"] / log_n_fact,
            "o-", label=f"N={N}"
        )
    ax2.set_xlabel("k (nominal layers)")
    ax2.set_ylabel("log H / log(N!)")
    ax2.set_title("Normalized Entropy vs Layer Count")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "dose_response_normalized.png", dpi=150)
    print(f"  Saved: dose_response_normalized.png")

    # ── Conclusion ──────────────────────────────────────────
    print("\n" + "=" * 60)
    all_slopes_negative = all(r["slope"] < 0 for _, r in reg_df.iterrows())
    all_significant = all(r["p_value"] < 0.05 for _, r in reg_df.iterrows())
    if all_slopes_negative and all_significant:
        print("  ✓ Dose-response CONFIRMED: log_H decreases with k at all N")
    elif all_slopes_negative:
        print("  ~ Direction confirmed but not all slopes significant")
    else:
        print("  ✗ Dose-response NOT confirmed")
    print("=" * 60)


if __name__ == "__main__":
    run()
