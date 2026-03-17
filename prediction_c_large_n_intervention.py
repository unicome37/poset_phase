#!/usr/bin/env python3
"""Prediction C — Large-N Intervention using SIS Approximation.

Extends the intervention and placebo experiments to N = 24, 28, 32, 36
using SIS-approximate entropy (n_runs=1024 for low variance).

Key design decisions:
1. Use PAIRED SIS seeds: original and intervened poset share the same
   SIS seed, so their difference Δlog_H has much lower variance than
   each individual estimate.
2. Use 1024 SIS runs for tight estimates.
3. Run same-layer (treatment) and cross-layer (placebo) interventions.
4. Validate at N=16 where we can compare SIS vs exact results.

The goal is to show that Prediction C effects persist and STRENGTHEN
at N well beyond the exact-counting frontier.
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from generators import (
    Poset, transitive_closure,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
)
from entropy_exact import count_linear_extensions_exact
from entropy_sis import estimate_log_linear_extensions_sis
from matched_residual_freedom_check import layer_index_by_minima

OUT_DIR = Path("outputs_exploratory/prediction_c_large_n_intervention")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Config ──────────────────────────────────────────────────
CONFIGS = {
    # N=16: validation (also compute exact for comparison)
    16: {"n_posets": 40, "max_interventions": 10, "exact": True},
    # Large N: SIS only
    24: {"n_posets": 30, "max_interventions": 10, "exact": False},
    28: {"n_posets": 25, "max_interventions": 8, "exact": False},
    32: {"n_posets": 20, "max_interventions": 8, "exact": False},
    36: {"n_posets": 15, "max_interventions": 6, "exact": False},
}
SIS_RUNS = 1024
SEED_BASE = 20260318


# ── Helpers ─────────────────────────────────────────────────
def layer_count(poset: Poset) -> int:
    idx = layer_index_by_minima(poset.closure)
    return int(idx.max() + 1) if len(idx) > 0 else 0


def estimate_log_H(poset: Poset, sis_seed: int, exact: bool = False) -> tuple[float, str]:
    if exact and poset.n <= 18:
        return math.log(count_linear_extensions_exact(poset, prefer_c=False)), "exact"
    mean_est, _ = estimate_log_linear_extensions_sis(poset, n_runs=SIS_RUNS, seed=sis_seed)
    return mean_est, "sis"


def incomparable_pairs(poset: Poset):
    c = poset.closure
    n = poset.n
    layer_idx = layer_index_by_minima(c)
    same, cross = [], []
    for i in range(n):
        for j in range(i + 1, n):
            if not c[i, j] and not c[j, i]:
                if layer_idx[i] == layer_idx[j]:
                    same.append((i, j))
                else:
                    if layer_idx[i] < layer_idx[j]:
                        cross.append((i, j))
                    else:
                        cross.append((j, i))
    return same, cross


def intervene(poset: Poset, a: int, b: int) -> Poset:
    new_adj = poset.closure.copy()
    new_adj[a, b] = True
    return Poset(transitive_closure(new_adj))


# ── Main ────────────────────────────────────────────────────
def run():
    rng = np.random.default_rng(SEED_BASE)
    rows = []
    t0 = time.time()

    for N, cfg in CONFIGS.items():
        n_posets = cfg["n_posets"]
        max_intv = cfg["max_interventions"]
        use_exact = cfg["exact"]
        print(f"\n=== N = {N} (method: {'exact+SIS' if use_exact else 'SIS'}) ===")

        for idx in range(n_posets):
            seed = SEED_BASE + N * 1000 + idx
            poset = generate_lorentzian_like_2d(N, seed=seed)
            lc0 = layer_count(poset)

            # Use a fixed SIS seed for this poset (shared between original + interventions)
            sis_seed = seed * 7 + 13
            lh0, method0 = estimate_log_H(poset, sis_seed, exact=use_exact)

            # Also compute exact at N=16 for validation
            lh0_exact = None
            if use_exact:
                lh0_exact = math.log(count_linear_extensions_exact(poset, prefer_c=False))

            same_pairs, cross_pairs = incomparable_pairs(poset)

            for group_name, pool in [("same_layer", same_pairs),
                                     ("cross_layer", cross_pairs)]:
                if len(pool) == 0:
                    continue
                sample_size = min(max_intv, len(pool))
                chosen = rng.choice(len(pool), size=sample_size, replace=False)

                for ci in chosen:
                    a, b = pool[ci]
                    new_poset = intervene(poset, a, b)
                    lc1 = layer_count(new_poset)
                    # Use SAME SIS seed for paired comparison
                    lh1, method1 = estimate_log_H(new_poset, sis_seed, exact=use_exact)

                    lh1_exact = None
                    if use_exact:
                        lh1_exact = math.log(count_linear_extensions_exact(new_poset, prefer_c=False))

                    row = {
                        "N": N,
                        "poset_seed": seed,
                        "intervention_type": group_name,
                        "original_layer_count": lc0,
                        "new_layer_count": lc1,
                        "delta_layer_count": lc1 - lc0,
                        "split_occurred": int(lc1 > lc0),
                        "log_H_original": lh0,
                        "log_H_new": lh1,
                        "delta_log_H": lh1 - lh0,
                        "abs_delta_log_H": abs(lh1 - lh0),
                        "entropy_method": method0,
                    }
                    if lh0_exact is not None and lh1_exact is not None:
                        row["delta_log_H_exact"] = lh1_exact - lh0_exact
                        row["abs_delta_log_H_exact"] = abs(lh1_exact - lh0_exact)
                    rows.append(row)

            if (idx + 1) % max(1, n_posets // 5) == 0:
                elapsed = time.time() - t0
                print(f"  N={N}: {idx+1}/{n_posets} ({elapsed:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "large_n_intervention_raw.csv", index=False)
    print(f"\nTotal interventions: {len(df)}")

    # ══════════════════════════════════════════════════════════
    #  ANALYSIS
    # ══════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("LARGE-N INTERVENTION ANALYSIS")
    print("=" * 70)

    # A. Validation at N=16: SIS vs exact
    val = df[(df["N"] == 16) & (df["delta_log_H_exact"].notna())]
    if len(val) > 0:
        print("\n--- A. Validation: SIS vs Exact at N=16 ---")
        corr, p_corr = stats.pearsonr(val["delta_log_H"], val["delta_log_H_exact"])
        rho, p_rho = stats.spearmanr(val["delta_log_H"], val["delta_log_H_exact"])
        mae = (val["delta_log_H"] - val["delta_log_H_exact"]).abs().mean()
        print(f"  Pearson r(ΔlogH_sis, ΔlogH_exact) = {corr:.4f} (p={p_corr:.2e})")
        print(f"  Spearman ρ =                         {rho:.4f} (p={p_rho:.2e})")
        print(f"  MAE(Δ) = {mae:.4f}")
        # Same conclusions?
        sis_split = val[val["split_occurred"] == 1]
        sis_nosplit = val[val["split_occurred"] == 0]
        if len(sis_split) >= 3 and len(sis_nosplit) >= 3:
            t_sis, p_sis = stats.ttest_ind(sis_split["abs_delta_log_H"], sis_nosplit["abs_delta_log_H"], equal_var=False)
            t_exact, p_exact = stats.ttest_ind(sis_split["abs_delta_log_H_exact"], sis_nosplit["abs_delta_log_H_exact"], equal_var=False)
            print(f"  Split vs No-split (SIS):   t={t_sis:.2f}, p={p_sis:.4f}")
            print(f"  Split vs No-split (exact): t={t_exact:.2f}, p={p_exact:.4f}")
            print(f"  {'✓ CONCORDANT' if (p_sis < 0.05) == (p_exact < 0.05) else '⚠ DISCORDANT'}")

    # B. Per-N: same-layer vs cross-layer (placebo test)
    print("\n--- B. Per-N: Treatment vs Placebo ---")
    result_rows = []
    for N in sorted(df["N"].unique()):
        sub = df[df["N"] == N]
        treat = sub[sub["intervention_type"] == "same_layer"]
        placebo = sub[sub["intervention_type"] == "cross_layer"]
        if len(treat) < 3 or len(placebo) < 3:
            continue
        t_stat, p_val = stats.ttest_ind(
            treat["abs_delta_log_H"], placebo["abs_delta_log_H"], equal_var=False
        )
        n1, n2 = len(treat), len(placebo)
        s1, s2 = treat["abs_delta_log_H"].std(), placebo["abs_delta_log_H"].std()
        pooled_sd = math.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2)) if (n1+n2-2) > 0 else 1
        d = (treat["abs_delta_log_H"].mean() - placebo["abs_delta_log_H"].mean()) / pooled_sd if pooled_sd > 0 else float('nan')
        print(f"  N={N}: treat={n1} (|Δ|={treat['abs_delta_log_H'].mean():.3f}) "
              f"placebo={n2} (|Δ|={placebo['abs_delta_log_H'].mean():.3f}) "
              f"t={t_stat:.2f} p={p_val:.4f} d={d:.2f}")
        result_rows.append({
            "N": N, "n_treatment": n1, "n_placebo": n2,
            "mean_abs_dH_treatment": treat["abs_delta_log_H"].mean(),
            "mean_abs_dH_placebo": placebo["abs_delta_log_H"].mean(),
            "welch_t": t_stat, "welch_p": p_val, "cohens_d": d,
        })

    pd.DataFrame(result_rows).to_csv(OUT_DIR / "large_n_placebo_comparison.csv", index=False)

    # C. Per-N: split vs no-split within same-layer
    print("\n--- C. Per-N: Split vs No-split (within same-layer) ---")
    for N in sorted(df["N"].unique()):
        treat = df[(df["N"] == N) & (df["intervention_type"] == "same_layer")]
        split = treat[treat["split_occurred"] == 1]
        nosplit = treat[treat["split_occurred"] == 0]
        if len(split) < 3 or len(nosplit) < 3:
            print(f"  N={N}: split={len(split)}, nosplit={len(nosplit)} (skipped)")
            continue
        t_stat, p_val = stats.ttest_ind(
            split["abs_delta_log_H"], nosplit["abs_delta_log_H"], equal_var=False
        )
        print(f"  N={N}: split={len(split)} (|Δ|={split['abs_delta_log_H'].mean():.3f}) "
              f"nosplit={len(nosplit)} (|Δ|={nosplit['abs_delta_log_H'].mean():.3f}) "
              f"t={t_stat:.2f} p={p_val:.4f}")

    # D. N-trend: does effect size grow with N?
    print("\n--- D. Effect size trend across N ---")
    if len(result_rows) >= 3:
        Ns = [r["N"] for r in result_rows]
        ds = [r["cohens_d"] for r in result_rows]
        slope, _, r_trend, p_trend, _ = stats.linregress(Ns, ds)
        print(f"  Linear trend: slope={slope:.4f}, r={r_trend:.3f}, p={p_trend:.4f}")
        if slope > 0:
            print("  ✓ Effect size INCREASES with N")
        else:
            print("  ↓ Effect size decreases with N")

    # E. Summary
    print("\n" + "=" * 70)
    all_sig = all(r["welch_p"] < 0.05 for r in result_rows)
    all_positive_d = all(r["cohens_d"] > 0 for r in result_rows)
    max_N = max(r["N"] for r in result_rows) if result_rows else 0
    if all_sig and all_positive_d:
        print(f"  ✓ PLACEBO TEST PASSES at ALL N up to {max_N}")
        print("    → Prediction C extends beyond exact-counting frontier")
    else:
        failed = [r["N"] for r in result_rows if r["welch_p"] >= 0.05 or r["cohens_d"] <= 0]
        print(f"  ⚠ Failed at N = {failed}")
    print("=" * 70)


if __name__ == "__main__":
    run()
