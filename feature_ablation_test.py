"""
Feature Ablation + Minimal Complete Basis Test
===============================================
Test whether (d_eff, C₁/C₀, width) is the minimal complete basis for
Lor4D identification by systematic feature ablation:

  1. Full triple: (d, c, w)  — baseline
  2. Drop d_eff: (c, w)
  3. Drop C₁/C₀: (d, w)
  4. Drop width: (d, c)
  5. Single features: d only, c only, w only

Also test adding candidate 4th features:
  6. + longest chain (height ratio h/N)
  7. + total relations R (order fraction)

For each configuration: rank Lor4D among 17 families at each N.
"""
from __future__ import annotations

import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from bd_action import count_intervals_fast
from generators import (
    Poset,
    generate_kr_like, generate_kr_2layer, generate_kr_4layer,
    generate_lorentzian_like_2d, generate_lorentzian_like_3d,
    generate_lorentzian_like_4d, generate_lorentzian_like_5d,
    generate_transitive_percolation, generate_interval_order,
    generate_absolute_layered, generate_multi_layer_random,
    generate_random_layered_k4_uniform, generate_random_layered_k6_uniform,
    generate_random_layered_k8_uniform, generate_random_layered_k6_tapered,
    generate_random_layered_k6_middle_heavy, generate_random_layered_k6_longjump,
)
from unified_functional import compute_xi_dim


def longest_chain_length(poset: Poset) -> int:
    n = poset.n
    if n <= 1:
        return n
    c = poset.closure
    dp = [1] * n
    for i in range(n):
        for j in range(i + 1, n):
            if c[i, j]:
                dp[j] = max(dp[j], dp[i] + 1)
    return max(dp)


def max_antichain_width(poset: Poset) -> int:
    n = poset.n
    if n <= 1:
        return n
    c = poset.closure
    remaining = set(range(n))
    max_w = 0
    while remaining:
        minimals = []
        for i in remaining:
            is_min = True
            for j in remaining:
                if j != i and c[j, i] and not c[i, j]:
                    is_min = False
                    break
            if is_min:
                minimals.append(i)
        if not minimals:
            break
        max_w = max(max_w, len(minimals))
        for m in minimals:
            remaining.discard(m)
    return max_w


def compute_features(poset: Poset, N: int) -> dict:
    counts = count_intervals_fast(poset, k_max=5)
    C0 = counts.get(0)
    C1 = counts.get(1)
    total_rel = counts.total_relations
    R = 1.0 - float(C0) / float(total_rel) if total_rel > 0 else 0.0
    c1_c0 = C1 / max(1, C0)
    xi_val, d_eff = compute_xi_dim(poset)
    aw = max_antichain_width(poset)
    width_ratio = aw / max(1, N)
    lc = longest_chain_length(poset)
    height_ratio = lc / max(1, N)
    order_frac = R
    return {
        "d_eff": d_eff,
        "c1_c0": c1_c0,
        "width_ratio": width_ratio,
        "height_ratio": height_ratio,
        "order_frac": order_frac,
    }


FAMILIES = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
    "KR_2layer": generate_kr_2layer,
    "KR_4layer": generate_kr_4layer,
    "AbsLayer": generate_absolute_layered,
    "MLR": generate_multi_layer_random,
    "RLk4": generate_random_layered_k4_uniform,
    "RLk6": generate_random_layered_k6_uniform,
    "RLk8": generate_random_layered_k8_uniform,
    "RLk6_tap": generate_random_layered_k6_tapered,
    "RLk6_mid": generate_random_layered_k6_middle_heavy,
    "RLk6_lj": generate_random_layered_k6_longjump,
    "TransPerc": generate_transitive_percolation,
    "IntOrder": generate_interval_order,
}

CATEGORY = {}
for f in ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]:
    CATEGORY[f] = "Lorentzian"
for f in ["KR_like", "KR_2layer", "KR_4layer"]:
    CATEGORY[f] = "KR-family"
for f in ["AbsLayer", "MLR", "RLk4", "RLk6", "RLk8",
          "RLk6_tap", "RLk6_mid", "RLk6_lj"]:
    CATEGORY[f] = "Layered"
for f in ["TransPerc", "IntOrder"]:
    CATEGORY[f] = "Other"


def evaluate_config(all_feats, N_VALUES, feature_keys, targets_fn):
    """Evaluate Mahalanobis distance for a given feature subset."""
    results = {}
    for N in N_VALUES:
        # Lor4D stats
        lor_rows = [f for f in all_feats if f["family"] == "Lor4D" and f["N"] == N]
        if len(lor_rows) < 3:
            continue
        data = np.array([[r[k] for k in feature_keys] for r in lor_rows])
        mu = np.mean(data, axis=0)
        cov = np.cov(data.T) if data.shape[1] > 1 else np.array([[np.var(data, ddof=1)]])
        if data.shape[1] == 1:
            cov = cov.reshape(1, 1)
        try:
            inv_cov = np.linalg.inv(cov + 1e-10 * np.eye(len(feature_keys)))
        except np.linalg.LinAlgError:
            inv_cov = np.eye(len(feature_keys))

        targets = targets_fn(N, mu)

        # Compute Mahalanobis for each family
        means = {}
        for fam in FAMILIES:
            rows = [f for f in all_feats if f["family"] == fam and f["N"] == N]
            if not rows:
                continue
            vals = []
            for r in rows:
                x = np.array([r[k] for k in feature_keys])
                delta = x - targets
                vals.append(float(delta @ inv_cov @ delta))
            means[fam] = np.mean(vals)

        ranked = sorted(means, key=means.get)
        r4d = ranked.index("Lor4D") + 1 if "Lor4D" in ranked else 99
        lor4d_f = means.get("Lor4D", 0)
        nonlor = [f for f in ranked if CATEGORY[f] != "Lorentzian"]
        if nonlor:
            ru = nonlor[0]
            margin = means[ru] - lor4d_f
        else:
            ru = "N/A"
            margin = 0
        results[N] = {"rank": r4d, "runner_up": ru, "margin": margin}
    return results


def main():
    N_VALUES = [16, 20, 28, 36, 48, 64, 96, 128]
    REPS = 25
    SEED_BASE = 42

    print("=" * 80)
    print("Feature Ablation + Minimal Complete Basis Test")
    print("=" * 80)

    # Phase 1: Generate data
    all_feats = []
    total = len(FAMILIES) * len(N_VALUES) * REPS
    done = 0
    t0 = time.time()
    for fam_name, gen_fn in FAMILIES.items():
        for N in N_VALUES:
            for rep in range(REPS):
                seed = (SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep) % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    feat = compute_features(poset, N)
                    feat["family"] = fam_name
                    feat["N"] = N
                    all_feats.append(feat)
                except Exception:
                    pass
                done += 1
                if done % 500 == 0:
                    elapsed = time.time() - t0
                    rate = done / max(elapsed, 0.01)
                    eta = (total - done) / max(rate, 0.01)
                    print(f"  [{done}/{total}] {elapsed:.0f}s, ETA {eta:.0f}s", flush=True)
    elapsed = time.time() - t0
    print(f"\n  Generated {len(all_feats)} samples in {elapsed:.1f}s\n")

    # Lor4D centroids for target construction
    lor4d_centroids = {}
    for N in N_VALUES:
        rows = [f for f in all_feats if f["family"] == "Lor4D" and f["N"] == N]
        if rows:
            lor4d_centroids[N] = {
                k: np.mean([r[k] for r in rows])
                for k in ["d_eff", "c1_c0", "width_ratio", "height_ratio", "order_frac"]
            }

    report = []
    report.append("# Feature Ablation + Minimal Complete Basis Test\n")

    # Define configurations
    configs = {
        "Full (d,c,w)": {
            "keys": ["d_eff", "c1_c0", "width_ratio"],
            "targets": lambda N, mu: np.array([4.0, mu[1], mu[2]]),
        },
        "Drop d: (c,w)": {
            "keys": ["c1_c0", "width_ratio"],
            "targets": lambda N, mu: mu,
        },
        "Drop c: (d,w)": {
            "keys": ["d_eff", "width_ratio"],
            "targets": lambda N, mu: np.array([4.0, mu[1]]),
        },
        "Drop w: (d,c)": {
            "keys": ["d_eff", "c1_c0"],
            "targets": lambda N, mu: np.array([4.0, mu[1]]),
        },
        "d only": {
            "keys": ["d_eff"],
            "targets": lambda N, mu: np.array([4.0]),
        },
        "c only": {
            "keys": ["c1_c0"],
            "targets": lambda N, mu: mu,
        },
        "w only": {
            "keys": ["width_ratio"],
            "targets": lambda N, mu: mu,
        },
        "+height: (d,c,w,h)": {
            "keys": ["d_eff", "c1_c0", "width_ratio", "height_ratio"],
            "targets": lambda N, mu: np.array([4.0, mu[1], mu[2], mu[3]]),
        },
        "+order_frac: (d,c,w,R)": {
            "keys": ["d_eff", "c1_c0", "width_ratio", "order_frac"],
            "targets": lambda N, mu: np.array([4.0, mu[1], mu[2], mu[3]]),
        },
    }

    # Evaluate each configuration
    all_results = {}
    for name, cfg in configs.items():
        print(f"  Evaluating: {name}", flush=True)
        res = evaluate_config(all_feats, N_VALUES, cfg["keys"], cfg["targets"])
        all_results[name] = res

    # Report: rank table
    report.append("\n## 1. Lor4D Rank by Configuration\n")
    header = "| Config | " + " | ".join(f"N={N}" for N in N_VALUES) + " | All#1? |"
    sep = "|--------|" + "|".join(":---:" for _ in N_VALUES) + "|:-----:|"
    report.append(header)
    report.append(sep)

    for name in configs:
        res = all_results[name]
        ranks = [res.get(N, {}).get("rank", "?") for N in N_VALUES]
        all_one = all(r == 1 for r in ranks if isinstance(r, int))
        rank_strs = [f"#{r}" for r in ranks]
        report.append(f"| {name} | " + " | ".join(rank_strs) + f" | {'✅' if all_one else '❌'} |")

    # Report: margin table
    report.append("\n\n## 2. Margin (to nearest non-Lor) by Configuration\n")
    header2 = "| Config | " + " | ".join(f"N={N}" for N in N_VALUES) + " | Mean |"
    sep2 = "|--------|" + "|".join(":---:" for _ in N_VALUES) + "|:----:|"
    report.append(header2)
    report.append(sep2)

    for name in configs:
        res = all_results[name]
        margins = [res.get(N, {}).get("margin", 0) for N in N_VALUES]
        mean_m = np.mean(margins) if margins else 0
        m_strs = [f"{m:.1f}" for m in margins]
        report.append(f"| {name} | " + " | ".join(m_strs) + f" | {mean_m:.1f} |")

    # Report: What each dropped feature costs
    report.append("\n\n## 3. Ablation Impact Analysis\n")
    full_res = all_results["Full (d,c,w)"]
    for drop_name, dep_name in [
        ("Drop d: (c,w)", "d_eff"),
        ("Drop c: (d,w)", "C₁/C₀"),
        ("Drop w: (d,c)", "width"),
    ]:
        drop_res = all_results[drop_name]
        fails = []
        margin_losses = []
        for N in N_VALUES:
            full_r = full_res.get(N, {}).get("rank", 99)
            drop_r = drop_res.get(N, {}).get("rank", 99)
            if drop_r > 1:
                fails.append(N)
            full_m = full_res.get(N, {}).get("margin", 0)
            drop_m = drop_res.get(N, {}).get("margin", 0)
            margin_losses.append(full_m - drop_m)

        report.append(f"### Dropping {dep_name}\n")
        if fails:
            report.append(f"- **Failures** (Lor4D not #1): N = {fails}")
        else:
            report.append(f"- **No failures**: Lor4D still #1 at all N")
        avg_loss = np.mean(margin_losses)
        report.append(f"- **Mean margin loss**: {avg_loss:.1f}")
        report.append(f"- **Impact**: {'CRITICAL' if fails else 'Moderate' if avg_loss > 10 else 'Low'}\n")

    # Report: 4th feature benefit
    report.append("\n## 4. Fourth Feature Benefit\n")
    for add_name, feat_name in [
        ("+height: (d,c,w,h)", "height_ratio"),
        ("+order_frac: (d,c,w,R)", "order_frac"),
    ]:
        add_res = all_results[add_name]
        margin_gains = []
        for N in N_VALUES:
            full_m = full_res.get(N, {}).get("margin", 0)
            add_m = add_res.get(N, {}).get("margin", 0)
            margin_gains.append(add_m - full_m)
        avg_gain = np.mean(margin_gains)
        report.append(f"- **{feat_name}**: mean margin gain = {avg_gain:.1f}")

    report.append(f"\n**Conclusion**: If gains are small and no new failures are resolved, "
                  f"the triple is the minimal complete basis.")

    # Summary
    report.append(f"\n\n## 5. Verdict\n")
    # Check if any single drop causes failure
    critical_features = []
    for drop_name, dep_name in [
        ("Drop d: (c,w)", "d_eff"),
        ("Drop c: (d,w)", "C₁/C₀"),
        ("Drop w: (d,c)", "width"),
    ]:
        res = all_results[drop_name]
        if any(res.get(N, {}).get("rank", 99) > 1 for N in N_VALUES):
            critical_features.append(dep_name)

    if len(critical_features) == 3:
        report.append("✅ **All three features are individually necessary**: dropping any one "
                      "causes Lor4D to lose #1 at some N.")
        report.append("\n→ **(d_eff, C₁/C₀, width) is the minimal complete basis** for Lor4D identification.")
    elif len(critical_features) > 0:
        report.append(f"🟡 **Critical features**: {critical_features}")
        report.append(f"Other features are helpful but not strictly necessary for rank #1.")
    else:
        report.append("❌ No single feature is individually necessary — the basis may be redundant.")

    # Write report
    out_path = Path("outputs_carlip") / "feature_ablation_test.md"
    out_path.write_text("\n".join(report), encoding="utf-8")
    print(f"\nReport saved to {out_path}")

    # Console summary
    print(f"\n{'='*60}")
    print(f"ABLATION SUMMARY:")
    for name in configs:
        res = all_results[name]
        ranks = [res.get(N, {}).get("rank", 99) for N in N_VALUES]
        all1 = all(r == 1 for r in ranks)
        print(f"  {name:30s}  {'ALL#1' if all1 else 'FAILS at ' + str([N for N,r in zip(N_VALUES,ranks) if r>1])}")
    print(f"\nCritical (dropping causes failure): {critical_features}")


if __name__ == "__main__":
    main()
