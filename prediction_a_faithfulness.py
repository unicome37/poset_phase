"""Prediction A — Faithfulness / Manifold-Likeness Filter

Instead of adding S_BD as an energy term (requires huge α), test whether
2D sprinklings systematically FAIL manifold-likeness diagnostics that 4D passes.

Three diagnostics:
  1. Dimension Stability: variance of multi-method d_eff estimates
  2. Interval Regularity: how well the interval distribution matches a
     smooth-manifold prediction (Myrheim-Meyer type)
  3. Link Dominance: fraction of causal relations that are links (C_0/total)

If 2D fails one or more of these as a HARD FILTER while 4D passes,
then Prediction A's "why not 2D" is answered by a constraint, not an energy term.

Usage:
  python prediction_a_faithfulness.py
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from generators import Poset
from observables import comparable_fraction, layer_profile
from observables_geo import (
    dimension_proxy_views,
    dimension_consistency_penalty,
    estimate_dimension_proxy_from_order_fraction,
)
from prediction_a_bd_bridge import (
    regenerate_poset,
    count_intervals_fast,
    compute_f5,
    CALIBRATED_WEIGHTS,
)


# ---------------------------------------------------------------------------
# Faithfulness diagnostics
# ---------------------------------------------------------------------------

def diagnostic_dim_stability(poset: Poset) -> dict:
    """Dimension Stability: variance across 3 dimension proxies.

    A manifold-faithful causal set should give consistent d_eff estimates
    from different methods. High variance = poor manifold-likeness.
    """
    views = dimension_proxy_views(poset)
    d_vals = [views["d_order"], views["d_chain"], views["d_width"]]
    mean_d = np.mean(d_vals)
    var_d = np.var(d_vals)
    spread = max(d_vals) - min(d_vals)

    # Also get local-vs-global consistency from existing function
    total_pen, global_d, mean_local_d, var_local, n_intervals = \
        dimension_consistency_penalty(poset)

    return {
        "d_order": views["d_order"],
        "d_chain": views["d_chain"],
        "d_width": views["d_width"],
        "d_mean": float(mean_d),
        "d_var": float(var_d),
        "d_spread": float(spread),
        "d_consistency_penalty": float(total_pen),
        "d_global": float(global_d),
        "d_local_mean": float(mean_local_d),
        "d_local_var": float(var_local),
        "n_local_intervals": int(n_intervals),
    }


def diagnostic_interval_regularity(poset: Poset) -> dict:
    """Interval Regularity: how well the interval distribution fits
    a smooth manifold.

    For d-dim Minkowski sprinkling, the expected fraction of k-element
    intervals follows a specific distribution. Deviations indicate
    non-manifold-like structure.

    Key metric: the "interval depth excess" — how much deeper are the
    intervals than expected for a manifold of the measured dimension?
    """
    counts = count_intervals_fast(poset)
    total = sum(counts.values())
    if total == 0:
        return {"total_relations": 0, "link_fraction": 0.0,
                "mean_interval": 0.0, "interval_excess": 0.0,
                "interval_entropy": 0.0, "max_interval": 0}

    C_0 = counts.get(0, 0)
    link_frac = C_0 / total

    weighted_sum = sum(k * v for k, v in counts.items())
    mean_interval = weighted_sum / total

    max_interval = max(counts.keys())

    # Interval distribution entropy (higher = more spread out = richer)
    probs = np.array([counts.get(k, 0) / total for k in range(max_interval + 1)])
    probs = probs[probs > 0]
    interval_entropy = -float(np.sum(probs * np.log(probs)))

    # Interval depth excess: compare mean interval to what's expected
    # for the measured dimension.
    # For d-dim Minkowski, link fraction scales as ~1 - c/N^(2/d)
    # At our small N, a rough expectation for d-dim:
    #   d=2: link_frac ~ 0.5-0.7 for N=16-36
    #   d=4: link_frac ~ 0.8-0.95
    # The fact that our 2D gives link_frac=0.35 is ANOMALOUS —
    # it means 2D sprinklings are "too connected" at these N values.
    #
    # Expected link fraction for faithful d-dim at density ρ:
    #   f_link(d) ~ 1 / (1 + ρ * V_d * τ^d / d!)
    # We don't compute this precisely; instead use the deviation as a signal.

    # Simple excess: how much deeper than "pure links" is the structure?
    interval_excess = mean_interval  # 0 = all links, higher = deeper

    return {
        "total_relations": int(total),
        "link_fraction": float(link_frac),
        "mean_interval": float(mean_interval),
        "interval_excess": float(interval_excess),
        "interval_entropy": float(interval_entropy),
        "max_interval": int(max_interval),
    }


def diagnostic_link_dominance(poset: Poset) -> dict:
    """Link Dominance: fraction of relations that are Hasse links.

    In a faithful d-dim causal set, most relations at moderate N should
    be links (direct covers). A low link fraction indicates the causal
    structure is "too connected" — either the dimension is too low for
    the density, or the structure isn't manifold-like.

    This is closely related to the Benincasa-Dowker action's physical
    content: the BD action weights links differently from longer paths.
    """
    c = poset.closure
    n = poset.n

    # Total causal relations
    total_relations = int(c.sum())
    if total_relations == 0:
        return {"total_relations": 0, "n_links": 0, "link_fraction": 1.0,
                "non_link_fraction": 0.0, "relations_per_element": 0.0}

    # Links = relations with no intermediate element
    # A link i→j means: c[i,j]=True AND there is no k with c[i,k] AND c[k,j]
    # This equals: pairs where interval_size = 0
    interval_matrix = c.astype(np.int32) @ c.astype(np.int32)

    n_links = 0
    for i in range(n):
        for j in range(n):
            if c[i, j] and interval_matrix[i, j] == 0:
                n_links += 1

    link_frac = n_links / total_relations
    relations_per_element = total_relations / n

    return {
        "total_relations": total_relations,
        "n_links": n_links,
        "link_fraction": float(link_frac),
        "non_link_fraction": 1.0 - float(link_frac),
        "relations_per_element": float(relations_per_element),
    }


# ---------------------------------------------------------------------------
# Composite faithfulness score
# ---------------------------------------------------------------------------

def faithfulness_score(poset: Poset) -> dict:
    """Compute all three faithfulness diagnostics and a composite score."""
    dim = diagnostic_dim_stability(poset)
    interval = diagnostic_interval_regularity(poset)
    link = diagnostic_link_dominance(poset)

    # Composite: combine three failure modes
    # Each diagnostic produces a "failure score" (higher = less manifold-like)

    # 1. Dimension instability: d_var + d_consistency_penalty
    dim_fail = dim["d_var"] + dim["d_consistency_penalty"]

    # 2. Interval depth excess: mean_interval (0 = pure links = good)
    interval_fail = interval["interval_excess"]

    # 3. Non-link excess: 1 - link_fraction (0 = pure links = good)
    link_fail = link["non_link_fraction"]

    # Relations per element (high = dense causal structure)
    density = link["relations_per_element"]

    return {
        **{f"dim_{k}": v for k, v in dim.items()},
        **{f"int_{k}": v for k, v in interval.items()},
        **{f"link_{k}": v for k, v in link.items()},
        "fail_dim": float(dim_fail),
        "fail_interval": float(interval_fail),
        "fail_link": float(link_fail),
        "fail_density": float(density),
    }


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def main() -> None:
    raw_path = "outputs_unified_functional/raw_features.csv"
    rows = list(csv.DictReader(open(raw_path)))
    lor_families = {"Lor2D", "Lor3D", "Lor4D", "Lor5D"}
    rows = [r for r in rows if r["family"] in lor_families]
    families_order = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]

    print("=" * 80)
    print("PREDICTION A — FAITHFULNESS / MANIFOLD-LIKENESS DIAGNOSTICS")
    print("=" * 80)

    # Compute diagnostics for all posets
    all_data = []
    for i, r in enumerate(rows):
        fam, n, rep = r["family"], int(r["N"]), int(r["rep"])
        poset = regenerate_poset(fam, n, rep, 42)
        f5 = compute_f5(r, CALIBRATED_WEIGHTS)
        scores = faithfulness_score(poset)
        all_data.append({"family": fam, "N": n, "rep": rep, "f5": f5, **scores})
        if (i + 1) % 16 == 0:
            print(f"  [{i+1}/{len(rows)}] computed")

    print(f"  Total: {len(all_data)} posets\n")

    # --- Report 1: Family statistics ---
    print("=" * 80)
    print("FAMILY STATISTICS")
    print("=" * 80)

    metrics = [
        ("fail_dim", "Dim Instability"),
        ("fail_interval", "Interval Depth"),
        ("fail_link", "Non-Link Frac"),
        ("fail_density", "Relations/Element"),
        ("dim_d_mean", "Mean d_eff"),
        ("dim_d_var", "d_eff Variance"),
        ("dim_d_spread", "d_eff Spread"),
        ("dim_d_consistency_penalty", "Dim Consistency Pen"),
        ("int_link_fraction", "Link Fraction"),
        ("int_mean_interval", "Mean Interval"),
        ("int_interval_entropy", "Interval Entropy"),
        ("int_max_interval", "Max Interval"),
        ("link_relations_per_element", "Relations/Element"),
    ]

    for metric_key, metric_label in metrics:
        print(f"\n  {metric_label} ({metric_key}):")
        print(f"    {'Family':<8s} {'Mean':>10s} {'Std':>10s} {'Min':>10s} {'Max':>10s}")
        for fam in families_order:
            vals = [d[metric_key] for d in all_data if d["family"] == fam]
            if vals:
                print(f"    {fam:<8s} {np.mean(vals):>10.4f} {np.std(vals):>10.4f} "
                      f"{np.min(vals):>10.4f} {np.max(vals):>10.4f}")

    # --- Report 2: Filter tests ---
    print("\n" + "=" * 80)
    print("FILTER TESTS: Can faithfulness filters eliminate 2D while preserving 4D?")
    print("=" * 80)

    # For each metric, find thresholds that separate families
    filter_metrics = [
        ("fail_interval", "Interval Depth", ">"),      # reject if too deep
        ("fail_link", "Non-Link Fraction", ">"),        # reject if too many non-links
        ("fail_density", "Relations/Element", ">"),     # reject if too dense
        ("int_link_fraction", "Link Fraction", "<"),    # reject if too few links
        ("int_interval_entropy", "Interval Entropy", ">"),  # reject if too spread
        ("dim_d_var", "d_eff Variance", ">"),           # reject if dimension unstable
    ]

    for metric_key, metric_label, direction in filter_metrics:
        print(f"\n--- Filter: {metric_label} ({metric_key}) ---")
        print(f"    Direction: reject if {direction} threshold")

        # Get per-family distributions
        fam_vals = {}
        for fam in families_order:
            fam_vals[fam] = [d[metric_key] for d in all_data if d["family"] == fam]

        # Try percentile thresholds
        all_vals = [d[metric_key] for d in all_data]
        percentiles = [50, 60, 70, 75, 80, 85, 90, 95]

        print(f"    {'Pctl':>5s} {'Thresh':>8s} | {'2D pass':>8s} {'3D pass':>8s} "
              f"{'4D pass':>8s} {'5D pass':>8s} | Note")
        print(f"    {'-'*5} {'-'*8} | {'-'*8} {'-'*8} {'-'*8} {'-'*8} | {'-'*30}")

        for pctl in percentiles:
            threshold = np.percentile(all_vals, pctl)
            pass_rates = {}
            for fam in families_order:
                if direction == ">":
                    passed = sum(1 for v in fam_vals[fam] if v <= threshold)
                else:
                    passed = sum(1 for v in fam_vals[fam] if v >= threshold)
                pass_rates[fam] = passed / len(fam_vals[fam]) * 100

            note = ""
            if pass_rates["Lor2D"] < 30 and pass_rates["Lor4D"] > 70:
                note = "*** 2D eliminated, 4D preserved ***"
            elif pass_rates["Lor2D"] < 50 and pass_rates["Lor4D"] > 80:
                note = "* partial separation *"

            print(f"    {pctl:5d} {threshold:8.4f} | {pass_rates['Lor2D']:7.0f}% "
                  f"{pass_rates['Lor3D']:7.0f}% {pass_rates['Lor4D']:7.0f}% "
                  f"{pass_rates['Lor5D']:7.0f}% | {note}")

    # --- Report 3: Combined filter ---
    print("\n" + "=" * 80)
    print("COMBINED FILTER: intersection of multiple criteria")
    print("=" * 80)

    # Find best single filters and combine them
    # Strategy: use link_fraction > threshold AND interval_depth < threshold
    print("\n  Combined: link_fraction >= T_link AND mean_interval <= T_interval")
    print()

    for link_thresh in [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80]:
        for int_thresh in [1.0, 0.8, 0.5, 0.3]:
            pass_rates = {}
            for fam in families_order:
                fam_data = [d for d in all_data if d["family"] == fam]
                passed = sum(1 for d in fam_data
                             if d["int_link_fraction"] >= link_thresh
                             and d["int_mean_interval"] <= int_thresh)
                pass_rates[fam] = passed / len(fam_data) * 100

            if pass_rates["Lor2D"] <= 10 and pass_rates["Lor4D"] >= 80:
                print(f"  *** link>={link_thresh:.2f} AND int<={int_thresh:.1f}: "
                      f"2D={pass_rates['Lor2D']:.0f}% 3D={pass_rates['Lor3D']:.0f}% "
                      f"4D={pass_rates['Lor4D']:.0f}% 5D={pass_rates['Lor5D']:.0f}%")

    # Also try: density-based filter
    print("\n  Density filter: relations_per_element <= threshold")
    for dens_thresh in np.arange(2.0, 12.0, 1.0):
        pass_rates = {}
        for fam in families_order:
            fam_data = [d for d in all_data if d["family"] == fam]
            passed = sum(1 for d in fam_data if d["fail_density"] <= dens_thresh)
            pass_rates[fam] = passed / len(fam_data) * 100

        if pass_rates["Lor2D"] <= 20 and pass_rates["Lor4D"] >= 80:
            print(f"  *** dens<={dens_thresh:.1f}: "
                  f"2D={pass_rates['Lor2D']:.0f}% 3D={pass_rates['Lor3D']:.0f}% "
                  f"4D={pass_rates['Lor4D']:.0f}% 5D={pass_rates['Lor5D']:.0f}%")

    # --- Report 4: Per-N filter stability ---
    print("\n" + "=" * 80)
    print("PER-N FILTER STABILITY (link_fraction >= 0.60)")
    print("=" * 80)

    n_values = sorted(set(d["N"] for d in all_data))
    for n_val in n_values:
        n_data = [d for d in all_data if d["N"] == n_val]
        print(f"\n  N = {n_val}:")
        for fam in families_order:
            fam_data = [d for d in n_data if d["family"] == fam]
            link_fracs = [d["int_link_fraction"] for d in fam_data]
            mean_ints = [d["int_mean_interval"] for d in fam_data]
            pass_link60 = sum(1 for v in link_fracs if v >= 0.60)
            pass_link50 = sum(1 for v in link_fracs if v >= 0.50)
            total = len(fam_data)
            print(f"    {fam}: link_frac={np.mean(link_fracs):.3f}+/-{np.std(link_fracs):.3f} "
                  f"mean_int={np.mean(mean_ints):.3f} "
                  f"pass(>=0.50)={pass_link50}/{total} pass(>=0.60)={pass_link60}/{total}")

    # --- Report 5: F5 ranking AFTER filter ---
    print("\n" + "=" * 80)
    print("F5 RANKING AFTER FAITHFULNESS FILTER")
    print("=" * 80)

    for link_thresh in [0.50, 0.55, 0.60]:
        print(f"\n  Filter: link_fraction >= {link_thresh}")
        filtered = [d for d in all_data if d["int_link_fraction"] >= link_thresh]

        for fam in families_order:
            fam_data = [d for d in filtered if d["family"] == fam]
            if fam_data:
                f5_mean = np.mean([d["f5"] for d in fam_data])
                print(f"    {fam}: {len(fam_data)} surviving, mean F5 = {f5_mean:.2f}")
            else:
                print(f"    {fam}: 0 surviving (ELIMINATED)")

        # What's the ranking among survivors?
        surviving_fams = [fam for fam in families_order
                          if any(d["family"] == fam for d in filtered)]
        if surviving_fams:
            fam_means = {}
            for fam in surviving_fams:
                fam_data = [d for d in filtered if d["family"] == fam]
                fam_means[fam] = np.mean([d["f5"] for d in fam_data])
            rank = sorted(surviving_fams, key=lambda f: fam_means[f])
            print(f"    Ranking among survivors: {' < '.join(rank)}")
            if "Lor4D" in rank and rank[0] != "Lor4D":
                print(f"    4D is not minimum (position {rank.index('Lor4D')+1}/{len(rank)})")
            elif "Lor4D" in rank:
                print(f"    *** 4D IS MINIMUM ***")


if __name__ == "__main__":
    main()
