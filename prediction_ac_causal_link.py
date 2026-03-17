"""
Analysis: Causal structure diagnostics across dimensions (2D-5D)
Tests the hypothesis that higher-dimensional sprinklings have sparser
causal structure, linking Prediction C's hierarchy-entropy relationship
to Prediction A's dimensional dominance question.
"""
import numpy as np
import pandas as pd
from pathlib import Path

from generators import (
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_kr_like,
)
from observables import neutral_penalty
from observables_geo import geometric_components
from entropy_exact import log_linear_extensions_exact
from entropy_sis import estimate_log_linear_extensions_sis

OUT_DIR = Path("outputs_exploratory/prediction_ac_causal_link")
OUT_DIR.mkdir(parents=True, exist_ok=True)

FAMILIES = {
    "lorentzian_like_2d": generate_lorentzian_like_2d,
    "lorentzian_like_3d": generate_lorentzian_like_3d,
    "lorentzian_like_4d": generate_lorentzian_like_4d,
    "lorentzian_like_5d": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}
N_VALUES = [20, 28, 36, 44, 52]
SAMPLES = 4

rows = []
for n in N_VALUES:
    print(f"N={n}")
    for fname, gen in FAMILIES.items():
        for sid in range(SAMPLES):
            seed = 950000 + 1000 * n + sid
            poset = gen(n=n, seed=seed)
            adj = poset.closure

            # Basic causal structure diagnostics
            n_edges = int(np.sum(adj))
            max_possible = n * (n - 1) // 2
            edge_density = n_edges / max_possible if max_possible > 0 else 0

            # Layer structure (topological sort into layers)
            # Use the same method as in HII computation
            remaining = set(range(n))
            layers = []
            while remaining:
                # Minimal elements: no predecessors in remaining
                minimals = []
                for v in remaining:
                    has_pred = False
                    for u in remaining:
                        if u != v and adj[u, v]:
                            has_pred = True
                            break
                    if not has_pred:
                        minimals.append(v)
                if not minimals:
                    break
                layers.append(minimals)
                remaining -= set(minimals)

            n_layers = len(layers)
            layer_sizes = [len(l) for l in layers]
            max_antichain_approx = max(layer_sizes) if layer_sizes else 0

            # Longest chain (height)
            # Simple DP
            topo = []
            for layer in layers:
                topo.extend(layer)
            chain_len = np.zeros(n, dtype=int)
            for v in topo:
                for u in range(n):
                    if adj[u, v] and chain_len[u] + 1 > chain_len[v]:
                        chain_len[v] = chain_len[u] + 1
            height = int(chain_len.max()) + 1

            # Comparable fraction
            comp_pairs = n_edges
            total_pairs = max_possible
            comparable_fraction = comp_pairs / total_pairs if total_pairs > 0 else 0

            # Entropy (exact for small, SIS for large)
            if n <= 24 or (fname == "lorentzian_like_2d" and n <= 60):
                log_h = log_linear_extensions_exact(poset)
                method = "exact"
            else:
                log_h, _ = estimate_log_linear_extensions_sis(poset, n_runs=2048, seed=seed)
                method = "sis"

            rows.append({
                "n": n,
                "family": fname,
                "sample_id": sid,
                "n_edges": n_edges,
                "edge_density": edge_density,
                "n_layers": n_layers,
                "height": height,
                "max_antichain": max_antichain_approx,
                "comparable_fraction": comparable_fraction,
                "log_H": log_h,
                "method": method,
            })
        print(f"  {fname}: edges={rows[-1]['n_edges']}, layers={rows[-1]['n_layers']}, "
              f"height={rows[-1]['height']}, logH={rows[-1]['log_H']:.2f}")

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "causal_diagnostics.csv", index=False)

# Summary table
print("\n" + "="*80)
print("CAUSAL STRUCTURE × ENTROPY SUMMARY")
print("="*80)

summary = df.groupby(["n", "family"]).agg(
    mean_edges=("n_edges", "mean"),
    mean_density=("edge_density", "mean"),
    mean_layers=("n_layers", "mean"),
    mean_height=("height", "mean"),
    mean_max_antichain=("max_antichain", "mean"),
    mean_comp_frac=("comparable_fraction", "mean"),
    mean_logH=("log_H", "mean"),
).reset_index()

summary.to_csv(OUT_DIR / "causal_summary.csv", index=False)

# Print dimension comparison at each N
for n in N_VALUES:
    sub = summary[summary["n"] == n].sort_values("family")
    print(f"\n--- N = {n} ---")
    print(f"{'Family':25s} {'Edges':>7s} {'Density':>8s} {'Layers':>7s} {'Height':>7s} {'MaxAnti':>8s} {'CompFrac':>9s} {'logH':>8s}")
    for _, r in sub.iterrows():
        print(f"{r['family']:25s} {r['mean_edges']:7.1f} {r['mean_density']:8.4f} "
              f"{r['mean_layers']:7.1f} {r['mean_height']:7.1f} {r['mean_max_antichain']:8.1f} "
              f"{r['mean_comp_frac']:9.4f} {r['mean_logH']:8.2f}")

# Cross-dimension trend
print("\n" + "="*80)
print("DIMENSIONAL TREND: Edge density & entropy vs dimension")
print("="*80)
lor_families = ["lorentzian_like_2d", "lorentzian_like_3d", "lorentzian_like_4d", "lorentzian_like_5d"]
for n in N_VALUES:
    print(f"\nN={n}:")
    sub = summary[(summary["n"] == n) & (summary["family"].isin(lor_families))].copy()
    sub["dim"] = sub["family"].map({"lorentzian_like_2d": 2, "lorentzian_like_3d": 3,
                                     "lorentzian_like_4d": 4, "lorentzian_like_5d": 5})
    sub = sub.sort_values("dim")
    for _, r in sub.iterrows():
        dim = int(r["dim"])
        print(f"  {dim}D: density={r['mean_density']:.4f}, layers={r['mean_layers']:.1f}, "
              f"height={r['mean_height']:.1f}, logH={r['mean_logH']:.2f}")

# Correlation: density vs logH across all Lorentzian families
from scipy import stats
lor_data = df[df["family"].isin(lor_families)]
r_dens_ent, p_dens_ent = stats.pearsonr(lor_data["edge_density"], lor_data["log_H"])
r_layer_ent, p_layer_ent = stats.pearsonr(lor_data["n_layers"], lor_data["log_H"])
r_height_ent, p_height_ent = stats.pearsonr(lor_data["height"], lor_data["log_H"])

print(f"\n=== CROSS-DIMENSIONAL CORRELATIONS (Lorentzian families only) ===")
print(f"  edge_density vs logH:  r = {r_dens_ent:+.4f}, p = {p_dens_ent:.2e}")
print(f"  n_layers vs logH:     r = {r_layer_ent:+.4f}, p = {p_layer_ent:.2e}")
print(f"  height vs logH:       r = {r_height_ent:+.4f}, p = {p_height_ent:.2e}")

# Within fixed N
print(f"\n=== WITHIN-N CORRELATIONS (Lorentzian, density vs logH) ===")
for n in N_VALUES:
    sub = lor_data[lor_data["n"] == n]
    if len(sub) >= 4:
        r, p = stats.pearsonr(sub["edge_density"], sub["log_H"])
        print(f"  N={n}: r = {r:+.4f}, p = {p:.4f}")

print(f"\nSaved to {OUT_DIR}")
