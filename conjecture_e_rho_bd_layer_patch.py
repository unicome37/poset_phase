"""Conjecture E — ρ_BD local density via layer block patches.

Replaces the Alexandrov interval patch scheme (additivity ratio = 0.12)
with a layer-based partition:

  1. Assign each node to its layer (topological sort level)
  2. Group consecutive layers into blocks of size `block_size` (default 2-3)
  3. For each block, extract the induced sub-poset
  4. Compute bd_ratio on each block
  5. Test additivity: weighted mean of local bd_ratios ≈ global bd_ratio

The layer block scheme guarantees:
  - Full coverage (every node belongs to exactly one block)
  - Non-trivial sub-posets (each block has O(N/L * block_size) nodes)
  - Preservation of local causal structure within each block

Outputs:
  - outputs_unified_functional/conjecture_e_rho_bd_layer_patch.csv
  - outputs_unified_functional/conjecture_e_rho_bd_layer_patch.md
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from bd_action import (
    IntervalCounts,
    bd_ratio_metric,
    bdg_action_d2_corrected,
    count_intervals_fast,
)
from generators import (
    Poset,
    generate_kr_like,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)


# ── Layer assignment (node → layer index) ────────────────────────────────

def layer_assignment(poset: Poset) -> np.ndarray:
    """Return array of shape (n,) where layer_assignment[i] = layer index of node i.

    Layer 0 = minimal elements, layer 1 = next level, etc.
    Uses the same algorithm as observables.layer_profile but returns per-node labels.
    """
    c = poset.closure.copy()
    n = poset.n
    indeg = c.sum(axis=0).astype(int)
    remaining = set(range(n))
    labels = np.full(n, -1, dtype=int)
    layer_idx = 0

    while remaining:
        mins = [i for i in remaining if indeg[i] == 0]
        if not mins:
            raise ValueError("Cycle detected — not a valid poset")
        for u in mins:
            labels[u] = layer_idx
            remaining.remove(u)
            for v in np.where(c[u])[0]:
                indeg[v] -= 1
        layer_idx += 1

    return labels


def layer_nodes_list(poset: Poset) -> list[list[int]]:
    """Return list of lists: layer_nodes[k] = [node indices in layer k]."""
    labels = layer_assignment(poset)
    n_layers = int(labels.max()) + 1
    result = [[] for _ in range(n_layers)]
    for i, lbl in enumerate(labels):
        result[lbl].append(i)
    return result


# ── Layer block patch ────────────────────────────────────────────────────

@dataclass(frozen=True)
class PatchResult:
    block_id: int
    layers_start: int
    layers_end: int
    n_nodes: int
    n_relations: int
    bd_ratio_local: float
    bdg_d2c_local: float
    C0: int
    C1: int


def compute_layer_block_patches(
    poset: Poset,
    block_size: int = 2,
) -> dict:
    """Partition poset into consecutive layer blocks and compute local BD observables.

    Args:
        poset: the input poset
        block_size: number of consecutive layers per block (2 or 3)

    Returns dict with:
        - patches: list of PatchResult
        - global_bd_ratio: global bd_ratio
        - global_bdg_d2c: global bdg_d2_corrected_norm
        - weighted_mean_bd: node-weighted mean of local bd_ratios
        - additivity_ratio: weighted_mean_bd / global_bd_ratio
        - coverage: fraction of nodes covered by non-degenerate patches
        - n_blocks: number of blocks
        - n_layers: total layers
    """
    layers = layer_nodes_list(poset)
    n_layers = len(layers)
    n = poset.n
    closure = poset.closure

    # Create blocks of consecutive layers
    blocks: list[tuple[int, int, list[int]]] = []  # (start, end, node_indices)
    start = 0
    while start < n_layers:
        end = min(start + block_size, n_layers)
        nodes = []
        for k in range(start, end):
            nodes.extend(layers[k])
        blocks.append((start, end - 1, nodes))
        start = end

    # Compute local BD on each block
    patches = []
    for bid, (lstart, lend, nodes) in enumerate(blocks):
        if len(nodes) < 2:
            patches.append(PatchResult(
                block_id=bid, layers_start=lstart, layers_end=lend,
                n_nodes=len(nodes), n_relations=0,
                bd_ratio_local=0.0, bdg_d2c_local=0.0, C0=0, C1=0,
            ))
            continue

        node_arr = np.array(nodes, dtype=int)
        sub_closure = closure[np.ix_(node_arr, node_arr)]
        sub_poset = Poset(closure=sub_closure)
        sub_counts = count_intervals_fast(sub_poset)

        patches.append(PatchResult(
            block_id=bid,
            layers_start=lstart,
            layers_end=lend,
            n_nodes=len(nodes),
            n_relations=sub_counts.total_relations,
            bd_ratio_local=bd_ratio_metric(sub_counts),
            bdg_d2c_local=bdg_action_d2_corrected(sub_counts, len(nodes), normalized=True),
            C0=sub_counts.get(0),
            C1=sub_counts.get(1),
        ))

    # Global reference
    global_counts = count_intervals_fast(poset)
    global_br = bd_ratio_metric(global_counts)
    global_bdg = bdg_action_d2_corrected(global_counts, n, normalized=True)

    # Node-weighted mean
    total_nodes_in_patches = sum(p.n_nodes for p in patches if p.n_nodes >= 2)
    if total_nodes_in_patches > 0 and global_br > 0:
        weighted_sum = sum(
            p.bd_ratio_local * p.n_nodes
            for p in patches if p.n_nodes >= 2
        )
        weighted_mean = weighted_sum / total_nodes_in_patches
        additivity = weighted_mean / global_br
    else:
        weighted_mean = float("nan")
        additivity = float("nan")

    # Unweighted mean (for comparison)
    valid_patches = [p for p in patches if p.n_nodes >= 2 and p.n_relations > 0]
    if valid_patches:
        unweighted_mean = float(np.mean([p.bd_ratio_local for p in valid_patches]))
    else:
        unweighted_mean = float("nan")

    coverage = total_nodes_in_patches / n if n > 0 else 0.0

    return {
        "patches": patches,
        "global_bd_ratio": global_br,
        "global_bdg_d2c": global_bdg,
        "weighted_mean_bd": weighted_mean,
        "unweighted_mean_bd": unweighted_mean,
        "additivity_ratio": additivity,
        "coverage": coverage,
        "n_blocks": len(blocks),
        "n_layers": n_layers,
    }


# ── Helpers ──────────────────────────────────────────────────────────────

FAMILY_GENS = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}


def seed_for(family: str, n: int, rep: int, base_seed: int) -> int:
    return int(base_seed + rep * 1000 + n * 100)


@dataclass(frozen=True)
class SummaryRow:
    family: str
    N: int
    rep: int
    block_size: int
    n_layers: int
    n_blocks: int
    global_bd_ratio: float
    weighted_mean_bd: float
    unweighted_mean_bd: float
    additivity_ratio: float
    coverage: float
    patch_std: float
    patch_cv: float
    min_patch_nodes: int
    max_patch_nodes: int


# ── Main ─────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(description="Conjecture E: ρ_BD layer block patch locality test")
    ap.add_argument("--raw-features", default="outputs_unified_functional/raw_features.csv")
    ap.add_argument("--out-csv", default="outputs_unified_functional/conjecture_e_rho_bd_layer_patch.csv")
    ap.add_argument("--out-md", default="outputs_unified_functional/conjecture_e_rho_bd_layer_patch.md")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--block-sizes", nargs="*", type=int, default=[2, 3, 4],
                     help="Layer block sizes to test")
    args = ap.parse_args()

    raw_path = Path(args.raw_features)
    if not raw_path.exists():
        raise FileNotFoundError(str(raw_path))

    raw_rows = list(csv.DictReader(raw_path.open("r", encoding="utf-8")))
    families = {"Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"}
    raw_rows = [r for r in raw_rows if r.get("family") in families]
    print(f"Loaded {len(raw_rows)} rows from {raw_path}")

    all_summary: list[SummaryRow] = []
    all_patches: list[dict] = []

    for bs in args.block_sizes:
        print(f"\n=== Block size = {bs} ===")
        for i, r in enumerate(raw_rows, start=1):
            fam = r["family"]
            n = int(r["N"])
            rep = int(r["rep"])
            s = seed_for(fam, n, rep, args.seed)
            poset = FAMILY_GENS[fam](n, seed=s)

            result = compute_layer_block_patches(poset, block_size=bs)

            valid = [p for p in result["patches"] if p.n_nodes >= 2 and p.n_relations > 0]
            if valid:
                ratios = np.array([p.bd_ratio_local for p in valid])
                pstd = float(np.std(ratios))
                pmean = float(np.mean(ratios))
                pcv = pstd / pmean if pmean > 0 else float("nan")
            else:
                pstd = float("nan")
                pcv = float("nan")

            all_summary.append(SummaryRow(
                family=fam, N=n, rep=rep, block_size=bs,
                n_layers=result["n_layers"], n_blocks=result["n_blocks"],
                global_bd_ratio=result["global_bd_ratio"],
                weighted_mean_bd=result["weighted_mean_bd"],
                unweighted_mean_bd=result["unweighted_mean_bd"],
                additivity_ratio=result["additivity_ratio"],
                coverage=result["coverage"],
                patch_std=pstd, patch_cv=pcv,
                min_patch_nodes=min(p.n_nodes for p in result["patches"]),
                max_patch_nodes=max(p.n_nodes for p in result["patches"]),
            ))

            # Save individual patch data for detailed analysis
            for p in result["patches"]:
                all_patches.append({
                    "family": fam, "N": n, "rep": rep, "block_size": bs,
                    **asdict(p),
                })

            if i % 40 == 0:
                sr = all_summary[-1]
                print(f"  [{i}/{len(raw_rows)}] {fam} N={n} bs={bs}: "
                      f"add={sr.additivity_ratio:.3f} blocks={sr.n_blocks} "
                      f"cov={sr.coverage:.2f}")

    # ── Save CSV ──
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(all_summary[0]).keys()))
        w.writeheader()
        for sr in all_summary:
            w.writerow(asdict(sr))
    print(f"\nSaved: {out_csv}")

    # Save patch details CSV
    patch_csv = out_csv.with_name("conjecture_e_rho_bd_layer_patch_details.csv")
    if all_patches:
        with patch_csv.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(all_patches[0].keys()))
            w.writeheader()
            w.writerows(all_patches)
        print(f"Saved: {patch_csv}")

    # ── Analysis ──
    print(f"\n{'='*70}")
    print("LAYER BLOCK PATCH ρ_BD ADDITIVITY ANALYSIS")
    print(f"{'='*70}")

    for bs in args.block_sizes:
        sub = [s for s in all_summary if s.block_size == bs]
        valid = [s for s in sub if not math.isnan(s.additivity_ratio)]
        if not valid:
            print(f"\nBlock size {bs}: no valid results")
            continue

        add_arr = np.array([s.additivity_ratio for s in valid])
        print(f"\n--- Block size = {bs} ---")
        print(f"Valid samples: {len(valid)}/{len(sub)}")
        print(f"Additivity ratio: mean={np.mean(add_arr):.4f}, median={np.median(add_arr):.4f}, std={np.std(add_arr):.4f}")

        # Per-family breakdown
        print(f"\n{'Family':>7s} | {'add_ratio':>9s} | {'std':>6s} | {'n':>3s} | {'blocks':>6s} | {'coverage':>8s}")
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
            fsub = [s for s in valid if s.family == fam]
            if not fsub:
                continue
            fa = np.array([s.additivity_ratio for s in fsub])
            fb = np.mean([s.n_blocks for s in fsub])
            fc = np.mean([s.coverage for s in fsub])
            print(f"{fam:>7s} | {np.mean(fa):9.4f} | {np.std(fa):6.4f} | {len(fsub):3d} | {fb:6.1f} | {fc:8.3f}")

        # Per-N breakdown
        ns = sorted({s.N for s in valid})
        print(f"\n{'N':>4s} | {'add_ratio':>9s} | {'std':>6s} | {'n':>3s}")
        for n in ns:
            nsub = [s for s in valid if s.N == n]
            na = np.array([s.additivity_ratio for s in nsub])
            print(f"{n:>4d} | {np.mean(na):9.4f} | {np.std(na):6.4f} | {len(nsub):3d}")

    # ── Comparison with Alexandrov baseline ──
    print(f"\n{'='*70}")
    print("COMPARISON: Layer block (best) vs Alexandrov interval (baseline)")
    print(f"{'='*70}")
    # Pick the block size with highest mean additivity
    best_bs = None
    best_mean = -1.0
    for bs in args.block_sizes:
        sub = [s for s in all_summary if s.block_size == bs and not math.isnan(s.additivity_ratio)]
        if sub:
            m = float(np.mean([s.additivity_ratio for s in sub]))
            if m > best_mean:
                best_mean = m
                best_bs = bs
    print(f"Best block size: {best_bs} (mean additivity = {best_mean:.4f})")
    print(f"Alexandrov interval baseline: mean additivity = 0.117")
    if best_mean > 0:
        improvement = best_mean / 0.117
        print(f"Improvement factor: {improvement:.1f}x")

    # ── Save markdown report ──
    md = Path(args.out_md)
    lines = []
    lines.append("# Conjecture E: ρ_BD Layer Block Patch Locality Test")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append("Partition poset into consecutive layer blocks of size `k` (=2,3,4 layers),")
    lines.append("compute bd_ratio on each block's induced sub-poset, test additivity:")
    lines.append("")
    lines.append("  additivity_ratio = (node-weighted mean of local bd_ratios) / global bd_ratio")
    lines.append("")
    lines.append("A ratio near 1.0 means bd_ratio decomposes into local densities ρ_BD(x).")
    lines.append("")

    for bs in args.block_sizes:
        sub = [s for s in all_summary if s.block_size == bs]
        valid = [s for s in sub if not math.isnan(s.additivity_ratio)]
        if not valid:
            continue

        add_arr = np.array([s.additivity_ratio for s in valid])
        lines.append(f"## Block size = {bs}")
        lines.append("")
        lines.append(f"Valid: {len(valid)}/{len(sub)}, mean={np.mean(add_arr):.4f}, median={np.median(add_arr):.4f}, std={np.std(add_arr):.4f}")
        lines.append("")
        lines.append("| family | add_ratio | std | n | blocks | coverage |")
        lines.append("|---|---:|---:|---:|---:|---:|")
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
            fsub = [s for s in valid if s.family == fam]
            if not fsub:
                continue
            fa = np.array([s.additivity_ratio for s in fsub])
            fb = np.mean([s.n_blocks for s in fsub])
            fc = np.mean([s.coverage for s in fsub])
            lines.append(f"| {fam} | {np.mean(fa):.4f} | {np.std(fa):.4f} | {len(fsub)} | {fb:.1f} | {fc:.3f} |")
        lines.append("")

        ns = sorted({s.N for s in valid})
        lines.append("| N | add_ratio | std | n |")
        lines.append("|---:|---:|---:|---:|")
        for n in ns:
            nsub = [s for s in valid if s.N == n]
            na = np.array([s.additivity_ratio for s in nsub])
            lines.append(f"| {n} | {np.mean(na):.4f} | {np.std(na):.4f} | {len(nsub)} |")
        lines.append("")

    lines.append("## Comparison with Alexandrov interval baseline")
    lines.append("")
    lines.append(f"| method | mean add_ratio |")
    lines.append(f"|---|---:|")
    lines.append(f"| Alexandrov interval | 0.117 |")
    if best_bs is not None:
        lines.append(f"| **Layer block (bs={best_bs})** | **{best_mean:.4f}** |")
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("If layer blocks achieve additivity_ratio >> 0.12 (the Alexandrov baseline),")
    lines.append("this validates that bd_ratio can be decomposed into a local density ρ_BD(x)")
    lines.append("defined on layer patches — a necessary step for the continuum limit.")
    lines.append("")
    lines.append("Perfect additivity (ratio=1) is not expected because boundary effects")
    lines.append("between blocks create extra/missing pairs. The key test is whether the")
    lines.append("ratio is stable across N (finite-size scaling) and across families.")

    md.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"\nSaved: {md}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
