"""Conjecture E — ρ_BD vs R(p) curvature cross-check.

Dual strategy:

  Strategy A — Cross-poset global correlation:
    For many poset instances, compute global bd_ratio and global bdg_d4.
    Correlate across instances within each family and across families.
    This tests: do posets with richer interval structure also have
    different discrete curvature?

  Strategy B — Within-poset sliding-window (Lor2D only):
    Lor2D has ~9 layers at N=48, enough for within-poset layer analysis.
    Use a sliding window (center ± 1 layer) to compute local bd_ratio
    and local bdg_d4, then correlate across windows within each poset.
    High-d families (Lor4D, KR_like) have only ~3 layers where
    boundary layers have bd_ratio=0, making within-poset correlation
    impossible.

Outputs:
  - outputs_unified_functional/conjecture_e_rho_bd_vs_curvature.csv
  - outputs_unified_functional/conjecture_e_rho_bd_vs_curvature.md
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from bd_action import (
    IntervalCounts,
    bd_ratio_metric,
    bdg_action_d2_corrected,
    bdg_action_d4_standard,
    count_intervals_fast,
)
from generators import (
    Poset,
    generate_kr_like,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_transitive_percolation,
)


# ── Layer assignment (reuse from rho_bd_layer_patch) ─────────────────────

def _layer_assignment(poset: Poset) -> np.ndarray:
    """Per-node layer index (0 = minimal elements)."""
    c = poset.closure
    n = poset.n
    indeg = c.sum(axis=0).astype(int)
    remaining = set(range(n))
    labels = np.full(n, -1, dtype=int)
    layer_idx = 0
    while remaining:
        mins = [i for i in remaining if indeg[i] == 0]
        if not mins:
            break
        for u in mins:
            labels[u] = layer_idx
            remaining.remove(u)
            for v in np.where(c[u])[0]:
                indeg[v] -= 1
        layer_idx += 1
    return labels


FAMILY_GENS = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
    "TransPerc": lambda n, seed=None: generate_transitive_percolation(n, p=0.08, seed=seed),
}


# ── Strategy A: Cross-poset global correlation ───────────────────────

@dataclass(frozen=True)
class GlobalRow:
    """Global observables for one poset instance."""
    family: str
    N: int
    rep: int
    seed: int
    n_layers: int
    bd_ratio: float
    bdg_d2c: float
    bdg_d4: float
    link_frac: float
    total_relations: int
    C0: int
    C1: int
    C2: int
    C3: int


def run_strategy_a(
    families: list[str],
    n_values: list[int],
    reps: int,
    base_seed: int,
) -> list[GlobalRow]:
    """Compute global BD + curvature for each poset instance."""
    results: list[GlobalRow] = []
    total = len(families) * len(n_values) * reps
    count = 0

    for fam in families:
        gen = FAMILY_GENS[fam]
        for n0 in n_values:
            for rep in range(reps):
                count += 1
                s = base_seed + rep * 1000 + n0 * 100
                poset = gen(n0, seed=s)
                labels = _layer_assignment(poset)
                n_layers = int(labels.max()) + 1
                counts = count_intervals_fast(poset, k_max=3)
                tr = counts.total_relations

                results.append(GlobalRow(
                    family=fam, N=n0, rep=rep, seed=s,
                    n_layers=n_layers,
                    bd_ratio=bd_ratio_metric(counts),
                    bdg_d2c=bdg_action_d2_corrected(counts, n0, normalized=True),
                    bdg_d4=bdg_action_d4_standard(counts, n0, normalized=True),
                    link_frac=float(counts.get(0)) / tr if tr > 0 else 0.0,
                    total_relations=tr,
                    C0=counts.get(0),
                    C1=counts.get(1),
                    C2=counts.get(2),
                    C3=counts.get(3),
                ))

                if count % 100 == 0 or count == total:
                    print(f"  [A: {count}/{total}] {fam} N={n0} rep={rep}")

    return results


# ── Strategy B: Within-poset sliding window (deep-layer families) ──────

@dataclass(frozen=True)
class WindowRow:
    """Observables for a sliding-window layer neighbourhood."""
    family: str
    N: int
    rep: int
    center_layer: int
    n_layers_total: int
    n_nodes: int
    total_relations: int
    bd_ratio: float
    bdg_d2c: float
    bdg_d4: float
    link_frac: float


def run_strategy_b(
    families: list[str],
    n_values: list[int],
    reps: int,
    base_seed: int,
    window: int = 1,
) -> list[WindowRow]:
    """Sliding-window layer analysis for families with enough layers."""
    results: list[WindowRow] = []
    total = len(families) * len(n_values) * reps
    count = 0

    for fam in families:
        gen = FAMILY_GENS[fam]
        for n0 in n_values:
            for rep in range(reps):
                count += 1
                s = base_seed + rep * 1000 + n0 * 100
                poset = gen(n0, seed=s)
                c = poset.closure
                labels = _layer_assignment(poset)
                n_layers = int(labels.max()) + 1

                for center in range(n_layers):
                    lo = max(0, center - window)
                    hi = min(n_layers - 1, center + window)
                    nodes = np.where((labels >= lo) & (labels <= hi))[0]
                    if len(nodes) < 3:
                        continue
                    sub_c = c[np.ix_(nodes, nodes)]
                    sub_p = Poset(closure=sub_c)
                    sub_counts = count_intervals_fast(sub_p, k_max=3)
                    tr = sub_counts.total_relations
                    if tr < 1:
                        continue
                    m = len(nodes)
                    results.append(WindowRow(
                        family=fam, N=n0, rep=rep,
                        center_layer=center,
                        n_layers_total=n_layers,
                        n_nodes=m,
                        total_relations=tr,
                        bd_ratio=bd_ratio_metric(sub_counts),
                        bdg_d2c=bdg_action_d2_corrected(sub_counts, m, normalized=True),
                        bdg_d4=bdg_action_d4_standard(sub_counts, m, normalized=True),
                        link_frac=float(sub_counts.get(0)) / tr if tr > 0 else 0.0,
                    ))

                if count % 100 == 0 or count == total:
                    print(f"  [B: {count}/{total}] {fam} N={n0} rep={rep}: "
                          f"{n_layers} layers")

    return results


# ── Main ─────────────────────────────────────────────────────────────

def safe_spearman(x, y):
    if np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return 0.0, 1.0
    r, p = sp_stats.spearmanr(x, y)
    return float(r), float(p)


def main() -> int:
    ap = argparse.ArgumentParser(description="Conjecture E: ρ_BD vs R(p) curvature cross-check")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--out-csv", default="outputs_unified_functional/conjecture_e_rho_bd_vs_curvature.csv")
    ap.add_argument("--out-md", default="outputs_unified_functional/conjecture_e_rho_bd_vs_curvature.md")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    if args.quick:
        families_a = ["Lor2D", "Lor3D", "Lor4D", "KR_like"]
        n_values_a = [20, 28, 36]
        reps = 10
        families_b = ["Lor2D"]
        n_values_b = [36, 48]
    else:
        families_a = ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like", "TransPerc"]
        n_values_a = [20, 28, 36, 48]
        reps = args.reps
        families_b = ["Lor2D", "TransPerc"]
        n_values_b = [36, 48]

    print(f"=== Conjecture E: ρ_BD vs R(p) Cross-Check ===")
    print(f"Strategy A families: {families_a}, N: {n_values_a}")
    print(f"Strategy B families: {families_b}, N: {n_values_b}")
    print(f"Reps: {reps}")

    # ── Strategy A ──
    print("\n--- Strategy A: Cross-poset global correlation ---")
    global_rows = run_strategy_a(families_a, n_values_a, reps, args.seed)

    # ── Strategy B ──
    print("\n--- Strategy B: Within-poset sliding-window ---")
    window_rows = run_strategy_b(families_b, n_values_b, reps, args.seed, window=1)

    # ── Save CSVs ──
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(global_rows[0]).keys()))
        w.writeheader()
        for r in global_rows:
            w.writerow(asdict(r))
    print(f"\nSaved: {out_csv}")

    out_csv_b = out_csv.with_name("conjecture_e_rho_bd_vs_curvature_windows.csv")
    if window_rows:
        with out_csv_b.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(asdict(window_rows[0]).keys()))
            w.writeheader()
            for r in window_rows:
                w.writerow(asdict(r))
        print(f"Saved: {out_csv_b}")

    # ── Analysis A: Cross-poset ──
    print(f"\n{'='*72}")
    print("STRATEGY A: CROSS-POSET GLOBAL CORRELATION")
    print(f"{'='*72}")

    ga = [asdict(r) for r in global_rows]

    # Within-family Spearman(bd_ratio, bdg_d4)
    print(f"\n--- Spearman ρ(bd_ratio, bdg_d4) within each family ---")
    print(f"{'Family':>10s} | {'ρ':>6s} | {'p-val':>7s} | {'n':>4s} | {'mean_bd':>8s} | {'mean_d4':>8s}")
    for fam in families_a:
        sub = [d for d in ga if d["family"] == fam]
        if len(sub) < 5:
            continue
        bd_arr = np.array([d["bd_ratio"] for d in sub])
        d4_arr = np.array([d["bdg_d4"] for d in sub])
        r, p = safe_spearman(bd_arr, d4_arr)
        print(f"{fam:>10s} | {r:+6.3f} | {p:7.4f} | {len(sub):4d} | {np.mean(bd_arr):8.4f} | {np.mean(d4_arr):+8.3f}")

    # Within-family Spearman(bd_ratio, link_frac)
    print(f"\n--- Spearman ρ(bd_ratio, link_frac) within each family ---")
    print(f"{'Family':>10s} | {'ρ':>6s} | {'p-val':>7s} | {'n':>4s}")
    for fam in families_a:
        sub = [d for d in ga if d["family"] == fam]
        if len(sub) < 5:
            continue
        bd_arr = np.array([d["bd_ratio"] for d in sub])
        lf_arr = np.array([d["link_frac"] for d in sub])
        r, p = safe_spearman(bd_arr, lf_arr)
        print(f"{fam:>10s} | {r:+6.3f} | {p:7.4f} | {len(sub):4d}")

    # Cross-family: pool all
    print(f"\n--- Cross-family: all posets pooled ---")
    all_bd = np.array([d["bd_ratio"] for d in ga])
    all_d4 = np.array([d["bdg_d4"] for d in ga])
    all_lf = np.array([d["link_frac"] for d in ga])
    all_d2c = np.array([d["bdg_d2c"] for d in ga])
    r1, p1 = safe_spearman(all_bd, all_d4)
    r2, p2 = safe_spearman(all_bd, all_lf)
    r3, p3 = safe_spearman(all_bd, all_d2c)
    print(f"  ρ(bd_ratio, bdg_d4)  = {r1:+.3f}  (p = {p1:.2e}, n = {len(ga)})")
    print(f"  ρ(bd_ratio, link_frac) = {r2:+.3f}  (p = {p2:.2e})")
    print(f"  ρ(bd_ratio, bdg_d2c)  = {r3:+.3f}  (p = {p3:.2e})")

    # Per-family mean bd_ratio vs bdg_d4 summary
    print(f"\n--- Family-level means ---")
    print(f"{'Family':>10s} | {'bd_ratio':>9s} | {'bdg_d4':>8s} | {'bdg_d2c':>8s} | {'link_frac':>9s} | {'n_layers':>8s}")
    for fam in families_a:
        sub = [d for d in ga if d["family"] == fam]
        if not sub:
            continue
        print(f"{fam:>10s} | {np.mean([d['bd_ratio'] for d in sub]):9.4f} "
              f"| {np.mean([d['bdg_d4'] for d in sub]):+8.3f} "
              f"| {np.mean([d['bdg_d2c'] for d in sub]):+8.3f} "
              f"| {np.mean([d['link_frac'] for d in sub]):9.4f} "
              f"| {np.mean([d['n_layers'] for d in sub]):8.1f}")

    # ── Analysis B: Within-poset ──
    if window_rows:
        print(f"\n{'='*72}")
        print("STRATEGY B: WITHIN-POSET SLIDING WINDOW")
        print(f"{'='*72}")

        wb = [asdict(r) for r in window_rows]
        # Group by (family, N, rep) and compute within-poset Spearman
        from itertools import groupby
        key_fn = lambda d: (d["family"], d["N"], d["rep"])
        wb_sorted = sorted(wb, key=key_fn)

        within_results = []
        for key, group in groupby(wb_sorted, key=key_fn):
            items = list(group)
            # Filter: only interior windows (not boundary layers with bd=0)
            interior = [d for d in items if d["bd_ratio"] > 0]
            if len(interior) < 3:
                continue
            bd_arr = np.array([d["bd_ratio"] for d in interior])
            d4_arr = np.array([d["bdg_d4"] for d in interior])
            lf_arr = np.array([d["link_frac"] for d in interior])
            r_d4, p_d4 = safe_spearman(bd_arr, d4_arr)
            r_lf, p_lf = safe_spearman(bd_arr, lf_arr)
            within_results.append({
                "family": key[0], "N": key[1], "rep": key[2],
                "n_windows": len(interior),
                "rho_bd_d4": r_d4, "p_bd_d4": p_d4,
                "rho_bd_lf": r_lf, "p_bd_lf": p_lf,
            })

        if within_results:
            print(f"\n--- Within-poset Spearman (interior windows only, bd>0) ---")
            print(f"{'Family':>10s} | {'N':>3s} | {'mean ρ(bd,d4)':>13s} | {'mean ρ(bd,lf)':>13s} | {'n_posets':>8s}")
            for fam in families_b:
                for n in n_values_b:
                    sub = [d for d in within_results if d["family"] == fam and d["N"] == n]
                    if not sub:
                        continue
                    r_d4_arr = np.array([d["rho_bd_d4"] for d in sub])
                    r_lf_arr = np.array([d["rho_bd_lf"] for d in sub])
                    print(f"{fam:>10s} | {n:>3d} | {np.mean(r_d4_arr):+13.3f} | {np.mean(r_lf_arr):+13.3f} | {len(sub):>8d}")

        # All interior windows pooled
        interior_all = [d for d in wb if d["bd_ratio"] > 0]
        if len(interior_all) >= 5:
            bd_all = np.array([d["bd_ratio"] for d in interior_all])
            d4_all = np.array([d["bdg_d4"] for d in interior_all])
            lf_all = np.array([d["link_frac"] for d in interior_all])
            r_d4, _ = safe_spearman(bd_all, d4_all)
            r_lf, _ = safe_spearman(bd_all, lf_all)
            print(f"\n  Pooled interior windows: ρ(bd,d4)={r_d4:+.3f}, ρ(bd,lf)={r_lf:+.3f} (n={len(interior_all)})")

    # ── Markdown report ──
    md = Path(args.out_md)
    lines = []
    lines.append("# Conjecture E: ρ_BD vs R(p) Curvature Cross-Check")
    lines.append("")
    lines.append("## Question")
    lines.append("")
    lines.append("Does the BD action density (ρ_BD, Layer 3 of Conjecture E) co-vary with")
    lines.append("the discrete curvature proxy (bdg_d4, Layer 2) — both globally across")
    lines.append("poset instances and locally within individual posets?")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append("**Strategy A (cross-poset):** For each poset instance, compute global")
    lines.append("bd_ratio and global bdg_d4. Spearman-correlate across instances.")
    lines.append("")
    lines.append("**Strategy B (within-poset):** For deep-layer families (Lor2D, TransPerc),")
    lines.append("use a sliding window (center ± 1 layer) to compute local observables.")
    lines.append("Correlate across windows within each poset.")
    lines.append("")
    lines.append("Note: High-d families (Lor4D/5D, KR_like) have only ~3 layers, with")
    lines.append("boundary layers having bd_ratio=0. Within-poset correlation is impossible.")
    lines.append("")

    # Strategy A results
    lines.append("## Strategy A: Cross-poset results")
    lines.append("")
    lines.append("### Within-family Spearman ρ(bd_ratio, bdg_d4)")
    lines.append("")
    lines.append("| family | ρ | p-value | n | mean bd | mean d4 |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for fam in families_a:
        sub = [d for d in ga if d["family"] == fam]
        if len(sub) < 5:
            continue
        bd_arr = np.array([d["bd_ratio"] for d in sub])
        d4_arr = np.array([d["bdg_d4"] for d in sub])
        r, p = safe_spearman(bd_arr, d4_arr)
        lines.append(f"| {fam} | {r:+.3f} | {p:.4f} | {len(sub)} | {np.mean(bd_arr):.4f} | {np.mean(d4_arr):+.3f} |")

    lines.append("")
    lines.append(f"### Cross-family pooled: ρ(bd_ratio, bdg_d4) = {r1:+.3f} (p = {p1:.2e}, n = {len(ga)})")
    lines.append(f"### Cross-family pooled: ρ(bd_ratio, link_frac) = {r2:+.3f} (p = {p2:.2e})")
    lines.append("")

    # Strategy B results
    if window_rows and within_results:
        lines.append("## Strategy B: Within-poset results")
        lines.append("")
        lines.append("| family | N | mean ρ(bd,d4) | mean ρ(bd,lf) | n_posets |")
        lines.append("|---|---:|---:|---:|---:|")
        for fam in families_b:
            for n in n_values_b:
                sub = [d for d in within_results if d["family"] == fam and d["N"] == n]
                if not sub:
                    continue
                r_d4_arr = np.array([d["rho_bd_d4"] for d in sub])
                r_lf_arr = np.array([d["rho_bd_lf"] for d in sub])
                lines.append(f"| {fam} | {n} | {np.mean(r_d4_arr):+.3f} | {np.mean(r_lf_arr):+.3f} | {len(sub)} |")
        lines.append("")

    lines.append("## Interpretation")
    lines.append("")
    lines.append("- Positive ρ(bd_ratio, bdg_d4): regions with richer interval structure")
    lines.append("  also have higher BDG d=4 action — the BD-to-EH bridge is locally consistent.")
    lines.append("- Negative ρ(bd_ratio, link_frac): more non-trivial intervals (fewer links)")
    lines.append("  → richer BD structure — bd_ratio measures 'causal thickness'.")
    lines.append("- Cross-family correlations reflect that different families have")
    lines.append("  systematically different bd_ratio and curvature profiles.")

    md.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"\nSaved: {md}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
