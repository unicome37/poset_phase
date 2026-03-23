"""Prediction E — Time Arrow from Structural Asymmetry.

E1: Augmentation entropy asymmetry (CST-style partial augmentation).
  From a base poset X₀, add k nodes:
  - Forward: new node is in the FUTURE of a random subset of top-layer elements
  - Backward: new node is in the PAST of a random subset of bottom-layer elements
  - Random: control baseline
  Measure: A_entropy = Δlog H(forward) - Δlog H(backward).
  Prediction: A_entropy > 0 (forward augmentation increases entropy more).

  CRITICAL FIX (v2): The v1 operators connected to ALL maxima/minima, which
  trivially gives H(X ∪ {top}) = H(X) — a standard poset identity.
  The v2 operators use CST-style PARTIAL augmentation: new node connects
  to a random subset of elements in the top/bottom half of the layer
  structure, breaking the symmetry.

E2: Causal depth directional growth.
  Multi-step augmentation (T steps), track layer depth ℓ(X_t).
  Measure: R_depth = (ℓ_forward - ℓ_0) / (ℓ_backward - ℓ_0).
  Prediction: R_depth > 1 (forward growth deepens structure faster).

E3: Coarse-graining robustness.
  After augmentation, coarse-grain and check if A_entropy preserves sign.

Outputs:
  - outputs_unified_functional/prediction_e_e1.csv
  - outputs_unified_functional/prediction_e_e1.md
"""

from __future__ import annotations

import argparse
import csv
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from coarse_grain import coarse_grain_delete_nodes
from entropy_exact import log_linear_extensions_exact
from generators import Poset, transitive_closure
from generators import (
    generate_kr_like,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_transitive_percolation,
)
from observables import layer_profile


# ── Layer assignment helper ──────────────────────────────────────────────

def _layer_assignment(closure: np.ndarray) -> np.ndarray:
    """Return per-node layer index (0 = minimal elements)."""
    n = closure.shape[0]
    indeg = closure.sum(axis=0).astype(int)
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
            for v in np.where(closure[u])[0]:
                indeg[v] -= 1
        layer_idx += 1
    return labels


# ── CST-style augmentation operators (v2) ────────────────────────────────

def augment_forward(poset: Poset, k: int = 1, seed: int = 0, p_connect: float = 0.6) -> Poset:
    """CST-style forward augmentation: add k nodes to the FUTURE.

    Each new node connects to a random subset of elements in the TOP half
    of the layer structure (probability p_connect per candidate).
    Then transitive closure propagates all predecessor relationships.

    This is NOT a top-element addition (which trivially preserves H).
    The partial connectivity creates genuine structural change.
    """
    rng = np.random.RandomState(seed)
    closure = poset.closure.copy()

    for step in range(k):
        cur_n = closure.shape[0]
        layers = _layer_assignment(closure)
        n_layers = int(layers.max()) + 1
        threshold = max(n_layers // 2, 1)

        # Candidates: elements in top half of layers
        candidates = np.where(layers >= threshold)[0]
        if len(candidates) == 0:
            candidates = np.where(layers == layers.max())[0]

        new_c = np.zeros((cur_n + 1, cur_n + 1), dtype=bool)
        new_c[:cur_n, :cur_n] = closure
        new_idx = cur_n

        # Each candidate becomes a predecessor of new node with prob p_connect
        selected = [c for c in candidates if rng.random() < p_connect]
        if not selected:
            selected = [candidates[rng.randint(len(candidates))]]

        for c in selected:
            new_c[c, new_idx] = True
            # Transitivity: predecessors of c also precede new_idx
            preds = np.where(closure[:, c])[0]
            new_c[preds, new_idx] = True

        np.fill_diagonal(new_c, False)
        closure = new_c

    return Poset(closure)


def augment_backward(poset: Poset, k: int = 1, seed: int = 0, p_connect: float = 0.6) -> Poset:
    """CST-style backward augmentation: add k nodes to the PAST.

    Each new node connects to a random subset of elements in the BOTTOM half
    of the layer structure. Then transitive closure propagates successor
    relationships.
    """
    rng = np.random.RandomState(seed)
    closure = poset.closure.copy()

    for step in range(k):
        cur_n = closure.shape[0]
        layers = _layer_assignment(closure)
        n_layers = int(layers.max()) + 1
        threshold = n_layers // 2

        # Candidates: elements in bottom half of layers
        candidates = np.where(layers < threshold)[0]
        if len(candidates) == 0:
            candidates = np.where(layers == 0)[0]

        new_c = np.zeros((cur_n + 1, cur_n + 1), dtype=bool)
        new_c[:cur_n, :cur_n] = closure
        new_idx = cur_n

        selected = [c for c in candidates if rng.random() < p_connect]
        if not selected:
            selected = [candidates[rng.randint(len(candidates))]]

        for c in selected:
            new_c[new_idx, c] = True
            # Transitivity: successors of c also succeed new_idx
            succs = np.where(closure[c, :])[0]
            new_c[new_idx, succs] = True

        np.fill_diagonal(new_c, False)
        closure = new_c

    return Poset(closure)


def augment_random(poset: Poset, k: int = 1, seed: int = 0) -> Poset:
    """Random augmentation (control): add k nodes at random layer positions.

    Each new node connects to a random subset of ALL elements (not biased
    toward top or bottom), with probability scaled by comparable fraction.
    """
    rng = np.random.RandomState(seed)
    closure = poset.closure.copy()

    for step in range(k):
        cur_n = closure.shape[0]
        layers = _layer_assignment(closure)
        n_layers = int(layers.max()) + 1

        # Pick a random layer threshold to split predecessors/successors
        split_layer = rng.randint(0, max(n_layers, 1))

        new_c = np.zeros((cur_n + 1, cur_n + 1), dtype=bool)
        new_c[:cur_n, :cur_n] = closure
        new_idx = cur_n

        # Elements below split become predecessors, above become successors
        below = np.where(layers < split_layer)[0]
        for b in below:
            if rng.random() < 0.4:
                new_c[b, new_idx] = True
                preds = np.where(closure[:, b])[0]
                new_c[preds, new_idx] = True

        above = np.where(layers > split_layer)[0]
        for a in above:
            if rng.random() < 0.4:
                new_c[new_idx, a] = True
                succs = np.where(closure[a, :])[0]
                new_c[new_idx, succs] = True

        np.fill_diagonal(new_c, False)
        closure = new_c

    return Poset(closure)


# ── Entropy computation ──────────────────────────────────────────────────

def compute_log_H(poset: Poset, use_exact: bool = True, n_sis: int = 256) -> float:
    """Compute log H using exact DP (if N <= max_exact) or SIS estimator."""
    if use_exact and poset.n <= 22:
        return log_linear_extensions_exact(poset)
    else:
        from unified_functional import compute_log_H as sis_log_H
        return sis_log_H(poset, n_runs=n_sis)


# ── E1 experiment ────────────────────────────────────────────────────────

@dataclass(frozen=True)
class E1Row:
    family: str
    N0: int
    k: int
    rep: int
    seed: int
    log_H_base: float
    log_H_forward: float
    log_H_backward: float
    log_H_random: float
    delta_forward: float
    delta_backward: float
    delta_random: float
    A_entropy: float
    n_layers_base: int
    n_layers_forward: int
    n_layers_backward: int
    n_maxima: int
    n_minima: int


FAMILY_GENS = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "KR_like": generate_kr_like,
    "TransPerc": lambda n, seed=None: generate_transitive_percolation(n, p=0.08, seed=seed),
}


def run_e1(
    families: list[str],
    n_values: list[int],
    k_values: list[int],
    reps: int,
    base_seed: int,
    use_exact: bool,
) -> list[E1Row]:
    """Run the E1 augmentation asymmetry experiment."""
    results = []
    total = len(families) * len(n_values) * len(k_values) * reps
    count = 0

    for fam in families:
        gen = FAMILY_GENS[fam]
        for n0 in n_values:
            for k in k_values:
                for rep in range(reps):
                    count += 1
                    s = base_seed + rep * 1000 + n0 * 100 + k * 10
                    poset = gen(n0, seed=s)

                    # Base entropy
                    log_H_base = compute_log_H(poset, use_exact=use_exact)

                    # Forward augmentation
                    p_fwd = augment_forward(poset, k=k, seed=s)
                    log_H_fwd = compute_log_H(p_fwd, use_exact=use_exact)

                    # Backward augmentation
                    p_bwd = augment_backward(poset, k=k, seed=s)
                    log_H_bwd = compute_log_H(p_bwd, use_exact=use_exact)

                    # Random augmentation (control)
                    p_rnd = augment_random(poset, k=k, seed=s)
                    log_H_rnd = compute_log_H(p_rnd, use_exact=use_exact)

                    d_fwd = log_H_fwd - log_H_base
                    d_bwd = log_H_bwd - log_H_base
                    d_rnd = log_H_rnd - log_H_base
                    A = d_fwd - d_bwd

                    # Layer counts
                    nl_base = len(layer_profile(poset))
                    nl_fwd = len(layer_profile(p_fwd))
                    nl_bwd = len(layer_profile(p_bwd))

                    # Maxima / minima counts
                    has_succ = poset.closure.any(axis=1)
                    has_pred = poset.closure.any(axis=0)
                    n_max = int((~has_succ).sum())
                    n_min = int((~has_pred).sum())

                    results.append(E1Row(
                        family=fam, N0=n0, k=k, rep=rep, seed=s,
                        log_H_base=log_H_base,
                        log_H_forward=log_H_fwd,
                        log_H_backward=log_H_bwd,
                        log_H_random=log_H_rnd,
                        delta_forward=d_fwd,
                        delta_backward=d_bwd,
                        delta_random=d_rnd,
                        A_entropy=A,
                        n_layers_base=nl_base,
                        n_layers_forward=nl_fwd,
                        n_layers_backward=nl_bwd,
                        n_maxima=n_max,
                        n_minima=n_min,
                    ))

                    if count % 50 == 0:
                        print(f"  [{count}/{total}] {fam} N={n0} k={k} rep={rep}: "
                              f"A={A:+.3f} (Δfwd={d_fwd:+.3f}, Δbwd={d_bwd:+.3f})")

    return results


# ── E2 experiment (multi-step depth tracking) ────────────────────────────

@dataclass(frozen=True)
class E2Row:
    family: str
    N0: int
    T: int
    rep: int
    seed: int
    layers_base: int
    layers_forward_T: int
    layers_backward_T: int
    depth_gain_forward: int
    depth_gain_backward: int
    R_depth: float


def run_e2(
    families: list[str],
    n_values: list[int],
    T_values: list[int],
    reps: int,
    base_seed: int,
) -> list[E2Row]:
    """Run E2: multi-step directional growth, track layer depth."""
    results = []
    total = len(families) * len(n_values) * len(T_values) * reps
    count = 0

    for fam in families:
        gen = FAMILY_GENS[fam]
        for n0 in n_values:
            for T in T_values:
                for rep in range(reps):
                    count += 1
                    s = base_seed + rep * 1000 + n0 * 100 + T * 10
                    poset = gen(n0, seed=s)

                    nl_base = len(layer_profile(poset))

                    # Forward T steps (1 node per step)
                    p_fwd = poset
                    for t in range(T):
                        p_fwd = augment_forward(p_fwd, k=1, seed=s + t)
                    nl_fwd = len(layer_profile(p_fwd))

                    # Backward T steps
                    p_bwd = poset
                    for t in range(T):
                        p_bwd = augment_backward(p_bwd, k=1, seed=s + t)
                    nl_bwd = len(layer_profile(p_bwd))

                    dg_fwd = nl_fwd - nl_base
                    dg_bwd = nl_bwd - nl_base
                    R = dg_fwd / dg_bwd if dg_bwd > 0 else (
                        float("inf") if dg_fwd > 0 else float("nan")
                    )

                    results.append(E2Row(
                        family=fam, N0=n0, T=T, rep=rep, seed=s,
                        layers_base=nl_base,
                        layers_forward_T=nl_fwd,
                        layers_backward_T=nl_bwd,
                        depth_gain_forward=dg_fwd,
                        depth_gain_backward=dg_bwd,
                        R_depth=R,
                    ))

                    if count % 50 == 0:
                        print(f"  [{count}/{total}] {fam} N={n0} T={T} rep={rep}: "
                              f"R_depth={R:.3f} (fwd={dg_fwd}, bwd={dg_bwd})")

    return results


# ── E3: Coarse-graining robustness ──────────────────────────────────────

def e3_cg_robustness(e1_rows: list[E1Row], base_seed: int, use_exact: bool) -> list[dict]:
    """Check if A_entropy preserves sign after coarse-graining."""
    results = []
    # Only test a subset (k=2, all families, first 5 reps)
    subset = [r for r in e1_rows if r.k == 2 and r.rep < 5]

    for i, r in enumerate(subset):
        s = r.seed
        poset = FAMILY_GENS[r.family](r.N0, seed=s)

        p_fwd = augment_forward(poset, k=r.k, seed=s)
        p_bwd = augment_backward(poset, k=r.k, seed=s)

        # Coarse-grain (keep 70%)
        p_fwd_cg = coarse_grain_delete_nodes(p_fwd, keep_ratio=0.7, seed=s + 9999)
        p_bwd_cg = coarse_grain_delete_nodes(p_bwd, keep_ratio=0.7, seed=s + 9999)
        p_base_cg = coarse_grain_delete_nodes(poset, keep_ratio=0.7, seed=s + 9999)

        log_H_base_cg = compute_log_H(p_base_cg, use_exact=use_exact)
        log_H_fwd_cg = compute_log_H(p_fwd_cg, use_exact=use_exact)
        log_H_bwd_cg = compute_log_H(p_bwd_cg, use_exact=use_exact)

        d_fwd_cg = log_H_fwd_cg - log_H_base_cg
        d_bwd_cg = log_H_bwd_cg - log_H_base_cg
        A_cg = d_fwd_cg - d_bwd_cg

        results.append({
            "family": r.family, "N0": r.N0, "k": r.k, "rep": r.rep,
            "A_entropy_raw": r.A_entropy,
            "A_entropy_cg": A_cg,
            "sign_preserved": (r.A_entropy > 0) == (A_cg > 0) if abs(r.A_entropy) > 0.01 else True,
        })

    return results


# ── Main ─────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(description="Prediction E: Time Arrow from Structural Asymmetry")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--out-csv", default="outputs_unified_functional/prediction_e_e1.csv")
    ap.add_argument("--out-md", default="outputs_unified_functional/prediction_e_e1.md")
    ap.add_argument("--quick", action="store_true", help="Reduced run for testing")
    args = ap.parse_args()

    if args.quick:
        families = ["Lor2D", "Lor4D", "KR_like"]
        n_values = [12, 16]
        k_values = [1, 2]
        T_values = [2, 4]
        reps = 5
    else:
        families = ["Lor2D", "Lor3D", "Lor4D", "KR_like", "TransPerc"]
        n_values = [12, 16, 20]
        k_values = [1, 2, 4]
        T_values = [2, 4, 8]
        reps = args.reps

    print(f"=== Prediction E: Time Arrow Experiment ===")
    print(f"Families: {families}")
    print(f"N values: {n_values}")
    print(f"k values (E1): {k_values}, T values (E2): {T_values}")
    print(f"Reps: {reps}")

    t0 = time.time()

    # ── E1 ──
    print(f"\n--- E1: Augmentation entropy asymmetry ---")
    e1_rows = run_e1(families, n_values, k_values, reps, args.seed, use_exact=True)

    # ── E2 ──
    print(f"\n--- E2: Causal depth directional growth ---")
    e2_rows = run_e2(families, n_values, T_values, reps, args.seed)

    # ── E3 ──
    print(f"\n--- E3: Coarse-graining robustness ---")
    e3_rows = e3_cg_robustness(e1_rows, args.seed, use_exact=True)

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed:.1f}s")

    # ── Save E1 CSV ──
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(e1_rows[0]).keys()))
        w.writeheader()
        for r in e1_rows:
            w.writerow(asdict(r))
    print(f"Saved: {out_csv}")

    # Save E2 CSV
    e2_csv = out_csv.with_name("prediction_e_e2.csv")
    with e2_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(e2_rows[0]).keys()))
        w.writeheader()
        for r in e2_rows:
            w.writerow(asdict(r))
    print(f"Saved: {e2_csv}")

    # ── Analysis ──
    data1 = [asdict(r) for r in e1_rows]
    data2 = [asdict(r) for r in e2_rows]

    print(f"\n{'='*70}")
    print("E1: AUGMENTATION ENTROPY ASYMMETRY")
    print(f"{'='*70}")

    # Per-family A_entropy
    print(f"\n{'Family':>10s} | {'mean A':>8s} | {'std A':>7s} | {'A>0 %':>6s} | {'n':>4s}")
    for fam in families:
        sub = [d for d in data1 if d["family"] == fam]
        A_arr = np.array([d["A_entropy"] for d in sub])
        pct_pos = 100.0 * np.mean(A_arr > 0)
        print(f"{fam:>10s} | {np.mean(A_arr):+8.4f} | {np.std(A_arr):7.4f} | {pct_pos:5.1f}% | {len(sub):4d}")

    # Per-k breakdown
    print(f"\n{'k':>4s} | {'mean A':>8s} | {'std A':>7s} | {'A>0 %':>6s}")
    for k in k_values:
        sub = [d for d in data1 if d["k"] == k]
        A_arr = np.array([d["A_entropy"] for d in sub])
        pct_pos = 100.0 * np.mean(A_arr > 0)
        print(f"{k:>4d} | {np.mean(A_arr):+8.4f} | {np.std(A_arr):7.4f} | {pct_pos:5.1f}%")

    # Per-N breakdown
    print(f"\n{'N':>4s} | {'mean A':>8s} | {'std A':>7s} | {'A>0 %':>6s}")
    for n in n_values:
        sub = [d for d in data1 if d["N0"] == n]
        A_arr = np.array([d["A_entropy"] for d in sub])
        pct_pos = 100.0 * np.mean(A_arr > 0)
        print(f"{n:>4d} | {np.mean(A_arr):+8.4f} | {np.std(A_arr):7.4f} | {pct_pos:5.1f}%")

    # Detailed family × k table
    print(f"\n--- Family × k table (mean A_entropy) ---")
    header = f"{'Family':>10s}"
    for k in k_values:
        header += f" | k={k:>2d}"
    print(header)
    for fam in families:
        row = f"{fam:>10s}"
        for k in k_values:
            sub = [d for d in data1 if d["family"] == fam and d["k"] == k]
            A_arr = np.array([d["A_entropy"] for d in sub])
            row += f" | {np.mean(A_arr):+6.3f}"
        print(row)

    print(f"\n{'='*70}")
    print("E2: CAUSAL DEPTH DIRECTIONAL GROWTH")
    print(f"{'='*70}")

    # Per-family R_depth
    print(f"\n{'Family':>10s} | {'mean R':>8s} | {'std R':>7s} | {'R>1 %':>6s} | {'n':>4s}")
    for fam in families:
        sub = [d for d in data2 if d["family"] == fam and math.isfinite(d["R_depth"])]
        if not sub:
            continue
        R_arr = np.array([d["R_depth"] for d in sub])
        pct_gt1 = 100.0 * np.mean(R_arr > 1.0)
        print(f"{fam:>10s} | {np.mean(R_arr):8.3f} | {np.std(R_arr):7.3f} | {pct_gt1:5.1f}% | {len(sub):4d}")

    # Per-T breakdown
    print(f"\n{'T':>4s} | {'mean R':>8s} | {'R>1 %':>6s}")
    for T in T_values:
        sub = [d for d in data2 if d["T"] == T and math.isfinite(d["R_depth"])]
        if not sub:
            continue
        R_arr = np.array([d["R_depth"] for d in sub])
        pct_gt1 = 100.0 * np.mean(R_arr > 1.0)
        print(f"{T:>4d} | {np.mean(R_arr):8.3f} | {pct_gt1:5.1f}%")

    print(f"\n{'='*70}")
    print("E3: COARSE-GRAINING ROBUSTNESS")
    print(f"{'='*70}")

    if e3_rows:
        n_preserved = sum(1 for r in e3_rows if r["sign_preserved"])
        n_total = len(e3_rows)
        print(f"Sign preserved: {n_preserved}/{n_total} ({100*n_preserved/n_total:.1f}%)")
        for fam in families:
            sub = [r for r in e3_rows if r["family"] == fam]
            if sub:
                n_p = sum(1 for r in sub if r["sign_preserved"])
                print(f"  {fam:>10s}: {n_p}/{len(sub)} preserved")

    # ── Save markdown report ──
    md = Path(args.out_md)
    lines = []
    lines.append("# Prediction E: Time Arrow from Structural Asymmetry")
    lines.append("")
    lines.append("## E1: Augmentation Entropy Asymmetry")
    lines.append("")
    lines.append("A_entropy = Δlog H(forward) - Δlog H(backward)")
    lines.append("")
    lines.append("Prediction: A_entropy > 0 (forward augmentation increases entropy more)")
    lines.append("")

    lines.append("### Per-family results")
    lines.append("")
    lines.append("| family | mean A | std | A>0 % | n |")
    lines.append("|---|---:|---:|---:|---:|")
    for fam in families:
        sub = [d for d in data1 if d["family"] == fam]
        A_arr = np.array([d["A_entropy"] for d in sub])
        pct_pos = 100.0 * np.mean(A_arr > 0)
        lines.append(f"| {fam} | {np.mean(A_arr):+.4f} | {np.std(A_arr):.4f} | {pct_pos:.1f}% | {len(sub)} |")

    lines.append("")
    lines.append("### Family × k table (mean A_entropy)")
    lines.append("")
    hdr = "| family |"
    for k in k_values:
        hdr += f" k={k} |"
    lines.append(hdr)
    lines.append("|---|" + "---:|" * len(k_values))
    for fam in families:
        row = f"| {fam} |"
        for k in k_values:
            sub = [d for d in data1 if d["family"] == fam and d["k"] == k]
            A_arr = np.array([d["A_entropy"] for d in sub])
            row += f" {np.mean(A_arr):+.3f} |"
        lines.append(row)

    lines.append("")
    lines.append("## E2: Causal Depth Directional Growth")
    lines.append("")
    lines.append("R_depth = (layers_forward - layers_base) / (layers_backward - layers_base)")
    lines.append("")
    lines.append("Prediction: R_depth > 1")
    lines.append("")
    lines.append("| family | mean R | std | R>1 % | n |")
    lines.append("|---|---:|---:|---:|---:|")
    for fam in families:
        sub = [d for d in data2 if d["family"] == fam and math.isfinite(d["R_depth"])]
        if not sub:
            continue
        R_arr = np.array([d["R_depth"] for d in sub])
        pct_gt1 = 100.0 * np.mean(R_arr > 1.0)
        lines.append(f"| {fam} | {np.mean(R_arr):.3f} | {np.std(R_arr):.3f} | {pct_gt1:.1f}% | {len(sub)} |")

    lines.append("")
    lines.append("## E3: Coarse-graining Robustness")
    lines.append("")
    if e3_rows:
        n_preserved = sum(1 for r in e3_rows if r["sign_preserved"])
        n_total = len(e3_rows)
        lines.append(f"Sign preservation rate: {n_preserved}/{n_total} ({100*n_preserved/n_total:.1f}%)")
        lines.append("")
        lines.append("| family | preserved | total |")
        lines.append("|---|---:|---:|")
        for fam in families:
            sub = [r for r in e3_rows if r["family"] == fam]
            if sub:
                n_p = sum(1 for r in sub if r["sign_preserved"])
                lines.append(f"| {fam} | {n_p} | {len(sub)} |")

    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("If A_entropy > 0 consistently across families and N values, this supports")
    lines.append("the claim that time direction emerges from structural growth asymmetry:")
    lines.append("augmenting toward the future (adding successors of maxima) increases")
    lines.append("combinatorial entropy more than augmenting toward the past.")
    lines.append("")
    lines.append("R_depth > 1 would confirm that forward growth deepens causal structure")
    lines.append("faster than backward growth — a complementary structural arrow of time.")

    md.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"\nSaved: {md}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
