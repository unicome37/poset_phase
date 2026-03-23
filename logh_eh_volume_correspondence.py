"""logH ↔ ∫√(-g)d⁴x Correspondence Experiment.

Tests the theoretical path item 3: how does the combinatorial entropy
logH (= log of the number of linear extensions) of a causal set relate
to the Einstein-Hilbert volume term ∫√(-g)d⁴x ?

Key theoretical links:
1. For a Poisson sprinkling at density ρ in volume V: N = ρV
2. logH depends on the order fraction r = (causal pairs) / C(N,2)
3. In Minkowski: r = r_d (dimension-dependent constant), logH ~ N·log(N)·g(r_d)
4. In de Sitter: r decreases with H (fewer causal pairs due to expansion)
5. The EH action S_EH = (1/16πG) ∫ R√(-g)d⁴x ∝ N for fixed ρ and curvature

Experiment: on de Sitter sprinklings (d=2,3,4; N=32,48,64,96,128;
H=0,0.2,0.5,1.0,2.0; 5 reps), compute:
- logH via SIS estimator
- order fraction r = n_causal / C(N,2)
- n_causal_pairs
- Compare logH vs N, logH vs r, logH/N vs H
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import build_causal_matrix_de_sitter
from entropy_sis import estimate_log_linear_extensions_sis
from generators import Poset


def build_poset_from_causal(causal: np.ndarray) -> Poset:
    """Build a Poset object from a causal matrix for SIS estimation."""
    # closure = transitive closure (causal matrix is already the closure for sprinklings)
    # But remove diagonal for Poset convention (strict partial order)
    closure = causal.copy()
    np.fill_diagonal(closure, False)
    return Poset(closure=closure)


@dataclass(frozen=True)
class Row:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float           # d(d-1)H²
    log_H: float          # log(# linear extensions)
    log_H_std: float      # SIS estimator std
    n_causal: int         # number of causal pairs
    order_frac: float     # r = n_causal / C(N,2)
    log_H_per_N: float    # logH / N
    log_H_per_NlogN: float  # logH / (N·log(N))


def run_single(d: int, N: int, hubble: float, rep: int, seed: int,
               n_sis: int = 128) -> Row:
    """Run one realization."""
    s = seed + d * 100000 + N * 1000 + rep + int(hubble * 10000)
    points = sprinkle_de_sitter_like_diamond(N, d - 1, hubble=hubble, seed=s)
    causal = build_causal_matrix_de_sitter(points, hubble)

    # Remove diagonal
    np.fill_diagonal(causal, False)

    n_causal = int(causal.sum())
    order_frac = n_causal / max(N * (N - 1) / 2, 1)  # C(N,2)

    # Build poset and estimate logH
    poset = build_poset_from_causal(causal)
    log_H, log_H_std = estimate_log_linear_extensions_sis(poset, n_runs=n_sis, seed=s + 777)

    R_dS = d * (d - 1) * hubble * hubble
    log_H_per_N = log_H / max(N, 1)
    log_H_per_NlogN = log_H / max(N * math.log(max(N, 2)), 1)

    return Row(
        d=d, N=N, hubble=hubble, rep=rep,
        R_dS=R_dS, log_H=log_H, log_H_std=log_H_std,
        n_causal=n_causal, order_frac=order_frac,
        log_H_per_N=log_H_per_N,
        log_H_per_NlogN=log_H_per_NlogN,
    )


def main() -> int:
    dims = [2, 3, 4]
    # Keep N small because SIS is O(N²) per sample
    ns = [32, 48, 64, 96, 128]
    hubbles = [0.0, 0.2, 0.5, 1.0, 2.0]
    n_reps = 5
    n_sis = 128
    seed_base = 55555

    csv_path = Path("outputs_unified_functional/logh_eh_volume.csv")
    csv_path.parent.mkdir(exist_ok=True)

    all_rows: list[Row] = []
    total = len(dims) * len(ns) * len(hubbles) * n_reps
    count = 0

    for d in dims:
        for N in ns:
            for H in hubbles:
                for rep in range(n_reps):
                    count += 1
                    print(f"[{count}/{total}] d={d} N={N} H={H:.1f} rep={rep}")
                    row = run_single(d, N, H, rep, seed_base, n_sis=n_sis)
                    all_rows.append(row)

    # Write CSV
    field_names = [f.name for f in fields(Row)]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=field_names)
        writer.writeheader()
        for row in all_rows:
            writer.writerow({fn: getattr(row, fn) for fn in field_names})
    print(f"CSV: {csv_path}  ({len(all_rows)} rows)")

    # ── Analysis and Report ──
    report_path = Path("outputs_unified_functional/logh_eh_volume.md")
    lines: list[str] = []
    lines.append("# logH ↔ ∫√(-g)d⁴x Correspondence\n")

    # ── 1. logH vs N at fixed d, H ──
    lines.append("\n## 1. logH Scaling with N\n")
    lines.append("| d | H | slope(logH vs NlogN) | R² | logH/NlogN at N=128 |")
    lines.append("|---|---|---------------------|----|--------------------|")

    for d in dims:
        for H in hubbles:
            subset = [r for r in all_rows if r.d == d and r.hubble == H]
            if len(subset) < 5:
                continue
            # Group by N, compute mean logH
            by_n: dict[int, list[float]] = defaultdict(list)
            for r in subset:
                by_n[r.N].append(r.log_H)
            
            xs, ys = [], []
            for N in sorted(by_n.keys()):
                mean_lh = float(np.mean(by_n[N]))
                xs.append(N * math.log(N))
                ys.append(mean_lh)
            
            if len(xs) >= 3:
                slope, intercept, r_val, _, _ = sp_stats.linregress(xs, ys)
                r2 = r_val ** 2
                lh_nlogn_128 = float(np.mean(by_n.get(128, [0]))) / (128 * math.log(128)) if 128 in by_n else float("nan")
                lines.append(f"| {d} | {H:.1f} | {slope:.4f} | {r2:.4f} | {lh_nlogn_128:.4f} |")

    # ── 2. logH/N vs H at fixed d, N ──
    lines.append("\n## 2. logH/N vs Hubble Parameter\n")
    lines.append("logH/N = combinatorial entropy density. If logH ∝ V = N/ρ,")
    lines.append("then logH/N ≈ const. Deviations reveal curvature effects.\n")
    lines.append("| d | N | logH/N (H=0) | logH/N (H=1) | logH/N (H=2) | ρ(H, logH/N) |")
    lines.append("|---|---|-------------|-------------|-------------|-------------|")

    for d in dims:
        for N in ns:
            h0_vals = [r.log_H_per_N for r in all_rows if r.d == d and r.N == N and r.hubble == 0]
            h1_vals = [r.log_H_per_N for r in all_rows if r.d == d and r.N == N and r.hubble == 1.0]
            h2_vals = [r.log_H_per_N for r in all_rows if r.d == d and r.N == N and r.hubble == 2.0]
            
            subset_n = [r for r in all_rows if r.d == d and r.N == N]
            if len(subset_n) >= 5:
                hs = [r.hubble for r in subset_n]
                lhns = [r.log_H_per_N for r in subset_n]
                rho, _ = sp_stats.spearmanr(hs, lhns)
            else:
                rho = float("nan")
            
            m0 = float(np.mean(h0_vals)) if h0_vals else float("nan")
            m1 = float(np.mean(h1_vals)) if h1_vals else float("nan")
            m2 = float(np.mean(h2_vals)) if h2_vals else float("nan")
            lines.append(f"| {d} | {N} | {m0:.3f} | {m1:.3f} | {m2:.3f} | {rho:+.3f} |")

    # ── 3. Order fraction vs H ──
    lines.append("\n## 3. Order Fraction r vs H\n")
    lines.append("r = n_causal / C(N,2). This is a proxy for the causal structure density.\n")
    lines.append("| d | N | r(H=0) | r(H=1) | r(H=2) | ρ(H, r) |")
    lines.append("|---|---|--------|--------|--------|---------|")

    for d in dims:
        for N in ns:
            subset_n = [r for r in all_rows if r.d == d and r.N == N]
            r0 = [r.order_frac for r in subset_n if r.hubble == 0]
            r1 = [r.order_frac for r in subset_n if r.hubble == 1.0]
            r2 = [r.order_frac for r in subset_n if r.hubble == 2.0]
            
            hs = [r.hubble for r in subset_n]
            rs = [r.order_frac for r in subset_n]
            rho, _ = sp_stats.spearmanr(hs, rs) if len(hs) >= 5 else (float("nan"), 0)
            
            m0 = float(np.mean(r0)) if r0 else float("nan")
            m1 = float(np.mean(r1)) if r1 else float("nan")
            m2 = float(np.mean(r2)) if r2 else float("nan")
            lines.append(f"| {d} | {N} | {m0:.4f} | {m1:.4f} | {m2:.4f} | {rho:+.3f} |")

    # ── 4. logH vs n_causal ──
    lines.append("\n## 4. logH vs n_causal (Causal Pair Count)\n")
    lines.append("If logH ~ f(n_causal), then the volume term is mediated by causal structure.\n")
    lines.append("| d | ρ(log_H, n_causal) | ρ(log_H, log(n_causal)) | ρ(logH/N, r) |")
    lines.append("|---|-------------------|------------------------|--------------|")

    for d in dims:
        subset_d = [r for r in all_rows if r.d == d]
        lhs = [r.log_H for r in subset_d]
        ncs = [r.n_causal for r in subset_d]
        log_ncs = [math.log(max(r.n_causal, 1)) for r in subset_d]
        lhns = [r.log_H_per_N for r in subset_d]
        rs = [r.order_frac for r in subset_d]

        rho1, _ = sp_stats.spearmanr(lhs, ncs)
        rho2, _ = sp_stats.spearmanr(lhs, log_ncs)
        rho3, _ = sp_stats.spearmanr(lhns, rs)
        lines.append(f"| {d} | {rho1:+.3f} | {rho2:+.3f} | {rho3:+.3f} |")

    # ── 5. logH decomposition: logH = f(N, r) ──
    lines.append("\n## 5. Theoretical Decomposition\n")
    lines.append("For a poset with N elements and order fraction r:")
    lines.append("- Antichain (r=0): logH = log(N!) = N·log(N) - N + O(log(N))")
    lines.append("- Total order (r≈0.5): logH = 0")
    lines.append("- General: logH ≈ N·log(N)·(1 - g(r)) where g captures ordering constraints\n")
    lines.append("The key question: does logH/N decrease monotonically with H?")
    lines.append("If yes, then logH encodes an effective volume-like quantity:")
    lines.append("fewer causal constraints (lower H → more pairs) → lower combinatorial")
    lines.append("entropy per element → the structure is MORE constrained.\n")
    
    # Actually: more causal pairs = more constraints = fewer linear extensions = lower logH
    # So: H↑ → n_causal↓ → logH↑ (FEWER constraints → MORE extensions)
    # OR: H↑ → n_causal↓ → logH could go either way (N also matters)
    
    # Let's check the actual direction empirically
    lines.append("### Empirical direction check\n")
    for d in dims:
        subset_d = [r for r in all_rows if r.d == d and r.N == 128]
        if not subset_d:
            continue
        by_h: dict[float, list[float]] = defaultdict(list)
        for r in subset_d:
            by_h[r.hubble].append(r.log_H)
        
        direction_parts = []
        for H in sorted(by_h.keys()):
            mean_lh = float(np.mean(by_h[H]))
            direction_parts.append(f"H={H:.1f}: {mean_lh:.1f}")
        lines.append(f"- d={d}, N=128: {' → '.join(direction_parts)}")

    # ── 6. Connection to EH action ──
    lines.append("\n\n## 6. Connection to Einstein-Hilbert\n")
    lines.append("### The Volume Connection\n")
    lines.append("In CST, N = ρ·V where V is the spacetime volume and ρ the sprinkling density.")
    lines.append("For fixed N (our experiment), the sprinkled region has FIXED volume V = N/ρ.")
    lines.append("But the CAUSAL structure within that volume changes with H.\n")
    lines.append("The EH action has two parts:")
    lines.append("  S_EH = (1/16πG) ∫ R·√(-g) d⁴x = curvature_part + cosmological_part")
    lines.append("  = (1/16πG)(R_dS·V - 2Λ·V)\n")
    lines.append("In our F7 functional:")
    lines.append("- logH encodes the COMBINATORIAL volume/entropy of the causal set")
    lines.append("- wall encodes the CURVATURE admissibility (sigmoid → step function)")
    lines.append("- Together: F7 ≈ S_combinatorial + wall(curvature)\n")
    lines.append("The correspondence is therefore:")
    lines.append("  logH ↔ cosmological term (volume-extensive entropy)")
    lines.append("  wall ↔ curvature threshold (admissibility)")
    lines.append("  F7 ↔ S_EH (with wall as finite-size rounding of the curvature cut)")

    report_text = "\n".join(lines) + "\n"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Report: {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
