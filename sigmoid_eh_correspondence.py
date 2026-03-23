"""Sigmoid → EH Correspondence Experiment.

Tests whether the F7 sigmoid wall σ((R - Rc)/w) converges to a
hard threshold Θ(S_EH - S_min) as N → ∞.

On de Sitter sprinklings with known R_dS = d(d-1)H², we:
1. Compute R (interval occupancy) for each sprinkling
2. Compute wall = σ((R - Rc)/w) with F7's parameters
3. Compute R_hat (curvature proxy from interval regression)
4. Test whether wall is a monotone function of R_hat
5. Track wall sharpness: as N grows, does R cluster → Rc becomes sharper separator?

Key predictions:
- wall(H=0) ≈ 0 or small (flat → low R → below Rc)
- wall(H>0) ≈ 1 (curved → high R → above Rc)
- The transition sharpens with N (R distribution narrows)
- Monotonicity: ρ(wall, R_hat) → +1
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from curvature_layer2_baseline import fit_curvature_proxy, fit_flat_volume_law
from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import build_causal_matrix_de_sitter


def sigmoid(x: float) -> float:
    if x > 500:
        return 1.0
    if x < -500:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


# F7 sigmoid wall parameters
ALPHA0 = 16.0
Q = -0.5
N0 = 20.0
RC = 0.25
W = 0.015


def compute_wall(R: float, N: int) -> float:
    alpha_N = ALPHA0 * (N0 / max(N, 1)) ** abs(Q)
    return alpha_N * sigmoid((R - RC) / W)


def compute_sigma_only(R: float) -> float:
    """The pure sigmoid σ((R - Rc)/w), without α(N) scaling."""
    return sigmoid((R - RC) / W)


@dataclass(frozen=True)
class Row:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    R_occupancy: float
    sigma_val: float
    wall_val: float
    alpha_N: float
    R_hat: float
    n_causal_pairs: int
    n_links: int


def count_links(causal: np.ndarray) -> int:
    """Count links (immediate predecessor relations, not mediated by others)."""
    N = causal.shape[0]
    c_int = causal.astype(np.int32)
    # interval_matrix[i,j] = number of z with i≺z≺j
    # A pair (i,j) is a link iff causal[i,j] and interval_matrix[i,j] == 0
    # For large N, compute interval sizes only for causal pairs
    ii, jj = np.where(causal)
    if len(ii) == 0:
        return 0
    
    # Batch computation
    batch_size = 5000
    n_links = 0
    c_uint8 = causal.astype(np.uint8)
    
    for start in range(0, len(ii), batch_size):
        end = min(start + batch_size, len(ii))
        bi = ii[start:end]
        bj = jj[start:end]
        rows = c_uint8[bi]
        cols = c_uint8[:, bj]
        k = np.einsum('pz,zp->p', rows, cols)
        n_links += int(np.sum(k == 0))
    
    return n_links


def run_single(d: int, N: int, hubble: float, rep: int, seed: int,
               tau_quantile: float) -> Row:
    """Run one realization."""
    s = seed + d * 100000 + N * 100 + rep + int(hubble * 1000)
    points = sprinkle_de_sitter_like_diamond(N, d - 1, hubble=hubble, seed=s)
    causal = build_causal_matrix_de_sitter(points, hubble)

    n_pairs = int(causal.sum())
    
    # Count links and compute occupancy
    n_lnk = count_links(causal)
    R_occ = 1.0 - n_lnk / max(n_pairs, 1) if n_pairs > 0 else 0.0

    # Sigmoid and wall
    sigma_val = compute_sigma_only(R_occ)
    wall_val = compute_wall(R_occ, N)
    alpha_N = ALPHA0 * (N0 / max(N, 1)) ** abs(Q)

    # Also compute R_hat for comparison
    t = points[:, 0]
    spatial = points[:, 1:]
    dt = t[None, :] - t[:, None]
    if hubble == 0.0:
        tau_matrix = np.sqrt(np.clip(dt * dt - np.sum(
            (spatial[None, :, :] - spatial[:, None, :]) ** 2, axis=2
        ), 0.0, None))
    else:
        ti = t[:, None]
        tj = t[None, :]
        chi = (np.exp(-hubble * ti) - np.exp(-hubble * tj)) / hubble
        tau_matrix = np.clip(chi, 0.0, None)

    c_int = causal.astype(np.int32)
    
    # For R_hat, use sampled pairs if N > 2048
    if N <= 2048:
        interval_matrix = c_int @ c_int
        tau_arr = tau_matrix[causal]
        k_arr = interval_matrix[causal].astype(float)
    else:
        ii, jj = np.where(causal)
        rng = np.random.RandomState(s + 999)
        max_pairs = 30000
        if len(ii) > max_pairs:
            idx = rng.choice(len(ii), max_pairs, replace=False)
            ii_s, jj_s = ii[idx], jj[idx]
        else:
            ii_s, jj_s = ii, jj
        
        tau_arr = tau_matrix[ii_s, jj_s]
        c_uint8 = causal.astype(np.uint8)
        k_arr = np.zeros(len(ii_s), dtype=float)
        bs = 2000
        for start in range(0, len(ii_s), bs):
            end = min(start + bs, len(ii_s))
            rows = c_uint8[ii_s[start:end]]
            cols = c_uint8[:, jj_s[start:end]]
            k_arr[start:end] = np.einsum('pz,zp->p', rows, cols).astype(float)

    if tau_arr.size < 10:
        r_hat = float('nan')
    else:
        y = k_arr / max(N - 2, 1)
        x = np.power(tau_arr, d)
        ratio_arr = y / np.clip(x, 1e-12, None)
        fit_mask = tau_arr <= np.quantile(tau_arr, tau_quantile)
        if fit_mask.sum() < 10:
            fit_mask = np.ones_like(tau_arr, dtype=bool)
        _, r_hat = fit_curvature_proxy(tau_arr[fit_mask], ratio_arr[fit_mask], d)

    R_dS = d * (d - 1) * hubble ** 2

    return Row(
        d=d, N=N, hubble=hubble, rep=rep, R_dS=R_dS,
        R_occupancy=R_occ, sigma_val=sigma_val, wall_val=wall_val,
        alpha_N=alpha_N, R_hat=r_hat,
        n_causal_pairs=n_pairs, n_links=n_lnk,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="Sigmoid → EH Correspondence")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[64, 128, 256, 512, 1024])
    ap.add_argument("--hubbles", nargs="*", type=float, default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--tau-quantile", type=float, default=0.75)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default="outputs_unified_functional/sigmoid_eh_correspondence.csv")
    ap.add_argument("--report", default="outputs_unified_functional/sigmoid_eh_correspondence.md")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[Row] = []
    total = len(args.dims) * len(args.ns) * len(args.hubbles) * args.reps
    done = 0

    for d in args.dims:
        for N in args.ns:
            for hubble in args.hubbles:
                for rep in range(args.reps):
                    row = run_single(d, N, hubble, rep, args.seed, args.tau_quantile)
                    rows.append(row)
                    done += 1
                    if done % 10 == 0 or done == total:
                        print(f"  [{done}/{total}] d={d} N={N} H={hubble:.1f} "
                              f"R={row.R_occupancy:.3f} σ={row.sigma_val:.3f} "
                              f"wall={row.wall_val:.2f} R_hat={row.R_hat:+.1f}",
                              flush=True)

    # Save CSV
    fieldnames = [f.name for f in fields(Row)]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({fn: getattr(row, fn) for fn in fieldnames})
    print(f"\nSaved: {out_path}")

    # ── Report ──
    report_path = Path(args.report)
    lines: list[str] = []
    lines.append("# Sigmoid → EH Correspondence\n")
    lines.append(f"Total: {len(rows)} realizations\n")
    lines.append(f"F7 wall: α(N)·σ((R−Rc)/w), α₀={ALPHA0}, q={Q}, Rc={RC}, w={W}\n")

    ds = sorted(set(r.d for r in rows))
    ns = sorted(set(r.N for r in rows))
    hs = sorted(set(r.hubble for r in rows))

    # ── 1. R(occupancy) vs H ──
    lines.append("\n## 1. Occupancy R vs Hubble Parameter\n")
    lines.append("Does R increase with curvature?\n")
    lines.append("| d | N | H=0 | H=0.25 | H=0.5 | H=1.0 | H=2.0 |")
    lines.append("|---|---|-----|--------|-------|-------|-------|")

    for d in ds:
        for N in ns:
            vals = []
            for H in hs:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == H]
                if subset:
                    mean_R = float(np.mean([r.R_occupancy for r in subset]))
                    vals.append(f"{mean_R:.3f}")
                else:
                    vals.append("—")
            lines.append(f"| {d} | {N} | " + " | ".join(vals) + " |")

    # ── 2. σ((R-Rc)/w) vs H — the pure sigmoid ──
    lines.append("\n## 2. Pure Sigmoid σ((R−Rc)/w) vs H\n")
    lines.append(f"Threshold: Rc={RC}, width: w={W}\n")
    lines.append("| d | N | H=0 | H=0.25 | H=0.5 | H=1.0 | H=2.0 |")
    lines.append("|---|---|-----|--------|-------|-------|-------|")

    for d in ds:
        for N in ns:
            vals = []
            for H in hs:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == H]
                if subset:
                    mean_s = float(np.mean([r.sigma_val for r in subset]))
                    vals.append(f"{mean_s:.4f}")
                else:
                    vals.append("—")
            lines.append(f"| {d} | {N} | " + " | ".join(vals) + " |")

    # ── 3. Wall sharpening: std(R) per (d, N, H) ──
    lines.append("\n## 3. Occupancy Variance Shrinking (Wall Sharpening)\n")
    lines.append("If std(R|d,N,H) → 0, the sigmoid becomes a step function.\n")
    lines.append("| d | N | std(R|H=0) | std(R|H=0.5) | std(R|H=1.0) | std(R|H=2.0) |")
    lines.append("|---|---|-----------|-------------|-------------|-------------|")

    for d in ds:
        for N in ns:
            vals = []
            for H in [0.0, 0.5, 1.0, 2.0]:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == H]
                if subset:
                    std_R = float(np.std([r.R_occupancy for r in subset]))
                    vals.append(f"{std_R:.4f}")
                else:
                    vals.append("—")
            lines.append(f"| {d} | {N} | " + " | ".join(vals) + " |")

    # ── 4. Monotonicity: σ vs R_hat ──
    lines.append("\n## 4. Monotonicity: σ((R−Rc)/w) vs R_hat\n")
    lines.append("| d | N | Spearman(σ, R_hat) | p-value |")
    lines.append("|---|---|-------------------|---------|")

    for d in ds:
        for N in ns:
            subset = [r for r in rows if r.d == d and r.N == N
                      and not math.isnan(r.R_hat) and r.hubble > 0]
            if len(subset) < 5:
                continue
            s_arr = np.array([r.sigma_val for r in subset])
            rh_arr = np.array([r.R_hat for r in subset])
            if np.std(s_arr) < 1e-12:
                lines.append(f"| {d} | {N} | σ constant | — |")
                continue
            rho, p = sp_stats.spearmanr(s_arr, rh_arr)
            lines.append(f"| {d} | {N} | {rho:+.3f} | {p:.2e} |")

    # ── 5. Separation: does σ separate H=0 from H>0? ──
    lines.append("\n## 5. Hard Threshold Test: σ(H=0) vs σ(H>0)\n")
    lines.append("If sigmoid → Θ, then σ(H=0) → 0 and σ(H>0) → 1 for large N.\n")
    lines.append("| d | N | mean σ(H=0) | mean σ(H>0) | gap | AUC |")
    lines.append("|---|---|-------------|-------------|-----|-----|")

    for d in ds:
        for N in ns:
            s0 = [r.sigma_val for r in rows if r.d == d and r.N == N and r.hubble == 0]
            s1 = [r.sigma_val for r in rows if r.d == d and r.N == N and r.hubble > 0]
            if not s0 or not s1:
                continue
            m0 = float(np.mean(s0))
            m1 = float(np.mean(s1))
            gap = m1 - m0

            # Simple AUC: fraction of (s1_i > s0_j) pairs
            n_correct = sum(1 for a in s1 for b in s0 if a > b)
            n_total_pairs = len(s1) * len(s0)
            auc = n_correct / max(n_total_pairs, 1)

            lines.append(f"| {d} | {N} | {m0:.4f} | {m1:.4f} | {gap:+.4f} | {auc:.3f} |")

    # ── 6. N-scaling of separation ──
    lines.append("\n## 6. N-Scaling of the Sigmoid-to-Step Convergence\n")
    lines.append("Does the gap σ(H>0) - σ(H=0) increase with N?\n")

    for d in ds:
        d_ns = sorted([N for N in ns])
        gaps = []
        for N in d_ns:
            s0 = [r.sigma_val for r in rows if r.d == d and r.N == N and r.hubble == 0]
            s1 = [r.sigma_val for r in rows if r.d == d and r.N == N and r.hubble > 0]
            if s0 and s1:
                gaps.append(float(np.mean(s1)) - float(np.mean(s0)))
            else:
                gaps.append(float('nan'))

        valid = [(n, g) for n, g in zip(d_ns, gaps) if not math.isnan(g)]
        if len(valid) >= 3:
            ns_v, gs_v = zip(*valid)
            rho, p = sp_stats.spearmanr(ns_v, gs_v)
            lines.append(f"\n### d={d}\n")
            lines.append(f"Gap trajectory: " + " → ".join(f"{g:.4f}" for g in gaps))
            lines.append(f"- Spearman(N, gap): ρ={rho:+.3f} (p={p:.3e})")
            if rho > 0.7:
                lines.append(f"- → **Gap increasing with N** — sigmoid sharpening toward step")

    # ── Summary ──
    lines.append("\n## Summary\n")
    lines.append("The sigmoid-to-step convergence requires:\n")
    lines.append("1. **R increases with curvature**: R(H>0) > R(H=0)")
    lines.append("2. **σ separates flat from curved**: σ(H>0) ≈ 1, σ(H=0) ≈ 0")
    lines.append("3. **Separation sharpens with N**: gap and AUC increase")
    lines.append("4. **σ is monotone with R_hat**: Spearman(σ, R_hat) → +1")
    lines.append("")
    lines.append("If all four hold, then σ((R−Rc)/w) → Θ(R_dS − R_dS,min) = Θ(S_EH − S_min).")

    report_text = "\n".join(lines) + "\n"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved: {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
