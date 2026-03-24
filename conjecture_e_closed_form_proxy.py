"""Conjecture E — Closed-Form Curvature Proxy Coefficients (§4.1.26b).

Following §4.1.26's PC-space proxy design, this script extracts the
**per-N closed-form C_k coefficients** for d=3 and d=4 optimal proxies,
analyzes their stability across N, and tests whether a simplified
2–3 coefficient formula can match full PC-optimal performance.

CRITICAL: bdg_d2c = N − 2C₀ + 2C₁ operates on RAW C_k counts, NOT on
normalized p_k = C_k/ΣC_k.  The sign flips under normalization because
ΣC_k itself anti-correlates with H.  This script therefore optimizes
in BOTH spaces and uses bdg_d2c values from the CSV (raw-C_k-based)
for fair comparison.

Key questions:
  Q1. What are the optimal C_k weights at each (d, N)?
  Q2. Do these weights converge with N?  (stability → physical meaning)
  Q3. How do they compare structurally with bdg_d2c's [-2, +2, 0, 0, 0]?
  Q4. Can a sparse approximation (2–3 nonzero C_k) match the full proxy?
  Q5. What is the physical interpretation of the dominant coefficients?

Data: loads §4.1.24 CSV (360 realizations with C0–C4).
"""

from __future__ import annotations

import argparse
import csv
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats


# ──────────────────────────────────────────────────────────────────
# Utility
# ──────────────────────────────────────────────────────────────────

def load_data(csv_path: str) -> list[dict]:
    rows = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            d = {}
            for k, v in row.items():
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
            rows.append(d)
    return rows


def build_matrices(rows: list[dict], max_k: int = 4):
    """Return (P, C_raw) where P is normalized p_k and C_raw is raw counts."""
    n = len(rows)
    C = np.zeros((n, max_k + 1))
    for i, r in enumerate(rows):
        for k in range(max_k + 1):
            C[i, k] = float(r.get(f"C{k}", 0))
    row_sums = C.sum(axis=1, keepdims=True)
    row_sums = np.clip(row_sums, 1.0, None)
    P = C / row_sums
    return P, C


def eval_proxy(X, coef, h2):
    """Return (Spearman ρ, p-value) for proxy = X @ coef vs h2."""
    proxy = X @ coef
    if np.std(proxy) < 1e-15:
        return 0.0, 1.0
    return sp_stats.spearmanr(h2, proxy)


def optimize_direction(X, h2, n_restarts=100, seed=42):
    """Find unit vector w maximizing |Spearman(X@w, H²)| via
    random-restart coordinate descent."""
    n_feat = X.shape[1]
    if n_feat == 0:
        return np.array([]), 0.0

    best_rho = 0.0
    best_w = np.zeros(n_feat)
    rng = np.random.RandomState(seed)

    for _ in range(n_restarts):
        w = rng.randn(n_feat)
        w /= np.linalg.norm(w)

        for _ in range(40):
            rho, _ = sp_stats.spearmanr(h2, X @ w)
            if abs(rho) > best_rho:
                best_rho = abs(rho)
                best_w = w.copy() * np.sign(rho)

            improved = False
            for dim in range(n_feat):
                for delta in [0.20, -0.20, 0.08, -0.08, 0.03, -0.03]:
                    w2 = w.copy()
                    w2[dim] += delta
                    nrm = np.linalg.norm(w2)
                    if nrm < 1e-10:
                        continue
                    w2 /= nrm
                    r2, _ = sp_stats.spearmanr(h2, X @ w2)
                    if abs(r2) > abs(rho):
                        w, rho = w2, r2
                        improved = True
                        if abs(r2) > best_rho:
                            best_rho = abs(r2)
                            best_w = w2.copy() * np.sign(r2)
            if not improved:
                break

    return best_w, best_rho


def normalize_coef(w):
    m = np.max(np.abs(w))
    return w / m if m > 0 else w


def sparse_search(X, h2, max_nz, n_restarts=40):
    """Best proxy with ≤ max_nz nonzero coefficients."""
    n_feat = X.shape[1]
    best_rho, best_w, best_idx = 0.0, np.zeros(n_feat), ()
    for idx in combinations(range(n_feat), max_nz):
        w_sel, rho_sel = optimize_direction(
            X[:, list(idx)], h2, n_restarts=n_restarts, seed=77)
        if abs(rho_sel) > abs(best_rho):
            best_rho = rho_sel
            full_w = np.zeros(n_feat)
            for j, k in enumerate(idx):
                full_w[k] = w_sel[j]
            best_w, best_idx = full_w, idx
    return best_w, best_rho, best_idx


# ──────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="outputs_unified_functional/conjecture_e_shape_ratio_family.csv")
    ap.add_argument("--report",
                    default="outputs_unified_functional/conjecture_e_closed_form_proxy.md")
    ap.add_argument("--max-k", type=int, default=4)
    args = ap.parse_args()

    rows = load_data(args.csv)
    print(f"Loaded {len(rows)} rows")

    P, C_raw = build_matrices(rows, args.max_k)
    n_feat = P.shape[1]

    dims    = np.array([r["d"] for r in rows])
    ns      = np.array([r["N"] for r in rows])
    hubbles = np.array([r["hubble"] for r in rows])
    h2      = hubbles ** 2
    bdg_csv = np.array([r.get("S8_bdg_d2c_shape", 0.0) for r in rows])
    kl_csv  = np.array([r.get("S10_kl_vs_flat", float("nan")) for r in rows])

    unique_ns = sorted(set(ns))
    L: list[str] = []          # report lines

    L.append("# Conjecture E — Closed-Form Curvature Proxy (§4.1.26b)\n")
    L.append(f"Data: {len(rows)} realizations, C_0 … C_{args.max_k}\n")

    # ── bdg_d2c reference ────────────────────────────────────────
    bdg_Ck = np.array([-2.0, +2.0, 0.0, 0.0, 0.0])  # in raw C_k
    L.append("## Reference\n")
    L.append("bdg_d2c = N − 2C₀ + 2C₁  →  raw-C_k coef = [-2, +2, 0, 0, 0]")
    L.append("NOTE: bdg_d2c is evaluated on raw C_k (from CSV), NOT on p_k.\n")

    # ═══════════════════════════════════════════════════════════════
    # Per-dimension analysis:  d = 3  and  d = 4
    # ═══════════════════════════════════════════════════════════════
    all_res: dict[tuple[int,int], dict] = {}

    for d_val in (3.0, 4.0):
        dm = dims == d_val
        C_d, P_d = C_raw[dm], P[dm]
        h2_d, ns_d = h2[dm], ns[dm]
        bdg_d, kl_d = bdg_csv[dm], kl_csv[dm]
        di = int(d_val)

        L.append(f"\n{'='*60}")
        L.append(f"# d = {di}")
        L.append(f"{'='*60}\n")

        # ── A. Full 5-coef optimization in raw C_k space ─────────
        L.append("## A. Full 5-coefficient optimal proxy (raw C_k space)\n")
        L.append("| N | w₀ | w₁ | w₂ | w₃ | w₄ | |ρ| | bdg ρ | KL ρ |")
        L.append("|---|----|----|----|----|----|----|-------|------|")

        for Nv in unique_ns:
            nm = ns_d == Nv
            if nm.sum() < 5:
                continue
            C_s, h2_s = C_d[nm], h2_d[nm]
            bdg_s, kl_s = bdg_d[nm], kl_d[nm]

            w_opt, rho_opt = optimize_direction(C_s, h2_s, n_restarts=120)
            wn = normalize_coef(w_opt)

            # bdg_d2c from CSV (raw)
            rho_bdg, _ = sp_stats.spearmanr(h2_s, bdg_s)

            # KL from CSV
            valid = ~np.isnan(kl_s)
            rho_kl = sp_stats.spearmanr(h2_s[valid], kl_s[valid])[0] if valid.sum() >= 5 else float("nan")

            all_res[(di, int(Nv))] = dict(
                w_Ck=wn.copy(), rho_Ck=rho_opt,
                rho_bdg=rho_bdg, rho_kl=rho_kl)

            L.append(f"| {int(Nv)} | "
                     f"{wn[0]:+.3f} | {wn[1]:+.3f} | {wn[2]:+.3f} | "
                     f"{wn[3]:+.3f} | {wn[4]:+.3f} | "
                     f"{rho_opt:.3f} | {rho_bdg:+.3f} | {rho_kl:+.3f} |")

        # ── B. Full 5-coef optimization in normalized p_k space ──
        L.append(f"\n## B. Full 5-coefficient optimal proxy (normalized p_k space)\n")
        L.append("| N | w₀ | w₁ | w₂ | w₃ | w₄ | |ρ| |")
        L.append("|---|----|----|----|----|----|-----|")

        for Nv in unique_ns:
            nm = ns_d == Nv
            if nm.sum() < 5:
                continue
            P_s, h2_s = P_d[nm], h2_d[nm]

            w_pk, rho_pk = optimize_direction(P_s, h2_s, n_restarts=120)
            wn_pk = normalize_coef(w_pk)

            key = (di, int(Nv))
            all_res[key]["w_pk"] = wn_pk.copy()
            all_res[key]["rho_pk"] = rho_pk

            L.append(f"| {int(Nv)} | "
                     f"{wn_pk[0]:+.3f} | {wn_pk[1]:+.3f} | {wn_pk[2]:+.3f} | "
                     f"{wn_pk[3]:+.3f} | {wn_pk[4]:+.3f} | {rho_pk:.3f} |")

        # ── C. Sparse search in raw C_k space ────────────────────
        L.append(f"\n## C. Sparse proxy search (raw C_k)\n")
        for max_nz in (2, 3):
            L.append(f"\n### Best {max_nz}-coefficient proxy\n")
            L.append("| N | Indices | Coefs (norm) | |ρ| sparse | |ρ| full | bdg ρ |")
            L.append("|---|---------|-------------|-----------|---------|-------|")

            for Nv in unique_ns:
                nm = ns_d == Nv
                if nm.sum() < 5:
                    continue
                C_s, h2_s = C_d[nm], h2_d[nm]

                sw, sr, si = sparse_search(C_s, h2_s, max_nz)
                sn = normalize_coef(sw)
                key = (di, int(Nv))
                full_rho = all_res[key]["rho_Ck"]
                bdg_rho = all_res[key]["rho_bdg"]

                idx_s = ", ".join(f"C{k}" for k in si)
                coef_s = ", ".join(f"{sn[k]:+.3f}" for k in si)

                tag = ""
                if key not in all_res:
                    tag = ""
                all_res[key][f"sparse{max_nz}_idx"] = si
                all_res[key][f"sparse{max_nz}_rho"] = sr
                all_res[key][f"sparse{max_nz}_w"] = sn.copy()

                L.append(f"| {int(Nv)} | {idx_s} | [{coef_s}] | "
                         f"{abs(sr):.3f} | {full_rho:.3f} | {bdg_rho:+.3f} |")

    # ═══════════════════════════════════════════════════════════════
    # Coefficient convergence
    # ═══════════════════════════════════════════════════════════════
    L.append(f"\n{'='*60}")
    L.append("# Coefficient Convergence")
    L.append(f"{'='*60}\n")

    for di in (3, 4):
        L.append(f"\n## d = {di}\n")

        # Cosine similarity across N for raw-C_k optimal
        coefs_N = [(Nv, all_res[(di,Nv)]["w_Ck"])
                   for Nv in [128,256,512] if (di,Nv) in all_res]
        for i in range(len(coefs_N)-1):
            N1, w1 = coefs_N[i]
            N2, w2 = coefs_N[i+1]
            cos = np.dot(w1, w2) / (np.linalg.norm(w1)*np.linalg.norm(w2)+1e-15)
            L.append(f"cos(w@N={N1}, w@N={N2}) = {cos:.3f}")

        # Cosine with bdg_d2c (in C_k space)
        bdg_n = normalize_coef(bdg_Ck)
        for Nv in [128,256,512]:
            if (di,Nv) not in all_res:
                continue
            w = all_res[(di,Nv)]["w_Ck"]
            cos = np.dot(w, bdg_n) / (np.linalg.norm(w)*np.linalg.norm(bdg_n)+1e-15)
            L.append(f"cos(opt@N={Nv}, bdg_d2c) = {cos:.3f}")

    # ═══════════════════════════════════════════════════════════════
    # Interpretation
    # ═══════════════════════════════════════════════════════════════
    L.append(f"\n{'='*60}")
    L.append("# Interpretation")
    L.append(f"{'='*60}\n")
    L.append("""### C_k meaning (raw counts):
- C₀ = number of links (order-1 intervals)
- C₁ = number of order-2 intervals
- C₂ = number of order-3 intervals
- C₃ = number of order-4 intervals
- C₄ = number of order-5 intervals

### bdg_d2c = N − 2C₀ + 2C₁  (raw C_k, NOT normalized):
- Positive correlation with H² at d≥3 (from CSV: ρ ≈ +0.82 to +0.93)
- Operates only on (C₀, C₁) — ignores higher intervals

### Key questions answered by this analysis:
Q3: Does the optimal proxy structurally resemble bdg_d2c?
Q4: Can a 2-coefficient formula match the full 5-coefficient optimal?
Q5: Do higher-order C_k carry additional curvature information?
""")

    # Save
    rp = Path(args.report)
    rp.parent.mkdir(parents=True, exist_ok=True)
    rp.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\nSaved: {rp}")

    # Terminal summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for line in L:
        print(line)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
