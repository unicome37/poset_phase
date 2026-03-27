"""
Expanded Family Robustness Test — 25 Families
===============================================

Test Prediction B robustness by adding 8 adversarial families to the
original 17, for a total of 25.  The core claim: Lor4D uniquely minimises
both LSD-Well and Mahalanobis LSD at every tested N.

New families (constructed from existing primitives):
  KR_8layer        — 8-layer fully-bipartite layered poset
  RandomDAG        — Erdős–Rényi DAG, p=0.3
  RandomDAG_sp     — Erdős–Rényi DAG, p=0.1 (sparse)
  RandomDAG_dn     — Erdős–Rényi DAG, p=0.6 (dense)
  ChainAnti4       — 4 antichains connected sequentially
  ChainAnti12      — 12 antichains connected sequentially (deep)
  SparseChain      — long chain + random side branches
  MixedLor4D       — Lor4D with 20% random relation corruption

Metrics:
  LSD-Well:  0.5·(d-4)² + 1.0·(c-c*(N))² + 5.0·(w-w*(N))²
  Mahalanobis: (I-μ)ᵀ Σ⁻¹ (I-μ),  μ and Σ from Lor4D samples
"""
from __future__ import annotations

import time
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from bd_action import count_intervals_fast
from generators import (
    Poset,
    transitive_closure,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_kr_like,
    generate_kr_2layer,
    generate_kr_4layer,
    generate_absolute_layered,
    generate_multi_layer_random,
    generate_random_layered_k4_uniform,
    generate_random_layered_k6_uniform,
    generate_random_layered_k8_uniform,
    generate_random_layered_k6_tapered,
    generate_random_layered_k6_middle_heavy,
    generate_random_layered_k6_longjump,
    generate_transitive_percolation,
    generate_interval_order,
)
from unified_functional import compute_xi_dim


# ══════════════════════════════════════════════════════════════════════════
#  New adversarial generators
# ══════════════════════════════════════════════════════════════════════════

def generate_deep_layered(n: int, n_layers: int = 8, seed: int | None = None) -> Poset:
    """Deep layered poset: n_layers groups with full bipartite edges between
    adjacent layers."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    sizes = [n // n_layers] * n_layers
    for i in range(n % n_layers):
        sizes[i] += 1
    adj = np.zeros((n, n), dtype=bool)
    offset = 0
    layer_starts = []
    for s in sizes:
        layer_starts.append(offset)
        offset += s
    for L in range(n_layers - 1):
        s_cur = sizes[L]
        s_nxt = sizes[L + 1]
        for i in perm[layer_starts[L]:layer_starts[L] + s_cur]:
            for j in perm[layer_starts[L + 1]:layer_starts[L + 1] + s_nxt]:
                adj[i, j] = True
    return Poset(transitive_closure(adj))


def generate_random_dag(n: int, p: float = 0.3, seed: int | None = None) -> Poset:
    """Erdős–Rényi random DAG on a fixed topological order."""
    rng = np.random.default_rng(seed)
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p:
                adj[i, j] = True
    return Poset(transitive_closure(adj))


def generate_chain_of_antichains(n: int, k: int = 4, seed: int | None = None) -> Poset:
    """k antichains of ~n/k elements, connected sequentially (full bipartite
    between adjacent antichains)."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    sizes = [n // k] * k
    for i in range(n % k):
        sizes[i] += 1
    adj = np.zeros((n, n), dtype=bool)
    offset = 0
    for layer_idx in range(k - 1):
        next_start = offset + sizes[layer_idx]
        next_end = next_start + sizes[layer_idx + 1]
        for i in perm[offset:offset + sizes[layer_idx]]:
            for j in perm[next_start:next_end]:
                adj[i, j] = True
        offset += sizes[layer_idx]
    return Poset(transitive_closure(adj))


def generate_sparse_chain(n: int, branch_prob: float = 0.2, seed: int | None = None) -> Poset:
    """Long chain (n//2 nodes) with remaining nodes attached as random
    side-branches."""
    rng = np.random.default_rng(seed)
    adj = np.zeros((n, n), dtype=bool)
    chain_len = n // 2
    for i in range(chain_len - 1):
        adj[i, i + 1] = True
    for i in range(chain_len, n):
        parent = rng.integers(0, chain_len)
        adj[parent, i] = True
        if parent + 1 < chain_len:
            adj[i, parent + 1] = True
    return Poset(transitive_closure(adj))


def generate_mixed_lor4d(n: int, corruption: float = 0.2, seed: int | None = None) -> Poset:
    """Lor4D with ~corruption fraction of random relation changes."""
    p = generate_lorentzian_like_4d(n, seed=seed)
    rng = np.random.default_rng(seed + 999 if seed is not None else None)
    adj = p.closure.copy()
    n_changes = int(corruption * n * (n - 1) / 2)
    for _ in range(n_changes):
        pair = rng.choice(n, 2, replace=False)
        i, j = int(min(pair)), int(max(pair))
        adj[i, j] = not adj[i, j]
    # Enforce upper-triangular (DAG with natural order)
    for i in range(n):
        for j in range(i):
            adj[i, j] = False
    return Poset(transitive_closure(adj))


# ══════════════════════════════════════════════════════════════════════════
#  Family registry
# ══════════════════════════════════════════════════════════════════════════

ORIGINAL_FAMILIES = {
    "Lor2D":     generate_lorentzian_like_2d,
    "Lor3D":     generate_lorentzian_like_3d,
    "Lor4D":     generate_lorentzian_like_4d,
    "Lor5D":     generate_lorentzian_like_5d,
    "KR_like":   generate_kr_like,
    "KR_2layer": generate_kr_2layer,
    "KR_4layer": generate_kr_4layer,
    "AbsLayer":  generate_absolute_layered,
    "MLR":       generate_multi_layer_random,
    "RLk4":      generate_random_layered_k4_uniform,
    "RLk6":      generate_random_layered_k6_uniform,
    "RLk8":      generate_random_layered_k8_uniform,
    "RLk6_tap":  generate_random_layered_k6_tapered,
    "RLk6_mid":  generate_random_layered_k6_middle_heavy,
    "RLk6_lj":   generate_random_layered_k6_longjump,
    "TransPerc": generate_transitive_percolation,
    "IntOrder":  generate_interval_order,
}

NEW_FAMILIES = {
    "KR_8layer":     lambda n, seed=None: generate_deep_layered(n, n_layers=8, seed=seed),
    "RandomDAG":     lambda n, seed=None: generate_random_dag(n, p=0.3, seed=seed),
    "RandomDAG_sp":  lambda n, seed=None: generate_random_dag(n, p=0.1, seed=seed),
    "RandomDAG_dn":  lambda n, seed=None: generate_random_dag(n, p=0.6, seed=seed),
    "ChainAnti4":    lambda n, seed=None: generate_chain_of_antichains(n, k=4, seed=seed),
    "ChainAnti12":   lambda n, seed=None: generate_chain_of_antichains(n, k=12, seed=seed),
    "SparseChain":   lambda n, seed=None: generate_sparse_chain(n, seed=seed),
    "MixedLor4D":    lambda n, seed=None: generate_mixed_lor4d(n, seed=seed),
}

ALL_FAMILIES = {**ORIGINAL_FAMILIES, **NEW_FAMILIES}


# ══════════════════════════════════════════════════════════════════════════
#  Feature extraction & scoring
# ══════════════════════════════════════════════════════════════════════════

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


def compute_features(poset: Poset, N: int) -> np.ndarray:
    counts = count_intervals_fast(poset, k_max=5)
    C0 = counts.get(0)
    C1 = counts.get(1)
    c1_c0 = C1 / max(1, C0)
    _, d_eff = compute_xi_dim(poset)
    aw = max_antichain_width(poset)
    width_ratio = aw / max(1, N)
    return np.array([d_eff, c1_c0, width_ratio])


def cstar(N: int) -> float:
    return 0.2485 - 2.33 / N


def wstar(N: int) -> float:
    return 0.3255 + 3.80 / N


def lsd_well_score(feat: np.ndarray, N: int) -> float:
    d_eff, c1_c0, width_ratio = feat
    return (0.5 * (d_eff - 4.0) ** 2
            + 1.0 * (c1_c0 - cstar(N)) ** 2
            + 5.0 * (width_ratio - wstar(N)) ** 2)


def mahalanobis_score(feat: np.ndarray, mu: np.ndarray, cov_inv: np.ndarray) -> float:
    delta = feat - mu
    return float(delta @ cov_inv @ delta)


# ══════════════════════════════════════════════════════════════════════════
#  Main experiment
# ══════════════════════════════════════════════════════════════════════════

def main() -> None:
    # Avoid Windows console encoding crashes when printing non-ASCII symbols.
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    N_VALUES = [16, 28, 48, 64, 96]
    REPS = 20
    SEED_BASE = 42
    # High-precision Lor4D reference ensemble at small N (no tunable model knob;
    # this is purely estimation precision for the reference manifold).
    LOR4D_REF_REPS_N16 = 80

    n_fam = len(ALL_FAMILIES)
    total_samples = n_fam * len(N_VALUES) * REPS
    print("=" * 80)
    print("EXPANDED FAMILY ROBUSTNESS TEST — 25 Families")
    print(f"  Original 17 + 8 new adversarial = {n_fam} families")
    print(f"  N values: {N_VALUES}")
    print(f"  Reps per (family, N): {REPS}")
    print(f"  Total samples: {total_samples}")
    print("=" * 80)

    t0 = time.time()

    # Collect evaluation features: eval_data[N][family_name] = list of feature vectors
    # Collect reference features for Lor4D: ref_lor4d[N] = list of feature vectors
    eval_data: dict[int, dict[str, list[np.ndarray]]] = defaultdict(lambda: defaultdict(list))
    ref_lor4d: dict[int, list[np.ndarray]] = defaultdict(list)

    for fi, (fam_name, gen_fn) in enumerate(ALL_FAMILIES.items()):
        print(f"  [{fi+1:2d}/{n_fam}] {fam_name:15s}", end="", flush=True)
        for N in N_VALUES:
            ok = 0
            for rep in range(REPS):
                seed = SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep
                seed = seed % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    feat = compute_features(poset, N)
                    eval_data[N][fam_name].append(feat)
                    if fam_name == "Lor4D":
                        ref_lor4d[N].append(feat)
                    ok += 1
                except Exception:
                    pass
            # Extra Lor4D reference samples at N=16 only (separate seeds)
            if fam_name == "Lor4D" and N == 16 and LOR4D_REF_REPS_N16 > 0:
                ok_ref = 0
                for rep in range(LOR4D_REF_REPS_N16):
                    seed = SEED_BASE + 999999 + N * 1000 + rep
                    seed = seed % (2**31)
                    try:
                        poset = gen_fn(N, seed=seed)
                        feat = compute_features(poset, N)
                        ref_lor4d[N].append(feat)
                        ok_ref += 1
                    except Exception:
                        pass
                print(f"+ref:{ok_ref}", end="", flush=True)
            print(f"  N={N}:{ok}", end="", flush=True)
        print()

    elapsed_gen = time.time() - t0
    print(f"\nGeneration done in {elapsed_gen:.1f}s")

    # ── Scoring ──────────────────────────────────────────────────────────
    # For each N: Lor4D statistics, then score every family
    results_lsd: dict[int, list[tuple[str, float]]] = {}
    results_mahal: dict[int, list[tuple[str, float]]] = {}

    for N in N_VALUES:
        # Use the high-precision reference ensemble if present; otherwise fall back
        # to the Lor4D evaluation set.
        lor4d_ref = np.array(ref_lor4d.get(N, []))
        lor4d_feats = lor4d_ref if len(lor4d_ref) >= 5 else np.array(eval_data[N]["Lor4D"])
        if len(lor4d_feats) == 0:
            print(f"  WARNING: no Lor4D data at N={N}")
            continue
        mu = np.mean(lor4d_feats, axis=0)
        cov = np.cov(lor4d_feats.T)
        cov_inv = np.linalg.inv(cov + 1e-12 * np.eye(3))

        fam_lsd_scores = []
        fam_mahal_scores = []
        for fam_name in ALL_FAMILIES:
            feats = eval_data[N][fam_name]
            if not feats:
                continue
            feats_arr = np.array(feats)
            avg_lsd = np.mean([lsd_well_score(f, N) for f in feats_arr])
            avg_mahal = np.mean([mahalanobis_score(f, mu, cov_inv) for f in feats_arr])
            fam_lsd_scores.append((fam_name, avg_lsd))
            fam_mahal_scores.append((fam_name, avg_mahal))

        results_lsd[N] = sorted(fam_lsd_scores, key=lambda x: x[1])
        results_mahal[N] = sorted(fam_mahal_scores, key=lambda x: x[1])

    # ── Report ───────────────────────────────────────────────────────────
    report: list[str] = []
    report.append("# Expanded Family Robustness Test — 25 Families\n")
    report.append(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    report.append(f"Families: 17 original + 8 adversarial = {n_fam}")
    report.append(f"N values: {N_VALUES}")
    report.append(f"Reps: {REPS}\n")
    report.append("Reference manifold precision:\n")
    report.append(f"- Lor4D reference reps at N=16: {REPS} + {LOR4D_REF_REPS_N16} (extra ref)")
    report.append(f"- Lor4D reference reps at other N: {REPS}\n")

    new_family_names = set(NEW_FAMILIES.keys())

    # 1. Lor4D rank at each N
    report.append("## 1. Lor4D Rank Summary\n")
    report.append("| N | LSD-Well rank | LSD-Well winner | Mahalanobis rank | Mahalanobis winner |")
    report.append("|---|---|---|---|---|")
    lor4d_lsd_ranks = []
    lor4d_mahal_ranks = []
    for N in N_VALUES:
        if N not in results_lsd:
            continue
        lsd_ranked = results_lsd[N]
        mahal_ranked = results_mahal[N]
        lsd_rank = next(i + 1 for i, (f, _) in enumerate(lsd_ranked) if f == "Lor4D")
        mahal_rank = next(i + 1 for i, (f, _) in enumerate(mahal_ranked) if f == "Lor4D")
        lor4d_lsd_ranks.append(lsd_rank)
        lor4d_mahal_ranks.append(mahal_rank)
        report.append(f"| {N} | {lsd_rank} | {lsd_ranked[0][0]} | {mahal_rank} | {mahal_ranked[0][0]} |")

    lsd_all_one = all(r == 1 for r in lor4d_lsd_ranks)
    mahal_all_one = all(r == 1 for r in lor4d_mahal_ranks)
    report.append(f"\n**Lor4D #1 at ALL N (LSD-Well)?** {'✅ YES' if lsd_all_one else '❌ NO'}")
    report.append(f"**Lor4D #1 at ALL N (Mahalanobis)?** {'✅ YES' if mahal_all_one else '❌ NO'}\n")

    # 2. New challenger analysis
    report.append("## 2. New Challenger Analysis\n")
    report.append("Does any of the 8 new families beat Lor4D?\n")
    report.append("| N | Metric | New family challengers (score < Lor4D) |")
    report.append("|---|---|---|")
    any_challenger = False
    for N in N_VALUES:
        if N not in results_lsd:
            continue
        # LSD-Well
        lor4d_lsd = next(s for f, s in results_lsd[N] if f == "Lor4D")
        challengers_lsd = [(f, s) for f, s in results_lsd[N]
                           if f in new_family_names and s < lor4d_lsd]
        if challengers_lsd:
            any_challenger = True
            ch_str = ", ".join(f"{f}({s:.4f})" for f, s in challengers_lsd)
        else:
            ch_str = "None"
        report.append(f"| {N} | LSD-Well | {ch_str} |")

        # Mahalanobis
        lor4d_mahal = next(s for f, s in results_mahal[N] if f == "Lor4D")
        challengers_mahal = [(f, s) for f, s in results_mahal[N]
                             if f in new_family_names and s < lor4d_mahal]
        if challengers_mahal:
            any_challenger = True
            ch_str = ", ".join(f"{f}({s:.2f})" for f, s in challengers_mahal)
        else:
            ch_str = "None"
        report.append(f"| {N} | Mahalanobis | {ch_str} |")

    if not any_challenger:
        report.append("\n**No new adversarial family beats Lor4D at any N.** ✅\n")
    else:
        report.append("\n**⚠ Some new families challenge Lor4D — see details above.**\n")

    # 3. Top-5 at each N
    report.append("## 3. Top-5 Rankings at Each N\n")
    for N in N_VALUES:
        if N not in results_lsd:
            continue
        report.append(f"### N = {N}\n")
        report.append("**LSD-Well Top 5:**\n")
        report.append("| Rank | Family | Score | New? |")
        report.append("|---|---|---|---|")
        for rank, (fam, score) in enumerate(results_lsd[N][:5], 1):
            is_new = "★" if fam in new_family_names else ""
            report.append(f"| {rank} | {fam} | {score:.6f} | {is_new} |")

        report.append("\n**Mahalanobis Top 5:**\n")
        report.append("| Rank | Family | Score | New? |")
        report.append("|---|---|---|---|")
        for rank, (fam, score) in enumerate(results_mahal[N][:5], 1):
            is_new = "★" if fam in new_family_names else ""
            report.append(f"| {rank} | {fam} | {score:.4f} | {is_new} |")
        report.append("")

    # 4. Margin analysis
    report.append("## 4. Margin: Lor4D vs Best Non-Lorentzian\n")
    report.append("| N | LSD-Well margin | Mahalanobis margin |")
    report.append("|---|---|---|")
    lorentzian_fams = {"Lor2D", "Lor3D", "Lor4D", "Lor5D"}
    for N in N_VALUES:
        if N not in results_lsd:
            continue
        lor4d_lsd = next(s for f, s in results_lsd[N] if f == "Lor4D")
        best_nonlor_lsd = min(
            (s for f, s in results_lsd[N] if f not in lorentzian_fams),
            default=float("inf"),
        )
        margin_lsd = best_nonlor_lsd - lor4d_lsd

        lor4d_mahal = next(s for f, s in results_mahal[N] if f == "Lor4D")
        best_nonlor_mahal = min(
            (s for f, s in results_mahal[N] if f not in lorentzian_fams),
            default=float("inf"),
        )
        margin_mahal = best_nonlor_mahal - lor4d_mahal

        report.append(f"| {N} | {margin_lsd:+.6f} | {margin_mahal:+.4f} |")

    # 5. Summary
    report.append("\n## 5. Conclusion\n")
    report.append(f"- Total families tested: {n_fam} (17 original + 8 adversarial)")
    report.append(f"- Lor4D LSD-Well #1 at all N: {'YES' if lsd_all_one else 'NO'}")
    report.append(f"- Lor4D Mahalanobis #1 at all N: {'YES' if mahal_all_one else 'NO'}")
    report.append(f"- New family challenger found: {'YES' if any_challenger else 'NO'}")
    report.append(f"- Runtime: {time.time() - t0:.1f}s")

    # ── Console output ───────────────────────────────────────────────────
    full_report = "\n".join(report)
    print("\n" + full_report)

    # ── Save ─────────────────────────────────────────────────────────────
    outdir = Path(__file__).parent / "outputs_carlip"
    outdir.mkdir(exist_ok=True)
    outpath = outdir / "expanded_family_results.txt"
    outpath.write_text(full_report, encoding="utf-8")
    print(f"\nSaved to {outpath}")


if __name__ == "__main__":
    main()
