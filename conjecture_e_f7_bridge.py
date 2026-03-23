"""Conjecture E — F7↔S_BD bridge: upgrade from old F5 to definitive F7 (sigmoid wall).

This script is the next step after the initial F5-based bridge in
conjecture_e_bridge_fit.py.  It:

1) Regenerates posets from raw_features.csv seeds (same pipeline)
2) Computes F7 (the definitive main model with sigmoid wall, §5.10.7)
3) Computes BD/BDG actions (reusing bd_action.py)
4) Computes ρ_BD(x) local density via layer-block patches
5) Runs the bridge analysis: F7 ~ N + family + bd_ratio
6) Tests patch additivity of ρ_BD(x)

Key upgrade: F7 = logH + γΠ_geo - λΣ_hist + ηΞ_d + α(N)σ((R-Rc)/w)
             where R = interval occupancy ratio (same quantity as bd_ratio uses)

Outputs:
  - outputs_unified_functional/conjecture_e_f7_bridge.csv
  - outputs_unified_functional/conjecture_e_f7_bridge.md
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
    bd_action_d4_truncated,
    bd_ratio_metric,
    bdg_action_d2_corrected,
    bdg_action_d2_link,
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
)
from unified_functional import (
    compute_log_H,
    compute_pi_geo,
    compute_sigma_hist,
    compute_xi_dim,
)


# ── F7 computation ──────────────────────────────────────────────────────

def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_R_from_counts(counts: IntervalCounts) -> float:
    """Interval occupancy ratio R = 1 - C0/total (= 1 - f_link)."""
    if counts.total_relations <= 0:
        return 0.0
    return 1.0 - float(counts.get(0)) / float(counts.total_relations)


def compute_F7_from_components(
    log_H: float,
    pi_geo: float,
    sigma_hist: float,
    xi_dim: float,
    R: float,
    N: int,
    *,
    alpha0: float = 16.0,
    q: float = -0.5,
    lam: float = 10.0,
    eta: float = 0.6,
    Rc: float = 0.25,
    w: float = 0.015,
    N0: float = 20.0,
) -> tuple[float, float]:
    """Compute F7 from pre-computed components. Returns (F7, wall_value)."""
    alpha_N = alpha0 * (N0 / max(N, 1)) ** abs(q)
    wall = alpha_N * sigmoid((R - Rc) / w)
    F7 = log_H + 0.0004 * pi_geo - lam * sigma_hist + eta * xi_dim + wall
    return F7, wall


# ── ρ_BD local density ──────────────────────────────────────────────────

def layer_partition(poset: Poset) -> list[np.ndarray]:
    """Return a topological layer partition as arrays of vertex indices."""
    closure = poset.closure
    indeg = closure.sum(axis=0).astype(int)
    remaining = set(range(poset.n))
    layers: list[np.ndarray] = []

    while remaining:
        mins = [i for i in remaining if indeg[i] == 0]
        if not mins:
            raise ValueError("Input relation is not a poset: cycle detected.")
        layer = np.array(sorted(mins), dtype=int)
        layers.append(layer)
        for u in mins:
            remaining.remove(u)
            for v in np.where(closure[u])[0]:
                indeg[v] -= 1
    return layers


def compute_rho_bd_local(
    poset: Poset,
    *,
    block_sizes: tuple[int, ...] = (2, 3),
    seed: int = 42,
) -> dict:
    """Compute local BD density by splitting poset into layer-block patches.

    Each patch is the induced sub-poset on a contiguous block of 2–3 layers.
    This is intentionally closer to the main layered structure than Alexandrov
    interval patches, and is therefore the preferred ρ_BD carrier in the F7
    bridge stage.

    Returns dict with:
      - rho_bd_patches: list of local bd_ratio values
      - rho_bd_mean: mean local density
      - rho_bd_std: std of local density
      - global_bd_ratio: global bd_ratio for comparison
      - additivity_ratio: mean(local) / global
      - n_patches: number of patches used
    """
    layers = layer_partition(poset)
    if len(layers) < min(block_sizes):
        global_counts = count_intervals_fast(poset)
        return {
            "rho_bd_patches": [],
            "rho_bd_mean": float("nan"),
            "rho_bd_std": float("nan"),
            "global_bd_ratio": bd_ratio_metric(global_counts),
            "additivity_ratio": float("nan"),
            "n_patches": 0,
            "patch_mode": "layer-block",
        }

    candidates: list[np.ndarray] = []
    for block_size in block_sizes:
        if block_size <= 0:
            continue
        for start in range(0, len(layers) - block_size + 1):
            block = np.concatenate(layers[start : start + block_size])
            if block.size >= 3:
                candidates.append(block)

    if not candidates:
        global_counts = count_intervals_fast(poset)
        return {
            "rho_bd_patches": [],
            "rho_bd_mean": float("nan"),
            "rho_bd_std": float("nan"),
            "global_bd_ratio": bd_ratio_metric(global_counts),
            "additivity_ratio": float("nan"),
            "n_patches": 0,
            "patch_mode": "layer-block",
        }

    rng = np.random.RandomState(seed)
    rng.shuffle(candidates)

    local_ratios = []
    for block in candidates:
        sub_closure = poset.closure[np.ix_(block, block)]
        sub_poset = Poset(closure=sub_closure)
        sub_counts = count_intervals_fast(sub_poset)
        local_ratios.append(bd_ratio_metric(sub_counts))

    global_counts = count_intervals_fast(poset)
    global_ratio = bd_ratio_metric(global_counts)
    mean_local = float(np.mean(local_ratios))
    std_local = float(np.std(local_ratios))
    additivity = mean_local / global_ratio if global_ratio > 0 else float("nan")

    return {
        "rho_bd_patches": local_ratios,
        "rho_bd_mean": mean_local,
        "rho_bd_std": std_local,
        "global_bd_ratio": global_ratio,
        "additivity_ratio": additivity,
        "n_patches": len(local_ratios),
        "patch_mode": "layer-block",
    }


# ── Helpers ──────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class F7BridgeRow:
    family: str
    N: int
    rep: int
    seed: int
    log_H: float
    pi_geo: float
    sigma_hist: float
    xi_dim: float
    R: float
    wall: float
    F7: float
    f5_calibrated: float
    bd_ratio: float
    bdg_d2_corrected_norm: float
    bdg_d4_standard_norm: float
    rho_bd_mean: float
    rho_bd_std: float
    additivity_ratio: float
    n_patches: int
    patch_mode: str


def seed_for(family: str, n: int, rep: int, base_seed: int) -> int:
    return int(base_seed + rep * 1000 + n * 100)


FAMILY_GENS = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}


def spearman_rank(x: np.ndarray, y: np.ndarray) -> float:
    if x.size != y.size or x.size == 0:
        return float("nan")

    def rankdata(a: np.ndarray) -> np.ndarray:
        order = np.argsort(a, kind="mergesort")
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, a.size + 1, dtype=float)
        sorted_a = a[order]
        i = 0
        while i < a.size:
            j = i + 1
            while j < a.size and sorted_a[j] == sorted_a[i]:
                j += 1
            if j - i > 1:
                ranks[order[i:j]] = (i + 1 + j) / 2.0
            i = j
        return ranks

    rx = rankdata(x)
    ry = rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = float(np.sqrt(np.sum(rx * rx) * np.sum(ry * ry)))
    return float(np.sum(rx * ry) / denom) if denom > 0 else float("nan")


def ols_r2(y: np.ndarray, X: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return r2, beta, pred


def design_matrix(rows: list[dict], extra_keys: list[str] | None = None) -> np.ndarray:
    fams = ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]
    n_values = sorted({r["N"] for r in rows})
    cols = []
    for r in rows:
        row = [1.0]
        for nv in n_values[:-1]:
            row.append(1.0 if r["N"] == nv else 0.0)
        for ff in fams[:-1]:
            row.append(1.0 if r["family"] == ff else 0.0)
        if extra_keys:
            for k in extra_keys:
                row.append(float(r[k]))
        cols.append(row)
    return np.array(cols, dtype=float)


def hit_rate_matched(focus: list[float], control: list[float]) -> float:
    m = min(len(focus), len(control))
    if m == 0:
        return float("nan")
    return sum(1 for i in range(m) if focus[i] < control[i]) / float(m)


F5_WEIGHTS = {"beta": 2.0, "gamma": 0.5, "lam": 1.5, "eta": 0.1, "kappa": 0.05}


def compute_f5_calibrated(raw_row: dict) -> float:
    return (
        F5_WEIGHTS["beta"] * float(raw_row["log_H"])
        + F5_WEIGHTS["gamma"] * float(raw_row["pi_geo"])
        - F5_WEIGHTS["lam"] * float(raw_row["sigma_hist"])
        + F5_WEIGHTS["eta"] * float(raw_row["xi_dim"])
        + F5_WEIGHTS["kappa"] * float(raw_row["pi_cg"])
    )


# ── Main ─────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(description="Conjecture E: F7↔S_BD bridge with sigmoid wall")
    ap.add_argument("--raw-features", default="outputs_unified_functional/raw_features.csv")
    ap.add_argument("--out-csv", default="outputs_unified_functional/conjecture_e_f7_bridge.csv")
    ap.add_argument("--out-md", default="outputs_unified_functional/conjecture_e_f7_bridge.md")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--block-sizes", default="2,3", help="comma-separated layer-block sizes for ρ_BD")
    args = ap.parse_args()

    raw_path = Path(args.raw_features)
    if not raw_path.exists():
        raise FileNotFoundError(str(raw_path))

    raw_rows = list(csv.DictReader(raw_path.open("r", encoding="utf-8")))
    families = {"Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"}
    raw_rows = [r for r in raw_rows if r.get("family") in families]
    print(f"Loaded {len(raw_rows)} rows from {raw_path}")

    block_sizes = tuple(int(x) for x in args.block_sizes.split(",") if x.strip())

    out_rows: list[F7BridgeRow] = []
    for i, r in enumerate(raw_rows, start=1):
        fam = r["family"]
        n = int(r["N"])
        rep = int(r["rep"])
        s = seed_for(fam, n, rep, args.seed)
        poset = FAMILY_GENS[fam](n, seed=s)

        # F7 components
        log_H = float(r["log_H"])
        pi_geo = float(r["pi_geo"])
        sigma_hist = float(r["sigma_hist"])
        xi_dim = float(r["xi_dim"])

        # BD/BDG from interval counts
        counts = count_intervals_fast(poset)
        R = compute_R_from_counts(counts)
        br = bd_ratio_metric(counts)
        bdg2c = bdg_action_d2_corrected(counts, n, normalized=True)
        bdg4s = bdg_action_d4_standard(counts, n, normalized=True)

        # F7 with sigmoid wall
        F7, wall_val = compute_F7_from_components(log_H, pi_geo, sigma_hist, xi_dim, R, n)
        f5_cal = compute_f5_calibrated(r)

        # ρ_BD local density
        rho_info = compute_rho_bd_local(poset, block_sizes=block_sizes, seed=s)

        out_rows.append(F7BridgeRow(
            family=fam, N=n, rep=rep, seed=s,
            log_H=log_H, pi_geo=pi_geo, sigma_hist=sigma_hist, xi_dim=xi_dim,
            R=R, wall=wall_val, F7=F7, f5_calibrated=f5_cal,
            bd_ratio=br, bdg_d2_corrected_norm=bdg2c, bdg_d4_standard_norm=bdg4s,
            rho_bd_mean=rho_info["rho_bd_mean"],
            rho_bd_std=rho_info["rho_bd_std"],
            additivity_ratio=rho_info["additivity_ratio"],
            n_patches=rho_info["n_patches"],
            patch_mode=rho_info["patch_mode"],
        ))
        if i % 10 == 0:
            print(f"  [{i}/{len(raw_rows)}] {fam} N={n} rep={rep}  F7={F7:.2f}  R={R:.4f}  wall={wall_val:.2f}  bd_ratio={br:.4f}  patches={rho_info['n_patches']}")

    # ── Save CSV ──
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(out_rows[0]).keys()))
        w.writeheader()
        for rr in out_rows:
            w.writerow(asdict(rr))
    print(f"\nSaved: {out_csv}")

    # ── Analysis ──
    data = [asdict(rr) for rr in out_rows]

    # 1. F7 vs F5 comparison
    f7_arr = np.array([d["F7"] for d in data], dtype=float)
    f5_arr = np.array([d["f5_calibrated"] for d in data], dtype=float)
    rho_f7_f5 = spearman_rank(f7_arr, f5_arr)
    print(f"\n=== F7 vs F5 ===")
    print(f"Spearman(F7, F5) = {rho_f7_f5:+.4f}")

    # 2. F7 vs BD observables
    br_arr = np.array([d["bd_ratio"] for d in data], dtype=float)
    r_arr = np.array([d["R"] for d in data], dtype=float)
    wall_arr = np.array([d["wall"] for d in data], dtype=float)
    print(f"\n=== F7 vs BD observables ===")
    print(f"Spearman(F7, bd_ratio)  = {spearman_rank(f7_arr, br_arr):+.4f}")
    print(f"Spearman(F7, R)         = {spearman_rank(f7_arr, r_arr):+.4f}")
    print(f"Spearman(wall, R)       = {spearman_rank(wall_arr, r_arr):+.4f}")
    print(f"Spearman(wall, bd_ratio)= {spearman_rank(wall_arr, br_arr):+.4f}")

    # 3. Residual bridge: F7 ~ N + family + bd_ratio
    y_f7 = f7_arr.copy()
    X_base = design_matrix(data, extra_keys=None)
    r2_base, _, pred_base = ols_r2(y_f7, X_base)
    resid = y_f7 - pred_base

    bridge_results = []
    for key in ["bd_ratio", "bdg_d2_corrected_norm", "bdg_d4_standard_norm", "R"]:
        X_ext = design_matrix(data, extra_keys=[key])
        r2_ext, beta_ext, _ = ols_r2(y_f7, X_ext)
        xs = np.array([d[key] for d in data], dtype=float)
        coef = float(beta_ext[-1])
        bridge_results.append({
            "metric": key,
            "r2_base": r2_base,
            "r2_ext": r2_ext,
            "delta_r2": r2_ext - r2_base,
            "coef": coef,
            "spearman_resid": spearman_rank(xs, resid),
        })

    print(f"\n=== Residual bridge: F7 ~ N + family ===")
    print(f"Baseline R2 = {r2_base:.4f}")
    for br in bridge_results:
        print(f"  + {br['metric']:30s}: dR2={br['delta_r2']:+.4f}, coef={br['coef']:+.4f}, rho(resid)={br['spearman_resid']:+.4f}")

    # 4. F7 Lor4D vs KR win rates
    ns = sorted({d["N"] for d in data})
    print(f"\n=== F7 ordering: Lor4D vs KR_like ===")
    print(f"{'N':>4} | F7_win | F5_win")
    for n in ns:
        lor = sorted([d for d in data if d["family"] == "Lor4D" and d["N"] == n], key=lambda x: x["rep"])
        kr = sorted([d for d in data if d["family"] == "KR_like" and d["N"] == n], key=lambda x: x["rep"])
        m = min(len(lor), len(kr))
        if m == 0:
            continue
        f7_lor = [d["F7"] for d in lor[:m]]
        f7_kr = [d["F7"] for d in kr[:m]]
        f5_lor = [d["f5_calibrated"] for d in lor[:m]]
        f5_kr = [d["f5_calibrated"] for d in kr[:m]]
        wr7 = hit_rate_matched(f7_lor, f7_kr)
        wr5 = hit_rate_matched(f5_lor, f5_kr)
        print(f"{n:>4} | {wr7:.3f}  | {wr5:.3f}")

    # 5. ρ_BD local density additivity check
    valid_rho = [d for d in data if not math.isnan(d["additivity_ratio"])]
    if valid_rho:
        add_arr = np.array([d["additivity_ratio"] for d in valid_rho], dtype=float)
        print(f"\n=== ρ_BD local density additivity ===")
        print(f"Samples with patches: {len(valid_rho)}/{len(data)}")
        print(f"Additivity ratio: mean={np.mean(add_arr):.3f}, std={np.std(add_arr):.3f}")
        print(f"  (=1 means perfectly additive, <1 means local underestimates global)")
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
            sub = [d for d in valid_rho if d["family"] == fam]
            if sub:
                add_sub = np.array([d["additivity_ratio"] for d in sub], dtype=float)
                print(f"  {fam:>7s}: add_ratio={np.mean(add_sub):.3f}±{np.std(add_sub):.3f}  (n={len(sub)})")

    # 6. Per-family F7 means
    print(f"\n=== Per-family means ===")
    print(f"{'Family':>7s} | {'F7':>8s} | {'F5':>8s} | {'R':>6s} | {'wall':>6s} | {'bd_ratio':>8s}")
    for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
        sub = [d for d in data if d["family"] == fam]
        if not sub:
            continue
        mf7 = np.mean([d["F7"] for d in sub])
        mf5 = np.mean([d["f5_calibrated"] for d in sub])
        mr = np.mean([d["R"] for d in sub])
        mw = np.mean([d["wall"] for d in sub])
        mbr = np.mean([d["bd_ratio"] for d in sub])
        print(f"{fam:>7s} | {mf7:8.2f} | {mf5:8.2f} | {mr:6.4f} | {mw:6.2f} | {mbr:8.4f}")

    # ── Save markdown report ──
    md = Path(args.out_md)
    lines = []
    lines.append("# Conjecture E: F7↔S_BD Bridge (sigmoid wall upgrade)")
    lines.append("")
    lines.append("## F7 main model (§5.10.7)")
    lines.append("")
    lines.append("F7 = logH + 0.0004·Π_geo - 10·Σ_hist + 0.6·Ξ_d + α(N)·σ((R-Rc)/w)")
    lines.append("")
    lines.append("α₀=16, q=-0.5, Rc=0.25, w=0.015, N₀=20")
    lines.append("")
    lines.append(f"## F7 vs F5 correlation")
    lines.append(f"")
    lines.append(f"Spearman(F7, F5) = {rho_f7_f5:+.4f}")
    lines.append(f"")
    lines.append(f"## F7 vs BD observables")
    lines.append("")
    lines.append(f"| pair | Spearman |")
    lines.append(f"|---|---:|")
    lines.append(f"| F7 vs bd_ratio | {spearman_rank(f7_arr, br_arr):+.4f} |")
    lines.append(f"| F7 vs R | {spearman_rank(f7_arr, r_arr):+.4f} |")
    lines.append(f"| wall vs R | {spearman_rank(wall_arr, r_arr):+.4f} |")
    lines.append(f"| wall vs bd_ratio | {spearman_rank(wall_arr, br_arr):+.4f} |")
    lines.append("")
    lines.append(f"## Residual bridge: F7 ~ N + family (+BD)")
    lines.append("")
    lines.append(f"Baseline R² = {r2_base:.4f}")
    lines.append("")
    lines.append("| metric | R² ext | ΔR² | coef | ρ(resid) |")
    lines.append("|---|---:|---:|---:|---:|")
    for br in bridge_results:
        lines.append(f"| {br['metric']} | {br['r2_ext']:.4f} | {br['delta_r2']:+.4f} | {br['coef']:+.4f} | {br['spearman_resid']:+.4f} |")
    lines.append("")
    lines.append("## F7 ordering: Lor4D vs KR_like")
    lines.append("")
    lines.append("| N | F7 win | F5 win |")
    lines.append("|---:|---:|---:|")
    for n in ns:
        lor = sorted([d for d in data if d["family"] == "Lor4D" and d["N"] == n], key=lambda x: x["rep"])
        kr = sorted([d for d in data if d["family"] == "KR_like" and d["N"] == n], key=lambda x: x["rep"])
        m = min(len(lor), len(kr))
        if m == 0:
            continue
        wr7 = hit_rate_matched([d["F7"] for d in lor[:m]], [d["F7"] for d in kr[:m]])
        wr5 = hit_rate_matched([d["f5_calibrated"] for d in lor[:m]], [d["f5_calibrated"] for d in kr[:m]])
        lines.append(f"| {n} | {wr7:.3f} | {wr5:.3f} |")
    lines.append("")
    lines.append("## ρ_BD local density additivity")
    lines.append("")
    if valid_rho:
        lines.append(f"Samples with patches: {len(valid_rho)}/{len(data)}")
        lines.append(f"")
        lines.append("| family | add_ratio | std | n |")
        lines.append("|---|---:|---:|---:|")
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
            sub = [d for d in valid_rho if d["family"] == fam]
            if sub:
                add_sub = np.array([d["additivity_ratio"] for d in sub], dtype=float)
                lines.append(f"| {fam} | {np.mean(add_sub):.3f} | {np.std(add_sub):.3f} | {len(sub)} |")
    else:
        lines.append("No patches found (N too small or posets too sparse).")
    lines.append("")
    lines.append("## Per-family means")
    lines.append("")
    lines.append("| family | F7 | F5 | R | wall | bd_ratio |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
        sub = [d for d in data if d["family"] == fam]
        if not sub:
            continue
        lines.append(f"| {fam} | {np.mean([d['F7'] for d in sub]):.2f} | {np.mean([d['f5_calibrated'] for d in sub]):.2f} | {np.mean([d['R'] for d in sub]):.4f} | {np.mean([d['wall'] for d in sub]):.2f} | {np.mean([d['bd_ratio'] for d in sub]):.4f} |")
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("The F7↔S_BD bridge tests whether the sigmoid wall (which uses R = 1 - f_link) ")
    lines.append("is already capturing the same information as the BD action's interval statistics. ")
    lines.append("A high Spearman(wall, bd_ratio) confirms that the sigmoid wall in F7 is a ")
    lines.append("monotone proxy for the same causal interval richness that BD actions encode.")
    lines.append("")
    lines.append("The ρ_BD local density additivity check tests whether bd_ratio can be ")
    lines.append("decomposed into local layer-block patches — a prerequisite for writing it as ")
    lines.append("ρ_BD(x) in the continuum limit.")

    md.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"Saved: {md}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
