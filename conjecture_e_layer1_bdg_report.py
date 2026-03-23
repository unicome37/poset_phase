"""Conjecture E — Layer 1 baseline: BD/BDG action reproduction on existing ensembles.

This script targets the "first layer" deliverable described in
`理论体系/推论核查与验证优先级分析.md`:

1) Implement a BD/BDG action calculator (see `bd_action.py`).
2) Recompute BD/BDG observables on the exact same dataset as `raw_features.csv`,
   using the same seed rule, and export a compact table for follow-up analysis.

Focus: numerical reproduction + interface hooks, not a continuum-limit proof.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from bd_action import (
    bd_action_d4_truncated,
    bd_ratio_metric,
    bdg_action_d2_corrected,
    bdg_action_d2_link,
    bdg_action_d4_standard,
    count_intervals_fast,
)
from generators import (
    generate_kr_like,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)


@dataclass(frozen=True)
class BdRow:
    family: str
    N: int
    rep: int
    seed: int
    log_H: float
    f5_calibrated: float
    bdg_d2_link_norm: float
    bdg_d2_corrected_norm: float
    bdg_d4_standard_norm: float
    bd_d4_trunc_norm: float
    bd_ratio: float
    C0: int
    C1: int
    C2: int
    C3: int
    total_relations: int


def seed_for(family: str, n: int, rep: int, base_seed: int) -> int:
    return int(base_seed + rep * 1000 + n * 100)


def generate_family(family: str, n: int, seed: int):
    gens = {
        "Lor2D": generate_lorentzian_like_2d,
        "Lor3D": generate_lorentzian_like_3d,
        "Lor4D": generate_lorentzian_like_4d,
        "Lor5D": generate_lorentzian_like_5d,
        "KR_like": generate_kr_like,
    }
    if family not in gens:
        raise ValueError(f"Unsupported family: {family}")
    return gens[family](n, seed=seed)


def spearmanr(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation without SciPy (ties handled by average ranks)."""
    if x.size != y.size or x.size == 0:
        return float("nan")

    def rankdata(a: np.ndarray) -> np.ndarray:
        order = np.argsort(a, kind="mergesort")
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, a.size + 1, dtype=float)
        # average ranks for ties
        sorted_a = a[order]
        i = 0
        while i < a.size:
            j = i + 1
            while j < a.size and sorted_a[j] == sorted_a[i]:
                j += 1
            if j - i > 1:
                avg = (i + 1 + j) / 2.0
                ranks[order[i:j]] = avg
            i = j
        return ranks

    rx = rankdata(x)
    ry = rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = float(np.sqrt(np.sum(rx * rx) * np.sum(ry * ry)))
    return float(np.sum(rx * ry) / denom) if denom > 0 else float("nan")


CALIBRATED_WEIGHTS = {
    "beta": 2.0,
    "gamma": 0.5,
    "lam": 1.5,
    "eta": 0.1,
    "kappa": 0.05,
}


def compute_f5_calibrated(raw_row: dict) -> float:
    return (
        CALIBRATED_WEIGHTS["beta"] * float(raw_row["log_H"])
        + CALIBRATED_WEIGHTS["gamma"] * float(raw_row["pi_geo"])
        - CALIBRATED_WEIGHTS["lam"] * float(raw_row["sigma_hist"])
        + CALIBRATED_WEIGHTS["eta"] * float(raw_row["xi_dim"])
        + CALIBRATED_WEIGHTS["kappa"] * float(raw_row["pi_cg"])
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="Conjecture E layer-1: recompute BD/BDG actions on raw_features.csv")
    ap.add_argument("--raw-features", default="outputs_unified_functional/raw_features.csv")
    ap.add_argument("--out", default="outputs_unified_functional/bd_actions.csv")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--families", nargs="*", default=["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"])
    ap.add_argument("--k-max", type=int, default=3, help="max interval size k to retain in output (counts still computed exactly)")
    args = ap.parse_args()

    raw_path = Path(args.raw_features)
    if not raw_path.exists():
        raise FileNotFoundError(str(raw_path))

    rows = list(csv.DictReader(raw_path.open("r", encoding="utf-8")))
    rows = [r for r in rows if r.get("family") in set(args.families)]
    if not rows:
        raise RuntimeError("No rows selected from raw_features.csv (check --families)")

    out_rows: list[BdRow] = []
    for i, r in enumerate(rows, start=1):
        fam = r["family"]
        n = int(r["N"])
        rep = int(r["rep"])
        s = seed_for(fam, n, rep, args.seed)
        poset = generate_family(fam, n, s)

        counts = count_intervals_fast(poset, k_max=None)
        c0, c1, c2, c3 = counts.get(0), counts.get(1), counts.get(2), counts.get(3)

        out_rows.append(
            BdRow(
                family=fam,
                N=n,
                rep=rep,
                seed=s,
                log_H=float(r["log_H"]),
                f5_calibrated=compute_f5_calibrated(r),
                bdg_d2_link_norm=bdg_action_d2_link(counts, n, normalized=True),
                bdg_d2_corrected_norm=bdg_action_d2_corrected(counts, n, normalized=True),
                bdg_d4_standard_norm=bdg_action_d4_standard(counts, n, normalized=True),
                bd_d4_trunc_norm=bd_action_d4_truncated(counts, n, normalized=True),
                bd_ratio=bd_ratio_metric(counts),
                C0=c0,
                C1=c1,
                C2=c2,
                C3=c3,
                total_relations=int(counts.total_relations),
            )
        )
        if i % 25 == 0:
            print(f"  [{i}/{len(rows)}] computed")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(asdict(out_rows[0]).keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for rr in out_rows:
            w.writerow(asdict(rr))
    print(f"Saved: {out_path}")

    # Quick summary: per-family means and rank correlations with log_H
    print("\nPer-family mean (pooled N):")
    fams = sorted(set(r.family for r in out_rows))
    for fam in fams:
        sub = [r for r in out_rows if r.family == fam]
        mean_logh = float(np.mean([r.log_H for r in sub]))
        mean_f5 = float(np.mean([r.f5_calibrated for r in sub]))
        mean_bdg4 = float(np.mean([r.bdg_d4_standard_norm for r in sub]))
        mean_link2 = float(np.mean([r.bdg_d2_link_norm for r in sub]))
        mean_ratio = float(np.mean([r.bd_ratio for r in sub]))
        print(
            f"  {fam:>7s}: log_H={mean_logh:8.3f} | F5={mean_f5:8.3f} | "
            f"BDG4={mean_bdg4:+8.4f} | d2_link={mean_link2:+8.4f} | bd_ratio={mean_ratio:8.4f}"
        )

    x = np.array([r.log_H for r in out_rows], dtype=float)
    f5 = np.array([r.f5_calibrated for r in out_rows], dtype=float)
    y = np.array([r.bdg_d4_standard_norm for r in out_rows], dtype=float)
    z = np.array([r.bdg_d2_link_norm for r in out_rows], dtype=float)
    print("\nSpearman correlations (pooled):")
    print(f"  Spearman(log_H, BDG4_norm)   = {spearmanr(x, y):+.4f}")
    print(f"  Spearman(log_H, d2_link_norm)= {spearmanr(x, z):+.4f}")
    print(f"  Spearman(F5, BDG4_norm)      = {spearmanr(f5, y):+.4f}")
    print(f"  Spearman(F5, d2_link_norm)   = {spearmanr(f5, z):+.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
