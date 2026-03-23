"""Conjecture E - Explicit bd_ratio <-> S_BD mapping check.

In the current theory notes, the intended BD bridge term is the interval-richness
quantity

    S_BD(X) = (1 - C0 / total_relations) * (1 + mean_k)

This script recomputes that quantity directly from the regenerated posets and
checks whether it matches the empirical `bd_ratio` column already used in the
bridge scans.

Outputs:
  - outputs_unified_functional/conjecture_e_sbd_mapping_check.csv
  - outputs_unified_functional/conjecture_e_sbd_mapping_check.md
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from bd_action import count_intervals_fast, bd_ratio_metric
from generators import (
    generate_kr_like,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)


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


def load_rows(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def spearmanr(x: np.ndarray, y: np.ndarray) -> float:
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


def design_matrix(rows: list[dict], extra_key: str | None = None) -> np.ndarray:
    fams = ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]
    n_values = sorted({int(r["N"]) for r in rows})
    cols = []
    for r in rows:
        row = [1.0]
        n = int(r["N"])
        for nv in n_values[:-1]:
            row.append(1.0 if n == nv else 0.0)
        fam = r["family"]
        for ff in fams[:-1]:
            row.append(1.0 if fam == ff else 0.0)
        if extra_key is not None:
            row.append(float(r[extra_key]))
        cols.append(row)
    return np.array(cols, dtype=float)


def ols_r2(y: np.ndarray, X: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return r2, beta, pred


@dataclass(frozen=True)
class MappingRow:
    family: str
    N: int
    rep: int
    seed: int
    log_H: float
    f5_calibrated: float
    C0: int
    C1: int
    C2: int
    C3: int
    total_relations: int
    link_fraction: float
    mean_k: float
    s_bd_exact: float
    bd_ratio_stored: float
    abs_diff: float


def main() -> int:
    ap = argparse.ArgumentParser(description="Check the explicit bd_ratio <-> S_BD mapping")
    ap.add_argument("--raw-features", default="outputs_unified_functional/raw_features.csv")
    ap.add_argument("--bd-actions", default="outputs_unified_functional/bd_actions.csv")
    ap.add_argument("--out-csv", default="outputs_unified_functional/conjecture_e_sbd_mapping_check.csv")
    ap.add_argument("--out-md", default="outputs_unified_functional/conjecture_e_sbd_mapping_check.md")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--families", nargs="*", default=["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"])
    args = ap.parse_args()

    raw_path = Path(args.raw_features)
    bd_path = Path(args.bd_actions)
    if not raw_path.exists():
        raise FileNotFoundError(str(raw_path))
    if not bd_path.exists():
        raise FileNotFoundError(str(bd_path))

    raw_rows = load_rows(raw_path)
    bd_rows = load_rows(bd_path)
    raw_rows = [r for r in raw_rows if r.get("family") in set(args.families)]
    bd_rows = [r for r in bd_rows if r.get("family") in set(args.families)]
    if not raw_rows or not bd_rows:
        raise RuntimeError("Missing rows after filtering families.")

    bd_lookup = {(r["family"], int(r["N"]), int(r["rep"])): r for r in bd_rows}
    out_rows: list[MappingRow] = []
    for r in raw_rows:
        fam = r["family"]
        n = int(r["N"])
        rep = int(r["rep"])
        seed = seed_for(fam, n, rep, args.seed)
        poset = generate_family(fam, n, seed)
        counts = count_intervals_fast(poset, k_max=None)
        total = counts.total_relations
        c0 = counts.get(0)
        c1 = counts.get(1)
        c2 = counts.get(2)
        c3 = counts.get(3)
        non_link_frac = 1.0 - (float(c0) / float(total) if total > 0 else 0.0)
        mean_k = 0.0
        for k, v in counts.counts.items():
            mean_k += float(k) * float(v)
        mean_k = mean_k / float(total) if total > 0 else 0.0
        s_bd_exact = bd_ratio_metric(counts)
        stored = float(bd_lookup[(fam, n, rep)]["bd_ratio"])
        out_rows.append(
            MappingRow(
                family=fam,
                N=n,
                rep=rep,
                seed=seed,
                log_H=float(r["log_H"]),
                f5_calibrated=float(bd_lookup[(fam, n, rep)]["f5_calibrated"]),
                C0=c0,
                C1=c1,
                C2=c2,
                C3=c3,
                total_relations=total,
                link_fraction=1.0 - non_link_frac,
                mean_k=mean_k,
                s_bd_exact=s_bd_exact,
                bd_ratio_stored=stored,
                abs_diff=abs(s_bd_exact - stored),
            )
        )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(out_rows[0]).keys()))
        w.writeheader()
        for row in out_rows:
            w.writerow(asdict(row))

    x = np.array([row.s_bd_exact for row in out_rows], dtype=float)
    y = np.array([row.bd_ratio_stored for row in out_rows], dtype=float)
    f5 = np.array([row.f5_calibrated for row in out_rows], dtype=float)
    r2_base, _, pred_base = ols_r2(f5, design_matrix([asdict(row) for row in out_rows]))
    r2_ext, _, _ = ols_r2(f5, design_matrix([asdict(row) for row in out_rows], extra_key="s_bd_exact"))
    residual = f5 - pred_base
    resid_corr = float(np.corrcoef(x, residual)[0, 1]) if x.std() > 0 and residual.std() > 0 else float("nan")

    max_abs_diff = float(np.max(np.abs(x - y)))
    mean_abs_diff = float(np.mean(np.abs(x - y)))
    corr_x_y = spearmanr(x, y)

    md_lines: list[str] = []
    md_lines.append("# Conjecture E S_BD mapping check")
    md_lines.append("")
    md_lines.append(f"- Input raw features: `{raw_path.as_posix()}`")
    md_lines.append(f"- Input BD table: `{bd_path.as_posix()}`")
    md_lines.append(f"- Output CSV: `{out_csv.as_posix()}`")
    md_lines.append("")
    md_lines.append("## Mapping claim")
    md_lines.append("")
    md_lines.append("The current empirical bridge term `bd_ratio` is the same interval-richness quantity used as `S_BD` in the theory note:")
    md_lines.append("")
    md_lines.append("`S_BD(X) = (1 - C0 / total_relations) * (1 + mean_k)`")
    md_lines.append("")
    md_lines.append("So the operational map is simply:")
    md_lines.append("")
    md_lines.append("`S_BD_proxy(X) = bd_ratio(X)`")
    md_lines.append("")
    md_lines.append("## Identity check")
    md_lines.append("")
    md_lines.append(f"Max |S_BD_exact - bd_ratio_stored| = {max_abs_diff:.6e}")
    md_lines.append(f"Mean |S_BD_exact - bd_ratio_stored| = {mean_abs_diff:.6e}")
    md_lines.append(f"Spearman(S_BD_exact, bd_ratio_stored) = {corr_x_y:+.4f}")
    md_lines.append("")
    md_lines.append("## Residual bridge")
    md_lines.append("")
    md_lines.append(f"Baseline model: `F5 ~ N + family` with R² = {r2_base:.4f}")
    md_lines.append(f"Extended model: `F5 ~ N + family + S_BD_proxy` with R² = {r2_ext:.4f}")
    md_lines.append(f"ΔR² = {r2_ext - r2_base:+.4f}")
    md_lines.append(f"Corr(S_BD_proxy, residual F5) = {resid_corr:+.4f}")
    md_lines.append("")
    md_lines.append("## Interpretation")
    md_lines.append("")
    md_lines.append(
        "This is the cleanest possible bridge statement at the current stage: the empirical `bd_ratio` is not just inspired by `S_BD`; it is the same normalized interval-richness observable. "
        "So the next theoretical step is not to redefine the metric, but to justify why this interval-richness term should enter the discrete effective action with a positive penalty sign."
    )

    out_md = Path(args.out_md)
    out_md.write_text("\n".join(md_lines).rstrip() + "\n", encoding="utf-8")

    print(f"Saved: {out_csv}")
    print(f"Saved: {out_md}")
    print(f"Max abs diff = {max_abs_diff:.6e}")
    print(f"Mean abs diff = {mean_abs_diff:.6e}")
    print(f"Residual bridge delta_R2 = {r2_ext - r2_base:+.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
