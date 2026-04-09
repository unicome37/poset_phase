#!/usr/bin/env python3
"""Prototype: curvature-adaptive reference manifold μ(N,κ).

目标：
1) 构建平坦参考 μ(N,0), Σ^{-1}(N,0)
2) 构建曲率自适应参考 μ(N,κ), Σ^{-1}(N,κ)
3) 对同一目标族分别用两套参考打分，比较 rank/score 稳定性

说明：
- 这是原型脚本，不修改主实验管线。
- 默认参数偏轻量，可快速验证可运行性。
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from curvature_backgrounds import sprinkle_flrw_matter_diamond, poset_from_flrw_matter_points
from expanded_family_robustness import ALL_FAMILIES, compute_features, mahalanobis_score


def _make_flrw_gen(kappa: float) -> Callable[[int, int], object]:
    def _gen(n: int, seed: int):
        pts = sprinkle_flrw_matter_diamond(n=n, d_spatial=3, kappa=kappa, seed=seed)
        return poset_from_flrw_matter_points(pts, kappa=kappa)

    return _gen


def build_reference(
    gen_fn: Callable[[int, int], object],
    n: int,
    reps: int,
    seed_base: int,
) -> tuple[np.ndarray, np.ndarray]:
    feats = []
    for i in range(reps):
        p = gen_fn(n, seed_base + i)
        feats.append(compute_features(p, n))
    mat = np.array(feats)
    mu = mat.mean(axis=0)
    cov = np.cov(mat.T) + 1e-8 * np.eye(3)
    cov_inv = np.linalg.inv(cov)
    return mu, cov_inv


def _target_generators(kappa: float) -> dict[str, Callable[[int, int], object]]:
    return {
        "Lor4D": ALL_FAMILIES["Lor4D"],
        "FLRW_kappa": _make_flrw_gen(kappa),
        "Lor3D": ALL_FAMILIES["Lor3D"],
        "KR_like": ALL_FAMILIES["KR_like"],
    }


def evaluate_once(
    n: int,
    kappa: float,
    reps_eval: int,
    mu_flat: np.ndarray,
    cov_inv_flat: np.ndarray,
    mu_adapt: np.ndarray,
    cov_inv_adapt: np.ndarray,
    seed_base: int,
) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    gens = _target_generators(kappa)

    for fam, gen in gens.items():
        f_list = []
        for r in range(reps_eval):
            p = gen(n, seed_base + 1000 * hash(fam) % 100000 + r)
            f_list.append(compute_features(p, n))
        feat = np.array(f_list).mean(axis=0)

        s_flat = mahalanobis_score(feat, mu_flat, cov_inv_flat)
        s_adapt = mahalanobis_score(feat, mu_adapt, cov_inv_adapt)

        rows.append(
            {
                "n": n,
                "kappa": kappa,
                "family": fam,
                "d_eff": float(feat[0]),
                "c1_c0": float(feat[1]),
                "width_ratio": float(feat[2]),
                "score_flat_ref": float(s_flat),
                "score_adapt_ref": float(s_adapt),
            }
        )

    df = pd.DataFrame(rows)
    df["rank_flat_ref"] = df["score_flat_ref"].rank(method="min").astype(int)
    df["rank_adapt_ref"] = df["score_adapt_ref"].rank(method="min").astype(int)
    return df.sort_values("score_flat_ref").reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prototype μ(N,κ) adaptive reference manifold")
    p.add_argument("--n-values", nargs="+", type=int, default=[256, 512])
    p.add_argument("--kappa", type=float, default=1.0)
    p.add_argument("--ref-reps", type=int, default=12)
    p.add_argument("--eval-reps", type=int, default=6)
    p.add_argument("--seed-base", type=int, default=20260409)
    p.add_argument("--out-dir", type=str, default="outputs_mu_kappa_prototype")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    summary_rows = []

    for n in args.n_values:
        # Flat reference: Lor4D
        mu_flat, cov_inv_flat = build_reference(
            gen_fn=ALL_FAMILIES["Lor4D"],
            n=n,
            reps=args.ref_reps,
            seed_base=args.seed_base + n * 10,
        )

        # Adaptive reference: FLRW(kappa)
        mu_adapt, cov_inv_adapt = build_reference(
            gen_fn=_make_flrw_gen(args.kappa),
            n=n,
            reps=args.ref_reps,
            seed_base=args.seed_base + n * 10 + 500000,
        )

        df_n = evaluate_once(
            n=n,
            kappa=args.kappa,
            reps_eval=args.eval_reps,
            mu_flat=mu_flat,
            cov_inv_flat=cov_inv_flat,
            mu_adapt=mu_adapt,
            cov_inv_adapt=cov_inv_adapt,
            seed_base=args.seed_base + n * 100,
        )
        all_rows.append(df_n)

        row_flrw = df_n[df_n["family"] == "FLRW_kappa"].iloc[0]
        row_lor4d = df_n[df_n["family"] == "Lor4D"].iloc[0]
        summary_rows.append(
            {
                "n": n,
                "kappa": args.kappa,
                "flrw_rank_flat_ref": int(row_flrw["rank_flat_ref"]),
                "flrw_rank_adapt_ref": int(row_flrw["rank_adapt_ref"]),
                "lor4d_rank_flat_ref": int(row_lor4d["rank_flat_ref"]),
                "lor4d_rank_adapt_ref": int(row_lor4d["rank_adapt_ref"]),
                "delta_mu_norm": float(np.linalg.norm(mu_adapt - mu_flat)),
            }
        )

    all_df = pd.concat(all_rows, ignore_index=True)
    sum_df = pd.DataFrame(summary_rows)

    raw_path = out_dir / "mu_kappa_raw.csv"
    sum_path = out_dir / "mu_kappa_summary.csv"
    all_df.to_csv(raw_path, index=False, encoding="utf-8-sig")
    sum_df.to_csv(sum_path, index=False, encoding="utf-8-sig")

    md_lines = [
        "# μ(N,κ) Adaptive Reference Prototype",
        "",
        f"- kappa: {args.kappa}",
        f"- n_values: {args.n_values}",
        f"- ref_reps: {args.ref_reps}",
        f"- eval_reps: {args.eval_reps}",
        "",
        "## Summary",
        "",
        sum_df.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "- 若 `flrw_rank_adapt_ref < flrw_rank_flat_ref`，说明曲率自适应参考对 FLRW κ 分支更友好。",
        "- 若 `lor4d_rank_flat_ref` 维持低位，说明平坦参考仍能稳定识别 Lor4D。",
        "- `delta_mu_norm` 衡量 μ(N,κ) 与 μ(N,0) 的几何偏移规模。",
    ]
    (out_dir / "mu_kappa_report.md").write_text("\n".join(md_lines), encoding="utf-8")

    print("Done.")
    print(f"raw: {raw_path}")
    print(f"summary: {sum_path}")
    print(f"report: {out_dir / 'mu_kappa_report.md'}")


if __name__ == "__main__":
    main()
