#!/usr/bin/env python3
"""β-1 增强版：独立 γ 序贯筛选 — Step1(A2, γ₁=soft) → TOP-K → Step2(A4, γ₂=hard)。

核心思想：
  几何惩罚在低 γ 下做"宽容预筛"（保留 Lor4D），
  信息论惩罚在高 γ 下做"精确选择"（Lor4D #1 in small pool）。
"""

import pandas as pd
import numpy as np
from pathlib import Path

OUT = Path("outputs_info_hybrid")
RAW = OUT / "hybrid_raw.csv"


def load_family_scores(raw_csv: Path) -> pd.DataFrame:
    """汇总为 (n, gamma, action_mode, family) → mean_score."""
    df = pd.read_csv(raw_csv)
    grouped = df.groupby(["n", "gamma", "action_mode", "family"])["score"].mean()
    return grouped.reset_index()


def get_ranking(df: pd.DataFrame, n: int, gamma: float, mode: str,
                subset: list | None = None) -> dict:
    """返回 {family: rank}。"""
    mask = (df["n"] == n) & (df["gamma"] == gamma) & (df["action_mode"] == mode)
    sub = df[mask].copy()
    if subset is not None:
        sub = sub[sub["family"].isin(subset)]
    sub = sub.sort_values("score")
    return dict(zip(sub["family"], range(1, len(sub) + 1)))


def main():
    df = load_family_scores(RAW)

    print("=" * 110)
    print("β-1 增强版：独立 γ 序贯筛选 — Step1(A2, γ₁) → TOP-K → Step2(A4, γ₂)")
    print("=" * 110)

    gamma1_candidates = [0.0, 0.1, 0.2, 0.5]   # Step 1: soft geo pre-filter
    gamma2_candidates = [0.2, 0.5, 1.0, 2.0]    # Step 2: hard info select
    K_candidates = [3, 5, 8]

    lor4d = "lorentzian_like_4d"
    results = []

    for n in [10, 20, 40]:
        for g1 in gamma1_candidates:
            # Step 1 ranking
            r1 = get_ranking(df, n, g1, "A2")
            if lor4d not in r1:
                continue
            lor4d_r1 = r1[lor4d]

            for K in K_candidates:
                top_k = [f for f, r in sorted(r1.items(), key=lambda x: x[1]) if r <= K]
                in_topk = lor4d in top_k

                for g2 in gamma2_candidates:
                    if not in_topk:
                        seq_rank = None
                    else:
                        r2 = get_ranking(df, n, g2, "A4", subset=top_k)
                        seq_rank = r2.get(lor4d)

                    # 对照基线
                    a4_full = get_ranking(df, n, g2, "A4")
                    a2_full = get_ranking(df, n, g1, "A2")
                    a4r = a4_full.get(lor4d, 99)
                    a2r = a2_full.get(lor4d, 99)

                    results.append({
                        "n": n, "g1": g1, "g2": g2, "K": K,
                        "lor4d_a2_rank(g1)": lor4d_r1,
                        "in_topk": in_topk,
                        "lor4d_seq_rank": seq_rank,
                        "lor4d_a4_rank(g2)": a4r,
                        "lor4d_a2_rank": a2r,
                        "improvement": (a4r - seq_rank) if (seq_rank is not None) else None,
                        "top_k": top_k,
                    })

    rdf = pd.DataFrame(results)

    # 找所有序贯胜出条件
    wins = rdf[(rdf["lor4d_seq_rank"].notna()) &
               (rdf["lor4d_seq_rank"] == 1)].copy()
    wins = wins.sort_values(["n", "g1", "g2", "K"])

    print(f"\n总条件数: {len(rdf)}, 其中 Lor4D 进入 TOP-K: {rdf['in_topk'].sum()}")
    print(f"\n序贯筛选使 Lor4D = #1 的条件数: {len(wins)}")
    print()

    if len(wins) > 0:
        print(f"{'n':>3} | {'γ₁':>4} | {'γ₂':>4} | {'K':>3} | "
              f"{'A2(γ₁)':>6} | {'Seq':>3} | {'A4(γ₂)':>6} | {'Δ':>3} | TOP-K")
        print("-" * 110)
        for _, r in wins.iterrows():
            a4r = r["lor4d_a4_rank(g2)"]
            delta = int(r["improvement"])
            tk = r["top_k"]
            # 缩短 top_k 显示
            short = [f.replace("lorentzian_like_", "Lor").replace("_2layer", "2L")
                     .replace("_4layer", "4L").replace("_like", "~")
                     .replace("multi_layer_random", "MLR").replace("transitive_percolation", "TP")
                     .replace("random_layered_k", "RL").replace("_uniform", "")
                     .replace("absolute_layered", "AbsL").replace("interval_order", "IntO")
                     for f in tk]
            mark = " ★" if delta > 0 else ""
            print(f"{int(r['n']):>3} | {r['g1']:>4.1f} | {r['g2']:>4.1f} | {int(r['K']):>3} | "
                  f"{int(r['lor4d_a2_rank(g1)']):>6} | {int(r['lor4d_seq_rank']):>3} | "
                  f"{int(a4r):>6} | {delta:>+3} | {short}{mark}")

    # 汇总：最佳序贯配置 vs 最佳单模态
    print("\n" + "=" * 110)
    print("各 n 的最佳策略比较")
    print("=" * 110)
    for n in [10, 20, 40]:
        sub = rdf[rdf["n"] == n]
        # 最佳 A4-only
        a4_rows = sub.drop_duplicates(subset=["g2"])
        best_a4 = a4_rows.loc[a4_rows["lor4d_a4_rank(g2)"].idxmin()]
        # 最佳 sequential
        seq_valid = sub[sub["lor4d_seq_rank"].notna()]
        if len(seq_valid) > 0:
            best_seq = seq_valid.loc[seq_valid["lor4d_seq_rank"].idxmin()]
            print(f"n={n}: Best A4-only = #{int(best_a4['lor4d_a4_rank(g2)'])} "
                  f"(γ₂={best_a4['g2']}), "
                  f"Best Sequential = #{int(best_seq['lor4d_seq_rank'])} "
                  f"(γ₁={best_seq['g1']}, γ₂={best_seq['g2']}, K={int(best_seq['K'])})")
        else:
            print(f"n={n}: Best A4-only = #{int(best_a4['lor4d_a4_rank(g2)'])} "
                  f"(γ₂={best_a4['g2']}), Sequential: N/A (Lor4D never in TOP-K)")

    # 保存
    out_path = OUT / "sequential_screening_enhanced.csv"
    rdf.to_csv(out_path, index=False)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
