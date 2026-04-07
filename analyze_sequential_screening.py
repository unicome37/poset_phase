#!/usr/bin/env python3
"""β-1: 序贯筛选模拟 — 从已有 hybrid 数据分析 Step1(geo)→TOP-K→Step2(info) 效果。

原理：
  Step 1 (A2: neutral+geo): 排除非几何族，保留 TOP-K
  Step 2 (A4: pure info):   在 TOP-K 中精选 Lor4D
  比较：序贯 vs A2-only vs A4-only vs A6/A7-hybrid
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


def rank_families(df_scores: pd.DataFrame, n: int, gamma: float,
                  mode: str, family_subset: list | None = None) -> pd.DataFrame:
    """给定条件下对族进行排名（score 越低越好）。"""
    mask = (df_scores["n"] == n) & (df_scores["gamma"] == gamma) & (df_scores["action_mode"] == mode)
    sub = df_scores[mask].copy()
    if family_subset is not None:
        sub = sub[sub["family"].isin(family_subset)]
    sub = sub.sort_values("score")
    sub["rank"] = range(1, len(sub) + 1)
    return sub


def sequential_screening(df_scores: pd.DataFrame, n: int, gamma: float,
                         K: int) -> dict:
    """执行序贯筛选并返回结果。"""
    # Step 1: A2 排名
    a2_ranked = rank_families(df_scores, n, gamma, "A2")
    top_k_families = a2_ranked.head(K)["family"].tolist()

    # 检查 Lor4D 是否进入 TOP-K
    lor4d_in_topk = "lorentzian_like_4d" in top_k_families
    lor4d_a2_rank = int(a2_ranked[a2_ranked["family"] == "lorentzian_like_4d"]["rank"].iloc[0])

    # Step 2: A4 在 TOP-K 中重排名
    if lor4d_in_topk:
        a4_reranked = rank_families(df_scores, n, gamma, "A4",
                                    family_subset=top_k_families)
        lor4d_seq_rank = int(a4_reranked[a4_reranked["family"] == "lorentzian_like_4d"]["rank"].iloc[0])
        seq_top1 = a4_reranked.iloc[0]["family"]
    else:
        lor4d_seq_rank = None  # 未进入 TOP-K，直接淘汰
        a4_reranked = rank_families(df_scores, n, gamma, "A4",
                                    family_subset=top_k_families)
        seq_top1 = a4_reranked.iloc[0]["family"]

    # 对照：A4-only 排名
    a4_full = rank_families(df_scores, n, gamma, "A4")
    lor4d_a4_rank = int(a4_full[a4_full["family"] == "lorentzian_like_4d"]["rank"].iloc[0])

    # 对照：A6 排名
    a6_full = rank_families(df_scores, n, gamma, "A6")
    lor4d_a6_rank = int(a6_full[a6_full["family"] == "lorentzian_like_4d"]["rank"].iloc[0])

    return {
        "n": n, "gamma": gamma, "K": K,
        "lor4d_a2_rank": lor4d_a2_rank,
        "lor4d_in_topk": lor4d_in_topk,
        "lor4d_sequential_rank": lor4d_seq_rank,
        "lor4d_a4_only_rank": lor4d_a4_rank,
        "lor4d_a6_hybrid_rank": lor4d_a6_rank,
        "sequential_top1": seq_top1,
        "topk_families": top_k_families,
    }


def main():
    df_scores = load_family_scores(RAW)
    families_list = sorted(df_scores["family"].unique())
    n_families = len(families_list)
    print(f"Total families: {n_families}")
    print(f"Families: {families_list[:5]}...{families_list[-3:]}\n")

    results = []
    for n in [10, 20, 40]:
        for gamma in [0.0, 0.1, 0.2, 0.5, 1.0, 2.0]:
            for K in [3, 5, 8, 10]:
                r = sequential_screening(df_scores, n, gamma, K)
                results.append(r)

    df_results = pd.DataFrame(results)

    # 汇总表
    print("=" * 100)
    print("β-1: 序贯筛选 vs 单一模态 vs 混合 — Lor4D 排名对比")
    print("=" * 100)

    for n in [10, 20, 40]:
        print(f"\n--- n = {n} ---")
        print(f"{'γ':>5} | {'A2':>3} | {'A4':>3} | {'A6':>3} |"
              f" {'Seq K=3':>7} | {'Seq K=5':>7} | {'Seq K=8':>7} | {'Seq K=10':>8} |"
              f" Best")
        print("-" * 90)
        for gamma in [0.0, 0.1, 0.2, 0.5, 1.0, 2.0]:
            row = df_results[(df_results["n"] == n) & (df_results["gamma"] == gamma)]
            a2 = row[row["K"] == 3].iloc[0]["lor4d_a2_rank"]
            a4 = row[row["K"] == 3].iloc[0]["lor4d_a4_only_rank"]
            a6 = row[row["K"] == 3].iloc[0]["lor4d_a6_hybrid_rank"]

            seq_ranks = []
            for K in [3, 5, 8, 10]:
                r = row[row["K"] == K].iloc[0]
                sr = r["lor4d_sequential_rank"]
                seq_ranks.append(sr if sr is not None else "OUT")

            # 找最佳方法
            all_ranks = {"A2": a2, "A4": a4, "A6": a6}
            for i, K in enumerate([3, 5, 8, 10]):
                sr = seq_ranks[i]
                if sr != "OUT":
                    all_ranks[f"Seq-{K}"] = sr
            best = min(all_ranks.items(), key=lambda x: x[1] if isinstance(x[1], (int, float)) else 99)

            seq_strs = [f"{s:>7}" if isinstance(s, (int, float)) else f"{'OUT':>7}" for s in seq_ranks]
            print(f"{gamma:>5.1f} | {a2:>3} | {a4:>3} | {a6:>3} | "
                  f"{seq_strs[0]} | {seq_strs[1]} | {seq_strs[2]} | {seq_strs[3]:>8} | "
                  f"← {best[0]}={best[1]}")

    # 序贯筛选胜率统计
    print("\n" + "=" * 100)
    print("序贯筛选胜率统计（非零γ条件）")
    print("=" * 100)

    for K in [3, 5, 8, 10]:
        sub = df_results[(df_results["K"] == K) & (df_results["gamma"] > 0)]
        n_cond = len(sub)
        n_seq_best = 0
        n_seq_eq_best = 0
        n_lor4d_eliminated = 0
        for _, row in sub.iterrows():
            sr = row["lor4d_sequential_rank"]
            a4 = row["lor4d_a4_only_rank"]
            a6 = row["lor4d_a6_hybrid_rank"]
            a2 = row["lor4d_a2_rank"]
            if sr is None:
                n_lor4d_eliminated += 1
                continue
            best_single = min(a2, a4, a6)
            if sr < best_single:
                n_seq_best += 1
            elif sr == best_single:
                n_seq_eq_best += 1
        print(f"K={K:>2}: seq wins {n_seq_best}/{n_cond}, "
              f"ties {n_seq_eq_best}/{n_cond}, "
              f"Lor4D eliminated {n_lor4d_eliminated}/{n_cond}")

    # 详细：哪些条件下序贯最优
    print("\n" + "=" * 100)
    print("序贯筛选独胜条件明细（Seq < min(A2,A4,A6)）")
    print("=" * 100)
    for K in [5, 8]:
        sub = df_results[(df_results["K"] == K) & (df_results["gamma"] > 0)]
        for _, row in sub.iterrows():
            sr = row["lor4d_sequential_rank"]
            if sr is None:
                continue
            best_single = min(row["lor4d_a2_rank"], row["lor4d_a4_only_rank"],
                              row["lor4d_a6_hybrid_rank"])
            if sr < best_single:
                print(f"  K={K}, n={row['n']}, γ={row['gamma']}: "
                      f"Seq={sr}, A2={row['lor4d_a2_rank']}, "
                      f"A4={row['lor4d_a4_only_rank']}, A6={row['lor4d_a6_hybrid_rank']}"
                      f"  — TOP-K: {row['topk_families']}")

    # 保存完整结果
    out_path = OUT / "sequential_screening_analysis.csv"
    df_results.to_csv(out_path, index=False)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
