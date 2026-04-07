"""信息论+几何混合惩罚实验 (Direction β)

核心假设：
- 几何惩罚选维度(d=4)但不选族类型
- 信息论惩罚选族类型(Lorentzian)但不选维度
- 二者联合应在更宽 (n, γ) 范围内唯一选出 Lor4D

实验设计：
- A6 模式：geo + info（无 neutral）
- A7 模式：neutral + geo + info
- 对比基线：A2(neutral+geo), A4(pure info)
- n ∈ {10, 20, 40}
- γ ∈ {0, 0.1, 0.2, 0.5, 1.0, 2.0}（加入 γ=0.1 精细探测低惩罚区）
- 10 samples/family, 48 SIS runs
"""
from __future__ import annotations
from pathlib import Path

import pandas as pd
import numpy as np

from experiment import run_experiment
from normalization import add_normalized_columns, add_size_scaled_columns

OUTPUT_DIR = Path("outputs_info_hybrid")


def run_hybrid_experiment():
    """Run A2, A4, A6, A7 comparison experiment."""
    
    n_values = (10, 20, 40)
    gammas = (0.0, 0.1, 0.2, 0.5, 1.0, 2.0)
    modes = ("A2", "A4", "A6", "A7")
    
    print(f"Running hybrid experiment: {len(n_values)} N × {len(gammas)} γ × {len(modes)} modes")
    
    df = run_experiment(
        n_values=n_values,
        gammas=gammas,
        beta=1.0,
        samples_per_family=10,
        sis_runs=48,
        action_modes=modes,
    )
    
    df = add_size_scaled_columns(df)
    df = add_normalized_columns(
        df, method="robust_zscore",
        group_cols=("n", "gamma", "action_mode"),
    )
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_DIR / "hybrid_raw.csv", index=False, encoding="utf-8-sig")
    
    # ── Analysis: Lor4D rank per (mode, n, gamma) ──
    summaries = []
    for (mode, n, gamma), group in df.groupby(["action_mode", "n", "gamma"]):
        ranked = group.groupby("family")["score"].mean().sort_values()
        families_ranked = list(ranked.index)
        lor4d_rank = (
            families_ranked.index("lorentzian_like_4d") + 1
            if "lorentzian_like_4d" in families_ranked else -1
        )
        top1 = families_ranked[0] if families_ranked else "N/A"
        top1_score = ranked.iloc[0] if len(ranked) > 0 else np.nan
        lor4d_score = ranked.get("lorentzian_like_4d", np.nan)
        gap = lor4d_score - top1_score if lor4d_rank > 1 else 0.0
        
        summaries.append({
            "action_mode": mode,
            "n": n,
            "gamma": gamma,
            "lor4d_rank": lor4d_rank,
            "lor4d_score": lor4d_score,
            "top1_family": top1,
            "top1_score": top1_score,
            "gap_to_top1": gap,
            "total_families": len(families_ranked),
        })
    
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUTPUT_DIR / "hybrid_summary.csv", index=False, encoding="utf-8-sig")
    
    # ── Print results ──
    print("\n" + "="*100)
    print("HYBRID EXPERIMENT: Lor4D rank comparison across action modes")
    print("="*100)
    
    for n in sorted(summary['n'].unique()):
        pivot = summary[summary['n']==n].pivot_table(
            index="action_mode", columns="gamma",
            values="lor4d_rank", aggfunc="first",
        )
        pivot["mean"] = pivot.mean(axis=1)
        pivot = pivot.sort_values("mean")
        print(f"\nn={n}:")
        print(pivot.to_string(float_format=lambda x: f"{x:.0f}"))
    
    # ── Best mode per (n, gamma) ──
    print("\n" + "="*100)
    print("BEST MODE per (n, gamma)")
    print("="*100)
    
    for (n, gamma), grp in summary.groupby(["n", "gamma"]):
        best = grp.loc[grp["lor4d_rank"].idxmin()]
        print(f"  n={n:3d} γ={gamma:.1f}: {best['action_mode']} rank={best['lor4d_rank']:.0f} "
              f"top1={best['top1_family']}")
    
    # ── Hybrid advantage ──
    print("\n" + "="*100)
    print("HYBRID ADVANTAGE: A6/A7 vs max(A2, A4)")
    print("="*100)
    
    for (n, gamma), grp in summary.groupby(["n", "gamma"]):
        if gamma == 0.0:
            continue
        rows = {r['action_mode']: r['lor4d_rank'] for _, r in grp.iterrows()}
        a2 = rows.get('A2', 99)
        a4 = rows.get('A4', 99)
        a6 = rows.get('A6', 99)
        a7 = rows.get('A7', 99)
        best_single = min(a2, a4)
        best_hybrid = min(a6, a7)
        delta = best_single - best_hybrid  # positive = hybrid is better
        marker = "★" if delta > 0 else ("=" if delta == 0 else "▼")
        print(f"  n={n:3d} γ={gamma:.1f}: A2=#{int(a2)} A4=#{int(a4)} | "
              f"A6=#{int(a6)} A7=#{int(a7)} | "
              f"Δ={delta:+.0f} {marker}")


if __name__ == "__main__":
    run_hybrid_experiment()
