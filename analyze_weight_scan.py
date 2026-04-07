"""权重扫描结果深度分析

读取 weight_scan_summary.csv，生成：
1. Lor4D 排名热力图（config × gamma）
2. 最优权重配比识别
3. 权重敏感性梯度
4. 与理论预测的交叉验证
"""
import pandas as pd
import numpy as np

CSV = "outputs_info_weight_scan/weight_scan_summary.csv"

def main():
    df = pd.read_csv(CSV)
    print(f"Total records: {len(df)}, configs: {df['config'].nunique()}")
    print(f"N values: {sorted(df['n'].unique())}")
    print(f"Gamma values: {sorted(df['gamma'].unique())}")
    
    # ── 1. 排名热力图: config × (n, gamma) ──
    print("\n" + "="*100)
    print("1. Lor4D RANK HEATMAP (lower = better)")
    print("="*100)
    
    for n in sorted(df['n'].unique()):
        sub = df[df['n']==n]
        pivot = sub.pivot_table(index="config", columns="gamma", 
                                values="lor4d_rank", aggfunc="first")
        pivot["mean"] = pivot.mean(axis=1)
        pivot = pivot.sort_values("mean")
        print(f"\nn={n}:")
        print(pivot.to_string(float_format=lambda x: f"{x:.0f}"))
    
    # ── 2. 按 gamma 区间统计最优配比 ──
    print("\n" + "="*100)
    print("2. BEST CONFIG PER (n, gamma)")
    print("="*100)
    
    for (n, gamma), grp in df.groupby(["n", "gamma"]):
        best = grp.loc[grp["lor4d_rank"].idxmin()]
        worst = grp.loc[grp["lor4d_rank"].idxmax()]
        print(f"  n={n:3d} γ={gamma:.1f}: best={best['config']:30s} rank={best['lor4d_rank']:.0f} "
              f"| worst={worst['config']:30s} rank={worst['lor4d_rank']:.0f}")
    
    # ── 3. Interval weight 扫描曲线 ──
    print("\n" + "="*100)
    print("3. INTERVAL WEIGHT SWEEP (rest=5.0)")
    print("="*100)
    
    interval_configs = df[df['config'].str.startswith('interval_') | (df['config']=='baseline_5x5')]
    interval_configs = interval_configs.copy()
    # baseline has interval_weight=5.0
    
    for n in sorted(interval_configs['n'].unique()):
        for gamma in [0.2, 0.5, 1.0, 2.0]:
            sub = interval_configs[(interval_configs['n']==n) & (interval_configs['gamma']==gamma)]
            sub = sub.sort_values('interval_weight')
            pairs = [(r['interval_weight'], r['lor4d_rank']) for _, r in sub.iterrows()]
            curve = " → ".join([f"w={w:.1f}:#{int(r)}" for w, r in pairs])
            print(f"  n={n:3d} γ={gamma:.1f}: {curve}")
    
    # ── 4. 统一缩放效应 ──
    print("\n" + "="*100)
    print("4. UNIFORM SCALING EFFECT")
    print("="*100)
    
    uniform = df[df['config'].str.startswith('uniform_') | (df['config']=='baseline_5x5')]
    for n in sorted(uniform['n'].unique()):
        for gamma in [0.2, 0.5, 1.0]:
            sub = uniform[(uniform['n']==n) & (uniform['gamma']==gamma)]
            sub = sub.sort_values('interval_weight')
            for _, r in sub.iterrows():
                print(f"  n={n:3d} γ={gamma:.1f}: {r['config']:20s} rank={r['lor4d_rank']:.0f} "
                      f"gap={r['gap_to_top1']:.2f} top1={r['top1_family']}")
    
    # ── 5. 关键发现总结 ──
    print("\n" + "="*100)
    print("5. KEY FINDINGS")
    print("="*100)
    
    # 5a. Which config has Lor4D #1 at most (n, gamma) combos?
    best_counts = {}
    for config in df['config'].unique():
        sub = df[df['config']==config]
        n_best = (sub['lor4d_rank'] == 1).sum()
        best_counts[config] = n_best
    best_sorted = sorted(best_counts.items(), key=lambda x: -x[1])
    print("\nConfigs ranked by # of (n,gamma) combos where Lor4D = #1:")
    for config, count in best_sorted[:10]:
        total = df[df['config']==config].shape[0]
        print(f"  {config:35s}: {count}/{total}")
    
    # 5b. Robust ranking: mean rank across all (n, gamma)
    mean_ranks = df.groupby("config")["lor4d_rank"].mean().sort_values()
    print("\nMean Lor4D rank across all (n, gamma):")
    for config, mr in mean_ranks.items():
        print(f"  {config:35s}: {mr:.2f}")
    
    # 5c. Theoretical prediction validation
    print("\n--- Theoretical Prediction Validation ---")
    print("Theory predicted: baseline is optimal at n=20 γ=1.0, Lor4D=#1")
    pred = df[(df['config']=='baseline_5x5') & (df['n']==20) & (df['gamma']==1.0)]
    if len(pred) > 0:
        actual = pred.iloc[0]['lor4d_rank']
        print(f"Actual: Lor4D rank = #{actual:.0f}  {'✓ CONFIRMED' if actual <= 2 else '✗ REFUTED'}")
    
    print("\nTheory predicted: interval=0 hurts Lor4D (KR benefits more)")
    pred2 = df[(df['config']=='interval_0.0') & (df['n']==20) & (df['gamma']==1.0)]
    base2 = df[(df['config']=='baseline_5x5') & (df['n']==20) & (df['gamma']==1.0)]
    if len(pred2) > 0 and len(base2) > 0:
        rank_int0 = pred2.iloc[0]['lor4d_rank']
        rank_base = base2.iloc[0]['lor4d_rank']
        print(f"Actual: interval=0 rank={rank_int0:.0f}, baseline rank={rank_base:.0f}  "
              f"{'✓ CONFIRMED' if rank_int0 > rank_base else '✗ REFUTED'}")

if __name__ == "__main__":
    main()
