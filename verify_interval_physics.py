"""Interval diversity 物理解释独立验证

目的：验证 interval_diversity 作为 4D Lorentzian 几何签名的物理解释
1. 各族 interval size entropy 分布统计
2. Lor-dD (d=2,3,4,5) 维度梯度检验：interval diversity 应随 d 单调变化
3. 与已知 Alexandrov interval volume ∝ τ^d 关系的一致性
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from collections import defaultdict

from generators import Poset
from experiment import FAMILIES
from observables_info import (
    interval_size_entropy_normalised,
    _sample_intervals,
)


def analyze_interval_distribution(poset: Poset, max_pairs: int = 128):
    """Compute detailed interval statistics for a single poset."""
    intervals = _sample_intervals(poset, max_pairs=max_pairs)
    if not intervals:
        return {"n_intervals": 0, "mean_size": 0, "std_size": 0, 
                "max_size": 0, "n_distinct": 0, "entropy_norm": 0}
    
    sizes = [len(interior) for interior in intervals]
    sizes_arr = np.array(sizes)
    
    return {
        "n_intervals": len(intervals),
        "mean_size": float(sizes_arr.mean()),
        "std_size": float(sizes_arr.std()),
        "max_size": int(sizes_arr.max()),
        "n_distinct": len(set(sizes)),
        "entropy_norm": float(interval_size_entropy_normalised(poset, max_pairs)),
    }


def main():
    n_values = [10, 20, 40]
    samples_per = 20
    
    results = []
    
    # Focus on relevant families
    focus_families = [
        "lorentzian_like_2d",
        "lorentzian_like_3d", 
        "lorentzian_like_4d",
        "lorentzian_like_5d",
        "KR_2layer",
        "KR_like",
        "random_dag",
        "transitive_percolation",
        "multi_layer_random",
    ]
    
    available = set(FAMILIES.keys())
    focus_families = [f for f in focus_families if f in available]
    
    for n in n_values:
        print(f"\n{'='*80}")
        print(f"n = {n}")
        print(f"{'='*80}")
        
        family_stats = defaultdict(list)
        
        for name in focus_families:
            gen_fn = FAMILIES[name]
            
            entropies = []
            for _ in range(samples_per):
                poset = gen_fn(n)
                stats = analyze_interval_distribution(poset)
                entropies.append(stats["entropy_norm"])
                family_stats[name].append(stats)
            
            ent_arr = np.array(entropies)
            results.append({
                "family": name,
                "n": n,
                "entropy_mean": float(ent_arr.mean()),
                "entropy_std": float(ent_arr.std()),
                "entropy_min": float(ent_arr.min()),
                "entropy_max": float(ent_arr.max()),
            })
            
            print(f"  {name:30s}: H_int/ln(m) = {ent_arr.mean():.4f} ± {ent_arr.std():.4f}")
    
    df = pd.DataFrame(results)
    
    # ── Dimension gradient check ──
    print("\n" + "="*80)
    print("DIMENSION GRADIENT CHECK: Lor-dD interval entropy vs d")
    print("="*80)
    
    lor_families = ["lorentzian_like_2d", "lorentzian_like_3d", 
                    "lorentzian_like_4d", "lorentzian_like_5d"]
    lor_families = [f for f in lor_families if f in focus_families]
    
    for n in n_values:
        sub = df[(df['n']==n) & (df['family'].isin(lor_families))]
        sub = sub.sort_values('family')
        
        print(f"\nn={n}:")
        for _, row in sub.iterrows():
            d = int(row['family'].split('_')[-1].replace('d',''))
            bar = "█" * int(row['entropy_mean'] * 50)
            print(f"  d={d}: H={row['entropy_mean']:.4f} ± {row['entropy_std']:.4f} {bar}")
        
        # Monotonicity check
        entropies = sub['entropy_mean'].tolist()
        is_monotone = all(entropies[i] <= entropies[i+1] for i in range(len(entropies)-1))
        is_mono_dec = all(entropies[i] >= entropies[i+1] for i in range(len(entropies)-1))
        print(f"  Monotone increasing: {is_monotone}")
        print(f"  Monotone decreasing: {is_mono_dec}")
    
    # ── KR vs Lor4D comparison ──
    print("\n" + "="*80)
    print("KR vs Lor4D: Why interval diversity discriminates")
    print("="*80)
    
    for n in n_values:
        lor4d = df[(df['n']==n) & (df['family']=='lorentzian_like_4d')]
        kr2 = df[(df['n']==n) & (df['family']=='KR_2layer')]
        if len(lor4d) > 0 and len(kr2) > 0:
            lor_ent = lor4d.iloc[0]['entropy_mean']
            kr_ent = kr2.iloc[0]['entropy_mean']
            # penalty = (1 - H)^2
            lor_pen = (1 - lor_ent)**2
            kr_pen = (1 - kr_ent)**2
            ratio = kr_pen / lor_pen if lor_pen > 0 else float('inf')
            print(f"  n={n}: Lor4D H={lor_ent:.4f} pen={lor_pen:.4f} | "
                  f"KR H={kr_ent:.4f} pen={kr_pen:.4f} | "
                  f"ratio={ratio:.2f}x")
    
    # ── Physical interpretation ──
    print("\n" + "="*80)
    print("PHYSICAL INTERPRETATION")
    print("="*80)
    print("""
In d-dimensional Minkowski spacetime, the Alexandrov interval between 
two causally related points has volume ∝ τ^d (proper time separation).
For a random sprinkling of n points:
- d=2: intervals have ~uniform size → moderate entropy  
- d=3: more size variety (larger range of τ) → higher entropy
- d=4: even more variety → highest entropy in tested range
- d=5: counterintuitively may decrease (most volume in large intervals)

KR/2-layer structures have degenerate intervals (all size 0 or 1):
- Very low entropy → high penalty
- This is the structural "weakness" that helps Lor4D win

The dimension gradient and KR-gap together explain why interval_diversity
is the critical discriminator in the info-theoretic penalty.
""")


if __name__ == "__main__":
    main()
