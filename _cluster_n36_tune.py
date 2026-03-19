"""Tuning cluster moves specifically for N=36 — the critical size."""
import numpy as np
from generators import generate_lorentzian_like_2d
from unified_functional_metropolis import (
    microcanonical_metropolis_sample, CALIBRATED_WEIGHTS,
    compute_F_fast, compute_observables_snapshot,
)
from _cluster_move import cluster_metropolis_sample

def tune_n36():
    N = 36
    n_reps = 5
    T = 2.0
    seed_base = 42
    
    configs = [
        {"label": "single-300",   "method": "single",  "steps": 300, "cf": 0.0, "k": 0},
        {"label": "single-600",   "method": "single",  "steps": 600, "cf": 0.0, "k": 0},
        {"label": "clust30-300",  "method": "cluster", "steps": 300, "cf": 0.3, "k": 3},
        {"label": "clust30-600",  "method": "cluster", "steps": 600, "cf": 0.3, "k": 3},
        {"label": "clust50-300",  "method": "cluster", "steps": 300, "cf": 0.5, "k": 3},
        {"label": "clust50-600",  "method": "cluster", "steps": 600, "cf": 0.5, "k": 3},
        {"label": "clust50-600-k5", "method": "cluster", "steps": 600, "cf": 0.5, "k": 5},
    ]
    
    print("="*80)
    print(f"N=36 CLUSTER MOVE TUNING (T={T}, reps={n_reps})")
    print("="*80)
    
    results = {}
    
    for cfg in configs:
        label = cfg["label"]
        wins = []
        comps = []
        dFs = []
        accepts = []
        
        for rep in range(n_reps):
            s = seed_base + rep * 777
            init = generate_lorentzian_like_2d(N, seed=s)
            init_F = compute_F_fast(init, CALIBRATED_WEIGHTS)
            
            if cfg["method"] == "single":
                res = microcanonical_metropolis_sample(
                    init, CALIBRATED_WEIGHTS, n_steps=cfg["steps"],
                    temperature=T, seed=s,
                )
            else:
                res = cluster_metropolis_sample(
                    init, CALIBRATED_WEIGHTS, n_steps=cfg["steps"],
                    temperature=T, seed=s,
                    cluster_fraction=cfg["cf"], multi_k=cfg["k"],
                )
            
            f = res.trajectory[-1]
            in_win = 0.30 <= f['comp_frac'] <= 0.55
            wins.append(in_win)
            comps.append(f['comp_frac'])
            dFs.append(f['F'] - init_F)
            accepts.append(res.acceptance_rate)
        
        results[label] = {
            "win_rate": np.mean(wins),
            "comp_mean": np.mean(comps),
            "comp_std": np.std(comps),
            "dF_mean": np.mean(dFs),
            "accept": np.mean(accepts),
        }
        
        r = results[label]
        print(f"  {label:>20s}: win={r['win_rate']:4.0%}, "
              f"comp={r['comp_mean']:.3f}±{r['comp_std']:.3f}, "
              f"ΔF={r['dF_mean']:+.1f}, accept={r['accept']:.3f}")
    
    print(f"\n  BEST config for N=36 window hit rate:")
    best = max(results.items(), key=lambda x: (x[1]['win_rate'], -x[1]['comp_std']))
    print(f"    {best[0]}: {best[1]['win_rate']:.0%} window hit")

if __name__ == "__main__":
    tune_n36()
