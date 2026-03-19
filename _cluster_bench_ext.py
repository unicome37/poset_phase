"""Extended cluster move benchmark: more reps, include N=64."""
import numpy as np
from generators import Poset, generate_lorentzian_like_2d
from unified_functional_metropolis import (
    microcanonical_metropolis_sample, CALIBRATED_WEIGHTS,
    compute_F_fast, compute_observables_snapshot,
)
from _cluster_move import cluster_metropolis_sample

def run_extended():
    n_values = [16, 28, 36]
    n_reps = 5
    steps = 300
    T = 2.0
    seed_base = 42
    
    print("="*75)
    print("EXTENDED CLUSTER MOVE BENCHMARK")
    print(f"  Steps={steps}, T={T}, Reps={n_reps}")
    print("="*75)
    
    all_results = []
    
    for n in n_values:
        for rep in range(n_reps):
            s = seed_base + rep * 777 + n * 13
            init = generate_lorentzian_like_2d(n, seed=s)
            init_F = compute_F_fast(init, CALIBRATED_WEIGHTS)
            init_obs = compute_observables_snapshot(init)
            
            # Single-swap
            res_s = microcanonical_metropolis_sample(
                init, CALIBRATED_WEIGHTS, n_steps=steps,
                temperature=T, seed=s, label=f"s_N{n}_r{rep}",
            )
            fs = res_s.trajectory[-1]
            
            # Cluster
            res_c = cluster_metropolis_sample(
                init, CALIBRATED_WEIGHTS, n_steps=steps,
                temperature=T, seed=s, label=f"c_N{n}_r{rep}",
                cluster_fraction=0.3, multi_k=3,
            )
            fc = res_c.trajectory[-1]
            
            all_results.append({
                "N": n, "rep": rep, "method": "single",
                "accept": res_s.acceptance_rate,
                "comp": fs['comp_frac'], "F": fs['F'],
                "dF": fs['F'] - init_F,
                "in_win": 0.30 <= fs['comp_frac'] <= 0.55,
                "layers": fs['n_layers'], "d_eff": fs['d_eff'],
            })
            all_results.append({
                "N": n, "rep": rep, "method": "cluster",
                "accept": res_c.acceptance_rate,
                "comp": fc['comp_frac'], "F": fc['F'],
                "dF": fc['F'] - init_F,
                "in_win": 0.30 <= fc['comp_frac'] <= 0.55,
                "layers": fc['n_layers'], "d_eff": fc['d_eff'],
            })
            
            print(f"  N={n:2d} rep={rep}: single comp={fs['comp_frac']:.3f} "
                  f"win={'✓' if 0.30<=fs['comp_frac']<=0.55 else '✗'} | "
                  f"cluster comp={fc['comp_frac']:.3f} "
                  f"win={'✓' if 0.30<=fc['comp_frac']<=0.55 else '✗'}")
    
    # Summary table
    print(f"\n{'='*75}")
    print("SUMMARY (mean ± std over reps)")
    print(f"{'='*75}")
    print(f"  {'N':>3s} {'method':>8s} {'accept':>10s} {'comp':>12s} "
          f"{'ΔF':>12s} {'win%':>6s} {'layers':>10s}")
    print(f"  {'─'*3} {'─'*8} {'─'*10} {'─'*12} {'─'*12} {'─'*6} {'─'*10}")
    
    for n in n_values:
        for method in ["single", "cluster"]:
            rows = [r for r in all_results if r["N"]==n and r["method"]==method]
            acc = np.mean([r["accept"] for r in rows])
            comp_m = np.mean([r["comp"] for r in rows])
            comp_s = np.std([r["comp"] for r in rows])
            dF_m = np.mean([r["dF"] for r in rows])
            dF_s = np.std([r["dF"] for r in rows])
            win = np.mean([r["in_win"] for r in rows]) * 100
            lay_m = np.mean([r["layers"] for r in rows])
            lay_s = np.std([r["layers"] for r in rows])
            
            print(f"  {n:3d} {method:>8s} {acc:10.3f} "
                  f"{comp_m:.3f}±{comp_s:.3f} "
                  f"{dF_m:+.1f}±{dF_s:.1f} {win:5.0f}% "
                  f"{lay_m:.1f}±{lay_s:.1f}")
    
    # Improvement summary
    print(f"\n  KEY METRIC: Window hit rate improvement at each N:")
    for n in n_values:
        s_win = np.mean([r["in_win"] for r in all_results if r["N"]==n and r["method"]=="single"]) * 100
        c_win = np.mean([r["in_win"] for r in all_results if r["N"]==n and r["method"]=="cluster"]) * 100
        delta = c_win - s_win
        print(f"    N={n:2d}: single={s_win:.0f}% → cluster={c_win:.0f}% (Δ={delta:+.0f}%)")

if __name__ == "__main__":
    run_extended()
