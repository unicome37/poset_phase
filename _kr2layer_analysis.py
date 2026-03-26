"""Deep analysis of KR_2layer rank inversion at large gamma.

Decomposes the score into logH vs penalty contributions to explain
why KR_2layer rises to rank #2 at gamma=0.8 despite having the worst
geometry among all 17 families.
"""
import pandas as pd
import numpy as np

df = pd.read_csv("outputs/raw_samples.csv")

# Focus on A2
sub = df[df["action_mode"] == "A2"]

# 1. Score decomposition: score = logH - gamma * penalty_effective
print("=" * 90)
print("ANALYSIS 1: Score decomposition for KR_2layer vs Lor2D vs Lor5D (N=80)")
print("=" * 90)
focus = ["KR_2layer", "lorentzian_like_2d", "lorentzian_like_5d", "KR_like"]
for fam in focus:
    print(f"\n  {fam}:")
    fam_sub = sub[(sub["family"] == fam) & (sub["n"] == 80)]
    agg = fam_sub.groupby("gamma").agg(
        mean_logH=("log_H_mean", "mean"),
        mean_penalty_eff=("penalty_effective", "mean"),
        mean_geo=("penalty_geometric", "mean"),
        mean_neutral=("penalty_neutral", "mean"),
        mean_score=("score", "mean"),
        mean_score_norm=("score_norm", "mean"),
    ).reset_index()
    for _, r in agg.iterrows():
        logH_contrib = r["mean_logH"]
        pen_contrib = r["gamma"] * r["mean_penalty_eff"]
        print(f"    gamma={r['gamma']:.1f}: logH={r['mean_logH']:.1f}, "
              f"penalty_eff={r['mean_penalty_eff']:.1f}, "
              f"gamma*pen={pen_contrib:.1f}, "
              f"score={r['mean_score']:.1f}, "
              f"score_norm={r['mean_score_norm']:+.3f}")

# 2. Why normalization amplifies KR_2layer at large gamma
print("\n" + "=" * 90)
print("ANALYSIS 2: Score distribution statistics per gamma (N=80)")
print("=" * 90)
for g in [0.0, 0.2, 0.4, 0.8]:
    g_sub = sub[(sub["gamma"] == g) & (sub["n"] == 80)]
    agg = g_sub.groupby("family")["score"].mean()
    print(f"\n  gamma={g}: score range [{agg.min():.1f}, {agg.max():.1f}], "
          f"spread={agg.max()-agg.min():.1f}, "
          f"median={agg.median():.1f}, "
          f"std={agg.std():.1f}")
    # KR_2layer and Lor2D positions
    kr2 = agg.get("KR_2layer", np.nan)
    lor2 = agg.get("lorentzian_like_2d", np.nan)
    print(f"    KR_2layer score={kr2:.1f}, Lor2D score={lor2:.1f}, "
          f"KR_2layer-Lor2D={kr2-lor2:.1f}")

# 3. Geometric penalty breakdown for KR_2layer
print("\n" + "=" * 90)
print("ANALYSIS 3: KR_2layer geometric sub-penalties vs Lor2D (N=80, gamma=0)")
print("=" * 90)
from generators import generate_kr_2layer, generate_lorentzian_like_2d
from observables_geo import geometric_components, DEFAULT_GEOMETRIC_WEIGHTS

weighted_keys = list(DEFAULT_GEOMETRIC_WEIGHTS.keys())

for fam_name, gen_fn in [("KR_2layer", generate_kr_2layer),
                          ("Lor2D", generate_lorentzian_like_2d)]:
    comp_vals = {k: [] for k in weighted_keys}
    for s in range(16):
        p = gen_fn(n=80, seed=80000 + s)
        comps = geometric_components(p)
        for k in weighted_keys:
            comp_vals[k].append(comps[k])
    print(f"\n  {fam_name} (N=80, 16 samples):")
    total_weighted = 0.0
    for k in weighted_keys:
        arr = np.array(comp_vals[k])
        w = DEFAULT_GEOMETRIC_WEIGHTS[k]
        contrib = arr.mean() * w
        total_weighted += contrib
        print(f"    {k:30s}: raw={arr.mean():.4f} (std={arr.std():.4f}), "
              f"weight={w:.0f}, weighted={contrib:.3f}")
    print(f"    {'TOTAL':30s}: {total_weighted:.3f}")

# 4. N-scaling of the crossover gamma
print("\n" + "=" * 90)
print("ANALYSIS 4: At what gamma does KR_2layer first beat Lor2D? (per N)")
print("=" * 90)
for n_val in [20, 40, 60, 80]:
    n_sub = sub[sub["n"] == n_val]
    agg = n_sub.groupby(["gamma", "family"])["score_norm"].mean().reset_index()
    for g in sorted(agg["gamma"].unique()):
        kr2 = agg[(agg["gamma"] == g) & (agg["family"] == "KR_2layer")]["score_norm"].values[0]
        lor2 = agg[(agg["gamma"] == g) & (agg["family"] == "lorentzian_like_2d")]["score_norm"].values[0]
        if kr2 > lor2:
            print(f"  N={n_val}: KR_2layer beats Lor2D at gamma={g} "
                  f"(KR2={kr2:+.3f} vs Lor2D={lor2:+.3f})")
            break
    else:
        print(f"  N={n_val}: KR_2layer never beats Lor2D in the tested gamma range")
