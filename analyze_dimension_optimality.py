"""Dimensional analysis: why d=4 is info-theoretically optimal."""
import numpy as np
from generators import (
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)
from entropy_sis import estimate_log_linear_extensions_sis
from observables_info import info_components

gens = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
}

for n in [10, 20, 40]:
    print(f"\n=== n={n} ===")
    header = f"{'Family':>8s}  {'logH':>8s}  {'I_info':>8s}  {'A_0.5':>10s}  {'P1_spec':>8s}  {'P2_deg':>8s}  {'P3_lay':>8s}  {'P4_edge':>8s}  {'P5_int':>8s}"
    print(header)
    for name, gen in gens.items():
        scores = []
        components_list = []
        for seed in range(4):
            p = gen(n=n, seed=seed * 100)
            lh, _ = estimate_log_linear_extensions_sis(p, n_runs=32, seed=seed)
            c = info_components(p)
            scores.append(lh)
            components_list.append(c)
        lh_mean = np.mean(scores)
        info_mean = np.mean([c["info_total"] for c in components_list])
        action_05 = -1.0 * lh_mean + 0.5 * info_mean
        p1 = np.mean([c["info_spectral_entropy_deficit"] for c in components_list])
        p2 = np.mean([c["info_degree_heterogeneity"] for c in components_list])
        p3 = np.mean([c["info_layer_concentration"] for c in components_list])
        p4 = np.mean([c["info_edge_density_extremity"] for c in components_list])
        p5 = np.mean([c["info_interval_diversity_deficit"] for c in components_list])
        print(
            f"{name:>8s}  {lh_mean:8.2f}  {info_mean:8.3f}  {action_05:10.3f}"
            f"  {p1:8.4f}  {p2:8.4f}  {p3:8.4f}  {p4:8.4f}  {p5:8.4f}"
        )
