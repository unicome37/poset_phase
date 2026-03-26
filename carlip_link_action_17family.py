"""
Post-Carlip Experiment: BDG Link Action on All 17 Families
===========================================================
Replace logH with S_link = (N - 2C₀)/N as the primary ordering variable.

Motivation (Carlip critique):
- logH (linear extensions) has no established physical meaning in CST
- S_link = N - 2C₀ is the d=2 BD action, with established CST physics:
  * Loomis-Carlip 2018: link term dominates KR suppression
  * Carlip 2024: C₀ controls leading-order entropic suppression
  * Mathur-Singh-Surya 2020: link action suppresses non-manifoldlike sets

Three functional candidates tested:
  F_link1 = S_link/N                           (pure link action)
  F_link2 = S_link/N - λ·Σ_hist + η·ξ_d       (link + geometric terms)
  F_link3 = S_link/N + wall                     (link + wall)

Lower F = more favored. We want Lor4D to rank LOW (favored).
"""
from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from bd_action import count_intervals_fast, bdg_action_d2_link
from generators import (
    Poset,
    generate_kr_like,
    generate_kr_2layer,
    generate_kr_4layer,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_transitive_percolation,
    generate_interval_order,
    generate_absolute_layered,
    generate_multi_layer_random,
    generate_random_layered_k4_uniform,
    generate_random_layered_k6_uniform,
    generate_random_layered_k8_uniform,
    generate_random_layered_k6_tapered,
    generate_random_layered_k6_middle_heavy,
    generate_random_layered_k6_longjump,
)
from unified_functional import (
    compute_sigma_hist,
    compute_xi_dim,
)


def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_metrics(poset: Poset, N: int) -> dict:
    counts = count_intervals_fast(poset, k_max=3)
    C0 = counts.get(0)
    total_rel = counts.total_relations

    s_link = bdg_action_d2_link(counts, N, normalized=True)
    link_density = C0 / max(1, N * (N - 1) // 2)  # f_link = C₀ / C(N,2)

    R = 1.0 - float(C0) / float(total_rel) if total_rel > 0 else 0.0
    sigma_hist = compute_sigma_hist(poset)
    xi_dim_val, d_eff = compute_xi_dim(poset)

    # F_link1: pure link action (lower = more links = more constrained)
    F_link1 = s_link

    # F_link2: link action + geometric terms
    lam, eta = 0.1, 0.05
    F_link2 = s_link - lam * sigma_hist + eta * xi_dim_val

    # F_link3: link action + wall (penalize high-R = sparse structures)
    alpha_N = 2.0 * (20.0 / max(N, 1)) ** 0.5
    wall = alpha_N * sigmoid((R - 0.25) / 0.015)
    F_link3 = s_link + wall

    return {
        "C0": C0,
        "total_rel": total_rel,
        "s_link": s_link,
        "link_density": link_density,
        "R": R,
        "d_eff": d_eff,
        "sigma_hist": sigma_hist,
        "xi_dim": xi_dim_val,
        "F_link1": F_link1,
        "F_link2": F_link2,
        "F_link3": F_link3,
    }


# All 17 families
FAMILIES = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
    "KR_2layer": generate_kr_2layer,
    "KR_4layer": generate_kr_4layer,
    "AbsLayer": generate_absolute_layered,
    "MLR": generate_multi_layer_random,
    "RLk4": generate_random_layered_k4_uniform,
    "RLk6": generate_random_layered_k6_uniform,
    "RLk8": generate_random_layered_k8_uniform,
    "RLk6_tap": generate_random_layered_k6_tapered,
    "RLk6_mid": generate_random_layered_k6_middle_heavy,
    "RLk6_lj": generate_random_layered_k6_longjump,
    "TransPerc": generate_transitive_percolation,
    "IntOrder": generate_interval_order,
}

CATEGORY = {}
for f in ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]:
    CATEGORY[f] = "Lorentzian"
for f in ["KR_like", "KR_2layer", "KR_4layer"]:
    CATEGORY[f] = "KR-family"
for f in ["AbsLayer", "MLR", "RLk4", "RLk6", "RLk8", "RLk6_tap", "RLk6_mid", "RLk6_lj"]:
    CATEGORY[f] = "Layered"
for f in ["TransPerc", "IntOrder"]:
    CATEGORY[f] = "Other"


def main():
    N_VALUES = [16, 20, 28, 36, 48]
    REPS = 10
    SEED_BASE = 42

    print("=" * 80)
    print("POST-CARLIP: Link Action on All 17 Families")
    print("=" * 80)
    print(f"N = {N_VALUES}, reps = {REPS}")
    print()

    all_rows = []
    total = len(FAMILIES) * len(N_VALUES) * REPS
    done = 0

    for fam_name, gen_fn in FAMILIES.items():
        for N in N_VALUES:
            for rep in range(REPS):
                seed = SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep
                seed = seed % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    result = compute_metrics(poset, N)
                    result["family"] = fam_name
                    result["category"] = CATEGORY[fam_name]
                    result["N"] = N
                    result["rep"] = rep
                    all_rows.append(result)
                except Exception as e:
                    print(f"  ERROR: {fam_name} N={N} rep={rep}: {e}")

                done += 1
                if done % 100 == 0:
                    print(f"  [{done}/{total}] {fam_name} N={N}", flush=True)

    by_nf = defaultdict(list)
    for r in all_rows:
        by_nf[(r["N"], r["family"])].append(r)

    report = []
    report.append("# Post-Carlip: Link Action (S_link = N-2C₀) on All 17 Families\n")
    report.append(f"N = {N_VALUES}, reps = {REPS}\n")
    report.append("**Key question**: Does S_link/N (physically grounded in BDG) rank")
    report.append("Lorentzian structures above layered/KR controls?\n")

    for func_name, func_label in [
        ("F_link1", "F_link1 = S_link/N (pure link action)"),
        ("F_link2", "F_link2 = S_link/N - 0.1·Σ_hist + 0.05·ξ_d"),
        ("F_link3", "F_link3 = S_link/N + wall"),
    ]:
        report.append(f"\n## {func_label}\n")

        for N in N_VALUES:
            report.append(f"\n### N = {N}\n")
            report.append("| Rank | Family | Category | F | S_link/N | link_den | d_eff | R |")
            report.append("|------|--------|----------|---|---------|----------|-------|---|")

            means = {}
            details = {}
            for fam in FAMILIES:
                rows = by_nf.get((N, fam), [])
                if rows:
                    f_vals = [r[func_name] for r in rows]
                    means[fam] = np.mean(f_vals)
                    details[fam] = {
                        func_name: np.mean(f_vals),
                        "s_link": np.mean([r["s_link"] for r in rows]),
                        "link_density": np.mean([r["link_density"] for r in rows]),
                        "d_eff": np.mean([r["d_eff"] for r in rows]),
                        "R": np.mean([r["R"] for r in rows]),
                    }

            ranked = sorted(means, key=means.get)
            for rank, fam in enumerate(ranked, 1):
                d = details[fam]
                cat = CATEGORY[fam]
                tag = ""
                if cat == "Lorentzian":
                    tag = " ◆"
                elif cat == "KR-family":
                    tag = " ●"
                report.append(
                    f"| {rank} | {fam}{tag} | {cat} | {d[func_name]:+.4f} | "
                    f"{d['s_link']:.4f} | {d['link_density']:.4f} | "
                    f"{d['d_eff']:.2f} | {d['R']:.3f} |"
                )

        # Summary: does Lor4D beat all non-Lor?
        report.append(f"\n### {func_label} — Lor4D vs All Summary\n")
        for N in N_VALUES:
            lor4d_rows = by_nf.get((N, "Lor4D"), [])
            if not lor4d_rows:
                continue
            lor4d_vals = [r[func_name] for r in lor4d_rows]
            lor4d_mean = np.mean(lor4d_vals)

            non_lor = [f for f in FAMILIES if CATEGORY[f] != "Lorentzian"]
            beats_all = True
            losers = []
            for nf in non_lor:
                nr = by_nf.get((N, nf), [])
                if nr:
                    nf_mean = np.mean([r[func_name] for r in nr])
                    if lor4d_mean >= nf_mean:
                        beats_all = False
                        losers.append(f"{nf}({nf_mean:.4f})")

            if beats_all:
                report.append(f"- N={N}: **Lor4D beats ALL non-Lorentzian** ✅")
            else:
                report.append(f"- N={N}: Lor4D({lor4d_mean:.4f}) LOSES to: {', '.join(losers[:5])}")

    # Key diagnostic: link density by family
    report.append("\n## Diagnostic: Link Density (C₀/C(N,2)) by Family\n")
    report.append("Higher link density = more causal links = more 'manifold-like'\n")
    report.append("| Family | Category | N=20 | N=36 | N=48 |")
    report.append("|--------|----------|------|------|------|")
    for fam in sorted(FAMILIES.keys()):
        vals = []
        for N in [20, 36, 48]:
            rows = by_nf.get((N, fam), [])
            if rows:
                vals.append(f"{np.mean([r['link_density'] for r in rows]):.4f}")
            else:
                vals.append("—")
        report.append(f"| {fam} | {CATEGORY[fam]} | {' | '.join(vals)} |")

    # Conclusion
    report.append("\n## Conclusion\n")
    report.append("### Does replacing logH with S_link address Carlip's critique?\n")
    report.append("- **C1 (physical basis)**: S_link = N - 2C₀ is the d=2 BD action,")
    report.append("  established in CST literature as physically meaningful.")
    report.append("- **C2 (family selection)**: Test whether S_link ranks Lor above ALL 17 families.")
    report.append("- **C3 (literature)**: S_link directly connects to Loomis-Carlip 2018,")
    report.append("  Carlip 2024, Mathur-Singh-Surya 2020.\n")

    out = "\n".join(report)
    outdir = Path("outputs_carlip")
    outdir.mkdir(exist_ok=True)
    md_path = outdir / "link_action_17family_test.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"\nSaved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
