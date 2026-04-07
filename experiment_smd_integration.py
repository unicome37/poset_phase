#!/usr/bin/env python3
"""β-2: S_MD 路径积分整合实验 — exp(−β·logH + γ·P + α·S_MD) 联合作用量测试。

测试：将 Mahalanobis identity 层直接纳入 Boltzmann 权重是否改善 Lor4D 选择。
对比：
  - A4-only (info): A = −logH + γ·I_info
  - A4+S_MD:        A = −logH + γ·I_info + α·S_MD
  - S_MD-only:      A = −logH + α·S_MD
  - A2+S_MD:        A = −logH + γ·P_geo + α·S_MD
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Ensure project root on path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from experiment import FAMILIES
from entropy_sis import estimate_log_linear_extensions_sis
from action import get_action_penalty
from bd_action import count_intervals_fast
from expanded_family_robustness import compute_features, mahalanobis_score
from observables_info import weighted_info_total

OUT = Path("outputs_smd_integration")
OUT.mkdir(exist_ok=True)

# ───────────────────────────
# Configuration
# ───────────────────────────
N_VALUES = (10, 20, 40)
SAMPLES = 10
SIS_RUNS = 48
GAMMAS = (0.0, 0.2, 0.5, 1.0, 2.0)
ALPHAS = (0.0, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0)
MODES = ("info_only", "info_plus_smd", "smd_only", "geo_plus_smd")
LOR4D_REF_REPS = 40


def build_lor4d_reference(n: int, seed_base: int = 99000) -> tuple[np.ndarray, np.ndarray]:
    """Build Lor4D reference manifold (μ, Σ⁻¹) at given n."""
    gen = FAMILIES["lorentzian_like_4d"]
    feats = []
    for i in range(LOR4D_REF_REPS):
        p = gen(n=n, seed=seed_base + i)
        feats.append(compute_features(p, n))
    feat_mat = np.array(feats)
    mu = feat_mat.mean(axis=0)
    cov = np.cov(feat_mat.T) + 1e-8 * np.eye(3)
    cov_inv = np.linalg.inv(cov)
    return mu, cov_inv


def main():
    print("β-2: S_MD Path Integral Integration Experiment")
    print("=" * 60)

    rows = []
    total = len(N_VALUES) * len(FAMILIES) * SAMPLES
    count = 0

    for n in N_VALUES:
        # Pre-compute reference manifold for this n
        mu_ref, cov_inv_ref = build_lor4d_reference(n)
        print(f"\nn={n}: μ={mu_ref}, det(Σ⁻¹)={np.linalg.det(cov_inv_ref):.2e}")

        for family_name, generator in FAMILIES.items():
            for sid in range(SAMPLES):
                count += 1
                if count % 50 == 0:
                    print(f"  [{count}/{total}]")

                seed = 1000 * n + sid
                poset = generator(n=n, seed=seed)

                # Entropy
                log_h_mean, log_h_std = estimate_log_linear_extensions_sis(
                    poset, n_runs=SIS_RUNS, seed=seed
                )

                # Penalties
                penalty_info = get_action_penalty(poset, "A4")
                penalty_geo = get_action_penalty(poset, "A3")

                # S_MD
                feat = compute_features(poset, n)
                smd = mahalanobis_score(feat, mu_ref, cov_inv_ref)

                for gamma in GAMMAS:
                    for alpha in ALPHAS:
                        for mode in MODES:
                            if mode == "info_only":
                                score = -log_h_mean + gamma * penalty_info
                            elif mode == "info_plus_smd":
                                score = -log_h_mean + gamma * penalty_info + alpha * smd
                            elif mode == "smd_only":
                                score = -log_h_mean + alpha * smd
                            elif mode == "geo_plus_smd":
                                score = -log_h_mean + gamma * penalty_geo + alpha * smd
                            else:
                                raise ValueError(mode)

                            rows.append({
                                "n": n, "family": family_name,
                                "sample_id": sid, "gamma": gamma,
                                "alpha": alpha, "mode": mode,
                                "log_H_mean": log_h_mean,
                                "penalty_info": penalty_info,
                                "penalty_geo": penalty_geo,
                                "smd": smd, "score": score,
                            })

    df = pd.DataFrame(rows)
    raw_path = OUT / "smd_integration_raw.csv"
    df.to_csv(raw_path, index=False)
    print(f"\nRaw data: {raw_path} ({len(df)} rows)")

    # ─── Summary: family-average ranking ───
    summary_rows = []
    for n in N_VALUES:
        for gamma in GAMMAS:
            for alpha in ALPHAS:
                for mode in MODES:
                    sub = df[(df["n"] == n) & (df["gamma"] == gamma) &
                             (df["alpha"] == alpha) & (df["mode"] == mode)]
                    if sub.empty:
                        continue
                    fam_means = sub.groupby("family")["score"].mean().sort_values()
                    rank_list = list(fam_means.index)
                    lor4d_rank = rank_list.index("lorentzian_like_4d") + 1
                    top1 = rank_list[0]
                    gap = fam_means[top1] - fam_means.get("lorentzian_like_4d", 0)
                    summary_rows.append({
                        "n": n, "gamma": gamma, "alpha": alpha,
                        "mode": mode, "lor4d_rank": lor4d_rank,
                        "top1_family": top1,
                        "gap_to_top1": fam_means["lorentzian_like_4d"] - fam_means.iloc[0],
                        "total_families": len(rank_list),
                    })

    sdf = pd.DataFrame(summary_rows)
    sum_path = OUT / "smd_integration_summary.csv"
    sdf.to_csv(sum_path, index=False)
    print(f"Summary: {sum_path} ({len(sdf)} rows)")

    # ─── Quick print: best configs per mode ───
    print("\n" + "=" * 80)
    print("Best Lor4D rank per mode (across all γ, α)")
    print("=" * 80)
    for mode in MODES:
        msub = sdf[sdf["mode"] == mode]
        for n in N_VALUES:
            nsub = msub[msub["n"] == n]
            if nsub.empty:
                continue
            best = nsub.loc[nsub["lor4d_rank"].idxmin()]
            print(f"  {mode:>18}, n={n}: rank={int(best['lor4d_rank'])}, "
                  f"γ={best['gamma']}, α={best['alpha']}, top1={best['top1_family']}")

    # ─── Comparison table at key conditions ───
    print("\n" + "=" * 80)
    print("n=20: info_only vs info+smd — Lor4D rank by (γ, α)")
    print("=" * 80)
    for gamma in GAMMAS:
        row_parts = [f"γ={gamma:>3.1f}"]
        for alpha in ALPHAS:
            info_r = sdf[(sdf["n"] == 20) & (sdf["gamma"] == gamma) &
                         (sdf["alpha"] == alpha) & (sdf["mode"] == "info_only")]
            ismd_r = sdf[(sdf["n"] == 20) & (sdf["gamma"] == gamma) &
                         (sdf["alpha"] == alpha) & (sdf["mode"] == "info_plus_smd")]
            ir = int(info_r["lor4d_rank"].iloc[0]) if len(info_r) else "-"
            isr = int(ismd_r["lor4d_rank"].iloc[0]) if len(ismd_r) else "-"
            row_parts.append(f"α={alpha}: {ir}/{isr}")
        print("  " + " | ".join(row_parts))

    print("\nDone.")


if __name__ == "__main__":
    main()
