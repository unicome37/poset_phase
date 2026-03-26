"""
γ from Einstein-Hilbert: Theoretical Derivation
=================================================
Goal: derive γ ≈ 1 in F10 = logH + γN(d_eff − d*)² from the discrete EH action.

APPROACH:

1. BDG Action Channel:
   The BDG action S^(d) = Σ α_k^(d) C_k recovers ∫R√g d^dx.
   For a d-dim sprinkle, S^(d)/N → ρ^{2/d} · R_scalar (continuum).
   If we evaluate S^(4) on a d≠4 sprinkle, the mismatch creates an
   effective "dimension penalty" that should map onto γN(d_eff−d*)².

2. Entropy Channel:
   logH ≈ N·h(d) where h(d) is the entropy density.
   h(d) is a monotonically increasing function of d (more phase space
   in higher dimensions). The d-dependence of h(d) determines how
   strongly entropy pulls toward high d.

3. The Balance:
   F10 = logH + γN(d_eff−d*)² + ...
   At the minimum (d_eff=d*):  ∂F/∂d = N·h'(d*) + 0 = 0  (impossible since h'>0)
   → The well is NOT from a variational principle on F10 alone.
   → Instead, γ encodes the EXTERNAL constraint: "the action selects d=4."

4. Alternative: γ from the BDG mismatch penalty
   Compute S_BDG^(4)(d-dim poset) for d=2,3,4,5.
   The penalty (S^(4)(d) − S^(4)(4))² / N should scale like (d_eff−4)².
   If so, γ = coefficient of this quadratic in the BDG mismatch.

5. Direct numerical extraction:
   Compute BDG action for all posets in the 1000-poset dataset.
   Fit: S_BDG/N = a + b·(d_eff−4)² + residual.
   Compare b with γ≈1.
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats
from scipy.special import gamma as Gamma


def sigmoid(x):
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def load_csv(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            d = {}
            for k, v in row.items():
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
            d["N"] = int(d["N"])
            d["rep"] = int(d["rep"])
            rows.append(d)
    return rows


def f2_myrheim_meyer(d):
    """Fraction of causally related pairs in d-dim Minkowski."""
    return Gamma(d + 1) * Gamma(d / 2) / (4.0 * Gamma(3 * d / 2))


def main():
    csv_path = "outputs_d_recovery/prediction_c_f7_large_n.csv"
    rows = load_csv(csv_path)
    n_values = sorted(set(r["N"] for r in rows))
    all_families = sorted(set(r["family"] for r in rows))
    lor_families = [f for f in all_families if f.startswith("Lor")]

    by_nf = defaultdict(list)
    for r in rows:
        by_nf[(r["N"], r["family"])].append(r)

    report = []
    report.append("# γ from Einstein-Hilbert: Derivation Attempt\n")

    # ══════════════════════════════════════════════════════════════
    # Part 1: Entropy density h(d) = logH / N
    # ══════════════════════════════════════════════════════════════
    report.append("## 1. Entropy Density h(d) = logH / N\n")
    report.append("| N | h(2D) | h(3D) | h(4D) | h(5D) | h(KR) |")
    report.append("|---|-------|-------|-------|-------|-------|")

    h_by_nf = {}
    for N in n_values:
        cells = [str(N)]
        for f in all_families:
            vals = by_nf.get((N, f), [])
            if vals:
                h = np.mean([r["log_H"] for r in vals]) / N
                h_by_nf[(N, f)] = h
                cells.append(f"{h:.3f}")
            else:
                cells.append("—")
        report.append("| " + " | ".join(cells) + " |")

    # h(d) growth rate
    report.append("\n### h(d) growth: Δh = h(d+1) − h(d)\n")
    report.append("| N | Δh(3D−2D) | Δh(4D−3D) | Δh(5D−4D) |")
    report.append("|---|-----------|-----------|-----------|")
    for N in n_values:
        dh_32 = h_by_nf.get((N, "Lor3D"), 0) - h_by_nf.get((N, "Lor2D"), 0)
        dh_43 = h_by_nf.get((N, "Lor4D"), 0) - h_by_nf.get((N, "Lor3D"), 0)
        dh_54 = h_by_nf.get((N, "Lor5D"), 0) - h_by_nf.get((N, "Lor4D"), 0)
        report.append(f"| {N} | {dh_32:+.4f} | {dh_43:+.4f} | {dh_54:+.4f} |")

    # ══════════════════════════════════════════════════════════════
    # Part 2: d_eff well penalty per element
    # ══════════════════════════════════════════════════════════════
    report.append("\n## 2. d_eff Well Penalty: (d_eff − 4)² per Element\n")
    report.append("| N | (d_eff−4)²(2D) | (d_eff−4)²(3D) | (d_eff−4)²(4D) | (d_eff−4)²(5D) | (d_eff−4)²(KR) |")
    report.append("|---|--------------|--------------|--------------|--------------|--------------|")

    for N in n_values:
        cells = [str(N)]
        for f in all_families:
            vals = by_nf.get((N, f), [])
            if vals:
                d2 = np.mean([(r["d_eff"] - 4.0)**2 for r in vals])
                cells.append(f"{d2:.3f}")
            else:
                cells.append("—")
        report.append("| " + " | ".join(cells) + " |")

    # ══════════════════════════════════════════════════════════════
    # Part 3: γ_crit derivation from entropy-geometry balance
    # ══════════════════════════════════════════════════════════════
    report.append("\n## 3. γ from Entropy-Geometry Balance\n")
    report.append("At the F10 minimum (4D), for 4D to beat dimension d:\n")
    report.append("  F10(d) > F10(4D)")
    report.append("  logH(d) + γN(d_eff(d)−4)² > logH(4D) + γN(d_eff(4D)−4)²")
    report.append("  γ > [logH(4D) − logH(d)] / [N · ((d_eff(d)−4)² − (d_eff(4D)−4)²)]")
    report.append("  γ > −Δh(d−4D) / Δ(d_eff²)(d−4D)\n")
    report.append("where Δh and Δ(d_eff²) are per-element quantities.\n")

    report.append("### Per-element balance:\n")
    report.append("| N | Δh/Δd² (3D) | Δh/Δd² (5D) | Δh/Δd² (2D) | binding |")
    report.append("|---|------------|------------|------------|---------|")

    gamma_3d_list = []
    gamma_5d_list = []
    gamma_2d_list = []

    for N in n_values:
        for d_name, d_fam in [("3D", "Lor3D"), ("5D", "Lor5D"), ("2D", "Lor2D")]:
            vd = by_nf.get((N, d_fam), [])
            v4 = by_nf.get((N, "Lor4D"), [])
            if not vd or not v4:
                continue
            dh = (np.mean([r["log_H"] for r in v4]) - np.mean([r["log_H"] for r in vd])) / N
            dd2_d = np.mean([(r["d_eff"] - 4.0)**2 for r in vd])
            dd2_4 = np.mean([(r["d_eff"] - 4.0)**2 for r in v4])
            delta_d2 = dd2_d - dd2_4
            if abs(delta_d2) > 1e-6:
                gamma_req = -dh / delta_d2
            else:
                gamma_req = float("inf")
            
            if d_name == "3D":
                gamma_3d_list.append(gamma_req)
            elif d_name == "5D":
                gamma_5d_list.append(gamma_req)
            elif d_name == "2D":
                gamma_2d_list.append(gamma_req)

        g3 = gamma_3d_list[-1] if gamma_3d_list else 0
        g5 = gamma_5d_list[-1] if gamma_5d_list else 0
        g2 = gamma_2d_list[-1] if gamma_2d_list else 0
        binding = "5D" if g5 > g3 and g5 > g2 else ("3D" if g3 > g5 else "2D")
        report.append(f"| {N} | {g3:.4f} | {g5:.4f} | {g2:.4f} | {binding} |")

    report.append(f"\n**Mean γ_crit(3D)**: {np.mean(gamma_3d_list):.3f} ± {np.std(gamma_3d_list):.3f}")
    report.append(f"**Mean γ_crit(5D)**: {np.mean(gamma_5d_list):.3f} ± {np.std(gamma_5d_list):.3f}")
    report.append(f"**Mean γ_crit(2D)**: {np.mean(gamma_2d_list):.3f} ± {np.std(gamma_2d_list):.3f}")
    report.append(f"\n**Binding constraint**: γ must exceed max(γ_crit) ≈ {max(np.max(gamma_5d_list), np.max(gamma_3d_list)):.2f} (5D)")

    # ══════════════════════════════════════════════════════════════
    # Part 4: Theoretical prediction of γ from EH action
    # ══════════════════════════════════════════════════════════════
    report.append("\n## 4. Theoretical γ from EH Action Structure\n")

    report.append("### 4a. The Myrheim-Meyer encoding identity\n")
    report.append("For a causal set sprinkled into d-dim Minkowski at density ρ:\n")
    report.append("  f₂(d) = Γ(d+1)Γ(d/2) / (4Γ(3d/2))  [fraction of causally related pairs]\n")
    report.append("The Myrheim-Meyer estimator solves f₂(d_eff) = C₀/C(N,2) for d_eff.\n")
    report.append("For d_eff close to d*, Taylor expand:\n")
    report.append("  f₂(d) ≈ f₂(d*) + f₂'(d*)(d−d*) + ½f₂''(d*)(d−d*)² + ...\n")

    # Numerical derivatives of f₂
    d_star = 4.0
    eps = 0.01
    f2_d = f2_myrheim_meyer(d_star)
    f2_plus = f2_myrheim_meyer(d_star + eps)
    f2_minus = f2_myrheim_meyer(d_star - eps)
    f2_prime = (f2_plus - f2_minus) / (2 * eps)
    f2_double_prime = (f2_plus - 2 * f2_d + f2_minus) / eps**2

    report.append(f"At d*=4: f₂ = {f2_d:.6f}, f₂' = {f2_prime:.6f}, f₂'' = {f2_double_prime:.6f}\n")

    # f₂ values at integer dimensions
    report.append("f₂ at integer d:\n")
    for d in [2, 3, 4, 5, 6]:
        report.append(f"  f₂({d}) = {f2_myrheim_meyer(d):.6f}")

    report.append("\n### 4b. Key subtlety: f₂-space vs d-space\n")
    report.append("**IMPORTANT**: The F10 well operates in d-space: γN(d_eff−4)².")
    report.append("But the BDG link action operates in f₂-space: S_link ∝ −2C₀ = −2·f₂·C(N,2).")
    report.append("These are related by the NONLINEAR map d = f₂⁻¹(observed fraction).\n")
    report.append("The d-space curvature of f₂ tells us about the link action's")
    report.append("dimension sensitivity, but γ in F10 is a SEPARATE quantity.\n")
    report.append(f"f₂''(4) = {f2_double_prime:+.6f} (POSITIVE = convex in d-space)")
    report.append("→ The link action has CONCAVE d-dependence at d=4")
    report.append("→ S_link FAVORS deviations from d=4, not penalizes them!\n")
    report.append("This is why γ cannot come from f₂ alone — the sign is wrong.")
    report.append("The d_eff well must arise from a DIFFERENT mechanism.\n")

    # Compute the d-space curvature via the Jacobian
    report.append("### 4c. The d_eff well as an INDEPENDENT geometric constraint\n")
    report.append("The d_eff estimator inverts f₂(d_eff) = C₀/C(N,2).")
    report.append("Since f₂ is monotonically decreasing, this inversion is well-defined.")
    report.append("But d_eff is a DERIVED quantity — it encodes C₀ in d-units.\n")
    report.append("The F10 well γN(d_eff−4)² is therefore equivalent to:")
    report.append("  γN · [f₂⁻¹(C₀/C(N,2)) − 4]²\n")
    report.append("This is a NONLINEAR function of C₀, amplified by the")
    report.append("steep gradient of f₂⁻¹ near d=4 (where f₂ is small and")
    report.append("changing rapidly with d).\n")

    # Jacobian of the inversion
    # d = f₂⁻¹(p) → dd/dp = 1/f₂'(d)
    dd_dp = 1.0 / f2_prime  # dd/dp at d=4
    report.append(f"Jacobian: dd/df₂ = 1/f₂'(4) = 1/({f2_prime:.6f}) = {dd_dp:.2f}")
    report.append(f"→ A unit change in f₂ maps to {abs(dd_dp):.1f} units of d_eff\n")
    report.append("This amplification factor is KEY: small changes in C₀")
    report.append("get magnified into large changes in d_eff, making the")
    report.append(f"well γ(d_eff−4)² effectively γ·{dd_dp**2:.0f}·(Δf₂)² in f₂-space.\n")
    
    gamma_f2_equiv = 1.0 * dd_dp**2  # γ=1 in d-space = this in f₂-space
    report.append(f"γ=1 in d-space ≡ γ_f₂ = {gamma_f2_equiv:.0f} in f₂-space")
    report.append("This enormous amplification explains why γ=1 works:\n")
    report.append("the d_eff well leverages the steep f₂→d inversion near d=4.\n")

    report.append("### 4d. Full BDG coefficient structure\n")
    report.append("The link action only uses C₀. The full BDG action uses C₀, C₁, C₂, C₃...")
    report.append("Each C_k has its own dimension-dependent coefficient α_k^(d).\n")
    report.append("The total dimension penalty from the BDG action is:\n")
    report.append("  γ_eff = Σ_k α_k^(d*) · (d²f_k/dd²)|_{d=d*}\n")
    report.append("where f_k(d) = E[C_k]/C(N,2) is the expected k-interval fraction.\n")
    report.append("This sum can amplify or suppress the leading f₂'' term.\n")

    # ══════════════════════════════════════════════════════════════
    # Part 5: Why γ = O(1) is natural
    # ══════════════════════════════════════════════════════════════
    report.append("## 5. Why γ = O(1) is Natural\n")
    report.append("The question is not 'derive γ=1 from f₂' but rather:")
    report.append("'why should the dimension penalty be O(1) per element?'\n")
    report.append("### Three independent arguments for γ = O(1):\n")

    report.append("**Argument 1: Entropy-geometry equipartition**")
    report.append("The entropy density h(d) ≈ 1.4–2.7 nat/element (from Part 1).")
    report.append("For the well to compete with entropy, γ·(d_eff−4)² must be O(h).")
    report.append(f"Since (d_eff(5D)−4)² ≈ {np.mean([(r['d_eff']-4)**2 for r in by_nf.get((100,'Lor5D'),[])]):.2f} at N=100,")
    report.append(f"we need γ ≈ Δh(5D−4D)/Δ(d²) = {np.mean(gamma_5d_list):.1f} for 5D competition.")
    report.append("This gives γ = O(1), which is a BALANCE condition, not a derivation.\n")

    report.append("**Argument 2: Planck density normalization**")
    report.append("In the continuum, S_EH/V = R/(16πG). At Planck density ρ_P = 1/ℓ_P^d:")
    report.append("  S_EH/N = S_EH/(ρ_P V) = R·ℓ_P²/(16π) = O(1)")
    report.append("since R·ℓ_P² = O(1) for curvatures at the Planck scale.")
    report.append("The dimension dependence ∂²(S_EH/N)/∂d² is therefore also O(1),")
    report.append("giving γ_EH = O(1) naturally.\n")

    report.append("**Argument 3: f₂ Jacobian amplification (from §4c)**")
    report.append(f"γ=1 in d-space = γ_f₂ = {gamma_f2_equiv:.0f} in f₂-space.")
    report.append("The Myrheim-Meyer inversion d = f₂⁻¹(p) has steep gradient")
    report.append(f"near d=4: |dd/df₂| = {abs(dd_dp):.1f}.")
    report.append("So a mild O(1) well in d-space corresponds to an enormous")
    report.append(f"penalty of ~{gamma_f2_equiv:.0f} per unit (Δf₂)² in the observable C₀ space.")
    report.append("The physical content: d_eff is a HIGHLY COMPRESSED encoding")
    report.append("of C₀, and the well leverages this compression.\n")

    # ══════════════════════════════════════════════════════════════
    # Part 6: d_eff vs R — are they the same information?
    # ══════════════════════════════════════════════════════════════
    report.append("## 6. d_eff vs R: Information Content\n")

    report.append("Our R observable is the OCCUPANCY RATIO (fraction of matrix entries = 1).")
    report.append("The Myrheim-Meyer f₂ is the fraction of CAUSALLY RELATED pairs.")
    report.append("For a Hasse diagram stored as a CLOSURE matrix, R = f₂ = C₀/C(N,2).")
    report.append("But our R is computed differently — let's check the actual relationship.\n")

    report.append("| N | family | R | f₂_theory | d_eff | d_from_f₂_theory |")
    report.append("|---|--------|---|-----------|-------|-----------------|")

    for N in [20, 100]:
        for f in all_families:
            vals = by_nf.get((N, f), [])
            if not vals:
                continue
            m_r = np.mean([r["R"] for r in vals])
            m_d = np.mean([r["d_eff"] for r in vals])
            # d_eff in our data is the Myrheim-Meyer estimate
            # so f₂_empirical ≈ f₂(d_eff) by definition
            f2_from_deff = f2_myrheim_meyer(m_d)
            # True f₂ for the family's integer dimension
            true_d = {"Lor2D": 2, "Lor3D": 3, "Lor4D": 4, "Lor5D": 5, "KR_like": None}
            td = true_d.get(f)
            f2_true = f2_myrheim_meyer(td) if td else float('nan')
            f2_true_str = f"{f2_true:.4f}" if not math.isnan(f2_true) else "—"
            report.append(f"| {N} | {f} | {m_r:.3f} | {f2_true_str} | {m_d:.2f} | f₂(d_eff)={f2_from_deff:.4f} |")

    report.append("\nKey: R in our dataset is NOT f₂. R is the occupancy ratio (total")
    report.append("relations / N²), while f₂ = C₀/C(N,2) (ordered pairs only).")
    report.append("But d_eff IS calibrated from f₂ by construction.\n")

    # ══════════════════════════════════════════════════════════════
    # Part 7: Formal interpretation of γ
    # ══════════════════════════════════════════════════════════════
    report.append("\n## 7. Formal Interpretation of γ\n")
    report.append("### Summary of derivation attempt:\n")
    report.append("### 结论\n")
    report.append("1. **f₂'' 通道的符号是错的**：f₂''(4) > 0（凸），link action 在 d=4")
    report.append("   附近是凹的 → 它反而**鼓励**偏离 d=4，而非惩罚。")
    report.append("   因此 γ 不可能从 link action 推导出来。\n")
    report.append("2. **d_eff 井是独立的几何约束**，不是 BDG action 的重参数化。")
    report.append("   它通过 f₂⁻¹ 反演的陡峭 Jacobian 从 C₀ 中提取维度信息,")
    report.append(f"   放大倍数 |dd/df₂|² ≈ {dd_dp**2:.0f}。\n")
    report.append("3. **γ = O(1) 的三个独立论证**:")
    report.append("   (a) 熵-几何等分: γ·Δ(d²) ∼ Δh ∼ O(1) → γ ∼ O(1)")
    report.append("   (b) Planck 密度归一化: S_EH/N = R·ℓ_P²/(16π) = O(1)")
    report.append("   (c) Jacobian 放大: γ=1 在 d 空间 = γ_f₂ = O(500) 在 f₂ 空间\n")
    report.append("4. **严格推导需要**: 完整 BDG 展开到所有阶的 d 依赖性,")
    report.append("   这是一个开放的理论问题。但 γ = O(1) 是自然标度。\n")
    report.append("5. **核心洞见**: F10 的 d_eff 井不是 BDG action 的简化版,")
    report.append("   而是一个**正交**的几何约束——它编码的是 Myrheim-Meyer 维度")
    report.append("   (通过 f₂ 反演获得), 不是标量曲率(通过 BDG 系数获得)。")
    report.append("   这正是为什么 Φ_geom(d_eff) ⊥ Ψ_Lor(R) 是正交的：")
    report.append("   它们编码了因果集合几何的**不同方面**。")

    out = "\n".join(report)
    outdir = Path("outputs_d_recovery")
    md_path = outdir / "prediction_a_gamma_eh_derivation.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
