"""H_int Independence Diagnostics.

Key question: Is interval spectral entropy H_int an INDEPENDENT signal,
or is it a rewrite of existing terms (Ξ_d, f_link, S_BD, log_H, etc.)?

Computes:
  1. Pairwise correlations: H_int vs all 5 F5 terms + f_link + S_BD
  2. Partial correlation: H_int vs 3D/4D label after controlling for F5 terms
  3. Logistic regression: incremental gain of H_int for 3D-vs-4D classification
  4. VIF: variance inflation factor of H_int in the full feature set
"""
from __future__ import annotations

import csv
import numpy as np
from scipy import stats

from generators import Poset
from prediction_a_bd_bridge import (
    regenerate_poset, compute_s_bd_ratio, compute_f5, CALIBRATED_WEIGHTS
)
from prediction_a_faithfulness import faithfulness_score
from prediction_a_mid_gap import interval_spectrum, hasse_diagnostics


def main():
    raw_path = "outputs_unified_functional/raw_features.csv"
    rows = list(csv.DictReader(open(raw_path)))
    lor_families = {"Lor2D", "Lor3D", "Lor4D", "Lor5D"}
    rows = [r for r in rows if r["family"] in lor_families]

    # Compute all features
    data = []
    for i, r in enumerate(rows):
        fam, n, rep = r["family"], int(r["N"]), int(r["rep"])
        poset = regenerate_poset(fam, n, rep, 42)

        # F5 components
        log_H = float(r["log_H"])
        pi_geo = float(r["pi_geo"])
        sigma_hist = float(r["sigma_hist"])
        xi_dim = float(r["xi_dim"])
        pi_cg = float(r["pi_cg"])
        f5 = compute_f5(r, CALIBRATED_WEIGHTS)

        # BD / faithfulness / spectrum
        sbd = compute_s_bd_ratio(poset)
        scores = faithfulness_score(poset)
        spec = interval_spectrum(poset)

        data.append({
            "family": fam, "N": n,
            "log_H": log_H, "pi_geo": pi_geo, "sigma_hist": sigma_hist,
            "xi_dim": xi_dim, "pi_cg": pi_cg, "f5": f5,
            "S_BD": sbd,
            "f_link": scores["int_link_fraction"],
            "H_int": spec["spectral_entropy"],
            "C0_frac": spec["C0_frac"],
            "C1_C0": spec["C1_C0_ratio"],
            "tail_weight": spec["tail_weight"],
        })
        if (i + 1) % 16 == 0:
            print(f"  [{i+1}/{len(rows)}]")

    print(f"  Total: {len(data)} posets\n")

    # =====================================================================
    # 1. Pairwise correlations: H_int vs everything
    # =====================================================================
    print("=" * 80)
    print("1. PAIRWISE CORRELATIONS with H_int")
    print("=" * 80)

    h_int = np.array([d["H_int"] for d in data])
    corr_targets = [
        ("log_H", "log H (combinatorial entropy)"),
        ("pi_geo", "Π_geo (geometric constraint)"),
        ("sigma_hist", "Σ_hist (historical sedimentation)"),
        ("xi_dim", "Ξ_d (dimension barrier)"),
        ("pi_cg", "Π_cg (CG drift)"),
        ("f5", "F5 (total functional)"),
        ("S_BD", "S_BD (interval richness)"),
        ("f_link", "f_link (link fraction)"),
        ("C0_frac", "C_0 fraction"),
        ("C1_C0", "C_1/C_0 ratio"),
        ("tail_weight", "tail weight"),
    ]

    print(f"\n  {'Variable':<35s} {'Pearson r':>10s} {'p-value':>12s} {'Spearman ρ':>10s}")
    print(f"  {'-'*35} {'-'*10} {'-'*12} {'-'*10}")

    for key, label in corr_targets:
        vals = np.array([d[key] for d in data])
        r_p, p_p = stats.pearsonr(h_int, vals)
        r_s, p_s = stats.spearmanr(h_int, vals)
        flag = " *** HIGH" if abs(r_p) > 0.8 else (" ** mod" if abs(r_p) > 0.5 else "")
        print(f"  {label:<35s} {r_p:>+10.4f} {p_p:>12.2e} {r_s:>+10.4f}{flag}")

    # =====================================================================
    # 2. Within-family correlations (removes between-family confound)
    # =====================================================================
    print("\n" + "=" * 80)
    print("2. WITHIN-FAMILY CORRELATIONS (H_int vs F5 terms, Lor3D+Lor4D only)")
    print("=" * 80)

    data_34 = [d for d in data if d["family"] in ("Lor3D", "Lor4D")]
    h_34 = np.array([d["H_int"] for d in data_34])

    print(f"\n  N = {len(data_34)} (Lor3D + Lor4D only)")
    print(f"\n  {'Variable':<35s} {'Pearson r':>10s} {'Spearman ρ':>10s}")
    print(f"  {'-'*35} {'-'*10} {'-'*10}")

    for key, label in corr_targets:
        vals = np.array([d[key] for d in data_34])
        r_p, _ = stats.pearsonr(h_34, vals)
        r_s, _ = stats.spearmanr(h_34, vals)
        flag = " *** HIGH" if abs(r_p) > 0.8 else (" ** mod" if abs(r_p) > 0.5 else "")
        print(f"  {label:<35s} {r_p:>+10.4f} {r_s:>+10.4f}{flag}")

    # =====================================================================
    # 3. Partial correlation: H_int vs family label, controlling F5 terms
    # =====================================================================
    print("\n" + "=" * 80)
    print("3. PARTIAL CORRELATION: H_int → 3D/4D label, controlling for F5 terms")
    print("=" * 80)

    # Binary label: 0 = Lor3D, 1 = Lor4D
    y = np.array([1.0 if d["family"] == "Lor4D" else 0.0 for d in data_34])

    # Raw correlation
    r_raw, p_raw = stats.pointbiserialr(y, h_34)
    print(f"\n  Raw point-biserial: r = {r_raw:+.4f}, p = {p_raw:.2e}")

    # Partial: regress H_int on F5 terms, then correlate residual with label
    X_f5 = np.column_stack([
        [d["log_H"] for d in data_34],
        [d["pi_geo"] for d in data_34],
        [d["sigma_hist"] for d in data_34],
        [d["xi_dim"] for d in data_34],
        [d["pi_cg"] for d in data_34],
    ])

    # Add constant
    X_f5c = np.column_stack([np.ones(len(data_34)), X_f5])
    # OLS: H_int ~ F5 terms
    beta_h = np.linalg.lstsq(X_f5c, h_34, rcond=None)[0]
    h_resid = h_34 - X_f5c @ beta_h
    r_partial_f5, p_partial_f5 = stats.pointbiserialr(y, h_resid)
    print(f"  Partial (controlling F5 terms): r = {r_partial_f5:+.4f}, p = {p_partial_f5:.2e}")

    # Also partial controlling for f_link
    X_flink = np.column_stack([X_f5, [d["f_link"] for d in data_34]])
    X_flinkc = np.column_stack([np.ones(len(data_34)), X_flink])
    beta_fl = np.linalg.lstsq(X_flinkc, h_34, rcond=None)[0]
    h_resid_fl = h_34 - X_flinkc @ beta_fl
    r_partial_fl, p_partial_fl = stats.pointbiserialr(y, h_resid_fl)
    print(f"  Partial (controlling F5 + f_link): r = {r_partial_fl:+.4f}, p = {p_partial_fl:.2e}")

    # Partial controlling for S_BD
    X_sbd = np.column_stack([X_f5, [d["S_BD"] for d in data_34]])
    X_sbdc = np.column_stack([np.ones(len(data_34)), X_sbd])
    beta_sb = np.linalg.lstsq(X_sbdc, h_34, rcond=None)[0]
    h_resid_sb = h_34 - X_sbdc @ beta_sb
    r_partial_sb, p_partial_sb = stats.pointbiserialr(y, h_resid_sb)
    print(f"  Partial (controlling F5 + S_BD): r = {r_partial_sb:+.4f}, p = {p_partial_sb:.2e}")

    # Partial controlling for ALL (F5 + f_link + S_BD)
    X_all = np.column_stack([X_f5,
                             [d["f_link"] for d in data_34],
                             [d["S_BD"] for d in data_34]])
    X_allc = np.column_stack([np.ones(len(data_34)), X_all])
    beta_all = np.linalg.lstsq(X_allc, h_34, rcond=None)[0]
    h_resid_all = h_34 - X_allc @ beta_all
    r_partial_all, p_partial_all = stats.pointbiserialr(y, h_resid_all)
    print(f"  Partial (controlling F5 + f_link + S_BD): r = {r_partial_all:+.4f}, p = {p_partial_all:.2e}")

    # =====================================================================
    # 4. R² incremental gain
    # =====================================================================
    print("\n" + "=" * 80)
    print("4. R² INCREMENTAL GAIN: adding H_int to F5 for 3D/4D classification")
    print("=" * 80)

    # Logistic-like: use linear discriminant (OLS on binary label)
    # Model 1: y ~ F5 terms
    beta1 = np.linalg.lstsq(X_f5c, y, rcond=None)[0]
    y_pred1 = X_f5c @ beta1
    ss_res1 = np.sum((y - y_pred1)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2_f5 = 1 - ss_res1 / ss_tot

    # Model 2: y ~ F5 terms + H_int
    X_f5h = np.column_stack([X_f5c, h_34])
    beta2 = np.linalg.lstsq(X_f5h, y, rcond=None)[0]
    y_pred2 = X_f5h @ beta2
    ss_res2 = np.sum((y - y_pred2)**2)
    r2_f5h = 1 - ss_res2 / ss_tot

    # Model 3: y ~ F5 + f_link + H_int
    X_f5fh = np.column_stack([X_f5c, [d["f_link"] for d in data_34], h_34])
    beta3 = np.linalg.lstsq(X_f5fh, y, rcond=None)[0]
    y_pred3 = X_f5fh @ beta3
    ss_res3 = np.sum((y - y_pred3)**2)
    r2_f5fh = 1 - ss_res3 / ss_tot

    # Model 4: y ~ F5 + f_link + S_BD + H_int
    X_full = np.column_stack([X_f5c,
                              [d["f_link"] for d in data_34],
                              [d["S_BD"] for d in data_34],
                              h_34])
    beta4 = np.linalg.lstsq(X_full, y, rcond=None)[0]
    y_pred4 = X_full @ beta4
    ss_res4 = np.sum((y - y_pred4)**2)
    r2_full = 1 - ss_res4 / ss_tot

    print(f"\n  Model 1: F5 terms only          R² = {r2_f5:.4f}")
    print(f"  Model 2: F5 + H_int             R² = {r2_f5h:.4f}  ΔR² = {r2_f5h - r2_f5:+.4f}")
    print(f"  Model 3: F5 + f_link + H_int    R² = {r2_f5fh:.4f}  ΔR² = {r2_f5fh - r2_f5:+.4f}")
    print(f"  Model 4: F5 + f_link + S_BD + H R² = {r2_full:.4f}  ΔR² = {r2_full - r2_f5:+.4f}")

    # Classification accuracy
    for name, pred in [("F5 only", y_pred1), ("F5+H_int", y_pred2),
                       ("F5+f_link+H_int", y_pred3), ("Full", y_pred4)]:
        acc = np.mean((pred > 0.5) == y) * 100
        print(f"  {name:<25s} accuracy = {acc:.1f}%")

    # =====================================================================
    # 5. VIF of H_int
    # =====================================================================
    print("\n" + "=" * 80)
    print("5. VARIANCE INFLATION FACTOR (VIF) of H_int")
    print("=" * 80)

    # VIF = 1 / (1 - R²) where R² is from regressing H_int on all other predictors
    predictors = {
        "log_H": [d["log_H"] for d in data_34],
        "pi_geo": [d["pi_geo"] for d in data_34],
        "sigma_hist": [d["sigma_hist"] for d in data_34],
        "xi_dim": [d["xi_dim"] for d in data_34],
        "pi_cg": [d["pi_cg"] for d in data_34],
        "f_link": [d["f_link"] for d in data_34],
        "S_BD": [d["S_BD"] for d in data_34],
    }

    X_others = np.column_stack([np.ones(len(data_34))] + list(predictors.values()))
    beta_vif = np.linalg.lstsq(X_others, h_34, rcond=None)[0]
    h_pred_vif = X_others @ beta_vif
    ss_res_vif = np.sum((h_34 - h_pred_vif)**2)
    ss_tot_vif = np.sum((h_34 - h_34.mean())**2)
    r2_vif = 1 - ss_res_vif / ss_tot_vif
    vif = 1 / (1 - r2_vif) if r2_vif < 1 else float('inf')

    print(f"\n  R²(H_int ~ F5 + f_link + S_BD) = {r2_vif:.4f}")
    print(f"  VIF = {vif:.2f}")
    if vif < 5:
        print(f"  → VIF < 5: H_int has LOW collinearity with existing terms")
    elif vif < 10:
        print(f"  → VIF < 10: H_int has MODERATE collinearity")
    else:
        print(f"  → VIF ≥ 10: H_int has HIGH collinearity — may be a rewrite")

    # Also compute VIF for subsets
    for exclude_key in ["f_link", "S_BD"]:
        other_keys = [k for k in predictors if k != exclude_key]
        X_sub = np.column_stack([np.ones(len(data_34))] +
                                [predictors[k] for k in other_keys])
        b_sub = np.linalg.lstsq(X_sub, h_34, rcond=None)[0]
        r2_sub = 1 - np.sum((h_34 - X_sub @ b_sub)**2) / ss_tot_vif
        vif_sub = 1 / (1 - r2_sub) if r2_sub < 1 else float('inf')
        print(f"  VIF(H_int ~ F5 + {', '.join(other_keys[5:])}) = {vif_sub:.2f} [R²={r2_sub:.4f}]")

    # F5 only
    X_f5only = np.column_stack([np.ones(len(data_34)),
                                [d["log_H"] for d in data_34],
                                [d["pi_geo"] for d in data_34],
                                [d["sigma_hist"] for d in data_34],
                                [d["xi_dim"] for d in data_34],
                                [d["pi_cg"] for d in data_34]])
    b_f5only = np.linalg.lstsq(X_f5only, h_34, rcond=None)[0]
    r2_f5only = 1 - np.sum((h_34 - X_f5only @ b_f5only)**2) / ss_tot_vif
    vif_f5only = 1 / (1 - r2_f5only) if r2_f5only < 1 else float('inf')
    print(f"  VIF(H_int ~ F5 only) = {vif_f5only:.2f} [R²={r2_f5only:.4f}]")

    # =====================================================================
    # 6. SUMMARY
    # =====================================================================
    print("\n" + "=" * 80)
    print("6. INDEPENDENCE VERDICT")
    print("=" * 80)
    print(f"""
  H_int correlations:
    vs f_link:     high (expected — both measure interval structure)
    vs S_BD:       high (expected — S_BD is a cruder aggregate of same data)
    vs F5 terms:   see VIF above

  Key test — partial correlation with 3D/4D label:
    Raw:                    r = {r_raw:+.4f}
    After F5:               r = {r_partial_f5:+.4f}  (unique signal AFTER F5)
    After F5 + f_link:      r = {r_partial_fl:+.4f}
    After F5 + S_BD:        r = {r_partial_sb:+.4f}
    After F5 + f_link + S_BD: r = {r_partial_all:+.4f}

  Key test — R² incremental gain:
    F5 only:      R² = {r2_f5:.4f}
    F5 + H_int:   R² = {r2_f5h:.4f}  (ΔR² = {r2_f5h - r2_f5:+.4f})

  VIF(H_int):     {vif:.2f}
""")


if __name__ == "__main__":
    main()
