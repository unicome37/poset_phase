"""Focused one-sided Ξ_d scan: penalties that only affect d > 4."""
import math
import sys
from collections import defaultdict
sys.path.insert(0, ".")
from xi_d_redesign import load_data, make_model, evaluate, sigmoid

data = load_data()

def make_onesided_model(*, eta_orig: float = 0.6, lam: float = 10.0,
                         alpha0: float = 16.0, q: float = 0.5,
                         Rc: float = 0.25, w: float = 0.015,
                         Rfloor: float = 0.18, alpha0_f: float = 12.0, w_f: float = 0.02,
                         mode: str = "none",
                         kappa: float = 0.0,
                         d_cut: float = 4.0, steepness: float = 1.0,
                         N0: float = 20.0):
    def model(row):
        N = row["N"]
        nlogn = N * math.log(N) if N > 1 else 1.0
        h = row["log_H"] / nlogn
        R = row["R"]
        alpha_N = alpha0 * (N0 / max(N, 1)) ** q
        wall_up = alpha_N * sigmoid((R - Rc) / w)
        alpha_Nf = alpha0_f * (N0 / max(N, 1)) ** q
        wall_lo = alpha_Nf * sigmoid((Rfloor - R) / w_f)
        xi_orig = row["xi_dim"]
        d = row["d_eff"]

        if mode == "none":
            xi_extra = 0.0
        elif mode == "onesided_quad":
            # max(0, d - d_cut)^2 * steepness
            xi_extra = kappa * max(0, d - d_cut) ** 2 * steepness
        elif mode == "onesided_cube":
            # max(0, d - d_cut)^3 * steepness  (steeper)
            xi_extra = kappa * max(0, d - d_cut) ** 3 * steepness
        elif mode == "onesided_exp":
            # kappa * exp(steepness * max(0, d - d_cut)) - kappa
            xi_extra = kappa * (math.exp(steepness * max(0, d - d_cut)) - 1)
        elif mode == "onesided_sig":
            # sigmoid barrier at d_cut
            xi_extra = kappa * sigmoid((d - d_cut) / (0.1 * steepness))
        elif mode == "onesided_pow4":
            # max(0, d - d_cut)^4
            xi_extra = kappa * max(0, d - d_cut) ** 4 * steepness
        else:
            xi_extra = 0.0

        return (h + 0.0004 * row["pi_geo"]
                - lam * row["sigma_hist"]
                + eta_orig * xi_orig
                + xi_extra
                + wall_up + wall_lo)
    return model

configs = [
    ("Baseline", dict(mode="none")),
    # Fine-tune Quad around k=20-60
    ("Quad k=25", dict(mode="onesided_quad", kappa=25)),
    ("Quad k=30", dict(mode="onesided_quad", kappa=30)),
    ("Quad k=35", dict(mode="onesided_quad", kappa=35)),
    ("Quad k=40", dict(mode="onesided_quad", kappa=40)),
    ("Quad k=45", dict(mode="onesided_quad", kappa=45)),
    ("Quad k=50", dict(mode="onesided_quad", kappa=50)),
    # Fine-tune Cube around k=60-150
    ("Cube k=60", dict(mode="onesided_cube", kappa=60)),
    ("Cube k=80", dict(mode="onesided_cube", kappa=80)),
    ("Cube k=100", dict(mode="onesided_cube", kappa=100)),
    ("Cube k=120", dict(mode="onesided_cube", kappa=120)),
    ("Cube k=150", dict(mode="onesided_cube", kappa=150)),
    ("Cube k=200", dict(mode="onesided_cube", kappa=200)),
    # Also vary floor wall with best Quad/Cube
    ("Quad k=35 fl=0.15", dict(mode="onesided_quad", kappa=35, Rfloor=0.15)),
    ("Quad k=35 fl=0.20", dict(mode="onesided_quad", kappa=35, Rfloor=0.20)),
    ("Quad k=35 af=16", dict(mode="onesided_quad", kappa=35, alpha0_f=16)),
    ("Cube k=100 fl=0.15", dict(mode="onesided_cube", kappa=100, Rfloor=0.15)),
    ("Cube k=100 fl=0.20", dict(mode="onesided_cube", kappa=100, Rfloor=0.20)),
    ("Cube k=100 af=16", dict(mode="onesided_cube", kappa=100, alpha0_f=16)),
    ("Cube k=120 af=16", dict(mode="onesided_cube", kappa=120, alpha0_f=16)),
    # Try d_cut=3.9 (tighter than 4.0) with quad
    ("Quad k=30 dc=3.9", dict(mode="onesided_quad", kappa=30, d_cut=3.9)),
    ("Quad k=40 dc=3.9", dict(mode="onesided_quad", kappa=40, d_cut=3.9)),
    ("Quad k=50 dc=3.9", dict(mode="onesided_quad", kappa=50, d_cut=3.9)),
]

print(f"{'Label':<25s} B2   B3   B4   A     4v5a 4v5h 4v2")
print("-" * 78)
for label, kw in configs:
    m = make_onesided_model(**kw)
    r = evaluate(data, m)
    ok = " ***" if r["B4"] >= 99 and r["B2"] >= 99 and r["A5h"] >= 70 else ""
    print(f"{label:<25s} {r['B2']:3.0f}% {r['B3']:3.0f}% {r['B4']:3.0f}% "
          f"{r['A']:5.1f}% {r['A5']:3.0f}% {r['A5h']:3.0f}% {r['A2']:3.0f}%{ok}")
