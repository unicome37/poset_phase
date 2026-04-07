"""Extract key Lor4D rankings from large-N experiment."""
import pandas as pd

bs = pd.read_csv("outputs_info_largeN/bootstrap_summary.csv")
sm = pd.read_csv("outputs_info_largeN/summary.csv")

print("=" * 90)
print("大 N 扩展实验结果 (A4 模态)")
print("=" * 90)

for n_val in [60, 80, 100]:
    print(f"\n{'='*90}\nn = {n_val}")
    for g in [0.0, 0.2, 0.5, 1.0]:
        sub = sm[(sm["n"] == n_val) & (sm["gamma"] == g) & (sm["action_mode"] == "A4")].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("mean_score_norm")
        ranking = list(sub["family"].values)
        lor4d_row = sub[sub["family"] == "lorentzian_like_4d"]
        kr_row = sub[sub["family"] == "KR_like"]
        if lor4d_row.empty:
            continue
        lor_rank = ranking.index("lorentzian_like_4d") + 1
        kr_rank = ranking.index("KR_like") + 1 if "KR_like" in ranking else "N/A"
        lor_score = lor4d_row["mean_score_norm"].values[0]
        kr_score = kr_row["mean_score_norm"].values[0] if not kr_row.empty else float("nan")
        gap = lor_score - kr_score if not kr_row.empty else float("nan")

        # Bootstrap CI for Lor4D
        b_lor = bs[(bs["n"] == n_val) & (bs["gamma"] == g) & (bs["action_mode"] == "A4") &
                    (bs["family"] == "lorentzian_like_4d")]
        ci_lo = b_lor["score_norm_ci_low"].values[0] if not b_lor.empty else float("nan")
        ci_hi = b_lor["score_norm_ci_high"].values[0] if not b_lor.empty else float("nan")
        
        top3 = ", ".join([f"{ranking[i]}({sub['mean_score_norm'].values[i]:.3f})" for i in range(min(3, len(ranking)))])
        
        print(f"  γ={g:.1f}: Lor4D rank #{lor_rank}, KR_like rank #{kr_rank}, "
              f"gap={gap:+.4f}, Lor4D CI=[{ci_lo:.3f}, {ci_hi:.3f}]")
        print(f"         Top 3: {top3}")

# Cross-N trend for Lor4D at gamma=1.0
print(f"\n{'='*90}")
print("Lor4D rank trend at γ=1.0:")
for n_val in [60, 80, 100]:
    sub = sm[(sm["n"] == n_val) & (sm["gamma"] == 1.0) & (sm["action_mode"] == "A4")].copy()
    sub = sub.sort_values("mean_score_norm")
    if sub.empty:
        continue
    ranking = list(sub["family"].values)
    lor_rank = ranking.index("lorentzian_like_4d") + 1 if "lorentzian_like_4d" in ranking else "N/A"
    lor_score = sub[sub["family"]=="lorentzian_like_4d"]["mean_score_norm"].values[0]
    print(f"  n={n_val:3d}: rank #{lor_rank}, score_norm={lor_score:+.4f}")
