"""Diagnose the normalization reversal at n=60."""
import pandas as pd

df = pd.read_csv("outputs_info_full/raw_samples.csv")
sub = df[(df["n"] == 60) & (df["gamma"] == 1.0) & (df["action_mode"] == "A4")]

# Raw score ranking
grp = sub.groupby("family")["score"].agg(["mean", "std"]).sort_values("mean")
print("RAW score ranking (n=60, g=1.0, A4):")
for i, (fam, row) in enumerate(grp.iterrows()):
    tag = " *" if "lorentzian" in fam else ""
    print(f"  {i+1:2d}. {fam:40s}  mean={row['mean']:+.2f}  std={row['std']:.2f}{tag}")

# Check if normalized columns exist
if "score_norm" in df.columns:
    sub2 = df[(df["n"] == 60) & (df["gamma"] == 1.0) & (df["action_mode"] == "A4")]
    grp2 = sub2.groupby("family")["score_norm"].agg(["mean", "std"]).sort_values("mean")
    print("\nNORMALIZED score ranking:")
    for i, (fam, row) in enumerate(grp2.iterrows()):
        tag = " *" if "lorentzian" in fam else ""
        print(f"  {i+1:2d}. {fam:40s}  mean={row['mean']:+.4f}  std={row['std']:.4f}{tag}")
else:
    print("\nNo score_norm column in raw output (normalization happens post-hoc)")
    # Compute robust z-score manually
    import numpy as np
    scores = sub["score"].values
    median = np.median(scores)
    mad = np.median(np.abs(scores - median))
    if mad == 0:
        mad = 1.0
    print(f"Global median={median:.2f}, MAD={mad:.2f}")
    sub = sub.copy()
    sub["z"] = (sub["score"] - median) / (mad * 1.4826)
    grp3 = sub.groupby("family")["z"].agg(["mean", "std"]).sort_values("mean")
    print("\nManual robust z-score ranking:")
    for i, (fam, row) in enumerate(grp3.iterrows()):
        tag = " *" if "lorentzian" in fam else ""
        print(f"  {i+1:2d}. {fam:40s}  mean={row['mean']:+.4f}  std={row['std']:.4f}{tag}")
