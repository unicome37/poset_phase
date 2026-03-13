"""w_con sensitivity scan for Prediction A.

Tests 4D dominance robustness across a range of consistency weights.
Uses the exact computation range (N=20-36) for speed and accuracy.
"""
import subprocess
import sys
import yaml
import pandas as pd
from pathlib import Path

WEIGHTS = [0.3, 0.4, 0.5, 0.6, 0.682882, 0.8, 1.0, 1.2]
TEMPLATE = Path("config_wcon_sensitivity_template.yaml")
OUT_DIR = Path("outputs_exploratory/wcon_sensitivity")
OUT_DIR.mkdir(parents=True, exist_ok=True)

results = []

for w in WEIGHTS:
    tag = f"w{w:.4f}".replace(".", "p")
    cfg_path = OUT_DIR / f"config_{tag}.yaml"
    run_dir = OUT_DIR / tag

    # Read winner counts (check if already completed)
    winners_csv = run_dir / "prediction_a_winners.csv"
    if winners_csv.exists() and winners_csv.stat().st_size > 0:
        print(f"\n{'='*60}")
        print(f"SKIP w_con = {w:.4f} (already completed)")
        print(f"{'='*60}")
    else:
        # Write config
        with open(TEMPLATE) as f:
            cfg = yaml.safe_load(f)
        cfg["experiment"]["consistency_weight"] = w
        cfg["experiment"]["multi_consistency_weight"] = w
        cfg["output"]["directory"] = str(run_dir)
        with open(cfg_path, "w") as f:
            yaml.dump(cfg, f)

        # Run experiment
        print(f"\n{'='*60}")
        print(f"Running w_con = {w:.4f} ...")
        print(f"{'='*60}")
        subprocess.run([sys.executable, "prediction_a_dimension_scan.py", "--config", str(cfg_path)],
                       check=True)

    if winners_csv.exists():
        df = pd.read_csv(winners_csv)
        counts = df["winner_family"].value_counts()
        total = len(df)
        lor4d_wins = int(counts.get("lorentzian_like_4d", 0))
        kr_wins = int(counts.get("KR_like", 0))
        lor2d_wins = int(counts.get("lorentzian_like_2d", 0))
        lor3d_wins = int(counts.get("lorentzian_like_3d", 0))
        results.append({
            "w_con": w,
            "total_configs": total,
            "lor4d_wins": lor4d_wins,
            "kr_wins": kr_wins,
            "lor2d_wins": lor2d_wins,
            "lor3d_wins": lor3d_wins,
            "lor4d_pct": f"{100*lor4d_wins/total:.1f}%"
        })
        print(f"  w_con={w:.4f}: Lor4D={lor4d_wins}/{total} ({100*lor4d_wins/total:.1f}%)")

# Summary
print("\n" + "="*60)
print("w_con SENSITIVITY SUMMARY")
print("="*60)
summary_df = pd.DataFrame(results)
print(summary_df.to_string(index=False))
summary_df.to_csv(OUT_DIR / "wcon_sensitivity_summary.csv", index=False)
print(f"\nSaved to {OUT_DIR / 'wcon_sensitivity_summary.csv'}")
