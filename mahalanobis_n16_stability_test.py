"""
N=16 Stability Test for Mahalanobis Variants
============================================

Goal:
  Evaluate whether a "principled" small-N modification can make Lor4D rank #1
  at N=16 consistently across seeds, without introducing hand-tuned knobs.

Methods compared (all parameter-free given the Lor4D ensemble at that N):
  1) Baseline(full):   (I - mu)^T Sigma^{-1} (I - mu)
  2) Anchored(full):   (I - mu_anchor)^T Sigma^{-1} (I - mu_anchor),
                       mu_anchor = (4.0, mu_c, mu_w)
  3) Anchored(diag):   sum_i (I_i - mu_anchor_i)^2 / var_i

We test over multiple seed bases, keeping the generator/seed protocol
consistent with other scripts in this project.
"""
from __future__ import annotations

import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from expanded_family_robustness import ALL_FAMILIES, compute_features


def mahalanobis_score(feat: np.ndarray, mu: np.ndarray, cov_inv: np.ndarray) -> float:
    delta = feat - mu
    return float(delta @ cov_inv @ delta)


def main() -> None:
    N = 16
    REPS = 80
    SEED_BASES = [42, 137, 271, 500, 777, 1001, 2023, 3141, 5000, 8888]

    print("=" * 80)
    print("N=16 Stability Test — Mahalanobis Variants")
    print(f"  Families: {len(ALL_FAMILIES)}")
    print(f"  Seeds: {SEED_BASES}")
    print(f"  Reps per family: {REPS}")
    print("=" * 80)

    t0 = time.time()

    winners = {"baseline_full": Counter(), "anchored_full": Counter(), "anchored_diag": Counter()}
    lor_ranks = {"baseline_full": [], "anchored_full": [], "anchored_diag": []}
    lor_margins = {"baseline_full": [], "anchored_full": [], "anchored_diag": []}

    per_seed_rows: list[dict] = []

    for seed_base in SEED_BASES:
        data: dict[str, list[np.ndarray]] = defaultdict(list)
        for fam_name, gen_fn in ALL_FAMILIES.items():
            for rep in range(REPS):
                seed = seed_base + hash(fam_name) % 10000 + N * 100 + rep
                seed = seed % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    feat = compute_features(poset, N)
                    data[fam_name].append(feat)
                except Exception:
                    pass

        lor4d_feats = np.array(data["Lor4D"])
        if len(lor4d_feats) < 5:
            print(f"  WARNING: seed_base={seed_base} has too few Lor4D samples ({len(lor4d_feats)})")
            continue

        mu = np.mean(lor4d_feats, axis=0)
        cov = np.cov(lor4d_feats.T)
        cov_inv = np.linalg.inv(cov + 1e-12 * np.eye(3))

        mu_anchor = mu.copy()
        mu_anchor[0] = 4.0
        var = np.var(lor4d_feats, axis=0, ddof=1)
        inv_var = 1.0 / np.maximum(var, 1e-12)

        mean_scores = {"baseline_full": {}, "anchored_full": {}, "anchored_diag": {}}
        for fam_name, feats in data.items():
            if not feats:
                continue
            feats_arr = np.array(feats)
            mean_scores["baseline_full"][fam_name] = float(
                np.mean([mahalanobis_score(f, mu, cov_inv) for f in feats_arr])
            )
            mean_scores["anchored_full"][fam_name] = float(
                np.mean([mahalanobis_score(f, mu_anchor, cov_inv) for f in feats_arr])
            )
            deltas = feats_arr - mu_anchor
            mean_scores["anchored_diag"][fam_name] = float(np.mean(np.sum(inv_var * (deltas ** 2), axis=1)))

        for key in ["baseline_full", "anchored_full", "anchored_diag"]:
            ranked = sorted(mean_scores[key], key=mean_scores[key].get)
            w = ranked[0]
            winners[key][w] += 1
            r_lor = ranked.index("Lor4D") + 1 if "Lor4D" in ranked else 99
            lor_ranks[key].append(r_lor)
            # margin vs best non-Lor4D
            best_nonlor = min(mean_scores[key][f] for f in ranked if f != "Lor4D")
            margin = best_nonlor - mean_scores[key]["Lor4D"]
            lor_margins[key].append(margin)

        per_seed_rows.append(
            {
                "seed_base": seed_base,
                "mu_d": float(mu[0]),
                "mu_c": float(mu[1]),
                "mu_w": float(mu[2]),
                "winner_baseline": min(mean_scores["baseline_full"], key=mean_scores["baseline_full"].get),
                "winner_anchor_full": min(mean_scores["anchored_full"], key=mean_scores["anchored_full"].get),
                "winner_anchor_diag": min(mean_scores["anchored_diag"], key=mean_scores["anchored_diag"].get),
                "rank_baseline": sorted(mean_scores["baseline_full"], key=mean_scores["baseline_full"].get).index("Lor4D") + 1,
                "rank_anchor_full": sorted(mean_scores["anchored_full"], key=mean_scores["anchored_full"].get).index("Lor4D") + 1,
                "rank_anchor_diag": sorted(mean_scores["anchored_diag"], key=mean_scores["anchored_diag"].get).index("Lor4D") + 1,
            }
        )

        print(
            f"  seed={seed_base:4d}: winners "
            f"baseline={per_seed_rows[-1]['winner_baseline']}, "
            f"anchor_full={per_seed_rows[-1]['winner_anchor_full']}, "
            f"anchor_diag={per_seed_rows[-1]['winner_anchor_diag']}"
        )

    # ── Report ───────────────────────────────────────────────────────────
    report: list[str] = []
    report.append("# N=16 Stability Test — Mahalanobis Variants\n")
    report.append(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    report.append(f"N: {N}")
    report.append(f"Families: {len(ALL_FAMILIES)}")
    report.append(f"Seeds tested: {len(SEED_BASES)}")
    report.append(f"Reps per family per seed: {REPS}\n")

    report.append("## 1. Winner Census (by seed)\n")
    for key, title in [
        ("baseline_full", "Baseline (full)"),
        ("anchored_full", "Anchored (full)"),
        ("anchored_diag", "Anchored (diag)"),
    ]:
        report.append(f"### {title}\n")
        report.append("| Winner | Count |")
        report.append("|---|---:|")
        for fam, cnt in winners[key].most_common():
            report.append(f"| {fam} | {cnt} |")
        report.append("")

    report.append("## 2. Lor4D Rank and Margin Summary\n")
    report.append("| Method | Mean rank | #1 rate | Mean margin |")
    report.append("|---|---:|---:|---:|")
    for key, title in [
        ("baseline_full", "baseline_full"),
        ("anchored_full", "anchored_full"),
        ("anchored_diag", "anchored_diag"),
    ]:
        ranks = np.array(lor_ranks[key], dtype=float)
        margins = np.array(lor_margins[key], dtype=float)
        mean_rank = float(np.mean(ranks)) if len(ranks) else float("nan")
        rate1 = float(np.mean(ranks == 1.0)) if len(ranks) else float("nan")
        mean_margin = float(np.mean(margins)) if len(margins) else float("nan")
        report.append(f"| {title} | {mean_rank:.2f} | {rate1*100:.0f}% | {mean_margin:.4f} |")

    report.append("\n## 3. Per-Seed Detail\n")
    report.append("| seed | mu_d | mu_c | mu_w | winner baseline | winner anchor_full | winner anchor_diag | Lor4D rank baseline | rank anchor_full | rank anchor_diag |")
    report.append("|---:|---:|---:|---:|---|---|---|---:|---:|---:|")
    for row in per_seed_rows:
        report.append(
            f"| {row['seed_base']} | {row['mu_d']:.4f} | {row['mu_c']:.4f} | {row['mu_w']:.4f} | "
            f"{row['winner_baseline']} | {row['winner_anchor_full']} | {row['winner_anchor_diag']} | "
            f"{row['rank_baseline']} | {row['rank_anchor_full']} | {row['rank_anchor_diag']} |"
        )

    report.append(f"\nRuntime: {time.time() - t0:.1f}s\n")

    outdir = Path(__file__).parent / "outputs_carlip"
    outdir.mkdir(exist_ok=True)
    outpath = outdir / "mahalanobis_n16_stability_results.md"
    outpath.write_text("\n".join(report), encoding="utf-8")
    print(f"\nSaved to {outpath}")


if __name__ == "__main__":
    main()
