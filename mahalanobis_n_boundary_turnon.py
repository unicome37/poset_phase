"""
Mahalanobis S_MD Turn-On Experiment across N Boundary
=====================================================

Purpose:
  Produce a definitive "turn-on table" for the Letter:
  At each N in {12, 14, 16, 18, 20, 24, 28, 32}, with REPS=80 reference
  ensemble, does Lor4D rank #1 across all seeds?

  This answers: what is the minimum N at which S_MD unconditionally identifies
  Lor4D as the unique minimiser?

Protocol:
  - 25 families (17 original + 8 adversarial)
  - 10 independent seed bases
  - REPS=80 samples per family per (N, seed) → Lor4D reference μ(N), Σ(N)
  - Score = (I - μ)^T Σ^{-1} (I - μ)  (zero free parameters)
  - Record: Lor4D rank, margin to runner-up, runner-up identity, cond(Σ)

Output:
  outputs_carlip/mahalanobis_n_boundary_turnon.md
"""
from __future__ import annotations

import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from expanded_family_robustness import ALL_FAMILIES, compute_features


def mahalanobis_score(feat: np.ndarray, mu: np.ndarray, cov_inv: np.ndarray) -> float:
    delta = feat - mu
    return float(delta @ cov_inv @ delta)


def run_single(N: int, seed_base: int, reps: int) -> dict:
    """Run one (N, seed_base) experiment. Returns result dict."""
    data: dict[str, list[np.ndarray]] = defaultdict(list)
    fail_count = 0

    for fam_name, gen_fn in ALL_FAMILIES.items():
        for rep in range(reps):
            seed = (seed_base + hash(fam_name) % 10000 + N * 100 + rep) % (2**31)
            try:
                poset = gen_fn(N, seed=seed)
                feat = compute_features(poset, N)
                data[fam_name].append(feat)
            except Exception:
                fail_count += 1

    lor4d_feats = np.array(data.get("Lor4D", []))
    if len(lor4d_feats) < 5:
        return {"valid": False, "N": N, "seed_base": seed_base}

    mu = np.mean(lor4d_feats, axis=0)
    cov = np.cov(lor4d_feats.T)
    reg = 1e-12 * np.eye(3)
    cov_inv = np.linalg.inv(cov + reg)
    cond_num = float(np.linalg.cond(cov))

    # Score every family
    mean_scores: dict[str, float] = {}
    for fam_name, feats in data.items():
        if not feats:
            continue
        feats_arr = np.array(feats)
        scores = [mahalanobis_score(f, mu, cov_inv) for f in feats_arr]
        mean_scores[fam_name] = float(np.mean(scores))

    ranked = sorted(mean_scores, key=mean_scores.get)
    lor_rank = ranked.index("Lor4D") + 1 if "Lor4D" in ranked else 99

    # Runner-up: best non-Lor4D family
    runner_up = next((f for f in ranked if f != "Lor4D"), "N/A")
    runner_up_score = mean_scores.get(runner_up, float("inf"))
    lor_score = mean_scores.get("Lor4D", float("inf"))
    margin = runner_up_score - lor_score

    return {
        "valid": True,
        "N": N,
        "seed_base": seed_base,
        "lor_rank": lor_rank,
        "winner": ranked[0],
        "runner_up": runner_up,
        "margin": margin,
        "lor_score": lor_score,
        "runner_up_score": runner_up_score,
        "cond_sigma": cond_num,
        "mu_d": float(mu[0]),
        "mu_c": float(mu[1]),
        "mu_w": float(mu[2]),
        "n_lor4d_samples": len(lor4d_feats),
        "n_families_tested": len(mean_scores),
        "fail_count": fail_count,
    }


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    N_GRID = [12, 14, 16, 18, 20, 24, 28, 32]
    REPS = 80
    SEED_BASES = [42, 137, 271, 500, 777, 1001, 2023, 3141, 5000, 8888]

    print("=" * 80)
    print("Mahalanobis S_MD Turn-On Experiment")
    print(f"  N grid: {N_GRID}")
    print(f"  Families: {len(ALL_FAMILIES)}")
    print(f"  Seeds: {len(SEED_BASES)}")
    print(f"  REPS: {REPS}")
    print(f"  Total runs: {len(N_GRID) * len(SEED_BASES)}")
    print("=" * 80)

    t0 = time.time()
    all_results: list[dict] = []

    for N in N_GRID:
        print(f"\n--- N = {N} ---")
        for seed_base in SEED_BASES:
            t1 = time.time()
            res = run_single(N, seed_base, REPS)
            dt = time.time() - t1
            all_results.append(res)
            if res["valid"]:
                tag = "OK" if res["lor_rank"] == 1 else f"RANK={res['lor_rank']}"
                print(
                    f"  seed={seed_base:5d}  rank={res['lor_rank']}  "
                    f"margin={res['margin']:+.4f}  runner={res['runner_up']:12s}  "
                    f"cond={res['cond_sigma']:.1f}  [{tag}]  {dt:.1f}s"
                )
            else:
                print(f"  seed={seed_base:5d}  INVALID (too few Lor4D samples)")

    elapsed = time.time() - t0

    # ── Aggregate by N ──────────────────────────────────────────────────
    summary_by_n: dict[int, dict] = {}
    for N in N_GRID:
        rows = [r for r in all_results if r["N"] == N and r["valid"]]
        if not rows:
            summary_by_n[N] = {"n_valid": 0}
            continue
        ranks = [r["lor_rank"] for r in rows]
        margins = [r["margin"] for r in rows]
        conds = [r["cond_sigma"] for r in rows]
        runners = Counter(r["runner_up"] for r in rows if r["lor_rank"] != 1)
        winners = Counter(r["winner"] for r in rows)
        summary_by_n[N] = {
            "n_valid": len(rows),
            "n_rank1": sum(1 for r in ranks if r == 1),
            "mean_rank": float(np.mean(ranks)),
            "max_rank": max(ranks),
            "mean_margin": float(np.mean(margins)),
            "min_margin": float(np.min(margins)),
            "max_margin": float(np.max(margins)),
            "std_margin": float(np.std(margins)),
            "mean_cond": float(np.mean(conds)),
            "max_cond": float(np.max(conds)),
            "intruders": dict(runners.most_common(3)),
            "winner_census": dict(winners.most_common(3)),
        }

    # ── Report ──────────────────────────────────────────────────────────
    report: list[str] = []
    report.append("# Mahalanobis S_MD — N-Boundary Turn-On Table\n")
    report.append(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    report.append(f"Families: {len(ALL_FAMILIES)} (17 standard + 8 adversarial)")
    report.append(f"Seeds: {len(SEED_BASES)}")
    report.append(f"REPS: {REPS}")
    report.append(f"Runtime: {elapsed:.0f}s\n")

    report.append("## 1. Turn-On Summary Table\n")
    report.append("| N | #1 rate | Mean rank | Mean margin | Min margin | Max cond(Σ) | Top intruder |")
    report.append("|---:|---:|---:|---:|---:|---:|---|")
    for N in N_GRID:
        s = summary_by_n[N]
        if s["n_valid"] == 0:
            report.append(f"| {N} | — | — | — | — | — | — |")
            continue
        rate = f"{s['n_rank1']}/{s['n_valid']}"
        intruder = list(s["intruders"].keys())[0] if s["intruders"] else "—"
        report.append(
            f"| {N} | {rate} | {s['mean_rank']:.2f} | "
            f"{s['mean_margin']:+.4f} | {s['min_margin']:+.4f} | "
            f"{s['max_cond']:.1f} | {intruder} |"
        )

    report.append("\n## 2. Interpretation\n")
    # Auto-detect turn-on boundary
    turn_on_n = None
    for N in N_GRID:
        s = summary_by_n[N]
        if s["n_valid"] > 0 and s["n_rank1"] == s["n_valid"]:
            turn_on_n = N
            break
    if turn_on_n:
        report.append(f"**Turn-on boundary**: N ≥ {turn_on_n} — Lor4D is unconditionally rank #1 "
                       f"across all {len(SEED_BASES)} seeds with REPS={REPS}.\n")
    else:
        report.append("**No clean turn-on found** within the tested N range.\n")

    # Check if all N achieve 100%
    all_perfect = all(
        summary_by_n[N]["n_valid"] > 0 and summary_by_n[N]["n_rank1"] == summary_by_n[N]["n_valid"]
        for N in N_GRID
    )
    if all_perfect:
        report.append("**Remarkable**: With REPS=80, Lor4D holds rank #1 at every N in the grid, "
                       "including the previously problematic N=12-16 range.\n")

    report.append("## 3. Reference Manifold Parameters μ(N)\n")
    report.append("| N | mean μ_d | mean μ_c | mean μ_w |")
    report.append("|---:|---:|---:|---:|")
    for N in N_GRID:
        rows = [r for r in all_results if r["N"] == N and r["valid"]]
        if not rows:
            report.append(f"| {N} | — | — | — |")
            continue
        mu_d = np.mean([r["mu_d"] for r in rows])
        mu_c = np.mean([r["mu_c"] for r in rows])
        mu_w = np.mean([r["mu_w"] for r in rows])
        report.append(f"| {N} | {mu_d:.4f} | {mu_c:.4f} | {mu_w:.4f} |")

    report.append("\n## 4. Margin Growth with N\n")
    report.append("| N | Mean margin | Std margin | Min margin | Max margin |")
    report.append("|---:|---:|---:|---:|---:|")
    for N in N_GRID:
        s = summary_by_n[N]
        if s["n_valid"] == 0:
            report.append(f"| {N} | — | — | — | — |")
            continue
        report.append(
            f"| {N} | {s['mean_margin']:.4f} | {s['std_margin']:.4f} | "
            f"{s['min_margin']:.4f} | {s['max_margin']:.4f} |"
        )

    report.append("\n## 5. Winner Census by N\n")
    for N in N_GRID:
        s = summary_by_n[N]
        if s["n_valid"] == 0:
            continue
        report.append(f"\n### N = {N}\n")
        report.append("| Family | Wins |")
        report.append("|---|---:|")
        for fam, cnt in sorted(s["winner_census"].items(), key=lambda x: -x[1]):
            report.append(f"| {fam} | {cnt} |")

    report.append("\n## 6. Per-Run Detail\n")
    report.append("| N | seed | Lor4D rank | margin | runner-up | cond(Σ) | μ_d | μ_c | μ_w |")
    report.append("|---:|---:|---:|---:|---|---:|---:|---:|---:|")
    for r in all_results:
        if not r["valid"]:
            report.append(f"| {r['N']} | {r['seed_base']} | INVALID | — | — | — | — | — | — |")
            continue
        report.append(
            f"| {r['N']} | {r['seed_base']} | {r['lor_rank']} | "
            f"{r['margin']:+.4f} | {r['runner_up']} | {r['cond_sigma']:.1f} | "
            f"{r['mu_d']:.4f} | {r['mu_c']:.4f} | {r['mu_w']:.4f} |"
        )

    report.append(f"\n---\n*Generated by mahalanobis_n_boundary_turnon.py — {time.strftime('%Y-%m-%d %H:%M')}*\n")

    # Save
    outdir = Path(__file__).parent / "outputs_carlip"
    outdir.mkdir(exist_ok=True)
    outpath = outdir / "mahalanobis_n_boundary_turnon.md"
    outpath.write_text("\n".join(report), encoding="utf-8")
    print(f"\n{'=' * 80}")
    print(f"Saved: {outpath}")
    print(f"Total runtime: {elapsed:.0f}s")


if __name__ == "__main__":
    main()
