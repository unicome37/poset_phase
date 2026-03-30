"""
F2 — Margin-aware turn-on refit runner

Goal:
- Upgrade the current winner-only low-N refit into a margin-aware turn-on
  stability test that can be compared directly to manuscript-safe onset claims.

Outputs:
- JSON + markdown report
- seed-level checkpoint logs
- partial snapshots for long runs
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from expanded_family_robustness import ALL_FAMILIES, compute_features, mahalanobis_score


DEFAULT_ONSET_LEVELS: dict[str, Any] = {
    "consecutive_n": 3,
    "manuscript_safe": {
        "min_rank1_count": None,
        "require_min_margin_positive": True,
        "require_ci95_lower_positive": True,
        "max_cond_sigma": 60.0,
    },
    "operational": {
        "min_rank1_count": None,
        "require_mean_margin_positive": True,
    },
    "winner_only": {
        "min_rank1_count": None,
    },
}


def _load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "PyYAML is required for F2 configs. Install with pip install pyyaml"
        ) from e
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _stable_name_offset(name: str) -> int:
    digest = hashlib.sha256(name.encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % 10000


def _ensemble_features(gen_fn, n: int, reps: int, seed_base: int) -> np.ndarray:
    feats = []
    for r in range(reps):
        p = gen_fn(n, seed=seed_base + r)
        feats.append(compute_features(p, n))
    return np.array(feats, dtype=float)


def _rank(scores: dict[str, float]) -> list[str]:
    return [k for k, _ in sorted(scores.items(), key=lambda kv: kv[1])]


def _find_blocks_in_order(n_grid: list[int], hit_map: dict[int, bool], min_len: int) -> list[list[int]]:
    blocks: list[list[int]] = []
    cur: list[int] = []
    for n in n_grid:
        if hit_map.get(n, False):
            cur.append(n)
        else:
            if len(cur) >= min_len:
                blocks.append(cur[:])
            cur = []
    if len(cur) >= min_len:
        blocks.append(cur[:])
    return blocks


def _ci95_bounds(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"mean": None, "std": None, "half_width": None, "lower": None, "upper": None}
    arr = np.array(values, dtype=float)
    mean = float(np.mean(arr))
    if len(arr) == 1:
        return {"mean": mean, "std": 0.0, "half_width": 0.0, "lower": mean, "upper": mean}
    std = float(np.std(arr, ddof=1))
    half_width = float(1.96 * std / math.sqrt(len(arr)))
    return {
        "mean": mean,
        "std": std,
        "half_width": half_width,
        "lower": mean - half_width,
        "upper": mean + half_width,
    }


def _resolve_family_names(cfg: dict[str, Any]) -> list[str]:
    registry = cfg.get("family_registry", "all")
    focus = list(cfg.get("focus_families", []))

    if registry == "all":
        family_names = sorted(ALL_FAMILIES.keys())
    elif registry == "focus_only":
        family_names = sorted(set(["Lor4D", *focus]))
    elif isinstance(registry, list):
        family_names = sorted(set(registry + ["Lor4D"]))
    else:
        raise RuntimeError(f"Unsupported family_registry={registry!r}")

    missing = [name for name in family_names if name not in ALL_FAMILIES]
    if missing:
        raise RuntimeError(f"Unknown families in registry: {missing}")
    if "Lor4D" not in family_names:
        raise RuntimeError("Lor4D missing in family registry")
    return family_names


def _merge_onset_levels(cfg: dict[str, Any], seed_runs: int) -> dict[str, Any]:
    onset_cfg = dict(DEFAULT_ONSET_LEVELS)
    user_cfg = cfg.get("onset_levels", {})
    onset_cfg.update(user_cfg)

    consecutive_n = int(onset_cfg.get("consecutive_n", 3))
    result = {"consecutive_n": consecutive_n}
    for level_name in ("manuscript_safe", "operational", "winner_only"):
        merged = dict(DEFAULT_ONSET_LEVELS[level_name])
        merged.update(onset_cfg.get(level_name, {}))
        if merged.get("min_rank1_count") is None:
            merged["min_rank1_count"] = seed_runs
        result[level_name] = merged
    return result


def _passes_level(stats: dict[str, Any], criteria: dict[str, Any]) -> bool:
    if stats["rank1_count"] < int(criteria.get("min_rank1_count", 0)):
        return False
    if criteria.get("require_mean_margin_positive", False) and not (stats["mean_margin"] is not None and stats["mean_margin"] > 0.0):
        return False
    if criteria.get("require_min_margin_positive", False) and not (stats["min_margin"] is not None and stats["min_margin"] > 0.0):
        return False
    if criteria.get("require_ci95_lower_positive", False):
        ci_lower = stats["margin_ci95"].get("lower")
        if not (ci_lower is not None and ci_lower > 0.0):
            return False
    max_cond_limit = criteria.get("max_cond_sigma")
    if max_cond_limit is not None:
        max_cond_sigma = stats.get("max_cond_sigma")
        if max_cond_sigma is None or max_cond_sigma >= float(max_cond_limit):
            return False
    return True


def _scores_for_seed(
    *,
    n: int,
    reps: int,
    reference_reps: int,
    seed_base: int,
    family_names: list[str],
    focus_families: list[str],
) -> dict[str, Any]:
    ref = _ensemble_features(ALL_FAMILIES["Lor4D"], n, reference_reps, seed_base + 100000)
    mu = ref.mean(axis=0)
    cov = np.cov(ref.T) + 1e-8 * np.eye(3)
    cov_inv = np.linalg.inv(cov)
    cond_sigma = float(np.linalg.cond(cov))

    scores: dict[str, float] = {}
    family_mean_features: dict[str, list[float]] = {}
    for name in family_names:
        feats = _ensemble_features(
            ALL_FAMILIES[name],
            n,
            reps,
            seed_base + 1000 + _stable_name_offset(name),
        )
        mean_feat = feats.mean(axis=0)
        family_mean_features[name] = [float(x) for x in mean_feat]
        scores[name] = float(mahalanobis_score(mean_feat, mu, cov_inv))

    ordered = _rank(scores)
    lor4d_rank = ordered.index("Lor4D") + 1
    lor4d_score = float(scores["Lor4D"])

    competitor_order = [name for name in ordered if name != "Lor4D"]
    best_competitor = competitor_order[0]
    best_competitor_score = float(scores[best_competitor])

    focus_pool = [name for name in focus_families if name in scores and name != "Lor4D"]
    if focus_pool:
        focus_best = min(focus_pool, key=lambda name: scores[name])
        focus_best_score = float(scores[focus_best])
    else:
        focus_best = None
        focus_best_score = None

    return {
        "winner": ordered[0],
        "winner_score": float(scores[ordered[0]]),
        "lor4d_rank": lor4d_rank,
        "lor4d_score": lor4d_score,
        "runner_up_family": best_competitor,
        "runner_up_score": best_competitor_score,
        "margin": float(best_competitor_score - lor4d_score),
        "cond_sigma": cond_sigma,
        "mu_d": float(mu[0]),
        "mu_c": float(mu[1]),
        "mu_w": float(mu[2]),
        "focus_runner_up_family": focus_best,
        "focus_runner_up_score": focus_best_score,
        "focus_margin": None if focus_best_score is None else float(focus_best_score - lor4d_score),
        "family_mean_features": family_mean_features,
    }


def run(cfg: dict[str, Any]) -> dict[str, Any]:
    out_dir = Path(cfg.get("output_dir", "outputs_carlip"))
    out_dir.mkdir(parents=True, exist_ok=True)
    output_basename = str(cfg.get("output_basename", "f2_turnon_margin_refit"))
    checkpoint_every_seed = bool(cfg.get("checkpoint_every_seed", True))

    experiment_id = str(cfg.get("experiment_id", "f2_turnon_margin_refit"))
    n_grid: list[int] = list(cfg["n_grid"])
    reps = int(cfg.get("reps", 80))
    reference_reps = int(cfg.get("reference_reps", reps))
    seed_base = int(cfg.get("seed_base", 42))
    seed_runs = int(cfg.get("seed_runs", cfg.get("hard_fail", {}).get("seed_runs", 10)))
    family_names = _resolve_family_names(cfg)
    focus_families = list(cfg.get("focus_families", []))
    onset_levels = _merge_onset_levels(cfg, seed_runs)

    seedlog_path = out_dir / f"{output_basename}.seedlog.jsonl"
    snapshot_path = out_dir / f"{output_basename}.partial.json"
    if checkpoint_every_seed:
        seedlog_path.write_text("", encoding="utf-8")

    per_seed_records: list[dict[str, Any]] = []
    per_n_records: dict[int, list[dict[str, Any]]] = defaultdict(list)

    for s in range(seed_runs):
        for n in n_grid:
            seed_stats = _scores_for_seed(
                n=n,
                reps=reps,
                reference_reps=reference_reps,
                seed_base=seed_base + 1000000 * s + 10000 * n,
                family_names=family_names,
                focus_families=focus_families,
            )
            record = {
                "seed_run": s,
                "N": n,
                **seed_stats,
            }
            per_seed_records.append(record)
            per_n_records[n].append(record)

        if checkpoint_every_seed:
            event = {
                "event": "seed_complete",
                "seed_run": s,
                "completed_seed_runs": s + 1,
                "total_seed_runs": seed_runs,
                "records_so_far": len(per_seed_records),
            }
            with seedlog_path.open("a", encoding="utf-8") as f:
                f.write(json.dumps(event, ensure_ascii=False) + "\n")

            partial = {
                "experiment_id": experiment_id,
                "output_basename": output_basename,
                "n_grid": n_grid,
                "reps": reps,
                "reference_reps": reference_reps,
                "completed_seed_runs": s + 1,
                "total_seed_runs": seed_runs,
                "focus_families": focus_families,
                "family_registry": cfg.get("family_registry", "all"),
                "per_seed_records": per_seed_records,
            }
            snapshot_path.write_text(json.dumps(partial, ensure_ascii=False, indent=2), encoding="utf-8")

    per_n_summary: dict[str, Any] = {}
    for n in n_grid:
        records = per_n_records[n]
        margins = [float(r["margin"]) for r in records]
        focus_margins = [float(r["focus_margin"]) for r in records if r["focus_margin"] is not None]
        conds = [float(r["cond_sigma"]) for r in records]
        rank1_count = sum(1 for r in records if int(r["lor4d_rank"]) == 1)
        winner_census = Counter(str(r["winner"]) for r in records)
        top_intruder = Counter(str(r["runner_up_family"]) for r in records)
        focus_intruder = Counter(str(r["focus_runner_up_family"]) for r in records if r["focus_runner_up_family"] is not None)
        margin_ci95 = _ci95_bounds(margins)
        focus_margin_ci95 = _ci95_bounds(focus_margins)

        per_n_summary[str(n)] = {
            "seed_runs": len(records),
            "rank1_count": rank1_count,
            "rank1_rate": float(rank1_count / len(records)) if records else None,
            "mean_margin": float(np.mean(margins)) if margins else None,
            "min_margin": float(np.min(margins)) if margins else None,
            "margin_ci95": margin_ci95,
            "mean_focus_margin": float(np.mean(focus_margins)) if focus_margins else None,
            "min_focus_margin": float(np.min(focus_margins)) if focus_margins else None,
            "focus_margin_ci95": focus_margin_ci95,
            "top_intruder_census": dict(sorted(top_intruder.items())),
            "focus_intruder_census": dict(sorted(focus_intruder.items())),
            "winner_census": dict(sorted(winner_census.items())),
            "max_cond_sigma": float(np.max(conds)) if conds else None,
            "mean_cond_sigma": float(np.mean(conds)) if conds else None,
            "mu_d_mean": float(np.mean([float(r["mu_d"]) for r in records])) if records else None,
            "mu_c_mean": float(np.mean([float(r["mu_c"]) for r in records])) if records else None,
            "mu_w_mean": float(np.mean([float(r["mu_w"]) for r in records])) if records else None,
        }

    onset_summary: dict[str, Any] = {"consecutive_n": onset_levels["consecutive_n"]}
    for level_name in ("manuscript_safe", "operational", "winner_only"):
        criteria = onset_levels[level_name]
        pass_map = {
            n: _passes_level(per_n_summary[str(n)], criteria)
            for n in n_grid
        }
        blocks = _find_blocks_in_order(n_grid, pass_map, int(onset_levels["consecutive_n"]))
        onset_summary[level_name] = {
            "criteria": criteria,
            "pass_map": {str(k): bool(v) for k, v in pass_map.items()},
            "blocks": blocks,
            "first_block_start": None if not blocks else blocks[0][0],
        }

    result = {
        "experiment_id": experiment_id,
        "output_basename": output_basename,
        "family_registry": cfg.get("family_registry", "all"),
        "family_names": family_names,
        "focus_families": focus_families,
        "n_grid": n_grid,
        "reps": reps,
        "reference_reps": reference_reps,
        "seed_runs": seed_runs,
        "onset_levels": onset_levels,
        "per_n_summary": per_n_summary,
        "onset_summary": onset_summary,
        "per_seed_records": per_seed_records,
    }

    json_path = out_dir / f"{output_basename}.json"
    md_path = out_dir / f"{output_basename}.md"

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    if checkpoint_every_seed:
        event = {
            "event": "run_complete",
            "completed_seed_runs": seed_runs,
            "winner_only_start": onset_summary["winner_only"]["first_block_start"],
            "operational_start": onset_summary["operational"]["first_block_start"],
            "manuscript_safe_start": onset_summary["manuscript_safe"]["first_block_start"],
        }
        with seedlog_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(event, ensure_ascii=False) + "\n")

    lines = [
        "# F2 Margin-Aware Turn-on Refit",
        "",
        f"- Winner-only onset start: **{onset_summary['winner_only']['first_block_start']}**",
        f"- Operational onset start: **{onset_summary['operational']['first_block_start']}**",
        f"- Manuscript-safe onset start: **{onset_summary['manuscript_safe']['first_block_start']}**",
        f"- N grid: {n_grid}",
        f"- Seeds: {seed_runs}",
        f"- Reps: {reps}",
        f"- Reference reps: {reference_reps}",
        "",
        "## Per-N summary",
        "",
        "| N | rank1 | rate | mean margin | min margin | CI95 lower | CI95 upper | max cond(Σ) | top intruder |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]

    for n in n_grid:
        stats = per_n_summary[str(n)]
        top_intruder = "-"
        if stats["top_intruder_census"]:
            fam, count = max(stats["top_intruder_census"].items(), key=lambda kv: kv[1])
            top_intruder = f"{fam} ({count})"
        lines.append(
            "| {n} | {rank1} | {rate:.3f} | {mean_margin:+.6f} | {min_margin:+.6f} | {ci_lower:+.6f} | {ci_upper:+.6f} | {max_cond:.3f} | {top_intruder} |".format(
                n=n,
                rank1=stats["rank1_count"],
                rate=stats["rank1_rate"],
                mean_margin=stats["mean_margin"],
                min_margin=stats["min_margin"],
                ci_lower=stats["margin_ci95"]["lower"],
                ci_upper=stats["margin_ci95"]["upper"],
                max_cond=stats["max_cond_sigma"],
                top_intruder=top_intruder,
            )
        )

    lines += ["", "## Onset blocks", ""]
    for level_name in ("winner_only", "operational", "manuscript_safe"):
        info = onset_summary[level_name]
        lines.append(f"- {level_name}: start={info['first_block_start']}, blocks={info['blocks']}")

    md_path.write_text("\n".join(lines), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_f2_turnon_margin_refit.yaml")
    args = parser.parse_args()

    cfg = _load_yaml(Path(args.config))
    result = run(cfg)
    print("DONE")
    print(f"outputs_carlip/{result['output_basename']}.json")


if __name__ == "__main__":
    main()
