"""
F1 — Family pressure falsification runner (skeleton)

Goal:
- Test C1 hard-fail condition on expanded family library.
- Produce machine-readable JSON + markdown summary.

Hard-fail (pre-registered):
- A non-Lor4D family outranks Lor4D in >=3 consecutive N bins,
  with >=8/10 seed-runs per bin.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from expanded_family_robustness import ALL_FAMILIES, compute_features, mahalanobis_score


@dataclass
class C1Rule:
    consecutive_n: int
    min_seed_success: int
    seed_runs: int


def _load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "PyYAML is required for config_falsify_c1.yaml. Install with pip install pyyaml"
        ) from e
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _ensemble_features(gen_fn, n: int, reps: int, seed_base: int) -> np.ndarray:
    feats = []
    for r in range(reps):
        p = gen_fn(n, seed=seed_base + r)
        feats.append(compute_features(p, n))
    return np.array(feats)


def _scores_for_seed(n: int, reps: int, seed_base: int, family_names: list[str]) -> dict[str, float]:
    # Lor4D reference for this seed-run
    ref = _ensemble_features(ALL_FAMILIES["Lor4D"], n, reps, seed_base + 100000)
    mu = ref.mean(axis=0)
    cov = np.cov(ref.T) + 1e-8 * np.eye(3)
    cov_inv = np.linalg.inv(cov)

    scores: dict[str, float] = {}
    for name in family_names:
        feats = _ensemble_features(ALL_FAMILIES[name], n, reps, seed_base + 1000 + hash(name) % 10000)
        scores[name] = mahalanobis_score(feats.mean(axis=0), mu, cov_inv)
    return scores


def _rank(scores: dict[str, float]) -> list[str]:
    return [k for k, _ in sorted(scores.items(), key=lambda kv: kv[1])]


def _find_consecutive(ns: list[int], hit_map: dict[int, int], min_hits: int, k: int) -> list[list[int]]:
    blocks: list[list[int]] = []
    cur: list[int] = []
    for n in sorted(ns):
        if hit_map.get(n, 0) >= min_hits:
            if not cur or n > cur[-1]:
                cur.append(n)
            else:
                # should not happen for unique sorted N, but keep robust
                if len(cur) >= k:
                    blocks.append(cur[:])
                cur = [n]
        else:
            if len(cur) >= k:
                blocks.append(cur[:])
            cur = []
    if len(cur) >= k:
        blocks.append(cur[:])
    return blocks


def run(cfg: dict[str, Any]) -> dict[str, Any]:
    out_dir = Path(cfg.get("output_dir", "outputs_carlip"))
    out_dir.mkdir(parents=True, exist_ok=True)
    output_basename = str(cfg.get("output_basename", "falsify_c1_family_pressure"))
    checkpoint_every_seed = bool(cfg.get("checkpoint_every_seed", True))

    n_grid: list[int] = list(cfg["n_grid"])
    reps: int = int(cfg.get("reps", 80))
    seed_base: int = int(cfg.get("seed_base", 42))

    hf = cfg["hard_fail"]
    rule = C1Rule(
        consecutive_n=int(hf["consecutive_n"]),
        min_seed_success=int(hf["min_seed_success"]),
        seed_runs=int(hf["seed_runs"]),
    )

    family_names = sorted(list(ALL_FAMILIES.keys()))
    if "Lor4D" not in family_names:
        raise RuntimeError("Lor4D missing in ALL_FAMILIES")

    seedlog_path = out_dir / f"{output_basename}.seedlog.jsonl"
    snapshot_path = out_dir / f"{output_basename}.partial.json"
    if checkpoint_every_seed:
        seedlog_path.write_text("", encoding="utf-8")

    # per family, per N: count seed-runs where family outranks Lor4D
    outrank_counts: dict[str, dict[int, int]] = {
        f: {n: 0 for n in n_grid} for f in family_names if f != "Lor4D"
    }

    # run seed-runs
    per_seed_records = []
    for s in range(rule.seed_runs):
        for n in n_grid:
            scores = _scores_for_seed(n=n, reps=reps, seed_base=seed_base + 1000000 * s + 10000 * n, family_names=family_names)
            ordered = _rank(scores)
            lor_rank = ordered.index("Lor4D") + 1
            per_seed_records.append({"seed_run": s, "N": n, "lor4d_rank": lor_rank, "winner": ordered[0]})
            for fam in outrank_counts:
                if ordered.index(fam) < ordered.index("Lor4D"):
                    outrank_counts[fam][n] += 1

        if checkpoint_every_seed:
            event = {
                "event": "seed_complete",
                "seed_run": s,
                "completed_seed_runs": s + 1,
                "total_seed_runs": rule.seed_runs,
                "records_so_far": len(per_seed_records),
            }
            with seedlog_path.open("a", encoding="utf-8") as f:
                f.write(json.dumps(event, ensure_ascii=False) + "\n")

            partial = {
                "experiment_id": cfg.get("experiment_id", "falsify_c1_family_pressure"),
                "output_basename": output_basename,
                "rule": {
                    "consecutive_n": rule.consecutive_n,
                    "min_seed_success": rule.min_seed_success,
                    "seed_runs": rule.seed_runs,
                },
                "n_grid": n_grid,
                "reps": reps,
                "completed_seed_runs": s + 1,
                "total_seed_runs": rule.seed_runs,
                "outrank_counts": outrank_counts,
                "per_seed_records": per_seed_records,
            }
            snapshot_path.write_text(json.dumps(partial, ensure_ascii=False, indent=2), encoding="utf-8")

    # evaluate hard-fail
    hard_fail_families = []
    for fam, hit_map in outrank_counts.items():
        blocks = _find_consecutive(n_grid, hit_map, rule.min_seed_success, rule.consecutive_n)
        if blocks:
            hard_fail_families.append({"family": fam, "blocks": blocks, "hit_map": hit_map})

    hard_fail = len(hard_fail_families) > 0

    result = {
        "experiment_id": cfg.get("experiment_id", "falsify_c1_family_pressure"),
        "output_basename": output_basename,
        "rule": {
            "consecutive_n": rule.consecutive_n,
            "min_seed_success": rule.min_seed_success,
            "seed_runs": rule.seed_runs,
        },
        "n_grid": n_grid,
        "reps": reps,
        "hard_fail": hard_fail,
        "hard_fail_families": hard_fail_families,
        "outrank_counts": outrank_counts,
        "per_seed_records": per_seed_records,
    }

    # write outputs
    json_path = out_dir / f"{output_basename}.json"
    md_path = out_dir / f"{output_basename}.md"

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    if checkpoint_every_seed:
        event = {
            "event": "run_complete",
            "completed_seed_runs": rule.seed_runs,
            "hard_fail": hard_fail,
        }
        with seedlog_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(event, ensure_ascii=False) + "\n")

    lines = [
        "# F1 Falsification Report — Family Pressure",
        "",
        f"- Hard fail: **{'YES' if hard_fail else 'NO'}**",
        f"- Rule: consecutive_N >= {rule.consecutive_n}, seed_hits >= {rule.min_seed_success}/{rule.seed_runs}",
        f"- N grid: {n_grid}",
        "",
        "## Families triggering hard fail",
    ]
    if hard_fail_families:
        for item in hard_fail_families:
            lines.append(f"- {item['family']}: blocks={item['blocks']}")
    else:
        lines.append("- None")

    lines += ["", "## Notes", "", "- This is a skeleton runner; extend family registry as new adversarial generators are added."]

    md_path.write_text("\n".join(lines), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_falsify_c1.yaml")
    args = parser.parse_args()

    cfg = _load_yaml(Path(args.config))
    result = run(cfg)
    print("HARD_FAIL" if result["hard_fail"] else "PASS")
    output_basename = str(cfg.get("output_basename", "falsify_c1_family_pressure"))
    print(f"outputs_carlip/{output_basename}.json")


if __name__ == "__main__":
    main()
