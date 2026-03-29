"""
F3 — Background-response falsification runner (skeleton)

Goal:
- Test C2 hard-fail condition on weak/moderate curvature windows.
- Use metric-faithful generators where available; keep proxy branch explicit.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from curvature_backgrounds import (
    sprinkle_flrw_matter_diamond,
    poset_from_flrw_matter_points,
    sprinkle_schwarzschild_diamond,
    poset_from_schwarzschild_points,
)
from curvature_layer2_de_sitter_scan import (
    sprinkle_de_sitter_like_diamond,
    poset_from_de_sitter_points,
)
from expanded_family_robustness import ALL_FAMILIES, compute_features, mahalanobis_score


@dataclass
class C2Rule:
    min_fail_ratio: float
    top_k_required: int
    seed_runs: int


def _load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "PyYAML is required for config_falsify_c3.yaml. Install with pip install pyyaml"
        ) from e
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _ensemble_score(gen_fn, n: int, reps: int, mu: np.ndarray, cov_inv: np.ndarray, seed_base: int) -> float:
    feats = []
    for r in range(reps):
        p = gen_fn(n, seed=seed_base + r)
        feats.append(compute_features(p, n))
    return mahalanobis_score(np.array(feats).mean(axis=0), mu, cov_inv)


def _ref_stats(n: int, reps: int, seed: int) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    ref_feats = []
    for r in range(reps):
        p = ALL_FAMILIES["Lor4D"](n, seed=seed + r)
        ref_feats.append(compute_features(p, n))
    ref = np.array(ref_feats)
    mu = ref.mean(axis=0)
    cov = np.cov(ref.T) + 1e-8 * np.eye(3)
    cov_inv = np.linalg.inv(cov)

    base_scores = {}
    for name, gen_fn in ALL_FAMILIES.items():
        base_scores[name] = _ensemble_score(gen_fn, n, reps, mu, cov_inv, seed + 10000 + hash(name) % 10000)
    return mu, cov_inv, base_scores


def _make_ds_gen(h: float):
    def _gen(n: int, seed: int):
        pts = sprinkle_de_sitter_like_diamond(n=n, d_spatial=3, hubble=h, seed=seed)
        return poset_from_de_sitter_points(pts, hubble=h)

    return _gen


def _make_flrw_gen(kappa: float):
    def _gen(n: int, seed: int):
        pts = sprinkle_flrw_matter_diamond(n=n, d_spatial=3, kappa=kappa, seed=seed)
        return poset_from_flrw_matter_points(pts, kappa=kappa)

    return _gen


def _make_schw_gen(phi0: float):
    def _gen(n: int, seed: int):
        pts = sprinkle_schwarzschild_diamond(n=n, d_spatial=3, phi0=phi0, seed=seed)
        return poset_from_schwarzschild_points(pts, phi0=phi0)

    return _gen


def _rank_of_background(label: str, score: float, base_scores: dict[str, float]) -> int:
    all_scores = dict(base_scores)
    all_scores[label] = score
    names = [n for n, _ in sorted(all_scores.items(), key=lambda kv: kv[1])]
    return names.index(label) + 1


def run(cfg: dict[str, Any]) -> dict[str, Any]:
    out_dir = Path(cfg.get("output_dir", "outputs_carlip"))
    out_dir.mkdir(parents=True, exist_ok=True)

    n_grid: list[int] = list(cfg["n_grid"])
    reps: int = int(cfg.get("reps", 20))
    seed_base: int = int(cfg.get("seed_base", 42))

    rule = C2Rule(
        min_fail_ratio=float(cfg["hard_fail"]["min_fail_ratio"]),
        top_k_required=int(cfg["hard_fail"]["top_k_required"]),
        seed_runs=int(cfg["hard_fail"]["seed_runs"]),
    )

    windows = cfg["windows"]
    include_families = set(cfg.get("include_families", ["de_sitter", "flrw", "schwarzschild"]))
    h_weak = list(windows["de_sitter_h"])
    k_weak = list(windows["flrw_kappa"])
    p_weak = list(windows["schwarz_phi0"])

    family_param_counts = {
        "de_sitter": len(h_weak) if "de_sitter" in include_families else 0,
        "flrw": len(k_weak) if "flrw" in include_families else 0,
        "schwarzschild": len(p_weak) if "schwarzschild" in include_families else 0,
    }

    output_basename = cfg.get("output_basename", "falsify_c2_background_response")
    seedlog_path = out_dir / f"{output_basename}.seedlog.jsonl"
    snapshot_path = out_dir / f"{output_basename}.partial.json"
    checkpoint_every_seed = bool(cfg.get("checkpoint_every_seed", True))

    if checkpoint_every_seed:
        # clear previous run log
        seedlog_path.write_text("", encoding="utf-8")

    records = []

    n_seeds = rule.seed_runs
    n_cells = n_seeds * len(n_grid)
    cell_done = 0

    for s in range(n_seeds):
        for n in n_grid:
            cell_done += 1
            print(f"[F3] seed={s+1}/{n_seeds} N={n} ({cell_done}/{n_cells}) include={sorted(include_families)}", flush=True)
            mu, cov_inv, base_scores = _ref_stats(n=n, reps=reps, seed=seed_base + 1000000 * s + 10000 * n)

            if "de_sitter" in include_families:
                for h in h_weak:
                    score = _ensemble_score(_make_ds_gen(h), n, reps, mu, cov_inv, seed_base + s * 1000 + int(h * 100))
                    rank = _rank_of_background(f"dS4D_H{h}", score, base_scores)
                    records.append({"seed_run": s, "N": n, "family": "de_sitter", "param": h, "rank": rank})

            if "flrw" in include_families:
                for k in k_weak:
                    score = _ensemble_score(_make_flrw_gen(k), n, reps, mu, cov_inv, seed_base + s * 2000 + int(k * 100))
                    rank = _rank_of_background(f"FLRW_k{k}", score, base_scores)
                    records.append({"seed_run": s, "N": n, "family": "flrw", "param": k, "rank": rank})

            if "schwarzschild" in include_families:
                for p in p_weak:
                    score = _ensemble_score(_make_schw_gen(p), n, reps, mu, cov_inv, seed_base + s * 3000 + int(p * 1000))
                    rank = _rank_of_background(f"Schw_p{p}", score, base_scores)
                    records.append({"seed_run": s, "N": n, "family": "schwarzschild", "param": p, "rank": rank})

        if checkpoint_every_seed:
            event = {
                "event": "seed_complete",
                "seed_run": s,
                "completed_cells": (s + 1) * len(n_grid),
                "total_cells": n_cells,
                "records_so_far": len(records),
            }
            with seedlog_path.open("a", encoding="utf-8") as f:
                f.write(json.dumps(event, ensure_ascii=False) + "\n")

            partial = {
                "experiment_id": cfg.get("experiment_id", "falsify_c2_background_response"),
                "output_basename": output_basename,
                "include_families": sorted(include_families),
                "family_param_counts": family_param_counts,
                "rule": {
                    "min_fail_ratio": rule.min_fail_ratio,
                    "top_k_required": rule.top_k_required,
                    "seed_runs": rule.seed_runs,
                },
                "n_grid": n_grid,
                "reps": reps,
                "completed_seed_runs": s + 1,
                "total_seed_runs": n_seeds,
                "records_so_far": len(records),
                "records": records,
            }
            snapshot_path.write_text(json.dumps(partial, ensure_ascii=False, indent=2), encoding="utf-8")

    # aggregate fail ratio by (family,param)
    by_key: dict[tuple[str, float], dict[int, int]] = {}
    total_by_key: dict[tuple[str, float], dict[int, int]] = {}

    for r in records:
        key = (r["family"], float(r["param"]))
        by_key.setdefault(key, {n: 0 for n in n_grid})
        total_by_key.setdefault(key, {n: 0 for n in n_grid})
        total_by_key[key][r["N"]] += 1
        if int(r["rank"]) > rule.top_k_required:
            by_key[key][r["N"]] += 1

    hard_fail_items = []
    for key, fail_map in by_key.items():
        fail_bins = 0
        for n in n_grid:
            total = max(1, total_by_key[key][n])
            fail_rate_n = fail_map[n] / total
            if fail_rate_n >= 0.5:
                fail_bins += 1
        ratio = fail_bins / len(n_grid)
        if ratio >= rule.min_fail_ratio:
            hard_fail_items.append({
                "family": key[0],
                "param": key[1],
                "fail_bins": fail_bins,
                "total_bins": len(n_grid),
                "fail_ratio": ratio,
            })

    hard_fail = len(hard_fail_items) > 0

    result = {
        "experiment_id": cfg.get("experiment_id", "falsify_c2_background_response"),
        "include_families": sorted(include_families),
        "family_param_counts": family_param_counts,
        "rule": {
            "min_fail_ratio": rule.min_fail_ratio,
            "top_k_required": rule.top_k_required,
            "seed_runs": rule.seed_runs,
        },
        "n_grid": n_grid,
        "reps": reps,
        "records": records,
        "hard_fail": hard_fail,
        "hard_fail_items": hard_fail_items,
    }

    json_path = out_dir / f"{output_basename}.json"
    md_path = out_dir / f"{output_basename}.md"

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    lines = [
        "# F3 Falsification Report — Background Response",
        "",
        f"- Hard fail: **{'YES' if hard_fail else 'NO'}**",
        f"- Included families: {sorted(include_families)}",
        f"- Rule: fail_ratio >= {rule.min_fail_ratio}, top-k requirement={rule.top_k_required}",
        f"- N grid: {n_grid}",
        "",
        "## Hard-fail items",
    ]
    if hard_fail_items:
        for it in hard_fail_items:
            lines.append(f"- {it['family']} param={it['param']}: fail_ratio={it['fail_ratio']:.2f} ({it['fail_bins']}/{it['total_bins']})")
    else:
        lines.append("- None")

    lines += [
        "",
        "## Notes",
        "",
        "- This is a skeleton runner; metric-faithful background generators can be swapped in-place.",
        "- Current script explicitly keeps weak/moderate windows pre-registered in config.",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config_falsify_c3.yaml")
    args = parser.parse_args()

    cfg = _load_yaml(Path(args.config))
    result = run(cfg)
    output_basename = cfg.get("output_basename", "falsify_c2_background_response")
    print("HARD_FAIL" if result["hard_fail"] else "PASS")
    print(f"outputs_carlip/{output_basename}.json")


if __name__ == "__main__":
    main()
