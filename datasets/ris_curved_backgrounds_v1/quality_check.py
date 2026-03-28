import json
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
REPORTS = ROOT / "reports"


def check_exists(p: Path):
    return p.exists() and p.is_file()


def load_json(p: Path):
    with p.open("r", encoding="utf-8") as f:
        return json.load(f)


def main():
    findings = {
        "timestamp": datetime.now().isoformat(),
        "dataset_root": str(ROOT),
        "checks": {},
        "summary": {"passed": 0, "failed": 0},
        "warnings": []
    }

    required_files = [
        ROOT / "README.md",
        ROOT / "CITATION.cff",
        DATA / "large_scale_robustness_results.json",
        DATA / "cost_analysis.json",
        REPORTS / "large_scale_summary.md",
        REPORTS / "cost_analysis.md",
        ROOT / "code" / "README.md",
    ]

    for fp in required_files:
        ok = check_exists(fp)
        findings["checks"][f"exists::{fp.relative_to(ROOT)}"] = ok
        findings["summary"]["passed" if ok else "failed"] += 1

    # JSON parse checks
    robust = None
    cost = None
    for name in ["large_scale_robustness_results.json", "cost_analysis.json"]:
        target = DATA / name
        key = f"json_parse::{name}"
        try:
            obj = load_json(target)
            findings["checks"][key] = True
            findings["summary"]["passed"] += 1
            if name == "large_scale_robustness_results.json":
                robust = obj
            else:
                cost = obj
        except Exception as e:
            findings["checks"][key] = False
            findings["summary"]["failed"] += 1
            findings["warnings"].append(f"Failed parsing {name}: {e}")

    # Semantic checks for robustness dataset
    if robust is not None:
        expected_partitions = {"N512", "N1024", "N2048"}
        got_partitions = set(robust.keys())
        ok_partitions = expected_partitions.issubset(got_partitions)
        findings["checks"]["robustness::partitions_512_1024_2048"] = ok_partitions
        findings["summary"]["passed" if ok_partitions else "failed"] += 1

        total_trials = 0
        all_success_100 = True
        for p in expected_partitions:
            cfgs = robust.get(p, {}).get("configurations", {})
            total_trials += sum(int(v.get("num_trials", 0)) for v in cfgs.values())
            for _, v in cfgs.items():
                if float(v.get("success_rate", -1)) != 1.0:
                    all_success_100 = False

        findings["checks"]["robustness::total_trials_is_180"] = (total_trials == 180)
        findings["summary"]["passed" if total_trials == 180 else "failed"] += 1

        findings["checks"]["robustness::all_success_rate_100"] = all_success_100
        findings["summary"]["passed" if all_success_100 else "failed"] += 1

    # Semantic checks for cost dataset
    if cost is not None:
        sizes = cost.get("sizes", [])
        feas = cost.get("feasibility", [])
        mem = cost.get("memory_mb", [])

        ok_sizes = sizes == [512, 768, 1024, 1536, 2048]
        findings["checks"]["cost::sizes_expected"] = ok_sizes
        findings["summary"]["passed" if ok_sizes else "failed"] += 1

        ok_feas = all(bool(x) for x in feas) and len(feas) == len(sizes)
        findings["checks"]["cost::all_feasible"] = ok_feas
        findings["summary"]["passed" if ok_feas else "failed"] += 1

        ok_mem_max = (len(mem) > 0 and max(mem) <= 41.6 + 1e-9)
        findings["checks"]["cost::max_memory_41_6_mb"] = ok_mem_max
        findings["summary"]["passed" if ok_mem_max else "failed"] += 1

    # CITATION sanity (lightweight)
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8") if (ROOT / "CITATION.cff").exists() else ""
    cff_checks = {
        "citation::has_cff_version": "cff-version:" in cff,
        "citation::has_title": "title:" in cff,
        "citation::has_type_dataset": "type: dataset" in cff,
        "citation::has_repository_code": "repository-code:" in cff,
    }
    for k, ok in cff_checks.items():
        findings["checks"][k] = ok
        findings["summary"]["passed" if ok else "failed"] += 1

    out_json = REPORTS / "quality_check_report.json"
    out_md = REPORTS / "quality_check_report.md"

    with out_json.open("w", encoding="utf-8") as f:
        json.dump(findings, f, indent=2, ensure_ascii=False)

    status = "PASS" if findings["summary"]["failed"] == 0 else "FAIL"
    lines = [
        "# Dataset Quality Check Report",
        "",
        f"- Time: {findings['timestamp']}",
        f"- Status: **{status}**",
        f"- Passed: {findings['summary']['passed']}",
        f"- Failed: {findings['summary']['failed']}",
        "",
        "## Detailed Checks",
        "",
    ]
    for k, v in findings["checks"].items():
        lines.append(f"- {'✅' if v else '❌'} `{k}`")

    if findings["warnings"]:
        lines += ["", "## Warnings", ""]
        lines += [f"- {w}" for w in findings["warnings"]]

    out_md.write_text("\n".join(lines), encoding="utf-8")
    print(status)
    print(f"JSON report: {out_json}")
    print(f"MD report: {out_md}")


if __name__ == "__main__":
    main()
