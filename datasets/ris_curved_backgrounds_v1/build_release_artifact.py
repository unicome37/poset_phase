import hashlib
import json
import zipfile
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DIST = ROOT / "dist"
DIST.mkdir(exist_ok=True)

DATASET_NAME = "ris_curved_backgrounds_v1"
ARCHIVE_NAME = f"{DATASET_NAME}.zip"

INCLUDE_FILES = [
    "README.md",
    "CITATION.cff",
    "data/large_scale_robustness_results.json",
    "data/cost_analysis.json",
    "reports/large_scale_summary.md",
    "reports/cost_analysis.md",
    "reports/quality_check_report.md",
    "reports/quality_check_report.json",
    "code/README.md",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def main():
    now = datetime.now().isoformat()
    archive_path = DIST / ARCHIVE_NAME

    manifest = {
        "name": DATASET_NAME,
        "built_at": now,
        "files": [],
        "quality_report": "reports/quality_check_report.md",
    }

    # Validate input files
    missing = []
    for rel in INCLUDE_FILES:
        p = ROOT / rel
        if not p.exists():
            missing.append(rel)
    if missing:
        raise FileNotFoundError(f"Missing required files: {missing}")

    # Create zip archive
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for rel in INCLUDE_FILES:
            p = ROOT / rel
            zf.write(p, arcname=f"{DATASET_NAME}/{rel}")
            manifest["files"].append({
                "path": rel,
                "size_bytes": p.stat().st_size,
                "sha256": sha256_file(p)
            })

    archive_sha256 = sha256_file(archive_path)
    checksum_txt = DIST / f"{ARCHIVE_NAME}.sha256.txt"
    checksum_txt.write_text(f"{archive_sha256}  {ARCHIVE_NAME}\n", encoding="utf-8")

    manifest["archive"] = {
        "path": f"dist/{ARCHIVE_NAME}",
        "size_bytes": archive_path.stat().st_size,
        "sha256": archive_sha256,
    }

    manifest_path = DIST / "release_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)

    summary_md = DIST / "RELEASE_ARTIFACT_SUMMARY.md"
    lines = [
        "# Release Artifact Summary",
        "",
        f"- Built at: {now}",
        f"- Archive: `{ARCHIVE_NAME}`",
        f"- Archive size: {archive_path.stat().st_size} bytes",
        f"- Archive SHA256: `{archive_sha256}`",
        "",
        "## Included Files",
        "",
    ]
    for fitem in manifest["files"]:
        lines.append(f"- `{fitem['path']}` ({fitem['size_bytes']} bytes)")

    summary_md.write_text("\n".join(lines), encoding="utf-8")

    print("OK")
    print(f"Archive: {archive_path}")
    print(f"SHA256: {archive_sha256}")
    print(f"Manifest: {manifest_path}")
    print(f"Summary: {summary_md}")


if __name__ == "__main__":
    main()
