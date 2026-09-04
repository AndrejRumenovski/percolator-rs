#!/usr/bin/env python3
"""Freeze a complete integrity inventory for the homology experiment."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


EXCLUDED_NAMES = {"ARTIFACT_MANIFEST.json", "SHA256SUMS.txt"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def command(args: list[str], cwd: Path) -> str:
    return subprocess.check_output(args, cwd=cwd, text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-root", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment_root.resolve()
    repo = args.repo.resolve()

    files = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name in EXCLUDED_NAMES or "__pycache__" in path.parts:
            continue
        files.append({
            "path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })

    production_paths = ["src", "Cargo.toml", "Cargo.lock", "build.rs"]
    diff = command(["git", "diff", "--", *production_paths], repo)
    manifest = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "experiment_root": str(root),
        "repository": {
            "path": str(repo),
            "commit": command(["git", "rev-parse", "HEAD"], repo),
            "production_paths_checked": production_paths,
            "production_source_diff_empty": diff == "",
        },
        "audited_binary": {
            "path": str(args.binary.resolve()),
            "bytes": args.binary.stat().st_size,
            "sha256": sha256(args.binary),
        },
        "platform": platform.platform(),
        "inventory_excludes": sorted(EXCLUDED_NAMES) + ["__pycache__/"],
        "files": files,
    }
    manifest_path = root / "ARTIFACT_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    sums = root / "SHA256SUMS.txt"
    with sums.open("w") as handle:
        for row in files:
            handle.write(f"{row['sha256']}  {row['path']}\n")
        handle.write(f"{sha256(manifest_path)}  ARTIFACT_MANIFEST.json\n")

    print(json.dumps({
        "files": len(files),
        "bytes": sum(row["bytes"] for row in files),
        "production_source_diff_empty": manifest["repository"]["production_source_diff_empty"],
        "binary_sha256": manifest["audited_binary"]["sha256"],
        "manifest_sha256": sha256(manifest_path),
    }, indent=2))


if __name__ == "__main__":
    main()
