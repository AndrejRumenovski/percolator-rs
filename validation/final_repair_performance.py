#!/usr/bin/env python3
"""Check end-to-end output equivalence across build and parallel paths."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(binary: Path, pin: Path, output: Path, threads: int) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    target = output / "target.tsv"
    decoy = output / "decoy.tsv"
    command = [
        str(binary), "--canonical", "--no-select-c", "--seed", "1",
        "--num-threads", str(threads), "--results-psms", str(target),
        "--decoy-results-psms", str(decoy), str(pin),
    ]
    before = time.perf_counter()
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    elapsed = time.perf_counter() - before
    (output / "stdout.log").write_text(result.stdout)
    (output / "stderr.log").write_text(result.stderr)
    if result.returncode:
        raise RuntimeError(f"failed ({result.returncode}): {' '.join(command)}")
    return {
        "binary": str(binary),
        "binary_sha256": digest(binary),
        "threads": threads,
        "wall_seconds": elapsed,
        "target_sha256": digest(target),
        "decoy_sha256": digest(decoy),
        "target_bytes": target.stat().st_size,
        "decoy_bytes": decoy.stat().st_size,
        "command": command,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", required=True, type=Path)
    parser.add_argument("--debug", required=True, type=Path)
    parser.add_argument("--pin", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    variants = [
        ("release-serial", args.release, 1),
        ("release-parallel", args.release, 3),
        ("debug-serial", args.debug, 1),
    ]
    runs = {
        name: run(binary.resolve(), args.pin.resolve(), args.out / name, threads)
        for name, binary, threads in variants
    }
    signatures = {(item["target_sha256"], item["decoy_sha256"]) for item in runs.values()}
    payload = {
        "input": str(args.pin.resolve()),
        "input_sha256": digest(args.pin),
        "output_byte_equivalent": len(signatures) == 1,
        "runs": runs,
    }
    (args.out / "manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
