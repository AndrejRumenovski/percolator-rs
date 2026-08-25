#!/usr/bin/env python3
"""Five-seed matched Rust/C++ signal-present entrapment validation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shlex
import statistics
import subprocess
import sys
import time
from pathlib import Path

import psm_agreement


THRESHOLDS = psm_agreement.THRESHOLDS


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def execute(command: list[str], stdout: Path, stderr: Path, environment: dict) -> dict:
    started = time.time()
    with stdout.open("w") as out, stderr.open("w") as err:
        result = subprocess.run(command, stdout=out, stderr=err, env=environment, check=False)
    return {
        "argv": command, "shell_display": shlex.join(command), "exit_code": result.returncode,
        "wall_seconds": time.time() - started, "stdout": str(stdout), "stderr": str(stderr),
    }


def members(row: dict, decoy: bool) -> list[str]:
    values = []
    for key, value in row.items():
        if key == "proteinIds" or key is None:
            values.extend(value if isinstance(value, list) else [value])
    proteins = [protein for value in values if value for protein in value.replace(";", "\t").split("\t") if protein]
    if decoy:
        proteins = [protein.removeprefix("DECOY_").removeprefix("decoy_") for protein in proteins]
    return proteins


def load(path: Path, decoy: bool) -> list[tuple[float, bool, bool]]:
    rows = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            q = float(row["q-value"])
            proteins = members(row, decoy)
            pure = bool(proteins) and all(protein.startswith("ENT_") for protein in proteins)
            mixed = any(protein.startswith("ENT_") for protein in proteins) and not pure
            rows.append((q, pure, mixed))
    return rows


def calibrate(targets: list, decoys: list, fallback_fraction: float) -> list[dict]:
    output = []
    for threshold in THRESHOLDS:
        accepted = [row for row in targets if row[0] < threshold]
        accepted_decoys = [row for row in decoys if row[0] < threshold and not row[2]]
        fraction = (
            sum(row[1] for row in accepted_decoys) / len(accepted_decoys)
            if accepted_decoys else fallback_fraction
        )
        entrapment = sum(row[1] for row in accepted)
        estimated_false = entrapment / fraction
        output.append({
            "q_threshold": threshold, "comparison": "strict_less_than",
            "accepted_targets": len(accepted), "pure_entrapment": entrapment,
            "mixed": sum(row[2] for row in accepted),
            "effective_entrapment_fraction": fraction,
            "estimated_false": estimated_false,
            "adjusted_fdp": estimated_false / len(accepted) if accepted else 0.0,
        })
    return output


def describe(values: list[float]) -> dict:
    return {
        "n": len(values), "mean": statistics.fmean(values), "median": statistics.median(values),
        "sd": statistics.stdev(values) if len(values) > 1 else 0.0,
        "minimum": min(values), "maximum": max(values), "range": max(values) - min(values),
    }


def first_line(command: list[str], environment: dict) -> str:
    result = subprocess.run(command, text=True, capture_output=True, env=environment, check=False)
    lines = (result.stdout + result.stderr).splitlines()
    return lines[0] if lines else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pins", type=Path, required=True, help="directory containing comet-*/comet.pin")
    parser.add_argument("--seeds", default="1,2,3,4,5")
    parser.add_argument("--entrapment-fraction", type=float, required=True)
    parser.add_argument("--rust", type=Path, required=True)
    parser.add_argument("--cpp", type=Path, required=True)
    parser.add_argument("--cpp-library-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    seeds = [int(value) for value in args.seeds.split(",")]
    pins = sorted(
        pin for pin in args.pins.glob("comet-*/comet.pin")
        if pin.parent.name != "comet-out"
    )
    if len(pins) != 6:
        parser.error(f"expected exactly six entrapment PINs, found {len(pins)}")
    if args.output.exists():
        parser.error(f"output already exists: {args.output}")
    args.output.mkdir(parents=True)
    environment = os.environ.copy()
    old = environment.get("LD_LIBRARY_PATH", "")
    environment["LD_LIBRARY_PATH"] = str(args.cpp_library_dir.resolve()) + (f":{old}" if old else "")
    repo = Path(__file__).resolve().parents[1]
    manifest = {
        "schema_version": 1,
        "study": "matched five-seed signal-present foreign-proteome entrapment",
        "threshold_policy": "q < threshold",
        "thresholds": THRESHOLDS,
        "seeds": seeds,
        "entrapment_fraction_fallback": args.entrapment_fraction,
        "environment": {
            "platform": platform.platform(), "python": sys.version,
            "commit": subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True, capture_output=True, check=True).stdout.strip(),
            "worktree_dirty": bool(subprocess.run(["git", "-C", str(repo), "status", "--porcelain"], text=True, capture_output=True, check=True).stdout),
            "rust": first_line([str(args.rust.resolve()), "--help"], environment),
            "cpp": first_line([str(args.cpp.resolve()), "--help"], environment),
            "cpp_search_input": "concatenated",
        },
        "evaluation_scripts": [
            {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
            {"path": str((Path(__file__).parent / "psm_agreement.py").resolve()), "sha256": sha256((Path(__file__).parent / "psm_agreement.py").resolve())},
        ],
        "inputs": [{"path": str(pin), "sha256": sha256(pin), "bytes": pin.stat().st_size} for pin in pins],
        "runs": [],
    }
    (args.output / "manifest.partial.json").write_text(json.dumps(manifest, indent=2) + "\n")

    for seed in seeds:
        seed_root = args.output / f"seed-{seed}"
        method_rows = {}
        for method, binary in (("rust", args.rust.resolve()), ("cpp", args.cpp.resolve())):
            target_rows, decoy_rows, executions = [], [], []
            for pin in pins:
                stem = pin.parent.name.removeprefix("comet-")
                destination = seed_root / method / stem
                destination.mkdir(parents=True)
                command = [str(binary), "--seed", str(seed), "--num-threads", "1"]
                if method == "rust":
                    command.extend(["--canonical", "--no-select-c"])
                else:
                    command.extend(["--search-input", "concatenated"])
                command.extend([
                    "--results-psms", str(destination / "target.psms.tsv"),
                    "--decoy-results-psms", str(destination / "decoy.psms.tsv"), str(pin),
                ])
                execution = execute(command, destination / "stdout.log", destination / "stderr.log", environment)
                execution["input"] = str(pin)
                executions.append(execution)
                if execution["exit_code"] != 0:
                    raise RuntimeError(f"seed {seed} {method} failed on {pin}")
                target_rows.extend(load(destination / "target.psms.tsv", False))
                decoy_rows.extend(load(destination / "decoy.psms.tsv", True))
            method_rows[method] = calibrate(target_rows, decoy_rows, args.entrapment_fraction)
            manifest["runs"].append({
                "seed": seed, "method": method, "executions": executions,
                "calibration": method_rows[method], "status": "complete",
            })
            row = next(row for row in method_rows[method] if row["q_threshold"] == 0.01)
            print(f"seed={seed} {method}: accepted={row['accepted_targets']} adjusted_FDP={row['adjusted_fdp']:.4%}", flush=True)
        agreement = psm_agreement.compare(seed_root / "rust", seed_root / "cpp")
        agreement_path = seed_root / "agreement.json"
        agreement_path.write_text(json.dumps(agreement, indent=2, sort_keys=True) + "\n")
        manifest["runs"].append({"seed": seed, "method": "agreement", "path": str(agreement_path), "sha256": sha256(agreement_path), "status": "complete"})
        (args.output / "manifest.partial.json").write_text(json.dumps(manifest, indent=2) + "\n")

    aggregate = []
    for method in ("rust", "cpp"):
        method_runs = [run for run in manifest["runs"] if run["method"] == method]
        for threshold in THRESHOLDS:
            rows = [next(row for row in run["calibration"] if row["q_threshold"] == threshold) for run in method_runs]
            aggregate.append({
                "method": method, "q_threshold": threshold,
                **{metric: describe([row[metric] for row in rows]) for metric in (
                    "accepted_targets", "pure_entrapment", "effective_entrapment_fraction", "estimated_false", "adjusted_fdp"
                )},
            })
    agreement_runs = [json.loads(Path(run["path"]).read_text()) for run in manifest["runs"] if run["method"] == "agreement"]
    agreement_aggregate = []
    for threshold in THRESHOLDS:
        rows = [next(row for row in run["thresholds"] if row["q_threshold"] == threshold) for run in agreement_runs]
        agreement_aggregate.append({
            "q_threshold": threshold,
            **{metric: describe([row[metric] for row in rows]) for metric in ("rust", "cpp", "intersection", "rust_only", "cpp_only", "jaccard")},
        })
    manifest["aggregate"] = aggregate
    manifest["agreement_aggregate"] = agreement_aggregate
    manifest["completed_unix"] = time.time()
    final = args.output / "manifest.json"
    final.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (args.output / "manifest.partial.json").unlink()
    print(f"machine-readable manifest: {final}")


if __name__ == "__main__":
    main()
