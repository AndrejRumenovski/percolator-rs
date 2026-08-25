#!/usr/bin/env python3
"""Repeated exchangeable-label pure-null FDR experiment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import random
import shlex
import statistics
import subprocess
import sys
import time
from pathlib import Path


THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.10)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_input(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("input must be ID=PIN")
    name, path = value.split("=", 1)
    return name, Path(path).resolve()


def make_null(source: Path, destination: Path, seed: int) -> dict:
    with source.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        label_column = next((i for i, name in enumerate(header) if name.lower() == "label"), None)
        if label_column is None:
            raise ValueError(f"{source}: no Label column")
        rows = []
        for row in reader:
            if row and row[0] == "DefaultDirection":
                continue
            if len(row) <= label_column:
                raise ValueError(f"{source}: short row in input")
            if int(row[label_column]) < 0:
                rows.append(row)
    if len(rows) < 100:
        raise ValueError(f"{source}: fewer than 100 source decoy rows")

    order = list(range(len(rows)))
    random.Random(seed).shuffle(order)
    target_indices = set(order[: len(rows) // 2])
    for index, row in enumerate(rows):
        row[label_column] = "1" if index in target_indices else "-1"
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return {
        "source_decoy_rows": len(rows),
        "pseudo_targets": len(target_indices),
        "pseudo_decoys": len(rows) - len(target_indices),
        "relabel_seed": seed,
        "sha256": sha256(destination),
    }


def q_counts(path: Path) -> dict:
    values = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            q = float(row["q-value"])
            if not math.isfinite(q) or not 0 <= q <= 1:
                raise ValueError(f"{path}: invalid q-value {q}")
            values.append(q)
    return {f"q_lt_{threshold:g}": sum(q < threshold for q in values) for threshold in THRESHOLDS}


def run(command: list[str], stdout: Path, stderr: Path, environment: dict) -> dict:
    started = time.time()
    with stdout.open("w") as out, stderr.open("w") as err:
        result = subprocess.run(command, stdout=out, stderr=err, env=environment, check=False)
    return {
        "argv": command,
        "shell_display": shlex.join(command),
        "exit_code": result.returncode,
        "wall_seconds": time.time() - started,
        "stdout": str(stdout),
        "stderr": str(stderr),
    }


def wilson(successes: int, trials: int) -> tuple[float, float]:
    if trials == 0:
        return 0.0, 0.0
    z = 1.959963984540054
    estimate = successes / trials
    denominator = 1 + z * z / trials
    center = (estimate + z * z / (2 * trials)) / denominator
    half = z * math.sqrt(estimate * (1 - estimate) / trials + z * z / (4 * trials * trials)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def describe(values: list[int]) -> dict:
    return {
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "sd": statistics.stdev(values) if len(values) > 1 else 0.0,
        "minimum": min(values),
        "maximum": max(values),
        "range": max(values) - min(values),
    }


def first_line(command: list[str], environment: dict) -> str:
    result = subprocess.run(command, text=True, capture_output=True, env=environment, check=False)
    lines = (result.stdout + result.stderr).splitlines()
    return lines[0] if lines else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", action="append", type=parse_input, required=True)
    parser.add_argument("--relabel-seeds", default="1001,1002,1003,1004,1005,1006,1007,1008,1009,1010")
    parser.add_argument("--model-seed", type=int, default=1)
    parser.add_argument("--rust", type=Path, required=True)
    parser.add_argument("--cpp", type=Path)
    parser.add_argument("--cpp-library-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    seeds = [int(value) for value in args.relabel_seeds.split(",")]
    if len(seeds) < 2 or len(seeds) != len(set(seeds)):
        parser.error("relabel seeds must contain at least two unique integers")
    if bool(args.cpp) != bool(args.cpp_library_dir):
        parser.error("--cpp and --cpp-library-dir must be supplied together")
    if args.output.exists():
        parser.error(f"output already exists: {args.output}")
    for name, path in args.input:
        if not path.is_file():
            parser.error(f"input {name} does not exist: {path}")
    args.output.mkdir(parents=True)

    environment = os.environ.copy()
    if args.cpp_library_dir:
        old = environment.get("LD_LIBRARY_PATH", "")
        environment["LD_LIBRARY_PATH"] = str(args.cpp_library_dir.resolve()) + (f":{old}" if old else "")
    repo = Path(__file__).resolve().parents[1]
    scripts = [Path(__file__).resolve()]
    manifest = {
        "schema_version": 1,
        "study": "exchangeable-label pure-null FDR validation",
        "null_interpretation": (
            "Every pseudo-target is false because all retained rows were decoy search matches and "
            "target/decoy labels were reassigned independently of features with exact class balance. "
            "Under the complete null FDP is 1 when any target is accepted and 0 otherwise, so empirical "
            "FDR across relabel replicates is the probability of one or more acceptances."
        ),
        "threshold_policy": "q < threshold",
        "thresholds": THRESHOLDS,
        "relabel_seeds": seeds,
        "model_seed": args.model_seed,
        "environment": {
            "platform": platform.platform(),
            "python": sys.version,
            "commit": subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True, capture_output=True, check=True).stdout.strip(),
            "worktree_dirty": bool(subprocess.run(["git", "-C", str(repo), "status", "--porcelain"], text=True, capture_output=True, check=True).stdout),
            "rust_version": first_line([str(args.rust.resolve()), "--help"], environment),
            "cpp_version": first_line([str(args.cpp.resolve()), "--help"], environment) if args.cpp else None,
        },
        "evaluation_scripts": [{"path": str(path), "sha256": sha256(path)} for path in scripts],
        "inputs": [{"id": name, "path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in args.input],
        "runs": [],
    }
    (args.output / "manifest.partial.json").write_text(json.dumps(manifest, indent=2) + "\n")

    methods = [("rust", args.rust.resolve())]
    if args.cpp:
        methods.append(("cpp", args.cpp.resolve()))
    for name, source in args.input:
        for seed in seeds:
            root = args.output / name / f"relabel-{seed}"
            root.mkdir(parents=True)
            null_pin = root / "null.pin"
            construction = make_null(source, null_pin, seed)
            for method, binary in methods:
                destination = root / method
                destination.mkdir()
                command = [str(binary), "--seed", str(args.model_seed), "--num-threads", "1"]
                if method == "rust":
                    command.extend(["--canonical", "--no-select-c"])
                else:
                    command.extend(["--search-input", "concatenated"])
                command.extend([
                    "--results-psms", str(destination / "target.psms.tsv"),
                    "--decoy-results-psms", str(destination / "decoy.psms.tsv"),
                    str(null_pin),
                ])
                execution = run(command, destination / "stdout.log", destination / "stderr.log", environment)
                record = {
                    "input": name, "relabel_seed": seed, "model_seed": args.model_seed,
                    "method": method, "construction": construction, "execution": execution,
                }
                if execution["exit_code"] == 0:
                    record["false_target_psms"] = q_counts(destination / "target.psms.tsv")
                    record["target_output_sha256"] = sha256(destination / "target.psms.tsv")
                    record["status"] = "complete"
                else:
                    record["status"] = "failed"
                manifest["runs"].append(record)
                (args.output / "manifest.partial.json").write_text(json.dumps(manifest, indent=2) + "\n")
                if record["status"] == "complete":
                    print(f"{name} relabel={seed} {method}: q<0.01 false targets={record['false_target_psms']['q_lt_0.01']}", flush=True)
                else:
                    print(f"{name} relabel={seed} {method}: FAILED exit={execution['exit_code']}", flush=True)

    aggregate = []
    for name, _ in args.input:
        for method, _ in methods:
            attempted_runs = [run for run in manifest["runs"] if run["input"] == name and run["method"] == method]
            selected_runs = [run for run in attempted_runs if run["status"] == "complete"]
            for threshold in THRESHOLDS:
                key = f"q_lt_{threshold:g}"
                values = [run["false_target_psms"][key] for run in selected_runs]
                positives = sum(value > 0 for value in values)
                low, high = wilson(positives, len(values))
                aggregate.append({
                    "input": name,
                    "method": method,
                    "q_threshold": threshold,
                    "replicates": len(values),
                    "attempted_replicates": len(attempted_runs),
                    "failed_replicates": len(attempted_runs) - len(values),
                    "replicates_with_any_false_discovery": positives,
                    "empirical_fdr": positives / len(values),
                    "empirical_fdr_wilson95_low": low,
                    "empirical_fdr_wilson95_high": high,
                    "nominal_fdr": threshold,
                    "false_target_count": describe(values),
                })
    manifest["aggregate"] = aggregate
    manifest["completed_unix"] = time.time()
    final = args.output / "manifest.json"
    final.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    with (args.output / "calibration.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=(
            "input", "method", "q_threshold", "replicates", "attempted_replicates", "failed_replicates", "replicates_with_any_false_discovery",
            "empirical_fdr", "empirical_fdr_wilson95_low", "empirical_fdr_wilson95_high", "nominal_fdr",
            "false_target_mean", "false_target_median", "false_target_sd", "false_target_minimum",
            "false_target_maximum", "false_target_range",
        ))
        writer.writeheader()
        for row in aggregate:
            counts = row["false_target_count"]
            writer.writerow({
                **{key: value for key, value in row.items() if key != "false_target_count"},
                **{f"false_target_{key}": value for key, value in counts.items()},
            })
    (args.output / "manifest.partial.json").unlink()
    print(f"machine-readable manifest: {final}")


if __name__ == "__main__":
    main()
