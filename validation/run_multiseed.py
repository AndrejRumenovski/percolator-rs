#!/usr/bin/env python3
"""Run matched Rust/C++ rescoring across predeclared datasets and seeds."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
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


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def version(command: list[str], environment: dict | None = None) -> str:
    result = subprocess.run(command, text=True, capture_output=True, check=False, env=environment)
    return (result.stdout + result.stderr).strip().splitlines()[0] if result.stdout or result.stderr else ""


def parse_dataset(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("dataset must be ID=PIN")
    name, path = value.split("=", 1)
    if not name or not path:
        raise argparse.ArgumentTypeError("dataset must be ID=PIN")
    return name, Path(path).resolve()


def load_q(path: Path) -> list[float]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "q-value" not in (reader.fieldnames or ()):
            raise ValueError(f"{path}: no q-value column")
        return [float(row["q-value"]) for row in reader]


def counts(path: Path) -> dict[str, int]:
    values = load_q(path)
    return {f"q_lt_{threshold:g}": sum(value < threshold for value in values) for threshold in THRESHOLDS}


def summary(values: list[float]) -> dict:
    return {
        "n": len(values),
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "sd": statistics.stdev(values) if len(values) > 1 else 0.0,
        "minimum": min(values),
        "maximum": max(values),
        "range": max(values) - min(values),
    }


def execute(command: list[str], stdout: Path, stderr: Path, environment: dict | None = None) -> dict:
    started = time.time()
    with stdout.open("w") as out, stderr.open("w") as err:
        result = subprocess.run(command, stdout=out, stderr=err, check=False, env=environment)
    return {
        "argv": command,
        "shell_display": shlex.join(command),
        "exit_code": result.returncode,
        "wall_seconds": time.time() - started,
        "stdout": str(stdout),
        "stderr": str(stderr),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", action="append", type=parse_dataset, required=True)
    parser.add_argument("--seeds", default="1,2,3,4,5")
    parser.add_argument("--rust", type=Path, default=Path("target/release/percolator-rs"))
    parser.add_argument("--cpp", type=Path, required=True)
    parser.add_argument("--cpp-library-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpp-search-input", choices=("concatenated", "separate"), default="concatenated")
    args = parser.parse_args()

    seeds = [int(value) for value in args.seeds.split(",")]
    if len(seeds) < 2 or len(set(seeds)) != len(seeds) or any(seed < 1 for seed in seeds):
        parser.error("--seeds must contain at least two unique positive integers")
    if args.output.exists():
        parser.error(f"output already exists: {args.output}")
    for name, pin in args.dataset:
        if not pin.is_file():
            parser.error(f"dataset {name}: PIN not found: {pin}")
    args.output.mkdir(parents=True)
    child_environment = os.environ.copy()
    prior_library_path = child_environment.get("LD_LIBRARY_PATH", "")
    child_environment["LD_LIBRARY_PATH"] = str(args.cpp_library_dir.resolve()) + (
        f":{prior_library_path}" if prior_library_path else ""
    )

    repository = Path(__file__).resolve().parents[1]
    commit = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=True,
    ).stdout.strip()
    dirty = bool(subprocess.run(
        ["git", "-C", str(repository), "status", "--porcelain"],
        text=True, capture_output=True, check=True,
    ).stdout)
    manifest = {
        "schema_version": 1,
        "study": "predeclared five-seed Rust/C++ PSM rescoring reproducibility",
        "created_unix": time.time(),
        "threshold_policy": "q < threshold",
        "thresholds": THRESHOLDS,
        "seeds": seeds,
        "environment": {
            "platform": platform.platform(),
            "python": sys.version,
            "rustc": version(["rustc", "--version"]),
            "percolator_rs_commit": commit,
            "percolator_rs_worktree_dirty": dirty,
            "reference_percolator": version([str(args.cpp), "--help"], child_environment),
            "reference_library_dir": str(args.cpp_library_dir.resolve()),
            "cpp_search_input": args.cpp_search_input,
            "cwd": str(repository),
        },
        "evaluation_scripts": [
            {"path": str(Path(__file__).resolve()), "sha256": digest(Path(__file__).resolve())},
            {
                "path": str(Path(psm_agreement.__file__).resolve()),
                "sha256": digest(Path(psm_agreement.__file__).resolve()),
            },
        ],
        "datasets": [
            {"id": name, "pin": str(pin), "sha256": digest(pin), "bytes": pin.stat().st_size}
            for name, pin in args.dataset
        ],
        "runs": [],
    }
    (args.output / "manifest.partial.json").write_text(json.dumps(manifest, indent=2) + "\n")

    for name, pin in args.dataset:
        for seed in seeds:
            run_root = args.output / name / f"seed-{seed}"
            implementations = {}
            for implementation, binary in (("rust", args.rust.resolve()), ("cpp", args.cpp.resolve())):
                destination = run_root / implementation
                destination.mkdir(parents=True)
                command = [str(binary), "--seed", str(seed), "--num-threads", "1"]
                if implementation == "rust":
                    command.extend(["--canonical", "--no-select-c"])
                else:
                    command.extend(["--search-input", args.cpp_search_input])
                command.extend([
                    "--results-psms", str(destination / "target.psms.tsv"),
                    "--decoy-results-psms", str(destination / "decoy.psms.tsv"),
                    "--results-peptides", str(destination / "target.peptides.tsv"),
                    "--decoy-results-peptides", str(destination / "decoy.peptides.tsv"),
                    str(pin),
                ])
                execution = execute(
                    command, destination / "stdout.log", destination / "stderr.log", child_environment
                )
                if execution["exit_code"] != 0:
                    manifest["runs"].append({
                        "dataset": name, "seed": seed, "implementation": implementation,
                        "execution": execution, "status": "failed",
                    })
                    (args.output / "manifest.partial.json").write_text(json.dumps(manifest, indent=2) + "\n")
                    raise RuntimeError(f"{name} seed {seed} {implementation} failed")
                artifacts = {}
                for artifact in destination.glob("*.tsv"):
                    artifacts[artifact.name] = {
                        "path": str(artifact), "sha256": digest(artifact), "bytes": artifact.stat().st_size,
                    }
                implementations[implementation] = {
                    "execution": execution,
                    "psms": counts(destination / "target.psms.tsv"),
                    "peptides": counts(destination / "target.peptides.tsv"),
                    "artifacts": artifacts,
                }

            agreement = psm_agreement.compare(run_root / "rust", run_root / "cpp")
            agreement_path = run_root / "agreement.json"
            agreement_path.write_text(json.dumps(agreement, indent=2, sort_keys=True) + "\n")
            manifest["runs"].append({
                "dataset": name,
                "seed": seed,
                "status": "complete",
                "implementations": implementations,
                "agreement": str(agreement_path),
                "agreement_sha256": digest(agreement_path),
            })
            (args.output / "manifest.partial.json").write_text(json.dumps(manifest, indent=2) + "\n")
            q01 = next(row for row in agreement["thresholds"] if row["q_threshold"] == 0.01)
            print(f"{name} seed={seed}: Rust/C++ q<0.01={q01['rust']}/{q01['cpp']}, J={q01['jaccard']:.4f}", flush=True)

    aggregate = {}
    for name, _ in args.dataset:
        dataset_runs = [run for run in manifest["runs"] if run["dataset"] == name and run["status"] == "complete"]
        dataset_summary = {}
        for implementation in ("rust", "cpp"):
            dataset_summary[implementation] = {}
            for level in ("psms", "peptides"):
                for threshold in THRESHOLDS:
                    key = f"q_lt_{threshold:g}"
                    values = [run["implementations"][implementation][level][key] for run in dataset_runs]
                    dataset_summary[implementation][f"{level}_{key}"] = summary(values)
        agreement_rows = [json.loads(Path(run["agreement"]).read_text()) for run in dataset_runs]
        dataset_summary["agreement"] = {}
        for threshold in THRESHOLDS:
            rows = [next(row for row in item["thresholds"] if row["q_threshold"] == threshold) for item in agreement_rows]
            for metric in ("intersection", "rust_only", "cpp_only", "jaccard"):
                dataset_summary["agreement"][f"{metric}_q_lt_{threshold:g}"] = summary([row[metric] for row in rows])
        for metric in ("score_spearman", "q_value_spearman", "pep_spearman"):
            values = [item["correlations_on_matching_psms"][metric] for item in agreement_rows]
            dataset_summary["agreement"][metric] = summary([value for value in values if value is not None])
        aggregate[name] = dataset_summary

    manifest["aggregate"] = aggregate
    manifest["completed_unix"] = time.time()
    final = args.output / "manifest.json"
    final.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (args.output / "manifest.partial.json").unlink()
    print(f"machine-readable manifest: {final}")


if __name__ == "__main__":
    main()
