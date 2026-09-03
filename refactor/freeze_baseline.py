#!/usr/bin/env python3
"""Freeze a behavior/performance baseline before the architecture refactor.

The harness records exact commands, stdout/stderr, canonical TSVs, adversarial
observations, and optional repeated PXD032157 timings.  It never modifies
production source and refuses to overwrite an artifact directory.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "target/release/percolator-rs"
FIXTURE = ROOT / "tests/fixtures/sample.pin"
PRODUCTION_PATHS = ("src", "tests", "Cargo.toml", "Cargo.lock")
SHELL_GATES = (
    "tests/regression.sh",
    "tests/ensemble_regression.sh",
    "tests/selection_regression.sh",
    "tests/model_regression.sh",
    "tests/protein_regression.sh",
    "tests/feature_report.sh",
)
OUTPUT_NAMES = (
    "target.psms.tsv",
    "decoy.psms.tsv",
    "target.peptides.tsv",
    "decoy.peptides.tsv",
)
Q01_PSM = re.compile(r"target PSMs q<0\.01: ([0-9]+)")
Q01_PEPTIDE = re.compile(r"target peptides q<0\.01: ([0-9]+)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, Any]:
    return {"bytes": path.stat().st_size, "sha256": sha256(path)}


def command_text(command: list[str]) -> str:
    return " ".join(command)


def execute(
    command: list[str],
    record_dir: Path,
    name: str,
    *,
    check: bool = True,
    cwd: Path = ROOT,
) -> dict[str, Any]:
    record_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    wall = time.perf_counter() - started
    stdout = record_dir / f"{name}.stdout.txt"
    stderr = record_dir / f"{name}.stderr.txt"
    stdout.write_text(result.stdout)
    stderr.write_text(result.stderr)
    record = {
        "argv": command,
        "command": command_text(command),
        "exit_code": result.returncode,
        "wall_seconds": wall,
        "stdout": str(stdout.relative_to(record_dir.parents[1])),
        "stdout_sha256": sha256(stdout),
        "stderr": str(stderr.relative_to(record_dir.parents[1])),
        "stderr_sha256": sha256(stderr),
    }
    if check and result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {command_text(command)}\n{result.stderr}"
        )
    return record


def git_text(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, text=True, capture_output=True, check=True
    ).stdout


def assert_production_clean() -> None:
    changed = git_text("status", "--porcelain=v1", "--", *PRODUCTION_PATHS)
    if changed:
        raise RuntimeError(
            "baseline requires clean production/build/test paths; found:\n" + changed
        )


def output_command(
    destination: Path,
    extra: list[str],
    inputs: list[str],
    *,
    proteins: bool,
) -> list[str]:
    destination.mkdir(parents=True, exist_ok=True)
    command = [
        str(BINARY),
        "--canonical",
        "--seed",
        "1",
        *extra,
        "--results-psms",
        str(destination / "target.psms.tsv"),
        "--decoy-results-psms",
        str(destination / "decoy.psms.tsv"),
        "--results-peptides",
        str(destination / "target.peptides.tsv"),
        "--decoy-results-peptides",
        str(destination / "decoy.peptides.tsv"),
    ]
    if proteins:
        command.extend(
            [
                "--results-proteins",
                str(destination / "target.proteins.tsv"),
                "--decoy-results-proteins",
                str(destination / "decoy.proteins.tsv"),
            ]
        )
    command.extend(inputs)
    return command


def capture_outputs(root: Path, commands: Path) -> dict[str, Any]:
    variants = {
        "fixed-serial": (["--no-select-c", "--num-threads", "1"], [str(FIXTURE)], True),
        "fixed-parallel": (["--no-select-c", "--num-threads", "3"], [str(FIXTURE)], True),
        "select-c": (["--select-c", "--num-threads", "1"], [str(FIXTURE)], True),
        "ensemble": (
            ["--no-select-c", "--num-threads", "1", "--ensemble"],
            [f"comet={FIXTURE}", f"tide={FIXTURE}"],
            False,
        ),
    }
    output: dict[str, Any] = {}
    for name, (extra, inputs, proteins) in variants.items():
        destination = root / name
        execution = execute(
            output_command(destination, extra, inputs, proteins=proteins),
            commands,
            f"output-{name}",
        )
        files = {
            path.name: file_record(path)
            for path in sorted(destination.iterdir())
            if path.is_file()
        }
        output[name] = {"execution": execution, "files": files}

    fixed = output["fixed-serial"]["files"]
    parallel = output["fixed-parallel"]["files"]
    output["fixed_serial_parallel_byte_identical"] = all(
        fixed[name]["sha256"] == parallel[name]["sha256"]
        for name in fixed.keys() & parallel.keys()
    )

    with tempfile.TemporaryDirectory(prefix="percolator-baseline-repeat-") as temporary:
        repeat = Path(temporary)
        execute(
            output_command(
                repeat,
                ["--no-select-c", "--num-threads", "1"],
                [str(FIXTURE)],
                proteins=True,
            ),
            commands,
            "output-fixed-repeat",
        )
        output["fixed_repeat_byte_identical"] = all(
            sha256(repeat / name) == fixed[name]["sha256"] for name in fixed
        )
    return output


def compile_and_run_probe(source: str, commands: Path) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="percolator-baseline-probe-") as temporary:
        binary = Path(temporary) / source.removesuffix(".rs")
        compile_record = execute(
            ["rustc", "--edition=2021", f"validation/{source}", "-O", "-o", str(binary)],
            commands,
            f"compile-{source.removesuffix('.rs')}",
        )
        run_record = execute(
            [str(binary)], commands, f"run-{source.removesuffix('.rs')}"
        )
    return {"compile": compile_record, "run": run_record}


def capture_adversarial(artifacts: Path, commands: Path) -> dict[str, Any]:
    root = artifacts / "adversarial"
    result_json = root / "final-repair-adversarial.json"
    record = execute(
        [
            "python3",
            "validation/final_repair_adversarial.py",
            "--binary",
            str(BINARY),
            "--work-dir",
            str(root / "work"),
            "--json",
            str(result_json),
        ],
        commands,
        "final-repair-adversarial",
    )
    payload = json.loads(result_json.read_text())
    cv = payload["cv_isolation"]
    summary = {
        "single_file_winner_identity_invariant": payload["single_file_ties"][
            "winner_identity_invariant"
        ],
        "single_file_score_and_q_invariant": payload["single_file_ties"][
            "score_and_q_invariant"
        ],
        "single_file_target_pep_invariant": payload["single_file_ties"][
            "target_pep_invariant"
        ],
        "single_file_decoy_pep_invariant": payload["single_file_ties"][
            "decoy_pep_invariant"
        ],
        "joined_fixed_name_permutation_invariant": payload["joined_inputs"][
            "file_row_target_decoy_permutation_invariant"
        ],
        "joined_path_alias_invariant": payload["joined_inputs"]["path_alias"][
            "invariant"
        ],
        "duplicate_false_discoveries": payload["candidate_multiplicity"][
            "false_discoveries_created_by_exact_duplicate_rows"
        ],
        "cv_leakage_detected": {
            name: values["leakage_detected"] for name, values in sorted(cv.items())
        },
        "protein_mapping_changes_with_seed": payload["protein_repaired_areas"][
            "seed_changes_protein_mapping"
        ],
    }
    probes = {
        source: compile_and_run_probe(source, commands)
        for source in ("final_repair_stats_probe.rs", "final_repair_protein_probe.rs")
    }
    return {
        "execution": record,
        "result": file_record(result_json),
        "summary": summary,
        "probes": probes,
    }


def run_one_file(binary: Path, pin: Path, destination: Path, threads: int) -> dict[str, Any]:
    destination.mkdir(parents=True, exist_ok=True)
    command = output_command(
        destination,
        ["--no-select-c", "--num-threads", str(threads)],
        [str(pin)],
        proteins=False,
    )
    started = time.perf_counter()
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    wall = time.perf_counter() - started
    if result.returncode:
        raise RuntimeError(f"benchmark failed: {command_text(command)}\n{result.stderr}")
    hashes = {name: sha256(destination / name) for name in OUTPUT_NAMES}
    return {
        "input": str(pin),
        "input_bytes": pin.stat().st_size,
        "threads": threads,
        "wall_seconds": wall,
        "output_sha256": hashes,
        "target_psms_q01": int(Q01_PSM.search(result.stderr).group(1)),
        "target_peptides_q01": int(Q01_PEPTIDE.search(result.stderr).group(1)),
    }


def balanced_order(inputs: list[Path]) -> list[Path]:
    ordered: list[Path] = []
    left, right = 0, len(inputs) - 1
    while left <= right:
        ordered.append(inputs[left])
        if left != right:
            ordered.append(inputs[right])
        left += 1
        right -= 1
    return ordered


def benchmark_full(binary: Path, inputs: list[Path], concurrency: int) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="percolator-baseline-full-") as temporary:
        root = Path(temporary)

        def one(pin: Path) -> tuple[str, dict[str, Any]]:
            return pin.name, run_one_file(binary, pin, root / pin.stem, 1)

        started = time.perf_counter()
        if concurrency == 1:
            pairs = [one(pin) for pin in inputs]
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as pool:
                pairs = list(pool.map(one, inputs))
        wall = time.perf_counter() - started
        by_file = dict(sorted(pairs))
        return {
            "processes": concurrency,
            "wall_seconds": wall,
            "files": len(by_file),
            "target_psms_q01": sum(item["target_psms_q01"] for item in by_file.values()),
            "target_peptides_q01": sum(
                item["target_peptides_q01"] for item in by_file.values()
            ),
            "output_sha256": {
                name: {
                    output: values["output_sha256"][output]
                    for output in OUTPUT_NAMES
                }
                for name, values in by_file.items()
            },
        }


def summarize_timings(runs: list[dict[str, Any]]) -> dict[str, Any]:
    walls = [run["wall_seconds"] for run in runs]
    return {
        "runs": runs,
        "wall_seconds": walls,
        "median_wall_seconds": statistics.median(walls),
    }


def capture_benchmarks(repeats: int) -> dict[str, Any]:
    input_root = ROOT / "data/PXD032157"
    inputs = sorted(input_root.glob("*.pin"), key=lambda path: (-path.stat().st_size, path.name))
    if len(inputs) != 65:
        raise RuntimeError(f"expected 65 benchmark PINs, found {len(inputs)}")
    largest = inputs[0]
    order = balanced_order(inputs)
    single_t1 = []
    single_t3 = []
    full_sequential = []
    full_n4 = []
    for repeat in range(1, repeats + 1):
        with tempfile.TemporaryDirectory(prefix="percolator-baseline-largest-") as temporary:
            root = Path(temporary)
            single_t1.append(run_one_file(BINARY, largest, root / "t1", 1))
            single_t1[-1]["repeat"] = repeat
            single_t3.append(run_one_file(BINARY, largest, root / "t3", 3))
            single_t3[-1]["repeat"] = repeat
        sequential = benchmark_full(BINARY, order, 1)
        sequential["repeat"] = repeat
        full_sequential.append(sequential)
        n4 = benchmark_full(BINARY, order, 4)
        n4["repeat"] = repeat
        full_n4.append(n4)
        print(
            f"benchmark repeat {repeat}/{repeats}: "
            f"largest t1={single_t1[-1]['wall_seconds']:.3f}s "
            f"t3={single_t3[-1]['wall_seconds']:.3f}s "
            f"full seq={sequential['wall_seconds']:.3f}s "
            f"n4={n4['wall_seconds']:.3f}s",
            flush=True,
        )
    return {
        "largest_input": {"path": str(largest), **file_record(largest)},
        "input_files": len(inputs),
        "input_bytes": sum(path.stat().st_size for path in inputs),
        "largest_pin_t1": summarize_timings(single_t1),
        "largest_pin_t3": summarize_timings(single_t3),
        "full_sequential": summarize_timings(full_sequential),
        "full_n4": summarize_timings(full_n4),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--full-benchmarks", action="store_true")
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite existing output: {output}")
    assert_production_clean()
    output.mkdir(parents=True)
    commands = output / "commands"

    status = git_text("status", "--porcelain=v1")
    (output / "git-status.txt").write_text(status)
    diff = git_text("diff", "--binary")
    (output / "pre-existing.patch").write_text(diff)

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "purpose": "behavior-preserving refactor baseline",
        "git": {
            "head": git_text("rev-parse", "HEAD").strip(),
            "status": "git-status.txt",
            "status_sha256": sha256(output / "git-status.txt"),
            "tracked_diff": "pre-existing.patch",
            "tracked_diff_sha256": sha256(output / "pre-existing.patch"),
        },
        "environment": {
            "rustc": subprocess.run(
                ["rustc", "--version", "--verbose"], text=True, capture_output=True, check=True
            ).stdout,
            "cargo": subprocess.run(
                ["cargo", "--version", "--verbose"], text=True, capture_output=True, check=True
            ).stdout,
            "platform": dict(zip(
                ("sysname", "nodename", "release", "version", "machine"),
                os.uname(),
            )),
        },
        "fixture": {"path": str(FIXTURE), **file_record(FIXTURE)},
    }

    manifest["build"] = execute(
        ["cargo", "build", "--release", "--locked"], commands, "cargo-build-release"
    )
    manifest["binary"] = {"path": str(BINARY), **file_record(BINARY)}
    manifest["tests"] = {
        "cargo_release_all_targets": execute(
            ["cargo", "test", "--release", "--all-targets", "--locked"],
            commands,
            "cargo-test-release-all-targets",
        ),
        "shell_gates": {
            gate: execute(["bash", gate], commands, Path(gate).stem)
            for gate in SHELL_GATES
        },
    }
    manifest["outputs"] = capture_outputs(output / "outputs", commands)
    manifest["adversarial"] = capture_adversarial(output, commands)
    if args.full_benchmarks:
        manifest["benchmarks"] = capture_benchmarks(args.repeats)

    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"baseline complete: {manifest_path}")
    print(f"manifest SHA-256: {sha256(manifest_path)}")


if __name__ == "__main__":
    main()
