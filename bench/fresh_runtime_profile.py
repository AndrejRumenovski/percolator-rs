#!/usr/bin/env python3
"""Collect a fresh, behavior-preserving runtime profile of PXD032157.

The production binary, the feature-gated stage-timed binary, and a sampled CPU
binary are built in isolated target directories.  Every run writes the normal
TSV outputs; their hashes are retained and the bulky duplicate TSV payloads are
then removed from the artifact directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent.parent
Q01_PSM = re.compile(r"target PSMs q<0\.01: (\d+)")
Q01_PEPTIDE = re.compile(r"target peptides q<0\.01: (\d+)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_capture(command: list[str], path: Path, env: dict[str, str] | None = None) -> None:
    result = subprocess.run(
        command,
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    path.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}")


def proc_sample(pid: int) -> tuple[int, int, int] | None:
    """Return current RSS bytes and cumulative minor/major faults."""
    try:
        raw = Path(f"/proc/{pid}/stat").read_text()
    except (FileNotFoundError, ProcessLookupError):
        return None
    fields = raw[raw.rfind(")") + 2 :].split()
    page_size = os.sysconf("SC_PAGE_SIZE")
    return int(fields[21]) * page_size, int(fields[7]), int(fields[9])


@dataclass
class Running:
    process: subprocess.Popen[bytes]
    stdout: Any
    stderr: Any
    destination: Path
    command: list[str]
    peak_rss_bytes: int = 0
    minor_faults: int = 0
    major_faults: int = 0


def output_command(
    binary: Path,
    pins: list[Path],
    threads: int,
    destination: Path,
    profile_json: Path | None,
    cpu_prefix: Path | None,
    allocations: bool,
    proteins: bool,
    extra_args: list[str],
) -> list[str]:
    command = [
        str(binary),
        "--canonical",
        "--seed",
        "1",
        "--num-threads",
        str(threads),
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
        command += [
            "--results-proteins",
            str(destination / "target.proteins.tsv"),
            "--decoy-results-proteins",
            str(destination / "decoy.proteins.tsv"),
        ]
    if profile_json is not None:
        command += ["--profile-json", str(profile_json)]
    if cpu_prefix is not None:
        command += ["--profile-cpu", str(cpu_prefix)]
    if allocations:
        command.append("--profile-allocations")
    command += extra_args
    command += [str(pin) for pin in pins]
    return command


def start_one(command: list[str], destination: Path) -> Running:
    destination.mkdir(parents=True, exist_ok=True)
    stdout = (destination / "stdout.log").open("wb")
    stderr = (destination / "stderr.log").open("wb")
    process = subprocess.Popen(command, cwd=ROOT, stdout=stdout, stderr=stderr)
    return Running(process, stdout, stderr, destination, command)


def finish_one(run: Running) -> dict[str, Any]:
    if run.process.returncode:
        stderr = (run.destination / "stderr.log").read_text(errors="replace")
        raise RuntimeError(
            f"profile command failed ({run.process.returncode}): {' '.join(run.command)}\n{stderr}"
        )
    stderr = (run.destination / "stderr.log").read_text(errors="replace")
    psm = Q01_PSM.search(stderr)
    peptide = Q01_PEPTIDE.search(stderr)
    hashes: dict[str, dict[str, Any]] = {}
    for path in sorted(run.destination.glob("*.tsv")):
        hashes[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
        path.unlink()
    return {
        "command": run.command,
        "peak_rss_bytes": run.peak_rss_bytes,
        "minor_faults": run.minor_faults,
        "major_faults": run.major_faults,
        "target_psms_q_lt_0_01": int(psm.group(1)) if psm else None,
        "target_peptides_q_lt_0_01": int(peptide.group(1)) if peptide else None,
        "outputs": hashes,
    }


def run_processes(specifications: list[tuple[list[str], Path]], concurrency: int) -> dict[str, Any]:
    pending = list(specifications)
    active: list[Running] = []
    completed: list[Running] = []
    peak_tree_rss = 0
    started = time.perf_counter_ns()
    while pending or active:
        while pending and len(active) < concurrency:
            command, destination = pending.pop(0)
            active.append(start_one(command, destination))
        current_tree_rss = 0
        for run in active:
            sample = proc_sample(run.process.pid)
            if sample is not None:
                rss, minor, major = sample
                current_tree_rss += rss
                run.peak_rss_bytes = max(run.peak_rss_bytes, rss)
                run.minor_faults = max(run.minor_faults, minor)
                run.major_faults = max(run.major_faults, major)
        peak_tree_rss = max(peak_tree_rss, current_tree_rss)
        still_active = []
        for run in active:
            status = run.process.poll()
            if status is None:
                still_active.append(run)
            else:
                run.stdout.close()
                run.stderr.close()
                completed.append(run)
        active = still_active
        if active:
            time.sleep(0.02)
    wall_ns = time.perf_counter_ns() - started
    # Hashing proves output identity but is benchmark-harness work, not scorer
    # runtime.  Finalize only after the measured interval has ended.
    finished = [finish_one(run) for run in completed]
    return {
        "wall_ns": wall_ns,
        "peak_tree_rss_bytes": peak_tree_rss,
        "minor_faults": sum(item["minor_faults"] for item in finished),
        "major_faults": sum(item["major_faults"] for item in finished),
        "processes": finished,
    }


class Campaign:
    def __init__(self, artifacts: Path, inputs: list[Path]) -> None:
        self.artifacts = artifacts
        self.inputs = inputs
        self.largest = inputs[0]
        self.binaries: dict[str, Path] = {}
        self.timing_path = artifacts / "timings.tsv"
        self.timing_path.write_text(
            "configuration\tbuild\trepetition\twall_ns\tprocesses\t"
            "intra_file_threads\tcpu_sampling\tpeak_tree_rss_bytes\t"
            "minor_faults\tmajor_faults\n"
        )

    def record(
        self,
        configuration: str,
        artifact_configuration: str,
        build: str,
        repetition: int,
        pins: list[Path],
        concurrency: int,
        threads: int,
        allocations: bool = False,
        proteins: bool = False,
        extra_args: list[str] | None = None,
    ) -> None:
        tag = f"{artifact_configuration}_{build}_r{repetition}"
        output_root = self.artifacts / "outputs" / tag
        profile_root = self.artifacts / "profiles" / tag
        manifest_path = self.artifacts / "manifests" / f"{tag}.json"
        output_root.mkdir(parents=True)
        profile_root.mkdir(parents=True)
        specifications = []
        joined = "--join" in (extra_args or [])
        work = [pins] if joined else [[pin] for pin in pins]
        for group in work:
            stem = "joined" if joined else group[0].stem
            destination = output_root / stem
            profile_json = None
            cpu_prefix = None
            if build in {"instrumented", "cpu"}:
                profile_json = profile_root / f"{stem}.json"
            if build == "cpu":
                cpu_prefix = profile_root / f"{stem}.cpu"
            command = output_command(
                self.binaries[build],
                group,
                threads,
                destination,
                profile_json,
                cpu_prefix,
                allocations,
                proteins,
                extra_args or [],
            )
            specifications.append((command, destination))
        print(
            f"{configuration} {build} r{repetition}: "
            f"{len(specifications)} process(es), N={concurrency}, threads={threads}",
            flush=True,
        )
        result = run_processes(specifications, concurrency)
        manifest_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        with self.timing_path.open("a") as handle:
            handle.write(
                f"{configuration}\t{build}\t{repetition}\t{result['wall_ns']}\t"
                f"{concurrency}\t{threads}\t{int(build == 'cpu')}\t"
                f"{result['peak_tree_rss_bytes']}\t{result['minor_faults']}\t"
                f"{result['major_faults']}\n"
            )


def build_binaries(artifacts: Path) -> dict[str, Path]:
    build = artifacts / "build"
    binaries = artifacts / "bin"
    binaries.mkdir()
    targets = {
        "normal": build / "target-normal",
        "instrumented": build / "target-instrumented",
        "cpu": build / "target-cpu",
    }
    commands = {
        "normal": ["cargo", "build", "--release", "--locked"],
        "instrumented": [
            "cargo",
            "build",
            "--release",
            "--locked",
            "--features",
            "profiling",
        ],
        "cpu": [
            "cargo",
            "build",
            "--release",
            "--locked",
            "--features",
            "profiling",
        ],
    }
    result = {}
    for name in ("normal", "instrumented", "cpu"):
        command = commands[name] + ["--target-dir", str(targets[name])]
        env = os.environ.copy()
        if name == "cpu":
            existing = env.get("RUSTFLAGS", "")
            env["RUSTFLAGS"] = (
                f"{existing} -C force-frame-pointers=yes -C debuginfo=1"
            ).strip()
        run_capture(command, build / f"cargo-build-{name}.txt", env)
        source = targets[name] / "release" / "percolator-rs"
        destination = binaries / f"percolator-rs-{name}"
        shutil.copy2(source, destination)
        result[name] = destination
    hashes = {name: sha256(path) for name, path in result.items()}
    (build / "binary-sha256.json").write_text(json.dumps(hashes, indent=2) + "\n")
    return result


def rotated(values: list[int], repetition: int) -> list[int]:
    shift = (repetition - 1) % len(values)
    ordered = values[shift:] + values[:shift]
    return ordered if repetition % 2 else list(reversed(ordered))


def build_order(repetition: int, include_cpu: bool) -> list[str]:
    values = ["normal", "instrumented", "cpu"] if include_cpu else ["normal", "instrumented"]
    shift = (repetition - 1) % len(values)
    return values[shift:] + values[:shift]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", required=True, type=Path)
    parser.add_argument("--input", type=Path, default=ROOT / "data/PXD032157")
    parser.add_argument("--single-repeats", type=int, default=5)
    parser.add_argument("--full-repeats", type=int, default=3)
    parser.add_argument("--cpu-repeats", type=int, default=3)
    parser.add_argument(
        "--reuse-binaries",
        type=Path,
        help="directory containing percolator-rs-{normal,instrumented,cpu}",
    )
    parser.add_argument(
        "--timings-only",
        action="store_true",
        help="skip allocation and conditional-path runs",
    )
    args = parser.parse_args()
    if min(args.single_repeats, args.full_repeats, args.cpu_repeats) < 1:
        parser.error("repeat counts must be positive")
    artifacts = args.artifacts.resolve()
    if artifacts.exists():
        parser.error(f"refusing to overwrite {artifacts}")
    artifacts.mkdir(parents=True)
    for name in ("build", "profiles", "outputs", "manifests"):
        (artifacts / name).mkdir(exist_ok=True)

    inputs = sorted(args.input.resolve().glob("*.pin"), key=lambda path: (-path.stat().st_size, path.name))
    if len(inputs) != 65:
        parser.error(f"expected 65 PIN files, found {len(inputs)}")
    balanced = []
    left, right = 0, len(inputs) - 1
    while left <= right:
        balanced.append(inputs[left])
        if left != right:
            balanced.append(inputs[right])
        left += 1
        right -= 1
    (artifacts / "input-order.txt").write_text("\n".join(map(str, balanced)) + "\n")

    build = artifacts / "build"
    run_capture(["git", "rev-parse", "HEAD"], build / "git-head.txt")
    run_capture(["git", "status", "--porcelain=v1"], build / "git-status.txt")
    run_capture(["git", "diff", "--binary"], build / "instrumentation.patch")
    run_capture(["rustc", "-Vv"], build / "rustc.txt")
    run_capture(["cargo", "-V"], build / "cargo.txt")
    run_capture(["uname", "-a"], build / "uname.txt")
    run_capture(["lscpu"], build / "lscpu.txt")
    (build / "cargo-config.toml").write_text((ROOT / ".cargo/config.toml").read_text())
    perf = subprocess.run(
        ["perf", "stat", "-e", "task-clock", "--", "true"],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (build / "perf-probe.txt").write_text(perf.stdout)
    (build / "perf-probe-exit-code.txt").write_text(f"{perf.returncode}\n")
    run_capture(
        [sys.executable, "refactor/verify_baseline.py"],
        build / "verify-baseline.txt",
    )

    campaign = Campaign(artifacts, balanced)
    if args.reuse_binaries is None:
        campaign.binaries = build_binaries(artifacts)
    else:
        source = args.reuse_binaries.resolve()
        campaign.binaries = {
            name: source / f"percolator-rs-{name}"
            for name in ("normal", "instrumented", "cpu")
        }
        missing = [str(path) for path in campaign.binaries.values() if not path.is_file()]
        if missing:
            parser.error(f"missing reused binaries: {', '.join(missing)}")
        (build / "binary-sha256.json").write_text(
            json.dumps(
                {name: sha256(path) for name, path in campaign.binaries.items()},
                indent=2,
            )
            + "\n"
        )

    # Largest-file scaling and paired stage/CPU overhead.  Thread counts are
    # rotated and reversed across repetitions; build order uses a Latin rotation.
    for repetition in range(1, args.single_repeats + 1):
        for threads in rotated([1, 2, 3, 6], repetition):
            config = f"single_file_t{threads}"
            if threads in {1, 3}:
                for binary in build_order(repetition, repetition <= args.cpu_repeats):
                    artifact_config = config if binary != "cpu" else f"{config}_cpu"
                    campaign.record(
                        config,
                        artifact_config,
                        binary,
                        repetition,
                        [campaign.largest],
                        1,
                        threads,
                    )
            else:
                campaign.record(config, config, "normal", repetition, [campaign.largest], 1, threads)

    # Full-corpus scaling.  N=1 and N=4 also provide balanced normal/timed
    # pairs; N=4 has repeated sampled builds for an overhead median.
    for repetition in range(1, args.full_repeats + 1):
        for concurrency in rotated([1, 2, 4, 6], repetition):
            config = "full_sequential" if concurrency == 1 else f"full_n{concurrency}"
            if concurrency in {1, 4}:
                include_cpu = concurrency == 4 and repetition <= args.cpu_repeats
                for binary in build_order(repetition, include_cpu):
                    artifact_config = config if binary != "cpu" else f"{config}_cpu"
                    campaign.record(
                        config,
                        artifact_config,
                        binary,
                        repetition,
                        balanced,
                        concurrency,
                        1,
                    )
                if concurrency == 1 and repetition == 1:
                    campaign.record(
                        config,
                        f"{config}_cpu",
                        "cpu",
                        repetition,
                        balanced,
                        1,
                        1,
                    )
            else:
                campaign.record(config, config, "normal", repetition, balanced, concurrency, 1)

    if args.timings_only:
        print(f"fresh runtime timing matrix complete: {artifacts}", flush=True)
        return

    # Allocation counting is isolated from timing conclusions.
    campaign.record(
        "single_file_t1_allocations",
        "single_file_t1_allocations",
        "instrumented",
        1,
        [campaign.largest],
        1,
        1,
        allocations=True,
    )
    campaign.record(
        "full_sequential_allocations",
        "full_sequential_allocations",
        "instrumented",
        1,
        balanced,
        1,
        1,
        allocations=True,
    )
    campaign.record(
        "full_n4_allocations",
        "full_n4_allocations",
        "instrumented",
        1,
        balanced,
        4,
        1,
        allocations=True,
    )

    # Conditional paths use workloads where the feature is genuinely enabled.
    campaign.record(
        "single_file_t1_rt",
        "single_file_t1_rt",
        "instrumented",
        1,
        [campaign.largest],
        1,
        1,
        extra_args=["--rt-features"],
    )
    campaign.record(
        "joined_two_file_t1",
        "joined_two_file_t1",
        "instrumented",
        1,
        list(reversed(inputs[-2:])),
        1,
        1,
        extra_args=["--join"],
    )
    protein_pin = ROOT / "data/F_3.pin"
    if protein_pin.exists():
        campaign.record(
            "protein_f3_t1",
            "protein_f3_t1",
            "instrumented",
            1,
            [protein_pin],
            1,
            1,
            proteins=True,
        )

    print(f"fresh runtime profile complete: {artifacts}", flush=True)


if __name__ == "__main__":
    main()
