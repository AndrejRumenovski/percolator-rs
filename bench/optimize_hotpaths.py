#!/usr/bin/env python3
"""Paired hot-path experiments using the fresh-profile collector and exact hashes.

Build/copy each revision before running; compilation must not overlap timings.
Artifacts are new directories. TSV hashing is outside the measured interval.
"""
from __future__ import annotations

import argparse
import filecmp
import json
import statistics
import sys
from pathlib import Path

from fresh_runtime_profile import Campaign, ROOT, sha256
from runtime_profile_report import aggregate_profiles

sys.path.insert(0, str(ROOT / "refactor"))
import freeze_baseline as frozen


def verify_frozen(binary: Path, artifacts: Path) -> None:
    frozen.BINARY = binary.resolve()
    frozen.capture_outputs(artifacts / "outputs", artifacts / "commands")
    expected = ROOT / "refactor/baseline/e8d83d1/outputs"
    for path in expected.glob("*/*.tsv"):
        observed = artifacts / "outputs" / path.relative_to(expected)
        assert filecmp.cmp(path, observed, shallow=False), f"frozen byte mismatch: {path}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--profile", type=Path)
    parser.add_argument("--artifacts", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--single", action="store_true")
    parser.add_argument("--concurrency", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--profile-concurrency", type=int, nargs="+")
    args = parser.parse_args()
    if args.repeats < 1 or any(n < 1 for n in args.concurrency + (args.profile_concurrency or [])):
        parser.error("repeat and concurrency counts must be positive")
    args.artifacts.mkdir(parents=True, exist_ok=False)
    for name in ("profiles", "outputs", "manifests"):
        (args.artifacts / name).mkdir()
    inputs = sorted((ROOT / "data/PXD032157").glob("*.pin"),
                    key=lambda p: (-p.stat().st_size, p.name))
    assert len(inputs) == 65
    balanced = []
    while inputs:
        balanced.append(inputs.pop(0))
        if inputs:
            balanced.append(inputs.pop())
    if args.single:
        balanced = balanced[:1]
    campaign = Campaign(args.artifacts, balanced)
    campaign.binaries = {"before": args.before.resolve(), "after": args.after.resolve()}
    if args.profile:
        campaign.binaries["instrumented"] = args.profile.resolve()
    summary = {"binaries": {k: {"path": str(v), "sha256": sha256(v)}
                            for k, v in campaign.binaries.items()}, "timings": {}}
    for build in ("before", "after"):
        verify_frozen(campaign.binaries[build], args.artifacts / f"frozen-{build}")
    summary["frozen_outputs_byte_identical"] = True
    oracle = None
    yield_oracle = None
    for repetition in range(1, args.repeats + 1):
        for n in (args.concurrency if repetition % 2 else args.concurrency[::-1]):
            config = f"{'single' if args.single else 'full'}_n{n}"
            for build in (["before", "after"] if repetition % 2 else ["after", "before"]):
                campaign.record(config, config, build, repetition, balanced, n, 1)
                manifest = json.loads((args.artifacts / "manifests" /
                                      f"{config}_{build}_r{repetition}.json").read_text())
                hashes = {Path(p["command"][-1]).name: p["outputs"] for p in manifest["processes"]}
                if oracle is None:
                    oracle = hashes
                assert hashes == oracle, f"output mismatch: {config} {build} {repetition}"
                yields = {key: sum(p[key] for p in manifest["processes"]) for key in
                          ("target_psms_q_lt_0_01", "target_peptides_q_lt_0_01")}
                if yield_oracle is None:
                    yield_oracle = yields
                assert yields == yield_oracle, f"yield mismatch: {config} {build} {repetition}"
                seconds = manifest["wall_ns"] / 1e9
                summary["timings"].setdefault(config, {}).setdefault(build, []).append(seconds)
                print(f"  {seconds:.6f}s; all output bytes/hashes match", flush=True)
    for config, builds in summary["timings"].items():
        medians = {b: statistics.median(v) for b, v in builds.items()}
        print(f"{config}: {medians}, improvement {100*(1-medians['after']/medians['before']):.2f}%", flush=True)
    if args.profile:
        for n in (args.profile_concurrency or args.concurrency):
            config = f"{'single' if args.single else 'full'}_n{n}"
            campaign.record(config, config, "instrumented", 1, balanced, n, 1)
            manifest = json.loads((args.artifacts / "manifests" /
                                  f"{config}_instrumented_r1.json").read_text())
            assert {Path(p["command"][-1]).name: p["outputs"] for p in manifest["processes"]} == oracle
            paths = list((args.artifacts / "profiles" / f"{config}_instrumented_r1").glob("*.json"))
            profile = aggregate_profiles(paths)
            (args.artifacts / f"{config}-profile.json").write_text(json.dumps(profile, indent=2) + "\n")
    summary["yield"] = yield_oracle
    (args.artifacts / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
