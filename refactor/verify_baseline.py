#!/usr/bin/env python3
"""Verify the current checkout against the frozen pre-refactor behavior."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Any

import freeze_baseline as frozen


DEFAULT_BASELINE = Path(__file__).with_name("baseline") / "e8d83d1"


def compare_files(
    observed: dict[str, Any], expected: dict[str, Any]
) -> list[str]:
    failures: list[str] = []
    for variant in ("fixed-serial", "fixed-parallel", "select-c", "ensemble"):
        observed_files = observed[variant]["files"]
        expected_files = expected[variant]["files"]
        if observed_files.keys() != expected_files.keys():
            failures.append(
                f"{variant}: file set {sorted(observed_files)} != "
                f"{sorted(expected_files)}"
            )
            continue
        for name, expected_record in expected_files.items():
            observed_record = observed_files[name]
            if observed_record != expected_record:
                failures.append(
                    f"{variant}/{name}: {observed_record} != {expected_record}"
                )
    for invariant in (
        "fixed_serial_parallel_byte_identical",
        "fixed_repeat_byte_identical",
    ):
        if observed[invariant] != expected[invariant]:
            failures.append(
                f"{invariant}: {observed[invariant]} != {expected[invariant]}"
            )
    return failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        type=Path,
        default=DEFAULT_BASELINE,
        help="frozen artifact directory (default: %(default)s)",
    )
    parser.add_argument(
        "--skip-tests",
        action="store_true",
        help="skip the release Cargo suite and portable shell gates",
    )
    parser.add_argument(
        "--skip-adversarial",
        action="store_true",
        help="skip the adversarial driver and standalone probes",
    )
    args = parser.parse_args()

    baseline = args.baseline.resolve()
    manifest = json.loads((baseline / "manifest.json").read_text())
    failures: list[str] = []

    with tempfile.TemporaryDirectory(prefix="percolator-refactor-verify-") as temporary:
        artifacts = Path(temporary)
        commands = artifacts / "commands"
        frozen.execute(
            ["cargo", "build", "--release", "--locked"],
            commands,
            "cargo-build-release",
        )

        if not args.skip_tests:
            frozen.execute(
                ["cargo", "test", "--release", "--all-targets", "--locked"],
                commands,
                "cargo-test-release-all-targets",
            )
            for gate in frozen.SHELL_GATES:
                frozen.execute(
                    ["bash", gate],
                    commands,
                    Path(gate).stem,
                )

        observed_outputs = frozen.capture_outputs(artifacts / "outputs", commands)
        failures.extend(compare_files(observed_outputs, manifest["outputs"]))

        if not args.skip_adversarial:
            observed_adversarial = frozen.capture_adversarial(artifacts, commands)
            if observed_adversarial["summary"] != manifest["adversarial"]["summary"]:
                failures.append(
                    "adversarial summary changed:\n"
                    + json.dumps(observed_adversarial["summary"], indent=2, sort_keys=True)
                )

    if failures:
        raise SystemExit("baseline verification failed:\n- " + "\n- ".join(failures))
    print("release tests, shell gates, exact outputs, and adversarial behavior match")


if __name__ == "__main__":
    main()
