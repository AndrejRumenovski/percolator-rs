#!/usr/bin/env python3
"""Adversarial complete-null variants, added after the 2026-08-25 repair.

These are **new experiments**, not part of the predeclared complete-null study in
`run_null.py`. They exist to try to falsify the repaired estimator on null
constructions it was never checked against:

* ``source`` — which rows the pseudo-labels are assigned to. ``decoy`` reproduces
  the predeclared construction (original decoy rows only). ``target`` uses the
  original *target* rows instead: they are heterogeneous in quality, but a label
  assigned at random is still exchangeable with the features, so the complete
  null holds and no learner can beat chance.
* ``ratio`` — pseudo-targets per pseudo-decoy. The predeclared study is exactly
  balanced. An imbalanced relabeling stresses the part of the estimator that the
  repair replaced: a k:1 relabeling means an incorrect target outranks its paired
  decoy with probability k/(k+1), which is what ``--null-target-win-prob`` is for.
  Running an imbalanced null at the default 0.5 *should* be anti-conservative;
  running it at the declared value should not. That contrast is the test.

Under a complete null every accepted pseudo-target is false, so FDP is 1 whenever
anything is accepted and 0 otherwise, and empirical FDR is the probability of any
acceptance across replicates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import random
import shlex
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


def make_null(source: Path, destination: Path, seed: int, keep: str, ratio: float) -> dict:
    """Relabel a fixed row set at random, keeping the requested target share."""
    with source.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        label_column = next(
            (i for i, name in enumerate(header) if name.lower() == "label"), None
        )
        if label_column is None:
            raise ValueError(f"{source}: no Label column")
        want_target = keep == "target"
        rows = []
        for row in reader:
            if row and row[0] == "DefaultDirection":
                continue
            if len(row) <= label_column:
                raise ValueError(f"{source}: short row in input")
            if (int(row[label_column]) > 0) == want_target:
                rows.append(row)
    if len(rows) < 100:
        raise ValueError(f"{source}: fewer than 100 source rows of class {keep}")

    order = list(range(len(rows)))
    random.Random(seed).shuffle(order)
    target_share = ratio / (1.0 + ratio)
    cut = int(round(len(rows) * target_share))
    targets = set(order[:cut])
    for index, row in enumerate(rows):
        row[label_column] = "1" if index in targets else "-1"
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return {
        "source_rows_kept": keep,
        "source_rows": len(rows),
        "pseudo_targets": len(targets),
        "pseudo_decoys": len(rows) - len(targets),
        "requested_target_decoy_ratio": ratio,
        "relabel_seed": seed,
        "sha256": sha256(destination),
    }


def q_counts(path: Path) -> dict:
    values = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            q = float(row["q-value"])
            if not math.isfinite(q) or not 0 <= q <= 1:
                raise ValueError(f"{path}: invalid q-value {q}")
            values.append(q)
    return {
        f"q_lt_{threshold:g}": sum(q < threshold for q in values)
        for threshold in THRESHOLDS
    }


def wilson(successes: int, trials: int) -> tuple[float, float]:
    if trials == 0:
        return 0.0, 0.0
    z = 1.959963984540054
    estimate = successes / trials
    denominator = 1 + z * z / trials
    centre = (estimate + z * z / (2 * trials)) / denominator
    half = (
        z
        * math.sqrt(
            estimate * (1 - estimate) / trials + z * z / (4 * trials * trials)
        )
        / denominator
    )
    return max(0.0, centre - half), min(1.0, centre + half)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--rust", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--relabel-seeds", default="2001,2002,2003,2004,2005,2006,2007,2008,2009,2010")
    parser.add_argument("--model-seed", type=int, default=1)
    args = parser.parse_args()

    seeds = [int(value) for value in args.relabel_seeds.split(",")]
    if args.output.exists():
        parser.error(f"output already exists: {args.output}")
    args.output.mkdir(parents=True)

    # (label, keep, ratio, extra argv) -- each arm is declared here before it runs.
    arms = [
        ("decoy_balanced", "decoy", 1.0, []),
        ("target_balanced", "target", 1.0, []),
        ("decoy_2to1_default_p", "decoy", 2.0, []),
        ("decoy_2to1_declared_p", "decoy", 2.0, ["--null-target-win-prob", str(2.0 / 3.0)]),
    ]

    repo = Path(__file__).resolve().parents[1]
    manifest = {
        "schema_version": 1,
        "study": "post-repair adversarial complete-null variants",
        "status": "NEW EXPERIMENT, added after the repair; not part of run_null.py",
        "null_interpretation": (
            "Pseudo-labels are assigned at random to a fixed row set, so labels are exchangeable "
            "with the features and every accepted pseudo-target is false. FDP is 1 if anything is "
            "accepted and 0 otherwise; empirical FDR is the probability of any acceptance."
        ),
        "threshold_policy": "q < threshold",
        "thresholds": THRESHOLDS,
        "relabel_seeds": seeds,
        "model_seed": args.model_seed,
        "arms": [
            {"id": name, "rows_kept": keep, "target_decoy_ratio": ratio, "extra_argv": extra}
            for name, keep, ratio, extra in arms
        ],
        "environment": {
            "platform": platform.platform(),
            "python": sys.version,
            "commit": subprocess.run(
                ["git", "-C", str(repo), "rev-parse", "HEAD"],
                text=True, capture_output=True, check=True,
            ).stdout.strip(),
            "rust_binary": str(args.rust.resolve()),
            "rust_binary_sha256": sha256(args.rust.resolve()),
        },
        "evaluation_script_sha256": sha256(Path(__file__).resolve()),
        "inputs": [
            {"path": str(p.resolve()), "sha256": sha256(p.resolve())} for p in args.input
        ],
        "runs": [],
    }

    for source in args.input:
        source = source.resolve()
        for name, keep, ratio, extra in arms:
            for seed in seeds:
                root = args.output / source.stem / name / f"relabel-{seed}"
                root.mkdir(parents=True)
                null_pin = root / "null.pin"
                construction = make_null(source, null_pin, seed, keep, ratio)
                command = [
                    str(args.rust.resolve()), "--canonical", "--no-select-c",
                    "--seed", str(args.model_seed), "--num-threads", "1",
                    *extra,
                    "--results-psms", str(root / "target.psms.tsv"),
                    "--decoy-results-psms", str(root / "decoy.psms.tsv"),
                    str(null_pin),
                ]
                started = time.time()
                with (root / "stdout.log").open("w") as out, (root / "stderr.log").open("w") as err:
                    result = subprocess.run(command, stdout=out, stderr=err, check=False)
                record = {
                    "input": source.name, "arm": name, "relabel_seed": seed,
                    "construction": construction,
                    "argv": command, "shell_display": shlex.join(command),
                    "exit_code": result.returncode,
                    "wall_seconds": time.time() - started,
                    "status": "complete" if result.returncode == 0 else "failed",
                }
                if record["status"] == "complete":
                    record["false_target_psms"] = q_counts(root / "target.psms.tsv")
                manifest["runs"].append(record)
                print(
                    f"{source.stem} {name} relabel={seed}: "
                    + (
                        f"q<0.01 false targets={record['false_target_psms']['q_lt_0.01']}"
                        if record["status"] == "complete"
                        else f"FAILED exit={result.returncode}"
                    ),
                    flush=True,
                )

    aggregate = []
    for name, _, _, _ in arms:
        runs = [r for r in manifest["runs"] if r["arm"] == name and r["status"] == "complete"]
        attempted = [r for r in manifest["runs"] if r["arm"] == name]
        for threshold in THRESHOLDS:
            key = f"q_lt_{threshold:g}"
            values = [r["false_target_psms"][key] for r in runs]
            positives = sum(value > 0 for value in values)
            low, high = wilson(positives, len(values))
            aggregate.append({
                "arm": name, "q_threshold": threshold,
                "replicates": len(values),
                "attempted_replicates": len(attempted),
                "replicates_with_any_false_discovery": positives,
                "empirical_fdr": positives / len(values) if values else None,
                "empirical_fdr_wilson95_low": low,
                "empirical_fdr_wilson95_high": high,
                "max_false_targets": max(values) if values else None,
            })
    manifest["aggregate"] = aggregate
    manifest["completed_unix"] = time.time()
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(f"machine-readable manifest: {args.output / 'manifest.json'}")


if __name__ == "__main__":
    main()
