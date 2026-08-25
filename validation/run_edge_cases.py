#!/usr/bin/env python3
"""Adversarial parser, numerical, tie, duplicate, and small-sample cases."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shlex
import subprocess
import time
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_pin(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = [row for row in reader if row and row[0] != "DefaultDirection"]
    return header, rows


def write_pin(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def output_checks(target: Path, decoy: Path) -> dict:
    rows = []
    for path in (target, decoy):
        if not path.exists():
            continue
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows.append({name: float(row[name]) for name in ("score", "q-value", "posterior_error_prob")})
    finite = all(math.isfinite(value) for row in rows for value in row.values())
    bounded_q = all(0 <= row["q-value"] <= 1 for row in rows)
    bounded_pep = all(0 <= row["posterior_error_prob"] <= 1 for row in rows)
    order = sorted(rows, key=lambda row: row["score"], reverse=True)
    monotone_q = all(order[index]["q-value"] <= order[index + 1]["q-value"] + 1e-12 for index in range(len(order) - 1))
    tie_values = {}
    for row in rows:
        tie_values.setdefault(row["score"], set()).add(row["q-value"])
    tie_invariant = all(len(values) == 1 for values in tie_values.values())
    return {
        "output_rows": len(rows), "finite_statistics": finite, "bounded_qvalues": bounded_q,
        "bounded_peps": bounded_pep, "q_monotone_by_printed_score": monotone_q,
        "equal_printed_scores_have_equal_q": tie_invariant,
        "exact_zero_qvalues": sum(row["q-value"] == 0 for row in rows),
        "exact_zero_peps": sum(row["posterior_error_prob"] == 0 for row in rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--rust", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"output exists: {args.output}")
    args.output.mkdir(parents=True)
    header, source = read_pin(args.fixture)
    label = next(i for i, name in enumerate(header) if name.lower() == "label")
    scan = next(i for i, name in enumerate(header) if name.lower() == "scannr")
    peptide = next(i for i, name in enumerate(header) if name.lower() == "peptide")
    feature_columns = [i for i in range(scan + 1, peptide) if header[i] not in ("ExpMass", "CalcMass")]
    first_feature = feature_columns[0]
    targets = [row[:] for row in source if int(row[label]) > 0]
    decoys = [row[:] for row in source if int(row[label]) < 0]

    cases = {}
    cases["very_small"] = targets[:3] + decoys[:3]
    cases["target_decoy_imbalance"] = targets[:100] + decoys[:10]
    cases["duplicate_psms"] = [row[:] for row in source[:1000]] + [row[:] for row in source[:1000]]
    tied = [row[:] for row in targets[:1000] + decoys[:1000]]
    for row in tied:
        for column in feature_columns:
            row[column] = "0"
    cases["all_features_tied"] = tied
    malformed = [row[:] for row in source[:2000]]
    malformed[0][first_feature] = "not-a-number"
    cases["malformed_feature"] = malformed
    missing = [row[:] for row in source[:2000]]
    missing[0][first_feature] = ""
    cases["missing_feature"] = missing
    nonfinite = [row[:] for row in source[:2000]]
    nonfinite[0][first_feature] = "NaN"
    cases["nan_feature"] = nonfinite
    extreme = [row[:] for row in source[:2000]]
    extreme[0][first_feature] = "1e308"
    extreme[1][first_feature] = "-1e308"
    cases["extreme_feature"] = extreme
    unusual = [row[:] for row in source[:2000]]
    unusual[0] = unusual[0][: peptide + 1]
    unusual[1][peptide + 1 :] = ["P1", "P2", "P3", "P4"]
    cases["unusual_protein_mapping"] = unusual

    results = []
    for name, rows in cases.items():
        root = args.output / name
        root.mkdir()
        pin = root / "input.pin"
        write_pin(pin, header, rows)
        target_out = root / "target.psms.tsv"
        decoy_out = root / "decoy.psms.tsv"
        command = [
            str(args.rust.resolve()), "--canonical", "--no-select-c", "--seed", "1",
            "--results-psms", str(target_out), "--decoy-results-psms", str(decoy_out), str(pin),
        ]
        started = time.time()
        with (root / "stdout.log").open("w") as stdout, (root / "stderr.log").open("w") as stderr:
            execution = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
        record = {
            "case": name, "input_rows": len(rows), "input_sha256": sha256(pin),
            "command": command, "shell_display": shlex.join(command), "exit_code": execution.returncode,
            "wall_seconds": time.time() - started,
        }
        if execution.returncode == 0:
            record["checks"] = output_checks(target_out, decoy_out)
        else:
            record["stderr_tail"] = (root / "stderr.log").read_text(errors="replace").splitlines()[-20:]
        results.append(record)
        print(name, execution.returncode, record.get("checks", record.get("stderr_tail", [])[-1:]))
    manifest = {
        "schema_version": 1,
        "study": "adversarial failure and edge-case behavior",
        "fixture": {"path": str(args.fixture.resolve()), "sha256": sha256(args.fixture)},
        "runner": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        "results": results,
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
