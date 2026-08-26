#!/usr/bin/env python3
"""Independent fail-closed probes for the PIN parser.

This is an audit harness, not part of the production implementation.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, default=Path("tests/fixtures/sample.pin"))
    args = parser.parse_args()

    with args.fixture.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = rows[0]
    label = header.index("Label")
    scan = header.index("ScanNr")
    feature = next(i for i, name in enumerate(header) if name not in {
        "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "Peptide", "Proteins"
    })

    cases = {
        "baseline": None,
        "label_zero": (label, "0"),
        "label_two": (label, "2"),
        "label_minus_two": (label, "-2"),
        "label_text": (label, "target"),
        "feature_nan": (feature, "NaN"),
        "feature_infinity": (feature, "inf"),
        "scan_text": (scan, "not-a-scan"),
        "extreme_finite_feature": (feature, "1.7976931348623157e308"),
    }
    output: dict[str, dict[str, object]] = {}
    with tempfile.TemporaryDirectory(prefix="percolator-parser-audit-") as temporary:
        root = Path(temporary)
        for name, mutation in cases.items():
            candidate = [row.copy() for row in rows]
            if mutation is not None:
                candidate[1][mutation[0]] = mutation[1]
            pin = root / f"{name}.pin"
            with pin.open("w", newline="") as handle:
                csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(candidate)
            command = [
                str(args.binary.resolve()), "--seed", "1", "--maxiter", "1",
                "--results-psms", str(root / f"{name}.target.tsv"),
                "--decoy-results-psms", str(root / f"{name}.decoy.tsv"), str(pin),
            ]
            run = subprocess.run(command, text=True, capture_output=True)
            output[name] = {
                "exit_code": run.returncode,
                "accepted_by_parser": "Reading" in run.stderr or run.returncode == 0,
                "last_diagnostic": (run.stderr + run.stdout).strip().splitlines()[-1]
                if (run.stderr + run.stdout).strip() else "",
            }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
