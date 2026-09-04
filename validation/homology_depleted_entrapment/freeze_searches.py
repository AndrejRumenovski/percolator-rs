#!/usr/bin/env python3
"""Freeze hashes and labeling checks for all from-raw Comet searches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path


CONDITIONS = (
    "original", "homology_depleted", "size_control_130363",
    "size_control_155921", "size_control_196613",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_pin_sha256(path: Path) -> str:
    """Ignore only the output-directory prefix embedded in SpecId."""
    digest = hashlib.sha256()
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row_number, row in enumerate(reader):
            if row_number and row:
                row[0] = Path(row[0]).name
            digest.update(("\t".join(row) + "\n").encode())
    return digest.hexdigest()


def compare_pins(old: Path, new: Path) -> dict:
    with old.open(newline="") as old_handle, new.open(newline="") as new_handle:
        left = csv.reader(old_handle, delimiter="\t")
        right = csv.reader(new_handle, delimiter="\t")
        old_header = next(left); new_header = next(right)
        if old_header != new_header:
            return {"headers_identical": False}
        field_differences = {name: 0 for name in old_header}
        differing_rows = 0
        old_extra = new_extra = 0
        rows = 0
        for old_row, new_row in itertools.zip_longest(left, right):
            if old_row is None:
                new_extra += 1
                continue
            if new_row is None:
                old_extra += 1
                continue
            rows += 1
            old_row[0] = Path(old_row[0]).name
            new_row[0] = Path(new_row[0]).name
            if old_row != new_row:
                differing_rows += 1
                for index, (old_value, new_value) in enumerate(zip(old_row, new_row)):
                    if old_value != new_value:
                        field_differences[old_header[index]] += 1
        return {
            "headers_identical": True, "rows_compared": rows,
            "differing_rows": differing_rows, "old_extra_rows": old_extra,
            "new_extra_rows": new_extra,
            "field_differences": {name: count for name, count in field_differences.items() if count},
            "psm_identity_columns_unchanged": not any(
                field_differences.get(name, 0) for name in ("Label", "ScanNr", "ExpMass", "Peptide", "Proteins")
            ),
        }


def label_check(path: Path) -> dict:
    counts = {"rows": 0, "targets": 0, "decoys": 0,
              "target_rows_without_target_protein": 0,
              "target_rows_with_mixed_target_decoy_mapping": 0,
              "decoy_rows_without_decoy_protein": 0,
              "entrapment_target_rows": 0,
              "entrapment_decoy_rows": 0}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        label_index = header.index("Label")
        protein_index = header.index("Proteins")
        for row in reader:
            counts["rows"] += 1
            decoy = row[label_index] == "-1"
            counts["decoys" if decoy else "targets"] += 1
            proteins = [value for field in row[protein_index:]
                        for value in field.replace(";", "\t").split("\t") if value]
            all_decoy = bool(proteins) and all(value.startswith("DECOY_") for value in proteins)
            has_target = any(not value.startswith("DECOY_") for value in proteins)
            has_decoy = any(value.startswith("DECOY_") for value in proteins)
            if decoy and not all_decoy:
                counts["decoy_rows_without_decoy_protein"] += 1
            if not decoy and not has_target:
                counts["target_rows_without_target_protein"] += 1
            if not decoy and has_target and has_decoy:
                counts["target_rows_with_mixed_target_decoy_mapping"] += 1
            stripped = [value.removeprefix("DECOY_") if decoy else value for value in proteins]
            pure_entrapment = bool(stripped) and all(value.startswith("ENT_") for value in stripped)
            if pure_entrapment:
                counts["entrapment_decoy_rows" if decoy else "entrapment_target_rows"] += 1
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-root", type=Path, required=True)
    parser.add_argument("--canonical-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    searches = []
    original_reproduction = []
    for condition in CONDITIONS:
        directories = sorted((args.experiment_root / "searches" / condition).glob("comet-*"))
        if len(directories) != 6:
            raise ValueError(f"{condition}: expected six search directories, found {len(directories)}")
        for directory in directories:
            pin = directory / "comet.pin"
            log = directory / "comet.log.txt"
            params = directory / "comet.params.txt"
            pepxml = directory / "comet.target.pep.xml"
            for required in (pin, log, params, pepxml):
                if not required.exists() or required.stat().st_size == 0:
                    raise FileNotFoundError(required)
            checks = label_check(pin)
            if checks["target_rows_without_target_protein"] or checks["decoy_rows_without_decoy_protein"]:
                raise ValueError(f"label integrity failed for {pin}: {checks}")
            row = {
                "condition": condition, "run": directory.name,
                "pin": {"path": str(pin), "sha256": sha256(pin),
                        "normalized_specid_sha256": normalized_pin_sha256(pin),
                        "bytes": pin.stat().st_size},
                "pepxml": {"path": str(pepxml), "sha256": sha256(pepxml),
                           "bytes": pepxml.stat().st_size},
                "log": {"path": str(log), "sha256": sha256(log)},
                "resolved_parameters": {"path": str(params), "sha256": sha256(params)},
                "time": (directory / "time.txt").read_text(),
                "label_check": checks,
            }
            searches.append(row)
            if condition == "original":
                old = args.canonical_root / directory.name / "comet.pin"
                original_reproduction.append({
                    "run": directory.name,
                    "old_pin": str(old), "old_sha256": sha256(old),
                    "old_normalized_specid_sha256": normalized_pin_sha256(old),
                    "new_normalized_specid_sha256": row["pin"]["normalized_specid_sha256"],
                    "identical_after_specid_path_normalization": (
                        normalized_pin_sha256(old) == row["pin"]["normalized_specid_sha256"]
                    ),
                    "field_comparison": compare_pins(old, pin),
                })
    databases = {}
    for condition in CONDITIONS:
        path = args.experiment_root / "databases" / f"{condition}.fasta"
        databases[condition] = {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}
    result = {
        "schema_version": 1,
        "from_raw_searches": 30,
        "decoy_generation": "Comet decoy_search=1; enzyme-aware target-peptide reversal with terminal residue retained",
        "databases": databases,
        "searches": searches,
        "original_vs_prior_pin_reproduction": original_reproduction,
        "no_post_search_filtering": True,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "searches": len(searches),
        "label_checks_passed": all(not row["label_check"]["target_rows_without_target_protein"] and
                                   not row["label_check"]["decoy_rows_without_decoy_protein"]
                                   for row in searches),
        "original_pins_reproduced_after_path_normalization": [
            row["identical_after_specid_path_normalization"] for row in original_reproduction
        ],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
