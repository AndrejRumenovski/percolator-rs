#!/usr/bin/env python3
"""Render machine-readable TSV tables from the final JSON analyses."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


CONDITIONS = (
    "original", "homology_depleted", "size_control_130363",
    "size_control_155921", "size_control_196613",
)


def write_rows(path: Path, rows: list[dict]) -> None:
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment_root
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    result = json.loads((root / "analysis/rescored_results.json").read_text())
    raw = json.loads((root / "analysis/raw_xcorr.json").read_text())
    pre = json.loads((root / "databases/presearch_summary.json").read_text())

    write_rows(tables / "primary_endpoints.tsv", [
        {"condition": condition, **result["endpoints"][condition]} for condition in CONDITIONS
    ])
    write_rows(tables / "q_thresholds.tsv", [
        {"condition": condition, **row}
        for condition in CONDITIONS for row in result["q_thresholds"][condition]
    ])
    write_rows(tables / "pep_calibration.tsv", [
        {"condition": condition, **row}
        for condition in CONDITIONS for row in result["pep_calibration"][condition]["bins"]
    ])
    write_rows(tables / "training_dose_response.tsv", [
        {"condition": condition, **row}
        for condition in CONDITIONS for row in result["training_dose_response"][condition]
    ])
    write_rows(tables / "seed_reproducibility.tsv", [
        {"condition": condition, **row}
        for condition, rows in result["seed_reproducibility"].items() for row in rows
    ])
    write_rows(tables / "enzN_enzC_interaction.tsv", [
        {"condition": condition, "feature_set": feature_set, **row}
        for condition, comparison in result["enzN_enzC_interaction"].items()
        for feature_set, row in comparison.items()
    ])
    write_rows(tables / "raw_xcorr_matched_depth.tsv", [
        {"condition": condition, "entrapment_decoy_depth": int(depth), **row}
        for condition in CONDITIONS
        for depth, row in raw["conditions"][condition]["matched_entrapment_decoy_depths"].items()
        if row is not None
    ])
    write_rows(tables / "presearch_opportunity.tsv", [
        {"condition": condition,
         **pre["conditions"][condition]["database"],
         **{f"entrapment_{key}": value for key, value in
            pre["conditions"][condition]["entrapment_component"].items()
            if not isinstance(value, list)},
         **{f"retained_{key}": value for key, value in
            pre["conditions"][condition]["entrapment_opportunity_retained"].items()},
         "primary_near_homolog_proteins": pre["conditions"][condition]["similarity_to_native"]["proteins_with_primary_near_homolog"],
         "exact_shared_tryptic_sequences": pre["conditions"][condition]["exact_shared_full_tryptic_unique_sequences"]}
        for condition in CONDITIONS
    ])
    print(f"wrote {len(list(tables.glob('*.tsv')))} tables to {tables}")


if __name__ == "__main__":
    main()
