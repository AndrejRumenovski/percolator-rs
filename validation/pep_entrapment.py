#!/usr/bin/env python3
"""Entrapment-based aggregate calibration check for reported PSM PEPs."""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import math
from pathlib import Path


EDGES = (0.0, 1e-12, 1e-6, 1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0000001)


def proteins(row: dict, decoy: bool) -> list[str]:
    fields = []
    for key, value in row.items():
        if key == "proteinIds" or key is None:
            fields.extend(value if isinstance(value, list) else [value])
    result = [protein for value in fields if value for protein in value.replace(";", "\t").split("\t") if protein]
    if decoy:
        result = [protein.removeprefix("DECOY_").removeprefix("decoy_") for protein in result]
    return result


def load(pattern: str, decoy: bool) -> tuple[list[dict], list[dict]]:
    rows = []
    files = []
    for name in sorted(glob.glob(pattern)):
        path = Path(name)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        count = 0
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                pep = float(row["posterior_error_prob"])
                if not math.isfinite(pep) or not 0 <= pep <= 1:
                    raise ValueError(f"{path}: invalid PEP {pep}")
                members = proteins(row, decoy)
                pure = bool(members) and all(member.startswith("ENT_") for member in members)
                mixed = any(member.startswith("ENT_") for member in members) and not pure
                rows.append({"pep": pep, "pure_entrapment": pure, "mixed": mixed})
                count += 1
        files.append({"path": str(path), "sha256": digest, "rows": count})
    if not rows:
        raise ValueError(f"no rows match {pattern}")
    return rows, files


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--decoys", required=True)
    parser.add_argument("--fallback-entrapment-fraction", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    args = parser.parse_args()
    targets, target_files = load(args.targets, False)
    decoys, decoy_files = load(args.decoys, True)
    usable_decoys = [row for row in decoys if not row["mixed"]]
    global_fraction = (
        sum(row["pure_entrapment"] for row in usable_decoys) / len(usable_decoys)
        if usable_decoys else args.fallback_entrapment_fraction
    )

    bins = []
    for lower, upper in zip(EDGES, EDGES[1:]):
        selected = [row for row in targets if lower <= row["pep"] < upper]
        selected_decoys = [row for row in decoys if lower <= row["pep"] < upper and not row["mixed"]]
        fraction = (
            sum(row["pure_entrapment"] for row in selected_decoys) / len(selected_decoys)
            if selected_decoys else global_fraction
        )
        entrapment = sum(row["pure_entrapment"] for row in selected)
        estimated_false = entrapment / fraction if fraction else None
        adjusted_error = estimated_false / len(selected) if selected and estimated_false is not None else None
        bins.append({
            "method": args.method,
            "pep_lower_inclusive": lower,
            "pep_upper_exclusive": upper,
            "targets": len(selected),
            "mean_reported_pep": sum(row["pep"] for row in selected) / len(selected) if selected else None,
            "pure_entrapment_targets": entrapment,
            "mixed_targets": sum(row["mixed"] for row in selected),
            "decoys_for_fraction": len(selected_decoys),
            "effective_entrapment_fraction": fraction,
            "adjusted_observed_error": adjusted_error,
            "calibration_error_observed_minus_reported": (
                adjusted_error - sum(row["pep"] for row in selected) / len(selected)
                if selected and adjusted_error is not None else None
            ),
        })
    nonempty = [row for row in bins if row["targets"] and row["adjusted_observed_error"] is not None]
    total = sum(row["targets"] for row in nonempty)
    ece = sum(row["targets"] * abs(row["calibration_error_observed_minus_reported"]) for row in nonempty) / total
    signed = sum(row["targets"] * row["calibration_error_observed_minus_reported"] for row in nonempty) / total
    result = {
        "schema_version": 1,
        "analysis": "entrapment-adjusted aggregate PEP calibration",
        "evaluator": {
            "path": str(Path(__file__).resolve()),
            "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "method": args.method,
        "interpretation": (
            "Entrapment matches are known false, but native matches are not individually labelled. "
            "Observed error is therefore a bin-level estimate obtained by dividing entrapment matches "
            "by the empirical foreign fraction among non-mixed decoys in the same PEP bin. It is not "
            "individual ground truth and bins with few decoys are unstable."
        ),
        "target_files": target_files,
        "decoy_files": decoy_files,
        "global_effective_entrapment_fraction": global_fraction,
        "targets": len(targets),
        "exact_zero_peps": sum(row["pep"] == 0 for row in targets),
        "pure_entrapment_with_exact_zero_pep": sum(row["pep"] == 0 and row["pure_entrapment"] for row in targets),
        "weighted_absolute_calibration_error": ece,
        "weighted_signed_observed_minus_reported": signed,
        "bins": bins,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    with args.table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=bins[0].keys())
        writer.writeheader()
        writer.writerows(bins)
    print(json.dumps({key: result[key] for key in (
        "method", "targets", "exact_zero_peps", "pure_entrapment_with_exact_zero_pep",
        "weighted_absolute_calibration_error", "weighted_signed_observed_minus_reported",
    )}, indent=2))


if __name__ == "__main__":
    main()
