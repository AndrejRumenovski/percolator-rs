#!/usr/bin/env python3
"""Compare reported q-values with foreign-proteome entrapment errors."""

import argparse
import csv
import math
from pathlib import Path


def load(path: Path, decoy=False):
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            q = float(row["q-value"])
            protein_fields = []
            for key, value in row.items():
                if key != "proteinIds" and key is not None:
                    continue
                protein_fields.extend(value if isinstance(value, list) else [value])
            protein_text = "\t".join(value for value in protein_fields if value)
            proteins = [p for p in protein_text.replace(";", "\t").split("\t") if p]
            if decoy:
                proteins = [
                    p.removeprefix("DECOY_").removeprefix("decoy_") for p in proteins
                ]
            pure_entrapment = bool(proteins) and all(p.startswith("ENT_") for p in proteins)
            mixed = any(p.startswith("ENT_") for p in proteins) and not pure_entrapment
            rows.append((q, pure_entrapment, mixed))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--entrapment-fraction", required=True, type=float)
    parser.add_argument("--output", type=Path)
    parser.add_argument("inputs", nargs="+", type=Path)
    args = parser.parse_args()
    thresholds = (0.001, 0.005, 0.01, 0.02, 0.05, 0.1)
    grouped = {}
    grouped_decoys = {}
    for path in args.inputs:
        method = path.name.split(".")[0]
        grouped.setdefault(method, []).extend(load(path))
        decoy_path = Path(str(path).replace(".target.psms.tsv", ".decoy.psms.tsv"))
        if decoy_path.exists():
            grouped_decoys.setdefault(method, []).extend(load(decoy_path, decoy=True))
    output_rows = []

    for method, rows in sorted(grouped.items()):
        for threshold in thresholds:
            accepted = [row for row in rows if row[0] <= threshold]
            entrapment = sum(row[1] for row in accepted)
            mixed = sum(row[2] for row in accepted)
            accepted_decoys = [row for row in grouped_decoys.get(method, []) if row[0] <= threshold]
            pure_decoys = [row for row in accepted_decoys if not row[2]]
            entrapment_decoys = sum(row[1] for row in pure_decoys)
            effective_fraction = (
                entrapment_decoys / len(pure_decoys) if pure_decoys else args.entrapment_fraction
            )
            estimated_false = entrapment / effective_fraction
            observed_fdp = estimated_false / len(accepted) if accepted else 0.0
            # Wilson interval for the observed entrapment proportion, scaled by
            # its fraction of target search space to obtain an FDP interval.
            n = len(accepted)
            if n:
                z = 1.959963984540054
                p = entrapment / n
                center = (p + z * z / (2 * n)) / (1 + z * z / n)
                half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
                low = max(0.0, center - half) / effective_fraction
                high = min(1.0, center + half) / effective_fraction
            else:
                low = high = 0.0
            output_rows.append(
                (
                    method, threshold, n, entrapment, mixed, effective_fraction,
                    observed_fdp, low, high,
                )
            )

    header = (
        "method", "q_threshold", "accepted_targets", "pure_entrapment", "mixed",
        "effective_entrapment_fraction", "adjusted_fdp", "fdp_ci95_low", "fdp_ci95_high",
    )
    print("\t".join(header))
    for method, threshold, accepted, entrapment, mixed, fraction, fdp, low, high in output_rows:
        print(
            f"{method}\t{threshold:.3g}\t{accepted}\t{entrapment}\t{mixed}"
            f"\t{fraction:.6f}\t{fdp:.6f}\t{low:.6f}\t{high:.6f}"
        )

    if args.output:
        with args.output.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(output_rows)


if __name__ == "__main__":
    main()
