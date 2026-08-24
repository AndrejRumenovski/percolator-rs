#!/usr/bin/env python3
"""Remove circular Sage features and canonicalize metadata in a Sage PIN.

Sage's parallel output order and its order-derived ``SpecId`` are not stable
between identical searches.  Sorting complete normalized records and assigning
stable IDs makes generated PIN checksums and downstream CV folds reproducible.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


DROP = {
    "posterior_error",
    "ion_mobility",
    "predicted_mobility",
    "sqrt(delta_mobility)",
}
INTEGER = re.compile(r"^\d+$")
SCAN = re.compile(r"(?:scan=|#)(\d+)", re.IGNORECASE)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open(encoding="utf-8", newline="") as src, args.output.open(
        "w", encoding="utf-8", newline=""
    ) as dst:
        reader = csv.reader(src, delimiter="\t")
        header = next(reader)
        keep = [index for index, name in enumerate(header) if name not in DROP]
        scan_index = header.index("ScanNr")
        output_header = [header[index] for index in keep]
        spec_index = output_header.index("SpecId")
        output_scan_index = output_header.index("ScanNr")
        rows = []
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"row {row_number} has {len(row)} fields; expected {len(header)}"
                )
            if not INTEGER.fullmatch(row[scan_index]):
                match = SCAN.search(row[scan_index])
                if not match:
                    raise ValueError(
                        f"row {row_number} has no integer scan: {row[scan_index]!r}"
                    )
                row[scan_index] = match.group(1)
            rows.append([row[index] for index in keep])

        key_indices = [i for i in range(len(output_header)) if i != spec_index]
        rows.sort(
            key=lambda row: (
                int(row[output_scan_index]),
                *(row[i] for i in key_indices),
            )
        )

        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(output_header)
        for stable_id, row in enumerate(rows, start=1):
            row[spec_index] = str(stable_id)
            writer.writerow(row)


if __name__ == "__main__":
    main()
