#!/usr/bin/env python3
"""Normalize Sage PIN metadata for a fair Rust/C++ Percolator comparison.

Sage preserves an MGF TITLE string in ScanNr, whereas C++ Percolator requires an
integer.  Extract the trailing ``#scan`` value.  Sage writes otherwise identical
results in nondeterministic order and assigns ``SpecId`` from that order, so sort
the normalized records and replace ``SpecId`` with a stable sequential ID.
Also remove Sage's already-trained posterior error feature to avoid circular
rescoring, along with constant mobility fields lost during MGF export.
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
SCAN = re.compile(r"#(\d+)(?:\D.*)?$")


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
        keep = [i for i, name in enumerate(header) if name not in DROP]
        scan_idx = header.index("ScanNr")
        output_header = [header[i] for i in keep]
        spec_idx = output_header.index("SpecId")
        output_scan_idx = output_header.index("ScanNr")
        rows = []
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"row {row_number} has {len(row)} fields; expected {len(header)}"
                )
            match = SCAN.search(row[scan_idx])
            if not match:
                raise ValueError(
                    f"row {row_number} has no integer scan suffix: {row[scan_idx]!r}"
                )
            row[scan_idx] = match.group(1)
            rows.append([row[i] for i in keep])

        # Sage's parallel writer does not promise record order.  Exclude its
        # order-derived SpecId from the key, then use the complete normalized
        # record as a deterministic tie breaker.
        key_indices = [i for i in range(len(output_header)) if i != spec_idx]
        rows.sort(
            key=lambda row: (
                int(row[output_scan_idx]),
                *(row[i] for i in key_indices),
            )
        )

        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(output_header)
        for stable_id, row in enumerate(rows, start=1):
            row[spec_idx] = str(stable_id)
            writer.writerow(row)


if __name__ == "__main__":
    main()
