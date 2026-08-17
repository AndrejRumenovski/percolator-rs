#!/usr/bin/env python3
"""Normalize Sage PIN metadata for a fair Rust/C++ Percolator comparison.

Sage preserves an MGF TITLE string in ScanNr, whereas C++ Percolator requires an
integer.  Extract the trailing ``#scan`` value (or use a deterministic row index).
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
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(header[i] for i in keep)
        for row_number, row in enumerate(reader, start=1):
            match = SCAN.search(row[scan_idx])
            row[scan_idx] = match.group(1) if match else str(row_number)
            writer.writerow(row[i] for i in keep)


if __name__ == "__main__":
    main()
