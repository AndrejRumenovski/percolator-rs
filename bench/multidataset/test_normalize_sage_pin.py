#!/usr/bin/env python3
"""Regression tests for deterministic Sage PIN normalization."""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
HEADER = (
    "SpecId\tLabel\tScanNr\tposterior_error\tPeptide\tProteins\n"
)
ROWS = [
    "91\t-1\trun title #20\t0.9\tDECOYSEQ\tDECOY_p\n",
    "17\t1\trun title #3\t0.1\tTARGETSEQ\ttarget_p\n",
]


class NormalizeSagePinTests(unittest.TestCase):
    def normalize(self, rows: list[str]) -> str:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pin"
            output = root / "normalized.pin"
            source.write_text(HEADER + "".join(rows), encoding="utf-8")
            subprocess.run(
                ["python3", str(HERE / "normalize_sage_pin.py"), source, output],
                check=True,
            )
            return output.read_text(encoding="utf-8")

    def test_output_is_independent_of_sage_row_order_and_spec_ids(self) -> None:
        first = self.normalize(ROWS)
        reversed_rows = [
            ROWS[1].replace("17", "900", 1),
            ROWS[0].replace("91", "901", 1),
        ]
        second = self.normalize(reversed_rows)
        self.assertEqual(first, second)
        self.assertEqual(
            first,
            "SpecId\tLabel\tScanNr\tPeptide\tProteins\n"
            "1\t1\t3\tTARGETSEQ\ttarget_p\n"
            "2\t-1\t20\tDECOYSEQ\tDECOY_p\n",
        )


if __name__ == "__main__":
    unittest.main()
