import csv
import tempfile
import unittest
from pathlib import Path

from validation import psm_agreement


HEADER = ("PSMId", "score", "q-value", "posterior_error_prob", "peptide", "proteinIds")


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)


class AgreementTest(unittest.TestCase):
    def test_threshold_overlap_and_qualified_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rust = root / "rust"
            cpp = root / "cpp"
            write(rust / "a" / "target.psms.tsv", [
                ("same", 3, 0.001, 0.01, "A.PEP.K", "P1"),
                ("rust", 2, 0.009, 0.02, "A.RUST.K", "P2"),
                ("cpp", 1, 0.02, 0.10, "A.CPP.K", "P3"),
            ])
            write(cpp / "a" / "target.psms.tsv", [
                ("same", 30, 0.002, 0.02, "A.PEP.K", "P1"),
                ("rust", 10, 0.02, 0.20, "A.RUST.K", "P2"),
                ("cpp", 20, 0.009, 0.05, "A.CPP.K", "P3"),
            ])
            result = psm_agreement.compare(rust, cpp)
            q01 = next(row for row in result["thresholds"] if row["q_threshold"] == 0.01)
            self.assertEqual(q01["rust"], 2)
            self.assertEqual(q01["cpp"], 2)
            self.assertEqual(q01["intersection"], 1)
            self.assertEqual(q01["rust_only"], 1)
            self.assertEqual(q01["cpp_only"], 1)
            self.assertAlmostEqual(q01["jaccard"], 1 / 3)
            self.assertEqual(result["row_counts"]["matching_unambiguous"], 3)

    def test_rejects_nonfinite_statistics(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write(root / "rust.tsv", [("x", "nan", 0, 0, "P", "P")])
            write(root / "cpp.tsv", [("x", 1, 0, 0, "P", "P")])
            with self.assertRaisesRegex(ValueError, "non-finite score"):
                psm_agreement.compare(root / "rust.tsv", root / "cpp.tsv")


if __name__ == "__main__":
    unittest.main()
