#!/usr/bin/env python3
"""Focused tests for protein-calibration labels and metrics."""

from __future__ import annotations

import math
from pathlib import Path
import tempfile
import unittest

import report


class ReportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.truth = {"a": "A", "b": "B", "r1": "RANDOM", "r2": "RANDOM"}

    def test_group_truth_uses_at_least_one_present_member(self) -> None:
        self.assertEqual(report.classify(("a",), "A", self.truth), (0, "pure_present"))
        self.assertEqual(
            report.classify(("a", "b"), "A", self.truth),
            (0, "mixed_present_absent"),
        )
        self.assertEqual(
            report.classify(("b",), "A", self.truth),
            (1, "pure_absent_paired_pool"),
        )
        self.assertEqual(
            report.classify(("DECOY_r1",), "A", self.truth),
            (1, "pure_random_entrapment"),
        )

    def test_perfect_ranking_has_unit_auc_and_partial_auc(self) -> None:
        group = report.Group
        groups = [
            group(0.001, 0.01, 4.0, ("a",), 0, "pure_present"),
            group(0.01, 0.1, 3.0, ("a",), 0, "pure_present"),
            group(0.1, 0.8, 2.0, ("r1",), 1, "pure_random_entrapment"),
            group(0.2, 0.9, 1.0, ("r2",), 1, "pure_random_entrapment"),
        ]
        self.assertEqual(report.auc(groups), 1.0)
        self.assertEqual(report.partial_auc(groups), 1.0)

    def test_tied_scores_get_half_credit(self) -> None:
        group = report.Group
        groups = [
            group(0.1, 0.5, 1.0, ("a",), 0, "pure_present"),
            group(0.1, 0.5, 1.0, ("r1",), 1, "pure_random_entrapment"),
        ]
        self.assertEqual(report.auc(groups), 0.5)

    def test_adjusted_fdp_and_probability_metrics(self) -> None:
        group = report.Group
        groups = [
            group(0.01, 0.1, 2.0, ("a",), 0, "pure_present"),
            group(0.01, 0.9, 1.0, ("r1",), 1, "pure_random_entrapment"),
        ]
        metrics = report.threshold_metrics(groups, 0.01, 0.5)
        self.assertEqual(metrics["accepted"], 2)
        self.assertEqual(metrics["false"], 1)
        self.assertEqual(metrics["raw_fdp"], 0.5)
        self.assertEqual(metrics["adjusted_fdp"], 1.0)
        brier, ece = report.probability_metrics(groups)
        self.assertTrue(math.isclose(brier, 0.01))
        self.assertTrue(math.isclose(ece, 0.1))

    def test_current_picked_output_keeps_missing_posteriors_missing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "picked.tsv"
            path.write_text(
                "ProteinGroupId\tq-value\tposterior_error_prob\tscore\tnumPeptides\tproteinIds\n"
                "group1\t0.01\tNA\t2.0\t1\ta\n"
                "group2\t0.01\tNA\t1.0\t1\tr1\n"
            )
            groups = report.load_groups(path, "A", self.truth)
        self.assertTrue(all(group.pep is None for group in groups))
        self.assertEqual(report.threshold_metrics(groups, 0.01, 0.5)["false"], 1)
        self.assertEqual(report.auc(groups), 1.0)
        self.assertTrue(all(math.isnan(value) for value in report.probability_metrics(groups)))

    def test_incomplete_posteriors_do_not_select_a_calibration_subset(self) -> None:
        groups = [
            report.Group(0.01, None, 2.0, ("a",), 0, "pure_present"),
            report.Group(0.01, 0.9, 1.0, ("r1",), 1, "pure_random_entrapment"),
        ]
        self.assertTrue(all(math.isnan(value) for value in report.probability_metrics(groups)))

    def test_invalid_numeric_output_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.tsv"
            for q, pep, score in [
                ("NaN", "0.1", "1"), ("1.1", "0.1", "1"),
                ("0.1", "inf", "1"), ("0.1", "-0.1", "1"),
                ("0.1", "0.1", "inf"),
            ]:
                with self.subTest(q=q, pep=pep, score=score):
                    path.write_text(
                        "q-value\tposterior_error_prob\tscore\tproteinIds\n"
                        f"{q}\t{pep}\t{score}\ta\n"
                    )
                    with self.assertRaises(ValueError):
                        report.load_groups(path, "A", self.truth)


if __name__ == "__main__":
    unittest.main()
