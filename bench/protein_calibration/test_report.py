#!/usr/bin/env python3
"""Focused tests for protein-calibration labels and metrics."""

from __future__ import annotations

import math
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


if __name__ == "__main__":
    unittest.main()
