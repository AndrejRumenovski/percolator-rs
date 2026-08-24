#!/usr/bin/env python3
"""Evaluate PrEST protein groups against explicit present/absent ground truth."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.1)
METHODS = ("picked", "bayes-fixed", "bayes-selected")


@dataclass(frozen=True)
class Group:
    q: float
    pep: float
    score: float
    proteins: tuple[str, ...]
    error: int
    composition: str


def read_truth(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        truth = {row["protein_id"]: row["pool"] for row in rows}
    expected = {"A": 192, "B": 191, "RANDOM": 1000}
    observed = {pool: sum(value == pool for value in truth.values()) for pool in expected}
    if observed != expected:
        raise ValueError(f"ground-truth counts drifted: {observed}, expected {expected}")
    return truth


def present_pools(vial: str) -> set[str]:
    return {
        "A": {"A"},
        "B": {"B"},
        "AB": {"A", "B"},
        "BLANK": set(),
    }[vial.upper()]


def classify(proteins: tuple[str, ...], vial: str, truth: dict[str, str]) -> tuple[int, str]:
    pools: list[str] = []
    for protein in proteins:
        target = protein.removeprefix("DECOY_")
        if target not in truth:
            raise ValueError(f"unknown PrEST protein identifier: {protein}")
        pools.append(truth[target])
    present = present_pools(vial)
    has_present = any(pool in present for pool in pools)
    has_absent = any(pool not in present for pool in pools)
    if has_present and has_absent:
        composition = "mixed_present_absent"
    elif has_present:
        composition = "pure_present"
    elif pools and all(pool == "RANDOM" for pool in pools):
        composition = "pure_random_entrapment"
    elif pools and all(pool in {"A", "B"} for pool in pools):
        composition = "pure_absent_paired_pool"
    else:
        composition = "mixed_absent_sources"
    return (0 if has_present else 1), composition


def load_groups(path: Path, vial: str, truth: dict[str, str]) -> list[Group]:
    groups: list[Group] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"q-value", "posterior_error_prob", "score", "proteinIds"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"missing required protein columns in {path}")
        for row in reader:
            proteins = tuple(value for value in row["proteinIds"].split(",") if value)
            if not proteins:
                raise ValueError(f"empty protein group in {path}")
            error, composition = classify(proteins, vial, truth)
            groups.append(
                Group(
                    q=float(row["q-value"]),
                    pep=float(row["posterior_error_prob"]),
                    score=float(row["score"]),
                    proteins=proteins,
                    error=error,
                    composition=composition,
                )
            )
    return groups


def absent_fraction(vial: str, truth: dict[str, str]) -> float:
    present = present_pools(vial)
    return sum(pool not in present for pool in truth.values()) / len(truth)


def threshold_metrics(groups: list[Group], threshold: float, fraction: float) -> dict[str, float]:
    accepted = [group for group in groups if group.q <= threshold]
    false = sum(group.error for group in accepted)
    raw_fdp = false / len(accepted) if accepted else 0.0
    adjusted_fdp = false / (fraction * len(accepted)) if accepted and fraction else 0.0
    if accepted and fraction:
        z = 1.959963984540054
        n = len(accepted)
        proportion = false / n
        center = (proportion + z * z / (2 * n)) / (1 + z * z / n)
        half = z * math.sqrt(
            proportion * (1 - proportion) / n + z * z / (4 * n * n)
        ) / (1 + z * z / n)
        ci_low = max(0.0, center - half) / fraction
        ci_high = min(1.0, center + half) / fraction
    else:
        ci_low = ci_high = 0.0
    return {
        "accepted": len(accepted),
        "true": len(accepted) - false,
        "false": false,
        "raw_fdp": raw_fdp,
        "adjusted_fdp": adjusted_fdp,
        "adjusted_ci95_low": ci_low,
        "adjusted_ci95_high": ci_high,
    }


def auc(groups: list[Group]) -> float:
    """Tie-aware ROC AUC, treating present groups as positives."""
    positives = sum(group.error == 0 for group in groups)
    negatives = len(groups) - positives
    if not positives or not negatives:
        return math.nan
    ordered = sorted(groups, key=lambda group: group.score)
    rank_sum = 0.0
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end].score == ordered[index].score:
            end += 1
        average_rank = ((index + 1) + end) / 2.0
        rank_sum += average_rank * sum(group.error == 0 for group in ordered[index:end])
        index = end
    return (rank_sum - positives * (positives + 1) / 2) / (positives * negatives)


def partial_auc(groups: list[Group], max_fpr: float = 0.05) -> float:
    """Normalized trapezoidal ROC area between FPR 0 and max_fpr."""
    positives = sum(group.error == 0 for group in groups)
    negatives = len(groups) - positives
    if not positives or not negatives:
        return math.nan
    ordered = sorted(groups, key=lambda group: group.score, reverse=True)
    points = [(0.0, 0.0)]
    tp = fp = index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end].score == ordered[index].score:
            end += 1
        tp += sum(group.error == 0 for group in ordered[index:end])
        fp += sum(group.error == 1 for group in ordered[index:end])
        points.append((fp / negatives, tp / positives))
        index = end
    area = 0.0
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 >= max_fpr:
            break
        stop = min(x1, max_fpr)
        if x1 == x0:
            y_stop = y1
        else:
            y_stop = y0 + (y1 - y0) * (stop - x0) / (x1 - x0)
        area += (stop - x0) * (y0 + y_stop) / 2
        if x1 >= max_fpr:
            break
    return area / max_fpr


def probability_metrics(groups: list[Group], bins: int = 10) -> tuple[float, float]:
    if not groups:
        return math.nan, math.nan
    brier = sum((group.pep - group.error) ** 2 for group in groups) / len(groups)
    buckets: list[list[Group]] = [[] for _ in range(bins)]
    for group in groups:
        index = min(int(max(0.0, min(1.0, group.pep)) * bins), bins - 1)
        buckets[index].append(group)
    ece = 0.0
    for bucket in buckets:
        if not bucket:
            continue
        predicted = sum(group.pep for group in bucket) / len(bucket)
        observed = sum(group.error for group in bucket) / len(bucket)
        ece += len(bucket) / len(groups) * abs(predicted - observed)
    return brier, ece


def format_number(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    if math.isnan(value):
        return "NA"
    return f"{value:.8g}"


def write_tsv(path: Path, header: tuple[str, ...], rows: list[tuple[object, ...]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(format_number(value) if isinstance(value, (float, int)) else value for value in row)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    truth = read_truth(args.truth)
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    threshold_rows: list[tuple[object, ...]] = []
    summary_rows: list[tuple[object, ...]] = []
    composition_rows: list[tuple[object, ...]] = []
    aggregates: dict[tuple[str, str], list[Group]] = defaultdict(list)
    aggregate_times: dict[tuple[str, str], list[tuple[float, int]]] = defaultdict(list)

    def report_one(
        sample: str,
        vial: str,
        split: str,
        method: str,
        groups: list[Group],
        wall: float,
        rss: int,
    ) -> None:
        fraction = absent_fraction(vial, truth) if sample != "ALL" else math.nan
        metrics_at_threshold: list[dict[str, float]] = []
        if sample != "ALL":
            for threshold in THRESHOLDS:
                metrics = threshold_metrics(groups, threshold, fraction)
                metrics_at_threshold.append(metrics)
                threshold_rows.append(
                    (
                        sample, vial, split, method, threshold, metrics["accepted"],
                        metrics["true"], metrics["false"], fraction,
                        metrics["raw_fdp"], metrics["adjusted_fdp"],
                        metrics["adjusted_ci95_low"], metrics["adjusted_ci95_high"],
                    )
                )
        else:
            # Aggregate rows use their already-computed error labels and report the
            # unadjusted known-absence FDP. Per-vial adjusted values remain above.
            for threshold in THRESHOLDS:
                metrics = threshold_metrics(groups, threshold, 1.0)
                metrics_at_threshold.append(metrics)
                threshold_rows.append(
                    (
                        sample, vial, split, method, threshold, metrics["accepted"],
                        metrics["true"], metrics["false"], math.nan,
                        metrics["raw_fdp"], math.nan, math.nan, math.nan,
                    )
                )
        # Picked inference exposes its best-peptide PEP in the shared output
        # schema; that value is not a protein-group posterior and must not be
        # evaluated as one. Threshold calibration remains meaningful.
        brier, ece = probability_metrics(groups) if method != "picked" else (math.nan, math.nan)
        calibration_values = [
            abs(metric["adjusted_fdp"] - threshold)
            for metric, threshold in zip(metrics_at_threshold, THRESHOLDS)
            if metric["accepted"] and sample != "ALL"
        ]
        calibration_mae = (
            sum(calibration_values) / len(calibration_values) if calibration_values else math.nan
        )
        q01 = metrics_at_threshold[2]
        summary_rows.append(
            (
                sample, vial, split, method, len(groups), q01["accepted"], q01["false"],
                auc(groups), partial_auc(groups), brier, ece, calibration_mae, wall, rss,
            )
        )
        accepted_q01 = [group for group in groups if group.q <= 0.01]
        counts = defaultdict(int)
        for group in accepted_q01:
            counts[group.composition] += 1
        for composition in (
            "pure_present", "mixed_present_absent", "pure_absent_paired_pool",
            "pure_random_entrapment", "mixed_absent_sources",
        ):
            composition_rows.append(
                (sample, vial, split, method, 0.01, composition, counts[composition])
            )

    for item in manifest:
        sample, vial, split = item["sample"], item["vial"], item["split"]
        for method in METHODS:
            run = args.runs / sample / method
            target = run / "target.tsv"
            timing = run / "time.tsv"
            if not target.exists() or not timing.exists():
                raise FileNotFoundError(f"missing completed run for {sample}/{method}")
            groups = load_groups(target, vial, truth)
            wall_text, rss_text = timing.read_text(encoding="utf-8").strip().split("\t")
            wall, rss = float(wall_text), int(rss_text)
            report_one(sample, vial, split, method, groups, wall, rss)
            aggregates[(split, method)].extend(groups)
            aggregate_times[(split, method)].append((wall, rss))

    for (split, method), groups in sorted(aggregates.items()):
        timings = aggregate_times[(split, method)]
        report_one(
            "ALL", "ALL", split, method, groups,
            sum(wall for wall, _rss in timings), max(rss for _wall, rss in timings),
        )

    write_tsv(
        args.output_dir / "thresholds.tsv",
        (
            "sample", "vial", "split", "method", "q_threshold", "accepted", "true",
            "known_absent", "absent_database_fraction", "entrapment_fdp_lower_bound",
            "search_space_adjusted_fdp", "adjusted_fdp_ci95_low", "adjusted_fdp_ci95_high",
        ),
        threshold_rows,
    )
    write_tsv(
        args.output_dir / "summary.tsv",
        (
            "sample", "vial", "split", "method", "protein_groups", "accepted_q_le_0.01",
            "known_absent_q_le_0.01", "roc_auc", "partial_auc_fpr_le_0.05", "brier",
            "ece_10_bin", "threshold_calibration_mae", "wall_seconds", "peak_rss_kb",
        ),
        summary_rows,
    )
    write_tsv(
        args.output_dir / "composition.tsv",
        ("sample", "vial", "split", "method", "q_threshold", "composition", "groups"),
        composition_rows,
    )
    print(f"wrote {args.output_dir / 'thresholds.tsv'}")
    print(f"wrote {args.output_dir / 'summary.tsv'}")
    print(f"wrote {args.output_dir / 'composition.tsv'}")


if __name__ == "__main__":
    main()
