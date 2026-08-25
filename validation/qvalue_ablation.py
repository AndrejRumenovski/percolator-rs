#!/usr/bin/env python3
"""Post-hoc q-value estimator ablations on fixed model scores."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.10)
VARIANTS = (
    ("D_over_T_rowwise", 0, False, False),
    ("D_plus_1_over_T_rowwise", 1, False, False),
    ("D_over_T_tie_grouped", 0, True, False),
    ("D_plus_1_over_T_tie_grouped", 1, True, False),
    ("rust_count_pi0_D_over_T_tie_grouped", 0, True, True),
    ("rust_count_pi0_D_plus_1_over_T_tie_grouped", 1, True, True),
)


def result_path(argv: list[str], option: str) -> Path:
    index = argv.index(option)
    return Path(argv[index + 1])


def load(target_path: Path, decoy_path: Path) -> list[tuple[float, bool]]:
    rows = []
    for path, target in ((target_path, True), (decoy_path, False)):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                score = float(row["score"])
                if not math.isfinite(score):
                    raise ValueError(f"{path}: non-finite score")
                rows.append((score, target))
    return rows


def qvalues(rows: list[tuple[float, bool]], pseudocount: int, tie_grouped: bool, count_pi0: bool) -> list[tuple[bool, float]]:
    ordered = sorted(rows, key=lambda row: row[0], reverse=True)
    total_targets = sum(target for _, target in ordered)
    total_decoys = len(ordered) - total_targets
    pi0 = min(1.0, total_decoys / total_targets) if count_pi0 and total_targets else 1.0
    targets = 0
    decoys = pseudocount
    raw = []
    if tie_grouped:
        start = 0
        while start < len(ordered):
            end = start + 1
            while end < len(ordered) and ordered[end][0] == ordered[start][0]:
                end += 1
            targets += sum(target for _, target in ordered[start:end])
            decoys += sum(not target for _, target in ordered[start:end])
            value = min(1.0, pi0 * decoys / max(1, targets))
            raw.extend([value] * (end - start))
            start = end
    else:
        for _, target in ordered:
            targets += target
            decoys += not target
            raw.append(min(1.0, pi0 * decoys / max(1, targets)))
    running = 1.0
    q = [1.0] * len(raw)
    for index in range(len(raw) - 1, -1, -1):
        running = min(running, raw[index])
        q[index] = running
    return [(ordered[index][1], q[index]) for index in range(len(q))]


def counts(scored: list[tuple[bool, float]]) -> dict:
    return {f"q_lt_{threshold:g}": sum(target and q < threshold for target, q in scored) for threshold in THRESHOLDS}


def describe(values: list[int]) -> dict:
    return {
        "mean": statistics.fmean(values), "median": statistics.median(values),
        "sd": statistics.stdev(values) if len(values) > 1 else 0.0,
        "minimum": min(values), "maximum": max(values), "range": max(values) - min(values),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--null-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.null_manifest.read_text())
    records = []
    for run in source["runs"]:
        if run["status"] != "complete":
            continue
        argv = run["execution"]["argv"]
        rows = load(result_path(argv, "--results-psms"), result_path(argv, "--decoy-results-psms"))
        for name, pseudocount, tie_grouped, count_pi0 in VARIANTS:
            records.append({
                "input": run["input"], "method": run["method"],
                "relabel_seed": run["relabel_seed"], "variant": name,
                "pseudocount": pseudocount, "tie_grouped": tie_grouped,
                "count_ratio_pi0": count_pi0,
                "false_target_psms": counts(qvalues(rows, pseudocount, tie_grouped, count_pi0)),
            })
    aggregate = []
    groups = sorted({(row["method"], row["variant"]) for row in records})
    for method, variant in groups:
        selected = [row for row in records if row["method"] == method and row["variant"] == variant]
        for threshold in THRESHOLDS:
            key = f"q_lt_{threshold:g}"
            values = [row["false_target_psms"][key] for row in selected]
            aggregate.append({
                "method": method, "variant": variant, "q_threshold": threshold,
                "replicates": len(values),
                "replicates_with_any_false_discovery": sum(value > 0 for value in values),
                "empirical_fdr": sum(value > 0 for value in values) / len(values),
                "false_target_count": describe(values),
            })
    result = {
        "schema_version": 1,
        "analysis": "fixed-score q-value estimator ablation under the complete null",
        "evaluator": {
            "path": str(Path(__file__).resolve()),
            "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "source_manifest": str(args.null_manifest),
        "source_manifest_sha256": hashlib.sha256(args.null_manifest.read_bytes()).hexdigest(),
        "threshold_policy": "q < threshold",
        "limitation": (
            "This is a post-hoc ablation on scores printed to six decimal places, not on the "
            "full-precision internal score vector. Printed-score ties and their order can differ "
            "from internal ties. The ablation isolates estimator behavior on the preserved output "
            "but is not an exact rerun of either implementation's internal statistic."
        ),
        "variants": [name for name, *_ in VARIANTS],
        "runs": records,
        "aggregate": aggregate,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    for row in aggregate:
        if row["q_threshold"] == 0.01:
            print(row["method"], row["variant"], row["empirical_fdr"], row["false_target_count"])


if __name__ == "__main__":
    main()
