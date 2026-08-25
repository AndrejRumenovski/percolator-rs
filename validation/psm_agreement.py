#!/usr/bin/env python3
"""PSM-level agreement between percolator-rs and C++ Percolator outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter
from pathlib import Path
from typing import Iterable


THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.10)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def result_files(root: Path) -> list[Path]:
    if root.is_file():
        return [root]
    files = sorted(root.rglob("target.psms.tsv")) + sorted(root.rglob("decoy.psms.tsv"))
    if not files:
        files = sorted(root.rglob("*.psms.tsv"))
    return sorted(set(files))


def relative_key(root: Path, path: Path) -> str:
    if root.is_file():
        return path.name
    return str(path.relative_to(root))


def parse_float(value: str, field: str, path: Path, line: int) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError(f"{path}:{line}: invalid {field} {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"{path}:{line}: non-finite {field} {value!r}")
    return result


def load(root: Path) -> tuple[dict, dict, list[dict], dict]:
    grouped_rows: dict[tuple[str, str, str], list[dict]] = {}
    provenance = []
    for path in result_files(root):
        rel = relative_key(root, path)
        label = "decoy" if path.name.startswith("decoy") else "target"
        count = 0
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"PSMId", "score", "q-value", "posterior_error_prob"}
            missing = required - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{path}: missing columns {sorted(missing)}")
            for line, row in enumerate(reader, 2):
                key = (
                    rel.rsplit("/", 1)[0] if "/" in rel else ".",
                    label,
                    row["PSMId"],
                    row.get("peptide", ""),
                    row.get("proteinIds", ""),
                )
                parsed = {
                    "key": key,
                    "label": label,
                    "score": parse_float(row["score"], "score", path, line),
                    "q": parse_float(row["q-value"], "q-value", path, line),
                    "pep": parse_float(
                        row["posterior_error_prob"], "posterior_error_prob", path, line
                    ),
                    "peptide": row.get("peptide", ""),
                    "proteins": row.get("proteinIds", ""),
                }
                grouped_rows.setdefault(key, []).append(parsed)
                count += 1
        provenance.append({"path": str(path), "relative_path": rel, "sha256": sha256(path), "rows": count})
    if not grouped_rows:
        raise ValueError(f"no PSM result rows found under {root}")
    # Some legacy PINs reuse PSMId even for distinct feature rows. Output is score-sorted,
    # so there is no defensible occurrence-wise cross-tool matching key. Exclude every
    # ambiguous key rather than manufacturing a pairing.
    rows = {key: values[0] for key, values in grouped_rows.items() if len(values) == 1}
    ambiguity = {
        "ambiguous_qualified_keys_excluded": sum(len(values) > 1 for values in grouped_rows.values()),
        "ambiguous_rows_excluded": sum(len(values) for values in grouped_rows.values() if len(values) > 1),
    }
    return rows, grouped_rows, provenance, ambiguity


def pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) < 2:
        return None
    lm = statistics.fmean(left)
    rm = statistics.fmean(right)
    numerator = sum((x - lm) * (y - rm) for x, y in zip(left, right))
    ld = sum((x - lm) ** 2 for x in left)
    rd = sum((y - rm) ** 2 for y in right)
    return numerator / math.sqrt(ld * rd) if ld > 0 and rd > 0 else None


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def spearman(left: list[float], right: list[float]) -> float | None:
    return pearson(average_ranks(left), average_ranks(right))


def quantile(values: list[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def distribution(values: Iterable[float]) -> dict:
    data = list(values)
    if not data:
        return {"n": 0, "mean": None, "median": None, "sd": None, "min": None, "max": None,
                "p90": None, "p95": None, "p99": None}
    return {
        "n": len(data),
        "mean": statistics.fmean(data),
        "median": statistics.median(data),
        "sd": statistics.stdev(data) if len(data) > 1 else 0.0,
        "min": min(data),
        "max": max(data),
        "p90": quantile(data, 0.90),
        "p95": quantile(data, 0.95),
        "p99": quantile(data, 0.99),
    }


def selected(groups: dict, threshold: float) -> Counter:
    return Counter({
        key: sum(row["q"] < threshold for row in values)
        for key, values in groups.items()
        if key[1] == "target"
    })


def representative(key, rust, cpp) -> dict:
    return {
        "qualified_psm_id": list(key),
        "peptide": rust["peptide"] or cpp["peptide"],
        "rust": {name: rust[name] for name in ("score", "q", "pep")},
        "cpp": {name: cpp[name] for name in ("score", "q", "pep")},
        "score_difference": rust["score"] - cpp["score"],
        "q_difference": rust["q"] - cpp["q"],
        "pep_difference": rust["pep"] - cpp["pep"],
    }


def compare(rust_root: Path, cpp_root: Path) -> dict:
    rust, rust_groups, rust_files, rust_ambiguity = load(rust_root)
    cpp, cpp_groups, cpp_files, cpp_ambiguity = load(cpp_root)
    common = sorted(rust.keys() & cpp.keys())
    rust_only_rows = sorted(rust.keys() - cpp.keys())
    cpp_only_rows = sorted(cpp.keys() - rust.keys())
    if not common:
        raise ValueError("the Rust and C++ outputs have no matching qualified PSM IDs")

    rs = [rust[key]["score"] for key in common]
    cs = [cpp[key]["score"] for key in common]
    rq = [rust[key]["q"] for key in common]
    cq = [cpp[key]["q"] for key in common]
    rp = [rust[key]["pep"] for key in common]
    cp = [cpp[key]["pep"] for key in common]
    rust_score_ranks = average_ranks(rs)
    cpp_score_ranks = average_ranks(cs)
    rank_differences = [a - b for a, b in zip(rust_score_ranks, cpp_score_ranks)]
    rank_difference_by_key = dict(zip(common, rank_differences))

    thresholds = []
    exclusive_characterization = []
    for threshold in THRESHOLDS:
        rset = selected(rust_groups, threshold)
        cset = selected(cpp_groups, threshold)
        intersection = rset & cset
        union = rset | cset
        rust_count = sum(rset.values())
        cpp_count = sum(cset.values())
        intersection_count = sum(intersection.values())
        union_count = sum(union.values())
        thresholds.append({
            "q_threshold": threshold,
            "comparison": "strict_less_than",
            "rust": rust_count,
            "cpp": cpp_count,
            "intersection": intersection_count,
            "rust_only": rust_count - intersection_count,
            "cpp_only": cpp_count - intersection_count,
            "union": union_count,
            "jaccard": intersection_count / union_count if union_count else None,
        })
        rust_only_common = [
            key for key in common
            if rust[key]["label"] == "target" and rust[key]["q"] < threshold <= cpp[key]["q"]
        ]
        cpp_only_common = [
            key for key in common
            if cpp[key]["label"] == "target" and cpp[key]["q"] < threshold <= rust[key]["q"]
        ]
        exclusive_characterization.append({
            "q_threshold": threshold,
            "rust_only_unambiguous_matching_psms": len(rust_only_common),
            "rust_only_cpp_q": distribution(cpp[key]["q"] for key in rust_only_common),
            "rust_only_cpp_q_lt_2x_threshold": sum(cpp[key]["q"] < 2 * threshold for key in rust_only_common),
            "rust_only_cpp_q_ge_0.05": sum(cpp[key]["q"] >= 0.05 for key in rust_only_common),
            "rust_only_absolute_normalized_rank_difference": distribution(
                abs(rank_difference_by_key[key]) / len(common) for key in rust_only_common
            ),
            "cpp_only_unambiguous_matching_psms": len(cpp_only_common),
            "cpp_only_rust_q": distribution(rust[key]["q"] for key in cpp_only_common),
            "cpp_only_rust_q_lt_2x_threshold": sum(rust[key]["q"] < 2 * threshold for key in cpp_only_common),
            "cpp_only_rust_q_ge_0.05": sum(rust[key]["q"] >= 0.05 for key in cpp_only_common),
            "cpp_only_absolute_normalized_rank_difference": distribution(
                abs(rank_difference_by_key[key]) / len(common) for key in cpp_only_common
            ),
        })

    bands = []
    lower = 0.0
    for upper in THRESHOLDS:
        in_band = [key for key in common if lower <= min(rust[key]["q"], cpp[key]["q"]) < upper]
        bands.append({
            "lower_inclusive": lower,
            "upper_exclusive": upper,
            "matching_psms": len(in_band),
            "absolute_q_difference": distribution(abs(rust[key]["q"] - cpp[key]["q"]) for key in in_band),
            "absolute_rank_difference": distribution(
                abs(rank_difference_by_key[key]) for key in in_band
            ),
        })
        lower = upper

    cutoff = 0.01
    cross_cutoff = [
        key for key in common
        if (rust[key]["q"] < cutoff) != (cpp[key]["q"] < cutoff)
    ]
    cross_cutoff.sort(key=lambda key: abs(rust[key]["q"] - cpp[key]["q"]), reverse=True)
    largest_q = sorted(common, key=lambda key: abs(rust[key]["q"] - cpp[key]["q"]), reverse=True)

    return {
        "schema_version": 1,
        "analysis": "PSM-level percolator-rs versus C++ Percolator agreement",
        "evaluator": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "threshold_policy": "q < threshold",
        "inputs": {"rust": rust_files, "cpp": cpp_files},
        "ambiguous_psm_ids": {"rust": rust_ambiguity, "cpp": cpp_ambiguity},
        "row_counts": {
            "rust": sum(map(len, rust_groups.values())),
            "cpp": sum(map(len, cpp_groups.values())),
            "matching_unambiguous": len(common),
            "matching_multiset": sum(
                min(len(rust_groups[key]), len(cpp_groups[key]))
                for key in rust_groups.keys() & cpp_groups.keys()
            ),
            "rust_only": len(rust_only_rows), "cpp_only": len(cpp_only_rows),
        },
        "correlations_on_matching_psms": {
            "score_pearson": pearson(rs, cs),
            "score_spearman": spearman(rs, cs),
            "q_value_pearson": pearson(rq, cq),
            "q_value_spearman": spearman(rq, cq),
            "pep_pearson": pearson(rp, cp),
            "pep_spearman": spearman(rp, cp),
        },
        "differences_on_matching_psms": {
            "score_rust_minus_cpp": distribution(a - b for a, b in zip(rs, cs)),
            "q_rust_minus_cpp": distribution(a - b for a, b in zip(rq, cq)),
            "absolute_q": distribution(abs(a - b) for a, b in zip(rq, cq)),
            "pep_rust_minus_cpp": distribution(a - b for a, b in zip(rp, cp)),
            "absolute_pep": distribution(abs(a - b) for a, b in zip(rp, cp)),
            "rank_rust_minus_cpp": distribution(rank_differences),
            "absolute_rank": distribution(abs(value) for value in rank_differences),
            "absolute_normalized_rank": distribution(abs(value) / len(common) for value in rank_differences),
        },
        "thresholds": thresholds,
        "exclusive_discovery_characterization": exclusive_characterization,
        "q_bands_by_minimum_method_q": bands,
        "representative_disagreements": {
            "crossing_q_0.01": [representative(key, rust[key], cpp[key]) for key in cross_cutoff[:20]],
            "largest_absolute_q_difference": [representative(key, rust[key], cpp[key]) for key in largest_q[:20]],
        },
        "unmatched_examples": {
            "rust_only": [list(key) for key in rust_only_rows[:20]],
            "cpp_only": [list(key) for key in cpp_only_rows[:20]],
        },
    }


def write_threshold_table(path: Path, result: dict) -> None:
    with path.open("w", newline="") as handle:
        fields = ("q_threshold", "comparison", "rust", "cpp", "intersection", "rust_only",
                  "cpp_only", "union", "jaccard")
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        writer.writerows(result["thresholds"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rust", type=Path, required=True)
    parser.add_argument("--cpp", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--table", type=Path)
    args = parser.parse_args()
    result = compare(args.rust, args.cpp)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if args.table:
        args.table.parent.mkdir(parents=True, exist_ok=True)
        write_threshold_table(args.table, result)
    print(json.dumps({
        "matching_unambiguous_psms": result["row_counts"]["matching_unambiguous"],
        "matching_multiset_psms": result["row_counts"]["matching_multiset"],
        "score_spearman": result["correlations_on_matching_psms"]["score_spearman"],
        "q_0.01": next(row for row in result["thresholds"] if row["q_threshold"] == 0.01),
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
