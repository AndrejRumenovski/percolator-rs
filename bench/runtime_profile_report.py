#!/usr/bin/env python3
"""Aggregate feature-gated per-process profiles into JSON and Markdown."""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import re
import statistics
from pathlib import Path
from typing import Any


TAG_RE = re.compile(r"_(instrumented|cpu)_r\d+$")


def ns_ms(value: float) -> float:
    return value / 1_000_000.0


def percent(value: float, total: float) -> float:
    return 0.0 if total <= 0 else 100.0 * value / total


def add_counter(target: dict[str, dict[str, int]], name: str, event: dict[str, Any]) -> None:
    aggregate = target.setdefault(
        name, {"calls": 0, "total_ns": 0, "elements": 0, "bytes": 0}
    )
    aggregate["calls"] += 1
    aggregate["total_ns"] += event["duration_ns"]
    aggregate["elements"] += event.get("elements") or 0
    aggregate["bytes"] += event.get("bytes") or 0


def read_timings(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    rows: list[dict[str, Any]] = []
    for line in lines[1:]:
        if not line:
            continue
        row = dict(zip(header, line.split("\t")))
        for key in ("repetition", "wall_ns", "processes", "intra_file_threads", "cpu_sampling"):
            row[key] = int(row[key])
        rows.append(row)
    grouped: dict[tuple[str, str], list[int]] = collections.defaultdict(list)
    for row in rows:
        grouped[(row["configuration"], row["build"])].append(row["wall_ns"])
    summary: dict[str, Any] = {}
    for (configuration, build), values in sorted(grouped.items()):
        summary.setdefault(configuration, {})[build] = {
            "runs": len(values),
            "wall_ns": values,
            "median_wall_ns": int(statistics.median(values)),
        }
    for configuration, builds in summary.items():
        if "normal" in builds and "instrumented" in builds:
            normal_by_repetition = {
                row["repetition"]: row["wall_ns"]
                for row in rows
                if row["configuration"] == configuration and row["build"] == "normal"
            }
            instrumented_by_repetition = {
                row["repetition"]: row["wall_ns"]
                for row in rows
                if row["configuration"] == configuration
                and row["build"] == "instrumented"
            }
            paired_overheads = [
                percent(instrumented_by_repetition[repetition] - normal, normal)
                for repetition, normal in sorted(normal_by_repetition.items())
                if repetition in instrumented_by_repetition
            ]
            builds["paired_overhead_percent"] = paired_overheads
            builds["instrumentation_overhead_percent"] = statistics.median(
                paired_overheads
            )
            builds["overhead_method"] = "median of same-repetition paired deltas"
    return rows, summary


def profile_group(path: Path, artifacts: Path) -> str:
    relative = path.relative_to(artifacts / "profiles")
    return TAG_RE.sub("", relative.parts[0])


def aggregate_profiles(paths: list[Path]) -> dict[str, Any]:
    event_totals: dict[str, dict[str, int]] = {}
    category_name_totals: dict[str, dict[str, int]] = {}
    allocations: dict[str, dict[str, int]] = {}
    allocator = collections.Counter()
    iterations: dict[int, dict[str, dict[str, int]]] = collections.defaultdict(dict)
    folds: dict[int, dict[str, dict[str, int]]] = collections.defaultdict(dict)
    fold_heldout_ns: collections.Counter[int] = collections.Counter()
    phase_events: dict[str, dict[str, dict[str, int]]] = collections.defaultdict(dict)
    slowest_fold_counts: collections.Counter[int] = collections.Counter()
    fold_max_min_ratios: list[float] = []
    elapsed_ns = 0
    psms = 0
    input_bytes = 0
    cpu_collapsed: list[str] = []
    scheduling_ns = 0
    metadata_examples: list[dict[str, Any]] = []

    for path in paths:
        profile = json.loads(path.read_text())
        elapsed_ns += profile["elapsed_ns"]
        metadata = profile.get("metadata", {})
        psms += int(metadata.get("psms", 0))
        input_bytes += int(metadata.get("input_bytes", 0))
        if len(metadata_examples) < 3:
            metadata_examples.append(metadata)
        for key, value in profile.get("allocator", {}).items():
            allocator[key] += value
        for site in profile.get("allocation_sites", []):
            current = allocations.setdefault(site["name"], {"calls": 0, "bytes": 0})
            current["calls"] += site["calls"]
            current["bytes"] += site["bytes"]

        events = profile.get("events", [])
        fold_times: dict[int, int] = {}
        dispatch = 0
        for event in events:
            add_counter(event_totals, event["name"], event)
            add_counter(
                category_name_totals,
                f'{event["category"]}:{event["name"]}',
                event,
            )
            if event.get("iteration") is not None:
                add_counter(iterations[int(event["iteration"])], event["name"], event)
            add_counter(phase_events[event.get("phase") or "unscoped"], event["name"], event)
            if event.get("fold") is not None:
                fold_number = int(event["fold"])
                add_counter(folds[fold_number], event["name"], event)
                if (
                    event.get("phase") == "final_heldout_scoring"
                    and event["name"] == "model_score_rows"
                ):
                    fold_heldout_ns[fold_number] += event["duration_ns"]
            if event["name"] == "fold_total":
                fold_times[int(event["fold"])] = event["duration_ns"]
            elif event["name"] == "fold_dispatch_and_join":
                dispatch += event["duration_ns"]
        threads = int(metadata.get("num_threads", 1))
        if fold_times:
            durations = list(fold_times.values())
            useful = max(durations) if threads > 1 else sum(durations)
            scheduling_ns += max(0, dispatch - useful)
            slowest_fold_counts[max(fold_times, key=fold_times.get)] += 1
            fold_max_min_ratios.append(max(durations) / min(durations))
        cpu = profile.get("cpu")
        if cpu and cpu.get("collapsed_path"):
            cpu_collapsed.append(cpu["collapsed_path"])

    def named(name: str) -> dict[str, int]:
        return event_totals.get(name, {"calls": 0, "total_ns": 0, "elements": 0, "bytes": 0})

    top_stage_names = [
        "input_loading",
        "rescoring",
        "psm_level_processing",
        "peptide_level_processing",
        "result_output",
        "protein_inference_and_output",
    ]
    stages = {name: dict(named(name)) for name in top_stage_names if named(name)["calls"]}
    accounted = sum(stage["total_ns"] for stage in stages.values())
    stages["miscellaneous_unaccounted"] = {
        "calls": len(paths),
        "total_ns": max(0, elapsed_ns - accounted),
        "elements": 0,
        "bytes": 0,
    }
    for stage in stages.values():
        stage["percent_process_time"] = percent(stage["total_ns"], elapsed_ns)

    operation_names = [
        "pin_parse_total",
        "normalization_total",
        "initial_direction_selection",
        "fold_creation_and_setup",
        "model_score_rows",
        "qvalues_total",
        "confident_positive_selection",
        "svm_training_total",
        "pep_pava_total",
        "peptide_level_processing",
        "picked_protein_inference",
        "result_format_and_buffer",
        "result_file_write",
    ]
    operations = {name: dict(named(name)) for name in operation_names if named(name)["calls"]}
    for operation in operations.values():
        operation["percent_process_time"] = percent(operation["total_ns"], elapsed_ns)

    svm_names = [
        "allocation_and_buffer_initialization",
        "active_set_and_margin_scoring",
        "gradient_computation",
        "hessian_construction",
        "cholesky_factorization",
        "linear_solve",
        "convergence_logic",
        "solver_buffer_update",
        "line_search_weight_update",
        "line_search_total",
    ]
    svm_total = named("svm_training_total")["total_ns"]
    newton = named("newton_iteration_total")
    svm = {
        "total_ns": svm_total,
        "training_calls": named("svm_training_total")["calls"],
        "newton_iterations": newton["calls"],
        "mean_ns_per_newton_iteration": (
            newton["total_ns"] / newton["calls"] if newton["calls"] else 0.0
        ),
        "components": {},
    }
    for name in svm_names:
        component = dict(named(name))
        component["percent_svm_time"] = percent(component["total_ns"], svm_total)
        component["mean_ns_per_call"] = (
            component["total_ns"] / component["calls"] if component["calls"] else 0.0
        )
        svm["components"][name] = component

    iteration_output: list[dict[str, Any]] = []
    iteration_fields = {
        "training": "svm_training_total",
        "scoring": "model_score_rows",
        "qvalues": "qvalues_total",
        "positive_selection": "confident_positive_selection",
        "total": "iteration_total",
    }
    for iteration, totals in sorted(iterations.items()):
        row: dict[str, Any] = {"iteration": iteration}
        iteration_calls = totals.get("iteration_total", {}).get("calls", 0)
        row["fold_iteration_calls"] = iteration_calls
        for label, event_name in iteration_fields.items():
            value = totals.get(event_name, {}).get("total_ns", 0)
            row[f"{label}_total_ns"] = value
            row[f"{label}_mean_ns_per_fold"] = value / iteration_calls if iteration_calls else 0.0
        iteration_output.append(row)

    fold_output: list[dict[str, Any]] = []
    for fold, totals in sorted(folds.items()):
        process_calls = totals.get("fold_total", {}).get("calls", 0)
        row: dict[str, Any] = {"fold": fold, "process_calls": process_calls}
        for label, event_name in {
            "setup": "fold_setup",
            "training": "fold_training_total",
            "scoring": "model_score_rows",
            "qvalues": "qvalues_total",
            "heldout_scoring": "model_score_rows",
            "total": "fold_total",
        }.items():
            if label == "heldout_scoring":
                value = fold_heldout_ns[fold]
            elif label == "scoring":
                value = totals.get(event_name, {}).get("total_ns", 0)
                value -= fold_heldout_ns[fold]
            else:
                value = totals.get(event_name, {}).get("total_ns", 0)
            row[f"{label}_total_ns"] = value
            row[f"{label}_mean_ns_per_process"] = value / process_calls if process_calls else 0.0
        fold_output.append(row)

    parser_total = named("pin_parse_total")["total_ns"]
    parser = {
        "input_bytes": input_bytes,
        "psms": psms,
        "throughput_mib_per_second": (
            input_bytes / (1024 * 1024) / (parser_total / 1e9) if parser_total else 0.0
        ),
        "components": {
            name: dict(named(name))
            for name in [
                "mmap_setup",
                "header_and_feature_names",
                "row_loading_total",
                "field_splitting",
                "numeric_and_float_parsing",
                "string_allocation_and_copy",
                "pin_parse_total",
            ]
        },
    }

    sorts = []
    for key, value in category_name_totals.items():
        category, name = key.split(":", 1)
        if category == "sort":
            row = {"name": name, **value}
            row["percent_process_time"] = percent(row["total_ns"], elapsed_ns)
            sorts.append(row)
    sorts.sort(key=lambda item: (-item["total_ns"], item["name"]))

    allocation_sites = [
        {"name": name, **values}
        for name, values in sorted(
            allocations.items(), key=lambda item: (-item[1]["bytes"], item[0])
        )
    ]
    cpu = aggregate_cpu(cpu_collapsed)

    return {
        "profile_files": len(paths),
        "summed_process_elapsed_ns": elapsed_ns,
        "metadata_examples": metadata_examples,
        "top_level_stages": stages,
        "nested_operations": operations,
        "semi_supervised_iterations": iteration_output,
        "folds": fold_output,
        "fold_balance": {
            "slowest_fold_counts": dict(slowest_fold_counts),
            "mean_max_to_min_ratio": statistics.mean(fold_max_min_ratios)
            if fold_max_min_ratios
            else 0.0,
            "median_max_to_min_ratio": statistics.median(fold_max_min_ratios)
            if fold_max_min_ratios
            else 0.0,
        },
        "phase_events": phase_events,
        "parallel_scheduling_overhead_ns": scheduling_ns,
        "svm": svm,
        "parser": parser,
        "sorts": sorts,
        "allocator_totals": dict(allocator),
        "allocation_sites": allocation_sites,
        "cpu": cpu,
    }


def aggregate_cpu(paths: list[str]) -> dict[str, Any]:
    leaf: collections.Counter[str] = collections.Counter()
    inclusive: collections.Counter[str] = collections.Counter()
    total = 0
    valid_paths = []
    for raw_path in paths:
        path = Path(raw_path)
        if not path.exists():
            continue
        valid_paths.append(str(path))
        for line in path.read_text(errors="replace").splitlines():
            try:
                stack, raw_count = line.rsplit(" ", 1)
                count = int(raw_count)
            except ValueError:
                continue
            symbols = stack.split(";") if stack else []
            if not symbols:
                continue
            total += count
            leaf[symbols[-1]] += count
            for symbol in set(symbols):
                inclusive[symbol] += count

    def top(counter: collections.Counter[str]) -> list[dict[str, Any]]:
        return [
            {"symbol": symbol, "samples": count, "percent": percent(count, total)}
            for symbol, count in counter.most_common(30)
        ]

    return {
        "collapsed_profiles": valid_paths,
        "total_samples": total,
        "leaf_functions": top(leaf),
        "inclusive_functions": top(inclusive),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def output_equivalence(artifacts: Path) -> list[dict[str, Any]]:
    output = artifacts / "outputs"
    pairs = []
    for normal in sorted(output.glob("*_normal_r*")):
        instrumented = Path(str(normal).replace("_normal_r", "_instrumented_r"))
        if instrumented.exists():
            pairs.append((normal.name, normal, instrumented))
    special = [
        ("single_file_t1_cpu", output / "single_file_t1_instrumented_r1", output / "single_file_t1_cpu_cpu_r1"),
        ("full_n4_cpu", output / "full_n4_instrumented_r1", output / "full_n4_cpu_cpu_r1"),
        ("single_file_t1_allocations", output / "single_file_t1_instrumented_r1", output / "single_file_t1_allocations_instrumented_r1"),
    ]
    pairs.extend((name, left, right) for name, left, right in special if left.exists() and right.exists())
    results = []
    for name, left, right in pairs:
        left_files = {
            path.relative_to(left): path
            for path in left.rglob("*.tsv")
            if path.name not in {"summary.tsv"}
        }
        right_files = {
            path.relative_to(right): path
            for path in right.rglob("*.tsv")
            if path.name not in {"summary.tsv"}
        }
        differences = []
        for relative in sorted(set(left_files) | set(right_files)):
            if relative not in left_files or relative not in right_files:
                differences.append(str(relative))
            elif sha256(left_files[relative]) != sha256(right_files[relative]):
                differences.append(str(relative))
        results.append(
            {
                "comparison": name,
                "files_compared": len(set(left_files) | set(right_files)),
                "byte_identical": not differences,
                "differences": differences,
            }
        )
    return results


def workload_outcomes(artifacts: Path) -> dict[str, dict[str, int]]:
    psm_pattern = re.compile(r"target PSMs q<0\.01: (\d+)")
    peptide_pattern = re.compile(r"target peptides q<0\.01: (\d+)")
    outcomes: dict[str, dict[str, int]] = {}
    for directory in sorted((artifacts / "outputs").glob("full_*_r*")):
        psms = 0
        peptides = 0
        files = 0
        for log in directory.glob("*/stderr.log"):
            text = log.read_text(errors="replace")
            psm_match = psm_pattern.search(text)
            peptide_match = peptide_pattern.search(text)
            if psm_match and peptide_match:
                files += 1
                psms += int(psm_match.group(1))
                peptides += int(peptide_match.group(1))
        if files:
            outcomes[directory.name] = {
                "valid_files": files,
                "target_psms_q_lt_0_01": psms,
                "target_peptides_q_lt_0_01": peptides,
            }
    return outcomes


def environment(artifacts: Path) -> dict[str, str]:
    build = artifacts / "build"
    output = {}
    for name in ("git-head", "rustc", "cargo", "uname", "perf-probe-exit-code"):
        path = build / f"{name}.txt"
        if path.exists():
            output[name.replace("-", "_")] = path.read_text(errors="replace").strip()
    lscpu = build / "lscpu.txt"
    if lscpu.exists():
        for line in lscpu.read_text(errors="replace").splitlines():
            if line.startswith("Model name:"):
                output["cpu_model"] = line.split(":", 1)[1].strip()
            elif line.startswith("CPU(s):") and "cpu_count" not in output:
                output["cpu_count"] = line.split(":", 1)[1].strip()
    return output


def format_ms(value: float) -> str:
    if value >= 1e9:
        return f"{value / 1e9:.3f} s"
    return f"{ns_ms(value):.3f} ms"


def markdown_report(result: dict[str, Any]) -> str:
    lines = [
        "# percolator-rs runtime profile",
        "",
        "This report is generated from feature-gated wall-clock events and userspace sampled CPU profiles. Nested-operation percentages intentionally overlap; the top-level stage table does not.",
        "",
        "## Instrumentation overhead",
        "",
        "| Configuration | Normal median | Instrumented median | Overhead |",
        "|---|---:|---:|---:|",
    ]
    for configuration, builds in result["timing_summary"].items():
        if "normal" not in builds or "instrumented" not in builds:
            continue
        normal = builds["normal"]["median_wall_ns"]
        instrumented = builds["instrumented"]["median_wall_ns"]
        overhead = builds.get("instrumentation_overhead_percent", math.nan)
        lines.append(
            f"| {configuration} | {format_ms(normal)} | {format_ms(instrumented)} | {overhead:+.2f}% |"
        )

    preferred = "full_sequential"
    profile = result["configurations"].get(preferred)
    if profile:
        total = profile["summed_process_elapsed_ns"]
        lines += [
            "",
            f"## Runtime breakdown: {preferred}",
            "",
            "| Top-level stage | Summed process time | % process time |",
            "|---|---:|---:|",
        ]
        for name, stage in profile["top_level_stages"].items():
            lines.append(
                f'| {name} | {format_ms(stage["total_ns"])} | {stage["percent_process_time"]:.2f}% |'
            )
        lines += [
            "",
            "| Nested operation | Time | % process time |",
            "|---|---:|---:|",
        ]
        for name, operation in sorted(
            profile["nested_operations"].items(), key=lambda item: -item[1]["total_ns"]
        ):
            lines.append(
                f'| {name} | {format_ms(operation["total_ns"])} | {percent(operation["total_ns"], total):.2f}% |'
            )
        lines += [
            "",
            "## Semi-supervised iterations",
            "",
            "Values are mean wall time per fold invocation.",
            "",
            "| Iteration | Training | Scoring | q-values | Positive selection | Total |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
        for row in profile["semi_supervised_iterations"]:
            lines.append(
                "| {iteration} | {training} | {scoring} | {qvalues} | {positive} | {total} |".format(
                    iteration=row["iteration"],
                    training=format_ms(row["training_mean_ns_per_fold"]),
                    scoring=format_ms(row["scoring_mean_ns_per_fold"]),
                    qvalues=format_ms(row["qvalues_mean_ns_per_fold"]),
                    positive=format_ms(row["positive_selection_mean_ns_per_fold"]),
                    total=format_ms(row["total_mean_ns_per_fold"]),
                )
            )
        lines += [
            "",
            "## Per-fold breakdown",
            "",
            "Values are means per input process. Training scoring excludes final held-out scoring.",
            "",
            "| Fold | Setup | Training | Training scoring | q-values | Held-out scoring | Total |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for row in profile["folds"]:
            lines.append(
                "| {fold} | {setup} | {training} | {scoring} | {qvalues} | {heldout} | {total} |".format(
                    fold=row["fold"],
                    setup=format_ms(row["setup_mean_ns_per_process"]),
                    training=format_ms(row["training_mean_ns_per_process"]),
                    scoring=format_ms(row["scoring_mean_ns_per_process"]),
                    qvalues=format_ms(row["qvalues_mean_ns_per_process"]),
                    heldout=format_ms(row["heldout_scoring_mean_ns_per_process"]),
                    total=format_ms(row["total_mean_ns_per_process"]),
                )
            )
        lines += [
            "",
            "## SVM solver",
            "",
            f'Newton iterations: {profile["svm"]["newton_iterations"]}; mean {format_ms(profile["svm"]["mean_ns_per_newton_iteration"])} per Newton iteration.',
            "",
            "| Component | Time | Calls | % SVM time |",
            "|---|---:|---:|---:|",
        ]
        for name, component in sorted(
            profile["svm"]["components"].items(), key=lambda item: -item[1]["total_ns"]
        ):
            lines.append(
                f'| {name} | {format_ms(component["total_ns"])} | {component["calls"]} | {component["percent_svm_time"]:.2f}% |'
            )
        lines += [
            "",
            "## Sorts",
            "",
            "| Sort | Calls | Elements | Time | % process time |",
            "|---|---:|---:|---:|---:|",
        ]
        for sort in profile["sorts"]:
            lines.append(
                f'| {sort["name"]} | {sort["calls"]} | {sort["elements"]} | {format_ms(sort["total_ns"])} | {sort["percent_process_time"]:.2f}% |'
            )
        parser = profile["parser"]
        lines += [
            "",
            "## Parser and allocations",
            "",
            f'Parser throughput: {parser["throughput_mib_per_second"]:.1f} MiB/s over {parser["input_bytes"]} bytes and {parser["psms"]} PSMs.',
            "",
            "| Approximate allocation site | Calls | Bytes |",
            "|---|---:|---:|",
        ]
        for site in profile["allocation_sites"][:20]:
            lines.append(f'| {site["name"]} | {site["calls"]} | {site["bytes"]} |')

    lines += [
        "",
        "## CPU profiles",
        "",
    ]
    for configuration in ("single_file_t1_cpu", "full_n4_cpu"):
        cpu_profile = result["configurations"].get(configuration, {}).get("cpu", {})
        lines += [
            f"### {configuration}",
            "",
            "| Inclusive function | Samples | Approx. CPU |",
            "|---|---:|---:|",
        ]
        wrappers = (
            "std::rt::",
            "std::sys::backtrace::",
            "__libc_start",
        )
        useful_symbols = [
            symbol
            for symbol in cpu_profile.get("inclusive_functions", [])
            if symbol["symbol"] not in {"_start", "main", "percolator_rs::main"}
            and not symbol["symbol"].startswith(wrappers)
        ]
        for symbol in useful_symbols[:15]:
            lines.append(f'| `{symbol["symbol"]}` | {symbol["samples"]} | {symbol["percent"]:.2f}% |')
        lines.append("")

    lines += [
        "## Determinism",
        "",
        "| Comparison | Files | Byte-identical |",
        "|---|---:|---:|",
    ]
    for check in result["output_equivalence"]:
        lines.append(
            f'| {check["comparison"]} | {check["files_compared"]} | {str(check["byte_identical"]).lower()} |'
        )
    lines += [
        "",
        "Raw per-process JSON, protobuf CPU profiles, collapsed stacks, SVG flamegraphs, build metadata, and the perf capability probe are retained beside this report.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifacts", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--markdown", required=True, type=Path)
    args = parser.parse_args()

    timing_rows, timing_summary = read_timings(args.artifacts / "timings.tsv")
    grouped: dict[str, list[Path]] = collections.defaultdict(list)
    for path in sorted((args.artifacts / "profiles").rglob("*.json")):
        grouped[profile_group(path, args.artifacts)].append(path)
    configurations = {
        name: aggregate_profiles(paths) for name, paths in sorted(grouped.items())
    }
    result = {
        "schema_version": 1,
        "artifacts": str(args.artifacts),
        "environment": environment(args.artifacts),
        "timings": timing_rows,
        "timing_summary": timing_summary,
        "workload_outcomes": workload_outcomes(args.artifacts),
        "output_equivalence": output_equivalence(args.artifacts),
        "configurations": configurations,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(markdown_report(result))


if __name__ == "__main__":
    main()
