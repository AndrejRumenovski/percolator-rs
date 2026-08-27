#!/usr/bin/env python3
"""Frozen-method empirical reruns for the final scientific audit.

This driver keeps the predefined seeds, thresholds, entrapment adjustment,
compact datasets, and PrEST parameters unchanged.  It runs only the audited
Rust binary; no reference implementation is treated as an oracle.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path


THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.10)
PEP_EDGES = (0.0, 1e-12, 1e-6, 1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0000001)
SEEDS = (1, 2, 3, 4, 5)
ENTRAPMENT_FALLBACK = 0.50000027

ENTRAPMENT_INPUTS = (
    Path("/home/andrej-rumenovski/percolator_rs_out/entrapment/comet-09Dec2015-QEHF1-Anopheles-5-atrium-P-12hpm-3rd-01/comet.pin"),
    Path("/home/andrej-rumenovski/percolator_rs_out/entrapment/comet-22Oct2014-Anopheles-8-MAGs-S-01/comet.pin"),
    Path("/home/andrej-rumenovski/percolator_rs_out/entrapment/comet-28May2015-QE-HF-Anopheles-22-atrium-S-24H-2nd-01/comet.pin"),
    Path("/home/andrej-rumenovski/percolator_rs_out/entrapment/comet-28May2015-QE-HF-Anopheles-23-atrium-P-24H-2nd-01/comet.pin"),
    Path("/home/andrej-rumenovski/percolator_rs_out/entrapment/comet-28May2015-QE-HF-Anopheles-38-MAGs-P-3rd-02/comet.pin"),
    Path("/home/andrej-rumenovski/percolator_rs_out/entrapment/comet-9March2015-29-MAGs-pellet-2ndRep-14N-male-02/comet.pin"),
)

COMPACT_INPUTS = {
    "tide": Path("/home/andrej-rumenovski/percolator_rs_out/multidataset/inputs/hogrebe_tide.pin"),
    "msfragger": Path("/home/andrej-rumenovski/percolator_rs_out/multidataset/inputs/PXD020243_msfragger.pin"),
    "sage": Path("/home/andrej-rumenovski/percolator_rs_out/multidataset/inputs/PXD060954_sage.pin"),
    "yeast": Path("/home/andrej-rumenovski/percolator_rs_out/multidataset/inputs/percolator_yeast.pin"),
}

PROTEIN_ROOT = Path("/home/andrej-rumenovski/percolator_rs_out/protein-calibration")
PROTEIN_MANIFEST = PROTEIN_ROOT / "manifest.tsv"
PROTEIN_TRUTH = PROTEIN_ROOT / "input/ground-truth.tsv"
PROTEIN_PARAMS = {"alpha": 0.1, "beta": 0.0001, "gamma": 0.001, "max_iter": 1000}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def execute(command: list[str], directory: Path) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    started = time.time()
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    wall = time.time() - started
    (directory / "stdout.log").write_text(result.stdout)
    (directory / "stderr.log").write_text(result.stderr)
    record = {
        "argv": command,
        "exit_code": result.returncode,
        "wall_seconds": wall,
        "stdout": str(directory / "stdout.log"),
        "stderr": str(directory / "stderr.log"),
    }
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stderr}")
    return record


def proteins(row: dict[str | None, str | list[str]], decoy: bool) -> list[str]:
    values: list[str] = []
    for key, value in row.items():
        if key == "proteinIds" or key is None:
            if isinstance(value, list):
                values.extend(value)
            elif value:
                values.append(value)
    output = [protein for value in values for protein in value.replace(";", "\t").split("\t") if protein]
    if decoy:
        output = [protein.removeprefix("DECOY_").removeprefix("decoy_") for protein in output]
    return output


def load_psms(path: Path, decoy: bool) -> list[dict]:
    output = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            qvalue = float(row["q-value"])
            pep = float(row["posterior_error_prob"])
            members = proteins(row, decoy)
            pure = bool(members) and all(member.startswith("ENT_") for member in members)
            mixed = any(member.startswith("ENT_") for member in members) and not pure
            output.append({"q": qvalue, "pep": pep, "pure": pure, "mixed": mixed})
    return output


def entrapment_curve(targets: list[dict], decoys: list[dict]) -> list[dict]:
    output = []
    for threshold in THRESHOLDS:
        selected = [row for row in targets if row["q"] < threshold]
        usable = [row for row in decoys if row["q"] < threshold and not row["mixed"]]
        fraction = sum(row["pure"] for row in usable) / len(usable) if usable else ENTRAPMENT_FALLBACK
        known_false = sum(row["pure"] for row in selected)
        estimated_false = known_false / fraction
        output.append({
            "q_threshold": threshold,
            "accepted_targets": len(selected),
            "pure_entrapment": known_false,
            "mixed_targets": sum(row["mixed"] for row in selected),
            "effective_entrapment_fraction": fraction,
            "estimated_false": estimated_false,
            "adjusted_fdp": estimated_false / len(selected) if selected else 0.0,
        })
    return output


def pep_calibration(targets: list[dict], decoys: list[dict]) -> dict:
    usable = [row for row in decoys if not row["mixed"]]
    global_fraction = sum(row["pure"] for row in usable) / len(usable) if usable else ENTRAPMENT_FALLBACK
    bins = []
    for lower, upper in zip(PEP_EDGES, PEP_EDGES[1:]):
        selected = [row for row in targets if lower <= row["pep"] < upper]
        selected_decoys = [row for row in decoys if lower <= row["pep"] < upper and not row["mixed"]]
        fraction = (
            sum(row["pure"] for row in selected_decoys) / len(selected_decoys)
            if selected_decoys else global_fraction
        )
        known_false = sum(row["pure"] for row in selected)
        observed = known_false / fraction / len(selected) if selected and fraction else None
        mean_pep = statistics.fmean(row["pep"] for row in selected) if selected else None
        bins.append({
            "lower_inclusive": lower,
            "upper_exclusive": upper,
            "targets": len(selected),
            "mean_reported_pep": mean_pep,
            "pure_entrapment_targets": known_false,
            "decoys_for_fraction": len(selected_decoys),
            "effective_entrapment_fraction": fraction,
            "adjusted_observed_error": observed,
            "observed_minus_reported": observed - mean_pep if observed is not None else None,
        })
    populated = [row for row in bins if row["targets"] and row["observed_minus_reported"] is not None]
    total = sum(row["targets"] for row in populated)
    return {
        "targets": len(targets),
        "global_effective_entrapment_fraction": global_fraction,
        "exact_zero_peps": sum(row["pep"] == 0.0 for row in targets),
        "known_false_pep_lt_0.001": sum(row["pure"] and row["pep"] < 0.001 for row in targets),
        "minimum_known_false_pep": min((row["pep"] for row in targets if row["pure"]), default=None),
        "weighted_absolute_calibration_error": sum(row["targets"] * abs(row["observed_minus_reported"]) for row in populated) / total,
        "weighted_signed_observed_minus_reported": sum(row["targets"] * row["observed_minus_reported"] for row in populated) / total,
        "bins": bins,
    }


def describe(values: list[float]) -> dict:
    return {
        "n": len(values),
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "sd": statistics.stdev(values) if len(values) > 1 else 0.0,
        "minimum": min(values),
        "maximum": max(values),
    }


def run_entrapment(binary: Path, root: Path) -> dict:
    all_targets: list[dict] = []
    all_decoys: list[dict] = []
    runs = []
    for seed in SEEDS:
        seed_targets: list[dict] = []
        seed_decoys: list[dict] = []
        executions = []
        for pin in ENTRAPMENT_INPUTS:
            destination = root / "entrapment" / f"seed-{seed}" / pin.parent.name
            target = destination / "target.tsv"
            decoy = destination / "decoy.tsv"
            command = [
                str(binary), "--canonical", "--no-select-c", "--seed", str(seed),
                "--num-threads", "1", "--results-psms", str(target),
                "--decoy-results-psms", str(decoy), str(pin),
            ]
            executions.append(execute(command, destination))
            seed_targets.extend(load_psms(target, False))
            seed_decoys.extend(load_psms(decoy, True))
        curve = entrapment_curve(seed_targets, seed_decoys)
        runs.append({"seed": seed, "curve": curve, "executions": executions})
        all_targets.extend(seed_targets)
        all_decoys.extend(seed_decoys)
        q01 = next(row for row in curve if row["q_threshold"] == 0.01)
        print(f"entrapment seed={seed} accepted={q01['accepted_targets']} adjusted_fdp={q01['adjusted_fdp']:.6f}", flush=True)
    aggregate = []
    for threshold in THRESHOLDS:
        rows = [next(row for row in run["curve"] if row["q_threshold"] == threshold) for run in runs]
        aggregate.append({
            "q_threshold": threshold,
            "accepted_targets": describe([row["accepted_targets"] for row in rows]),
            "pure_entrapment": describe([row["pure_entrapment"] for row in rows]),
            "adjusted_fdp": describe([row["adjusted_fdp"] for row in rows]),
        })
    return {
        "inputs": [{"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for path in ENTRAPMENT_INPUTS],
        "seeds": list(SEEDS),
        "threshold_policy": "q < threshold",
        "fallback_entrapment_fraction": ENTRAPMENT_FALLBACK,
        "runs": runs,
        "aggregate": aggregate,
        "pep_calibration": pep_calibration(all_targets, all_decoys),
    }


def q_counts(path: Path) -> dict:
    with path.open(newline="") as handle:
        values = [float(row["q-value"]) for row in csv.DictReader(handle, delimiter="\t")]
    return {f"q_lt_{threshold:g}": sum(value < threshold for value in values) for threshold in THRESHOLDS}


def run_compact(binary: Path, root: Path) -> dict:
    runs = []
    for dataset, pin in COMPACT_INPUTS.items():
        for seed in SEEDS:
            repetitions = []
            for repeat in (1, 2):
                destination = root / "compact" / dataset / f"seed-{seed}" / f"repeat-{repeat}"
                target, decoy = destination / "target.tsv", destination / "decoy.tsv"
                command = [
                    str(binary), "--canonical", "--no-select-c", "--seed", str(seed),
                    "--num-threads", "1", "--results-psms", str(target),
                    "--decoy-results-psms", str(decoy), str(pin),
                ]
                execution = execute(command, destination)
                repetitions.append({
                    "repeat": repeat,
                    "execution": execution,
                    "target_sha256": sha256(target),
                    "decoy_sha256": sha256(decoy),
                    "counts": q_counts(target),
                })
            runs.append({
                "dataset": dataset,
                "seed": seed,
                "repetitions": repetitions,
                "byte_reproducible": all(
                    repetition["target_sha256"] == repetitions[0]["target_sha256"]
                    and repetition["decoy_sha256"] == repetitions[0]["decoy_sha256"]
                    for repetition in repetitions
                ),
            })
            print(f"compact {dataset} seed={seed} q01={repetitions[0]['counts']['q_lt_0.01']}", flush=True)
    aggregate = []
    for dataset in COMPACT_INPUTS:
        selected = [run for run in runs if run["dataset"] == dataset]
        for threshold in THRESHOLDS:
            key = f"q_lt_{threshold:g}"
            aggregate.append({
                "dataset": dataset,
                "q_threshold": threshold,
                "target_psms": describe([run["repetitions"][0]["counts"][key] for run in selected]),
            })
    return {
        "inputs": [{"id": name, "path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in COMPACT_INPUTS.items()],
        "seeds": list(SEEDS),
        "runs": runs,
        "aggregate": aggregate,
        "all_repeats_byte_identical": all(run["byte_reproducible"] for run in runs),
    }


def load_truth() -> dict[str, str]:
    with PROTEIN_TRUTH.open(newline="") as handle:
        return {row["protein_id"]: row["pool"] for row in csv.DictReader(handle, delimiter="\t")}


def present_pools(vial: str) -> set[str]:
    return {"A": {"A"}, "B": {"B"}, "AB": {"A", "B"}, "BLANK": set()}[vial]


def protein_groups(path: Path, vial: str, truth: dict[str, str]) -> list[dict]:
    output = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            members = [value for value in row["proteinIds"].split(",") if value]
            pools = [truth[member.removeprefix("DECOY_")] for member in members]
            error = not any(pool in present_pools(vial) for pool in pools)
            raw_pep = row["posterior_error_prob"]
            output.append({
                "q": float(row["q-value"]),
                "pep": None if raw_pep == "NA" else float(raw_pep),
                "error": error,
                "members": members,
            })
    return output


def probability_bins(groups: list[dict]) -> dict:
    usable = [group for group in groups if group["pep"] is not None]
    bins = []
    for index in range(10):
        lower, upper = index / 10.0, (index + 1) / 10.0
        selected = [group for group in usable if lower <= group["pep"] < upper or (index == 9 and group["pep"] == 1.0)]
        bins.append({
            "lower": lower,
            "upper": upper,
            "groups": len(selected),
            "mean_reported_pep": statistics.fmean(group["pep"] for group in selected) if selected else None,
            "observed_error": statistics.fmean(group["error"] for group in selected) if selected else None,
        })
    brier = statistics.fmean((float(group["error"]) - group["pep"]) ** 2 for group in usable) if usable else None
    ece = sum(
        row["groups"] * abs(row["observed_error"] - row["mean_reported_pep"])
        for row in bins if row["groups"]
    ) / len(usable) if usable else None
    return {"groups": len(usable), "brier": brier, "ece_10_bin": ece, "bins": bins}


def run_protein(binary: Path, root: Path) -> dict:
    truth = load_truth()
    with PROTEIN_MANIFEST.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    runs = []
    pooled: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for item in manifest:
        sample, vial, split, pin = item["sample"], item["vial"], item["split"], Path(item["pin"])
        for method in ("picked", "bayes-fixed", "bayes-selected"):
            destination = root / "protein" / "runs" / sample / method
            target, decoy = destination / "target.tsv", destination / "decoy.tsv"
            command = [
                str(binary), "--canonical", "--seed", "1", "--results-proteins", str(target),
                "--decoy-results-proteins", str(decoy),
            ]
            if method == "picked":
                command.extend(["--protein-inference", "picked"])
            else:
                command.extend(["--protein-inference", "bayesian"])
                if method == "bayes-selected":
                    command.extend([
                        "--protein-alpha", str(PROTEIN_PARAMS["alpha"]),
                        "--protein-beta", str(PROTEIN_PARAMS["beta"]),
                        "--protein-gamma", str(PROTEIN_PARAMS["gamma"]),
                        "--protein-max-iter", str(PROTEIN_PARAMS["max_iter"]),
                    ])
            command.append(str(pin))
            execution = execute(command, destination)
            (destination / "time.tsv").write_text(f"{execution['wall_seconds']}\t0\n")
            groups = protein_groups(target, vial, truth)
            accepted = [group for group in groups if group["q"] <= 0.01]
            record = {
                "sample": sample,
                "vial": vial,
                "split": split,
                "method": method,
                "input": str(pin),
                "input_sha256": sha256(pin),
                "execution": execution,
                "groups": len(groups),
                "accepted_q_le_0.01": len(accepted),
                "known_absent_q_le_0.01": sum(group["error"] for group in accepted),
                "raw_known_absent_fdp": statistics.fmean(group["error"] for group in accepted) if accepted else 0.0,
                "pep_all_na": all(group["pep"] is None for group in groups),
                "pep_all_numeric": all(group["pep"] is not None for group in groups),
            }
            runs.append(record)
            pooled[(split, method)].extend(groups)
            print(f"protein {sample} {method} accepted={len(accepted)} false={record['known_absent_q_le_0.01']}", flush=True)

    probability = {
        f"{split}/{method}": probability_bins(groups)
        for (split, method), groups in sorted(pooled.items())
        if method != "picked"
    }

    # Audit the committed evaluator against the current picked output.  It is
    # expected to reject the now-required `NA` field; that is preserved as a
    # validation-suite failure rather than patched during the audit.
    evaluator_output = root / "protein" / "committed-report"
    evaluator = subprocess.run([
        sys.executable, str(Path(__file__).resolve().parents[1] / "bench/protein_calibration/report.py"),
        "--truth", str(PROTEIN_TRUTH), "--manifest", str(PROTEIN_MANIFEST),
        "--runs", str(root / "protein" / "runs"), "--output-dir", str(evaluator_output),
    ], text=True, capture_output=True, check=False)
    (root / "protein" / "committed-report.stdout").write_text(evaluator.stdout)
    (root / "protein" / "committed-report.stderr").write_text(evaluator.stderr)
    return {
        "truth": {"path": str(PROTEIN_TRUTH), "sha256": sha256(PROTEIN_TRUTH)},
        "manifest": {"path": str(PROTEIN_MANIFEST), "sha256": sha256(PROTEIN_MANIFEST)},
        "selected_parameters_frozen_from_calibration_split": PROTEIN_PARAMS,
        "runs": runs,
        "bayesian_probability_calibration": probability,
        "picked_all_na": all(run["pep_all_na"] for run in runs if run["method"] == "picked"),
        "bayesian_all_numeric": all(run["pep_all_numeric"] for run in runs if run["method"] != "picked"),
        "committed_evaluator": {
            "exit_code": evaluator.returncode,
            "stdout": str(root / "protein" / "committed-report.stdout"),
            "stderr": str(root / "protein" / "committed-report.stderr"),
            "accepted_current_na_schema": evaluator.returncode == 0,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-entrapment", action="store_true")
    parser.add_argument("--skip-compact", action="store_true")
    parser.add_argument("--skip-protein", action="store_true")
    args = parser.parse_args()
    binary = args.binary.resolve()
    if args.output.exists():
        parser.error(f"output exists: {args.output}")
    required = [*ENTRAPMENT_INPUTS, *COMPACT_INPUTS.values(), PROTEIN_MANIFEST, PROTEIN_TRUTH]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        parser.error("missing predefined inputs:\n" + "\n".join(missing))
    args.output.mkdir(parents=True)
    result = {
        "schema_version": 1,
        "audit": "frozen-method final empirical rerun",
        "environment": {
            "platform": platform.platform(),
            "python": sys.version,
            "binary": str(binary),
            "binary_sha256": sha256(binary),
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256(Path(__file__).resolve()),
        },
    }
    if not args.skip_entrapment:
        result["entrapment"] = run_entrapment(binary, args.output)
    if not args.skip_compact:
        result["compact"] = run_compact(binary, args.output)
    if not args.skip_protein:
        result["protein"] = run_protein(binary, args.output)
    manifest = args.output / "manifest.json"
    manifest.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"manifest: {manifest}")


if __name__ == "__main__":
    main()
