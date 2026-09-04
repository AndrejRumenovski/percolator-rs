#!/usr/bin/env python3
"""Reproducible real-data controls for the PEP root-cause investigation.

This file is deliberately outside the production crate.  It can:

* copy PIN files while removing named feature columns (the enzN/enzC ablation);
* summarize a set of percolator-rs PSM outputs without changing any PEP; and
* measure the internal entrapment-target/entrapment-decoy null in score tails.

Entrapment-adjusted error divides the pure-entrapment target fraction by the
global pure-entrapment fraction among non-mixed decoys.  The script always also
reports the adjustment-free ``f=1`` lower bound and the internal null counts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


PEP_EDGES = np.array(
    [0.0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0 + 1e-12]
)
PEP_NAMES = [
    "[0,1e-4)", "[1e-4,1e-3)", "[1e-3,5e-3)", "[5e-3,.01)",
    "[.01,.02)", "[.02,.05)", "[.05,.10)", "[.10,.20)",
    "[.20,.50)", "[.50,1]",
]
PEP_REGIONS = (1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1)
Q_REGIONS = (1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1)
TAIL_FRACTIONS = (0.001, 0.005, 0.01, 0.02, 0.05)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def classify(protein_fields: list[str], decoy: bool) -> tuple[bool, bool]:
    members = [
        member
        for field in protein_fields
        for member in field.replace(";", "\t").split("\t")
        if member
    ]
    if decoy:
        members = [
            member.removeprefix("DECOY_").removeprefix("decoy_")
            for member in members
        ]
    pure = bool(members) and all(member.startswith("ENT_") for member in members)
    mixed = any(member.startswith("ENT_") for member in members) and not pure
    return pure, mixed


def prepare_ablation(input_root: Path, output_root: Path, drop: list[str]) -> dict:
    pins = sorted(
        path for path in input_root.glob("comet-*/comet.pin")
        if path.parent.name != "comet-out"
    )
    if len(pins) != 6:
        raise ValueError(f"expected six canonical PINs, found {len(pins)}")
    manifest = {"operation": "drop PIN feature columns", "drop": drop, "inputs": []}
    for source in pins:
        destination = output_root / source.parent.name / "comet.pin"
        destination.parent.mkdir(parents=True, exist_ok=True)
        rows = 0
        with source.open(newline="") as inp, destination.open("w", newline="") as out:
            reader = csv.reader(inp, delimiter="\t")
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            header = next(reader)
            missing = sorted(set(drop) - set(header))
            if missing:
                raise ValueError(f"{source}: missing columns {missing}")
            keep = [index for index, name in enumerate(header) if name not in drop]
            writer.writerow([header[index] for index in keep])
            for record in reader:
                if len(record) < len(header):
                    raise ValueError(f"{source}: short row {rows + 2}")
                # Proteins may occupy several trailing tab fields.  The named
                # feature columns all precede Proteins, so preserve that tail.
                fixed = [record[index] for index in keep]
                if len(record) > len(header):
                    fixed.extend(record[len(header):])
                writer.writerow(fixed)
                rows += 1
        manifest["inputs"].append({
            "source": str(source), "source_sha256": sha256(source),
            "destination": str(destination), "destination_sha256": sha256(destination),
            "rows": rows,
        })
    return manifest


def protein_fields(row: dict) -> list[str]:
    fields: list[str] = []
    for key, value in row.items():
        if key == "proteinIds" or key is None:
            if isinstance(value, list):
                fields.extend(value)
            elif value:
                fields.append(value)
    return fields


def core_peptide(peptide: str) -> str:
    if len(peptide) > 4 and peptide[1] == "." and peptide[-2] == ".":
        peptide = peptide[2:-2]
    result = []
    bracket = 0
    for character in peptide:
        if character == "[":
            bracket += 1
        elif character == "]":
            bracket = max(0, bracket - 1)
        elif bracket == 0 and character.isalpha() and character.isupper():
            result.append(character)
    return "".join(result)


def find_result_files(root: Path, target_name: str, method_subdir: str | None) -> list[tuple[int, str, Path, Path]]:
    found = []
    for seed_dir in sorted(root.glob("seed-*")):
        try:
            seed = int(seed_dir.name.split("-", 1)[1])
        except (IndexError, ValueError):
            continue
        base = seed_dir / method_subdir if method_subdir else seed_dir
        for target in sorted(base.glob(f"*/{target_name}")):
            decoy_name = target_name.replace("target", "decoy", 1)
            decoy = target.with_name(decoy_name)
            if not decoy.exists():
                raise FileNotFoundError(decoy)
            found.append((seed, target.parent.name, target, decoy))
    if not found:
        raise ValueError(f"no results under {root}")
    return found


def load_outputs(root: Path, target_name: str, method_subdir: str | None) -> dict[str, np.ndarray]:
    columns: dict[str, list] = defaultdict(list)
    files = []
    dataset_code: dict[str, int] = {}
    for seed, dataset, target, decoy in find_result_files(root, target_name, method_subdir):
        dataset_code.setdefault(dataset, len(dataset_code))
        for is_decoy, path in ((False, target), (True, decoy)):
            count = 0
            with path.open(newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    pure, mixed = classify(protein_fields(row), is_decoy)
                    pep_text = row["posterior_error_prob"]
                    pep = float("nan") if pep_text == "NA" else float(pep_text)
                    columns["seed"].append(seed)
                    columns["dataset"].append(dataset_code[dataset])
                    columns["dataset_name"].append(dataset)
                    columns["decoy"].append(is_decoy)
                    columns["pure"].append(pure)
                    columns["mixed"].append(mixed)
                    columns["score"].append(float(row["score"]))
                    columns["q"].append(float(row["q-value"]))
                    columns["pep"].append(pep)
                    columns["psmid"].append(row["PSMId"])
                    columns["peptide"].append(core_peptide(row["peptide"]))
                    count += 1
            files.append({"seed": seed, "dataset": dataset, "decoy": is_decoy,
                          "path": str(path), "sha256": sha256(path), "rows": count})
    numeric = {
        "seed": np.int16, "dataset": np.int16, "decoy": bool, "pure": bool,
        "mixed": bool, "score": float, "q": float, "pep": float,
    }
    result = {
        key: np.asarray(values, dtype=numeric.get(key, object))
        for key, values in columns.items()
    }
    result["files"] = np.asarray(files, dtype=object)
    return result


def f_global(data: dict[str, np.ndarray], mask: np.ndarray) -> float:
    usable = mask & data["decoy"] & ~data["mixed"]
    return float(data["pure"][usable].sum()) / max(int(usable.sum()), 1)


def counts(data: dict[str, np.ndarray], target_mask: np.ndarray, decoy_mask: np.ndarray) -> dict:
    t = target_mask & ~data["decoy"]
    d = decoy_mask & data["decoy"]
    ent_t = int((t & data["pure"]).sum())
    ent_d = int((d & data["pure"]).sum())
    native_t = int((t & ~data["pure"] & ~data["mixed"]).sum())
    return {
        "targets": int(t.sum()), "normal_targets": native_t,
        "entrapment_targets": ent_t, "entrapment_decoys": ent_d,
        "mixed_targets": int((t & data["mixed"]).sum()),
        "ent_target_over_ent_decoy": ent_t / ent_d if ent_d else None,
    }


def calibration(data: dict[str, np.ndarray], mask: np.ndarray) -> dict:
    f = f_global(data, mask)
    target = mask & ~data["decoy"]
    decoy = mask & data["decoy"]
    index = np.digitize(data["pep"], PEP_EDGES) - 1
    rows = []
    for bin_index, name in enumerate(PEP_NAMES):
        mt = target & (index == bin_index)
        if not mt.any():
            continue
        md = decoy & (index == bin_index)
        n = int(mt.sum())
        predicted = float(data["pep"][mt].mean())
        ent_t = int((mt & data["pure"]).sum())
        ent_d = int((md & data["pure"]).sum())
        native_d = int((md & ~data["pure"] & ~data["mixed"]).sum())
        observed_f1 = ent_t / n
        observed_adjusted = observed_f1 / f if f else math.nan
        # Wilson interval for the directly observed entrapment fraction.  It is
        # not an interval for the extrapolated total error probability.
        z = 1.959963984540054
        center = (observed_f1 + z * z / (2 * n)) / (1 + z * z / n)
        half = z * math.sqrt(observed_f1 * (1 - observed_f1) / n + z * z / (4 * n * n)) / (1 + z * z / n)
        rows.append({
            "bin": name, "n": n, "mean_predicted_pep": predicted,
            "entrapment_targets": ent_t, "entrapment_decoys": ent_d,
            "native_decoys": native_d, "observed_f1": observed_f1,
            "observed_f1_wilson95": [max(0.0, center - half), min(1.0, center + half)],
            "observed_adjusted": observed_adjusted,
            "predicted_minus_observed_f1": predicted - observed_f1,
            "predicted_minus_observed_adjusted": predicted - observed_adjusted,
            "observed_over_predicted_adjusted": observed_adjusted / predicted if predicted else None,
            "ent_target_over_ent_decoy": ent_t / ent_d if ent_d else None,
        })
    n = sum(row["n"] for row in rows)
    return {
        "f_global": f,
        "bins": rows,
        "summary": {
            "n_targets": n,
            "target_pep_sum": float(data["pep"][target].sum()),
            "entrapment_targets": int((target & data["pure"]).sum()),
            "pure_entrapment_decoys": int((decoy & data["pure"] & ~data["mixed"]).sum()),
            "usable_decoys": int((decoy & ~data["mixed"]).sum()),
            "observed_minus_predicted_adjusted": sum(
                row["n"] * -row["predicted_minus_observed_adjusted"] for row in rows
            ) / n,
            "observed_minus_predicted_f1": sum(
                row["n"] * -row["predicted_minus_observed_f1"] for row in rows
            ) / n,
        },
    }


def cluster_bootstrap(data: dict[str, np.ndarray], unit_field: str, replicates: int = 4000) -> dict:
    """Bootstrap pooled calibration with multiplicity preserved.

    The earlier audit script formed a union mask after resampling units, which
    accidentally discarded repeated draws.  Here every unit is reduced to
    sufficient counts/sums and a multinomial draw weights those quantities.
    """
    units = np.unique(data[unit_field])
    bins = np.digitize(data["pep"], PEP_EDGES) - 1
    # Columns per unit: target n, target PEP sum, entrapment targets,
    # pure-entrapment decoys, all non-mixed decoys.  Row 0 is the pooled total;
    # following rows are PEP bins (the two decoy columns stay global there).
    sufficient = np.zeros((units.size, len(PEP_NAMES) + 1, 5), dtype=float)
    for ui, unit in enumerate(units):
        group = data[unit_field] == unit
        target = group & ~data["decoy"]
        decoy = group & data["decoy"]
        sufficient[ui, 0] = [
            target.sum(), data["pep"][target].sum(), (target & data["pure"]).sum(),
            (decoy & data["pure"] & ~data["mixed"]).sum(),
            (decoy & ~data["mixed"]).sum(),
        ]
        for bi in range(len(PEP_NAMES)):
            mt = target & (bins == bi)
            sufficient[ui, bi + 1, :3] = [
                mt.sum(), data["pep"][mt].sum(), (mt & data["pure"]).sum(),
            ]
            sufficient[ui, bi + 1, 3:] = sufficient[ui, 0, 3:]

    def statistics(weight: np.ndarray) -> tuple[float, np.ndarray]:
        total = np.tensordot(weight, sufficient, axes=(0, 0))
        global_row = total[0]
        f = global_row[3] / global_row[4] if global_row[4] else math.nan
        signed = (global_row[2] / f - global_row[1]) / global_row[0]
        gaps = np.full(len(PEP_NAMES), np.nan)
        for bi, row in enumerate(total[1:]):
            if row[0] and f:
                gaps[bi] = row[2] / f / row[0] - row[1] / row[0]
        return signed, gaps

    point, point_bins = statistics(np.ones(units.size))
    rng = np.random.default_rng(20260830)
    pooled = np.empty(replicates)
    by_bin = np.empty((replicates, len(PEP_NAMES)))
    for replicate in range(replicates):
        draw = rng.integers(0, units.size, size=units.size)
        weight = np.bincount(draw, minlength=units.size)
        pooled[replicate], by_bin[replicate] = statistics(weight)
    bin_rows = []
    for bi, name in enumerate(PEP_NAMES):
        values = by_bin[:, bi]
        values = values[np.isfinite(values)]
        bin_rows.append({
            "bin": name, "point_observed_minus_predicted": float(point_bins[bi]),
            "replicates_populated": int(values.size),
            "percentile95": [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]
            if values.size else None,
        })
    return {
        "unit": unit_field, "n_units": int(units.size), "replicates": replicates,
        "point_observed_minus_predicted": float(point),
        "percentile95": [float(np.percentile(pooled, 2.5)), float(np.percentile(pooled, 97.5))],
        "bins": bin_rows,
    }


def cumulative_regions(data: dict[str, np.ndarray], mask: np.ndarray, field: str, thresholds: tuple[float, ...]) -> list[dict]:
    f = f_global(data, mask)
    out = []
    for threshold in thresholds:
        selected = mask & (data[field] < threshold)
        row = counts(data, selected, selected)
        t = row["targets"]
        ent = row["entrapment_targets"]
        predicted = float(data["pep"][selected & ~data["decoy"]].mean()) if t else None
        row.update({
            "threshold": threshold,
            "mean_predicted_pep": predicted,
            "adjusted_fdp": ent / f / t if f and t else None,
            "f1_fdp": ent / t if t else None,
        })
        out.append(row)
    return out


def tail_analysis(data: dict[str, np.ndarray], mask: np.ndarray) -> list[dict]:
    out = []
    for fraction in TAIL_FRACTIONS:
        selected = np.zeros(mask.size, dtype=bool)
        # Scores are comparable within one trained file, not necessarily across
        # files, so define every percentile within each seed x dataset list.
        for seed in np.unique(data["seed"][mask]):
            for dataset in np.unique(data["dataset"][mask & (data["seed"] == seed)]):
                group = mask & (data["seed"] == seed) & (data["dataset"] == dataset)
                indices = np.flatnonzero(group)
                take = max(1, int(math.ceil(fraction * indices.size)))
                order = indices[np.argsort(-data["score"][indices], kind="stable")]
                selected[order[:take]] = True
        row = counts(data, selected, selected)
        row.update({"top_fraction": fraction, "rows": int(selected.sum())})
        out.append(row)
    return out


def score_distributions(data: dict[str, np.ndarray], mask: np.ndarray) -> dict:
    classes = {
        "normal_targets": mask & ~data["decoy"] & ~data["pure"] & ~data["mixed"],
        "entrapment_targets": mask & ~data["decoy"] & data["pure"],
        "entrapment_decoys": mask & data["decoy"] & data["pure"],
        "native_decoys": mask & data["decoy"] & ~data["pure"] & ~data["mixed"],
    }
    result = {}
    for name, selected in classes.items():
        values = data["score"][selected]
        result[name] = {
            "n": int(values.size),
            "mean": float(values.mean()) if values.size else None,
            "q01_q10_q50_q90_q99": [float(x) for x in np.quantile(values, [.01, .1, .5, .9, .99])] if values.size else None,
        }
    return result


def summarize(data: dict[str, np.ndarray]) -> dict:
    all_rows = np.ones(data["pep"].size, dtype=bool)
    result = {
        "rows": int(all_rows.size),
        "seeds": [int(x) for x in np.unique(data["seed"])],
        "datasets": sorted(set(str(x) for x in data["dataset_name"])),
        "calibration": calibration(data, all_rows),
        "pep_regions": cumulative_regions(data, all_rows, "pep", PEP_REGIONS),
        "q_regions": cumulative_regions(data, all_rows, "q", Q_REGIONS),
        "score_tails": tail_analysis(data, all_rows),
        "score_distributions": score_distributions(data, all_rows),
        "uncertainty": {
            "dataset_cluster_bootstrap": cluster_bootstrap(data, "dataset"),
            "seed_reproducibility_bootstrap": cluster_bootstrap(data, "seed"),
        },
        "by_seed": [],
        "by_dataset": [],
    }
    for seed in np.unique(data["seed"]):
        mask = data["seed"] == seed
        result["by_seed"].append({
            "seed": int(seed), "calibration": calibration(data, mask),
            "pep_regions": cumulative_regions(data, mask, "pep", PEP_REGIONS),
            "q_regions": cumulative_regions(data, mask, "q", Q_REGIONS),
        })
    for code in np.unique(data["dataset"]):
        mask = data["dataset"] == code
        name = str(data["dataset_name"][np.flatnonzero(mask)[0]])
        result["by_dataset"].append({
            "dataset": name, "calibration": calibration(data, mask),
            "pep_regions": cumulative_regions(data, mask, "pep", PEP_REGIONS),
            "q_regions": cumulative_regions(data, mask, "q", Q_REGIONS),
        })
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    prep = sub.add_parser("prepare-ablation")
    prep.add_argument("--input-root", type=Path, required=True)
    prep.add_argument("--output-root", type=Path, required=True)
    prep.add_argument("--drop", nargs="+", default=["enzN", "enzC"])
    prep.add_argument("--manifest", type=Path, required=True)

    report = sub.add_parser("summarize")
    report.add_argument("--root", type=Path, required=True)
    report.add_argument("--target-name", default="target.tsv")
    report.add_argument("--method-subdir")
    report.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()
    if args.command == "prepare-ablation":
        result = prepare_ablation(args.input_root, args.output_root, args.drop)
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    else:
        data = load_outputs(args.root, args.target_name, args.method_subdir)
        result = summarize(data)
        result["inputs"] = list(data["files"])
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
