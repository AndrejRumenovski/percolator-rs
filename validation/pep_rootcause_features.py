#!/usr/bin/env python3
"""Feature and local-density audit in the extreme entrapment score tail.

Joins reported PSMs back to their PIN rows by (dataset, SpecId, label).  The
analysis is descriptive: models are not refit and PEPs are not altered.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pep_rootcause_experiments import load_outputs


TAIL_BANDS = ((0.0, .001), (.001, .005), (.005, .01), (.01, .02), (.02, .05), (.05, .10))


def locate_pin(pin_root: Path, dataset: str) -> Path:
    candidates = [
        pin_root / dataset / "comet.pin",
        pin_root / f"comet-{dataset}" / "comet.pin",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"PIN for {dataset} under {pin_root}")


def read_pin(path: Path) -> tuple[list[str], dict[tuple[str, bool], np.ndarray]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        index = {name: i for i, name in enumerate(header)}
        excluded = {"SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "Peptide", "Proteins"}
        names = [name for name in header if name not in excluded]
        positions = [index[name] for name in names]
        values = {}
        for record in reader:
            key = (record[index["SpecId"]], int(record[index["Label"]]) < 0)
            values[key] = np.asarray([float(record[position]) for position in positions])
    return names, values


def describe(features: np.ndarray, names: list[str], target: np.ndarray, decoy: np.ndarray) -> dict:
    rows = []
    for column, name in enumerate(names):
        t = features[target, column]
        d = features[decoy, column]
        if not t.size or not d.size:
            rows.append({"feature": name, "n_target": int(t.size), "n_decoy": int(d.size)})
            continue
        variance = ((t.size - 1) * t.var(ddof=1) + (d.size - 1) * d.var(ddof=1)) / max(t.size + d.size - 2, 1)
        sd = math.sqrt(max(variance, 0.0))
        rows.append({
            "feature": name, "n_target": int(t.size), "n_decoy": int(d.size),
            "target_mean": float(t.mean()), "decoy_mean": float(d.mean()),
            "standardized_difference_target_minus_decoy": float((t.mean() - d.mean()) / sd) if sd else 0.0,
            "target_q10_q50_q90": [float(x) for x in np.quantile(t, [.1, .5, .9])],
            "decoy_q10_q50_q90": [float(x) for x in np.quantile(d, [.1, .5, .9])],
        })
    rows.sort(key=lambda row: abs(row.get("standardized_difference_target_minus_decoy", 0.0)), reverse=True)
    return {"all_features_by_absolute_smd": rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--pin-root", type=Path, required=True)
    parser.add_argument("--target-name", default="target.tsv")
    parser.add_argument("--method-subdir")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = load_outputs(args.results_root, args.target_name, args.method_subdir)
    matrices = []
    feature_names = None
    matched = np.zeros(data["pep"].size, dtype=bool)
    for dataset in sorted(set(str(x) for x in data["dataset_name"])):
        selected = data["dataset_name"] == dataset
        names, pin = read_pin(locate_pin(args.pin_root, dataset))
        if feature_names is None:
            feature_names = names
        elif names != feature_names:
            raise ValueError(f"feature layout differs for {dataset}")
        matrix = np.full((int(selected.sum()), len(names)), np.nan)
        for row, global_index in enumerate(np.flatnonzero(selected)):
            key = (str(data["psmid"][global_index]), bool(data["decoy"][global_index]))
            value = pin.get(key)
            if value is not None:
                matrix[row] = value
                matched[global_index] = True
        matrices.append((np.flatnonzero(selected), matrix))
    features = np.full((data["pep"].size, len(feature_names)), np.nan)
    for indices, matrix in matrices:
        features[indices] = matrix
    if not matched.all():
        raise ValueError(f"only {matched.sum()} / {matched.size} output rows matched PIN rows")

    pure_target = data["pure"] & ~data["decoy"]
    pure_decoy = data["pure"] & data["decoy"]
    regions = []
    masks = {
        "PEP < 0.001": data["pep"] < .001,
        "PEP < 0.01": data["pep"] < .01,
        "PEP < 0.05": data["pep"] < .05,
        "q < 0.01": data["q"] < .01,
    }
    for name, region in masks.items():
        target = region & pure_target
        decoy = region & pure_decoy
        entry = {
            "region": name, "entrapment_targets": int(target.sum()),
            "entrapment_decoys": int(decoy.sum()),
            "target_over_decoy": float(target.sum() / decoy.sum()) if decoy.any() else None,
        }
        entry.update(describe(features, feature_names, target, decoy))
        regions.append(entry)

    tail_bands = []
    percentile = np.full(data["score"].size, np.nan)
    # Best score has percentile zero.  Compute within each seed x dataset.
    for seed in np.unique(data["seed"]):
        for dataset in np.unique(data["dataset"]):
            group = (data["seed"] == seed) & (data["dataset"] == dataset)
            indices = np.flatnonzero(group)
            if not indices.size:
                continue
            order = indices[np.argsort(-data["score"][indices], kind="stable")]
            percentile[order] = np.arange(order.size) / order.size
    for lo, hi in TAIL_BANDS:
        region = (percentile >= lo) & (percentile < hi)
        t = int((region & pure_target).sum())
        d = int((region & pure_decoy).sum())
        tail_bands.append({
            "percentile_band": [lo, hi], "rows": int(region.sum()),
            "entrapment_targets": t, "entrapment_decoys": d,
            "local_target_over_decoy": t / d if d else None,
        })

    output = {
        "matched_rows": int(matched.sum()), "feature_names": feature_names,
        "regions": regions, "score_percentile_bands": tail_bands,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
