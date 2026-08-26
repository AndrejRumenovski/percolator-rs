#!/usr/bin/env python3
"""Root-cause probe for the anti-conservative signal-present entrapment curve.

The five-seed entrapment study reports a mean adjusted FDP well above nominal at
every threshold.  That is a validation failure, but it does not by itself say
*what* is anti-conservative.  This probe separates the candidates that can be
separated with the data at hand, using an estimator written here rather than
called out of the production crate.

Arms
----

``raw``
    Target-decoy competition and q-values computed directly on the best single
    PIN feature, with no semi-supervised rescoring at all.  If the raw search
    score is already anti-conservative by the same amount, semi-supervised
    training is not the cause.  If it is calibrated and the rescored arm is not,
    training is implicated.

``rescored``
    The same entrapment accounting applied to a percolator-rs result directory,
    for a like-for-like comparison against ``raw``.

Sensitivity
-----------

The adjusted FDP divides pure-entrapment target counts by an estimate of the
probability that an incorrect target lands in the foreign database.  That
probability is estimated from accepted non-mixed decoys, which assumes incorrect
targets and decoys place into the foreign database at the same rate.  The probe
reports the adjusted FDP under that estimate, under the declared database share,
and under the whole-file decoy fraction, so the size of the assumption is visible
rather than buried.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.1)


def mix64(z: int) -> int:
    mask = (1 << 64) - 1
    z &= mask
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & mask
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & mask
    return (z ^ (z >> 31)) & mask


def coin(seed: int, *values: int) -> int:
    state = mix64(seed ^ 0x243F6A8885A308D3)
    for value in values:
        state = mix64(state ^ mix64(value & ((1 << 64) - 1)))
    return state


def classify(proteins: str) -> tuple[bool, bool, bool]:
    """(is_decoy, pure_entrapment, mixed) for a raw Proteins field."""
    members = [p for p in proteins.replace(";", "\t").replace(" ", "\t").split("\t") if p]
    decoy = any(p.upper().startswith(("DECOY_", "REV_", "RANDOM_", "RANDOM-")) for p in members)
    stripped = [
        p[6:] if p.upper().startswith("DECOY_") else p for p in members
    ]
    foreign = [p.startswith("ENT_") for p in stripped]
    pure = bool(foreign) and all(foreign)
    mixed = any(foreign) and not pure
    return decoy, pure, mixed


def qvalues(scores, labels, null_target_win_prob=0.5):
    """Same estimator as the production scan, written out here independently."""
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    factor = null_target_win_prob / (1.0 - null_target_win_prob)
    fdp = [1.0] * len(order)
    targets = 0.0
    decoys = 1.0
    start = 0
    for rank, index in enumerate(order):
        if labels[index] > 0:
            targets += 1.0
        else:
            decoys += 1.0
        last = rank + 1 == len(order) or scores[order[rank + 1]] != scores[index]
        if last:
            value = min(1.0, decoys * factor / max(1.0, targets))
            for position in range(start, rank + 1):
                fdp[position] = value
            start = rank + 1
    out = [1.0] * len(order)
    running = 1.0
    for rank in range(len(order) - 1, -1, -1):
        running = min(running, fdp[rank])
        out[order[rank]] = running
    return out


def compete(rows, scores, seed):
    """One winner per (ScanNr, ExpMass), exact ties drawn by a fair coin."""
    best: dict[tuple[int, str], list[int]] = {}
    for index, row in enumerate(rows):
        key = (row["scan"], row["mass"])
        current = best.get(key)
        if current is None or scores[index] > scores[current[0]]:
            best[key] = [index]
        elif scores[index] == scores[current[0]]:
            current.append(index)
    winners = []
    for key, group in best.items():
        if len(group) == 1:
            winners.append(group[0])
            continue
        canonical = sorted(group, key=lambda i: (rows[i]["label"], rows[i]["peptide"], rows[i]["spec"]))
        draw = (coin(seed, key[0], hash(key[1]) & ((1 << 64) - 1)) * len(canonical)) >> 64
        winners.append(canonical[draw])
    return winners


def read_pin(path: Path):
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        # Proteins is the trailing free field; csv gives it as extra columns.
        index = {name: position for position, name in enumerate(header)}
        feature_names = [
            name
            for name in header
            if name
            not in {"SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "Peptide", "Proteins"}
        ]
        feature_positions = [index[name] for name in feature_names]
        protein_start = index["Proteins"]
        rows = []
        features = []
        for record in reader:
            if len(record) <= protein_start:
                continue
            decoy, pure, mixed = classify("\t".join(record[protein_start:]))
            rows.append(
                {
                    "spec": record[index["SpecId"]],
                    "label": int(record[index["Label"]]),
                    "scan": int(record[index["ScanNr"]]),
                    "mass": record[index["ExpMass"]],
                    "peptide": record[index["Peptide"]],
                    "pure": pure,
                    "mixed": mixed,
                    "decoy": decoy,
                }
            )
            features.append([float(record[position]) for position in feature_positions])
    return rows, features, feature_names


def calibrate(winners, rows, q, fallback):
    """Entrapment accounting under three foreign-fraction assumptions."""
    targets = [i for i in winners if rows[i]["label"] > 0]
    decoys = [i for i in winners if rows[i]["label"] < 0]
    usable = [i for i in decoys if not rows[i]["mixed"]]
    whole_file_fraction = (
        sum(rows[i]["pure"] for i in usable) / len(usable) if usable else fallback
    )
    out = []
    for threshold in THRESHOLDS:
        accepted = [i for i in targets if q[i] < threshold]
        accepted_decoys = [i for i in usable if q[i] < threshold]
        local = (
            sum(rows[i]["pure"] for i in accepted_decoys) / len(accepted_decoys)
            if accepted_decoys
            else whole_file_fraction
        )
        entrapment = sum(rows[i]["pure"] for i in accepted)
        row = {
            "q_threshold": threshold,
            "accepted_targets": len(accepted),
            "pure_entrapment": entrapment,
            "accepted_decoys": len(accepted_decoys),
            "effective_entrapment_fraction": local,
            "whole_file_entrapment_fraction": whole_file_fraction,
        }
        for name, fraction in (
            ("adjusted_fdp_accepted_decoy_fraction", local),
            ("adjusted_fdp_whole_file_fraction", whole_file_fraction),
            ("adjusted_fdp_declared_fraction", fallback),
        ):
            row[name] = (
                entrapment / fraction / len(accepted) if accepted and fraction else None
            )
        out.append(row)
    return out


def load_result_dir(directory: Path):
    """percolator-rs / percolator PSM output for the rescored arm."""
    rows = []
    q = []
    for name, decoy in (("target.psms.tsv", False), ("decoy.psms.tsv", True)):
        path = directory / name
        with path.open(newline="") as handle:
            for record in csv.DictReader(handle, delimiter="\t"):
                _, pure, mixed = classify(
                    "\t".join(
                        v if isinstance(v, str) else "\t".join(v)
                        for k, v in record.items()
                        if k == "proteinIds" or k is None
                    )
                )
                rows.append(
                    {"label": -1 if decoy else 1, "pure": pure, "mixed": mixed}
                )
                q.append(float(record["q-value"]))
    return rows, q


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pins", type=Path, required=True)
    parser.add_argument("--rescored", type=Path, help="seed directory of a run_entrapment study")
    parser.add_argument("--entrapment-fraction", type=float, required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    # Same set the five-seed entrapment study uses: comet-out is a duplicate of
    # one of the six and is excluded there too.
    datasets = sorted(
        pin for pin in args.pins.glob("comet-*/comet.pin") if pin.parent.name != "comet-out"
    )
    if not datasets:
        parser.error(f"no comet-*/comet.pin under {args.pins}")

    report = {"seed": args.seed, "thresholds": list(THRESHOLDS), "datasets": [], "pooled": {}}
    pooled_raw = []
    pooled_rows = []
    pooled_q = []
    for path in datasets:
        name = path.parent.name
        rows, features, feature_names = read_pin(path)
        # Best single feature, chosen the way an initial direction is chosen:
        # highest competed target yield at q<0.01, either sign.
        best = None
        for position, feature_name in enumerate(feature_names):
            for sign in (1.0, -1.0):
                scores = [sign * row[position] for row in features]
                winners = compete(rows, scores, args.seed)
                subset_scores = [scores[i] for i in winners]
                subset_labels = [rows[i]["label"] for i in winners]
                qs = qvalues(subset_scores, subset_labels)
                count = sum(
                    1
                    for k, i in enumerate(winners)
                    if rows[i]["label"] > 0 and qs[k] < 0.01
                )
                if best is None or count > best[0]:
                    best = (count, feature_name, sign, winners, qs)
        count, feature_name, sign, winners, qs = best
        q_by_row = {row_index: qs[k] for k, row_index in enumerate(winners)}
        calibration = calibrate(winners, rows, q_by_row, args.entrapment_fraction)
        report["datasets"].append(
            {
                "dataset": name,
                "arm": "raw",
                "rows": len(rows),
                "best_feature": feature_name,
                "best_feature_sign": sign,
                "competed_winners": len(winners),
                "calibration": calibration,
            }
        )
        print(
            f"{name}\traw best_feature={sign:+.0f}*{feature_name}\t"
            + "\t".join(
                f"q<{row['q_threshold']}:n={row['accepted_targets']},fdp="
                + (
                    "n/a"
                    if row["adjusted_fdp_accepted_decoy_fraction"] is None
                    else f"{row['adjusted_fdp_accepted_decoy_fraction']:.4%}"
                )
                for row in calibration
                if row["q_threshold"] in (0.001, 0.01, 0.1)
            ),
            flush=True,
        )
        pooled_raw.append((rows, winners, q_by_row))

    # Pooled raw arm: the entrapment study pools all datasets before thresholding.
    offset = 0
    all_rows = []
    all_winners = []
    all_q = {}
    for rows, winners, q_by_row in pooled_raw:
        for index, row in enumerate(rows):
            all_rows.append(row)
        for index in winners:
            all_winners.append(index + offset)
            all_q[index + offset] = q_by_row[index]
        offset += len(rows)
    report["pooled"]["raw"] = calibrate(all_winners, all_rows, all_q, args.entrapment_fraction)
    print("POOLED raw", json.dumps(report["pooled"]["raw"], indent=1), flush=True)

    if args.rescored:
        directories = sorted(p.parent for p in args.rescored.glob("*/target.psms.tsv"))
        rows = []
        q = []
        for directory in directories:
            part_rows, part_q = load_result_dir(directory)
            rows.extend(part_rows)
            q.extend(part_q)
        winners = list(range(len(rows)))
        report["pooled"]["rescored"] = calibrate(
            winners, rows, dict(enumerate(q)), args.entrapment_fraction
        )
        print("POOLED rescored", json.dumps(report["pooled"]["rescored"], indent=1), flush=True)

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()
