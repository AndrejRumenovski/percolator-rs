#!/usr/bin/env python3
"""Fresh end-to-end adversarial fixtures for the final repair audit.

The script deliberately does not reuse the repository's regression fixtures.
It preserves every generated PIN, command, stderr log, and result table under a
caller-supplied work directory and writes one machine-readable synopsis.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import re
import shutil
import subprocess
from dataclasses import dataclass, replace
from pathlib import Path


HEADER = (
    "SpecId", "Label", "ScanNr", "ExpMass", "f0", "f1", "f2", "f3",
    "Peptide", "Proteins",
)
MASK = (1 << 64) - 1


@dataclass(frozen=True)
class Row:
    psm_id: str
    label: int
    scan: int
    mass: float
    features: tuple[float, float, float, float]
    peptide: str
    proteins: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_pin(path: Path, rows: list[Row]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        for row in rows:
            writer.writerow((
                row.psm_id, row.label, row.scan, row.mass, *row.features,
                row.peptide, row.proteins,
            ))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(command: list[str], root: Path, tag: str) -> dict:
    logs = root / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    execution = subprocess.run(command, text=True, capture_output=True, check=False)
    (logs / f"{tag}.stdout").write_text(execution.stdout)
    (logs / f"{tag}.stderr").write_text(execution.stderr)
    record = {
        "tag": tag,
        "argv": command,
        "exit_code": execution.returncode,
        "stdout": str(logs / f"{tag}.stdout"),
        "stderr": str(logs / f"{tag}.stderr"),
    }
    if execution.returncode:
        raise RuntimeError(f"{tag} failed ({execution.returncode}):\n{execution.stderr}")
    return record


def shuffled(rows: list[Row], seed: int) -> list[Row]:
    output = rows.copy()
    random.Random(seed).shuffle(output)
    return output


def arrange(rows: list[Row], mode: str, seed: int = 0) -> list[Row]:
    if mode == "original":
        return rows.copy()
    if mode == "reversed":
        return list(reversed(rows))
    if mode == "target-first":
        return sorted(rows, key=lambda row: row.label < 0)
    if mode == "decoy-first":
        return sorted(rows, key=lambda row: row.label > 0)
    if mode == "shuffle":
        return shuffled(rows, seed)
    raise ValueError(mode)


def combined_result(targets: Path, decoys: Path) -> dict[str, tuple[str, ...]]:
    output: dict[str, tuple[str, ...]] = {}
    for label, path in ((1, targets), (-1, decoys)):
        for row in read_rows(path):
            output[row["PSMId"]] = (
                str(label), row["score"], row["q-value"],
                row["posterior_error_prob"], row["peptide"], row["proteinIds"],
            )
    return output


def psm_run(
    binary: Path,
    root: Path,
    tag: str,
    pins: list[Path],
    *,
    seed: int,
    join: bool = False,
    maxiter: int = 0,
    extra: tuple[str, ...] = (),
) -> dict:
    output = root / "outputs" / tag
    output.mkdir(parents=True, exist_ok=True)
    targets = output / "targets.tsv"
    decoys = output / "decoys.tsv"
    command = [
        str(binary), "--canonical", "--no-select-c", "--seed", str(seed),
        "--num-threads", "1", "--maxiter", str(maxiter),
        "--results-psms", str(targets), "--decoy-results-psms", str(decoys),
        *extra,
    ]
    if join:
        command.append("--join")
    command.extend(str(pin) for pin in pins)
    execution = run(command, root, tag)
    rows = combined_result(targets, decoys)
    target_rows = read_rows(targets)
    return {
        "execution": execution,
        "rows": rows,
        "target_winners": len(target_rows),
        "decoy_winners": len(read_rows(decoys)),
        "targets_q_lt_0.01": sum(float(row["q-value"]) < 0.01 for row in target_rows),
        "minimum_target_q": min((float(row["q-value"]) for row in target_rows), default=None),
        "target_sha256": sha256(targets),
        "decoy_sha256": sha256(decoys),
    }


def single_file_tie_attack(binary: Path, root: Path) -> dict:
    base: list[Row] = []
    for scan in range(1, 234):
        level = float((scan % 17) - 8)
        features = (level, float(scan % 7), float(scan % 11), float(scan % 13))
        base.extend([
            Row(f"S{scan}_T", 1, scan, 500.0, features, f"K.S{scan}T.R", f"P{scan}"),
            Row(f"S{scan}_D", -1, scan, 500.0, features, f"K.S{scan}D.R", f"DECOY_P{scan}"),
        ])
    live = root / "single-file" / "fixture.pin"
    arms = [
        ("original", "original", 0),
        ("reversed", "reversed", 0),
        ("target-first", "target-first", 0),
        ("decoy-first", "decoy-first", 0),
        ("shuffle-1701", "shuffle", 1701),
        ("shuffle-9907", "shuffle", 9907),
    ]
    results = {}
    for name, mode, seed in arms:
        rows = arrange(base, mode, seed)
        preserved = root / "fixtures" / "single-file" / f"{name}.pin"
        write_pin(preserved, rows)
        write_pin(live, rows)
        results[name] = psm_run(binary, root, f"single-{name}", [live], seed=19)
    reference = results["original"]["rows"]
    winner_invariant = all(set(result["rows"]) == set(reference) for result in results.values())
    score_q_invariant = all(
        all(result["rows"][psm_id][1:3] == reference[psm_id][1:3] for psm_id in reference)
        for result in results.values()
    )
    target_pep_invariant = all(
        all(
            result["rows"][psm_id][3] == reference[psm_id][3]
            for psm_id, values in reference.items() if values[0] == "1"
        )
        for result in results.values()
    )
    decoy_pep_changes = {
        name: [
            psm_id for psm_id, values in reference.items()
            if values[0] == "-1" and result["rows"][psm_id][3] != values[3]
        ]
        for name, result in results.items()
    }
    return {
        "spectra": 233,
        "arms": {name: {key: value for key, value in result.items() if key != "rows"}
                 for name, result in results.items()},
        "winner_and_statistics_invariant": all(result["rows"] == reference for result in results.values()),
        "winner_identity_invariant": winner_invariant,
        "score_and_q_invariant": score_q_invariant,
        "target_pep_invariant": target_pep_invariant,
        "decoy_pep_invariant": all(not changed for changed in decoy_pep_changes.values()),
        "decoy_pep_changed_rows_by_arm": decoy_pep_changes,
        "distinct_complete_results": len({tuple(sorted(result["rows"].items())) for result in results.values()}),
    }


def joined_base() -> dict[str, list[Row]]:
    names = ["orchid.pin", "quartz.pin", "raven.pin", "willow.pin"]
    output: dict[str, list[Row]] = {name: [] for name in names}
    for source, name in enumerate(names):
        for scan in range(1, 82):
            unit = source * 10_000 + scan
            label = 1 if (scan + source) % 3 else -1
            signal = 1.3 * label + ((scan * 7 + source) % 19 - 9) / 7.0
            features = (signal, float(scan % 5), float((scan + source) % 9), float(scan % 13))
            output[name].append(Row(
                f"J{source}_{scan}", label, scan, 500.0, features,
                f"K.J{source}_{scan}.R", ("" if label > 0 else "DECOY_") + f"JP{unit}",
            ))
        # Fresh exact target/decoy ties and strict representable near-ties in every source.
        exact = (2.75 + source, 1.0, 2.0, 3.0)
        output[name].extend([
            Row(f"J{source}_EX_T", 1, 9001, 500.0, exact, f"K.J{source}EXT.R", f"JEX{source}"),
            Row(f"J{source}_EX_D", -1, 9001, 500.0, exact, f"K.J{source}EXD.R", f"DECOY_JEX{source}"),
        ])
        low = 3.5 + source
        output[name].extend([
            Row(f"J{source}_NEAR_T", 1, 9002, 500.0, (low, 0.0, 1.0, 2.0), f"K.J{source}NT.R", f"JNT{source}"),
            Row(
                f"J{source}_NEAR_D", -1, 9002, 500.0,
                (math.nextafter(low, math.inf), 0.0, 1.0, 2.0),
                f"K.J{source}ND.R", f"DECOY_JNT{source}",
            ),
        ])
    return output


def joined_permutation_attack(binary: Path, root: Path) -> dict:
    base = joined_base()
    names = list(base)
    live_root = root / "joined" / "live"
    arms = [
        ("original", names, "original", 0),
        ("reverse-all", list(reversed(names)), "reversed", 0),
        ("target-first", [names[2], names[0], names[3], names[1]], "target-first", 0),
        ("decoy-first", [names[1], names[3], names[0], names[2]], "decoy-first", 0),
        ("shuffle-31415", [names[3], names[0], names[2], names[1]], "shuffle", 31415),
        ("shuffle-27182", [names[1], names[2], names[3], names[0]], "shuffle", 27182),
    ]
    results = {}
    for arm, order, mode, seed in arms:
        live_paths = []
        for offset, name in enumerate(names):
            rows = arrange(base[name], mode, seed + offset)
            preserved = root / "fixtures" / "joined" / arm / name
            write_pin(preserved, rows)
            live = live_root / name
            write_pin(live, rows)
        for name in order:
            live_paths.append(live_root / name)
        results[arm] = psm_run(
            binary, root, f"joined-{arm}", live_paths, seed=73, join=True, maxiter=3,
        )
    reference = results["original"]["rows"]

    # The same four actual files are then referenced through symlinks whose
    # lexical order reverses their content-to-source-index assignment.
    aliases = root / "joined" / "aliases"
    aliases.mkdir(parents=True, exist_ok=True)
    alias_names = ["z_orchid.pin", "a_quartz.pin", "y_raven.pin", "b_willow.pin"]
    alias_paths = []
    for name, alias_name in zip(names, alias_names):
        alias = aliases / alias_name
        if alias.exists() or alias.is_symlink():
            alias.unlink()
        os.symlink(live_root / name, alias)
        alias_paths.append(alias)
    aliased = psm_run(binary, root, "joined-path-aliases", alias_paths, seed=73, join=True, maxiter=3)
    path_changed = sorted(
        psm_id for psm_id in reference
        if aliased["rows"].get(psm_id) != reference[psm_id]
    )
    cutoff = joined_path_alias_cutoff_attack(binary, root)
    return {
        "files": names,
        "arms": {name: {key: value for key, value in result.items() if key != "rows"}
                 for name, result in results.items()},
        "file_row_target_decoy_permutation_invariant": all(result["rows"] == reference for result in results.values()),
        "distinct_permutation_results": len({tuple(sorted(result["rows"].items())) for result in results.values()}),
        "path_alias": {
            **{key: value for key, value in aliased.items() if key != "rows"},
            "changed_psms": len(path_changed),
            "changed_examples": path_changed[:20],
            "invariant": not path_changed,
        },
        "path_alias_cutoff": cutoff,
    }


def joined_path_alias_cutoff_attack(binary: Path, root: Path) -> dict:
    """Minimize path-dependent source numbering into a 0-versus-101 boundary."""
    base = root / "joined-alias-cutoff" / "real"
    alpha = base / "a-alpha.pin"
    beta = base / "b-beta.pin"
    high = [
        Row(f"AHIGH{scan}", 1, scan, 500.0, (10.0, 0.0, 1.0, 2.0),
            f"K.AHIGH{scan}.R", f"AHIGH{scan}")
        for scan in range(1, 101)
    ]
    high.extend([
        Row("ABOUND_T", 1, 10_002, 500.0, (0.0, 0.0, 1.0, 2.0), "K.ABOUNDT.R", "ABOUND"),
        Row("ABOUND_D", -1, 10_002, 500.0, (0.0, 0.0, 1.0, 2.0), "K.ABOUNDD.R", "DECOY_ABOUND"),
    ])
    low = [
        Row(f"BLOW{scan}", -1, 20_000 + scan, 500.0, (-10.0, 0.0, 1.0, 2.0),
            f"K.BLOW{scan}.R", f"DECOY_BLOW{scan}")
        for scan in range(1, 101)
    ]
    write_pin(alpha, high)
    write_pin(beta, low)
    original = psm_run(binary, root, "joined-alias-cutoff-original", [alpha, beta], seed=1, join=True)

    aliases = root / "joined-alias-cutoff" / "aliases"
    aliases.mkdir(parents=True, exist_ok=True)
    alpha_alias = aliases / "z-alpha.pin"
    beta_alias = aliases / "a-beta.pin"
    os.symlink(alpha, alpha_alias)
    os.symlink(beta, beta_alias)
    aliased = psm_run(binary, root, "joined-alias-cutoff-aliased", [alpha_alias, beta_alias], seed=1, join=True)
    return {
        "original_boundary_winner": "ABOUND_T" if "ABOUND_T" in original["rows"] else "ABOUND_D",
        "aliased_boundary_winner": "ABOUND_T" if "ABOUND_T" in aliased["rows"] else "ABOUND_D",
        "original_targets_q_lt_0.01": original["targets_q_lt_0.01"],
        "aliased_targets_q_lt_0.01": aliased["targets_q_lt_0.01"],
        "changed_winner": set(original["rows"]) != set(aliased["rows"]),
    }


def candidate_multiplicity_attack(binary: Path, root: Path) -> dict:
    balanced: list[Row] = []
    duplicated: list[Row] = []
    for scan in range(1, 102):
        features = (0.0, 0.0, 0.0, 0.0)
        target = Row(f"M{scan}_T", 1, scan, 500.0, features, f"K.M{scan}T.R", f"MP{scan}")
        decoy = Row(f"M{scan}_D", -1, scan, 500.0, features, f"K.M{scan}D.R", f"DECOY_MP{scan}")
        balanced.extend([target, decoy])
        duplicated.extend([target] * 94)
        duplicated.append(decoy)
    balanced_pin = root / "fixtures" / "multiplicity" / "balanced.pin"
    duplicated_pin = root / "fixtures" / "multiplicity" / "target-duplicated-94x.pin"
    write_pin(balanced_pin, balanced)
    write_pin(duplicated_pin, duplicated)
    balanced_result = psm_run(binary, root, "multiplicity-balanced", [balanced_pin], seed=1)
    duplicated_result = psm_run(binary, root, "multiplicity-duplicated", [duplicated_pin], seed=1)
    return {
        "spectra": 101,
        "target_copies_per_spectrum": 94,
        "balanced": {key: value for key, value in balanced_result.items() if key != "rows"},
        "duplicated": {key: value for key, value in duplicated_result.items() if key != "rows"},
        "false_discoveries_created_by_exact_duplicate_rows": (
            duplicated_result["targets_q_lt_0.01"] - balanced_result["targets_q_lt_0.01"]
        ),
    }


def rng_next(state: int) -> int:
    state ^= (state << 13) & MASK
    state ^= state >> 7
    state ^= (state << 17) & MASK
    return state & MASK


def fold_map(weights: dict[int, int], seed: int) -> dict[int, int]:
    scans = sorted(weights)
    state = max(seed, 1)
    for index in range(len(scans) - 1, 0, -1):
        state = rng_next(state)
        other = state % (index + 1)
        scans[index], scans[other] = scans[other], scans[index]
    sizes = [0, 0, 0]
    output = {}
    for scan in scans:
        fold = min(range(3), key=lambda candidate: sizes[candidate])
        output[scan] = fold
        sizes[fold] += weights[scan]
    return output


def cv_fixture() -> list[Row]:
    generator = random.Random(0x20260827)
    rows = []
    for scan in range(411):
        label = 1 if generator.random() < 0.64 else -1
        signal = float(label)
        features = (
            signal * 1.1 + generator.gauss(0.0, 1.7),
            signal * -0.4 + generator.gauss(0.0, 1.1),
            generator.gauss(0.0, 1.0),
            generator.gauss(0.0, 1.0),
        )
        rows.append(Row(
            f"CVX{scan}", label, scan, 500.0, features,
            f"K.CVX{scan}.R", ("" if label > 0 else "DECOY_") + f"CVP{scan}",
        ))
    return rows


def parse_selections(stderr_path: Path) -> list[str]:
    text = stderr_path.read_text()
    patterns = [
        r"fold (\d+): Cpos=([0-9.]+) Cneg=([0-9.]+)",
        r"fold (\d+): C=([0-9.]+), class-weights=([0-9.]+):([0-9.]+), features=(\d+), tolerance=([0-9.eE+-]+)",
    ]
    values = []
    for pattern in patterns:
        values.extend("|".join(match) for match in re.findall(pattern, text))
    return values


def cv_paths(root: Path, mode: str, tag: str) -> list[Path]:
    count = 3 if mode == "ensemble" else 1
    return [root / "cv" / f"{mode}-{tag}-engine{engine}.pin" for engine in range(count)]


def engine_keeps(scan: int, engine: int) -> bool:
    return engine == 0 or (engine == 1 and scan % 2 == 0) or (engine == 2 and scan % 3 == 0)


def write_cv_inputs(root: Path, mode: str, tag: str, rows: list[Row]) -> list[Path]:
    paths = cv_paths(root, mode, tag)
    for engine, path in enumerate(paths):
        selected = rows if mode != "ensemble" else [row for row in rows if engine_keeps(row.scan, engine)]
        write_pin(path, selected)
    return paths


def cv_run(binary: Path, root: Path, mode: str, tag: str, inputs: list[Path]) -> dict:
    output = root / "outputs" / f"cv-{mode}-{tag}"
    output.mkdir(parents=True, exist_ok=True)
    targets, decoys = output / "targets.tsv", output / "decoys.tsv"
    command = [
        str(binary), "--canonical", "--seed", "47", "--num-threads", "1",
        "--maxiter", "3", "--no-psm-competition",
        "--results-psms", str(targets), "--decoy-results-psms", str(decoys),
    ]
    if mode == "select-c":
        command.append("--select-c")
        command.append(str(inputs[0]))
    elif mode == "ensemble":
        command.extend(["--no-select-c", "--ensemble"])
        command.extend(f"fresh{engine}={path}" for engine, path in enumerate(inputs))
    else:
        command.extend(["--no-select-c", str(inputs[0])])
    execution = run(command, root, f"cv-{mode}-{tag}")
    scores = {psm_id.rsplit(":", 1)[-1]: values[1]
              for psm_id, values in combined_result(targets, decoys).items()}
    return {
        "execution": execution,
        "scores": scores,
        "selections": parse_selections(Path(execution["stderr"])),
    }


def cv_isolation_attack(binary: Path, root: Path) -> dict:
    base = cv_fixture()
    results = {}
    for mode in ("fixed-c", "select-c", "ensemble"):
        weights = {
            row.scan: (sum(engine_keeps(row.scan, engine) for engine in range(3)) if mode == "ensemble" else 1)
            for row in base
        }
        folds = fold_map(weights, 47)
        attacked = {scan for scan, fold in folds.items() if fold == 1}
        sentinel = set(sorted(attacked)[:11])
        variants = {
            "clean": base,
            "labels": [replace(row, label=-row.label) if row.scan in attacked else row for row in base],
            "outliers": [
                replace(row, features=(1e12, -1e12, 5e11, -5e11))
                if row.scan in attacked - sentinel else row
                for row in base
            ],
            "reversed": list(reversed(base)),
            "shuffle": shuffled(base, 884422),
        }
        runs = {}
        for tag, rows in variants.items():
            inputs = write_cv_inputs(root, mode, tag, rows)
            for engine, path in enumerate(inputs):
                preserved = root / "fixtures" / "cv" / mode / tag / f"engine{engine}.pin"
                preserved.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, preserved)
            runs[tag] = cv_run(binary, root, mode, tag, inputs)
        clean = runs["clean"]
        attacked_ids = [f"CVX{scan}" for scan in attacked]
        sentinel_ids = [f"CVX{scan}" for scan in sentinel]
        label_changed = [psm_id for psm_id in attacked_ids if runs["labels"]["scores"].get(psm_id) != clean["scores"].get(psm_id)]
        outlier_changed = [psm_id for psm_id in sentinel_ids if runs["outliers"]["scores"].get(psm_id) != clean["scores"].get(psm_id)]
        reversed_changed = [psm_id for psm_id, score in clean["scores"].items() if runs["reversed"]["scores"].get(psm_id) != score]
        shuffle_changed = [psm_id for psm_id, score in clean["scores"].items() if runs["shuffle"]["scores"].get(psm_id) != score]
        fold_selection = lambda run_result: [value for value in run_result["selections"] if value.startswith("1|")]
        results[mode] = {
            "attacked_fold": 1,
            "attacked_spectra": len(attacked),
            "sentinels": len(sentinel),
            "heldout_label_changed_scores": len(label_changed),
            "heldout_outlier_changed_sentinel_scores": len(outlier_changed),
            "reversed_row_changed_scores": len(reversed_changed),
            "shuffled_row_changed_scores": len(shuffle_changed),
            "selection_clean": fold_selection(clean),
            "selection_labels": fold_selection(runs["labels"]),
            "selection_outliers": fold_selection(runs["outliers"]),
            "selection_invariant": (
                fold_selection(clean) == fold_selection(runs["labels"])
                == fold_selection(runs["outliers"])
            ),
            "leakage_detected": bool(label_changed or outlier_changed),
            "row_permutation_invariant_at_printed_precision": not reversed_changed and not shuffle_changed,
        }
    return results


def protein_pipeline_run(
    binary: Path,
    root: Path,
    tag: str,
    pin: Path,
    seed: int,
    *,
    competition: bool,
    bayesian: bool = False,
) -> dict:
    output = root / "outputs" / tag
    output.mkdir(parents=True, exist_ok=True)
    targets, decoys = output / "targets.tsv", output / "decoys.tsv"
    proteins, decoy_proteins = output / "proteins.tsv", output / "decoy-proteins.tsv"
    command = [
        str(binary), "--canonical", "--no-select-c", "--maxiter", "0", "--seed", str(seed),
        "--num-threads", "1", "--results-psms", str(targets),
        "--decoy-results-psms", str(decoys), "--results-proteins", str(proteins),
        "--decoy-results-proteins", str(decoy_proteins),
    ]
    if not competition:
        command.append("--no-psm-competition")
    if bayesian:
        command.extend(["--protein-inference", "bayesian"])
    command.append(str(pin))
    execution = run(command, root, tag)
    protein_rows = read_rows(proteins) + read_rows(decoy_proteins)
    return {
        "execution": execution,
        "psm_ids": sorted(combined_result(targets, decoys)),
        "groups": sorted((row["ProteinGroupId"], row["proteinIds"], row["posterior_error_prob"], row["q-value"]) for row in protein_rows),
        "picked_pep_all_na": all(row["posterior_error_prob"] == "NA" for row in protein_rows),
        "bayesian_pep_all_numeric": all(row["posterior_error_prob"] != "NA" for row in protein_rows),
    }


def protein_repaired_area_attack(binary: Path, root: Path) -> dict:
    rows = [
        Row("LOSS_A", 1, 9001, 500.0, (5.0, 0.0, 1.0, 2.0), "K.LOSTMAP.R", "LOSS_A_PROT"),
        Row("LOSS_B", 1, 9001, 500.0, (5.0, 0.0, 1.0, 2.0), "K.LOSTMAP.R", "LOSS_B_PROT"),
        Row("MIXED", 1, 9002, 500.0, (4.0, 0.0, 1.0, 2.0), "K.MIXED.R", "PAIRX DECOY_PAIRX"),
    ]
    for scan in range(1, 241):
        label = 1 if scan % 3 else -1
        rows.append(Row(
            f"PBG{scan}", label, scan, 500.0,
            (label * 1.2 + (scan % 7) / 10.0, float(scan % 5), float(scan % 11), float(scan % 13)),
            f"K.PBG{scan}.R", ("" if label > 0 else "DECOY_") + f"PBG{scan}",
        ))
    original = root / "fixtures" / "protein" / "original.pin"
    reversed_pin = root / "fixtures" / "protein" / "reversed.pin"
    write_pin(original, rows)
    write_pin(reversed_pin, list(reversed(rows)))
    no_comp = protein_pipeline_run(binary, root, "protein-no-competition", original, 1, competition=False)
    seed1 = protein_pipeline_run(binary, root, "protein-competition-seed1", original, 1, competition=True)
    seed3 = protein_pipeline_run(binary, root, "protein-competition-seed3", original, 3, competition=True)
    reversed_seed1 = protein_pipeline_run(binary, root, "protein-reversed-seed1", reversed_pin, 1, competition=True)
    bayesian = protein_pipeline_run(binary, root, "protein-bayesian", original, 1, competition=False, bayesian=True)
    group_has = lambda result, protein: any(protein in row[1].split(",") for row in result["groups"])
    return {
        "no_competition_union_contains_both": group_has(no_comp, "LOSS_A_PROT") and group_has(no_comp, "LOSS_B_PROT"),
        "competition_seed1": {
            "loss_a_present": group_has(seed1, "LOSS_A_PROT"),
            "loss_b_present": group_has(seed1, "LOSS_B_PROT"),
            "picked_pep_all_na": seed1["picked_pep_all_na"],
        },
        "competition_seed3": {
            "loss_a_present": group_has(seed3, "LOSS_A_PROT"),
            "loss_b_present": group_has(seed3, "LOSS_B_PROT"),
            "picked_pep_all_na": seed3["picked_pep_all_na"],
        },
        "seed_changes_protein_mapping": seed1["groups"] != seed3["groups"],
        "reversed_insertion_invariant_for_seed1": seed1["groups"] == reversed_seed1["groups"],
        "mixed_target_decoy_not_co_grouped": all(not ("PAIRX" in row[1].split(",") and "DECOY_PAIRX" in row[1].split(",")) for row in no_comp["groups"]),
        "picked_pep_is_na": no_comp["picked_pep_all_na"],
        "bayesian_pep_is_numeric": bayesian["bayesian_pep_all_numeric"],
        "raw": {
            "no_competition_groups": no_comp["groups"],
            "seed1_groups": seed1["groups"],
            "seed3_groups": seed3["groups"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"binary not found: {binary}")
    if args.work_dir.exists():
        parser.error(f"work directory already exists: {args.work_dir}")
    if args.json.exists():
        parser.error(f"JSON already exists: {args.json}")
    args.work_dir.mkdir(parents=True)

    result = {
        "schema_version": 1,
        "audit": "fresh final joined-input/protein/CV adversarial audit",
        "binary": str(binary),
        "binary_sha256": sha256(binary),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "single_file_ties": single_file_tie_attack(binary, args.work_dir),
        "joined_inputs": joined_permutation_attack(binary, args.work_dir),
        "candidate_multiplicity": candidate_multiplicity_attack(binary, args.work_dir),
        "cv_isolation": cv_isolation_attack(binary, args.work_dir),
        "protein_repaired_areas": protein_repaired_area_attack(binary, args.work_dir),
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "single_file_invariant": result["single_file_ties"]["winner_and_statistics_invariant"],
        "joined_permutation_invariant": result["joined_inputs"]["file_row_target_decoy_permutation_invariant"],
        "joined_path_alias_invariant": result["joined_inputs"]["path_alias"]["invariant"],
        "duplicate_false_discoveries": result["candidate_multiplicity"]["false_discoveries_created_by_exact_duplicate_rows"],
        "cv": result["cv_isolation"],
        "protein_seed_changes_mapping": result["protein_repaired_areas"]["seed_changes_protein_mapping"],
    }, indent=2))


if __name__ == "__main__":
    main()
