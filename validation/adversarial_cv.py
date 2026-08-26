#!/usr/bin/env python3
"""End-to-end held-out-label attacks on every advertised cross-validation mode.

Three independent attacks, reported per mode.  Nothing is inferred about one
mode from another.

``fold``
    The original attack.  The fixture and fold assignment are deterministic; the
    second input changes only the labels of outer fold 0 and every feature is
    byte-identical.  The model that scores fold 0 must give it the same scores,
    and — where the mode selects hyperparameters per fold — fold 0 must keep the
    same selection.

``row``
    Flip one row's label, leave the file otherwise byte-identical, and read back
    that row's own score.  A row is held out by exactly one model, and that model
    trained without it, so its own label cannot reach its own score.  This attack
    needs no knowledge of the fold assignment.

``ensemble``
    The same per-row attack across two ENGINE=PIN inputs that overlap on scan
    numbers, which is where the cross-engine agreement features are built.

Run against two builds to see a repair:

    python3 validation/adversarial_cv.py --binary OLD --json old.json
    python3 validation/adversarial_cv.py --binary NEW --json new.json
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import re
import subprocess
import tempfile
from pathlib import Path


MASK = (1 << 64) - 1
ROWS = 600
# Predeclared, spread across the file so the deterministic fold assignment
# cannot place them all in one fold.
ATTACKED_ROWS = (3, 97, 198, 301, 444, 577)


def rng_next(state: int) -> int:
    state ^= (state << 13) & MASK
    state ^= state >> 7
    state ^= (state << 17) & MASK
    return state & MASK


def outer_folds(rows: int, seed: int = 1) -> list[int]:
    order = list(range(rows))
    state = max(seed, 1)
    for index in range(rows - 1, 0, -1):
        state = rng_next(state)
        other = state % (index + 1)
        order[index], order[other] = order[other], order[index]
    sizes = [0, 0, 0]
    folds = [0] * rows
    for row in order:
        fold = min(range(3), key=lambda candidate: sizes[candidate])
        sizes[fold] += 1
        folds[row] = fold
    return folds


def fixture(rows: int = ROWS) -> tuple[list[tuple[int, list[float]]], list[int]]:
    # This is a predeclared deterministic case, found by enumerating seeds 1..4.
    generator = random.Random(4)
    strength = generator.uniform(0.3, 2.5)
    noise = generator.uniform(0.5, 3.0)
    target_rate = generator.uniform(0.55, 0.8)
    output = []
    for _ in range(rows):
        label = 1 if generator.random() < target_rate else -1
        features = [
            label * strength + generator.gauss(0, noise),
            label * generator.uniform(-strength, strength) + generator.gauss(0, noise),
            generator.gauss(0, 1),
            generator.gauss(0, 1),
        ]
        output.append((label, features))
    return output, outer_folds(rows)


def write_pin(
    path: Path,
    rows: list[tuple[int, list[float]]],
    flip_rows: set[int],
    keep: set[int] | None = None,
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("SpecId", "Label", "ScanNr", "ExpMass", "f0", "f1", "f2", "f3", "Peptide", "Proteins")
        )
        for index, (label, features) in enumerate(rows):
            if keep is not None and index not in keep:
                continue
            if index in flip_rows:
                label = -label
            writer.writerow(
                (f"p{index}", label, index, 500.0, *features, f"K.P{index}.R", f"P{index}")
            )


def load_scores(*paths: Path) -> dict[str, str]:
    """Scores as printed, so the comparison is exact and format-independent."""
    scores: dict[str, str] = {}
    for path in paths:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = row["PSMId"].rsplit(":", 1)[-1]
                scores[key] = row["score"]
    return scores


def selected_weights(stderr: str) -> list[str]:
    """Every class-weight pair the run reported, in fold order where per-fold."""
    per_fold = re.findall(r"fold (\d+): Cpos=([0-9.]+) Cneg=([0-9.]+)", stderr)
    if per_fold:
        return [f"fold{fold}={pos}:{neg}" for fold, pos, neg in per_fold]
    nested = re.findall(r"fold (\d+): C=([0-9.]+), class-weights=([0-9.]+):([0-9.]+)", stderr)
    if nested:
        return [f"fold{fold}={c}*{pos}:{neg}" for fold, c, pos, neg in nested]
    single = re.search(r"Cpos=([0-9.]+) Cneg=([0-9.]+)", stderr)
    if single:
        return [f"all={single.group(1)}:{single.group(2)}"]
    return []


def execute(
    binary: Path, root: Path, mode: str, tag: str, inputs: list[Path]
) -> tuple[list[str], dict[str, str]]:
    targets = root / f"{tag}.targets.tsv"
    decoys = root / f"{tag}.decoys.tsv"
    command = [
        str(binary),
        "--canonical",
        "--seed",
        "1",
        "--num-threads",
        "1",
        "--no-psm-competition",
        "--results-psms",
        str(targets),
        "--decoy-results-psms",
        str(decoys),
    ]
    if mode == "select-c":
        command.append("--select-c")
    else:
        command.append("--no-select-c")
    if mode == "ensemble":
        command.append("--ensemble")
        command += [f"engine{index}={path}" for index, path in enumerate(inputs)]
    else:
        command.append(str(inputs[0]))
    execution = subprocess.run(command, text=True, capture_output=True, check=False)
    if execution.returncode:
        raise RuntimeError(f"{mode}/{tag} failed:\n{execution.stderr}")
    return selected_weights(execution.stderr), load_scores(targets, decoys)


def ensemble_halves() -> list[set[int]]:
    """Two engines whose reported spectra overlap on the middle third."""
    return [set(range(0, ROWS * 2 // 3)), set(range(ROWS // 3, ROWS))]


def ensemble_outer_folds(seed: int = 1) -> list[int]:
    """Fold per source row for the merged ensemble dataset.

    Ensemble input groups by ScanNr alone, because the same spectrum is
    deliberately reported by several engines, and each scan then carries one or
    two candidate rows.  Folds are dealt to whichever is smallest by candidate
    count, so a plain equal-size split is the wrong map here — using it makes an
    unrelated part of the file look attacked and reports leakage that is not
    there.
    """
    halves = ensemble_halves()
    spectra = sorted({row for half in halves for row in half})
    sizes_of = {row: sum(row in half for half in halves) for row in spectra}
    order = list(spectra)
    state = max(seed, 1)
    for index in range(len(order) - 1, 0, -1):
        state = rng_next(state)
        other = state % (index + 1)
        order[index], order[other] = order[other], order[index]
    sizes = [0, 0, 0]
    folds = [-1] * ROWS
    for row in order:
        fold = min(range(3), key=lambda candidate: sizes[candidate])
        sizes[fold] += sizes_of[row]
        folds[row] = fold
    return folds


def paths_for(root: Path, mode: str, tag: str) -> list[Path]:
    if mode == "ensemble":
        return [root / f"{mode}.{tag}.{half}.pin" for half in range(2)]
    return [root / f"{mode}.{tag}.pin"]


def write_inputs(
    root: Path,
    mode: str,
    tag: str,
    rows: list[tuple[int, list[float]]],
    flip_rows: set[int],
) -> list[Path]:
    paths = paths_for(root, mode, tag)
    if mode == "ensemble":
        for path, keep in zip(paths, ensemble_halves()):
            write_pin(path, rows, flip_rows, keep)
    else:
        write_pin(paths[0], rows, flip_rows)
    return paths


def fold_attack(binary: Path, root: Path, mode: str, rows, folds) -> dict:
    clean_paths = write_inputs(root, mode, "fold-clean", rows, set())
    dirty_rows = {index for index, fold in enumerate(folds) if fold == 0}
    dirty_paths = write_inputs(root, mode, "fold-dirty", rows, dirty_rows)
    clean_weights, clean_scores = execute(binary, root, mode, f"{mode}.fold-clean", clean_paths)
    dirty_weights, dirty_scores = execute(binary, root, mode, f"{mode}.fold-dirty", dirty_paths)
    heldout = [index for index in dirty_rows if f"p{index}" in clean_scores]
    changed = [
        index
        for index in heldout
        if clean_scores[f"p{index}"] != dirty_scores.get(f"p{index}")
    ]
    fold0_clean = [value for value in clean_weights if value.startswith(("fold0", "all"))]
    fold0_dirty = [value for value in dirty_weights if value.startswith(("fold0", "all"))]
    return {
        "attack": "fold",
        "mode": mode,
        "heldout_rows": len(heldout),
        "heldout_scores_changed": len(changed),
        "fold0_selection_clean": fold0_clean,
        "fold0_selection_dirty": fold0_dirty,
        "fold0_selection_changed": fold0_clean != fold0_dirty,
        "all_selection_clean": clean_weights,
        "all_selection_dirty": dirty_weights,
        "leakage_detected": bool(changed) or fold0_clean != fold0_dirty,
    }


def row_attack(binary: Path, root: Path, mode: str, rows) -> dict:
    clean_paths = write_inputs(root, mode, "row-clean", rows, set())
    _, clean_scores = execute(binary, root, mode, f"{mode}.row-clean", clean_paths)
    changed = []
    for attacked in ATTACKED_ROWS:
        key = f"p{attacked}"
        if key not in clean_scores:
            continue
        paths = write_inputs(root, mode, f"row-dirty-{attacked}", rows, {attacked})
        _, scores = execute(binary, root, mode, f"{mode}.row-dirty-{attacked}", paths)
        if scores.get(key) != clean_scores[key]:
            changed.append(
                {"row": attacked, "clean": clean_scores[key], "dirty": scores.get(key)}
            )
    return {
        "attack": "row",
        "mode": mode,
        "rows_attacked": len(ATTACKED_ROWS),
        "own_scores_changed": len(changed),
        "detail": changed,
        "leakage_detected": bool(changed),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=Path("target/release/percolator-rs"))
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"binary not found: {binary}")

    rows, folds = fixture()
    results = []
    with tempfile.TemporaryDirectory(prefix="percolator-cv-audit-") as temporary:
        root = Path(temporary)
        for mode in ("fixed-c", "select-c", "ensemble"):
            mode_folds = ensemble_outer_folds() if mode == "ensemble" else folds
            results.append(fold_attack(binary, root, mode, rows, mode_folds))
            results.append(row_attack(binary, root, mode, rows))

    for result in results:
        print(
            "\t".join(
                f"{key}={value}"
                for key, value in result.items()
                if key not in {"detail", "all_selection_clean", "all_selection_dirty"}
            )
        )
    if args.json:
        args.json.write_text(json.dumps({"binary": str(binary), "results": results}, indent=1))


if __name__ == "__main__":
    main()
