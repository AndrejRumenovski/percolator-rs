#!/usr/bin/env python3
"""Independent post-repair attacks not present in the repair author's suite.

The script exercises three surfaces end to end through the frozen CLI:

* nonconstant partial/wide PSM ties under scientifically equivalent row orders;
* joined-input order, which should not change the per-file tie lottery;
* whole-held-out-fold label, outlier, and row-order attacks for fixed-C,
  ``--select-c``, and ``--ensemble``;
* equal-score peptide representatives carrying different protein mappings.

It writes a single machine-readable JSON file.  Production source is never
modified.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.10)
HEADER = ("SpecId", "Label", "ScanNr", "ExpMass", "f0", "f1", "f2", "Peptide", "Proteins")
MASK = (1 << 64) - 1


def write_pin(path: Path, rows: list[tuple]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def execute(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}):\n{' '.join(command)}\n{result.stderr}")
    return result


def xorshift(state: int) -> int:
    state ^= (state << 13) & MASK
    state ^= state >> 7
    state ^= (state << 17) & MASK
    return state & MASK


def spectrum_folds(spectrum_sizes: dict[int, int], seed: int = 1) -> dict[int, int]:
    """Independent reimplementation of the documented greedy spectrum split."""
    spectra = sorted(spectrum_sizes)
    state = max(seed, 1)
    for index in range(len(spectra) - 1, 0, -1):
        state = xorshift(state)
        other = state % (index + 1)
        spectra[index], spectra[other] = spectra[other], spectra[index]
    sizes = [0, 0, 0]
    folds = {}
    for scan in spectra:
        target = min(range(3), key=lambda fold: sizes[fold])
        sizes[target] += spectrum_sizes[scan]
        folds[scan] = target
    return folds


def shuffle_rows(rows: list[tuple], seed: int) -> list[tuple]:
    output = list(rows)
    random.Random(seed).shuffle(output)
    return output


def psm_command(binary: Path, pin: Path, targets: Path, decoys: Path, *, seed: int = 17) -> list[str]:
    return [
        str(binary), "--canonical", "--no-select-c", "--seed", str(seed),
        "--num-threads", "1", "--results-psms", str(targets),
        "--decoy-results-psms", str(decoys), str(pin),
    ]


def summarize_psms(targets: Path, decoys: Path) -> dict:
    target_rows = read_rows(targets)
    decoy_rows = read_rows(decoys)
    combined = target_rows + decoy_rows
    return {
        "target_winners": len(target_rows),
        "decoy_winners": len(decoy_rows),
        "counts": {
            f"q_lt_{threshold:g}": sum(float(row["q-value"]) < threshold for row in target_rows)
            for threshold in THRESHOLDS
        },
        "winner_ids": sorted(row["PSMId"] for row in combined),
        "winner_peptides": sorted(row["peptide"] for row in combined),
        "score_by_id": {row["PSMId"]: row["score"] for row in combined},
        "q_by_id": {row["PSMId"]: row["q-value"] for row in combined},
        "pep_by_id": {row["PSMId"]: row["posterior_error_prob"] for row in combined},
    }


def variable_tie_rows(spectra: int = 311) -> list[tuple]:
    """New nonconstant fixture with two-, four-, and five-way partial top ties."""
    rows = []
    for scan in range(1, spectra + 1):
        level = float((scan * 37) % 29) / 7.0
        width = 2 + scan % 4
        for candidate in range(width):
            label = 1 if candidate % 2 == 0 else -1
            rows.append((
                f"V{scan}_{candidate}", label, scan, 400.0 + (scan % 3),
                level, float(scan % 11), float((scan * 13) % 17),
                f"K.VAR{scan}_{candidate}.R",
                f"{'P' if label > 0 else 'DECOY_P'}{scan}_{candidate}",
            ))
        # A strict lower candidate makes these partial rather than whole-spectrum ties.
        label = -1 if scan % 2 else 1
        rows.append((
            f"V{scan}_LOW", label, scan, 400.0 + (scan % 3),
            level - 3.0, float(scan % 11), float((scan * 13) % 17),
            f"K.LOW{scan}.R", f"{'P' if label > 0 else 'DECOY_P'}{scan}_LOW",
        ))
    return rows


def tied_fixture_attacks(binary: Path, root: Path) -> dict:
    rows = variable_tie_rows()
    arms = {
        "canonical": list(rows),
        "reversed": list(reversed(rows)),
        "pair_blocks_reversed": [row for scan in range(311, 0, -1) for row in rows if row[2] == scan],
        **{f"shuffle_{seed}": shuffle_rows(rows, seed) for seed in (101, 202, 303, 404, 505)},
    }
    summaries = {}
    for arm, arranged in arms.items():
        pin = root / f"ties-{arm}.pin"
        target = root / f"ties-{arm}.target.tsv"
        decoy = root / f"ties-{arm}.decoy.tsv"
        write_pin(pin, arranged)
        execute(psm_command(binary, pin, target, decoy))
        summaries[arm] = summarize_psms(target, decoy)
    reference = summaries["canonical"]
    return {
        "fixture": "311 spectra, variable score levels, 2-5 tied top candidates and one strict lower candidate",
        "arms": summaries,
        "winner_invariant": all(value["winner_ids"] == reference["winner_ids"] for value in summaries.values()),
        "q_invariant": all(value["q_by_id"] == reference["q_by_id"] for value in summaries.values()),
        "pep_invariant": all(value["pep_by_id"] == reference["pep_by_id"] for value in summaries.values()),
    }


def all_tied_rows(prefix: str, scan_start: int, spectra: int) -> list[tuple]:
    rows = []
    for offset in range(spectra):
        scan = scan_start + offset
        rows.extend([
            (f"{prefix}_T_{scan}", 1, scan, 500.0, 0.0, 0.0, 0.0, f"K.{prefix}T{scan}.R", f"{prefix}_P{scan}"),
            (f"{prefix}_D_{scan}", -1, scan, 500.0, 0.0, 0.0, 0.0, f"K.{prefix}D{scan}.R", f"DECOY_{prefix}_P{scan}"),
        ])
    return rows


def joined_order_attack(binary: Path, root: Path) -> dict:
    # Different scan ranges make swapping numeric source ids change the key, not
    # merely exchange two identical key sets between files.
    first = root / "joined-alpha.pin"
    second = root / "joined-beta.pin"
    write_pin(first, all_tied_rows("ALPHA", 1, 137))
    write_pin(second, all_tied_rows("BETA", 10_001, 149))

    results = {}
    for name, inputs in (("alpha_beta", (first, second)), ("beta_alpha", (second, first))):
        targets = root / f"join-{name}.target.tsv"
        decoys = root / f"join-{name}.decoy.tsv"
        command = [
            str(binary), "--canonical", "--no-select-c", "--join", "--seed", "17",
            "--num-threads", "1", "--results-psms", str(targets),
            "--decoy-results-psms", str(decoys), *(str(path) for path in inputs),
        ]
        execute(command)
        results[name] = summarize_psms(targets, decoys)
    # Boundary attack: 100 strict high-scoring targets give q=0.01 exactly and
    # therefore no strict-q<0.01 discoveries.  The next spectrum is an exact
    # target/decoy tie.  Its label is keyed by numeric source position, so swapping
    # input files can turn that next winner into a target (101 targets, q<.01) or
    # a decoy (still no discoveries).  Low-scoring decoys stabilize training and
    # cannot improve the prefix's q-value.
    boundary_alpha = root / "boundary-alpha.pin"
    boundary_beta = root / "boundary-beta.pin"
    boundary_a_rows = [
        (f"HIGH_T_{scan}", 1, scan, 500.0, 10.0, 0.0, 0.0, f"K.HIGHT{scan}.R", f"HIGH_P{scan}")
        for scan in range(1, 101)
    ] + [
        ("BOUNDARY_T", 1, 10_002, 500.0, 0.0, 0.0, 0.0, "K.BOUNDARYT.R", "BOUNDARY_P"),
        ("BOUNDARY_D", -1, 10_002, 500.0, 0.0, 0.0, 0.0, "K.BOUNDARYD.R", "DECOY_BOUNDARY_P"),
    ]
    boundary_b_rows = [
        (f"LOW_D_{scan}", -1, scan, 500.0, -10.0, 0.0, 0.0, f"K.LOWD{scan}.R", f"DECOY_LOW_P{scan}")
        for scan in range(20_001, 20_101)
    ]
    write_pin(boundary_alpha, boundary_a_rows)
    write_pin(boundary_beta, boundary_b_rows)
    boundary_results = {}
    for name, inputs in (("alpha_beta", (boundary_alpha, boundary_beta)), ("beta_alpha", (boundary_beta, boundary_alpha))):
        targets = root / f"boundary-{name}.target.tsv"
        decoys = root / f"boundary-{name}.decoy.tsv"
        execute([
            str(binary), "--canonical", "--no-select-c", "--maxiter", "0", "--join",
            "--seed", "17", "--num-threads", "1", "--results-psms", str(targets),
            "--decoy-results-psms", str(decoys), *(str(path) for path in inputs),
        ])
        boundary_results[name] = summarize_psms(targets, decoys)

    return {
        "arms": results,
        "winner_invariant": results["alpha_beta"]["winner_ids"] == results["beta_alpha"]["winner_ids"],
        "counts_invariant": results["alpha_beta"]["counts"] == results["beta_alpha"]["counts"],
        "changed_winners": len(set(results["alpha_beta"]["winner_ids"]) ^ set(results["beta_alpha"]["winner_ids"])),
        "boundary_arms": boundary_results,
        "boundary_q01_invariant": boundary_results["alpha_beta"]["counts"]["q_lt_0.01"] == boundary_results["beta_alpha"]["counts"]["q_lt_0.01"],
    }


def cv_rows(spectra: int = 360, engine: int | None = None) -> list[tuple]:
    generator = random.Random(91 + (engine or 0))
    rows = []
    for scan in range(spectra):
        labels = (1, -1) if engine is None else (1 if (scan + engine) % 3 else -1,)
        for candidate, label in enumerate(labels):
            signal = float(label)
            peptide = f"K.CV{scan}_{candidate if engine is None else engine}.R"
            rows.append((
                f"CV{scan}_{candidate if engine is None else engine}", label, scan, 700.0,
                signal * 1.1 + generator.gauss(0.0, 1.3),
                signal * 0.4 + generator.gauss(0.0, 0.9),
                generator.gauss(0.0, 1.0), peptide,
                f"{'P' if label > 0 else 'DECOY_P'}{scan}_{candidate}",
            ))
    return rows


def write_cv_variant(
    path: Path,
    rows: list[tuple],
    folds: dict[int, int],
    variant: str,
    sentinel_scans: set[int],
) -> None:
    output = []
    for row in rows:
        values = list(row)
        if folds[row[2]] == 0:
            if variant == "labels":
                values[1] = -values[1]
            elif variant == "outliers" and row[2] not in sentinel_scans:
                values[4] = 1e12 + row[2]
                values[5] = -1e12 - row[2]
                values[6] = 5e11
        output.append(tuple(values))
    if variant == "reordered":
        training = [row for row in output if folds[row[2]] != 0]
        heldout = [row for row in output if folds[row[2]] == 0]
        output = training + list(reversed(heldout))
    write_pin(path, output)


def parse_selections(stderr: str) -> list[str]:
    selected = re.findall(
        r"fold (\d+): C=([0-9.]+), class-weights=([0-9.]+):([0-9.]+), "
        r"features=(\d+), tolerance=([^,]+), inner-q01-yield=(\d+)",
        stderr,
    )
    return [f"{fold}:{c}:{positive}:{negative}:{features}:{tolerance}:{inner}" for fold, c, positive, negative, features, tolerance, inner in selected]


def read_scores(*paths: Path) -> dict[str, str]:
    output = {}
    for path in paths:
        for row in read_rows(path):
            output[row["PSMId"]] = row["score"]
    return output


def run_cv(binary: Path, root: Path, mode: str, tag: str, inputs: list[Path]) -> tuple[dict[str, str], list[str]]:
    targets = root / f"cv-{mode}-{tag}.target.tsv"
    decoys = root / f"cv-{mode}-{tag}.decoy.tsv"
    command = [
        str(binary), "--canonical", "--seed", "1", "--num-threads", "1",
        "--no-psm-competition", "--results-psms", str(targets),
        "--decoy-results-psms", str(decoys),
    ]
    command.append("--select-c" if mode == "select-c" else "--no-select-c")
    if mode == "ensemble":
        command.append("--ensemble")
        command.extend(f"engine{index}={path}" for index, path in enumerate(inputs))
    else:
        command.append(str(inputs[0]))
    result = execute(command)
    return read_scores(targets, decoys), parse_selections(result.stderr)


def cv_attacks(binary: Path, root: Path) -> dict:
    output = {}
    for mode in ("fixed-c", "select-c", "ensemble"):
        row_sets = [cv_rows(engine=0), cv_rows(engine=1)] if mode == "ensemble" else [cv_rows()]
        spectrum_sizes = defaultdict(int)
        for rows in row_sets:
            for row in rows:
                spectrum_sizes[row[2]] += 1
        folds = spectrum_folds(dict(spectrum_sizes))
        fold0 = sorted(scan for scan, fold in folds.items() if fold == 0)
        sentinel_scans = set(fold0[:12])
        variant_runs = {}
        for variant in ("clean", "labels", "outliers", "reordered"):
            paths = []
            for engine, rows in enumerate(row_sets):
                path = root / f"cv-{mode}-{variant}-{engine}.pin"
                write_cv_variant(path, rows, folds, variant, sentinel_scans)
                paths.append(path)
            scores, selections = run_cv(binary, root, mode, variant, paths)
            variant_runs[variant] = {"scores": scores, "selections": selections}

        clean = variant_runs["clean"]
        fold0_ids = {
            psm_id
            for psm_id in clean["scores"]
            if int(psm_id.rsplit("CV", 1)[1].split("_", 1)[0]) in set(fold0)
        }
        sentinel_ids = {
            psm_id
            for psm_id in clean["scores"]
            if int(psm_id.rsplit("CV", 1)[1].split("_", 1)[0]) in sentinel_scans
        }
        label_changed = sorted(psm_id for psm_id in fold0_ids if variant_runs["labels"]["scores"].get(psm_id) != clean["scores"][psm_id])
        outlier_sentinel_changed = sorted(psm_id for psm_id in sentinel_ids if variant_runs["outliers"]["scores"].get(psm_id) != clean["scores"][psm_id])
        reorder_changed = sorted(psm_id for psm_id in clean["scores"] if variant_runs["reordered"]["scores"].get(psm_id) != clean["scores"][psm_id])
        clean_fold0_selection = [item for item in clean["selections"] if item.startswith("0:")]
        output[mode] = {
            "fold0_spectra": len(fold0),
            "heldout_rows_checked_for_label_attack": len(fold0_ids),
            "sentinel_rows_checked_for_outlier_attack": len(sentinel_ids),
            "heldout_label_changed_scores": label_changed,
            "outlier_changed_untouched_sentinel_scores": outlier_sentinel_changed,
            "reorder_changed_scores": reorder_changed,
            "fold0_selection_clean": clean_fold0_selection,
            "fold0_selection_after_labels": [item for item in variant_runs["labels"]["selections"] if item.startswith("0:")],
            "fold0_selection_after_outliers": [item for item in variant_runs["outliers"]["selections"] if item.startswith("0:")],
            "fold0_selection_after_reorder": [item for item in variant_runs["reordered"]["selections"] if item.startswith("0:")],
        }
    return output


def protein_representative_attack(binary: Path, root: Path) -> dict:
    rows = []
    # Same target peptide, same scan and identical features, but incompatible
    # protein mappings.  Peptide statistics are identical; protein inference must
    # not acquire an arbitrary mapping from whichever equal row appeared first.
    rows.extend([
        ("AMB_A", 1, 1, 500.0, 5.0, 0.0, 1.0, "K.AMBIGUOUS.R", "PROT_A"),
        ("AMB_B", 1, 1, 500.0, 5.0, 0.0, 1.0, "K.AMBIGUOUS.R", "PROT_B"),
        ("UNIQUE_A", 1, 2, 500.0, 4.0, 0.0, 1.0, "K.UNIQUEA.R", "PROT_A"),
        ("UNIQUE_C", 1, 3, 500.0, 3.0, 0.0, 1.0, "K.UNIQUEC.R", "PROT_C"),
    ])
    # Stabilize three-fold training and supply target/decoy background.
    for scan in range(10, 250):
        label = 1 if scan % 2 else -1
        rows.append((
            f"BG{scan}", label, scan, 500.0, float(label) + (scan % 7) / 10.0,
            float(scan % 5), float(scan % 11), f"K.BG{scan}.R",
            f"{'PROT' if label > 0 else 'DECOY_PROT'}_BG{scan}",
        ))

    arms = {"a_first": rows, "b_first": [rows[1], rows[0], *rows[2:]]}
    results = {}
    for arm, arranged in arms.items():
        pin = root / f"protein-representative-{arm}.pin"
        targets = root / f"protein-representative-{arm}.target.tsv"
        decoys = root / f"protein-representative-{arm}.decoy.tsv"
        proteins = root / f"protein-representative-{arm}.proteins.tsv"
        decoy_proteins = root / f"protein-representative-{arm}.decoy-proteins.tsv"
        peptides = root / f"protein-representative-{arm}.peptides.tsv"
        decoy_peptides = root / f"protein-representative-{arm}.decoy-peptides.tsv"
        write_pin(pin, arranged)
        execute([
            str(binary), "--canonical", "--no-select-c", "--no-psm-competition",
            "--seed", "1", "--num-threads", "1", "--results-psms", str(targets),
            "--decoy-results-psms", str(decoys), "--results-proteins", str(proteins),
            "--decoy-results-proteins", str(decoy_proteins),
            "--results-peptides", str(peptides), "--decoy-results-peptides", str(decoy_peptides),
            str(pin),
        ])
        protein_rows = read_rows(proteins) + read_rows(decoy_proteins)
        ambiguous_peptide = next(
            row for row in read_rows(peptides) if row["peptide"] == "K.AMBIGUOUS.R"
        )
        results[arm] = {
            "groups": sorted((row["ProteinGroupId"], row["proteinIds"], row["posterior_error_prob"]) for row in protein_rows),
            "ambiguous_psm_scores": {
                row["PSMId"]: row["score"] for row in read_rows(targets) if row["PSMId"].startswith("AMB_")
            },
            "ambiguous_peptide_representative": {
                key: ambiguous_peptide[key]
                for key in ("PSMId", "score", "q-value", "posterior_error_prob", "proteinIds")
            },
            "all_picked_pep_na": all(row["posterior_error_prob"] == "NA" for row in protein_rows),
        }
    return {
        "arms": results,
        "ambiguous_scores_tied": len(set(results["a_first"]["ambiguous_psm_scores"].values())) == 1,
        "groups_invariant": results["a_first"]["groups"] == results["b_first"]["groups"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"binary does not exist: {binary}")
    if args.json.exists():
        parser.error(f"refusing to overwrite {args.json}")

    with tempfile.TemporaryDirectory(prefix="percolator-post-repair-audit-") as temporary:
        root = Path(temporary)
        result = {
            "binary": str(binary),
            "tied_fixture": tied_fixture_attacks(binary, root),
            "joined_input_order": joined_order_attack(binary, root),
            "cross_validation": cv_attacks(binary, root),
            "protein_representative": protein_representative_attack(binary, root),
        }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print(f"variable ties: winners={result['tied_fixture']['winner_invariant']} q={result['tied_fixture']['q_invariant']} PEP={result['tied_fixture']['pep_invariant']}")
    print(
        f"joined input order: winners={result['joined_input_order']['winner_invariant']} "
        f"changed={result['joined_input_order']['changed_winners']} "
        f"boundary-q01={result['joined_input_order']['boundary_q01_invariant']}"
    )
    for mode, attack in result["cross_validation"].items():
        print(
            f"{mode}: label-changed={len(attack['heldout_label_changed_scores'])} "
            f"outlier-sentinel-changed={len(attack['outlier_changed_untouched_sentinel_scores'])} "
            f"reorder-changed={len(attack['reorder_changed_scores'])}"
        )
    print(f"protein equal-score mapping invariant={result['protein_representative']['groups_invariant']}")
    print(f"machine-readable result: {args.json}")


if __name__ == "__main__":
    main()
