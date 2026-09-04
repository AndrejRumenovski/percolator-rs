#!/usr/bin/env python3
"""Freeze and validate the deterministic pre-intervention baseline."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pep_rootcause_controls import read_pin  # noqa: E402


EXPECTED = {
    "targets": 19545,
    "entrapment_targets": 258,
    "entrapment_decoys": 133,
    "ratio": 258 / 133,
    "adjusted_fdp": 0.01687315883053162,
    "direct_known_false_fraction": 0.013200306983883347,
    "seed1_pooled_pep_error": 0.018505443259603106,
    "raw_xcorr_depth": 6117,
    "raw_xcorr_entrapment_targets": 134,
    "raw_xcorr_entrapment_decoys": 133,
}


def digest(path: Path, algorithm: str = "sha256") -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def command(*args: str) -> str:
    return subprocess.check_output(args, text=True, stderr=subprocess.STDOUT).strip()


def raw_xcorr(pin_root: Path) -> dict:
    rows = []
    for directory in sorted(pin_root.glob("comet-*")):
        if directory.name == "comet-out" or not (directory / "comet.pin").exists():
            continue
        score, label, pure, _ = read_pin(directory / "comet.pin", "Xcorr", 1.0)
        rows.extend(zip(score.tolist(), label.tolist(), pure.tolist()))
    rows.sort(key=lambda row: -row[0])
    ent_t = ent_d = 0
    for depth, (_, label, pure) in enumerate(rows, 1):
        if pure and label > 0:
            ent_t += 1
        elif pure and label < 0:
            ent_d += 1
        if ent_d == EXPECTED["raw_xcorr_entrapment_decoys"]:
            return {"depth": depth, "entrapment_targets": ent_t,
                    "entrapment_decoys": ent_d, "ratio": ent_t / ent_d}
    raise RuntimeError("raw XCorr list never reached the frozen entrapment-decoy depth")


def split_entrapment(combined: Path, output: Path) -> None:
    keep = False
    with combined.open() as source, output.open("w") as target:
        for line in source:
            if line.startswith(">"):
                keep = line[1:].startswith("ENT_")
            if keep:
                target.write(line)


def close(left: float, right: float, tolerance: float = 1e-15) -> bool:
    return abs(left - right) <= tolerance


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    summary_path = args.output_root / "summary.json"
    summary = json.loads(summary_path.read_text())
    q01 = next(row for row in summary["q_regions"] if row["threshold"] == 0.01)
    observed = {
        "targets": q01["targets"],
        "entrapment_targets": q01["entrapment_targets"],
        "entrapment_decoys": q01["entrapment_decoys"],
        "ratio": q01["ent_target_over_ent_decoy"],
        "adjusted_fdp": q01["adjusted_fdp"],
        "direct_known_false_fraction": q01["f1_fdp"],
        "seed1_pooled_pep_error": summary["calibration"]["summary"]["observed_minus_predicted_adjusted"],
    }
    for key in ("targets", "entrapment_targets", "entrapment_decoys"):
        if observed[key] != EXPECTED[key]:
            raise RuntimeError(f"baseline gate failed for {key}: {observed[key]} != {EXPECTED[key]}")
    for key in ("ratio", "adjusted_fdp", "direct_known_false_fraction", "seed1_pooled_pep_error"):
        if not close(observed[key], EXPECTED[key]):
            raise RuntimeError(f"baseline gate failed for {key}: {observed[key]} != {EXPECTED[key]}")

    raw = raw_xcorr(args.source_root)
    for key, expected_key in (("depth", "raw_xcorr_depth"),
                              ("entrapment_targets", "raw_xcorr_entrapment_targets"),
                              ("entrapment_decoys", "raw_xcorr_entrapment_decoys")):
        if raw[key] != EXPECTED[expected_key]:
            raise RuntimeError(f"raw baseline gate failed for {key}: {raw[key]} != {EXPECTED[expected_key]}")

    parameter_dir = args.output_root / "search_parameters"
    parameter_dir.mkdir(exist_ok=True)
    parameter_files = []
    for source in sorted(args.source_root.glob("*-comet.params.txt")):
        destination = parameter_dir / source.name
        shutil.copyfile(source, destination)
        parameter_files.append({"source": str(source), "copy": str(destination),
                                "sha256": digest(destination), "bytes": destination.stat().st_size})

    entrapment_fasta = args.output_root / "entrapment_targets.fasta"
    split_entrapment(args.source_root / "combined.fasta", entrapment_fasta)
    binary = args.repo / "target/release/percolator-rs"
    crux = args.source_root / "crux-4.0.Linux.x86_64/bin/crux"
    pins = sorted(path for path in args.source_root.glob("comet-*/comet.pin")
                  if path.parent.name != "comet-out")
    spectra = sorted(args.source_root.glob("*.mzML"))
    result_files = sorted(args.results_root.glob("seed-1/comet-*/*.tsv"))
    git_status = command("git", "-C", str(args.repo), "status", "--short")
    production_diff = command("git", "-C", str(args.repo), "diff", "--", "src", "Cargo.toml", "Cargo.lock")
    manifest = {
        "schema_version": 1,
        "baseline_gate_passed": True,
        "repository": {
            "path": str(args.repo.resolve()),
            "commit": command("git", "-C", str(args.repo), "rev-parse", "HEAD"),
            "status_at_freeze": git_status.splitlines(),
            "production_source_diff_empty": production_diff == "",
        },
        "percolator_rs": {
            "path": str(binary), "sha256": digest(binary), "bytes": binary.stat().st_size,
            "canonical_parameters": ["--canonical", "--no-select-c", "--seed", "1",
                                     "--num-threads", "1", "default maxiter=10"],
        },
        "search_engine": {
            "path": str(crux), "sha256": digest(crux),
            "version": command(str(crux), "version"),
            "decoy_generation": "Comet internal concatenated search (decoy_search=1): reverse each target peptide while retaining the enzyme-terminal residue",
            "parameters": parameter_files,
        },
        "databases": {
            "native": {"path": str(args.source_root / "native.fasta"),
                       "sha256": digest(args.source_root / "native.fasta")},
            "entrapment_only": {"path": str(entrapment_fasta),
                                "sha256": digest(entrapment_fasta),
                                "bytes": entrapment_fasta.stat().st_size},
            "combined_target": {"path": str(args.source_root / "combined.fasta"),
                                "sha256": digest(args.source_root / "combined.fasta")},
        },
        "spectra": [{"path": str(path), "sha1": digest(path, "sha1"), "bytes": path.stat().st_size}
                    for path in spectra],
        "pins": [{"path": str(path), "sha256": digest(path), "bytes": path.stat().st_size}
                 for path in pins],
        "result_files": [{"path": str(path), "sha256": digest(path), "bytes": path.stat().st_size}
                         for path in result_files],
        "random_seeds": {"percolator_primary": 1, "prior_canonical_reproducibility": [1, 2, 3, 4, 5],
                         "comet_decoys": "deterministic peptide reversal; no random seed"},
        "fresh_seed1": {"q_regions": summary["q_regions"],
                        "pep_calibration": summary["calibration"],
                        "raw_xcorr_matched_depth": raw},
        "expected": EXPECTED,
        "commands": {
            "rescore": "percolator-rs --canonical --no-select-c --seed 1 --num-threads 1 --results-psms TARGET --decoy-results-psms DECOY PIN",
            "summarize": "python3 validation/pep_rootcause_experiments.py summarize --root RESULTS --output summary.json",
        },
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "platform": platform.platform()},
    }
    output = args.output_root / "baseline_manifest.json"
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"baseline_gate_passed": True, "observed": observed,
                      "raw_xcorr": raw, "manifest": str(output)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
