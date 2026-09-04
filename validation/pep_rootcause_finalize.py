#!/usr/bin/env python3
"""Assemble the audited PEP investigation evidence into one JSON artifact."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np


def load(path):
    return json.loads(path.read_text())


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    root = args.evidence_root
    baseline = load(root / "baseline-summary.json")
    bins = baseline["calibration"]["bins"]
    observed_ratio = np.array([row["observed_over_predicted_adjusted"] for row in bins])
    internal_ratio = np.array([row["ent_target_over_ent_decoy"] for row in bins])
    relation = {
        "log_ratio_pearson": float(np.corrcoef(np.log(observed_ratio), np.log(internal_ratio))[0, 1]),
        "observed_over_predicted_divided_by_ent_target_over_ent_decoy": [
            float(value) for value in observed_ratio / internal_ratio
        ],
        "median_absolute_log_fold_residual": float(np.median(np.abs(np.log(observed_ratio / internal_ratio)))),
        "interpretation": "under local entrapment adjustment, observed/predicted = (entT/entD) * (local decoy count / predicted false count)",
    }
    scripts = sorted(Path("validation").glob("pep_rootcause_*"))
    result = {
        "date": "2026-08-30",
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "audited_binary": str(args.binary.resolve()),
        "audited_binary_sha256": sha256(args.binary),
        "production_statistical_methodology_modified": False,
        "investigation_scripts": [{"path": str(path), "sha256": sha256(path)} for path in scripts if path.is_file()],
        "canonical_semi_tryptic": baseline,
        "enzN_enzC_ablation": load(root / "enz-ablation-summary.json"),
        "fully_tryptic": load(root / "fully-tryptic-summary.json"),
        "feature_distributions": {
            "canonical": load(root / "baseline-features.json"),
            "enzN_enzC_ablation": load(root / "enz-ablation-features.json"),
            "fully_tryptic": load(root / "fully-tryptic-features.json"),
        },
        "training_iteration_dose_response": load(root / "dose-summary.json"),
        "raw_comet_controls": load(root / "raw-controls.json"),
        "standalone_qvality_identical_scores": load(root / "qvality-identical-scores.json"),
        "cpp_percolator_default": {
            "calibration": load(root / "cpp-summary.json"),
            "individual_overlap_with_rust": load(root / "cpp-rust-overlap.json"),
        },
        "near_homology": load(root / "homology.json"),
        "synthetic": {
            "ideal_grid": load(root / "synthetic-ideal.json"),
            "edge_case_matrix": load(root / "synthetic-matrix.json"),
            "assumption_violations": load(root / "synthetic-violations.json"),
            "training_mechanisms": load(root / "synthetic-mechanism.json"),
            "production_oracle_selftest": (root / "probe-selftest.txt").read_text(),
        },
        "quantitative_shared_cause_test": relation,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
