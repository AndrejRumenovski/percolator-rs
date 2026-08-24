#!/usr/bin/env python3
"""Select Bayesian alpha/beta/gamma on replicate 1 of the PrEST standard."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import subprocess
import tempfile
from pathlib import Path

import report


ALPHAS = (0.01, 0.05, 0.1, 0.3, 0.5)
BETAS = (0.0001, 0.001, 0.01, 0.05)
GAMMAS = (0.001, 0.01, 0.1, 0.5)
MAX_ITER = 1000


def mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / len(finite) if finite else math.nan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    truth = report.read_truth(args.truth)
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        calibration = [
            row for row in csv.DictReader(handle, delimiter="\t") if row["split"] == "calibration"
        ]
    if {row["vial"] for row in calibration} != {"A", "B", "AB", "BLANK"}:
        raise ValueError("calibration split must contain A, B, AB, and BLANK")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    combinations = list(itertools.product(ALPHAS, BETAS, GAMMAS))
    scored: list[dict[str, float]] = []
    with tempfile.TemporaryDirectory(prefix="protein-grid-", dir=args.output_dir) as temporary:
        temp = Path(temporary)
        for number, (alpha, beta, gamma) in enumerate(combinations, start=1):
            briers: list[float] = []
            eces: list[float] = []
            aucs: list[float] = []
            paucs: list[float] = []
            converged = True
            for item in calibration:
                target = temp / f"{item['sample']}.target.tsv"
                decoy = temp / f"{item['sample']}.decoy.tsv"
                command = [
                    str(args.binary), "--canonical", "--seed", "1",
                    "--protein-inference", "bayesian",
                    "--protein-alpha", str(alpha),
                    "--protein-beta", str(beta),
                    "--protein-gamma", str(gamma),
                    "--protein-max-iter", str(MAX_ITER),
                    "--results-proteins", str(target),
                    "--decoy-results-proteins", str(decoy),
                    item["pin"],
                ]
                completed = subprocess.run(
                    command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                    text=True, check=False,
                )
                if completed.returncode:
                    raise RuntimeError(
                        f"parameter run failed for {item['sample']} "
                        f"alpha={alpha} beta={beta} gamma={gamma}:\n{completed.stderr}"
                    )
                if "converged: true" not in completed.stderr:
                    converged = False
                groups = report.load_groups(target, item["vial"], truth)
                brier, ece = report.probability_metrics(groups)
                briers.append(brier)
                eces.append(ece)
                aucs.append(report.auc(groups))
                paucs.append(report.partial_auc(groups))

            average_brier = mean(briers)
            average_ece = mean(eces)
            average_auc = mean(aucs)
            average_pauc = mean(paucs)
            # Ground-truth probability calibration is primary. A smaller ranking
            # term prevents a flat, well-calibrated prior from winning the grid.
            objective = average_brier + average_ece + 0.25 * (1.0 - average_auc)
            if not converged:
                objective = math.inf
            scored.append(
                {
                    "alpha": alpha, "beta": beta, "gamma": gamma,
                    "mean_brier": average_brier, "mean_ece": average_ece,
                    "mean_auc": average_auc, "mean_pauc": average_pauc,
                    "objective": objective, "converged": converged,
                }
            )
            if number % 10 == 0 or number == len(combinations):
                print(f"parameter grid: {number}/{len(combinations)}", flush=True)

    scored.sort(
        key=lambda row: (
            row["objective"], row["mean_brier"], row["mean_ece"],
            -row["mean_auc"], row["alpha"], row["beta"], row["gamma"],
        )
    )
    if not math.isfinite(scored[0]["objective"]):
        raise RuntimeError(f"no parameter combination converged within {MAX_ITER} iterations")
    grid_path = args.output_dir / "calibration-grid.tsv"
    fields = (
        "alpha", "beta", "gamma", "mean_brier", "mean_ece",
        "mean_auc", "mean_pauc", "objective", "converged",
    )
    with grid_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in scored:
            writer.writerow(
                {
                    key: (str(row[key]).lower() if key == "converged" else f"{row[key]:.10g}")
                    for key in fields
                }
            )

    best = scored[0]
    selection = {
        "alpha": best["alpha"],
        "beta": best["beta"],
        "gamma": best["gamma"],
        "peptide_prior": 0.1,
        "max_iter": MAX_ITER,
        "objective": best["objective"],
        "selection_split": "replicate 1 only",
        "objective_definition": "mean Brier + mean 10-bin ECE + 0.25*(1-mean ROC AUC)",
        "grid": {"alpha": ALPHAS, "beta": BETAS, "gamma": GAMMAS},
    }
    (args.output_dir / "selected-params.json").write_text(
        json.dumps(selection, indent=2) + "\n", encoding="utf-8"
    )
    (args.output_dir / "selected-params.env").write_text(
        f"PROTEIN_ALPHA={best['alpha']}\n"
        f"PROTEIN_BETA={best['beta']}\n"
        f"PROTEIN_GAMMA={best['gamma']}\n"
        f"PROTEIN_MAX_ITER={MAX_ITER}\n",
        encoding="utf-8",
    )
    print(
        f"selected alpha={best['alpha']} beta={best['beta']} gamma={best['gamma']} "
        f"objective={best['objective']:.8g}"
    )
    print(f"wrote {grid_path}")


if __name__ == "__main__":
    main()
