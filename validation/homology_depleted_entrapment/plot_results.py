#!/usr/bin/env python3
"""Render compact diagnostic figures from the frozen machine-readable results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CONDITIONS = (
    "original", "homology_depleted", "size_control_130363",
    "size_control_155921", "size_control_196613",
)
LABELS = ("Original", "Homology depleted", "Control 1", "Control 2", "Control 3")
COLORS = ("#4C78A8", "#E45756", "#72B7B2", "#54A24B", "#B279A2")


def save(fig, root: Path, name: str) -> None:
    fig.tight_layout()
    fig.savefig(root / f"{name}.png", dpi=180, bbox_inches="tight")
    fig.savefig(root / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-root", type=Path, required=True)
    args = parser.parse_args()
    analysis = args.experiment_root / "analysis"
    figures = args.experiment_root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    raw = json.loads((analysis / "raw_xcorr.json").read_text())
    result = json.loads((analysis / "rescored_results.json").read_text())

    # Raw versus rescored internal-null ratios at the frozen matched depth/q region.
    raw_ratio = [raw["conditions"][name]["matched_entrapment_decoy_depths"]["133"]["ratio"]
                 for name in CONDITIONS]
    rescored = [result["endpoints"][name]["ent_target_over_ent_decoy"] for name in CONDITIONS]
    intervals = [result["uncertainty"]["conditions"][name]["percentile95"]
                 for name in CONDITIONS]
    rescored_yerr = np.array([
        [estimate - interval[0] for estimate, interval in zip(rescored, intervals)],
        [interval[1] - estimate for estimate, interval in zip(rescored, intervals)],
    ])
    x = np.arange(len(CONDITIONS))
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    width = 0.36
    ax.bar(x - width / 2, raw_ratio, width, label="Raw XCorr at D_ent=133", color="#9DADBC")
    ax.bar(x + width / 2, rescored, width, label="Percolator q<0.01", color=COLORS,
           yerr=rescored_yerr, capsize=3,
           error_kw={"elinewidth": 1, "ecolor": "black"})
    ax.axhline(1.0, color="black", lw=1, ls="--")
    ax.set_xticks(x, LABELS, rotation=20, ha="right")
    ax.set_ylabel("Entrapment target / entrapment decoy")
    ax.set_title("Internal-null exchangeability before and after rescoring")
    ax.legend(frameon=False)
    save(fig, figures, "raw-vs-rescored-exchangeability")

    # Threshold curve.
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    for condition, label, color in zip(CONDITIONS, LABELS, COLORS):
        rows = result["q_thresholds"][condition]
        ax.plot([row["threshold"] for row in rows],
                [row["ent_target_over_ent_decoy"] for row in rows],
                marker="o", ms=4, label=label, color=color)
    ax.axhline(1.0, color="black", lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("Reported q-value threshold")
    ax.set_ylabel("Entrapment target / entrapment decoy")
    ax.set_title("Rescored internal-null ratio across predefined thresholds")
    ax.legend(frameon=False, ncol=2)
    save(fig, figures, "q-threshold-exchangeability")

    # PEP calibration curves.
    fig, ax = plt.subplots(figsize=(6.8, 5.4))
    for condition, label, color in zip(CONDITIONS, LABELS, COLORS):
        bins = result["pep_calibration"][condition]["bins"]
        xvalue = [row["mean_predicted_pep"] for row in bins]
        yvalue = [row["observed_adjusted"] for row in bins]
        ax.plot(xvalue, yvalue, marker="o", ms=3.5, label=label, color=color)
    grid = np.geomspace(1e-4, 1, 100)
    ax.plot(grid, grid, color="black", lw=1, ls="--", label="Calibrated")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e-4, 1); ax.set_ylim(1e-4, 1)
    ax.set_xlabel("Mean predicted PEP")
    ax.set_ylabel("Entrapment-adjusted observed error")
    ax.set_title("PEP calibration by preregistered bin")
    ax.legend(frameon=False, ncol=2)
    save(fig, figures, "pep-calibration")

    # Training dose response: ratio and PEP error share maxiter x values.
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), sharex=True)
    for condition, label, color in zip(CONDITIONS, LABELS, COLORS):
        rows = result["training_dose_response"][condition]
        iterations = [row["maxiter"] for row in rows]
        axes[0].plot(iterations, [row["ent_target_over_ent_decoy"] for row in rows],
                     marker="o", label=label, color=color)
        axes[1].plot(iterations,
                     [row["pooled_pep_observed_minus_predicted_adjusted"] for row in rows],
                     marker="o", label=label, color=color)
    axes[0].axhline(1.0, color="black", lw=1, ls="--")
    axes[0].set_ylabel("R_ent / D_ent at q<0.01")
    axes[1].axhline(0.0, color="black", lw=1, ls="--")
    axes[1].set_ylabel("Pooled observed − predicted PEP")
    for ax in axes:
        ax.set_xlabel("Maximum training iterations")
        ax.set_xticks([0, 1, 2, 3, 10])
    axes[0].set_title("Training-induced null imbalance")
    axes[1].set_title("Training-induced PEP error")
    axes[1].legend(frameon=False, fontsize=8)
    save(fig, figures, "training-dose-response")

    # Bin-level factorization, one facet per condition.
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.2))
    axes = axes.ravel()
    for axis, condition, label, color in zip(axes, CONDITIONS, LABELS, COLORS):
        bins = result["pep_calibration"][condition]["bins"]
        points = [(row["ent_target_over_ent_decoy"], row["observed_over_predicted_adjusted"])
                  for row in bins if row.get("ent_target_over_ent_decoy") and
                  row.get("observed_over_predicted_adjusted")]
        axis.scatter([point[0] for point in points], [point[1] for point in points],
                     s=28, color=color)
        if points:
            low = min(min(point) for point in points); high = max(max(point) for point in points)
            axis.plot([low, high], [low, high], color="black", lw=1, ls="--")
        axis.set_xscale("log"); axis.set_yscale("log")
        axis.set_title(label)
        axis.set_xlabel("E_T / E_D"); axis.set_ylabel("Observed / predicted")
    axes[-1].axis("off")
    fig.suptitle("PEP optimism tracks local entrapment-null imbalance", y=1.01)
    save(fig, figures, "bin-level-shared-cause")


if __name__ == "__main__":
    main()
