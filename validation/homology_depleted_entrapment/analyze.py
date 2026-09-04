#!/usr/bin/env python3
"""Analyze the preregistered homology-depletion experiment without corrections."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pep_rootcause_controls import read_pin  # noqa: E402
from pep_rootcause_experiments import (  # noqa: E402
    Q_REGIONS,
    calibration,
    cumulative_regions,
    load_outputs,
    summarize,
)
from pep_rootcause_lib import qvalues_and_peps  # noqa: E402


CONDITIONS = (
    "original", "homology_depleted", "size_control_130363",
    "size_control_155921", "size_control_196613",
)
CONTROL_CONDITIONS = CONDITIONS[2:]
MATCHED_DEPTHS = (25, 50, 100, 133, 250, 500)
PRIMARY_Q = 0.01


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value) -> bool:
    return value is not None and math.isfinite(value)


def internal_curve(score: np.ndarray, label: np.ndarray, pure: np.ndarray) -> dict:
    order = np.argsort(-score, kind="stable")
    ent_t = np.cumsum((label[order] > 0) & pure[order])
    ent_d = np.cumsum((label[order] < 0) & pure[order])
    result = {}
    for depth in MATCHED_DEPTHS:
        position = int(np.searchsorted(ent_d, depth))
        if position >= order.size:
            result[str(depth)] = None
        else:
            result[str(depth)] = {
                "rank_depth": position + 1,
                "entrapment_targets": int(ent_t[position]),
                "entrapment_decoys": depth,
                "ratio": float(ent_t[position] / depth),
            }
    return result


def raw_condition(search_root: Path, condition: str) -> dict:
    score_parts, label_parts, pure_parts, mixed_parts = [], [], [], []
    pin_manifest = []
    per_run = []
    for directory in sorted((search_root / condition).glob("comet-*")):
        pin = directory / "comet.pin"
        if not pin.exists():
            continue
        score, label, pure, mixed = read_pin(pin, "Xcorr", 1.0)
        qvalue, pep = qvalues_and_peps(score, label)
        score_parts.append(score); label_parts.append(label)
        pure_parts.append(pure); mixed_parts.append(mixed)
        per_run.append({"run": directory.name, "rows_after_spectrum_competition": int(score.size),
                        "targets": int((label > 0).sum()), "decoys": int((label < 0).sum()),
                        "pure_entrapment_targets": int(((label > 0) & pure).sum()),
                        "pure_entrapment_decoys": int(((label < 0) & pure).sum())})
        pin_manifest.append({"path": str(pin), "sha256": sha256(pin), "bytes": pin.stat().st_size})
    if len(pin_manifest) != 6:
        raise ValueError(f"{condition}: expected six PINs, found {len(pin_manifest)}")
    score = np.concatenate(score_parts); label = np.concatenate(label_parts)
    pure = np.concatenate(pure_parts); mixed = np.concatenate(mixed_parts)
    # Recompute per-run q-values above, then concatenate them in the same order.
    q_parts = []
    for directory in sorted((search_root / condition).glob("comet-*")):
        pin = directory / "comet.pin"
        if pin.exists():
            s, l, _, _ = read_pin(pin, "Xcorr", 1.0)
            q, _ = qvalues_and_peps(s, l)
            q_parts.append(q)
    qvalue = np.concatenate(q_parts)
    d = label < 0
    usable = d & ~mixed
    f_global = float(pure[usable].sum()) / int(usable.sum())
    thresholds = []
    for threshold in Q_REGIONS:
        selected = qvalue < threshold
        targets = selected & (label > 0)
        decoys = selected & (label < 0)
        et = int((targets & pure).sum()); ed = int((decoys & pure).sum())
        n = int(targets.sum())
        thresholds.append({
            "threshold": threshold, "targets": n,
            "entrapment_targets": et, "entrapment_decoys": ed,
            "ent_target_over_ent_decoy": et / ed if ed else None,
            "direct_known_false_fraction": et / n if n else None,
            "adjusted_fdp_global_f": et / f_global / n if n else None,
        })
    return {
        "condition": condition, "rows": int(score.size), "f_global": f_global,
        "matched_entrapment_decoy_depths": internal_curve(score, label, pure),
        "q_thresholds": thresholds, "per_run": per_run, "pins": pin_manifest,
    }


def analyze_raw(experiment_root: Path) -> dict:
    result = {
        "schema_version": 1,
        "stage": "raw Comet XCorr before Percolator interpretation",
        "competition_key": "ScanNr + ExpMass, highest Xcorr, per PIN",
        "matched_depths": list(MATCHED_DEPTHS),
        "conditions": {condition: raw_condition(experiment_root / "searches", condition)
                       for condition in CONDITIONS},
    }
    output = experiment_root / "analysis/raw_xcorr.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def qrow(summary: dict, threshold: float = PRIMARY_Q) -> dict:
    return next(row for row in summary["q_regions"] if row["threshold"] == threshold)


def compact_endpoint(summary: dict) -> dict:
    row = qrow(summary)
    return {
        "accepted_targets": row["targets"],
        "entrapment_targets": row["entrapment_targets"],
        "entrapment_decoys": row["entrapment_decoys"],
        "ent_target_over_ent_decoy": row["ent_target_over_ent_decoy"],
        "direct_known_false_fraction": row["f1_fdp"],
        "adjusted_fdp": row["adjusted_fdp"],
        "mean_predicted_pep": row["mean_predicted_pep"],
        "pooled_pep_observed_minus_predicted_adjusted": summary["calibration"]["summary"]["observed_minus_predicted_adjusted"],
        "pooled_pep_observed_minus_predicted_direct": summary["calibration"]["summary"]["observed_minus_predicted_f1"],
    }


def bin_relationship(summary: dict) -> dict:
    used = []
    for row in summary["calibration"]["bins"]:
        observed_ratio = row.get("observed_over_predicted_adjusted")
        internal_ratio = row.get("ent_target_over_ent_decoy")
        if finite(observed_ratio) and finite(internal_ratio) and observed_ratio > 0 and internal_ratio > 0:
            used.append((row["bin"], observed_ratio, internal_ratio))
    if len(used) < 2:
        correlation = None
        residual = None
    else:
        observed = np.log([row[1] for row in used])
        internal = np.log([row[2] for row in used])
        correlation = float(np.corrcoef(observed, internal)[0, 1])
        residual = float(np.median(np.abs(observed - internal)))
    return {
        "bins_used": [row[0] for row in used], "n_bins": len(used),
        "log_ratio_pearson": correlation,
        "median_absolute_log_fold_residual": residual,
        "median_fold_residual": math.exp(residual) if residual is not None else None,
    }


def dataset_counts(data: dict) -> dict[str, dict[str, int]]:
    result = {}
    selected = data["q"] < PRIMARY_Q
    for name in sorted(set(str(value) for value in data["dataset_name"])):
        dataset = data["dataset_name"] == name
        mask = selected & dataset
        target = dataset & ~data["decoy"]
        usable_decoy = dataset & data["decoy"] & ~data["mixed"]
        result[name] = {
            "targets": int((mask & ~data["decoy"]).sum()),
            "entrapment_targets": int((mask & ~data["decoy"] & data["pure"]).sum()),
            "entrapment_decoys": int((mask & data["decoy"] & data["pure"]).sum()),
            "calibration_targets": int(target.sum()),
            "calibration_target_pep_sum": float(data["pep"][target].sum()),
            "calibration_entrapment_targets": int((target & data["pure"]).sum()),
            "usable_decoys": int(usable_decoy.sum()),
            "pure_entrapment_decoys": int((usable_decoy & data["pure"]).sum()),
        }
    return result


def run_bootstrap(per_condition: dict[str, dict[str, dict[str, int]]], replicates: int = 10000) -> dict:
    runs = sorted(per_condition["original"])
    rng = np.random.default_rng(20260831)
    metrics = ("r_ent_over_d_ent", "pooled_pep_error", "q_lt_0_01_adjusted_fdp")
    values = {metric: {condition: [] for condition in CONDITIONS} for metric in metrics}
    contrasts = {
        "log_homology_minus_mean_log_controls": [],
        "homology_minus_original_pep_error": [],
        "homology_minus_mean_control_pep_error": [],
        "homology_minus_original_adjusted_fdp": [],
        "homology_minus_mean_control_adjusted_fdp": [],
    }
    for control in CONTROL_CONDITIONS:
        contrasts[f"log_homology_minus_{control}"] = []
    for _ in range(replicates):
        draw = rng.integers(0, len(runs), size=len(runs))
        estimates = {metric: {} for metric in metrics}
        for condition in CONDITIONS:
            def total(field: str) -> float:
                return sum(per_condition[condition][runs[index]][field] for index in draw)

            et = total("entrapment_targets")
            ed = total("entrapment_decoys")
            usable_decoys = total("usable_decoys")
            pure_decoys = total("pure_entrapment_decoys")
            f = pure_decoys / usable_decoys if usable_decoys else math.nan
            calibration_targets = total("calibration_targets")
            q_targets = total("targets")
            estimates["r_ent_over_d_ent"][condition] = et / ed if ed else math.nan
            estimates["pooled_pep_error"][condition] = (
                (total("calibration_entrapment_targets") / f -
                 total("calibration_target_pep_sum")) / calibration_targets
                if f and calibration_targets else math.nan
            )
            estimates["q_lt_0_01_adjusted_fdp"][condition] = (
                et / f / q_targets if f and q_targets else math.nan
            )
            for metric in metrics:
                values[metric][condition].append(estimates[metric][condition])
        ratios = estimates["r_ent_over_d_ent"]
        if all(math.isfinite(ratios[condition]) and ratios[condition] > 0 for condition in CONDITIONS):
            contrasts["log_homology_minus_mean_log_controls"].append(
                math.log(ratios["homology_depleted"]) -
                float(np.mean([math.log(ratios[name]) for name in CONTROL_CONDITIONS]))
            )
            for control in CONTROL_CONDITIONS:
                contrasts[f"log_homology_minus_{control}"].append(
                    math.log(ratios["homology_depleted"]) - math.log(ratios[control])
                )
        for metric, original_key, control_key in (
            ("pooled_pep_error", "homology_minus_original_pep_error",
             "homology_minus_mean_control_pep_error"),
            ("q_lt_0_01_adjusted_fdp", "homology_minus_original_adjusted_fdp",
             "homology_minus_mean_control_adjusted_fdp"),
        ):
            estimate = estimates[metric]
            if all(math.isfinite(estimate[condition]) for condition in CONDITIONS):
                contrasts[original_key].append(
                    estimate["homology_depleted"] - estimate["original"]
                )
                contrasts[control_key].append(
                    estimate["homology_depleted"] -
                    float(np.mean([estimate[name] for name in CONTROL_CONDITIONS]))
                )

    def interval(raw: list[float]) -> dict:
        valid = np.asarray([value for value in raw if math.isfinite(value)])
        return {
            "valid_replicates": int(valid.size),
            "percentile95": [float(np.percentile(valid, 2.5)), float(np.percentile(valid, 97.5))]
            if valid.size else None,
            "fraction_below_zero": float((valid < 0).mean()) if valid.size else None,
        }

    result = {"unit": "LC-MS/MS run", "n_units": len(runs), "replicates": replicates,
              "rng_seed": 20260831, "conditions": {}}
    for condition in CONDITIONS:
        result["conditions"][condition] = interval(values["r_ent_over_d_ent"][condition])
        result["conditions"][condition].pop("fraction_below_zero")
    result["log_homology_minus_mean_log_controls"] = interval(
        contrasts["log_homology_minus_mean_log_controls"]
    )
    result["metrics"] = {
        metric: {
            "conditions": {condition: interval(values[metric][condition])
                           for condition in CONDITIONS}
        }
        for metric in metrics
    }
    result["contrasts"] = {name: interval(raw) for name, raw in contrasts.items()}
    return result


def load_slim(root: Path) -> tuple[dict, dict]:
    data = load_outputs(root, "target.tsv", None)
    mask = np.ones(data["pep"].size, dtype=bool)
    slim = {
        "calibration": calibration(data, mask),
        "q_regions": cumulative_regions(data, mask, "q", Q_REGIONS),
    }
    return slim, data


def dose_analysis(experiment_root: Path, primary_summaries: dict) -> dict:
    result = {}
    for condition in CONDITIONS:
        rows = []
        for maxiter in (0, 1, 2, 3):
            summary, _ = load_slim(experiment_root / f"percolator/dose/mi-{maxiter}/{condition}")
            rows.append({"maxiter": maxiter, **compact_endpoint(summary)})
        rows.append({"maxiter": 10, **compact_endpoint(primary_summaries[condition])})
        result[condition] = rows
    return result


def enz_analysis(experiment_root: Path, primary_summaries: dict) -> dict:
    result = {}
    for condition in ("original", "homology_depleted"):
        ablated, _ = load_slim(experiment_root / f"percolator/enz_ablation/{condition}")
        result[condition] = {
            "all_features": compact_endpoint(primary_summaries[condition]),
            "without_enzN_enzC": compact_endpoint(ablated),
        }
    return result


def seed_analysis(experiment_root: Path, primary_summaries: dict) -> dict:
    result = {}
    for condition in ("original", "homology_depleted"):
        rows = [{"seed": 1, **compact_endpoint(primary_summaries[condition])}]
        data = load_outputs(experiment_root / f"percolator/seeds/{condition}", "target.tsv", None)
        for seed in (2, 3):
            mask = data["seed"] == seed
            summary = {
                "calibration": calibration(data, mask),
                "q_regions": cumulative_regions(data, mask, "q", Q_REGIONS),
            }
            rows.append({"seed": seed, **compact_endpoint(summary)})
        result[condition] = rows
    return result


def analyze_rescored(experiment_root: Path) -> dict:
    summaries = {}
    data = {}
    for condition in CONDITIONS:
        root = experiment_root / f"percolator/primary/{condition}"
        loaded = load_outputs(root, "target.tsv", None)
        data[condition] = loaded
        summaries[condition] = summarize(loaded)
        output = experiment_root / f"analysis/summary_{condition}.json"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(summaries[condition], indent=2, sort_keys=True) + "\n")

    endpoints = {condition: compact_endpoint(summary) for condition, summary in summaries.items()}
    original_ratio = endpoints["original"]["ent_target_over_ent_decoy"]
    homology_ratio = endpoints["homology_depleted"]["ent_target_over_ent_decoy"]
    control_ratios = [endpoints[name]["ent_target_over_ent_decoy"] for name in CONTROL_CONDITIONS]
    comparison = {
        "original_ratio": original_ratio,
        "homology_ratio": homology_ratio,
        "control_ratios": control_ratios,
        "mean_control_ratio": float(np.mean(control_ratios)),
        "homology_minus_original": homology_ratio - original_ratio,
        "mean_control_minus_original": float(np.mean(control_ratios)) - original_ratio,
        "log_homology_minus_mean_log_controls": (
            math.log(homology_ratio) - float(np.mean(np.log(control_ratios)))
            if homology_ratio > 0 and all(value > 0 for value in control_ratios) else None
        ),
        "fraction_original_excess_removed_total": (
            (original_ratio - homology_ratio) / (original_ratio - 1.0)
            if original_ratio != 1.0 else None
        ),
        "fraction_original_excess_attributable_beyond_size_controls": (
            (float(np.mean(control_ratios)) - homology_ratio) / (original_ratio - 1.0)
            if original_ratio != 1.0 else None
        ),
        "fraction_original_excess_attributable_to_size_control_change": (
            (original_ratio - float(np.mean(control_ratios))) / (original_ratio - 1.0)
            if original_ratio != 1.0 else None
        ),
    }
    per_condition_runs = {condition: dataset_counts(loaded) for condition, loaded in data.items()}
    pep_errors = {condition: endpoints[condition]["pooled_pep_observed_minus_predicted_adjusted"]
                  for condition in CONDITIONS}
    fdp_values = {condition: endpoints[condition]["adjusted_fdp"] for condition in CONDITIONS}
    calibration_comparison = {
        "pooled_pep_error": pep_errors,
        "mean_control_pep_error": float(np.mean([pep_errors[name] for name in CONTROL_CONDITIONS])),
        "homology_minus_original_pep_error": pep_errors["homology_depleted"] - pep_errors["original"],
        "homology_minus_mean_control_pep_error": (
            pep_errors["homology_depleted"] - float(np.mean([pep_errors[name] for name in CONTROL_CONDITIONS]))
        ),
        "q_lt_0_01_adjusted_fdp": fdp_values,
        "mean_control_adjusted_fdp": float(np.mean([fdp_values[name] for name in CONTROL_CONDITIONS])),
        "homology_minus_original_adjusted_fdp": fdp_values["homology_depleted"] - fdp_values["original"],
        "homology_minus_mean_control_adjusted_fdp": (
            fdp_values["homology_depleted"] - float(np.mean([fdp_values[name] for name in CONTROL_CONDITIONS]))
        ),
    }
    result = {
        "schema_version": 1,
        "primary_endpoint": "R_ent / D_ent at q<0.01, seed 1, canonical maxiter 10",
        "endpoints": endpoints,
        "primary_comparison": comparison,
        "calibration_comparison": calibration_comparison,
        "q_thresholds": {condition: summary["q_regions"] for condition, summary in summaries.items()},
        "pep_calibration": {condition: summary["calibration"] for condition, summary in summaries.items()},
        "bin_relationship": {condition: bin_relationship(summary) for condition, summary in summaries.items()},
        "per_run_primary_counts": per_condition_runs,
        "uncertainty": run_bootstrap(per_condition_runs),
        "training_dose_response": dose_analysis(experiment_root, summaries),
        "enzN_enzC_interaction": enz_analysis(experiment_root, summaries),
        "seed_reproducibility": seed_analysis(experiment_root, summaries),
    }
    output = experiment_root / "analysis/rescored_results.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("raw", "rescored"))
    parser.add_argument("--experiment-root", type=Path, default=HERE)
    args = parser.parse_args()
    if args.stage == "raw":
        result = analyze_raw(args.experiment_root)
        print(json.dumps({condition: row["matched_entrapment_decoy_depths"]["133"]
                          for condition, row in result["conditions"].items()}, indent=2, sort_keys=True))
    else:
        result = analyze_rescored(args.experiment_root)
        print(json.dumps({"endpoints": result["endpoints"],
                          "primary_comparison": result["primary_comparison"]},
                         indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
