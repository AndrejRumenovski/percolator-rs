#!/usr/bin/env python3
"""Render the complete calibration-bin appendix from audited JSON results."""

import argparse
import json
from pathlib import Path


def f(value, digits=6):
    return "—" if value is None else f"{value:.{digits}f}"


def standard_table(calibration):
    lines = [
        "| PEP interval | target PSMs | mean predicted | ent T | ent D | observed f=1 | f=1 Wilson 95% | observed adjusted | predicted−observed (adjusted) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in calibration["bins"]:
        wilson = row.get("observed_f1_wilson95")
        interval = "—" if not wilson else f"[{f(wilson[0])}, {f(wilson[1])}]"
        lines.append(
            f"| {row['bin']} | {row['n']:,} | {f(row['mean_predicted_pep'])} | "
            f"{row['entrapment_targets']:,} | {row['entrapment_decoys']:,} | "
            f"{f(row['observed_f1'])} | {interval} | {f(row['observed_adjusted'])} | "
            f"{f(row['predicted_minus_observed_adjusted'])} |"
        )
    return "\n".join(lines)


def add_summary(lines, title, result):
    lines += [f"## {title}", "", f"Global entrapment fraction among usable decoys: `{result['calibration']['f_global']:.6f}`.", "", standard_table(result["calibration"]), ""]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.evidence_root
    baseline = json.loads((root / "baseline-summary.json").read_text())
    ablation = json.loads((root / "enz-ablation-summary.json").read_text())
    fully = json.loads((root / "fully-tryptic-summary.json").read_text())
    cpp = json.loads((root / "cpp-summary.json").read_text())
    dose = json.loads((root / "dose-summary.json").read_text())
    raw = json.loads((root / "raw-controls.json").read_text())

    lines = [
        "# Complete PEP calibration tables", "",
        "Generated 2026-08-30 by `pep_rootcause_render_tables.py`. Positive calibration error means observed minus predicted; the final column below is the opposite sign requested in the audit (`predicted−observed`). The Wilson interval is only for the directly observed entrapment-target fraction and treats PSMs as independent; it is descriptive, not a cluster-robust interval and does not include uncertainty in the entrapment adjustment.", "",
    ]
    add_summary(lines, "Canonical semi-tryptic percolator-rs, pooled", baseline)
    for item in baseline["by_seed"]:
        lines += [f"## Canonical semi-tryptic, seed {item['seed']}", "", standard_table(item["calibration"]), ""]
    for item in baseline["by_dataset"]:
        lines += [f"## Canonical semi-tryptic, dataset `{item['dataset']}`", "", standard_table(item["calibration"]), ""]
    add_summary(lines, "enzN/enzC ablation, pooled", ablation)
    add_summary(lines, "Fully-tryptic search, pooled", fully)
    add_summary(lines, "C++ Percolator 3.09 default I-spline, pooled", cpp)

    lines += ["## Training-iteration dose response: all bin tables", ""]
    for design in ("semi", "fully"):
        for item in dose[design]:
            lines += [f"### {design}-tryptic, maxiter {item['maxiter']}", "", standard_table(item["calibration"]), ""]

    lines += ["## Raw Comet score controls", ""]
    for score_name in ("Xcorr", "lnExpect"):
        lines += [f"### {score_name}", "",
                  "These are local-error estimates made by applying the same analysis estimator to a raw score; the raw scores are not themselves PEPs.", "",
                  "| PEP interval | target PSMs | mean predicted | ent T | ent D | observed f=1 | observed adjusted | predicted−observed (adjusted) |",
                  "|---|---:|---:|---:|---:|---:|---:|---:|"]
        for row in raw[score_name]["bins"]:
            lines.append(f"| {row['bin']} | {row['n']:,} | {f(row['mean_pep'])} | {row['n_ent_target']:,} | {row['n_ent_decoy']:,} | {f(row['obs_f1'])} | {f(row['obs_adj'])} | {f(-row['gap_adj'])} |")
        lines.append("")
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
