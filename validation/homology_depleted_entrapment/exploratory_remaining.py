#!/usr/bin/env python3
"""Exploratory characterization of residual q<0.01 entrapment nulls."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pep_rootcause_experiments import classify, core_peptide, protein_fields  # noqa: E402
from pep_rootcause_homology import canonical, search  # noqa: E402


FAMILY_TERMS = (
    "actin", "tubulin", "histone", "ribosomal", "heat shock", "chaperon",
    "ATP synthase", "elongation factor", "ubiquitin", "myosin", "keratin",
)


def fasta_headers(path: Path) -> dict[str, str]:
    result = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                text = line[1:].rstrip()
                result[text.split()[0]] = text
    return result


def pin_rows(path: Path) -> dict[str, dict]:
    result = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["SpecId"]] = row
    return result


def result_rows(path: Path, decoy: bool, pins: dict[str, dict], headers: dict[str, str]) -> list[dict]:
    result = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pure, mixed = classify(protein_fields(row), decoy)
            if not pure or mixed or float(row["q-value"]) >= 0.01:
                continue
            feature = pins.get(row["PSMId"])
            if feature is None:
                raise KeyError(f"missing PIN row {row['PSMId']}")
            proteins = []
            for value in protein_fields(row):
                proteins.extend(item for item in value.replace(";", "\t").split("\t") if item)
            if decoy:
                proteins = [value.removeprefix("DECOY_") for value in proteins]
            family = []
            full_headers = [headers.get(value, value) for value in proteins]
            combined = " ".join(full_headers).lower()
            for term in FAMILY_TERMS:
                if term.lower() in combined:
                    family.append(term)
            result.append({
                "decoy": decoy,
                "psmid": row["PSMId"],
                "peptide": canonical(core_peptide(row["peptide"])),
                "score": float(row["score"]),
                "q": float(row["q-value"]),
                "pep": float(row["posterior_error_prob"]),
                "proteins": proteins,
                "protein_headers": full_headers,
                "families": family,
                "features": feature,
            })
    return result


def numeric_summary(rows: list[dict], field: str, nested: bool = False) -> dict:
    values = []
    for row in rows:
        raw = row["features"].get(field) if nested else row.get(field)
        if raw not in (None, "", "NA"):
            values.append(float(raw))
    array = np.asarray(values)
    return {"n": int(array.size), "mean": float(array.mean()) if array.size else None,
            "median": float(np.median(array)) if array.size else None,
            "q10_q90": [float(value) for value in np.quantile(array, [0.1, 0.9])]
            if array.size else None}


def summarize(rows: list[dict], similarity: dict[str, dict | None]) -> dict:
    peptides = sorted(set(row["peptide"] for row in rows))
    distances = [similarity[peptide]["distance"] for peptide in peptides
                 if similarity.get(peptide) is not None]
    feature_names = (
        "ExpMass", "CalcMass", "lnrSp", "deltLCn", "deltCn", "lnExpect",
        "Xcorr", "Sp", "IonFrac", "PepLen", "enzN", "enzC", "enzInt",
        "lnNumSP", "dM", "absdM",
    )
    return {
        "rows": len(rows), "distinct_peptides": len(peptides),
        "peptide_length": {"mean": float(np.mean([len(value) for value in peptides])) if peptides else None,
                           "counts": dict(sorted(Counter(len(value) for value in peptides).items()))},
        "score": numeric_summary(rows, "score"),
        "pep": numeric_summary(rows, "pep"),
        "features": {name: numeric_summary(rows, name, True) for name in feature_names},
        "charge_counts": {str(charge): sum(row["features"].get(f"Charge{charge}") == "1" for row in rows)
                          for charge in range(1, 7)},
        "family_row_counts": {term: sum(term in row["families"] for row in rows) for term in FAMILY_TERMS},
        "proteins": {"distinct": len(set(value for row in rows for value in row["proteins"])),
                     "top": Counter(value for row in rows for value in row["proteins"]).most_common(25)},
        "native_substring_distance_le_2": {
            "distinct_hits": len(distances),
            "distance_0_1_2": [distances.count(value) for value in range(3)],
            "matches": [{"peptide": peptide, **similarity[peptide]} for peptide in peptides
                        if similarity.get(peptide) is not None],
        },
        "top_psms": sorted([{key: row[key] for key in ("psmid", "peptide", "score", "pep", "proteins", "families")}
                            for row in rows], key=lambda row: (-row["score"], row["psmid"]))[:50],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-root", type=Path, default=HERE)
    parser.add_argument("--native-fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    condition = "homology_depleted"
    headers = fasta_headers(args.experiment_root / "databases/homology_depleted.fasta")
    rows = []
    for result_dir in sorted((args.experiment_root / f"percolator/primary/{condition}/seed-1").glob("comet-*")):
        run = result_dir.name
        pins = pin_rows(args.experiment_root / f"searches/{condition}/{run}/comet.pin")
        rows.extend(result_rows(result_dir / "target.tsv", False, pins, headers))
        rows.extend(result_rows(result_dir / "decoy.tsv", True, pins, headers))
    peptides = sorted(set(row["peptide"] for row in rows if row["peptide"]))
    found = search(args.native_fasta, peptides)
    similarity = dict(zip(peptides, found))
    targets = [row for row in rows if not row["decoy"]]
    decoys = [row for row in rows if row["decoy"]]
    result = {
        "status": "exploratory; no second filter or confirmatory threshold is constructed",
        "region": "q<0.01, seed 1, canonical maxiter 10, pure entrapment only",
        "entrapment_targets": summarize(targets, similarity),
        "entrapment_decoys": summarize(decoys, similarity),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"entrapment_targets": len(targets), "entrapment_decoys": len(decoys),
                      "target_distinct_peptides": result["entrapment_targets"]["distinct_peptides"],
                      "decoy_distinct_peptides": result["entrapment_decoys"]["distinct_peptides"]},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
