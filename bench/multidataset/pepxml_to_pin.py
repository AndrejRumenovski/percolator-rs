#!/usr/bin/env python3
"""Convert unrescored pepXML search hits to a deterministic PIN table.

This deliberately has no third-party dependencies.  It retains every search hit,
uses all pepXML search scores as numeric features, and labels a hit as target when
at least one mapped protein does not have the requested decoy prefix.
"""

from __future__ import annotations

import argparse
import math
import xml.etree.ElementTree as ET
from pathlib import Path


PROTON = 1.00727646677


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def modified_peptide(hit: ET.Element) -> str:
    peptide = hit.attrib["peptide"]
    mod_info = next(
        (child for child in hit if local_name(child.tag) == "modification_info"),
        None,
    )
    if mod_info is None:
        return peptide

    insertions: list[tuple[int, str]] = []
    nterm = mod_info.get("mod_nterm_mass")
    if nterm is not None:
        insertions.append((0, f"[{nterm}]"))
    cterm = mod_info.get("mod_cterm_mass")
    if cterm is not None:
        insertions.append((len(peptide), f"[{cterm}]"))
    for mod in mod_info:
        if local_name(mod.tag) == "mod_aminoacid_mass":
            insertions.append((int(mod.attrib["position"]), f"[{mod.attrib['mass']}]"))

    for position, text in sorted(insertions, reverse=True):
        peptide = peptide[:position] + text + peptide[position:]
    return peptide


def scan_score_names(path: Path) -> list[str]:
    names: set[str] = set()
    for _, elem in ET.iterparse(path, events=("end",)):
        if local_name(elem.tag) == "search_score":
            names.add(elem.attrib["name"])
        elem.clear()
    return sorted(names)


def rows(path: Path, decoy_prefix: str, score_names: list[str]):
    base_name = path.stem
    scan = charge = exp_mass = retention_time = None
    for event, elem in ET.iterparse(path, events=("start", "end")):
        name = local_name(elem.tag)
        if event == "start" and name == "msms_run_summary":
            base_name = elem.get("base_name", base_name)
        elif event == "start" and name == "spectrum_query":
            scan = int(elem.attrib["end_scan"])
            charge = int(elem.attrib["assumed_charge"])
            exp_mass = float(elem.attrib["precursor_neutral_mass"])
            retention_time = float(elem.get("retention_time_sec", "0"))
        elif event == "end" and name == "search_hit":
            assert scan is not None and charge is not None and exp_mass is not None
            proteins = [elem.attrib["protein"].split(" ", 1)[0]]
            proteins.extend(
                child.attrib["protein"].split(" ", 1)[0]
                for child in elem
                if local_name(child.tag) == "alternative_protein"
            )
            label = 1 if any(not p.startswith(decoy_prefix) for p in proteins) else -1
            calc_mass = float(elem.attrib["calc_neutral_pep_mass"])
            exp_mz = exp_mass / charge + PROTON
            calc_mz = calc_mass / charge + PROTON
            scores = {
                child.attrib["name"]: float(child.attrib["value"])
                for child in elem
                if local_name(child.tag) == "search_score"
            }
            # Match mokapot's useful treatment of expectation values while making
            # the favorable direction explicit: larger NegLog10Expect is better.
            expect = scores.get("expect", 1.0)
            neg_log_expect = -math.log10(max(expect, 1e-300))
            spec_id = f"{base_name}_{scan}_{charge}_{elem.get('hit_rank', '1')}"
            fields = [
                spec_id,
                str(label),
                str(scan),
                f"{exp_mass:.10g}",
                f"{calc_mass:.10g}",
                f"{retention_time or 0:.10g}",
                elem.get("hit_rank", "1"),
                elem.get("num_missed_cleavages", "0"),
                elem.get("num_tol_term", "0"),
                elem.get("num_matched_peptides", "1"),
                str(len(elem.attrib["peptide"])),
                f"{exp_mass - calc_mass:.10g}",
                f"{abs(exp_mz - calc_mz):.10g}",
                f"{neg_log_expect:.10g}",
            ]
            fields.extend(f"{scores.get(score, 0.0):.10g}" for score in score_names)
            fields.extend("1" if charge == z else "0" for z in range(1, 8))
            fields.extend([modified_peptide(elem), *proteins])
            yield fields
            elem.clear()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--decoy-prefix", default="rev_")
    args = parser.parse_args()

    score_names = scan_score_names(args.input)
    header = [
        "SpecId",
        "Label",
        "ScanNr",
        "ExpMass",
        "CalcMass",
        "retentiontime",
        "rank",
        "missed_cleavages",
        "ntt",
        "num_matched_peptides",
        "peptide_length",
        "mass_diff",
        "abs_mz_diff",
        "NegLog10Expect",
        *score_names,
        *(f"Charge{z}" for z in range(1, 8)),
        "Peptide",
        "Proteins",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as out:
        out.write("\t".join(header) + "\n")
        for row in rows(args.input, args.decoy_prefix, score_names):
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
