#!/usr/bin/env python3
"""Assemble result-blind database characterization and collision checks."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from construct_databases import canonical, read_fasta, theoretical_peptides


CONDITIONS = (
    "original", "homology_depleted", "size_control_130363",
    "size_control_155921", "size_control_196613",
)


def read_removed(path: Path) -> set[int]:
    with path.open() as handle:
        next(handle)
        return {int(line.split("\t", 1)[0]) for line in handle if line.strip()}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fraction_counts(values: list[int]) -> list[float]:
    total = sum(values)
    return [value / total if total else 0.0 for value in values]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    construction = json.loads((args.database_root / "construction_manifest.json").read_text())
    characterization = json.loads((args.database_root / "presearch_characterization.json").read_text())
    records = read_fasta(args.database_root / "original.fasta")
    native = [record for record in records if not record.header.startswith("ENT_")]
    entrapment = [record for record in records if record.header.startswith("ENT_")]
    removed = {"original": set()}
    for condition in CONDITIONS[1:]:
        removed[condition] = read_removed(args.database_root / f"{condition}.removed.tsv")

    native_protein_sequences = {canonical(record.sequence) for record in native}
    native_peptides: dict[int, set[str]] = defaultdict(set)
    for record in native:
        for peptide, *_ in theoretical_peptides(record.sequence):
            native_peptides[len(peptide)].add(peptide)

    exact_shared_instances = Counter()
    exact_shared_sequences = {condition: set() for condition in CONDITIONS}
    protein_collisions = Counter()
    for index, record in enumerate(entrapment):
        active = [condition for condition in CONDITIONS if index not in removed[condition]]
        if canonical(record.sequence) in native_protein_sequences:
            for condition in active:
                protein_collisions[condition] += 1
        for peptide, *_ in theoretical_peptides(record.sequence):
            if peptide in native_peptides.get(len(peptide), ()):
                for condition in active:
                    exact_shared_instances[condition] += 1
                    exact_shared_sequences[condition].add(peptide)

    witness_by_index = {row["protein_index"]: row for row in construction["homology_witnesses"]}
    similarity = {}
    for condition in CONDITIONS:
        kept_witnesses = [row for index, row in witness_by_index.items() if index not in removed[condition]]
        similarity[condition] = {
            "entrapment_proteins": len(entrapment) - len(removed[condition]),
            "proteins_with_primary_near_homolog": len(kept_witnesses),
            "fraction_with_primary_near_homolog": (
                len(kept_witnesses) / (len(entrapment) - len(removed[condition]))
                if len(entrapment) != len(removed[condition]) else None
            ),
            "first_witness_distance_0_1_2": [
                sum(row["distance"] == distance for row in kept_witnesses) for distance in range(3)
            ],
            "first_witness_length_counts": dict(sorted(Counter(
                len(row["entrapment_peptide"]) for row in kept_witnesses
            ).items())),
        }

    conditions = {}
    original_ent = characterization["conditions"]["original"]["entrapment_component"]
    for condition in CONDITIONS:
        database_path = args.database_root / f"{condition}.fasta"
        char = characterization["conditions"][condition]
        ent = char["entrapment_component"]
        total = char["complete_target_database"]
        header_seen: set[str] = set()
        duplicate_headers = 0
        for record in read_fasta(database_path):
            accession = record.header.split()[0]
            duplicate_headers += accession in header_seen
            header_seen.add(accession)
        amino = ent["amino_acid_counts_A_to_Z"]
        conditions[condition] = {
            "database": construction["databases"][condition],
            "duplicate_accession_headers": duplicate_headers,
            "il_canonical_full_protein_collisions_with_native": protein_collisions[condition],
            "exact_shared_full_tryptic_peptide_instances": exact_shared_instances[condition],
            "exact_shared_full_tryptic_unique_sequences": len(exact_shared_sequences[condition]),
            "total_component": total,
            "entrapment_component": ent,
            "entrapment_opportunity_retained": {
                "proteins": ent["proteins"] / original_ent["proteins"],
                "residues": ent["residues"] / original_ent["residues"],
                "fully_tryptic_instances": ent["fully_tryptic_peptide_instances"] / original_ent["fully_tryptic_peptide_instances"],
                "searchable_semi_tryptic_instances": ent["searchable_semi_tryptic_peptide_instances"] / original_ent["searchable_semi_tryptic_peptide_instances"],
                "searchable_unique_hll": ent["searchable_unique_hll"] / original_ent["searchable_unique_hll"],
            },
            "entrapment_amino_acid_fraction_A_to_Z": fraction_counts(amino),
            "similarity_to_native": similarity[condition],
        }

    result = {
        "schema_version": 1,
        "presearch_only": True,
        "reads_search_or_percolator_results": False,
        "preregistration_sha256": sha256(args.database_root.parent / "PREREGISTRATION.md"),
        "construction_script_sha256": sha256(args.database_root.parent / "construct_databases.py"),
        "characterization_source_sha256": sha256(args.database_root.parent / "characterize_databases.cpp"),
        "characterization_binary_sha256": sha256(args.database_root.parent / "characterize_databases"),
        "definitions": {
            "primary_homology": construction["homology_rule"],
            "searchable_peptide": characterization["searchable_definition"],
            "unique_estimator": characterization["unique_method"],
            "collision_check": "I/L-canonical full protein equality and exact I/L-canonical full-tryptic peptide equality in the primary theoretical universe",
        },
        "conditions": conditions,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        condition: {
            "proteins": row["entrapment_component"]["proteins"],
            "searchable_instances": row["entrapment_component"]["searchable_semi_tryptic_peptide_instances"],
            "searchable_unique_hll": row["entrapment_component"]["searchable_unique_hll"],
            "primary_near_homolog_proteins": row["similarity_to_native"]["proteins_with_primary_near_homolog"],
            "exact_shared_tryptic_sequences": row["exact_shared_full_tryptic_unique_sequences"],
            "protein_collisions": row["il_canonical_full_protein_collisions_with_native"],
            "duplicate_headers": row["duplicate_accession_headers"],
        } for condition, row in conditions.items()
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
