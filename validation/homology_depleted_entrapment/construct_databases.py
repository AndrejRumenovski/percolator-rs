#!/usr/bin/env python3
"""Construct preregistered protein-level homology-depleted target FASTAs.

This is validation-only code. It never reads search or Percolator results.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np


CONTROL_SEEDS = (130363, 155921, 196613)
STANDARD = frozenset("ACDEFGHIKLMNPQRSTVWY")
MIN_LENGTH = 8
MAX_LENGTH = 63
MIN_MH = 600.0
MAX_MH = 5000.0
MIN_IDENTITY = 0.85
MAX_DISTANCE = 2
STATIC_C = 57.021464
WATER_PLUS_PROTON = 19.018389715
MASS = {
    "A": 71.037113805, "R": 156.101111050, "N": 114.042927470,
    "D": 115.026943065, "C": 103.009184505 + STATIC_C,
    "E": 129.042593135, "Q": 128.058577540, "G": 57.021463735,
    "H": 137.058911875, "I": 113.084063975, "L": 113.084063975,
    "K": 128.094963015, "M": 131.040484645, "F": 147.068413945,
    "P": 97.052763875, "S": 87.032028435, "T": 101.047678505,
    "W": 186.079312980, "Y": 163.063328575, "V": 99.068413945,
}


@dataclass(frozen=True)
class Record:
    header: str
    sequence: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> list[Record]:
    records: list[Record] = []
    header: str | None = None
    sequence: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append(Record(header, "".join(sequence).upper()))
                header, sequence = line[1:], []
            elif line:
                sequence.append(line)
    if header is not None:
        records.append(Record(header, "".join(sequence).upper()))
    return records


def write_fasta(path: Path, records: list[Record]) -> None:
    with path.open("w") as handle:
        for record in records:
            handle.write(f">{record.header}\n{record.sequence}\n")


def canonical(peptide: str) -> str:
    return peptide.replace("I", "L")


def cleavage_boundaries(sequence: str) -> list[int]:
    boundaries = [0]
    for index, residue in enumerate(sequence):
        if residue in "KR" and (index + 1 == len(sequence) or sequence[index + 1] != "P"):
            boundaries.append(index + 1)
    if boundaries[-1] != len(sequence):
        boundaries.append(len(sequence))
    return boundaries


def theoretical_peptides(sequence: str):
    boundaries = cleavage_boundaries(sequence)
    for left in range(len(boundaries) - 1):
        for missed in range(3):
            right_index = left + missed + 1
            if right_index >= len(boundaries):
                break
            start, end = boundaries[left], boundaries[right_index]
            peptide = sequence[start:end]
            length = len(peptide)
            if length < MIN_LENGTH or length > MAX_LENGTH:
                continue
            if not set(peptide) <= STANDARD:
                continue
            mh = WATER_PLUS_PROTON + sum(MASS[residue] for residue in peptide)
            if MIN_MH <= mh <= MAX_MH:
                yield canonical(peptide), start, end, mh, missed


def allowed_distance(length: int) -> int:
    by_identity = math.floor((1.0 - MIN_IDENTITY) * length + 1e-12)
    return min(MAX_DISTANCE, by_identity)


def anchors(sequence: str, pieces: int) -> list[tuple[int, str]]:
    return [
        (piece, sequence[(piece * len(sequence)) // pieces:((piece + 1) * len(sequence)) // pieces])
        for piece in range(pieces)
    ]


def native_peptide_sets(records: list[Record]) -> dict[int, set[str]]:
    result: dict[int, set[str]] = defaultdict(set)
    for number, record in enumerate(records, 1):
        for peptide, *_ in theoretical_peptides(record.sequence):
            result[len(peptide)].add(peptide)
        if number % 10000 == 0:
            print(f"native digest {number}/{len(records)}", file=sys.stderr, flush=True)
    return dict(result)


def build_indices(native: dict[int, set[str]]):
    indices = {}
    sequences = {}
    for length in sorted(native):
        dmax = allowed_distance(length)
        if dmax < 0:
            continue
        pieces = dmax + 1
        seqs = list(native[length])
        index: dict[tuple[int, str], list[int]] = defaultdict(list)
        for sequence_id, sequence in enumerate(seqs):
            for key in anchors(sequence, pieces):
                index[key].append(sequence_id)
        indices[length] = (pieces, dict(index))
        sequences[length] = seqs
        print(
            f"native index length={length} unique={len(seqs)} anchors={len(index)} dmax={dmax}",
            file=sys.stderr,
            flush=True,
        )
    return sequences, indices


def find_witness(query: str, sequences, indices):
    length = len(query)
    if length not in indices:
        return None
    pieces, index = indices[length]
    seen: set[int] = set()
    dmax = allowed_distance(length)
    for key in anchors(query, pieces):
        for sequence_id in index.get(key, ()):
            if sequence_id in seen:
                continue
            seen.add(sequence_id)
            native = sequences[length][sequence_id]
            distance = sum(left != right for left, right in zip(query, native))
            if distance <= dmax:
                return {"entrapment_peptide": query, "native_peptide": native,
                        "distance": distance, "identity": 1.0 - distance / length}
    return None


def source_tag(header: str) -> str:
    first = header.split()[0]
    parts = first.split("_")
    return parts[1] if len(parts) > 2 and parts[0] == "ENT" else "UNKNOWN"


def length_stratum(sequence: str) -> int:
    return min(len(sequence) // 50, 40)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--combined", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    records = read_fasta(args.combined)
    native = [record for record in records if not record.header.startswith("ENT_")]
    entrapment = [record for record in records if record.header.startswith("ENT_")]
    if len(native) != 139191 or len(entrapment) != 389504:
        raise ValueError(f"unexpected source counts native={len(native)} entrapment={len(entrapment)}")

    native_sets = native_peptide_sets(native)
    native_sequences, native_indices = build_indices(native_sets)
    removed_homology: set[int] = set()
    witnesses = []
    for index, record in enumerate(entrapment):
        witness = None
        for peptide, start, end, mh, missed in theoretical_peptides(record.sequence):
            hit = find_witness(peptide, native_sequences, native_indices)
            if hit is not None:
                witness = {"protein_index": index, "protein": record.header,
                           "start_0based": start, "end_0based_exclusive": end,
                           "mh_plus": mh, "missed_cleavages": missed, **hit}
                break
        if witness is not None:
            removed_homology.add(index)
            witnesses.append(witness)
        if (index + 1) % 10000 == 0:
            print(
                f"entrapment scan {index + 1}/{len(entrapment)} removed={len(removed_homology)}",
                file=sys.stderr,
                flush=True,
            )

    strata: dict[tuple[str, int], list[int]] = defaultdict(list)
    removed_by_stratum: dict[tuple[str, int], int] = defaultdict(int)
    for index, record in enumerate(entrapment):
        stratum = (source_tag(record.header), length_stratum(record.sequence))
        strata[stratum].append(index)
        if index in removed_homology:
            removed_by_stratum[stratum] += 1

    removals: dict[str, set[int]] = {"homology_depleted": removed_homology}
    for seed in CONTROL_SEEDS:
        rng = np.random.default_rng(seed)
        selected: set[int] = set()
        for stratum in sorted(strata):
            candidates = np.asarray(strata[stratum], dtype=np.int64)
            count = removed_by_stratum.get(stratum, 0)
            if count:
                chosen = rng.choice(candidates, size=count, replace=False)
                selected.update(int(value) for value in chosen)
        removals[f"size_control_{seed}"] = selected

    conditions: dict[str, list[Record]] = {"original": records}
    for name, removed in removals.items():
        kept = [record for index, record in enumerate(entrapment) if index not in removed]
        conditions[name] = native + kept

    database_manifest = {}
    for name, condition_records in conditions.items():
        path = args.output_root / f"{name}.fasta"
        write_fasta(path, condition_records)
        ent_count = sum(record.header.startswith("ENT_") for record in condition_records)
        database_manifest[name] = {
            "path": str(path), "sha256": sha256(path),
            "native_proteins": len(condition_records) - ent_count,
            "entrapment_proteins": ent_count,
            "removed_entrapment_proteins": len(entrapment) - ent_count,
            "bytes": path.stat().st_size,
        }

    for name, removed in removals.items():
        path = args.output_root / f"{name}.removed.tsv"
        with path.open("w") as handle:
            handle.write("protein_index\theader\tlength\tsource\tlength_stratum\n")
            for index in sorted(removed):
                record = entrapment[index]
                handle.write(
                    f"{index}\t{record.header}\t{len(record.sequence)}\t"
                    f"{source_tag(record.header)}\t{length_stratum(record.sequence)}\n"
                )

    manifest = {
        "schema_version": 1,
        "reads_search_or_percolator_results": False,
        "source": {"path": str(args.combined), "sha256": sha256(args.combined),
                   "native_proteins": len(native), "entrapment_proteins": len(entrapment)},
        "homology_rule": {
            "digest": "trypsin K/R not before P; 0-2 missed cleavages",
            "peptide_length": [MIN_LENGTH, MAX_LENGTH],
            "mh_plus_da": [MIN_MH, MAX_MH],
            "static_cysteine_da": STATIC_C,
            "il_canonicalization": "I->L",
            "maximum_hamming_distance": MAX_DISTANCE,
            "minimum_identity": MIN_IDENTITY,
            "depletion_unit": "whole entrapment protein",
        },
        "control": {
            "rng": "NumPy Generator(PCG64)", "seeds": list(CONTROL_SEEDS),
            "strata": "source proteome x fixed 50-aa length bins; >=2000 final bin",
        },
        "native_unique_theoretical_peptides_by_length": {
            str(length): len(values) for length, values in sorted(native_sets.items())
        },
        "homology_removed_count": len(removed_homology),
        "homology_witnesses": witnesses,
        "databases": database_manifest,
    }
    manifest_path = args.output_root / "construction_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"manifest": str(manifest_path),
                      "homology_removed_count": len(removed_homology),
                      "databases": database_manifest}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
