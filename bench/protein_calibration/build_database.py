#!/usr/bin/env python3
"""Build the exact PrEST target/decoy database and its ground-truth table."""

from __future__ import annotations

import argparse
from pathlib import Path


EXPECTED = {"A": 192, "B": 191, "RANDOM": 1000}


def read_fasta(path: Path):
    name = None
    description = ""
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, description, "".join(sequence)
                fields = line[1:].split(maxsplit=1)
                name = fields[0]
                description = fields[1] if len(fields) == 2 else ""
                sequence = []
            else:
                if name is None:
                    raise ValueError(f"sequence before FASTA header in {path}")
                sequence.append(line.upper())
    if name is not None:
        yield name, description, "".join(sequence)


def write_entry(handle, name: str, description: str, sequence: str) -> None:
    suffix = f" {description}" if description else ""
    handle.write(f">{name}{suffix}\n")
    for start in range(0, len(sequence), 60):
        handle.write(sequence[start : start + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool-a", required=True, type=Path)
    parser.add_argument("--pool-b", required=True, type=Path)
    parser.add_argument("--random", required=True, type=Path)
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    args = parser.parse_args()

    records: list[tuple[str, str, str, str]] = []
    seen: set[str] = set()
    for pool, path in (("A", args.pool_a), ("B", args.pool_b), ("RANDOM", args.random)):
        pool_records = list(read_fasta(path))
        if len(pool_records) != EXPECTED[pool]:
            raise ValueError(
                f"{path} has {len(pool_records)} entries; expected {EXPECTED[pool]}"
            )
        for name, description, sequence in pool_records:
            if name in seen:
                raise ValueError(f"duplicate protein identifier: {name}")
            if not sequence or any(aa not in "ABCDEFGHIKLMNPQRSTVWXYZJUO*" for aa in sequence):
                raise ValueError(f"invalid sequence for {name}")
            seen.add(name)
            records.append((name, description, sequence, pool))

    args.database.parent.mkdir(parents=True, exist_ok=True)
    with args.database.open("w", encoding="utf-8") as database:
        for name, description, sequence, _pool in records:
            write_entry(database, name, description, sequence)
        # The original benchmark used reversed proteins. Explicit paired names make
        # picked target-decoy competition deterministic in both inference methods.
        for name, description, sequence, _pool in records:
            write_entry(database, f"DECOY_{name}", description, sequence[::-1])

    with args.truth.open("w", encoding="utf-8") as truth:
        truth.write("protein_id\tpool\n")
        for name, _description, _sequence, pool in records:
            truth.write(f"{name}\t{pool}\n")

    residues = sum(len(sequence) for _name, _desc, sequence, _pool in records)
    print(f"target_proteins={len(records)}")
    print(f"decoy_proteins={len(records)}")
    print(f"target_residues={residues}")
    for pool in ("A", "B", "RANDOM"):
        print(f"pool_{pool.lower()}={sum(record[3] == pool for record in records)}")


if __name__ == "__main__":
    main()
