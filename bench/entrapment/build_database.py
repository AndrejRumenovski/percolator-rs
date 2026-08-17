#!/usr/bin/env python3
"""Build a native + equally sized foreign-proteome entrapment FASTA."""

import argparse
import gzip
from pathlib import Path


def records(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        header = None
        sequence = []
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:], []
            elif line:
                sequence.append(line)
        if header is not None:
            yield header, "".join(sequence)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("native", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("foreign", nargs="+", type=Path)
    args = parser.parse_args()

    native = list(records(args.native))
    native_aa = sum(len(seq) for _, seq in native)
    entrapment_aa = 0
    entrapment_count = 0

    with args.output.open("w") as out:
        for header, sequence in native:
            print(f">{header}", file=out)
            print(sequence, file=out)

        for source in args.foreign:
            source_tag = source.name.split(".")[0].upper()
            for header, sequence in records(source):
                if entrapment_aa >= native_aa:
                    break
                accession = header.split()[0].replace("|", "_")
                print(f">ENT_{source_tag}_{accession}", file=out)
                print(sequence, file=out)
                entrapment_aa += len(sequence)
                entrapment_count += 1
            if entrapment_aa >= native_aa:
                break

    ratio = entrapment_aa / (native_aa + entrapment_aa)
    print(f"native_proteins={len(native)} native_aa={native_aa}")
    print(f"entrapment_proteins={entrapment_count} entrapment_aa={entrapment_aa}")
    print(f"entrapment_fraction={ratio:.8f}")


if __name__ == "__main__":
    main()
