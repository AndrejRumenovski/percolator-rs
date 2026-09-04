#!/usr/bin/env python3
"""Bounded near-homology audit for high-confidence entrapment PSMs.

This is a validation-only exact algorithm.  It asks whether an entrapment
peptide has a same-length native-proteome substring at Hamming distance <= 2,
with I/L treated as mass-spectrometrically indistinguishable.  Three disjoint
anchors make the search exhaustive for that distance bound: with at most two
mismatches at least one anchor must match exactly.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict, deque
from pathlib import Path

import numpy as np

from pep_rootcause_experiments import load_outputs


def canonical(sequence: str) -> str:
    return sequence.replace("I", "L")


def anchors(sequence: str) -> list[tuple[int, str]]:
    n = len(sequence)
    cuts = (0, n // 3, (2 * n) // 3, n)
    return [(cuts[i], sequence[cuts[i]:cuts[i + 1]]) for i in range(3)]


def build_automaton(queries: list[str]):
    goto = [{}]
    fail = [0]
    output = [[]]
    for qi, query in enumerate(queries):
        for offset, word in anchors(query):
            state = 0
            for character in word:
                if character not in goto[state]:
                    goto[state][character] = len(goto)
                    goto.append({}); fail.append(0); output.append([])
                state = goto[state][character]
            output[state].append((qi, offset, len(word)))
    queue = deque(goto[0].values())
    while queue:
        state = queue.popleft()
        for character, nxt in goto[state].items():
            queue.append(nxt)
            fallback = fail[state]
            while fallback and character not in goto[fallback]:
                fallback = fail[fallback]
            fail[nxt] = goto[fallback].get(character, 0)
            output[nxt].extend(output[fail[nxt]])
    return goto, fail, output


def native_sequences(path: Path):
    header = None
    sequence = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None and not header.startswith("ENT_") and not header.startswith("DECOY_"):
                    yield header, canonical("".join(sequence))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line.strip().upper())
    if header is not None and not header.startswith("ENT_") and not header.startswith("DECOY_"):
        yield header, canonical("".join(sequence))


def search(fasta: Path, queries: list[str]) -> list[dict | None]:
    goto, fail, output = build_automaton(queries)
    best: list[dict | None] = [None] * len(queries)
    for header, sequence in native_sequences(fasta):
        state = 0
        for end, character in enumerate(sequence):
            while state and character not in goto[state]:
                state = fail[state]
            state = goto[state].get(character, 0)
            for qi, offset, anchor_length in output[state]:
                start = end + 1 - anchor_length - offset
                query = queries[qi]
                if start < 0 or start + len(query) > len(sequence):
                    continue
                distance = sum(a != b for a, b in zip(query, sequence[start:start + len(query)]))
                if distance <= 2 and (best[qi] is None or distance < best[qi]["distance"]):
                    best[qi] = {
                        "distance": distance,
                        "native_substring": sequence[start:start + len(query)],
                        "native_header": header,
                        "actin_or_tubulin_header": "actin" in header.lower() or "tubulin" in header.lower(),
                    }
    return best


def summarize_group(data, mask, hits_by_peptide):
    peptides = sorted(set(str(x) for x in data["peptide"][mask] if x))
    hit = {p: hits_by_peptide[p] for p in peptides}
    row_hit = np.array([hit.get(str(p)) is not None for p in data["peptide"]], dtype=bool)
    selected_hit = mask & row_hit
    selected_other = mask & ~row_hit
    distances = [value["distance"] for value in hit.values() if value is not None]
    families = [p for p, value in hit.items() if value is not None and value["actin_or_tubulin_header"]]
    def stats(which):
        return {
            "rows": int(which.sum()),
            "mean_pep": float(data["pep"][which].mean()) if which.any() else None,
            "mean_score": float(data["score"][which].mean()) if which.any() else None,
            "mean_peptide_length": float(np.mean([len(str(x)) for x in data["peptide"][which]])) if which.any() else None,
        }
    return {
        "rows": int(mask.sum()), "distinct_peptides": len(peptides),
        "distinct_with_native_distance_le_2": len(distances),
        "fraction_with_native_distance_le_2": len(distances) / len(peptides) if peptides else None,
        "distance_0_1_2": [distances.count(i) for i in range(3)],
        "actin_or_tubulin_near_homolog_peptides": families,
        "near_homolog_rows": stats(selected_hit),
        "other_rows": stats(selected_other),
        "matches": [dict(peptide=p, **value) for p, value in hit.items() if value is not None],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = load_outputs(args.results_root, "target.tsv", None)
    masks = {}
    for threshold in (0.01, 0.05):
        masks[(threshold, "entrapment_target")] = (~data["decoy"] & data["pure"] & (data["pep"] < threshold))
        masks[(threshold, "entrapment_decoy")] = (data["decoy"] & data["pure"] & (data["pep"] < threshold))
    peptides = sorted(set(
        canonical(str(peptide)) for mask in masks.values() for peptide in data["peptide"][mask] if peptide
    ))
    found = search(args.fasta, peptides)
    hits = dict(zip(peptides, found))
    # Canonicalize the peptide array in the summary lookup, too.
    data = dict(data); data["peptide"] = np.asarray([canonical(str(x)) for x in data["peptide"]], dtype=object)
    result = {
        "definition": "exhaustive native substring Hamming distance <=2 after I/L canonicalization",
        "query_peptides": len(peptides),
        "groups": [
            {"pep_threshold": threshold, "class": label,
             **summarize_group(data, mask, hits)}
            for (threshold, label), mask in masks.items()
        ],
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
