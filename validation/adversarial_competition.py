#!/usr/bin/env python3
"""End-to-end falsification probe for exact ties in PSM competition.

Every synthetic spectrum has one target and one decoy candidate with identical
features.  The multiset of rows is the same in every arm; only the order of the
rows in the file changes.  A label-symmetric competition must not turn that
metadata choice into hundreds of confident target discoveries, and it must not
change the reported result at all.

Arms:

``target_first`` / ``decoy_first``
    The original two arms: within each spectrum, the target row precedes the
    decoy row or the other way round.

``file_reversed``
    The whole file reversed.

``targets_grouped`` / ``decoys_grouped``
    All targets before all decoys, and the reverse.

``shuffle_N``
    Deterministic whole-file shuffles.

``wide_ties``
    Three tied target candidates and one tied decoy candidate per spectrum, to
    check that a tie wider than two is drawn uniformly rather than being handed
    to whichever label the rule happens to favour.

Run against two builds to see a repair:

    python3 validation/adversarial_competition.py --binary OLD --json old.json
    python3 validation/adversarial_competition.py --binary NEW --json new.json
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path


HEADER = ("SpecId", "Label", "ScanNr", "ExpMass", "constant", "Peptide", "Proteins")


def candidates(spectra: int, targets_per_spectrum: int) -> list[tuple]:
    """One row per candidate, in canonical (scan, target-then-decoy) order."""
    rows = []
    for scan in range(1, spectra + 1):
        for index in range(targets_per_spectrum):
            suffix = "" if targets_per_spectrum == 1 else f"_{index}"
            rows.append(
                (f"scan{scan}_TARGET{suffix}", 1, scan, 500.0, 0.0,
                 f"K.TARGETPEP{scan}{suffix}.R", f"TARGET_P{scan}{suffix}")
            )
        rows.append(
            (f"scan{scan}_DECOY", -1, scan, 500.0, 0.0,
             f"K.DECOYPEP{scan}.R", f"DECOY_P{scan}")
        )
    return rows


def permute(rows: list[tuple], arm: str) -> list[tuple]:
    if arm == "target_first":
        return list(rows)
    if arm == "decoy_first":
        out = []
        by_scan: dict[int, list[tuple]] = {}
        for row in rows:
            by_scan.setdefault(row[2], []).append(row)
        for scan in sorted(by_scan):
            out.extend(reversed(by_scan[scan]))
        return out
    if arm == "file_reversed":
        return list(reversed(rows))
    if arm == "targets_grouped":
        return [r for r in rows if r[1] > 0] + [r for r in rows if r[1] < 0]
    if arm == "decoys_grouped":
        return [r for r in rows if r[1] < 0] + [r for r in rows if r[1] > 0]
    if arm.startswith("shuffle_"):
        seed = int(arm.split("_", 1)[1])
        out = list(rows)
        state = max(seed, 1)
        mask = (1 << 64) - 1
        for index in range(len(out) - 1, 0, -1):
            state ^= (state << 13) & mask
            state ^= state >> 7
            state ^= (state << 17) & mask
            state &= mask
            other = state % (index + 1)
            out[index], out[other] = out[other], out[index]
        return out
    raise ValueError(arm)


def write_pin(path: Path, rows: list[tuple]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)


def load(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(binary: Path, root: Path, name: str, rows: list[tuple], seed: int) -> dict[str, object]:
    pin = root / f"{name}.pin"
    targets = root / f"{name}.target.tsv"
    decoys = root / f"{name}.decoy.tsv"
    write_pin(pin, rows)
    command = [
        str(binary), "--canonical", "--no-select-c",
        "--seed", str(seed), "--num-threads", "1",
        "--results-psms", str(targets),
        "--decoy-results-psms", str(decoys),
        str(pin),
    ]
    execution = subprocess.run(command, text=True, capture_output=True, check=False)
    if execution.returncode:
        raise RuntimeError(f"{name} failed:\n{execution.stderr}")
    target_rows = load(targets)
    decoy_rows = load(decoys)
    return {
        "arm": name,
        "target_winners": len(target_rows),
        "decoy_winners": len(decoy_rows),
        "targets_q_lt_0.01": sum(float(row["q-value"]) < 0.01 for row in target_rows),
        "minimum_target_q": min((float(row["q-value"]) for row in target_rows), default=None),
        # Identity of the surviving set, so two arms can be compared exactly.
        "winner_signature": sorted(row["peptide"] for row in target_rows + decoy_rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=Path("target/release/percolator-rs"))
    parser.add_argument("--spectra", type=int, default=200)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"binary not found: {binary}")
    if args.spectra < 101:
        parser.error("use at least 101 spectra so the +1 correction can cross q<0.01")

    arms = [
        "target_first", "decoy_first", "file_reversed",
        "targets_grouped", "decoys_grouped",
        "shuffle_1", "shuffle_2", "shuffle_3", "shuffle_4",
    ]
    results = []
    with tempfile.TemporaryDirectory(prefix="percolator-competition-audit-") as temporary:
        root = Path(temporary)
        pairs = candidates(args.spectra, 1)
        for arm in arms:
            results.append(run(binary, root, arm, permute(pairs, arm), args.seed))
        wide = candidates(args.spectra * 5, 3)
        results.append(run(binary, root, "wide_ties", wide, args.seed))
        # Same fixture, different seed: a fair coin must re-flip.
        results.append(run(binary, root, "target_first_seed2", pairs, args.seed + 1))

    reference = results[0]["winner_signature"]
    permutation_arms = [r for r in results if r["arm"] in arms]
    invariant = all(r["winner_signature"] == reference for r in permutation_arms)
    distinct_winner_sets = len({tuple(r["winner_signature"]) for r in permutation_arms})

    for result in results:
        print(
            "\t".join(
                f"{key}={value}" for key, value in result.items() if key != "winner_signature"
            )
        )
    print(
        f"PERMUTATION_INVARIANT={invariant}\t"
        f"distinct_winner_sets_across_{len(permutation_arms)}_permutations={distinct_winner_sets}"
    )
    if args.json:
        args.json.write_text(
            json.dumps(
                {
                    "binary": str(binary),
                    "spectra": args.spectra,
                    "seed": args.seed,
                    "permutation_invariant": invariant,
                    "distinct_winner_sets": distinct_winner_sets,
                    "results": results,
                },
                indent=1,
            )
        )


if __name__ == "__main__":
    main()
