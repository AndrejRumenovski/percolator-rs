#!/usr/bin/env python3
"""Mutation audit for the two latest repairs, run in a throwaway worktree.

Each mutant removes one repaired property.  It is successful only when the
repository's ordinary validation suite fails.  The caller is responsible for
supplying a detached worktree; this script never touches the primary checkout.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path


MUTATIONS = [
    {
        "name": "joined_rows_not_canonicalized",
        "file": "src/pin.rs",
        "changes": [
            (
                """    for part in &mut parts {
        canonicalize_rows(part);
    }
""",
                """    // MUTATION: retain source row insertion order.
""",
            )
        ],
        "test": ["cargo", "test", "--release", "--test", "join_permutation"],
    },
    {
        "name": "joined_files_not_canonicalized",
        "file": "src/pin.rs",
        "changes": [
            (
                "    parts.sort_by(compare_parts);\n",
                "    // MUTATION: retain file-argument order.\n",
            )
        ],
        "test": ["cargo", "test", "--release", "--test", "join_permutation"],
    },
    {
        "name": "protein_representative_mapping_only",
        "file": "src/main.rs",
        "changes": [
            (
                "        let entries = protein_entries(&ds, &reported_indices, &pep_idx, &pscore, &ppep);\n",
                """        // MUTATION: keep only the selected representative PSM mapping.
        let entries: Vec<(f64, f64, String)> = pep_idx
            .iter()
            .enumerate()
            .map(|(peptide, &index)| {
                (pscore[peptide], ppep[peptide], ds.proteins[index].clone())
            })
            .collect();
""",
            )
        ],
        "test": ["cargo", "test", "--release", "--test", "protein_grouping"],
    },
    {
        "name": "protein_target_decoy_class_collapsed",
        "file": "src/protein.rs",
        "changes": [
            (
                "    let mut group_of: HashMap<(bool, Vec<u32>), usize> = HashMap::new();\n",
                "    let mut group_of: HashMap<Vec<u32>, usize> = HashMap::new();\n",
            ),
            (
                "    let mut group_evidence: Vec<(bool, Vec<u32>)> = Vec::new();\n",
                "    let mut group_evidence: Vec<Vec<u32>> = Vec::new();\n",
            ),
            (
                """        let key = (
            is_decoy_protein(names[protein]),
            std::mem::take(&mut evidence[protein]),
        );
""",
                """        // MUTATION: evidence alone defines indistinguishability.
        let key = std::mem::take(&mut evidence[protein]);
""",
            ),
            (
                "        .map(|(members, (is_decoy, peptides))| {\n",
                "        .map(|(members, peptides)| {\n",
            ),
            (
                "            proteins.sort();\n            let score = peptides\n",
                """            proteins.sort();
            // A mixed group can no longer be a decoy competitor.
            let is_decoy = proteins.iter().all(|protein| is_decoy_protein(protein));
            let score = peptides
""",
            ),
            ("                is_decoy: *is_decoy,\n", "                is_decoy,\n"),
        ],
        "test": [
            "cargo", "test", "--release", "--bin", "percolator-rs",
            "protein::tests",
        ],
    },
    {
        "name": "ensemble_agreement_uses_heldout_label",
        "file": "src/pin.rs",
        "changes": [
            (
                "    let mut psm_engines: BTreeMap<(i64, &str), BTreeSet<u32>> = BTreeMap::new();\n",
                "    let mut psm_engines: BTreeMap<(i64, i8, &str), BTreeSet<u32>> = BTreeMap::new();\n",
            ),
            (
                "            .entry((out.scan[row], out.peptide[row].as_str()))\n",
                "            .entry((out.scan[row], out.labels[row], out.peptide[row].as_str()))\n",
            ),
            (
                "                psm_engines[&(out.scan[row], out.peptide[row].as_str())].len() as f64,\n",
                "                psm_engines[&(out.scan[row], out.labels[row], out.peptide[row].as_str())].len() as f64,\n",
            ),
        ],
        "test": [
            "cargo", "test", "--release", "--bin", "percolator-rs",
            "pin::tests::ensemble_features_do_not_depend_on_any_label",
        ],
    },
]


def run(command: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)


def failures(output: str) -> list[str]:
    return sorted(set(re.findall(r"^test (\S+) \.\.\. FAILED", output, re.MULTILINE)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worktree", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    args = parser.parse_args()
    worktree = args.worktree.resolve()

    baseline = run(["cargo", "test", "--release", "--all-targets"], worktree)
    if baseline.returncode != 0:
        raise SystemExit("mutation baseline did not pass")

    results = []
    for mutation in MUTATIONS:
        path = worktree / mutation["file"]
        original = path.read_text()
        changed = original
        for old, new in mutation["changes"]:
            count = changed.count(old)
            if count != 1:
                raise RuntimeError(
                    f"{mutation['name']}: expected one fragment, observed {count}"
                )
            changed = changed.replace(old, new, 1)
        path.write_text(changed)
        try:
            result = run(mutation["test"], worktree)
            combined = result.stdout + result.stderr
            compile_failed = "could not compile" in combined or "error[E" in combined
            failed_tests = failures(combined)
            caught = result.returncode != 0 and bool(failed_tests) and not compile_failed
            record = {
                "mutation": mutation["name"],
                "command": mutation["test"],
                "returncode": result.returncode,
                "compile_failed": compile_failed,
                "failing_tests": failed_tests,
                "caught_by_behavioral_test": caught,
                "stdout_tail": result.stdout[-3000:],
                "stderr_tail": result.stderr[-3000:],
            }
            results.append(record)
            status = "CAUGHT" if caught else "MISSED"
            print(f"{mutation['name']}: {status} {failed_tests}", flush=True)
        finally:
            path.write_text(original)

    payload = {
        "baseline_passed": True,
        "mutations": results,
        "caught": sum(item["caught_by_behavioral_test"] for item in results),
        "total": len(results),
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"caught {payload['caught']}/{payload['total']} by behavioral failures")


if __name__ == "__main__":
    main()
