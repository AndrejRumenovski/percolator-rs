#!/usr/bin/env python3
"""Reintroduce known scientific defects and check that the suite catches them.

A test suite that passes is only evidence if it would have failed. Each mutation
below is a defect this repository actually shipped, or the exact inverse of a
correction it made. The suite must fail on every one of them; a mutation that
passes marks a property nothing constrains.

Mutations are applied in a throwaway git worktree and reverted afterwards.
Production source is never touched.

    python3 validation/mutation_test.py --worktree /path/to/worktree --json out.json
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path


# (name, file, old fragment, new fragment, what defect this recreates)
MUTATIONS: list[tuple[str, str, str, str, str]] = [
    (
        "input_order_dependent_psm_ties",
        "src/main.rs",
        """                } else if score[i] == current.score {
                    current.tied += 1;
                }""",
        """                } else if score[i] == current.score {
                    // MUTATION: first row wins every tie.
                }""",
        "PSM competition resolves exact target/decoy ties by row order",
    ),
    (
        "missing_finite_sample_decoy",
        "src/stats.rs",
        """    pub fn reported(null_target_win_prob: f64) -> Self {
        Tdc {
            pi0: 1.0,
            null_target_win_prob,
            skip_decoys_plus_one: false,
        }
    }""",
        """    pub fn reported(null_target_win_prob: f64) -> Self {
        Tdc {
            pi0: 1.0,
            null_target_win_prob,
            // MUTATION: drop the finite-sample safeguard on the reported path.
            skip_decoys_plus_one: true,
        }
    }""",
        "reported q-values lose the +1 decoy, so a leading target run reaches q = 0",
    ),
    (
        "broken_tie_grouping",
        "src/stats.rs",
        """fn ends_score_group(order: &[usize], scores: &[f64], rank: usize) -> bool {
    rank + 1 == order.len() || scores[order[rank + 1]] != scores[order[rank]]
}""",
        """fn ends_score_group(order: &[usize], scores: &[f64], rank: usize) -> bool {
    // MUTATION: every row is its own tie group.
    let _ = (order, scores);
    let _ = rank;
    true
}""",
        "equal scores stop sharing one rejection boundary",
    ),
    (
        "pep_prior_doubled",
        "src/stats.rs",
        """                let share =
                    ((estimated_false - assigned) / group_targets as f64).max(0.0);""",
        """                let share =
                    ((estimated_false - assigned) / group_targets as f64).max(0.0)
                        + 0.5 / target_hint(labels) as f64;""",
        "an arbitrary constant is added to every PEP after the estimate",
    ),
    (
        "pep_floor_raised",
        "src/stats.rs",
        "const PEP_FLOOR: f64 = 1e-12;",
        "const PEP_FLOOR: f64 = 0.01; // MUTATION",
        "small PEPs are clamped upward by a constant",
    ),
    (
        "leaked_normalization",
        "src/percolator.rs",
        "    let (x, dim) = build_matrix_fit(ds, &train_rows, p);\n    let w0 = initial_direction(&x, dim, &ds.labels, &train_rows, p);",
        "    // MUTATION: fit normalization on every row, held-out included.\n    let all_rows: Vec<usize> = (0..n).collect();\n    let (x, dim) = build_matrix_fit(ds, &all_rows, p);\n    let w0 = initial_direction(&x, dim, &ds.labels, &train_rows, p);",
        "per-fold normalization is fitted on the whole dataset",
    ),
    (
        "leaked_initial_direction",
        "src/percolator.rs",
        "    let w0 = initial_direction(&x, dim, &ds.labels, &train_rows, p);",
        "    // MUTATION: choose the direction using every label.\n    let every_row: Vec<usize> = (0..n).collect();\n    let w0 = initial_direction(&x, dim, &ds.labels, &every_row, p);",
        "the initial direction is chosen with held-out labels",
    ),
    (
        "leaked_c_selection",
        "src/percolator.rs",
        "        let (alpha, beta, inner_yield) = select_c_for_fold(ds, &setup.train_rows, test, p);",
        "        // MUTATION: select C on every row, held-out fold included.\n        let every_row: Vec<usize> = (0..ds.n_psm).collect();\n        let (alpha, beta, inner_yield) = select_c_for_fold(ds, &every_row, test, p);",
        "class-weight selection sees the fold it will score",
    ),
    (
        "leaked_ensemble_feature",
        "src/pin.rs",
        """    let mut psm_engines: BTreeMap<(i64, &str), BTreeSet<u32>> = BTreeMap::new();""",
        """    // MUTATION: key the agreement feature on the label again.
    let mut psm_engines: BTreeMap<(i64, i8, &str), BTreeSet<u32>> = BTreeMap::new();""",
        "the cross-engine agreement feature is built from labels",
    ),
    (
        "target_favouring_protein_ties",
        "src/protein.rs",
        """                let target_wins = if groups[ti].score == groups[di].score {
                    Coin::new(seed).bytes(key.as_bytes()).heads()
                } else {
                    groups[ti].score > groups[di].score
                };""",
        """                // MUTATION: the target wins every exact tie.
                let target_wins = groups[ti].score >= groups[di].score;""",
        "picked-protein ties always go to the target",
    ),
    (
        "connected_component_protein_grouping",
        "src/protein.rs",
        """    let mut group_of: HashMap<Vec<u32>, usize> = HashMap::new();""",
        """    // MUTATION: group by connected components of peptide sharing, which
    // merges proteins whose evidence differs.
    {
        let mut parent: Vec<usize> = (0..names.len()).collect();
        fn find(parent: &mut Vec<usize>, mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        for (_, _, raw) in entries {
            let members: Vec<usize> = split_proteins(raw).iter().map(|p| id_of[*p]).collect();
            for pair in members.windows(2) {
                let (a, b) = (find(&mut parent, pair[0]), find(&mut parent, pair[1]));
                if a != b {
                    parent[a] = b;
                }
            }
        }
        for protein in 0..names.len() {
            let root = find(&mut parent, protein);
            evidence[protein] = vec![root as u32];
        }
    }
    let mut group_of: HashMap<Vec<u32>, usize> = HashMap::new();""",
        "distinguishable proteins are merged by shared-peptide connectivity",
    ),
    (
        "peptide_pep_as_protein_pep",
        "src/protein.rs",
        """                // Picked-protein FDR estimates no protein-level posterior.
                pep: None,""",
        """                // MUTATION: report the best peptide's PEP as the protein PEP.
                pep: Some(
                    peptides
                        .iter()
                        .map(|&peptide| entries[peptide as usize].1)
                        .fold(1.0f64, f64::min),
                ),""",
        "a peptide-level PEP is emitted under a protein-level column name",
    ),
]


HELPER = """
#[cfg(test)]
fn target_hint(labels: &[i8]) -> usize {
    labels.iter().filter(|&&label| label > 0).count().max(1)
}
#[cfg(not(test))]
fn target_hint(labels: &[i8]) -> usize {
    labels.iter().filter(|&&label| label > 0).count().max(1)
}
"""


def run(command: list[str], cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)


def failing_tests(output: str) -> list[str]:
    return sorted({name for name in re.findall(r"^test (\S+) \.\.\. FAILED", output, re.M)})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worktree", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    worktree = args.worktree.resolve()

    baseline = run(["cargo", "test", "--release"], worktree)
    if baseline.returncode:
        print("baseline suite does not pass; aborting", file=sys.stderr)
        print(baseline.stdout[-4000:], file=sys.stderr)
        sys.exit(1)
    print("baseline: all tests pass", flush=True)

    results = []
    for name, relative, old, new, description in MUTATIONS:
        path = worktree / relative
        original = path.read_text()
        if old not in original:
            results.append(
                {"mutation": name, "file": relative, "status": "not_applicable",
                 "description": description}
            )
            print(f"{name}: FRAGMENT NOT FOUND in {relative}", flush=True)
            continue
        mutated = original.replace(old, new, 1)
        if name == "pep_prior_doubled":
            mutated += HELPER
        path.write_text(mutated)
        try:
            outcome = run(["cargo", "test", "--release"], worktree)
            combined = outcome.stdout + outcome.stderr
            compiled = "error[E" not in combined and "could not compile" not in combined
            failures = failing_tests(combined)
            caught = bool(failures) or (not compiled)
            results.append(
                {
                    "mutation": name,
                    "file": relative,
                    "description": description,
                    "compiled": compiled,
                    "failing_tests": failures,
                    "caught": caught,
                    "status": "caught" if caught else "NOT CAUGHT",
                }
            )
            print(
                f"{name}: {'caught by ' + str(len(failures)) + ' test(s)' if failures else ('rejected at compile time' if not compiled else 'NOT CAUGHT')}"
                + (f" -> {failures[:6]}" if failures else ""),
                flush=True,
            )
        finally:
            path.write_text(original)

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps({"mutations": results}, indent=1) + "\n")
    missed = [r for r in results if r.get("status") not in ("caught",)]
    print(f"\n{len(results) - len(missed)}/{len(results)} mutations caught")
    if missed:
        print("NOT CAUGHT: " + ", ".join(r["mutation"] for r in missed))


if __name__ == "__main__":
    main()
