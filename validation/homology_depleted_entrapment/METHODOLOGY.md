# Methodology and exact reproduction record

## Scope and invariants

This experiment intervenes on the entrapment target database before spectrum
searching. It does not modify `percolator-rs`, recalibrate PEPs, introduce an
empirical correction, tune the homology threshold after seeing results, or
filter an existing PSM set. C++ Percolator is not used as an oracle. Every
database condition is searched independently from the six raw mzML files, and
Comet regenerates decoys for each search.

The frozen repository commit is
`e8d83d1c76e4cf651fdfcf22d98b0b499c35943a`. The audited production binary is
`target/release/percolator-rs`, SHA-256
`be9bf670bfd69df4dc3ba3b8be6c4c164acaf56a91f4a2819d115f49604b2c45`.
The experiment used that binary without rebuilding it. A final source-integrity
check confirms no diff under `src/`, `Cargo.toml`, `Cargo.lock`, or `build.rs`.

## Temporal order and baseline gate

The following order was enforced:

1. Re-run the canonical seed-1 rescoring against the six previously frozen
   original PINs and reproduce the canonical endpoint.
2. Freeze input, software, parameter, database, spectrum, PIN, and result
   hashes in `baseline_reproduction/baseline_manifest.json`.
3. Write and hash `PREREGISTRATION.md` before constructing or searching an
   intervention database.
4. Construct and characterize all target databases without reading search or
   Percolator results.
5. Search all conditions from raw mzML and analyze raw XCorr.
6. Run the same audited `percolator-rs` binary for all confirmatory and
   diagnostic conditions.
7. Run the preregistered analysis, followed by the explicitly exploratory
   residual-tail characterization.

The baseline gate passed exactly for the discrete endpoints:

| metric | frozen canonical value | reproduced value |
|---|---:|---:|
| accepted targets, q<0.01 | 19,545 | 19,545 |
| entrapment targets, R_ent | 258 | 258 |
| entrapment decoys, D_ent | 133 | 133 |
| R_ent/D_ent | 1.939850 | 1.939850 |
| adjusted FDP | 0.01687316 | 0.01687316 |
| direct known-false lower bound | 0.01320031 | 0.01320031 |
| pooled adjusted PEP error | +0.01850544 | +0.01850544 |
| raw XCorr R_ent/D_ent at D_ent=133 | 134/133 | 134/133 |

The complete calibration bins are in `baseline_reproduction/summary.json`.
The fresh original search produced the same primary counts and ratio. Across
all six new original PINs, PSM identities and substantive feature values were
unchanged; 394 rows differed only in `lnrSp` tie ranking or last-decimal
`CalcMass`. The resulting fresh-search pooled PEP error was +0.01853124. These
minor search-engine tie differences are preserved in `search_manifest.json`
and are not corrected.

## Software, spectra, and search parameters

- Crux: 4.0-b854e9b-2021-06-01, SHA-256
  `bad553f44f8fadbe8bc276c42caeb292ff58dcb8cbbac6215c852501d0adf687`.
- Embedded Comet: 2019.01 rev. 5.
- ProteoWizard reported by Crux: 3.0.21153.
- Python: 3.14.4; NumPy: 2.5.1.
- Platform: Linux 7.0.0-30-generic, x86-64, glibc 2.43.
- Spectra: six PXD032157 mzML runs; paths, byte sizes, and SHA-1 hashes are in
  the baseline manifest.

The six deposited per-run Comet parameter files are copied and hashed under
`baseline_reproduction/search_parameters/`. The relevant resolved settings
are trypsin, semi-tryptic searching (`num_enzyme_termini=1`), two missed
cleavages, MH+ 600--5000 Da, peptide length 1--63, 10 ppm precursor tolerance,
static Cys +57.021464, five reported candidates per spectrum, and a
concatenated target/decoy search (`decoy_search=1`). Comet reverses each target
peptide while retaining its enzyme-terminal residue. Decoy generation is
deterministic and has no random seed.

## Frozen homology intervention

The primary comparison universe is the theoretical fully tryptic digest:
cleavage after K/R except before P, zero to two missed cleavages, peptide length
8--63, MH+ 600--5000 Da with static Cys +57.021464, and no ambiguous residues.
I and L are canonicalized to L.

An entrapment peptide is prohibited when a native peptide of the same length
has both I/L-aware Hamming distance at most two and I/L-aware identity at least
85%. This permits at most one difference for lengths 8--13 and at most two for
length 14 or greater. This one rule was chosen from the prior same-length,
I/L-aware distance-at-most-two evidence and frozen before intervention search
results existed.

The depletion unit is the whole entrapment protein. If a protein generates any
prohibited peptide, the complete protein is removed. This preserves coherent
FASTA records and avoids artificial termini. It also makes the causal estimand
the effect of removing *proteins capable of generating prohibited peptides*,
not the isolated effect of deleting only the witness peptides.

Three negative-control databases use NumPy PCG64 seeds 130363, 155921, and
196613. Within each source proteome and fixed 50-residue protein-length stratum
(with a final >=2000 stratum), each control removes exactly the homology arm's
protein count without using homology status. The confirmatory size comparison
is the mean of these three controls; every control is also reported separately.

## Pre-search characterization

The original entrapment component contains 389,504 proteins. Homology depletion
removes 200,141 (51.38%) and leaves 189,363; each random control leaves the same
number with exact source-by-length-stratum removal counts.

| condition | proteins | searchable peptide instances | unique searchable, HLL |
|---|---:|---:|---:|
| original | 389,504 | 598,438,171 | 418,094,318 |
| homology depleted | 189,363 | 201,087,250 | 156,459,902 |
| control 130363 | 189,363 | 202,659,292 | 165,829,534 |
| control 155921 | 189,363 | 202,581,927 | 165,724,726 |
| control 196613 | 189,363 | 202,578,864 | 165,554,557 |

The intervention removes 65.69% of entrapment residues, 66.40% of searchable
semi-tryptic instances, and an estimated 62.58% of unique searchable sequences.
It has 0.75% fewer searchable instances and 5.59% fewer HLL-estimated unique
sequences than the control mean. This fixed residual mismatch is reported as a
limitation and is not corrected post hoc. HLL uses p=18 (nominal relative SE
0.203%); instance counts are exact.

The original database has 200,141 qualifying entrapment proteins and 7,857
exact native-shared full-tryptic sequences. The homology arm has zero of both;
the controls retain 81,183--81,391 qualifying proteins and 4,668--4,698 exact
shared sequences. No condition has duplicate accessions or I/L-canonical
native/entrapment full-protein collisions. Full peptide-length, precursor-mass,
enzymatic-terminus, amino-acid, and similarity distributions are in
`databases/presearch_summary.json` and
`databases/presearch_characterization.json`.

The rule proved broad after it was frozen: 168,196 of 200,141 first witnesses
have length 8, and 189,792 are one-substitution witnesses. The rule was not
changed after observing this.

## Search, labels, and rescoring

Five complete target FASTAs were built: original, homology depleted, and the
three controls. Their hashes are in `databases/construction_manifest.json` and
`search_manifest.json`. Thirty independent Comet searches were run: five
conditions by six raw mzML files. Every search emitted a new PIN, pepXML,
resolved parameter file, log, and regenerated decoys.

PIN integrity checks require every `Label=1` row to map to at least one target
protein and every `Label=-1` row to map only to decoy proteins. There are 14,300
target rows across all searches whose peptide also maps to a decoy protein;
they remain valid target rows and are marked mixed for entrapment adjustment.
They are never relabeled or filtered. Pure-entrapment counts require all target
protein mappings, after stripping `DECOY_`, to have the `ENT_` prefix.

All primary runs use:

```text
percolator-rs --canonical --no-select-c --seed 1 --num-threads 1 \
  --results-psms TARGET.tsv --decoy-results-psms DECOY.tsv PIN
```

The program default `maxiter=10` is canonical. Diagnostic runs use maxiter 0,
1, 2, 3, and 10; seeds 1, 2, and 3 for original and homology-depleted; and PIN
copies lacking exactly `enzN` and `enzC` for the two specified ablation arms.

## Endpoint definitions

The primary endpoint is pooled `R_ent/D_ent` among seed-1 PSMs with reported
`q<0.01`, strict inequality. The direct known-false lower bound is
`R_ent / accepted targets`. The adjusted FDP divides R_ent by the condition's
global pure-entrapment fraction among usable non-mixed decoys, then by accepted
targets. Each database therefore uses its own measured entrapment opportunity;
the fractions are 0.782323 original, 0.576893 homology depleted, and
approximately 0.5956 in the controls.

PEP calibration uses the ten frozen bins in `PREREGISTRATION.md`; nine are
populated. For each bin the observed adjusted false fraction is compared with
the mean predicted PEP. Pooled signed error is total entrapment-implied false
mass minus total predicted false mass, divided by all target PSMs. Bin-level
shared-cause analysis compares
`log(observed adjusted / predicted)` with `log(R_ent/D_ent)` and reports Pearson
correlation and median absolute log-fold residual. No bins are merged.

Raw XCorr first applies the same spectrum-key competition as the prior audit,
then records R_ent when D_ent first reaches 25, 50, 100, 133, 250, and 500.
This analysis was frozen before inspecting rescored intervention endpoints.

## Uncertainty

- Percolator seeds 1--3 measure algorithmic reproducibility only.
- Three database-depletion seeds show the behavior of the frozen random control
  construction; three seeds are not a population-level confidence interval.
- Wilson intervals in the calibration table describe PSM-level binomial
  counting for the direct known-false fraction only.
- A six-run paired cluster bootstrap, seed 20260831 and 10,000 replicates,
  retains all PSM multiplicity within each LC-MS/MS run. It is the primary
  sampling-uncertainty sensitivity analysis, but only six clusters are
  available and repeated peptides across runs remain correlated.
- Only one dataset family was tested. Dataset-level and organism-level
  uncertainty are unquantified and limit generalization.

## Exact command sequence

All search and rescoring commands are implemented literally in
`run_pipeline.sh`. From the repository root, the end-to-end sequence is:

```bash
EXP="$PWD/validation/homology_depleted_entrapment"
SOURCE="/home/andrej-rumenovski/percolator_rs_out/entrapment"

# Baseline rescoring was written below baseline_reproduction/results/seed-1
# with the canonical command above, then summarized and gated:
python3 validation/pep_rootcause_experiments.py summarize \
  --root "$EXP/baseline_reproduction/results/seed-1" \
  --output "$EXP/baseline_reproduction/summary.json"
python3 "$EXP/freeze_baseline.py" --repo "$PWD" --source-root "$SOURCE" \
  --results-root "$EXP/baseline_reproduction/results/seed-1" \
  --output-root "$EXP/baseline_reproduction"

# PREREGISTRATION.md and its SHA-256 were frozen here.
python3 "$EXP/construct_databases.py" --combined "$SOURCE/combined.fasta" \
  --output-root "$EXP/databases" 2>"$EXP/databases.construction.log"
g++ -O3 -std=c++17 -fopenmp "$EXP/characterize_databases.cpp" \
  -o "$EXP/characterize_databases"
"$EXP/characterize_databases" "$EXP/databases/original.fasta" \
  "$EXP/databases/homology_depleted.removed.tsv" \
  "$EXP/databases/size_control_130363.removed.tsv" \
  "$EXP/databases/size_control_155921.removed.tsv" \
  "$EXP/databases/size_control_196613.removed.tsv" \
  "$EXP/databases/presearch_characterization.json" \
  2>"$EXP/databases/presearch_characterization.log"
python3 "$EXP/presearch_report.py" --database-root "$EXP/databases" \
  --output "$EXP/databases/presearch_summary.json"

bash "$EXP/run_pipeline.sh" search
python3 "$EXP/freeze_searches.py" --experiment-root "$EXP" \
  --canonical-root "$SOURCE" --output "$EXP/search_manifest.json"
python3 "$EXP/analyze.py" --experiment-root "$EXP" raw

bash "$EXP/run_pipeline.sh" primary
bash "$EXP/run_pipeline.sh" dose
bash "$EXP/run_pipeline.sh" seeds
bash "$EXP/run_pipeline.sh" enz-ablation
python3 "$EXP/analyze.py" --experiment-root "$EXP" rescored

python3 "$EXP/exploratory_remaining.py" --experiment-root "$EXP" \
  --native-fasta "$SOURCE/native.fasta" \
  --output "$EXP/analysis/exploratory_remaining.json"
python3 "$EXP/render_tables.py" --experiment-root "$EXP"
python3 "$EXP/plot_results.py" --experiment-root "$EXP"
python3 "$EXP/freeze_artifacts.py" --experiment-root "$EXP" \
  --repo "$PWD" --binary "$PWD/target/release/percolator-rs"
```

The original baseline PIN hashes and fresh PIN hashes are recorded separately.
The full file-by-file integrity inventory is `ARTIFACT_MANIFEST.json` with a
compact companion `SHA256SUMS.txt`.
