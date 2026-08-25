# Multi-dataset benchmark

> **These measurements predate the 2026-08-25 statistical repair** and describe commit `d83a7ba`,
> whose q-value and PEP estimators, cross-validation isolation and PIN feature selection were all
> subsequently found defective and replaced. They are kept as the record of what was measured then.
> For what the current implementation does and what it has been revalidated against, see
> [`../validation/REPAIR.md`](../validation/REPAIR.md).

The speed advantage generalizes across all five tested search configurations; the identification
delta does not. On the four compact extension cases, percolator-rs is **7.3–14.6x faster** and uses
**37–67% of the C++ peak RSS**, while its PSM count at reported q≤0.01 ranges from **−1.8% to
+12.0%** relative to C++ Percolator 3.09. PXD032157 remains the large-scale case at 18x faster and
+3.9% PSMs. This is stronger evidence than a single dataset, but it is not evidence that either
implementation's reported q-values are calibrated to a true 1% FDR.

## Dataset matrix

| Case | Organism | Instrument | Search / PIN source | Database scale | Rescoring input |
|---|---|---|---|---:|---:|
| [PXD032157](https://www.ebi.ac.uk/pride/archive/projects/PXD032157) | *Anopheles gambiae* | Q Exactive HF | Comet | 139,191 target proteins | 65 PINs; 8,639,746 PSMs; 2.30 GB |
| [PXD007145](https://www.ebi.ac.uk/pride/archive/projects/PXD007145) / Hogrebe | Human phosphoproteome | Q Exactive HF / Fusion Lumos study | Tide re-search distributed by mokapot | 20,416 target Swiss-Prot proteins | 55,398 PSMs; 21 features |
| [PXD020243](https://www.ebi.ac.uk/pride/archive/projects/PXD020243) | Human ALT cells | Q Exactive HF | MSFragger 20180316 pepXML | Human Swiss-Prot + spike-ins; 5,153 protein IDs represented in this run | 9,475 PSMs; 19 features |
| [PXD060954](https://www.ebi.ac.uk/pride/archive/projects/PXD060954) | *A. muciniphila* / *B. breve* search DB | timsTOF Pro | Sage 0.14.7 release asset (binary reports 0.14.6) | 4,647 targets + 9,294 deposited decoys | 35,093 PSMs; 28 features |
| [Percolator upstream fixture](https://github.com/percolator/percolator/blob/8c7412ea0e556dddbc2dfa26c0641f49c948378e/data/percolator/tab/percolatorTab) | *S. cerevisiae* | Not retained in fixture metadata | Legacy SQT/SEQUEST-style PIN | 15,368 protein IDs represented | 19,674 parsed PSMs; 19 features |

“Database scale” is the searched target FASTA size where that FASTA is public. For the MSFragger
and legacy fixtures the original FASTA is not deposited with the input, so the table reports the
number of distinct protein identifiers actually represented rather than inventing an exact search
database size. The PXD007145 PIN is a Tide reanalysis, not the study's original MaxQuant search.

## Results

Same machine (Ryzen 5 5600G), local ext4 input/output, seed 1, full ten-iteration/default training,
and identical normalized PIN input for both implementations. The archived production searches use
concatenated target+decoy databases. C++ is explicitly told `--search-input concatenated` on all
compact PINs to match percolator-rs's direct target-decoy post-processing; its count-based
auto-detector otherwise misclassifies the target-heavy MSFragger and Sage inputs as separate and
switches to mix-max. The legacy fixture no longer retains its original search mode. Counts use
q≤0.01 by reading the named `q-value` column, including C++ outputs that add a `filename` column.

| Case | Wall, Rust / C++ | Speedup | Peak RSS, Rust / C++ | PSMs, Rust / C++ | Rust PSM delta | Peptides, Rust / C++ |
|---|---:|---:|---:|---:|---:|---:|
| PXD032157 (N=4 processes) | 20.8 / 376.2 s | **18.1x** | 0.87 / 1.49 GiB | 107,046 / 103,038 | **+3.89%** | 37,469 / 35,852 |
| PXD007145 Tide | 0.53 / 5.99 s | **11.3x** | 53.1 / 131.4 MiB | 29,264 / 27,617 | **+5.96%** | 20,614 / 19,722 |
| PXD020243 MSFragger | 0.07 / 0.51 s | **7.3x** | 12.3 / 18.2 MiB | 1,554 / 1,388 | **+11.96%** | 1,177 / 1,062 |
| PXD060954 Sage | 0.31 / 4.54 s | **14.6x** | 32.0 / 86.4 MiB | 26,624 / 25,795 | **+3.21%** | 11,420 / 11,336 |
| Upstream yeast fixture | 0.10 / 0.90 s | **9.0x** | 22.0 / 33.5 MiB | 1,126 / 1,147 | **−1.83%** | 903 / 928 |

Compact wall/RSS figures are medians of three runs and include process startup, so they should be
read as order-of-magnitude portability checks, not sub-100-ms microbenchmarks. The full PXD032157
result is the throughput measurement. The key generalization result is that every new schema runs
and the speed/RSS direction is consistent; yield moves both ways and is dataset-dependent.

## Reproduction and input handling

Run:

```bash
bash bench/multidataset/run.sh
```

The harness downloads about 180 MB, verifies every source against
[`sources.sha256`](multidataset/sources.sha256), re-searches the deposited PXD060954 MGF with the
checked-in Sage configuration, converts the deposited MSFragger pepXML without third-party Python
packages, verifies the generated PINs against [`generated.sha256`](multidataset/generated.sha256),
and runs both rescorers. The table's machine-readable compact results are preserved in
[`recorded-results.tsv`](multidataset/recorded-results.tsv). All bulk data and outputs go to
`$HOME/percolator_rs_out/multidataset` (override with `MULTIDATASET_OUT`). The existing 2.3-GB
PXD032157 benchmark is not downloaded or rerun by this compact harness; `bench/regression.sh` is its
separate reproducible gate.

`REPEATS=3` controls the compact timing repetitions; three is the default and the reported table
uses their median.

Two deterministic adaptations are required:

- [`pepxml_to_pin.py`](multidataset/pepxml_to_pin.py) retains every PXD020243 MSFragger search hit,
  search score, protein mapping, mass/charge feature, and target/decoy label.
- Sage writes a full MGF title into `ScanNr`. percolator-rs accepts it, but C++ 3.09 aborts because
  it requires an integer. [`normalize_sage_pin.py`](multidataset/normalize_sage_pin.py) extracts the
  title's deposited `#scan` integer. It also removes Sage's already-trained `posterior_error` to
  avoid circular rescoring and three constant mobility fields absent from the MGF export. Because
  Sage's parallel writer can vary row order and its `SpecId` is order-derived, the normalizer sorts
  complete records and assigns stable sequential IDs. Both implementations then receive that exact
  same normalized PIN.

## Interpretation and limits

- These are reported-q yield comparisons, not ground-truth accuracy comparisons. The
  [signal-present entrapment experiment](ENTRAPMENT.md) shows that both implementations are
  anti-conservative on the tested PXD032157 search.
- The +12.0% MSFragger result is the largest Rust/C++ divergence and should motivate targeted
  calibration work; it must not be advertised as an accuracy gain without entrapment or a known
  mixture for that schema.
- The yeast case usefully falsifies a universal “Rust always yields more” claim: C++ finds 1.8%
  more PSMs and 2.7% more peptides there.
- The current matrix is DDA-focused. A future extension should add a valid DIA-derived PIN. Protein
  inference calibration outside PXD032157 is handled separately by the
  [PrEST standard benchmark](PROTEIN_CALIBRATION.md).
