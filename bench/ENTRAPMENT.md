# Signal-present FDR calibration by proteome entrapment

## Result

On six PXD032157 mzML runs searched against an amino-acid-balanced native + plant-entrapment
database, both implementations' nominal q-values are anti-conservative. At reported q ≤ 0.01,
percolator-rs accepts 19,666 target PSMs with an entrapment-estimated FDP of **2.78%**; C++
Percolator 3.09 accepts 19,126 with an estimated FDP of **2.62%**. The conditional 95% Wilson
intervals are 2.53–3.06% and 2.36–2.89%, respectively. These intervals condition on the empirical
entrapment search-space fraction and do not include uncertainty in that fraction or dependence
between repeated PSMs, so they should be read as descriptive rather than definitive inferential
bounds.

| reported q cutoff | C++ accepted / entrapment / FDP | Rust accepted / entrapment / FDP |
|---:|---:|---:|
| 0.001 | 13,076 / 143 / 1.41% | 14,202 / 142 / 1.18% |
| 0.005 | 17,600 / 277 / 2.00% | 17,999 / 280 / 1.92% |
| 0.010 | 19,126 / 369 / 2.62% | 19,666 / 408 / 2.78% |
| 0.020 | 20,619 / 537 / 3.53% | 21,129 / 587 / 3.88% |
| 0.050 | 23,255 / 1,132 / 6.66% | 23,821 / 1,268 / 7.29% |
| 0.100 | 26,188 / 2,237 / 11.54% | 26,911 / 2,472 / 12.27% |

The nominal-q yield lead remains (+540 PSMs, +2.82% at q ≤ 0.01), but it cannot be described as
validated at an actual 1% FDR. The extra yield is accompanied by about 47 more empirically estimated
false PSMs after method-specific search-space correction. This experiment therefore closes the
previous uncertainty with a negative result: the pure-null control is conservative, while
signal-present target-decoy q-values are not exactly calibrated for either implementation on this
search. Rust is slightly more anti-conservative at the 1% cutoff, although the absolute difference
between implementations is much smaller than their shared departure from nominal.

## Design

- Spectra, native FASTA, and per-run Comet parameter snapshots are the deposited PXD032157 files
  from [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD032157). All mzML and native-FASTA
  SHA-1 values are checked against the archive manifest.
- Six deposited mzML files were used (about 6.1 GiB and 168,000 searchable spectra in total).
- The native database has 139,191 proteins / 145,351,640 residues. Nine plant reference proteomes
  from UniProt release 2026_02 are prefixed `ENT_` and truncated at a protein boundary to
  145,351,799 residues. The nominal entrapment fraction is therefore 0.50000027.
- Comet is run through Crux 4.0, which embeds Comet 2019.01 rev. 5—the version recorded by the
  deposited workflow—with the deposited per-run parameters and automatic decoy generation.
- Native/plant shared-protein assignments are excluded. Only PSMs assigned entirely to `ENT_`
  proteins count as known errors. Among high-scoring decoys, plant matches make up about 72–85%,
  rather than 50%, because equal residue counts do not imply equal effective search space under
  semi-tryptic digestion and database redundancy. Each table row therefore uses the observed
  plant fraction among non-mixed decoys passing the same method and q cutoff.
- The plant FASTA hashes in `bench/entrapment/proteomes.sha256` pin the exact UniProt inputs. Two of
  194 distinct Rust entrapment peptide sequences at q ≤ 0.01 occur in the native FASTA (including
  I/L equivalence); excluding those two only makes a negligible change.

## Reproduce

```bash
bash bench/entrapment/run.sh
```

The script downloads about 6.3 GB to `$HOME/percolator_rs_out/entrapment` by default, constructs
the balanced database, searches all six mzML files, runs both rescorers on every identical PIN, and
writes `report.tsv`. Override the workspace with `ENTRAPMENT_WORK`, C++ Percolator with `PERC`, or
Crux with `CRUX`.

Measured on the benchmark machine, the six Comet searches took 4,243.85 s (70.7 min) total and
peaked at 3.42 GiB RSS. Rescoring took 6.30 s total in Rust versus 105.06 s in C++ (16.7× faster).
The search dominates the experiment and is not included in the implementation speed comparison.
