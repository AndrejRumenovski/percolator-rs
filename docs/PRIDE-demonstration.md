# Real PRIDE demonstration — 2026-09-06

This demonstration used public **PXD032157**, “A male-specific steroid hormone controls female reproductive biology in the malaria mosquito,” through the official [PRIDE Archive v3 API](https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD032157). The complete experiment's remote IDs, expected/actual digests, executable identity, parameters and result hashes are in [the compact execution record](pride-demo-summary.json). No proteomics source files were added to this repository.

A dedicated cache at `/tmp/percolator-pride-demo-20260906` used a **100 MB test ceiling**, overriding the normal **50 GB maximum**. No space was preallocated. Safety margin remained the default 1 GB of filesystem headroom.

## Metadata first

```sh
percolator-rs pride --cache-dir /tmp/percolator-pride-demo-20260906 \
  --cache-limit 100MB info PXD032157
```

| Inspection | Observed |
| --- | ---: |
| Public project status | PUBLIC |
| Indexed file records | 651 |
| Supplemental checksum-table records | 2 |
| Combined reported project size | 275,424,583,570 bytes |
| RAW data | 189,712,154,388 bytes |
| Processed data under the documented category/format filter | 3,185,516,792 bytes |
| Uncompressed PIN candidates | 65 |
| Source bytes downloaded by metadata inspection | 0 |

Organism: *Anopheles gambiae* str. pest. Instrument: Q Exactive HF. Experiment type: shotgun proteomics. Reported tissues include atrium, male reproductive gland and whole body. The entire project exceeds both the demonstration ceiling and the default 50 GB ceiling. One independent input was selected explicitly.

## Plan, download and verify

Selected file:

```text
28May2015-QE-HF-Anopheles-23-atrium-P-24H-2nd-01-comet.pin
8,254,268 bytes
```

```sh
percolator-rs pride --cache-dir /tmp/percolator-pride-demo-20260906 \
  --cache-limit 100MB fetch PXD032157 \
  --file 28May2015-QE-HF-Anopheles-23-atrium-P-24H-2nd-01-comet.pin \
  --max-download 16MB --max-working-space 2GB --dry-run
# The same selection was executed with --yes instead of --dry-run.
```

Both supplied authoritative checksums matched the streamed bytes:

| Evidence | Expected = actual |
| --- | --- |
| PRIDE file-record SHA-1 | `256a45345dee1766fc63092fc956b83934e9d91c` |
| PRIDE checksum-table MD5 | `467be2d99f77ba5d4a94ef6cd0df1d7d` |
| Locally calculated SHA-256 (not repository-published) | `c45cfe99282b5a2dbab5ad5e153302d2f8b6540a68eae646417ae6ecc111fe6b` |

No `.part` file was published as verified. Cache source usage was exactly **8,254,268 bytes** after completion.

## Pin protection

The source was pinned, `cache prune --all-evictable` was executed, and then it was unpinned:

```text
Storage before:       8,254,268 bytes
Objects deleted:      0
Space freed:          0 bytes
Pinned data remaining:8,254,268 bytes
```

## Existing analysis, followed by source eviction

```sh
percolator-rs pride --cache-dir /tmp/percolator-pride-demo-20260906 \
  --cache-limit 100MB run PXD032157 \
  --file 28May2015-QE-HF-Anopheles-23-atrium-P-24H-2nd-01-comet.pin \
  --ephemeral --max-download 0 --max-working-space 2GB --yes \
  -- --profile fast --seed 1
```

The cached input was rehashed and reused, so the analysis required **zero additional download bytes**. The existing PIN reader accepted **31,660 PSMs and 21 features** (16,527 targets, 15,133 decoys). The ordinary fast profile ran with its existing SVM behavior and parameters. Its reported target PSM count at q < 0.01 was 351; this is an observed output of that profile, not a new validation of experimental database design.

All four target/decoy PSM and peptide tables passed structural/numeric validation and received SHA-256 hashes before source eviction. The manifest records an experiment state of `verified`, its CLI parameters, the build's parent commit, and the exact executable hash. The source record retains both repository checksums, local SHA-256 and successful PIN validation, with availability changed to `evicted`.

## After cleanup

Both a dry preview and actual `cache prune --all-evictable` were run after the ephemeral experiment. There was already no large data left to delete:

| Storage class | Bytes remaining |
| --- | ---: |
| Downloaded source objects | **0** |
| Prepared artifacts | **0** |
| Partial/temporary files | **0** |
| Pinned large data | **0** |
| Final result tables, retained | **1,559,909** |

The manifest and processing provenance remained available locally. All four retained result hashes were independently recomputed after cleanup and matched. The target PSM result SHA-256 is:

```text
9f2f209332ced24332cea546b038b1ed82aeac6e51e0889d6ef9eb04a7c3ae78
```

The machine-readable execution record is retained in the repository; the original full manifest and four result files are in the dedicated external demonstration cache. Source/intermediate usage reached **zero**, independently of the configured maximum.

## Offline regression evidence

The normal test suite additionally processes two independent PIN fixtures whose combined source size exceeds its deliberately small cache ceiling, releasing each batch after verification. It tests nonzero-to-zero manual pruning, KEEP overrides with `purge-data --yes`, pinned-object preservation, interrupted transfers, invalid ranges, checksum mismatches, failed output-size limits, and external preparation lineage. The PRIDE path's result bytes are compared against the ordinary CLI for the same fixture and options.
