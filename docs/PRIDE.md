# PRIDE Archive: a disposable working cache

PRIDE is the durable public source. `percolator-rs pride` retains metadata and reproducibility records while treating downloaded source data and prepared PINs as replaceable. It does not change rescoring, statistics, cross-validation, competition, SVMs, or inference.

## Quick start

```sh
# Default ceiling: 50,000,000,000 bytes, with no preallocation or target occupancy.
# Put large data on a filesystem of your choice, outside any Git working tree.
export PERCOLATOR_RS_PRIDE_CACHE=/mnt/scratch/percolator-pride
export PERCOLATOR_RS_PRIDE_CACHE_LIMIT=50GB

percolator-rs pride info PXD032157
percolator-rs pride files PXD032157 --format PIN
percolator-rs pride manifest PXD032157 > experiment-manifest.json

# Exact inventory ID or filename; commands below use one real 8.3 MB input.
percolator-rs pride fetch PXD032157 \
  --file 28May2015-QE-HF-Anopheles-23-atrium-P-24H-2nd-01-comet.pin \
  --max-download 16MB --dry-run

# Add --yes to execute the printed plan.
percolator-rs pride run PXD032157 \
  --file 28May2015-QE-HF-Anopheles-23-atrium-P-24H-2nd-01-comet.pin \
  --max-download 16MB --max-working-space 2GB --ephemeral --yes \
  -- --profile fast --seed 1

percolator-rs pride cache status
percolator-rs pride cache prune --all-evictable --dry-run
percolator-rs pride cache prune --all-evictable
```

`fetch`, `run`, and `prepare` print a plan and require `--yes` for execution. The default transfer budget is **1 GB per operation**, independent of the **50 GB cache ceiling**. Retry bytes consume this budget too. Selecting a whole dataset requires `--all`; an empty selection never silently selects everything. `--file` is repeatable. Ambiguous filenames require exact inventory IDs. Format/category alternatives are ORed within each option; different option groups intersect.

`--cache-dir`, `--cache-limit`, and `--dry-run` are global options. Size arguments accept integer bytes and decimal `KB/MB/GB/TB` or binary `KiB/MiB/GiB/TiB`. A CLI ceiling applies to that invocation; use the environment variable consistently for a persistent preference. Without a location override, use `$XDG_CACHE_HOME/percolator-rs/pride` or `$HOME/.cache/percolator-rs/pride`.

## Official interface inspected on 2026-09-06

Implementation target: **PRIDE Archive REST API 3.0**, base `https://www.ebi.ac.uk/pride/ws/archive/v3`, whose live [OpenAPI document](https://www.ebi.ac.uk/pride/ws/archive/v3/v3/api-docs) advertises that exact server URL. The doubled `v3` in the OpenAPI URL is intentional. The older v2 prose guide still exists, but its HAL examples and older MSRun endpoints are not the implementation contract.

| Interface | Use |
| --- | --- |
| `GET /status/{accession}` | Public/private accessibility; only `PUBLIC` projects are accepted |
| `GET /projects/{projectAccession}` | Project details and controlled vocabulary metadata |
| `GET /projects/{projectAccession}/files/count` | Inventory count before and after pagination |
| `GET /projects/{projectAccession}/files?pageSize=100&page=N` | Zero-based file pagination; follow count even if server returns fewer than requested |
| `GET /projects/{projectAccession}/files/all` | Official unpaginated alternative, inspected but not used by the client |
| `GET /files/{fileAccession}` | Individual file metadata; checksum provenance references point here |
| `GET /files/checksum/{projectAccession}` | Separately labelled legacy MD5/size TSV |
| `GET /projects/files-path/{projectAccession}` | Official project download root for supplemental checksum entries |
| `GET /search/projects?keyword=…&filter=…&page=N&pageSize=N` | Public discovery, keyword/filter expression and explicit page selection |

Project fields used: accession, title, description, organisms, organism parts (tissues), instruments, experiment types, modifications, references/publication IDs, submission/publication dates, DOI, submission type, sample-processing protocol and data-processing protocol. Search returns a different representation for some terms; this is normalized inside the client. Absent fields serialize as `null`; empty remote arrays remain empty. `submissionType` (for example COMPLETE/PARTIAL) is not public/private status.

File fields used: accession, filename, category, format when supplied, `fileSizeBytes`, `publicFileLocations`, checksum, analysis accessions and additional attributes. The current OpenAPI does not advertise the old dedicated MSRun endpoints. Available analysis links/additional metadata are retained; missing run details are not invented. Format inferred from a filename is distinct from a remotely supplied format field.

The client checks counts before and after paging and rejects repeated IDs, empty premature pages, and changed counts. It merges unambiguous MD5-table records with the indexed inventory and adds supplemental entries using the returned project FTP root. This matters in practice: the demonstration project has **651 indexed records and two additional checksum-table records**. Archive directories can contain additional files not represented in either API source; the manifest explicitly records this inventory limitation. The implementation does not claim to crawl every FTP directory.

### Downloads and checksum evidence

[PRIDE's official download guide](https://github.com/orgs/PRIDE-Archive/discussions/33) documents HTTPS on `ftp.pride.ebi.ac.uk` with the same path as its FTP reference. The client makes that specific scheme conversion; it never guesses a year/month path. HTTPS references work directly. Aspera and other FTP hosts are retained as references but not automatically invoked.

The file-record `checksum` field and the MD5 table are **different evidence**. PRIDE's current [checksum generator](https://github.com/PRIDE-Archive/pride-checksum) documents SHA-1. The client interprets a 40-character hexadecimal file-record checksum as SHA-1; other untyped encodings remain opaque and block downloading rather than silently bypassing verification. The TSV's explicitly labelled MD5 values are separately checked. Expected and actual digests, algorithms, source endpoints and outcomes are recorded. All supplied supported checksums must match. Conflicting inventory/table sizes stop selection before transfer. The real demonstration verifies both SHA-1 and MD5 against the downloaded PIN.

A separate local SHA-256 is always computed by streaming. It is never described as a repository-published checksum. If PRIDE supplies no checksum, fetch records `downloaded_unverified`; analysis requires explicit `--allow-unverified`. This records weaker reproducibility evidence, rather than claiming repository verification.

### Search

```sh
percolator-rs pride search 'mosquito proteome' --page 0 --page-size 20
percolator-rs pride search 'cancer' --filter 'organisms==Homo sapiens' --page 1
```

The supported API parameters are `keyword`, `filter`, `page`, and `pageSize` (CLI range 1–100). The filter is the official `field==value` expression passed unchanged and URL-encoded. The CLI does not invent organism/instrument/title/date switches or promise that every possible field name is indexed. Consult the deployed API for available filter fields. Search displays one requested page and does not download project data. Metadata inspection reuses saved manifests, including offline after eviction; `--refresh` updates remote metadata and preserves local history and changed remote identities.

## Native input and scientific scope

A file extension alone does not establish valid input. Uncompressed `.pin` files are native **candidates**. They become directly compatible only after the existing PIN parser accepts the downloaded content. Saved validation is tied to the object identity and survives eviction. mzML/RAW require database search; mzIdentML/mzTab/search results require appropriate external preparation, with complete target/decoy candidates and numeric features. Compressed PINs require external decompression. No archives are extracted automatically.

`run` invokes the same `percolator-rs` executable with its existing scientific CLI. Options after `--` are validated against a bounded list of existing scientific flags; outputs and input paths are managed by PRIDE. The default invocation inherits the ordinary CLI defaults. The existing concatenated target/decoy input contract applies; neither a filename nor a successful parser can establish correct experimental database design.

Multiple inputs require `--independent-runs`, explicitly selecting **separate models and statistics**. `--batch-size N` controls how many inputs are held before processing and releasing them. Execution is serial; CPU/network concurrency does not bypass storage planning. This permits aggregate dataset size to exceed the cache ceiling in ephemeral mode, provided each batch plus temporary output allowance fits. Pinned/KEEP inputs are budgeted as retained across batches. For pooled training or engine ensembles, use the ordinary CLI with simultaneous required inputs; PRIDE never silently substitutes independent runs for those workflows.

Four result tables are written: target/decoy PSMs and target/decoy peptides. Automatic PRIDE runs do not currently request protein reports or feature reports; those remain available through the ordinary CLI. Parameters and final hashes are retained. Native PRIDE analysis currently requires Unix file-size limits (`RLIMIT_FSIZE`); Linux additionally kills the analysis child if its parent exits. Download/metadata/cache operations do not use that analysis-specific mechanism.

## External preparation and lineage

There is no new RAW converter or search-engine runner. Prepare data with an external tool and import its PIN with an explicit recipe:

```sh
percolator-rs pride prepare PXD032157 --pin /scratch/export.pin \
  --recipe recipe.json --retention keep-if-pinned --dry-run
# After inspecting the plan, replace --dry-run with --yes.
percolator-rs pride run PXD032157 --prepared exported-pin \
  --ephemeral --yes -- --profile fast
```

Example `recipe.json` (replace the input ID and tool details with the actual operation):

```json
{
  "steps": [
    {
      "id": "exported-pin",
      "inputs": ["ACTUAL-PRIDE-FILE-ID"],
      "output_sha256": null,
      "kind": "pin",
      "tool": "actual-export-tool",
      "tool_version": "actual-version",
      "parameters": ["actual arguments and configuration references"],
      "protein_database": null,
      "database_sha256": null,
      "decoy_generation": null
    }
  ]
}
```

Each stage names its parents, tool, version and reproduction parameters. Stages can describe RAW → mzML → database search → PIN. A `database_search` stage must include database identity, its SHA-256, and decoy generation. The last stage has kind `pin`; the importer supplies its actual SHA-256. Stage IDs must be unique and reference remote IDs or preceding stages. Missing details remain missing; the software cannot independently validate the user's account of an external computation. Recipes are provenance, never shell-executed commands.

Import streams a bounded copy, detects input changes, validates with the existing parser, and retains the recipe separately. It never deletes the external input; that remains under the researcher's control. Regenerate an evicted prepared PIN from the recipe and import it under a new stage ID, preserving earlier lineage.

## Storage rules and recovery

```text
CACHE/
  .percolator-pride-cache-v1   ownership marker
  cache.lock                 filesystem lock
  index.json                 versioned atomic artifact/pinning index
  objects/                   authoritative content identities or remote identities
  prepared/                  locally hashed derived PINs
  tmp/                       partial downloads/imports/analysis outputs
  manifests/                 metadata, remote history, recipes, experiments, hashes
  results/                   verified final tables
```

The **hard large-data ceiling covers `objects/`, `prepared/`, and `tmp/` together**. It excludes lightweight manifests/index and final results, whose separate sizes are reported. Results necessarily accumulate if retained; each analysis has a configurable output bound (`--max-results-per-input`, default 64 MiB split across four tables). Safety margin defaults to 1 GB. Free-space planning reserves all selected final-output allowances, the peak source batch, and safety margin. It includes existing selected inputs in the working-space budget. Conservative estimates may reject operations that could fit with a smaller user-specified output budget or batch.

There is no disk reservation or cache filling. Plans show selected identities/sizes, download bytes, current usage, ceiling, actual free space, temporary output allowance, final-output allowance, safety margin, expected LRU evictions and expected remaining source data. Conversion/search estimates are not fabricated: these external steps are not automatically executed. Actual output sizes are unknown until analysis.

Downloads use a 64 KiB application buffer and `.part` files. Request timeouts, bounded retries, cancellation, Content-Length and exact Content-Range checks protect transfers. Resume requires an authoritative checksum or a strong ETag; servers ignoring Range restart the file. Every cached reuse rechecks local SHA-256 and supplied authoritative checksums. Checksumming never trusts a filename or index state alone.

Cache operations use filesystem locking to serialize reservations, downloads, processing, pinning and eviction. Non-cooperating processes can still consume filesystem space: downloads/imports recheck free space while writing, analysis has a strict per-output size bound, and disk errors fail the operation. An HTTP request may remain blocked until its timeout after cancellation; `.part` state remains recoverable. A partial can be retried by rerunning fetch or removed by cleanup. Corruption is recorded and fails the operation; retry fetch replaces a corrupt source. Corrupt metadata fails closed and preserves all files; restore a trusted `index.json` backup before proceeding. No state is reconstructed in a way that could silently lose pin protection.

After a crash, cache opening recovers ownership of a completed rename without assuming validity and marks unfinished experiments interrupted. `cache clean-abandoned` removes tracked abandoned partial/temporary files, respecting pins and retention. Files not represented by the managed index are reported as untracked and block unsafe reclamation. Cache directories require an ownership marker and reject symlinks, special files, path traversal and use inside Git trees. Remote names are display/provenance only; generated local identities prevent path escape and filename collisions. No recursive deletion outside the owned cache occurs.

## Retention and returning large-data usage to zero

```sh
percolator-rs pride cache pin PXD032157
percolator-rs pride cache status
percolator-rs pride cache unpin PXD032157
percolator-rs pride cache prune --dry-run
percolator-rs pride cache prune --all-evictable
percolator-rs pride cache purge-data          # preview
percolator-rs pride cache purge-data --yes    # override KEEP, preserve pins
```

| Policy | Behavior for unpinned data |
| --- | --- |
| `keep` | Protected from automatic eviction and ordinary/all-evictable pruning |
| `evict` | Disposable whenever not actively used |
| `keep-if-pinned` | Default: disposable unless a referencing project is pinned |
| `until-result-verified` | Protected until a successful analysis commits verified final results |

Pinning is an independent protection shared across all references to an object. Any pinned referencing dataset protects the shared bytes. A KEEP reference is never weakened by another workflow. Ephemeral completion deletes selected disposable source/PIN objects only after result validation and provenance commit. Failed runs do not claim completion; valid sources remain and unsafe temporary outputs are prunable.

Ordinary `prune` and explicit `--all-evictable` both remove **all eligible objects**, with no retained-usage target. `purge-data --yes` also overrides KEEP and unfinished-PIN retention, while respecting pins. All preserve manifests, remote references, checksums, local hashes, recipes, configuration, lineage and final results. Reports give storage before, object identities deleted, bytes freed, storage remaining, pinned bytes remaining, and confirmation that provenance/results are retained. With no pins or retention protections, large-data usage can reach **zero**.

## Validation

Normal CI uses frozen official metadata plus local HTTP servers and small fixtures; it makes no PRIDE requests. Tests cover accessions, wire parsing, missing values, pagination/retries, filtering/classification, identity/lineage round trips, streaming, resume, interruption, checksum failures, deduplication, budgets, LRU/pins, relocation, dry runs, purge, symlink/traversal refusal, crash recovery, external PIN import, output limits, failure provenance and multi-batch ephemeral processing. An end-to-end test compares PRIDE-run target PSM output byte-for-byte with the ordinary CLI.

```sh
cargo test --release
# Explicitly opt into a separate metadata-only live test:
cargo test --release --test pride live_official_pride_metadata -- --ignored --exact
```

See [the real-project demonstration](PRIDE-demonstration.md) for the executed download/verify/analyze/cleanup sequence and retained hashes.
