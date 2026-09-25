# Concordance Testing

These scripts compare vep-rs output against the reference Perl VEP implementation to verify correctness. They are **optional**: you do not need Perl VEP to use vep-rs.

## Prerequisites

- **Docker** (for running Perl VEP via `ensemblorg/ensembl-vep` image)
- **Perl VEP cache** (downloaded via Perl VEP installer or Docker)
- **JSON cache** (built via `vep-cache-builder`, or provide pre-existing with `--json-cache-dir`)
- **Benchmark VCF** (e.g., ClinVar)
- **Python 3.8+** (for comparison scripts)

## Quick Start

```bash
# Smoke test (5,000 variants)
scripts/concordance/run_concordance.sh --mode smoke --smoke-variants 5000 \
  --perl-cache-dir <perl-cache> --json-cache-dir <json-cache> --fasta <genome.fa>

# Full concordance (all variants)
scripts/concordance/run_concordance.sh --mode full \
  --perl-cache-dir <perl-cache> --json-cache-dir <json-cache> --fasta <genome.fa>
```

Both commands read `--benchmark-dir` (default `./.vep/benchmark/`), a flat directory of
`.vcf`/`.vcf.gz` inputs; after `scripts/data/download_real_world_vcfs.sh`, pass
`--benchmark-dir .vep/vcf/grch37/all_variants`. The structural-variant comparator,
`scripts/validation/compare_sv_concordance.py`, lives beside the validation scripts and is
documented in `tests/sv_validation/README.md`.

## Scripts

| Script                       | Purpose                                                                                                                                                                                                                                                                                |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `run_clone_measurement.sh`   | Measurement entrypoint: annotates every benchmark dataset for one engine (`--engine {vep-rs,fastvep,perl}`) against a populated data tree, with a discarded cold warmup (its output flushed and deleted before the timed run) plus one timed run, emitting `wall_times.csv` and per-dataset concordance reports; the wall times and F1 values under `manuscript/data/` are its outputs. |
| `test_run_clone_measurement.sh` | Adversarial tests for `run_clone_measurement.sh`: stub engines inject faults to prove each failure-propagation gate fires, fork parity included                                                                                                                                       |
| `run_concordance.sh`         | Development-time single-VCF orchestrator: runs Perl VEP (Docker) + Rust on one input, with plugin testing and perf telemetry. It does not wrap `run_clone_measurement.sh`.                                                                                                              |
| `compare_vep_outputs.py`     | SNP/indel comparator: extracts `(Location, Allele, Feature, Feature_type, Consequence)` tuples from both outputs and computes precision, recall, F1; its `EXCLUSION_REGISTRY` is the adjusted-F1 mask                                                                                    |
| `compare_vep_fields.py`      | Per-column and per-`Extra`-key agreement over the shared consequence keys of two default-format outputs; writes `fields.json` and reports duplicate keys on each side |
| `test_compare_vep_outputs.py` | Unit tests for the SNP/indel comparator and its exclusion registry                                                                                                                                                                                                                    |
| `compare_output_formats.py`  | Output-format self-parity: one engine's VCF, tab, JSON and Parquet outputs must carry the same consequence calls                                                                                                                                                                       |
| `test_compare_output_formats.py` | Unit tests for the format-parity harness                                                                                                                                                                                                                                          |
| `test_compare_vep_fields.py` | Unit tests for the per-column comparator                                                                                                                                                                                                                                               |
| `fastvep_normalize.py`       | Rewrites fastVEP's raw tab output so both engines' rows key alike before comparison                                                                                                                                                                                                    |
| `count_warnings.py`          | Aggregates per-suite ERROR/WARN/INFO counts from both engines' stderr logs into `reports/warning_counts.json` and a Markdown section (`--out-md`)                                                                                                                                      |
| `cache_parity_audit.pl`      | Audits the Storable-to-JSON cache conversion for field coverage over a transcript sample                                                                                                                                                                                               |
| `compare_plugin_outputs.py`  | Plugin-specific comparison with allele-aware field matching                                                                                                                                                                                                                            |
| `test_compare_plugin_outputs.py` | Unit tests for the plugin comparator                                                                                                                                                                                                                                              |
| `analyze_all_discordants.py` | Classifies the `discordant.tsv` of several suites and builds a cross-suite pattern map. Reads `<out>/suites/<id>/report(s)/discordant.tsv`, the file a `run_concordance.sh` run with `--work-dir <out>/suites/<id>` writes; see "Suite severity" below.                                    |
| `compare_intermediates.py`   | Joins a `trace_perl_intermediates.pl` TSV with vep-rs output on (location, allele, transcript) and reports which fields diverge                                                                                                                                                      |
| `trace_perl_intermediates.pl` | Runs inside the Perl VEP container and dumps per-transcript intermediate values (CDS and translation coordinates, codons, peptides, consequence set) for a VCF of discordant variants                                                                                                                                  |
| `trace_perl_variant.sh`      | Run Perl VEP on a single variant with injected `warn()` tracing in VariationEffect.pm, TranscriptVariationAllele.pm, and BaseTranscriptVariation.pm. Produces per-transcript predicate trace.                                                                                          |
| `generate_reference_output.sh` | Writes the Perl VEP ground truth `run_clone_measurement.sh` scores against: one `output.txt` per SNP/indel suite and one `<vcf>.txt` per SV input, via Docker, with the harness's own invocation. See below.                                                                        |
| `finalize_sv_gt_set.py` | Finalizes a per-VCF SV ground-truth set directory after a file is regenerated: refreshes that file's canonical-inputs manifest entry from the file, and writes the set's `provenance.json` (per-file digests of every reference output and canonical input, both assemblies), `PROVENANCE` and `STATUS` from the per-file sidecars `generate_reference_output.sh` writes. |
| `test_generate_reference_output.sh` | Tests for `generate_reference_output.sh`: a stub `docker` exercises the dry run and every fail-closed path                                                                                                                                                                       |
| `extract_perl_reference.sh`  | Extracts five Perl VEP modules from the Docker image into a local directory for debugging and tracing.                                                                                                                                                                                 |
| `prepare_benchmark_vcfs.py`  | Prepares input VCFs: decompresses bgzipped files, normalizes chromosome names, truncates for smoke mode                                                                                                                                                                                |

## Generating the Perl ground truth

`run_clone_measurement.sh` only reads the ground-truth trees under `<data>/ground_truth/perl/`.
`generate_reference_output.sh` writes them, with the same Docker image and the same `vep`
flags the harness's Perl engine uses, so the reference and the timed runs cannot diverge:

```bash
# Preflight everything and print the docker command for every output, running nothing.
scripts/concordance/generate_reference_output.sh --data-dir <data> --dry-run

# Generate both assemblies' SNP/indel and SV references (hours: ClinVar and gnomAD are large).
scripts/concordance/generate_reference_output.sh --data-dir <data>

# One assembly's SV set only, into an explicit location.
scripts/concordance/generate_reference_output.sh --data-dir <data> --assembly GRCh37 --suites sv \
  --sv-gt-dir <data>/ground_truth/perl/sv_per_vcf
```

What it needs, all local (it fetches nothing and uses no cloud credentials):

- Docker and the `ensemblorg/ensembl-vep:release_115.2` image (`--perl-image` to override).
- The **full-genome** Ensembl VEP cache for each assembly, extracted under one directory so
  it holds `homo_sapiens/115_GRCh37/` and `homo_sapiens/115_GRCh38/` (`--perl-cache-dir`,
  default `<data>/caches/perl/vep-cache`). Ensembl FTP:
  `https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz`
  and `https://ftp.ensembl.org/pub/grch37/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh37.tar.gz`
  (about 24 GB each). A chromosome-subset cache is refused before anything runs: it does
  not fail a VEP run, it under-annotates and exits 0.
- An indexed reference FASTA per assembly (`--fasta-grch37` / `--fasta-grch38`, or
  `--fasta` with a single `--assembly`; default `<data>/reference/<asm>/genome.fa`). VEP
  receives `--fasta` on every call.
- The canonical inputs: `vcf/<asm>/all_variants_canonical/<base>.canonical.vcf.gz` for the
  six SNP/indel suites and `<sv-gt-dir>/canonical_inputs/<asm>/*.vcf[.gz]` for the SV set.

Fork and buffer are part of the ground truth's identity: SNP/indel runs `--fork 16
--buffer_size 5000` and the SV set runs `--fork 4 --buffer_size 5000`, the settings of the
released references. VEP's transcript selection above `--max_sv_size` is
batch-dependent, so a different SV fork count changes the row count. `--fork` and `--sv-fork`
override them, but a reference generated at other settings is a different reference.

Every run is gated on annotation volume (rows >= input records, per suite and per SV set,
the same floor the harness's `assert_perl_annotated` applies), refuses to overwrite an
existing reference without `--force`, and writes a `<name>.provenance.json` beside each
output recording the image digest, flags, input and output digests and counts.
`test_generate_reference_output.sh` exercises the dry run and every fail-closed path with
a stub `docker`.

## Output

After a concordance run, results are in `<work-dir>/reports/` (default `tmp/concordance/reports/`):

- `summary.json`: machine-readable metrics (precision, recall, F1, variant counts)
- `summary.md`: human-readable report
- `discordant.tsv`: all variants where Perl and Rust disagree
- `provenance.json`: exact versions, flags, and cache checksums used

## How It Works

1. **Input preparation**: `prepare_benchmark_vcfs.py` normalizes and optionally truncates the input VCF
2. **Perl VEP**: Runs via Docker with matching flags (offline, cache, same plugins)
3. **Rust VEP**: Runs the locally-built `target/release/vep` binary with equivalent flags
4. **Comparison**: `compare_vep_outputs.py` extracts semantic tuples from both outputs and computes set-based metrics
5. **Reporting**: Results written to `<work-dir>/reports/`

## Suite severity

`analyze_all_discordants.py` reads an optional `<out>/dashboard.json` of the form
`{"suites": [{"id": "<suite-id>", "name": "<display name>", "severity": "<severity>"}]}`.
The severity decides how many representative discordant rows the pattern map keeps
per category and suite:

| Severity   | Discordant samples retained |
| ---------- | --------------------------- |
| `blocking` | 10                          |
| `major`    | 5                           |
| `minor`    | 5                           |

A suite with no severity set, or no `dashboard.json` at all, is treated as `minor`.
