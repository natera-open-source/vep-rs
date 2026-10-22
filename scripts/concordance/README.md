# Concordance Testing

These scripts compare vep-rs output against the reference Ensembl VEP implementation to verify correctness. They are **optional**: you do not need Ensembl VEP to use vep-rs.

## Prerequisites

- **Docker** (for running Ensembl VEP via `ensemblorg/ensembl-vep` image)
- **Ensembl VEP cache** (downloaded via the Ensembl VEP installer or Docker)
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
| `run_clone_measurement.sh`   | Measurement entrypoint: annotates every benchmark dataset for one engine (`--engine {vep-rs,fastvep,perl}`) against a populated data tree, with a discarded cold warmup (its output flushed and deleted before the timed run) plus one timed run, emitting `wall_times.csv` and per-dataset concordance reports; the wall times and F1 values under `manuscript/data/` are its outputs. Each comparator call carries `--reference-release`, read per cell from the ground truth's provenance (or `--reference-release`); a ground truth naming no release fails the cell. |
| `test_run_clone_measurement.sh` | Adversarial tests for `run_clone_measurement.sh`: stub engines inject faults to prove each failure-propagation gate fires, fork parity included                                                                                                                                       |
| `run_concordance.sh`         | Development-time single-VCF orchestrator: runs Ensembl VEP (Docker) + Rust on one input, with plugin testing and perf telemetry. It does not wrap `run_clone_measurement.sh`. The comparator's `--reference-release` is derived from the Perl image tag, the fixtures directory's provenance or the local checkout's `Constants.pm`, or given with `--reference-release`.      |
| `compare_vep_outputs.py`     | SNP/indel comparator: extracts `(Location, Allele, Feature, Feature_type, Consequence)` tuples from both outputs and computes precision, recall, F1; its `EXCLUSION_REGISTRY` is the adjusted-F1 pair mask under both reference releases, and `--reference-release {115.2,116.2}` (default 116.2, `REFERENCE_REGISTRIES`) selects the one-sided classes: the paper's scope class under 115.2, the record Ensembl VEP 116 skips and the three breakend shapes of its chromosome-aware overlap under 116.2 ([docs/concordance-provenance/116.2-divergence-classes.md](../../docs/concordance-provenance/116.2-divergence-classes.md)) |
| `reference_release.py`       | Derives a comparator's `--reference-release` from a reference set's provenance (the Ensembl VEP image tag, else the cache version), an image reference or a checkout's `Constants.pm`; fails loudly when none names a release                                                             |
| `test_reference_release.py`  | Unit tests for the release derivation                                                                                                                                                                                                                                                  |
| `compare_vep_fields.py`      | Per-column and per-`Extra`-key agreement over the shared consequence keys of two default-format outputs; writes `fields.json` and reports duplicate keys on each side |
| `test_compare_vep_outputs.py` | Unit tests for the SNP/indel comparator and its two reference registries, including the report the 115.2 registry must reproduce byte for byte (`testdata/reference_release/`)                                                                                                            |
| `compare_vep_formats.py`     | Per-column comparator for every output format: one Ensembl VEP output against a candidate's in the default, tab, VCF, JSON or Parquet format, every row, CSQ entry and JSON object, every field verbatim. Intended-divergence `--masks` come from a corpus manifest plus two built-in classes: `regulatory_feature_rows` sets the reference's RegulatoryFeature and MotifFeature rows aside, and `order_only` lets the set-valued columns `DOMAINS`, `FREQS` and `ZYG` agree as multisets, every such agreement counted in the report. `--floors` names a default-format column and a minimum populated count on both sides; the JSON format counts the column under its JSON key, and a column it has no row count for (a record-level column, `Feature_type`, `Extra`) is reported there as not applicable. `--field-pairs` accepts named value pairs, and `--nearest-sites` settles `NEAREST` by transcription start site through the candidate's JSON cache. `--input` keys the VCF lines and JSON objects by input record through a line's first five columns, since Ensembl VEP writes a structural variant's empty or `0` QUAL, FILTER or INFO back as `.`; `--unordered` compares the shared entries without their order, for a flag set whose rows Ensembl VEP emits in hash order (`"unordered": true` in the corpus manifest). `--self-parity` holds one engine's four formats to each other column by column; `--warnings` compares two engines' `_warnings.txt` and skipped-variants files; a JSON report with a `verdict`, exit 1 on red |
| `test_compare_vep_formats.py` | Unit tests for the per-column comparator on the committed `tests/golden/116/GRCh38-release116` outputs                                                                                                                                                                             |
| `test_compare_vep_fields.py` | Unit tests for `compare_vep_fields.py`, the per-column comparator of two default-format outputs                                                                                                                                                                                       |
| `fastvep_normalize.py`       | Rewrites fastVEP's raw tab output so both engines' rows key alike before comparison                                                                                                                                                                                                    |
| `count_warnings.py`          | Aggregates per-suite ERROR/WARN/INFO counts from both engines' stderr logs into `reports/warning_counts.json` and a Markdown section (`--out-md`)                                                                                                                                      |
| `cache_parity_audit.pl`      | Audits the Storable-to-JSON cache conversion for field coverage over a transcript sample                                                                                                                                                                                               |
| `compare_plugin_outputs.py`  | Plugin-specific comparison with allele-aware field matching                                                                                                                                                                                                                            |
| `test_compare_plugin_outputs.py` | Unit tests for the plugin comparator                                                                                                                                                                                                                                              |
| `analyze_all_discordants.py` | Classifies the `discordant.tsv` of several suites and builds a cross-suite pattern map. Reads `<out>/suites/<id>/report(s)/discordant.tsv`, the file a `run_concordance.sh` run with `--work-dir <out>/suites/<id>` writes; see "Suite severity" below.                                    |
| `compare_intermediates.py`   | Joins a `trace_perl_intermediates.pl` TSV with vep-rs output on (location, allele, transcript) and reports which fields diverge                                                                                                                                                      |
| `trace_perl_intermediates.pl` | Runs inside the Ensembl VEP container and dumps per-transcript intermediate values (CDS and translation coordinates, codons, peptides, consequence set) for a VCF of discordant variants                                                                                                                               |
| `trace_perl_variant.sh`      | Run Ensembl VEP on a single variant with injected `warn()` tracing in VariationEffect.pm, TranscriptVariationAllele.pm, and BaseTranscriptVariation.pm. Produces per-transcript predicate trace.                                                                                       |
| `generate_reference_output.sh` | Writes the Ensembl VEP ground truth `run_clone_measurement.sh` scores against (one `output.txt` per SNP/indel suite and one `<vcf>.txt` per SV input, via Docker, with the harness's own invocation), a flag family's reference in any output format (`--flags`, `--format`, one container run per format), and a golden corpus's four expected outputs with its `provenance.json` (`--corpus`); the image digest is asserted before any run and `--stats` keeps Ensembl VEP's statistics report. See below. |
| `finalize_sv_gt_set.py` | Finalizes a per-VCF SV ground-truth set directory after a file is regenerated: refreshes that file's canonical-inputs manifest entry from the file, and writes the set's `provenance.json` (per-file digests of every reference output and canonical input, both assemblies), `PROVENANCE` and `STATUS` from the per-file sidecars `generate_reference_output.sh` writes. |
| `test_generate_reference_output.sh` | Tests for `generate_reference_output.sh`: a stub `docker` exercises the dry run and every fail-closed path                                                                                                                                                                       |
| `feature_oracle.sh`          | Runs the comparison cases of a case table (input, cache, FASTA, flags, formats, Ensembl VEP image and mode, masks, floors, field pairs per case): Ensembl VEP in Docker once per case and format, offline or against Ensembl's database, then each candidate `vep` with the same flags, every format compared by `compare_vep_formats.py`; writes `<out>/<candidate_pin>/<case>/report.json` with a `verdict` and `<out>/status.tsv` (`case`, `candidate_pin`, `verdict`, `report_digest`, `run_utc`); `--reference-dir` replaces the Ensembl VEP runs with stored outputs in the golden-corpus layout; exit 1 on any red. The header of the script documents the table's columns. |
| `test_feature_oracle.sh`     | Tests for `feature_oracle.sh`: a stub `docker` and a stub `vep` answer with the committed `tests/golden/116/GRCh38-release116` outputs, so the real comparator decides each verdict; covers the table contract, the command shapes, every fail-closed path and the verdicts                                                                 |
| `extract_perl_reference.sh`  | Extracts five Ensembl VEP modules from the Docker image into a local directory for debugging and tracing.                                                                                                                                                                              |
| `prepare_benchmark_vcfs.py`  | Prepares input VCFs: decompresses bgzipped files, normalizes chromosome names, truncates for smoke mode                                                                                                                                                                                |

## Generating the Perl ground truth

`run_clone_measurement.sh` only reads the ground-truth trees under `<data>/ground_truth/perl/`.
`generate_reference_output.sh` writes them with the same Docker image (`release_116.2`), cache
version (116) and `vep` flags the harness's Perl engine uses, so the reference and the timed
runs cannot diverge; `--perl-image` and `--cache-version` move both together:

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

- Docker and the `ensemblorg/ensembl-vep:release_116.2` image (`--perl-image` to override).
- The **full-genome** Ensembl VEP cache for each assembly, of the release the image pins (cache
  version 116; `--cache-version` to override), extracted under one directory so it holds
  `homo_sapiens/116_GRCh37/` and `homo_sapiens/116_GRCh38/` (`--perl-cache-dir`,
  default `<data>/caches/perl/vep-cache`). Ensembl FTP:
  `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz`
  and `https://ftp.ensembl.org/pub/grch37/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh37.tar.gz`
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
output recording the image digest, flags, input and output digests and counts. Before
anything runs, the pulled image's digest must match `--perl-image-digest` (default the
`release_116.2` tag's manifest-list digest, or `VEP_PERL_IMAGE_DIGEST`), so a moved tag
cannot produce a reference under this name; `--perl-image` and `--perl-image-digest` move
together. `test_generate_reference_output.sh` exercises the dry run and every fail-closed
path with a stub `docker`.

### Flag families and output formats

Ensembl VEP writes one output format per run (`Runner.pm:871` builds the run's single
`OutputFactory`), so a reference in several formats is one container run per format:

```bash
# The --everything family in all four formats, into its own root beside the harness's.
scripts/concordance/generate_reference_output.sh --data-dir <data> --format all \
  --flags "--everything" \
  --snp-indel-gt-dir <data>/ground_truth/perl/everything/snp_indel \
  --sv-gt-dir <data>/ground_truth/perl/everything/sv_per_vcf --sv-inputs-dir <data>/ground_truth/perl/sv_per_vcf/canonical_inputs
```

- `--flags "<words>"` appends the words, split on whitespace, after the harness's flags
  and before the format flag; a flag the script sets itself (`--fork`, `--fasta`,
  `--no_stats`, `--tab`, ...) is refused rather than passed twice, in the single-dash
  and abbreviated spellings Ensembl VEP's option parser accepts as well.
- `--format default|tab|vcf|json|all` names the output format (Ensembl VEP's own
  `--format vcf` on the command line is its input format). The default format keeps the
  names above; the others take their extension on the same stem (`output.tsv`,
  `output.vcf`, `output.jsonl`; `<vcf-basename>.tsv`, `.vcf`, `.jsonl`), so a `*.txt`
  glob still pairs only the default output. `--no_headers` goes with the default format
  alone: the tab column line and the VCF `##INFO=<ID=CSQ...>` line name the columns of
  those formats. Sidecars take the output's name less a default-format `.txt`
  (`output.provenance.json`, `output.vcf.provenance.json`).
- `--stats` keeps Ensembl VEP's statistics report as `<stem>.summary.log` instead of
  passing `--no_stats` (the run passes `--stats_text`; without it Ensembl VEP writes the
  HTML report alone).
- The volume gate is per format: the default and tab formats count consequence rows
  against input records; the VCF and JSON formats write one line per kept record, so
  their floor is that every written record carries a CSQ field or a consequences array.
  A flag family that drops records or rows (`--chr`, `--coding_only`, `--no_intergenic`,
  `--filter_common`, `--allow_non_variant`, `--individual`, `--individual_zyg`,
  `--dont_skip`), or one that
  writes one row per record (`--pick`, `--most_severe`, `--summary`) on an input with
  records Ensembl VEP skips before annotation, takes `--volume-floor any`, which keeps
  only the one-annotated-unit floor. The VCF gate reads the field by name, so under
  `--keep_csq` a record whose input already carried the field counts as annotated.

### Golden corpus expected outputs

`--corpus <dir>` writes a corpus's `expected/{default.txt,tab.txt,vcf.vcf,json.jsonl}.gz`
and its `provenance.json` from `<dir>/variants.vcf`, with the manifest's `assembly`,
`flags` and `fasta`, `--fork 1` and every format's headers kept, the layout
`crates/vep-cli/tests/golden.rs` and `format_parity.rs` read:

```bash
scripts/concordance/generate_reference_output.sh --corpus tests/golden/116/GRCh38-release116 \
  --perl-cache-dir <perl-cache> --stats
```

`provenance.json` carries the image digest, the vep command and the sha256 of every
expected file's uncompressed bytes; a key the script does not own (the vep-rs binary the
divergences were classified against, the pruned-cache summary) survives a regeneration,
which needs `--force`, and `input` and `fasta` are filled only when absent, so a record's
own description of them stands. The warnings and statistics reports land
beside the expected files as `<format>.warnings.log` and `<format>.summary.log`.

## Output

After a concordance run, results are in `<work-dir>/reports/` (default `tmp/concordance/reports/`):

- `summary.json`: machine-readable metrics (precision, recall, F1, variant counts)
- `summary.md`: human-readable report
- `discordant.tsv`: all variants where Perl and Rust disagree
- `discordant_open.tsv`: the rows of `discordant.tsv` that no excluding rule removes, in the same order and layout; its row count is the adjusted one-sided total
- `provenance.json`: exact versions, flags, and cache checksums used

## How It Works

1. **Input preparation**: `prepare_benchmark_vcfs.py` normalizes and optionally truncates the input VCF
2. **Ensembl VEP**: Runs via Docker with matching flags (offline, cache, same plugins)
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
