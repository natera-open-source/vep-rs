# Concordance Testing

These scripts compare vep-rs output against the reference Ensembl VEP implementation to verify correctness. They are **optional**: you do not need Ensembl VEP to use vep-rs.

## Prerequisites

- **Docker** (for running Ensembl VEP via `ensemblorg/ensembl-vep` image)
- **Ensembl VEP cache** (downloaded via the Ensembl VEP installer or Docker)
- **JSON cache** (built via `vep-cache-builder`, or provide pre-existing with `--json-cache-dir`)
- **Benchmark VCF** (e.g., ClinVar)
- **Python 3.10+** (for comparison scripts), with **`pyarrow`** for the Parquet format
  (`compare_vep_formats.py --format parquet` and `--parquet` read the writer's files through
  `scripts/adapters/parquet_to_vep_tab.py`)

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
| `run_clone_measurement.sh`   | Measurement entrypoint: annotates every benchmark dataset for one engine (`--engine {vep-rs,fastvep,perl}`) against a populated data tree, with a discarded cold warmup (its output flushed and deleted before the timed run) plus one timed run, emitting `wall_times.csv` and per-dataset concordance reports; the wall times and F1 values under `manuscript/data/` are its outputs. Each comparator call carries `--reference-release`, read per cell from the ground truth's provenance (or `--reference-release`); a ground truth naming no release fails the cell. `--fastvep-extra-args "<flags>"` (fastVEP only) appends words, after the harness's own arguments, to every fastVEP `annotate` run, the warmup and the timed run of each SNP/indel cell and each input of the SV per-VCF loop, through `<output-dir>/fastvep-wrapper.sh`: each word is an option (`-x`, `--name`, `--name=value`) or a value of `[A-Za-z0-9_=.,:/-]` that opens with no dash, and any other word stops the run before anything is written; a run given words writes `<output-dir>/provenance.json` (the binary, `fastvep_extra_args` as the words joined by single spaces, `fastvep_extra_argv`, `fastvep_wrapper`, `fastvep_executed_argv`), and a run given none writes no wrapper and no record. Words that carry `--sa-dir` take cells of one assembly, one invocation per assembly: see [fastVEP with its annotation databases](#fastvep-with-its-annotation-databases). |
| `test_run_clone_measurement.sh` | Adversarial tests for `run_clone_measurement.sh`: stub engines inject faults to prove each failure-propagation gate fires, fork parity included; the FASTVEP-EXTRA-ARGS case proves that `--fastvep-extra-args` reaches every fastVEP `annotate` run and the record, each refusal (another engine, a word outside the rule, an output directory whose path holds `=` and, without the words, a binary whose path does, `--sa-dir` over cells of two assemblies), and that a run without the words writes no wrapper and no record |
| `run_concordance.sh`         | Development-time single-VCF orchestrator: runs Ensembl VEP (Docker) + Rust on one input, with plugin testing and perf telemetry. It does not wrap `run_clone_measurement.sh`. The comparator's `--reference-release` is derived from the Perl image tag, the fixtures directory's provenance or the local checkout's `Constants.pm`, or given with `--reference-release`.      |
| `compare_vep_outputs.py`     | SNP/indel comparator: extracts `(Location, Allele, Feature, Feature_type, Consequence)` tuples from both outputs and computes precision, recall, F1; its `EXCLUSION_REGISTRY` is the adjusted-F1 pair mask under both reference releases, the 116.2 registry adds the pair rules that read the pair's two rows in the outputs (`row_pair_rules`), one of them with the record's alleles from the first column or the input records (`--input-vcf`, `--input-dir`), and `--reference-release {115.2,116.2}` (default 116.2, `REFERENCE_REGISTRIES`) selects the one-sided classes: the paper's scope class under 115.2, the record Ensembl VEP 116 skips and the three breakend shapes of its chromosome-aware overlap under 116.2 ([docs/concordance-provenance/116.2-divergence-classes.md](../../docs/concordance-provenance/116.2-divergence-classes.md)) |
| `../../tools/vep-outputs-compare/` | A Rust build of `compare_vep_outputs.py` with the same command line in its three modes and the same report files byte for byte (`cargo test` there holds the binary to the goldens the Python wrote over its fixtures); the Python stays the reference implementation. The harness scripts score with this tool unless `VEP_OUTPUTS_COMPARATOR=python` selects the Python. See its [README](../../tools/vep-outputs-compare/README.md) |
| `reference_release.py`       | Derives a comparator's `--reference-release` from a reference set's provenance (the Ensembl VEP image tag, else the cache version), an image reference or a checkout's `Constants.pm`; fails loudly when none names a release                                                             |
| `test_reference_release.py`  | Unit tests for the release derivation                                                                                                                                                                                                                                                  |
| `compare_vep_fields.py`      | Per-column and per-`Extra`-key agreement over the shared consequence keys of two default-format outputs; writes `fields.json` and reports duplicate keys on each side |
| `score_plugin_tuples.py`     | Plugin tuple scorer: one tuple per populated plugin field of an output row, `(Location, Allele, Feature, field, value)`, each field resolved by the registry of `compare_plugin_outputs.py` among the plugins the scored run loaded (`--plugins`, by registry name or alias) or by a `--custom-source` prefix, and every other Extra key outside Ensembl VEP's own listed as unclaimed, as is a field two of the run's plugins write where no custom source claims it; `--run-command` hands the registry the run's command line, which refuses a run whose flags make a loaded plugin's column name one of Ensembl VEP's own keys as well; scores each `--engine` (a default-format output, or the five-column tuple file a comparator's normaliser writes) against one Ensembl VEP default-format output, streaming through `sort`, and writes per source, engine and projection (with and without the Feature) the integer counts a raw and an adjusted F1 are derived from, `text_only_differences`, `residual_rows` and `shared_keys`, rows that share a (Location, Allele, Feature) never collapsed: a tuple counts as often as its side's rows write it; `--class-table <tsv>`, `plugin_class_table.tsv` where a run names none, declares the documented classes with the plugins and fields each belongs to, `--classes` names the ones applied, every one but the rows that read `not met` where a run names none, and each is applied to its own plugin's fields alone through the predicates of `compare_vep_formats.py`, its count per entry under `masked`; `--gwas-catalog <file>` hands the classes of GWAS the associations file the run's plugin read, without which they set nothing aside; `--utrannotator-sweep <path>` hands the classes of UTRAnnotator Ensembl VEP's outputs of the input under fixed Perl hash seeds, one file a run or a directory of them, likewise; `--json-cache <dir>` hands the classes of Downstream and of Blosum62 the candidate cache whose sequences they recompute from, the class of MaxEntScan the cache's introns, and `spliceai_multi_entry_cell` its transcripts' symbols, from which the class takes the gene of a row that writes no `SYMBOL`, likewise; `--spliceai-table <file>`, once per score file, hands the classes of SpliceAI the files the run's plugin read, likewise, and the summary's last lines then count the rows `spliceai_multi_entry_cell` was asked about and refused, by the ground of the refusal; a predicate that reads a record's other rows is handed them where the `Uploaded_variation` spells the record's alleles, and under an identifier where `--input <vcf>`, the run's input, names its record, the one line that carries the identifier or of several the one that writes the row at its Location with its allele: the reference's rows of the record, and for a sibling-allele class a default-format engine's rows of it as well; `plugin_value_follows_splice_term_past_the_record_span` takes from the same input the reference bases the allele of a row under an identifier stands in for, read by position from the lines that carry the identifier, and reads the row without its record where no input is given; `plugin_value_on_reference_skipped_record` is decided without a reference row, from both engines' rows at a value's Location and the first column of the reference, by the routines of `compare_vep_outputs.py`; `--residuals <dir>` writes the tuples outside the intersection that no applied class set aside, one file per source, engine and projection; `test_score_plugin_tuples.py` holds its unit tests |
| `plugin_class_table.tsv`     | The class table `score_plugin_tuples.py` reads unless `--class-table` names another: one row per documented class of plugin fields, which declares the class, the plugins it belongs to and the field names it reads on each, then `met` or `not met` and the columns `commit`, `genome build` and `note`, which a report carries, and a row that reads `not met` is declared and not applied unless `--classes` names it; `test_score_plugin_tuples.py` holds every row to the comparator's classes and the scorer's reader |
| `../../tools/vep-plugin-pairs/` | The scorer's extraction, sort and pairing in Rust, on every core, run by `score_plugin_tuples.py --pairing-tool <binary>`: it hands the scorer its counts and every value it did not pair as equal text, with both sides' rows at the Locations of those values, and the scorer settles those by its own value rule and classes, so the report and the residual files are the same byte for byte on either path (`test_score_plugin_tuples.py` holds the two to each other); the Python path stays the reference implementation. See its [README](../../tools/vep-plugin-pairs/README.md) |
| `plugin_gate_runs.toml`      | The runs a release scores its supported plugins on, as data: four whole runs on the genome-wide ClinVar datasets for the plugins verified on real data, and one run per set of every plugin corpus, each naming the derived rows it scores, the words of both engines, the image of Ensembl VEP's reference and the scorer's words, with every place outside the tree a named input |
| `reproduce_plugin_concordance.py` | The driver of `plugin_gate_runs.toml`: for a run, a corpus or a derived row it checks the data files by name, takes Ensembl VEP's output by path or writes it from the image `generate_reference_output.sh` pins, runs `vep`, scores with `score_plugin_tuples.py` (through `--pairing-tool` where given) and prints each row's tuple counts and raw and adjusted F1 as plugin record entries; `--dry-run` prints every command. [Reproducing the plugin concordance](../../docs/plugins/reproducing-concordance.md) is its runbook. |
| `test_plugin_gate_runs.py`, `test_reproduce_plugin_concordance.py` | Tests for the data file against the derived rows, the corpora and the data layout, and for the driver's dry runs, refusals and printed rows |
| `test_compare_vep_outputs.py` | Unit tests for the SNP/indel comparator and its two reference registries, including the report the 115.2 registry must reproduce byte for byte (`testdata/reference_release/`); one test runs `tools/vep-outputs-compare` beside it over fabricated output pairs and holds both to the same reports, standard output and exit status (it skips without the tool's binary)                                                                                                            |
| `compare_vep_formats.py`     | Per-column comparator for every output format: one Ensembl VEP output against a candidate's in the default, tab, VCF, JSON or Parquet format (a Parquet actual, a file or a directory of parts, is projected to tab rows by `scripts/adapters/parquet_to_vep_tab.py` and compared with a tab expected as a multiset), every row, CSQ entry and JSON object, every field verbatim. Intended-divergence `--masks` come from a corpus manifest plus the built-in classes: `regulatory_feature_rows` sets the reference's RegulatoryFeature and MotifFeature rows aside, `order_only` lets the set-valued columns `DOMAINS`, `FREQS`, `ZYG` and `NearestExonJB` agree as multisets and the `5UTR_annotation` groups by their pairs, `individual_allele_order` compares an individual's alleles as a set, `conservation_short_name_spelling` reads the reference's `phlyoP<N>way` key as `phyloP<N>way`, `pc_description_text` reads its `PC` header line with the candidate's description of the share the field reports, `motif_score_reverse_strand_mnv` accepts a `MOTIF_SCORE_CHANGE` the candidate writes on a reverse-strand motif row whose allele is two or more bases, the registry pair classes `start_cooccurrence_swap` and `splice_family_swap` set aside the pairs the SNP/indel comparator's excluding rules admit, and one field or warnings class per corrected behaviour of `docs/intended-divergences.md` accepts that behaviour's one shape of difference, every such agreement or rewrite counted in the report under `masks.applied`. `--floors` names a default-format column and a minimum populated count on both sides; the JSON format counts the column under its JSON key or, for a plugin's own column, inside the plugin object that nests it, and a column it has no row count for (a record-level column, `Feature_type`, `Extra`, a custom source's columns) is reported there as not applicable, as is a tab or VCF column neither header declares. `--field-pairs` accepts named value pairs, and `--nearest-sites` settles `NEAREST` by transcription start site through the candidate's JSON cache. `--input` keys the VCF lines and JSON objects by input record through a line's first five columns, since Ensembl VEP writes a structural variant's empty or `0` QUAL, FILTER or INFO back as `.`, and a default or tab comparison reads the input records where a class reads more than a row carries: under the three structural-variant classes (`reference_skipped_record`, which names a record from them and refuses a run without them, `breakend_mate_local_read` and `bnd_own_chromosome_rows_lost`), under a class that reads a record's alleles, which finds the one line that carries a row's identifier, and under `splice_term_past_the_record_span` with its plugin value class, which take the reference bases the allele of a row under an identifier stands in for, since such a row's first column spells no alleles; `--unordered` compares the shared entries without their order, for a flag set whose rows Ensembl VEP emits in hash order (`"unordered": true` in the corpus manifest). `--self-parity` holds one engine's four formats to each other column by column; `--warnings` compares two engines' `_warnings.txt` and skipped-variants files; a JSON report with a `verdict`, exit 1 on red |
| `test_compare_vep_formats.py` | Unit tests for the per-column comparator on the committed `tests/golden/116/GRCh38-release116` outputs; its hand-over tests compare planted outputs on the whole outputs and through `tools/vep-formats-remainder` and hold both runs to the same report, standard output and exit code (they skip without the tool's release build) |
| `../../tools/vep-formats-remainder/` | The front end of `compare_vep_formats.py` for the default and tab formats: a Rust tool that reads both outputs on every core, keeps back the record groups (the rows of one Location) written alike and writes the rest for the Python comparator, which holds every rule. `compare_vep_formats.py` hands a comparison over when `VEP_FORMATS_REMAINDER` names the tool's executable, with the same report and standard output byte for byte; unset, it compares the whole outputs. See its [README](../../tools/vep-formats-remainder/README.md) |
| `compare_vep_stats.py`       | Section-by-section comparator of two statistics reports (the `<output>_summary.txt` that `--stats_text` writes; a corpus keeps Ensembl VEP's as `expected/<format>.summary.log`): `[VEP run statistics]`, `[Data version]` and `[General statistics]` row by row, the command line, the two times, the run time and the input and output rows compared by key alone and the paths of `Annotation sources` masked; every chart section as a set of `(key, count)` pairs, since Ensembl VEP sorts ties in hash order; `--tolerance N` lets a chart count differ by up to N; `--mask-row` sets one more row aside; the same sections in the same order on both sides; a JSON report with a `verdict`, exit 1 on red. `test_compare_vep_stats.py` holds it to the committed `GRCh38-statistics` report; the golden harness applies the same rules to a manifest run's `stats_file` |
| `test_compare_vep_fields.py` | Unit tests for `compare_vep_fields.py`, the per-column comparator of two default-format outputs                                                                                                                                                                                       |
| `fastvep_normalize.py`       | Rewrites fastVEP's raw tab output so both engines' rows key alike before comparison                                                                                                                                                                                                    |
| `fastvep_plugin_tuples.py`   | Writes the supplementary columns of fastVEP's tab output (one per source loaded with `--sa-dir`) as plugin tuples (Location, Allele, Feature, Ensembl VEP field name, decoded value), the form in which a comparator's values are scored against Ensembl VEP's plugin and `--custom` fields. Its data file `fastvep_field_map.tsv` names the Ensembl VEP field of each part it maps and the expected set a fixture source is scored against; `test_fastvep_plugin_tuples.py` holds its tests, lines of the binary's own output among them |
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
| `feature_oracle.sh`          | The comparison oracle: the harness that runs Ensembl VEP and vep-rs on the same inputs and compares every column of every row in every output format. Runs the comparison cases of a case table (input, cache, species, FASTA, flags, formats, Ensembl VEP image and mode, masks as classes or `<class>=<expected applied count>`, floors, field pairs, unordered, reference and refusal per case; a `-` species is homo_sapiens, and the cache slice is read under `<species>/<version>_<assembly>`): Ensembl VEP in Docker once per case and format, offline or against Ensembl's database, then each candidate `vep` with the same flags (less `--genomes`, which selects Ensembl VEP's server) and `--species`, every format compared by `compare_vep_formats.py`; writes `<out>/<candidate_pin>/<case>/report.json` with a `verdict` and `<out>/status.tsv` (`case`, `candidate_pin`, `verdict`, `report_digest`, `run_utc`); a case's `reference` directory of stored outputs in the golden-corpus layout replaces its Ensembl VEP runs, as `--reference-dir` does for every case; a `refusal` directory of the two engines' expected error texts makes a case green when both refuse; `unordered` compares the rows as a multiset; exit 1 on any red. The header of the script documents the table's columns. |
| `test_feature_oracle.sh`     | Tests for `feature_oracle.sh`: a stub `docker` and a stub `vep` answer with the committed `tests/golden/116/GRCh38-release116` outputs, so the real comparator decides each verdict; covers the table contract, the command shapes, every fail-closed path, the stored-reference, refusal, layout and species cases, and the verdicts                                             |
| `test_species_harness.py`    | The species as an input of the harness: `feature_oracle.sh` plans both engines under the case's species and slice, `generate_reference_output.sh` under a corpus manifest's `species` and `release`, both defaulting to homo_sapiens; each of the four other species has one committed corpus under `tests/golden/116/` whose manifest, cache and reference record name its species, assembly and cache version, and whose reference generation is planned under them |
| `extract_perl_reference.sh`  | Extracts five Ensembl VEP modules from the Docker image into a local directory for debugging and tracing.                                                                                                                                                                              |
| `prepare_benchmark_vcfs.py`  | Prepares input VCFs: decompresses bgzipped files, normalizes chromosome names, truncates for smoke mode                                                                                                                                                                                |

## Comparing regulatory feature rows

Both field comparators accept `--keep-only regulatory_feature_rows` to compare only
`RegulatoryFeature` and `MotifFeature` consequence rows from **both** inputs:

```bash
python3 scripts/concordance/compare_vep_fields.py \
  --perl perl.txt --rust rust.txt --out regulatory-fields.json \
  --keep-only regulatory_feature_rows

python3 scripts/concordance/compare_vep_formats.py \
  --format json --expected perl.jsonl --actual rust.jsonl --input variants.vcf \
  --report regulatory-json.json --keep-only regulatory_feature_rows
```

The formats comparator supports this selector for `default`, `tab`, `vcf`, `json` and
`parquet` (a tab-format expected file against the Parquet actual). VCF and JSON retain
record-level entries only for records with retained regulatory or motif rows; on a JSON
record that loses non-regulatory rows, `most_severe_consequence` is recomputed from the
remaining rows on each side. Other record fields and every retained consequence field
still compare normally, including `HIGH_INF_POS` and `MOTIF_SCORE_CHANGE` if present.
Floors apply to the filtered entries; set an appropriate populated-count floor when an
empty regulatory comparison should fail. Both reports record the selector as `keep_only`.
The fields comparator still reports one-sided keys without comparing their fields and
exits zero after writing its report; the formats comparator exits by its green/red verdict.

The selector is available only in the formats comparator's comparison mode, not
`--self-parity` or `--warnings`, and cannot be combined with
`--masks regulatory_feature_rows`. That mask instead excludes the **reference's**
regulatory rows for caches without regulatory data. Its warnings comparison accepts
either the cache-unavailable notice or the notice of releases 0.1.0 to 0.5.0.

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
  default `<data>/caches/perl/vep-cache`; another species' cache extracts beside them as
  `<species>/<version>_<assembly>/`, the layout Ensembl VEP resolves from `--species`,
  `--cache_version` and `--assembly`). Ensembl FTP:
  `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz`
  and `https://ftp.ensembl.org/pub/grch37/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh37.tar.gz`
  (about 24 GB each). A chromosome-subset cache is refused before anything runs: it does
  not fail a VEP run, it under-annotates and exits 0; the primary set is known for
  homo_sapiens, and another species' slice must hold at least one chromosome directory.
- An indexed reference FASTA per genome (`--fasta-for <species>:<assembly>=<path>`,
  repeatable, or `--fasta` with a single `--assembly`; `homo_sapiens:GRCh37` and
  `homo_sapiens:GRCh38` default to `<data>/reference/<asm>/genome.fa`). VEP receives
  `--fasta` on every call.
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
  and abbreviated spellings Ensembl VEP's option parser accepts as well. The one
  exception is `--format <input format>`, which names the input the words read (an
  Ensembl-format, region or HGVS file) and stands in place of the script's own
  `--format vcf`.
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
and its `provenance.json` from `<dir>/variants.vcf`, with the manifest's `species`
(homo_sapiens when absent), `assembly`, `release` (the cache version), `flags` and `fasta`,
`--fork 1` and every format's headers kept, the layout `crates/vep-cli/tests/golden.rs` and
`format_parity.rs` read:

```bash
scripts/concordance/generate_reference_output.sh --corpus tests/golden/116/GRCh38-release116 \
  --perl-cache-dir <perl-cache>
```

`provenance.json` carries the image digest, the vep command and the sha256 of every
expected file's uncompressed bytes; a key the script does not own (the vep-rs binary the
divergences were classified against, the pruned-cache summary) survives a regeneration,
which needs `--force`, and `input` and `fasta` are filled only when absent, so a record's
own description of them stands. The warnings and statistics reports land
beside the expected files as `<format>.warnings.log` and `<format>.summary.log`: a corpus
run keeps the statistics report (it passes `--stats_text`) and records each report's
sha256 beside the outputs', and `--no-stats` is the opt-out that passes `--no_stats`.

| Placeholder | In a manifest flag word (`generate_reference_output.sh --corpus`, the golden tests) | In a case-table flags cell (`feature_oracle.sh`) |
| ----------- | ------------------------------------------------------------------------------------- | ------------------------------------------------- |
| `{corpus}`  | The corpus directory: Ensembl VEP reads a read-only mount of it at `/work/corpus`, vep-rs the directory itself, so `--custom file={corpus}/custom/sites.vcf.gz,...` names a file beside the manifest (`scripts/golden/make_custom_fixtures.py` writes those) | The input's directory, the same two ways; `{sources}` is `--sources-dir`, mounted read-only at `/work/sources` for Ensembl VEP; a `--dir` flag names its directory through one of the two, and one holding a `vep.ini` is refused |

## fastVEP with its annotation databases

`run_clone_measurement.sh --engine fastvep` measures fastVEP, a comparator, and
`--fastvep-extra-args "<flags>"` appends its words to every `fastvep annotate` run of the
measurement. fastVEP reads its annotation databases from the directory `--sa-dir` names and
loads every database in it, and a directory holds the databases of one assembly. A
measurement with databases is therefore one invocation per assembly, each naming that
assembly's cells and directory: words that carry `--sa-dir` over cells of two assemblies
stop the run with one `ERROR: [run_clone_measurement]` line before any cell runs. A table
measured apart from the others, as the pLI gene table is, sits in a directory of its own
and takes one further invocation:

```bash
# GRCh37 cells against the GRCh37 database directory.
scripts/concordance/run_clone_measurement.sh --engine fastvep --fastvep-binary <fastvep> \
  --data-dir <data> --output-dir <out>/grch37 --cells <GRCh37 cells> \
  --fastvep-extra-args "--exclude-readthrough --exclude-artifact --sa-dir <sa>/grch37"

# GRCh38 cells against the GRCh38 database directory.
scripts/concordance/run_clone_measurement.sh --engine fastvep --fastvep-binary <fastvep> \
  --data-dir <data> --output-dir <out>/grch38 --cells <GRCh38 cells> \
  --fastvep-extra-args "--exclude-readthrough --exclude-artifact --sa-dir <sa>/grch38"

# pLI: cells of one assembly against the directory that holds the pLI table alone.
scripts/concordance/run_clone_measurement.sh --engine fastvep --fastvep-binary <fastvep> \
  --data-dir <data> --output-dir <out>/pli --cells <cells> \
  --fastvep-extra-args "--exclude-readthrough --exclude-artifact --sa-dir <sa>/pli"
```

`--cells` takes a comma-separated list of cell ids; the script's cell table names each
cell's id and assembly, and `--help` shows the form. Each invocation writes its own
`<output-dir>/provenance.json`: `fastvep_extra_args` is the words joined by single spaces,
and `fastvep_executed_argv` lists every argument vector fastVEP ran with. A directory is
one word of `[A-Za-z0-9_=.,:/-]`, given after `--sa-dir` or as `--sa-dir=<directory>`, so
its path holds no space. With words, an `--output-dir` whose path holds `=` is refused
before anything is written: the wrapper runs from that directory under `env`, which reads
a path holding `=` as an assignment. Without words `env` starts the binary itself, so a
`--fastvep-binary` path that holds `=` is refused the same way.

## Output

After a concordance run, results are in `<work-dir>/reports/` (default `tmp/concordance/reports/`):

- `summary.json`: machine-readable metrics (precision, recall, F1, variant counts)
- `summary.md`: human-readable report
- `discordant.tsv`: all variants where Perl and Rust disagree
- `discordant_open.tsv`: the rows of `discordant.tsv` that no excluding rule removes, in the same order and layout; its row count is the adjusted one-sided total
- `discordant.tsv.partial`: where the tuple comparison stops on a key file that does not ascend, the rows written before the stop stand under this name in place of `discordant.tsv`, and the directory holding the two key files is kept and named in the error line
- `provenance.json`: exact versions, flags, plugin argument strings, and cache checksums used

## How It Works

1. **Input preparation**: `prepare_benchmark_vcfs.py` normalizes and optionally truncates the input VCF
2. **Ensembl VEP**: Runs via Docker with matching flags (offline, cache, same plugins)
3. **Rust VEP**: Runs the locally-built `target/release/vep` binary with equivalent flags
4. **Comparison**: the tuple comparator extracts semantic tuples from both outputs and computes set-based metrics. It is `tools/vep-outputs-compare`, which the run builds through cargo before any engine runs; `VEP_OUTPUTS_COMPARATOR=python` runs `compare_vep_outputs.py`, the reference implementation, which writes the same reports
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
