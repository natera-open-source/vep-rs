# Scripts

Utility scripts for testing, benchmarking, validation, and data preparation.

## Environment

These scripts read their data locations from environment variables, each with a repo-relative default under `./.vep/` so nothing depends on a specific machine layout. Override any of them to point at your own directories, or pass the equivalent `--flag` where a script exposes one. The optional S3 copy is off unless `VEP_BENCHMARK_S3` or `VEP_BENCHMARK_BUCKET` is set.

| Variable                  | Default                 | Used for                                                     |
| ------------------------- | ----------------------- | ------------------------------------------------------------ |
| `VEP_VCF_DIR`             | `./.vep/vcf`            | Downloaded input VCFs, laid out `<asm>/{sv,all_variants}/` by `data/download_real_world_vcfs.sh`; the canonical input location |
| `VEP_BENCHMARK_DIR`       | `./.vep/benchmark`      | `--benchmark-dir` of `concordance/run_concordance.sh` (a flat directory of VCFs); point it at `./.vep/vcf/<asm>/all_variants` |
| `VEP_PERL_CACHE_DIR`      | `./.vep/cache`          | Ensembl VEP cache (mounted into the Docker container)        |
| `VEP_PERL_REF_DIR`        | `./.vep/perl_reference` | Ensembl VEP modules extracted from Docker for tracing        |
| `VEP_WORK_DIR`            | `./.vep/work`           | Scratch output: concordance runs and Perl trace output       |
| `VEP_BENCHMARK_S3`        | (unset)                 | Optional S3 URI holding a copy of the VCFs `data/download_real_world_vcfs.sh` fetches |
| `VEP_BENCHMARK_BUCKET`    | (unset)                 | Bucket name for `data/setup_plugin_data.sh --source s3`; `data/download_real_world_vcfs.sh` accepts it too |

## Directory Layout

| Directory      | Purpose                                                                       | Prerequisites                 |
| -------------- | ----------------------------------------------------------------------------- | ----------------------------- |
| `concordance/` | Compare vep-rs output against Ensembl VEP                                     | Docker, Ensembl VEP cache     |
| `validation/`  | Standalone test data generation, SV concordance, and cache builder validation | Python 3.10+                  |
| `data/`        | Download input VCFs and stage plugin annotation data                          | Internet access, tabix        |
| `golden/`      | Build the per-release, per-assembly golden corpora under `tests/golden/`      | Python 3.10+; Ensembl VEP outputs and a JSON cache |
| `adapters/`    | Convert vep-rs outputs between formats                                        | `duckdb` CLI                  |
| `release/`     | Render the GitHub release body the release workflow publishes as a draft, and hold its Linux binaries to the glibc floor the page states | Python 3.10+                  |
| `docker/`      | Stage, tag and smoke-test the container image `.github/workflows/docker.yml` publishes | `gh`, Docker            |
| `parity/`      | Hold `tests/parity/matrix.tsv` to the tree: every flag has a row and every landed row's tests resolve | Python 3.10+ |
| `docs/`        | Build the documentation site, the release at its root and the preview of main under `main/`, and check its links, metadata and prose (`.github/workflows/pages.yml`); screenshot it for review | Python 3.10+, mdBook; Playwright for `site_screenshots.py` only |

## Reproducing the published concordance

The published F1 values come from `concordance/run_clone_measurement.sh`, which runs one engine over every benchmark dataset against a populated data tree, not from `run_concordance.sh`. Every published value also re-derives without a run: each raw F1 is `2I/(A+B)` over the tuple counts in `manuscript/data/f1_observations.csv`; the six SNP/indel adjusted values subtract that suite's `masked_pairs` (`manuscript/data/discordance_taxonomy.csv`) from both denominators; the two structural-variant rows carry their own `adj_*` counts ([manuscript/data/README.md](../manuscript/data/README.md)). To re-measure rather than re-derive:

1. **Populate the data tree.** The layout the harness expects under `--data-dir` (each path overridable by the `VEP_*` environment variable named beside its default in the script; there are no path flags):

   ```
   vcf/<asm>/all_variants/<base>.vcf.gz                        SNP/indel downloads (read by --engine perl)
   vcf/<asm>/all_variants_canonical/<base>.canonical.vcf.gz   SNP/indel inputs (read by --engine vep-rs and fastvep)
   ground_truth/perl/snp_indel/<suite>/output.txt             SNP/indel reference output
   ground_truth/perl/sv_per_vcf/canonical_inputs/<asm>/       SV inputs
   ground_truth/perl/sv_per_vcf/<asm>/                        SV reference output, one file per input VCF
   caches/vep-rs/<asm>/                                       vep-rs JSON cache
   caches/fastvep/<asm>.bin                                   fastVEP cache
   caches/perl/vep-cache/                                     Perl Storable cache
   reference/<asm>/genome.fa[.fai]                            reference FASTA
   reference/<asm>/annotation.gff3                            GFF3 annotation
   ```

   `<asm>` is `grch37` or `grch38`; `genome.fa` and `annotation.gff3` are the literal default names, changed only through `VEP_FASTA_GRCH37`/`VEP_FASTA_GRCH38` and `VEP_GFF3_GRCH37`/`VEP_GFF3_GRCH38`. The raw inputs come from `data/download_real_world_vcfs.sh` and the 13 synthetic SV files from `tests/sv_validation/<asm>/`. `all_variants/` is the download script's tree: `--engine perl` reads `all_variants/<base>.vcf.gz` and derives a canonical copy per run with `bcftools`. `all_variants_canonical/` holds the staged `<base>.canonical.vcf.gz` files that `--engine vep-rs`, `--engine fastvep` and `generate_reference_output.sh` (step 2) read; the harness does not derive them, so stage them before either runs. `<base>` is the `vcf_base` column of the suite table in the script.

   "Canonical" means Ensembl-style contig names (`21`, not `chr21`) on the canonical contigs only (`1`-`22`, `X`, `Y`, `MT`, `M`): raw GRCh38 chr21 inputs score F1 = 0 against a reference that names contigs the Ensembl way, so canonicalise every input symmetrically. A staged file is derived from its download the way the harness's own `canonicalize` function derives the Perl copy: strip a `chr` prefix from the records and the `##contig` lines, then `bcftools view -G -t 1,2,...,22,X,Y,MT,M -O z` and `tabix -p vcf`. (`concordance/prepare_benchmark_vcfs.py --canonical-contigs` applies the same contig rule to a directory of inputs; it does not drop genotypes.) `-G` drops per-sample genotypes, which change no consequence: VEP annotates per position and allele, and the tab output carries no genotype column. It is a no-op on the four ClinVar and gnomAD downloads, which are sites-only at source, and it is what the two 1000 Genomes timed cells need, because their downloads (`1kg_integrated_chr21.vcf.gz`, `1kg_highcov_chr21.vcf.gz`) carry genotypes. The suite table names four canonical 1000 Genomes files, two of them under names the download script does not write:

   | `<base>`                         | Assembly | Derived from                  | Genotypes | Measures                                                                              |
   | -------------------------------- | -------- | ----------------------------- | --------- | ------------------------------------------------------------------------------------- |
   | `1kg_integrated_chr21_sitesonly` | GRCh37   | `1kg_integrated_chr21.vcf.gz` | dropped   | the published F1 and wall-time cell                                                   |
   | `1kg_integrated_chr21`           | GRCh37   | `1kg_integrated_chr21.vcf.gz` | kept      | the sample-scaling wall time (`1kg_integrated_chr21_multisample` in `wall_times.csv`) |
   | `1kg_highcov_chr21`              | GRCh38   | `1kg_highcov_chr21.vcf.gz`    | dropped   | the published F1 and wall-time cell                                                   |
   | `1kg_highcov_chr21_multisample`  | GRCh38   | `1kg_highcov_chr21.vcf.gz`    | kept      | the sample-scaling wall time                                                          |

   The genotype-keeping copies take the same contig filter without `-G`. Because `--engine perl` reads `all_variants/<base>.vcf.gz`, that tree needs a `1kg_integrated_chr21_sitesonly.vcf.gz` (a `bcftools view -G` copy of the download) beside the downloads; the Perl branch drops genotypes on every input it canonicalises, so it never times a genotype-carrying form. The sample-scaling rows run only when named in `--suites` and are wall-time cells only: no reference output exists for them, their consequence output is that of the sites-only form, and the released rows for them are `vep-rs` and `fastvep`.

2. **Generate the reference output.** `concordance/generate_reference_output.sh` runs Ensembl VEP (`ensemblorg/ensembl-vep:release_116.2` with cache version 116 by default; `--perl-image` and `--cache-version` override) over the canonical inputs and writes the SNP/indel reference output per suite and the structural-variant reference output per input VCF into the `ground_truth/` layout above, passing `--fasta` at every call site; each of the 16 structural-variant VCFs per assembly is annotated in its own invocation, so no file's transcript-cache population can affect another's. The published values were scored against reference output at cache version 115, except ClinVar GRCh37, whose reference is the cache 113 set that a cache 115 regeneration reproduces key for key (`docs/concordance-provenance/2026-09-18-grch37-cache-113-vs-115.json`).

3. **Convert the cache both engines read.** vep-rs annotates from a JSON conversion of the Storable cache Ensembl VEP reads (`data/storable_to_json.pl`, usage below), which holds the transcript set fixed across the two engines. `concordance/cache_parity_audit.pl` audits the conversion: it reads a seeded sample of transcripts from both caches and classifies every Storable field against the converter's rules (identical, renamed, converted, dropped because undefined, dropped by a named rule, or an unexplained gap); every field the audit classifies is carried or dropped under a named rule, so a conversion made with `data/storable_to_json.pl` audits with no gap. The manifest for the GRCh37 pair the published figures were measured on, the cache-version-113 Storable VEP release 115 ships for that assembly and the JSON conversion of it those figures ran on, is `docs/concordance-provenance/2026-09-16-cache-parity-audit.json` (the manifest names neither release nor assembly; its `json_transcripts_indexed` of 195,232 is the GRCh37 cache's transcript count). It records two Storable fields absent from that JSON side under no rule (`_gene_version` and `_uniprot_isoform`, neither of which that conversion carries and neither of which vep-rs reads for consequence calling) and the SIFT and PolyPhen prediction matrices dropped by a named rule; the concordance tuple carries no prediction field, so none of this enters a published F1. The translation's sequence edits, the `Bio::EnsEMBL::SeqEdit` list VEP applies to the reference allele's peptide at annotation time, are carried as `variation_effect_feature_cache.seq_edits`, one object per edit (`code`, `name`, `start`, `end`, `alt_seq`, `length_diff`; always an array, `[]` when the translation has none), and vep-rs applies them as VEP does. The cached peptide cannot stand in for them: the release 115 Storable caches store `M` as the first residue over any start codon of the codon table whether or not an `initial_met` edit exists. A cache converted without the key loads; on such a cache vep-rs applies no edit, and a variant at codon 1 of a transcript with an alternative initiation codon is called from the cached peptide. The manifest above records the published GRCh37 conversion, which carries no `seq_edits`.

4. **Run the measurement.**

   ```bash
   scripts/concordance/run_clone_measurement.sh \
     --engine vep-rs --vep-binary target/release/vep \
     --data-dir /path/to/data --output-dir /path/to/data/work
   ```

   Each SNP/indel suite is scored by `concordance/compare_vep_outputs.py` and each structural-variant set by `validation/compare_sv_concordance.py`; the per-suite `report/` directories carry the raw and adjusted tuple counts, and `wall_times.csv` the timing. The published cells ran with `--fork 16` on 32-vCPU hosts (the host specification is recorded in `docs/concordance-provenance/2026-09-20-n20-parallel-clones.json`); `--fork` is overridable and F1 does not depend on it, which `crates/vep-cli/src/runner.rs` gates with a fork-1-versus-fork-4 equality test.

   Every vep-rs SNP/indel cell also runs the per-column field comparator, `concordance/compare_vep_fields.py`, which writes `report/fields.json` beside the F1 report (`--skip-fields` drops it; F1 is computed either way). On the genome-wide ClinVar cells that comparator is the longest step of the run, so [`tools/vep-fields-compare`](../tools/vep-fields-compare/README.md) builds the same comparator in Rust: the same command line and a byte-identical `fields.json`, held to the Python's output over committed fixture pairs in CI. Build it with `cargo build --release --locked --manifest-path tools/vep-fields-compare/Cargo.toml` and export `VEP_FIELDS_COMPARATOR=<path to the binary, or a command on PATH>`; the harness then runs it in the Python's place and, when the value names no executable, warns and runs the Python. The Python stays the reference implementation.

F1 is deterministic per binary and architecture-independent: a re-run on the same binary, cache, reference FASTA, canonical inputs and reference output reproduces the committed tuple counts exactly, and any difference points at one of those inputs rather than at the machine.

## Script Index

### Repository gates (this directory)

| Script                          | Description                                                                                                                                  |
| ------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `check_module_docs.sh`          | Every `crates/**/*.rs` file's module documentation survives the two-line licence header (CI job `Module Docs`)                                |
| `check_advisory_reachability.sh` | `cargo audit` advisories are a failure only when reachable from the production binary; `advisory-reachability-allow.txt` holds dated exceptions |

### concordance/

| Script                        | Description                                                                                                                                                             |
| ----------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `run_clone_measurement.sh`    | Measurement entrypoint: runs every benchmark dataset for one engine (`--engine vep-rs\|fastvep\|perl`) against a populated data tree, one discarded warmup (its output flushed and deleted before the timed run) plus one timed run, emitting `wall_times.csv` and per-dataset concordance reports; the wall times and F1 values under `manuscript/data/` are its outputs. Every comparator call scores against the Ensembl VEP release the ground truth's provenance names (`--reference-release` overrides). |
| `test_run_clone_measurement.sh` | Adversarial tests for `run_clone_measurement.sh`: stub engines inject faults to prove each failure-propagation gate fires, fork parity included                          |
| `run_concordance.sh`          | Development-time single-VCF harness: runs both Ensembl VEP (via Docker) and vep-rs on the same input, then compares outputs semantically under the release the Perl image names (`--reference-release` overrides). Independent of `run_clone_measurement.sh`. |
| `compare_vep_fields.py`       | Per-column and per-`Extra`-key agreement of two default-format outputs over their shared consequence keys (`fields.json`)                                                |
| `compare_vep_outputs.py`      | SNP/indel comparator: tuple matching on (Location, Allele, Feature, Feature_type, Consequence) with precision, recall, F1; `EXCLUSION_REGISTRY` carries the adjusted-F1 pair mask and `--reference-release {115.2,116.2}` (default 116.2) selects the one-sided classes of the Ensembl VEP release scored against ([docs/concordance-provenance/116.2-divergence-classes.md](../docs/concordance-provenance/116.2-divergence-classes.md)) |
| `reference_release.py`        | Derives `--reference-release` from a reference set's provenance (image tag, else cache version), an image reference or a checkout's `Constants.pm`; fails loudly when none names a release |
| `test_reference_release.py`   | Unit tests for the release derivation                                                                                                                                   |
| `test_compare_vep_outputs.py` | Unit tests for the SNP/indel comparator and its exclusion registry                                                                                                      |
| `compare_vep_formats.py`      | Per-column comparator for every output format: an Ensembl VEP output against a candidate's (default, tab, VCF, JSON, Parquet) with `--masks`, `--floors` and `--field-pairs`; `--self-parity` across one engine's formats; `--warnings` for the warnings and skipped-variants files |
| `test_compare_vep_formats.py` | Unit tests for the per-column comparator on the committed `tests/golden/116/GRCh38-release116` outputs                                                                  |
| `test_compare_vep_fields.py`  | Unit tests for the per-column comparator                                                                                                                                |
| `fastvep_normalize.py`        | Rewrites fastVEP's raw tab output so both engines' rows key alike before comparison                                                                                     |
| `count_warnings.py`           | Aggregates per-suite ERROR/WARN/INFO counts from both engines' stderr logs into `reports/warning_counts.json` and a Markdown section (`--out-md`)                       |
| `cache_parity_audit.pl`       | Audits the Storable-to-JSON cache conversion for field coverage over a transcript sample                                                                                |
| `compare_plugin_outputs.py`   | Plugin-specific output comparison with allele-aware matching                                                                                                            |
| `test_compare_plugin_outputs.py` | Unit tests for the plugin comparator                                                                                                                                |
| `trace_perl_intermediates.pl` | Trace Ensembl VEP intermediate values (codons, peptides, CDS bounds) for concordance debugging                                                                          |
| `compare_intermediates.py`    | Compare Perl and Rust intermediate values side-by-side for targeted discordant variants                                                                                 |
| `analyze_all_discordants.py`  | Classifies the `discordant.tsv` of several suites and builds a cross-suite pattern map; reads `<out>/suites/<id>/report(s)/discordant.tsv`, the file a `run_concordance.sh` run with `--work-dir <out>/suites/<id>` writes |
| `trace_perl_variant.sh`       | Run Ensembl VEP on a single variant with injected `warn()` tracing in 3 Perl modules. Produces per-transcript predicate trace. Requires `extract_perl_reference.sh` setup. |
| `extract_perl_reference.sh`   | Extracts five Ensembl VEP modules from the Docker image into `$VEP_PERL_REF_DIR` for debugging and tracing                                                              |
| `generate_reference_output.sh` | Writes the Ensembl VEP ground truth `run_clone_measurement.sh` scores against (one `output.txt` per SNP/indel suite, one `<vcf>.txt` per SV input) via Docker, with the harness's own invocation |
| `test_generate_reference_output.sh` | Tests for `generate_reference_output.sh`: a stub `docker` exercises the dry run and every fail-closed path                                                          |
| `prepare_benchmark_vcfs.py`   | Decompress, normalize, and optionally truncate VCF files for benchmarking                                                                                               |

### validation/

| Script                      | Description                                                                                                                                                                             |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `generate_sv_test_vcfs.py`  | Generate 13 synthetic VCF files with 41 variant categories (15,400 variants per assembly, seeded) for GRCh37/GRCh38                                                                     |
| `classify_discordants.py`   | Classify discordant tuples into root-cause categories (coordinate_only, consequence_swap, transcript_selection, true_missing/true_extra)                                                |
| `compare_sv_concordance.py` | Structural variant concordance analysis with per-type F1 scoring; `--reference-release {115.2,116.2}` (default 116.2, `MASKS_BY_RELEASE`) selects the masks of the Ensembl VEP release scored against, and under 116.2 the input records are held to the reference's own skip routes (`--reference-max-sv-size` when it ran with another cap) |
| `test_compare_sv_concordance.py` | Unit tests for the SV comparator and its masks                                                                                                                                     |
| `validate_test_vcfs.py`     | Validate all test VCFs: REF allele correctness, bcftools compliance, assembly consistency                                                                                               |
| `compare_gtf_tags_to_cache.py` | Compare the transcript attributes a GTF yields (`gencode_basic`, `gencode_primary`, `cds_start_NF`, `cds_end_NF`) with those a JSON cache carries, per code, listing every disagreement; exit 1 on any |
| `diff_caches.py`            | Structurally compare two JSON cache directories (transcripts, variations, info.json) with field-by-field diff and match-rate reporting                                                  |
| `capture_baseline.sh`       | Capture baseline annotation metrics from a JSON cache: runs vep-rs, computes per-consequence-type counts, frequency field stats, and cache checksums                                    |
| `compare_cache_outputs.sh`  | Full side-by-side comparison of baseline vs builder caches: runs vep-rs with each, diffs caches structurally, compares annotations, checks hard validation gates (F1, delta thresholds) |
| `validate_builder_cache.sh` | Cache builder validation: builds `vep-cache-builder`, generates a cache for a test region, and compares it against fixtures from an earlier run (`--fixture-mode`; exits 3 when none exist) or a Perl-derived baseline (`--baseline-cache`) |

### golden/

| Script                    | Purpose                                                                                                                                                                                                                                                                                                                                                                       |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `build_golden_corpus.py`  | `select`: from VEP reference outputs and their input VCFs, pick exemplar records for every observed consequence-set combination (several per combination, stratified by variant class and strand, narrow spans preferred, with `--cache` the cheapest in pruned-cache bytes first, every variant class kept, multi-allelic strata added; with `--keep-samples N` every record keeps `FORMAT` and the first N sample columns and every genotype class those samples show gets records; with `--include sid=name` a named record of suite `sid` is carried whatever the sampling chooses, and with `--keep-records` every record of the manifest being regenerated is kept) and write `variants.vcf` + `manifest.json`. `classify`: record which consequence keys differ between VEP and vep-rs, by `EXCLUSION_REGISTRY` class, with the corpus records they belong to; a record no VEP row names (one VEP dropped before annotation, its reason in the record's `reference_warning`) gets `reference_rows` 0 and the keys only vep-rs emits for it the class `reference_skipped_record`. `check`: fail when an observed combination has no exemplar, when a `covers` name is not a row of `--matrix`, when a `--family-columns` stratum has none, or, on a corpus built with `--keep-samples`, when `variants.vcf` lacks the manifest's sample columns or a genotype class the inputs show has no record. `covers`: list or extend the manifest's `covers` (or a flag set's), refusing a name that is not a matrix row. `flag-set`: write a `flag_sets` entry, a further run of the corpus under its own flags, input, stdin mode, formats, expected directory and `covers`, `--unordered` when Ensembl VEP emits a record's rows in hash order. |
| `prune_json_cache.py`     | Cut a JSON transcript cache down to the transcripts a corpus can reach (5 kb flank, breakend mates included), keeping a transcript in every shard it occupies; `--gzip`, `--drop <vefc key>` (`--drop-except <ids file>` names the transcripts that keep every dropped key), and `info.json` from a VEP cache's `info.txt` (`--perl-info`, or `--info-only`); `--variations <dir>` cuts the corpus's variation slice as bgzip with a `.csi` (`--variation-flank`, `--exclude-variation <regex>`, `--variation-info <info.json>` when the slice's columns are those of a cache other than the one the metadata describes); `CHECKSUMS.sha256` is written last, and `--manifest-only` writes it alone over an assembled cache directory. |
| `test_golden_scripts.py`  | Tests for both, plus the committed corpora's digests against their `provenance.json`.                                                                                                                                                                                                                                                                                       |

### adapters/

| Script                   | Purpose                                                                                                                              |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------ |
| `parquet_to_vep_tab.py`  | Project a vep-rs Parquet directory (flat or nested) back to the `--tab` columns; the Parquet round-trip test compares the result line for line. |

### data/

| Script                             | Description                                                                                                                                                                  |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `storable_to_json.pl`              | Convert an Ensembl VEP Storable cache to JSON for vep-rs: every chromosome directory, every field the cache dumper keeps, carried, flattened or dropped under a rule the parity audit names, `info.txt` to the `info.json` the output headers read (source versions, `contig_order`), and `chr_synonyms.txt` copied beside it; a Sereal cache is refused. Runs inside the `ensemblorg/ensembl-vep` Docker container. Arguments: `cache_dir output_dir [chr1,chr2,...]`. See usage below. |
| `setup_plugin_data.sh`             | Stage plugin data per assembly (`--assembly` REQUIRED). `--source upstream` (default) prints each database's URL and preparation steps for you to run (it downloads nothing); `--source s3` syncs an already-prepared copy from your own bucket (`--bucket` or `$VEP_BENCHMARK_BUCKET`; layout and `MANIFEST.json` schema in the script header). Skips databases the assembly does not publish. quick/full tiers |
| `plugin_data_layout.sh`            | The one table of staged plugin data file names, sourced by `setup_plugin_data.sh` and `concordance/run_concordance.sh`                                                        |
| `test_plugin_data_layout.sh`       | Proves every plugin data path `run_concordance.sh` passes is a file `setup_plugin_data.sh` stages, on both assemblies, and that both Perl `LoF` argument strings pin `run_splice_predictions:0` and the Perl output-cache key changes with the plugin argument strings |
| `download_real_world_vcfs.sh`      | Idempotent seeder for the 14 real-world source VCFs (SV + all-variants, both assemblies, chr21 except the genome-wide ClinVar files). Checks local, then an optional S3 copy, then the upstream host; copies back to S3 when `$VEP_BENCHMARK_S3` or `$VEP_BENCHMARK_BUCKET` is set. Does not produce the canonical, sites-only or multi-sample derivatives (see its header). |
| `compare_storable_cache_versions.pl` | Dump every transcript of a VEP Storable cache as `stable_id contig start end md5_full md5_content` (the content digest strips the run-specific fields), one line per transcript; runs inside `ensemblorg/ensembl-vep:release_116.2`. |
| `compare_storable_cache_versions.py` | Join two such dumps: ids in one cache only, identical, identical after normalisation, differing; writes the id lists and a summary JSON, and with `--parity-csv` the `storable_113_*` columns of `manuscript/data/transcript_set_parity.csv`. |
| `dump_storable_transcripts.pl`     | For a file of transcript ids, write one canonical dump per transcript from a VEP Storable cache (prediction matrices summarised as length and md5); runs inside `ensemblorg/ensembl-vep:release_116.2`. |
| `classify_storable_transcript_differences.py` | Read two such dumps for the same ids (optionally a vep-rs JSON cache, `.json` or `.json.gz` shards) and classify each transcript's difference (`peptide`, `predictions`, `protein_features`, `other`, `identical_after_normalisation`), writing `differing_<a>_<b>_classes.tsv` and a per-class summary JSON. |

#### storable_to_json.pl Usage

Converts Ensembl VEP Storable cache files to the JSON format expected by `json_cache.rs`. This is a **testing tool** for concordance validation: it produces a cache derived from the exact same transcript data Ensembl VEP uses, eliminating transcript set differences.

The script takes positional arguments only: `cache_dir output_dir [chr1,chr2,...]`, with the
optional chromosome list as ONE comma-separated argument. It lives at `scripts/data/`, so
mount that directory into the container.

```bash
# Convert specific chromosomes (one comma-separated list)
docker run --rm -v /path/to/cache:/input:ro -v /output/dir:/output \
  -v $(pwd)/scripts/data:/scripts:ro \
  ensemblorg/ensembl-vep:release_116.2 \
  perl /scripts/storable_to_json.pl /input /output 21,22

# Convert all chromosomes
docker run --rm -v /path/to/cache:/input:ro -v /output/dir:/output \
  -v $(pwd)/scripts/data:/scripts:ro \
  ensemblorg/ensembl-vep:release_116.2 \
  perl /scripts/storable_to_json.pl /input /output
```

To inspect a Storable object's raw structure, use Perl's `Data::Dumper` inside the same
container; the converter has no dump mode.

Output: `{output}/transcripts/{chr}/{start}-{end}.json` for every chromosome directory of the
cache (the assembled chromosomes and the patch, alternate-locus, scaffold and `LRG_`
directories alike; the optional list restricts it), `{output}/info.json` (the `info.txt`
metadata as the output headers read it, plus `contig_order`: the converted contigs as autosomes by
number, then X, Y, MT, then every other contig by name; the `sift` and `polyphen` markers only when
a matrix of that tool was written) and `{output}/chr_synonyms.txt`, the cache's own copy. A cache
whose `info.txt` names the `sereal` serialiser, or whose chromosome directories hold `.sereal`
shards, is refused with `ERROR: [storable_to_json] ...` and exit status 1, as is a cache with no
`<start>-<end>.gz` shard at all, a listed chromosome that is not a directory of the cache, an
unreadable chromosome directory, and a run in which any shard failed to read or write (each such
shard is named on a `WARN:` line; the run stops before `chr_synonyms.txt` and `info.json` are
written, so a partial output lacks both). Each transcript carries every field the cache dumper
keeps, carried, flattened or dropped under a rule `cache_parity_audit.pl` names (dropped: the Perl
object references, the per-exon and per-intron strand, which is the transcript's, a sequence edit's
description and the dump-time flags): the identifiers
(`uniprot_isoform` included), the protein features with the analysis label Ensembl VEP's `DOMAINS`
column prints, as `analysis`, the SIFT and PolyPhen matrices (`predictions_data`, the Storable's gzip blob
base64-encoded, beside `translation_md5`), and in `variation_effect_feature_cache` the
`seq_edits`, the translation's `Bio::EnsEMBL::SeqEdit` list (`code`, `name`, `start`, `end`,
`alt_seq`, `length_diff`, protein coordinates; `[]` when there are none), which vep-rs applies to
the reference peptide. A cache without a `seq_edits` key (the published conversions among them)
loads; vep-rs applies no edit and reads the reference peptide as cached, so convert again to get
VEP's calls at alternative initiation codons.

### release/

| Script                          | Description                                                                                                                                              |
| ------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `render_release_notes.py`       | Renders the release page from the repository's own files, every header `##`, never more than 400 lines: the opener and the optional Highlights of `.github/release-notes/vX.Y.Z.md` (the only hand-written lines), then either its `What changed, by area` section (one bullet per area with the changes nested under it, each breaking one beginning `Breaking:` and first in its area) with its optional `Security and robustness`, `Dependencies and tooling` and `Removed` sections, or, for a notes file without one, its optional Upgrade notes and the version's `CHANGELOG.md` section with its kind headers at `##` and one line per bullet, the whole-genome table (the six datasets and `All six`) and the structural-variant F1 table from the two provenance records under `docs/concordance-provenance/`, the Install bullets from `SHA256SUMS`, the glibc floors `check_glibc_floor.py --json` measured, the build and the workspace `rust-version`, the Cite line from `CITATION.cff`, and the full-changelog footer. `.github/workflows/release.yml` runs it on every tag and on a `workflow_dispatch` dry run. |
| `test_render_release_notes.py`  | Tests on the 0.3.2 inputs copied into a scratch repository and, for the six-dataset shape, on release and population records the tests build under the same version: the whole page against its expected text, each section's rule (Highlights bounds, optional Upgrade notes, the page by area and its refusals, kind headers, one line per bullet, the two tables, the Install bullets, the Cite line, the footer, the line cap) and the failure modes |
| `check_glibc_floor.py`          | Reads the glibc version each Linux binary requires from its ELF `.gnu.version_r` section (a binary or every ELF member of a `.tar.gz`), one line per binary; `--floor X.Y` fails the run naming the symbols a binary needs above the floor, and `--json FILE` writes the floor measured per archive for `render_release_notes.py --glibc-floors`. `.github/workflows/release.yml` runs it on every Linux binary of every archive. |
| `test_check_glibc_floor.py`     | Tests on ELF files the tests assemble (both classes and byte orders): the highest version per binary, the floor's failure line and its symbols, archives, the JSON floors, and the refused inputs. |

### parity/

| Script                          | Description                                                                                                                                              |
| ------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `check_parity_matrix.py`        | Holds `tests/parity/matrix.tsv` to the tree (CI job `Comparator Tests`): a flag in `crates/vep-cli/src/args.rs` without a row, a landed row without its tests, or a named test or corpus that does not resolve fails. `crates/vep-cli/tests/parity_matrix.rs` asks clap for the flag list and is the arbiter. |
| `test_check_parity_matrix.py`   | Tests on a fixture tree: one per rule, every clap field shape the reader follows, the test-attribute rule, the header and vocabulary checks, the exit codes |

### docs/

| Script                  | Description                                                                                                                                       |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `build_release_site.sh` | Builds one channel of the site and checks it. `build_release_site.sh v0.3.2 target/site` builds a release: the tag's README, CHANGELOG, CONTRIBUTING and `docs/` with this checkout's `book.toml`, `theme/` and `scripts/docs/`, the frame pages and redirects the tag lacks dropped, then the link check against the tag's tree and `site_postbuild.py`. `build_release_site.sh --main target/site` builds the development preview of HEAD, titled and addressed as the `main/` channel. The `Pages` workflow publishes both |
| `mdbook_repo_links.py`  | mdBook preprocessor (`book.toml`). Makes `CHANGELOG.md` the changelog page and `CONTRIBUTING.md` the contributing page, and `README.md` the home page when `VEP_RS_DOCS_README_HOME=1`, or when that variable is unset and `docs/index.md` is the three-line stub. Pastes a README section where a page asks for one by heading, fills `{{vep-rs-ref}}` and `{{vep-rs-release}}`, and sends links that leave `docs/` to their GitHub view at the ref built (`$VEP_RS_DOCS_REF`, default `main`) |
| `check_site_links.py`   | Fails on a link in the built site to a missing page, a missing heading id, or a repository path that does not exist in the tree built from (CI job `Pages`) |
| `site_postbuild.py`     | Finishes one built channel after `mdbook build`: a canonical link and `og:url` on every page, `versions.json` for the header's channel switcher, `sitemap.xml` for a release, the preview's channel and `noindex` meta tags. Fails on a redirect whose destination page or id the site lacks, and on a benchmark chart whose `data-release` differs from the README's `## Current release` heading |
| `site_metrics.py`       | Prose metrics of the pages `docs/SUMMARY.md` lists, as a Markdown table or `--json`: words, headings per 300 words, paragraphs over 120 words, sentences over 45 words, tables over 4 columns, code lines over 80 characters. `--gate` exits 1 on a page, paragraph, sentence, table or code line over its limit outside the script's allow-list |
| `site_screenshots.py`   | Screenshots every page of a built site at 14 widths from 320 to 2560 px in the light and dark themes. Probes each view for horizontal overflow, wrapped inline code, font size, line length and clipped code or tables, runs axe when `axe-playwright-python` is importable, and writes a contact sheet per page group. The one script here with a third-party runtime, Playwright with Chromium. Run by hand, never in CI |
| `render_benchmark_chart.py` | Renders `docs/assets/benchmark.svg`, the home page's wall-time chart, from the release records under `docs/concordance-provenance/`. One horizontal bar per dataset and engine (Ensembl VEP, fastVEP, vep-rs) on a logarithmic axis, every bar labelled, colours through CSS classes the theme fills. Stamps the record's release as `data-release`, which `site_postbuild.py` holds to the README. Run at each release |
| `test_site_tools.py`    | Tests on fixtures the tests write: the preprocessor's home-page rule, README sections, tokens and link rewriting through its stdin protocol; the metrics' section selection, per-page measures and `--gate` exit codes; the post-build head tags, `versions.json`, `sitemap.xml`, redirect and chart-release checks; the chart renderer's record selection, failure modes and the committed `benchmark.svg` byte for byte; the screenshot script's summary helpers and `--summary` entry point without Playwright (CI job `Comparator Tests`) |
