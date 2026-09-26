# Scripts

Utility scripts for testing, benchmarking, validation, and data preparation.

## Environment

These scripts read their data locations from environment variables, each with a repo-relative default under `./.vep/` so nothing depends on a specific machine layout. Override any of them to point at your own directories, or pass the equivalent `--flag` where a script exposes one. The optional S3 copy is off unless `VEP_BENCHMARK_S3` or `VEP_BENCHMARK_BUCKET` is set.

| Variable                  | Default                 | Used for                                                     |
| ------------------------- | ----------------------- | ------------------------------------------------------------ |
| `VEP_VCF_DIR`             | `./.vep/vcf`            | Downloaded input VCFs, laid out `<asm>/{sv,all_variants}/` by `data/download_real_world_vcfs.sh`; the canonical input location |
| `VEP_BENCHMARK_DIR`       | `./.vep/benchmark`      | `--benchmark-dir` of `concordance/run_concordance.sh` (a flat directory of VCFs); point it at `./.vep/vcf/<asm>/all_variants` |
| `VEP_PERL_CACHE_DIR`      | `./.vep/cache`          | Perl VEP cache (mounted into the Docker container)           |
| `VEP_PERL_REF_DIR`        | `./.vep/perl_reference` | Perl VEP modules extracted from Docker for tracing           |
| `VEP_WORK_DIR`            | `./.vep/work`           | Scratch output: concordance runs and Perl trace output       |
| `VEP_BENCHMARK_S3`        | (unset)                 | Optional S3 URI holding a copy of the VCFs `data/download_real_world_vcfs.sh` fetches |
| `VEP_BENCHMARK_BUCKET`    | (unset)                 | Bucket name for `data/setup_plugin_data.sh --source s3`; `data/download_real_world_vcfs.sh` accepts it too |

## Directory Layout

| Directory      | Purpose                                                                       | Prerequisites                 |
| -------------- | ----------------------------------------------------------------------------- | ----------------------------- |
| `concordance/` | Compare vep-rs output against Perl VEP                                        | Docker, Perl VEP cache        |
| `validation/`  | Standalone test data generation, SV concordance, and cache builder validation | Python 3.10+                  |
| `data/`        | Download input VCFs and stage plugin annotation data                          | Internet access, tabix        |
| `golden/`      | Build the per-release, per-assembly golden corpora under `tests/golden/`      | Python 3.10+; Ensembl VEP outputs and a JSON cache |
| `adapters/`    | Convert vep-rs outputs between formats                                        | `duckdb` CLI                  |

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

2. **Generate the reference output.** `concordance/generate_reference_output.sh` runs Perl VEP release 115.2 over the canonical inputs and writes the SNP/indel reference output per suite and the structural-variant reference output per input VCF into the `ground_truth/` layout above, passing `--fasta` at every call site; each of the 16 structural-variant VCFs per assembly is annotated in its own invocation, so no file's transcript-cache population can affect another's. The published values were scored against reference output at cache version 115, except ClinVar GRCh37, whose reference is the cache 113 set that a cache 115 regeneration reproduces key for key (`docs/concordance-provenance/2026-09-18-grch37-cache-113-vs-115.json`).

3. **Convert the cache both engines read.** vep-rs annotates from a JSON conversion of the Storable cache Perl VEP read (`data/storable_to_json.pl`, usage below), which holds the transcript set fixed across the two engines. `concordance/cache_parity_audit.pl` audits the conversion: it reads a seeded sample of transcripts from both caches and classifies every Storable field against the converter's rules (identical, renamed, converted, dropped because undefined, dropped by a named rule, or an unexplained gap). Its manifest for the GRCh37 pair, the cache-version-113 Storable VEP release 115 ships for that assembly and its JSON conversion, is `docs/concordance-provenance/2026-09-16-cache-parity-audit.json` (the manifest names neither release nor assembly; its `json_transcripts_indexed` of 195,232 is the GRCh37 cache's transcript count). Every field vep-rs reads for consequence calling is carried. Two Storable fields are absent from the JSON side under no rule (`_gene_version`, which the converter maps and the published cache lacks, and `_uniprot_isoform`, unmapped), and vep-rs reads neither. The SIFT and PolyPhen prediction matrices are dropped by a named rule: the converter reads a `predictions_data` field that only the native cache builder writes, and the concordance tuple carries no prediction field, so this enters no published F1. The per-transcript `seq_edits` are dropped by a named rule as well, but their effect is carried: the cached peptide is VEP's own translation with those edits applied, and vep-rs reads the edited residues from it.

4. **Run the measurement.**

   ```bash
   scripts/concordance/run_clone_measurement.sh \
     --engine vep-rs --vep-binary target/release/vep \
     --data-dir /path/to/data --output-dir /path/to/data/work
   ```

   Each SNP/indel suite is scored by `concordance/compare_vep_outputs.py` and each structural-variant set by `validation/compare_sv_concordance.py`; the per-suite `report/` directories carry the raw and adjusted tuple counts, and `wall_times.csv` the timing. The published cells ran with `--fork 16` on 32-vCPU hosts (the host specification is recorded in `docs/concordance-provenance/2026-09-20-n20-parallel-clones.json`); `--fork` is overridable and F1 does not depend on it, which `crates/vep-cli/src/runner.rs` gates with a fork-1-versus-fork-4 equality test.

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
| `run_clone_measurement.sh`    | Measurement entrypoint: runs every benchmark dataset for one engine (`--engine vep-rs\|fastvep\|perl`) against a populated data tree, one discarded warmup (its output flushed and deleted before the timed run) plus one timed run, emitting `wall_times.csv` and per-dataset concordance reports; the wall times and F1 values under `manuscript/data/` are its outputs. |
| `test_run_clone_measurement.sh` | Adversarial tests for `run_clone_measurement.sh`: stub engines inject faults to prove each failure-propagation gate fires, fork parity included                          |
| `run_concordance.sh`          | Development-time single-VCF harness: runs both Perl VEP (via Docker) and vep-rs on the same input, then compares outputs semantically. Independent of `run_clone_measurement.sh`. |
| `compare_vep_fields.py`       | Per-column and per-`Extra`-key agreement of two default-format outputs over their shared consequence keys (`fields.json`)                                                |
| `compare_vep_outputs.py`      | SNP/indel comparator: tuple matching on (Location, Allele, Feature, Feature_type, Consequence) with precision, recall, F1; `EXCLUSION_REGISTRY` carries the adjusted-F1 mask |
| `test_compare_vep_outputs.py` | Unit tests for the SNP/indel comparator and its exclusion registry                                                                                                      |
| `compare_output_formats.py`   | Output-format self-parity: the same engine must emit the same consequence calls in VCF, tab, JSON and Parquet                                                          |
| `test_compare_output_formats.py` | Unit tests for the format-parity harness                                                                                                                            |
| `test_compare_vep_fields.py`  | Unit tests for the per-column comparator                                                                                                                                |
| `fastvep_normalize.py`        | Rewrites fastVEP's raw tab output so both engines' rows key alike before comparison                                                                                     |
| `count_warnings.py`           | Aggregates per-suite ERROR/WARN/INFO counts from both engines' stderr logs into `reports/warning_counts.json` and a Markdown section (`--out-md`)                       |
| `cache_parity_audit.pl`       | Audits the Storable-to-JSON cache conversion for field coverage over a transcript sample                                                                                |
| `compare_plugin_outputs.py`   | Plugin-specific output comparison with allele-aware matching                                                                                                            |
| `test_compare_plugin_outputs.py` | Unit tests for the plugin comparator                                                                                                                                |
| `trace_perl_intermediates.pl` | Trace Perl VEP intermediate values (codons, peptides, CDS bounds) for concordance debugging                                                                             |
| `compare_intermediates.py`    | Compare Perl and Rust intermediate values side-by-side for targeted discordant variants                                                                                 |
| `analyze_all_discordants.py`  | Classifies the `discordant.tsv` of several suites and builds a cross-suite pattern map; reads `<out>/suites/<id>/report(s)/discordant.tsv`, the file a `run_concordance.sh` run with `--work-dir <out>/suites/<id>` writes |
| `trace_perl_variant.sh`       | Run Perl VEP on a single variant with injected `warn()` tracing in 3 Perl modules. Produces per-transcript predicate trace. Requires `extract_perl_reference.sh` setup. |
| `extract_perl_reference.sh`   | Extracts five Perl VEP modules from the Docker image into `$VEP_PERL_REF_DIR` for debugging and tracing                                                                 |
| `generate_reference_output.sh` | Writes the Perl VEP ground truth `run_clone_measurement.sh` scores against (one `output.txt` per SNP/indel suite, one `<vcf>.txt` per SV input) via Docker, with the harness's own invocation |
| `test_generate_reference_output.sh` | Tests for `generate_reference_output.sh`: a stub `docker` exercises the dry run and every fail-closed path                                                          |
| `prepare_benchmark_vcfs.py`   | Decompress, normalize, and optionally truncate VCF files for benchmarking                                                                                               |

### validation/

| Script                      | Description                                                                                                                                                                             |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `generate_sv_test_vcfs.py`  | Generate 13 synthetic VCF files with 41 variant categories (15,400 variants per assembly, seeded) for GRCh37/GRCh38                                                                     |
| `classify_discordants.py`   | Classify discordant tuples into root-cause categories (coordinate_only, consequence_swap, transcript_selection, true_missing/true_extra)                                                |
| `compare_sv_concordance.py` | Structural variant concordance analysis with per-type F1 scoring                                                                                                                        |
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
| `build_golden_corpus.py`  | `select`: from VEP reference outputs and their input VCFs, pick exemplar records for every observed consequence-set combination (several per combination, stratified by variant class and strand, narrow spans preferred, multi-allelic strata added) and write `variants.vcf` + `manifest.json`. `classify`: record which consequence keys differ between VEP and vep-rs, by `EXCLUSION_REGISTRY` class, with the corpus records they belong to. `check`: fail when an observed combination has no exemplar. |
| `prune_json_cache.py`     | Cut a JSON transcript cache down to the transcripts a corpus can reach (5 kb flank, breakend mates included), keeping a transcript in every shard it occupies; `--gzip`, `--drop <vefc key>`, and `info.json` from a VEP cache's `info.txt` (`--perl-info`, or `--info-only`).                                                                                                  |
| `test_golden_scripts.py`  | Tests for both, plus the committed corpora's digests against their `provenance.json`.                                                                                                                                                                                                                                                                                       |

### adapters/

| Script                   | Purpose                                                                                                                              |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------ |
| `parquet_to_vep_tab.py`  | Project a vep-rs Parquet directory (flat or nested) back to the `--tab` columns; the Parquet round-trip test compares the result line for line. |

### data/

| Script                             | Description                                                                                                                                                                  |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `storable_to_json.pl`              | Convert Perl VEP Storable cache to JSON for vep-rs, and its `info.txt` to the `info.json` the output headers read (source versions). Runs inside the `ensemblorg/ensembl-vep` Docker container. Arguments: `cache_dir output_dir [chr1,chr2,...]`. See usage below. |
| `setup_plugin_data.sh`             | Stage plugin data per assembly (`--assembly` REQUIRED). `--source upstream` (default) prints each database's URL and preparation steps for you to run (it downloads nothing); `--source s3` syncs an already-prepared copy from your own bucket (`--bucket` or `$VEP_BENCHMARK_BUCKET`; layout and `MANIFEST.json` schema in the script header). Skips databases the assembly does not publish. quick/full tiers |
| `plugin_data_layout.sh`            | The one table of staged plugin data file names, sourced by `setup_plugin_data.sh` and `concordance/run_concordance.sh`                                                        |
| `test_plugin_data_layout.sh`       | Proves every plugin data path `run_concordance.sh` passes is a file `setup_plugin_data.sh` stages, on both assemblies                                                        |
| `download_real_world_vcfs.sh`      | Idempotent seeder for the 14 real-world source VCFs (SV + all-variants, both assemblies, chr21 except the genome-wide ClinVar files). Checks local, then an optional S3 copy, then the upstream host; copies back to S3 when `$VEP_BENCHMARK_S3` or `$VEP_BENCHMARK_BUCKET` is set. Does not produce the canonical, sites-only or multi-sample derivatives (see its header). |
| `compare_storable_cache_versions.pl` | Dump every transcript of a VEP Storable cache as `stable_id contig start end md5_full md5_content` (the content digest strips the run-specific fields), one line per transcript; runs inside `ensemblorg/ensembl-vep:release_115.2`. |
| `compare_storable_cache_versions.py` | Join two such dumps: ids in one cache only, identical, identical after normalisation, differing; writes the id lists and a summary JSON, and with `--parity-csv` the `storable_113_*` columns of `manuscript/data/transcript_set_parity.csv`. |
| `dump_storable_transcripts.pl`     | For a file of transcript ids, write one canonical dump per transcript from a VEP Storable cache (prediction matrices summarised as length and md5); runs inside `ensemblorg/ensembl-vep:release_115.2`. |
| `classify_storable_transcript_differences.py` | Read two such dumps for the same ids (optionally a vep-rs JSON cache, `.json` or `.json.gz` shards) and classify each transcript's difference (`peptide`, `predictions`, `protein_features`, `other`, `identical_after_normalisation`), writing `differing_<a>_<b>_classes.tsv` and a per-class summary JSON. |

#### storable_to_json.pl Usage

Converts Perl VEP Storable cache files to the JSON format expected by `json_cache.rs`. This is a **testing tool** for concordance validation: it produces a cache derived from the exact same transcript data Perl VEP uses, eliminating transcript set differences.

The script takes positional arguments only: `cache_dir output_dir [chr1,chr2,...]`, with the
optional chromosome list as ONE comma-separated argument. It lives at `scripts/data/`, so
mount that directory into the container.

```bash
# Convert specific chromosomes (one comma-separated list)
docker run --rm -v /path/to/cache:/input:ro -v /output/dir:/output \
  -v $(pwd)/scripts/data:/scripts:ro \
  ensemblorg/ensembl-vep:release_115.2 \
  perl /scripts/storable_to_json.pl /input /output 21,22

# Convert all chromosomes
docker run --rm -v /path/to/cache:/input:ro -v /output/dir:/output \
  -v $(pwd)/scripts/data:/scripts:ro \
  ensemblorg/ensembl-vep:release_115.2 \
  perl /scripts/storable_to_json.pl /input /output
```

To inspect a Storable object's raw structure, use Perl's `Data::Dumper` inside the same
container; the converter has no dump mode.

Output: `{output}/transcripts/{chr}/{start}-{end}.json` + `{output}/info.json`
