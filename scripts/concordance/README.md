# Concordance Testing

These scripts compare vep-rs output against the reference Ensembl VEP implementation to verify correctness. They are **optional**: you do not need Ensembl VEP to use vep-rs.

## Prerequisites

- **Docker** (for running Ensembl VEP via `ensemblorg/ensembl-vep` image)
- **Ensembl VEP cache** (downloaded via the Ensembl VEP installer or Docker)
- **JSON cache** (built via `vep-cache-builder`, or provide pre-existing with `--json-cache-dir`)
- **Benchmark VCF** (e.g., ClinVar)
- **Python 3.10+** (for the helper scripts the launchers run)
- **A release build of `tools/vep-formats-compare`**, the format comparator `feature_oracle.sh`
  starts (`cargo build --release --locked --manifest-path tools/vep-formats-compare/Cargo.toml`)

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
`tools/vep-sv-compare`, is a Rust tool with a Cargo root of its own and is
documented in `tests/sv_validation/README.md`.

## Scripts

| Script                       | Purpose                                                                                                                                                                                                                                                                                |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `run_clone_measurement.sh`   | Measurement entrypoint: annotates every benchmark dataset for one engine (`--engine {vep-rs,fastvep,perl}`) against a populated data tree, with a discarded cold warmup (its output flushed and deleted before the timed run) plus one timed run, emitting `wall_times.csv` and per-dataset concordance reports; the wall times and F1 values under `manuscript/data/` are its outputs. Each comparator call carries `--reference-release`, read per cell from the ground truth's provenance (or `--reference-release`); a ground truth naming no release fails the cell. `--fastvep-extra-args "<flags>"` (fastVEP only) appends words, after the harness's own arguments, to every fastVEP `annotate` run, the warmup and the timed run of each SNP/indel cell and each input of the SV per-VCF loop, through `<output-dir>/fastvep-wrapper.sh`: each word is an option (`-x`, `--name`, `--name=value`) or a value of `[A-Za-z0-9_=.,:/-]` that opens with no dash, and any other word stops the run before anything is written; a run given words writes `<output-dir>/provenance.json` (the binary, `fastvep_extra_args` as the words joined by single spaces, `fastvep_extra_argv`, `fastvep_wrapper`, `fastvep_executed_argv`), and a run given none writes no wrapper and no record. A fastVEP cell whose transcript cache is not a file with bytes is refused before fastVEP starts, with one `ERROR` line and a row of run kind `fastvep_cache_missing` and exit code 1: at a `--transcript-cache` path that names no file fastVEP parses the `--gff3` file, saves a cache of its own there and exits 0, so the cell would be timed and scored on the GFF3's gene models, while beside a cache file that loads fastVEP does not parse the GFF3. Words that carry `--sa-dir` take cells of one assembly, one invocation per assembly: see [fastVEP with its annotation databases](#fastvep-with-its-annotation-databases). |
| `test_run_clone_measurement.sh` | Adversarial tests for `run_clone_measurement.sh`: stub engines inject faults to prove each failure-propagation gate fires, fork parity included; the FASTVEP-EXTRA-ARGS case proves that `--fastvep-extra-args` reaches every fastVEP `annotate` run and the record, each refusal (another engine, a word outside the rule, an output directory whose path holds `=` and, without the words, a binary whose path does, `--sa-dir` over cells of two assemblies), and that a run without the words writes no wrapper and no record; the FASTVEP-CACHE case proves that a fastVEP cell whose transcript cache path names no file, an empty file or a link to no file is refused before fastVEP starts, with a row of run kind `fastvep_cache_missing`, and that the same tree with its cache file runs every cell on that file |
| `run_concordance.sh`         | Development-time single-VCF orchestrator: runs Ensembl VEP (Docker) + vep-rs on one input, with plugin testing and perf telemetry. It does not wrap `run_clone_measurement.sh`. The comparator's `--reference-release` is derived from the Ensembl VEP image tag, the fixtures directory's provenance or the local checkout's `Constants.pm`, or given with `--reference-release`.    |
| `test_run_concordance.sh`    | Tests for `run_concordance.sh` on the script's own text, each block extracted and run with fabricated inputs: the canonical-contigs flag is bare under GRCh37 and GRCh38, carries the `contig_order` of the JSON cache's `info.json` under any other assembly (a cache that lists none is refused) and is `--no-canonical-contigs` when the filter is off; the messages for a missing `--json-cache-dir` and a missing `--fasta` name the human primary-assembly FASTA under homo_sapiens and the species' own toplevel FASTA under any other |
| `../../tools/vep-outputs-compare/` | SNP/indel comparator, a Rust tool that scores in three modes across the host's cores: extracts `(Location, Allele, Feature, Feature_type, Consequence)` tuples from both outputs and computes precision, recall, F1; its `EXCLUSION_REGISTRY` is the adjusted-F1 pair mask under both reference releases, the 116.2 registry adds the pair rules that read the pair's two rows in the outputs (`row_pair_rules`), two of them with the record's alleles from the first column or the input records (`--input-vcf`, `--input-dir`), and `--reference-release {115.2,116.2,116.3}` (default 116.3) selects the one-sided classes: the paper's scope class under 115.2, the fourteen of the linked page under 116.2, which are the record Ensembl VEP 116 skips, the three shapes of its chromosome-aware overlap and ten further record shapes ([docs/concordance-provenance/116.2-divergence-classes.md](../../docs/concordance-provenance/116.2-divergence-classes.md)). `cargo test` there holds the binary to the goldens over its fixtures. `run_concordance.sh`, `run_clone_measurement.sh` and `../validation/compare_cache_outputs.sh` score with this tree's build of the tool and with no other program. See its [README](../../tools/vep-outputs-compare/README.md) |
| `reference_release.py`       | Derives a comparator's `--reference-release` from a reference set's provenance (the Ensembl VEP image tag, else the cache version), an image reference or a checkout's `Constants.pm`; fails loudly when none names a release                                                             |
| `test_reference_release.py`  | Unit tests for the release derivation                                                                                                                                                                                                                                                  |
| `../../tools/vep-fields-compare/` | Per-column and per-`Extra`-key agreement over the shared consequence keys of two default-format outputs; writes `fields.json` and reports duplicate keys on each side: a Rust tool, run by `run_clone_measurement.sh` from the tree's release build |
| `../../tools/vep-plugin-pairs/` | A Rust tool with two commands, each on every core. `score` is the plugin tuple scorer: one tuple per populated plugin field of an output row, `(Location, Allele, Feature, field, value)`, each field resolved by the tool's plugin registry among the plugins the scored run loaded (`--plugins`, by registry name or alias) or by a `--custom-source` prefix, and every other Extra key outside Ensembl VEP's own listed as unclaimed, as is a field two of the run's plugins write where no custom source claims it; `--run-command` hands the registry the run's command line, which refuses a run whose flags make a loaded plugin's column name one of Ensembl VEP's own keys as well; scores each `--engine` (a default-format output, or the five-column tuple file a comparator's normaliser writes) against one Ensembl VEP default-format output and writes per source, engine and projection (with and without the Feature) the integer counts a raw and an adjusted F1 are derived from, `text_only_differences`, `residual_rows` and `shared_keys`, rows that share a (Location, Allele, Feature) never collapsed: a tuple counts as often as its side's rows write it; `--class-table <tsv>`, `plugin_class_table.tsv` where a run names none, declares the documented classes with the plugins and fields each belongs to, `--classes` names the ones applied, every one but the rows that read `not met` where a run names none, and each is applied to its own plugin's fields alone through the tool's class rules, its count per entry under `masked`; `--gwas-catalog <file>` hands the classes of GWAS the associations file the run's plugin read, without which they set nothing aside; `--utrannotator-sweep <path>` hands the classes of UTRAnnotator Ensembl VEP's outputs of the input under fixed Perl hash seeds, one file a run or a directory of them, without which `utrannotator_consequence_precedence` sets nothing aside; `--json-cache <dir>` hands the classes of Downstream and of Blosum62 the candidate cache whose sequences they recompute from, the class of MaxEntScan the cache's introns, and `spliceai_multi_entry_cell` its transcripts' symbols, from which the class takes the gene of a row that writes no `SYMBOL`: without the cache the classes of Downstream, of Blosum62 and of MaxEntScan set nothing aside, and `spliceai_multi_entry_cell` sets aside no such row; `--spliceai-table <file>`, once per score file, hands the classes of SpliceAI the files the run's plugin read, without which `spliceai_multi_entry_cell` sets nothing aside, and the summary's last lines then count the rows `spliceai_multi_entry_cell` was asked about and refused, by the ground of the refusal; a rule that reads a record's other rows is handed them where the `Uploaded_variation` spells the record's alleles, and under an identifier where `--input <vcf>`, the run's input, names its record, the one line that carries the identifier or of several the one that writes the row at its Location with its allele: the reference's rows of the record, and for a sibling-allele class a default-format engine's rows of it as well; `plugin_value_follows_splice_term_past_the_record_span` and `plugin_value_follows_stop_codon_term_pair` take from the same input the reference bases the allele of a row under an identifier stands in for, read by position from the lines that carry the identifier, and read the row without its record where no input is given; `plugin_value_on_reference_skipped_record` is decided without a reference row, from both engines' rows at a value's Location and the first column of the reference, by the SNP/indel comparator's rule; `--residuals <dir>` writes the tuples outside the intersection that no applied class set aside, one file per source, engine and projection. `compare` is the per-field plugin comparison `run_concordance.sh` runs under `--plugins`: per plugin and file pair, the records both engines annotate compared field by field, a record on one side alone and an engine that annotated nothing each a failure, with a JSON report. `cargo test` there holds both commands. See its [README](../../tools/vep-plugin-pairs/README.md) |
| `plugin_class_table.tsv`     | The class table `vep-plugin-pairs score` reads unless `--class-table` names another: one row per documented class of plugin fields, which declares the class, the plugins it belongs to and the field names it reads on each, then `met` or `not met` and the columns `commit`, `genome build` and `note`, which a report carries, and a row that reads `not met` is declared and not applied unless `--classes` names it; the tool's tests hold every row to its class rules and its reader |
| `plugin_gate_runs.toml`      | The runs a release scores its supported plugins on, as data: five whole runs on the genome-wide ClinVar datasets for the plugins verified on real data, and one run per set of every plugin corpus, each naming the derived rows it scores, the words of both engines, the image of Ensembl VEP's reference and the scorer's words, with every place outside the tree a named input |
| `reproduce_plugin_concordance.py` | The driver of `plugin_gate_runs.toml`: for a run, a corpus or a derived row it checks the data files by name, takes Ensembl VEP's output by path or writes it from the image `generate_reference_output.sh` pins, runs `vep`, scores with `vep-plugin-pairs score` (the tree's build, or the binary `--pairing-tool` names) and prints each row's tuple counts and raw and adjusted F1 as plugin record entries; `--dry-run` prints every command. [Reproducing the plugin concordance](../../docs/plugins/reproducing-concordance.md) is its runbook. |
| `test_plugin_gate_runs.py`, `test_reproduce_plugin_concordance.py` | Tests for the data file against the derived rows, the corpora and the data layout, for the driver's dry runs, refusals and printed rows, and for the workflow step that counts the driver's tests that score |
| `../../tools/vep-formats-compare/` | Per-column comparator for every output format: one Ensembl VEP output against a candidate's in the default, tab, VCF, JSON or Parquet format (a Parquet actual, a file or a directory of parts, is projected to tab rows and compared with a tab expected as a multiset), every row, CSQ entry and JSON object, every field verbatim. Intended-divergence `--masks` come from a corpus manifest plus the built-in classes: `regulatory_feature_rows` sets the reference's RegulatoryFeature and MotifFeature rows aside, `order_only` lets the set-valued columns `DOMAINS`, `FREQS`, `ZYG` and `NearestExonJB` agree as multisets, the `5UTR_annotation` groups by their pairs and the `CLIN_SIG` lists of a row's known variants as a multiset of lists, `individual_allele_order` compares an individual's alleles as a set, `conservation_short_name_spelling` reads the reference's `phlyoP<N>way` key as `phyloP<N>way`, `pc_description_text` reads its `PC` header line with the candidate's description of the share the field reports, `motif_score_reverse_strand_mnv` accepts a `MOTIF_SCORE_CHANGE` the candidate writes on a reverse-strand motif row whose allele is two or more bases, the registry pair classes `start_cooccurrence_swap` and `splice_family_swap` set aside the pairs the SNP/indel comparator's excluding rules admit, and one field or warnings class per corrected behaviour of `docs/intended-divergences.md` accepts that behaviour's one shape of difference, every such agreement or rewrite counted in the report under `masks.applied`. `--floors` names a default-format column and a minimum populated count on both sides; the JSON format counts the column under its JSON key or, for a plugin's own column, inside the plugin object that nests it, and a column it has no row count for (a record-level column, `Feature_type`, `Extra`, a custom source's columns) is reported there as not applicable, as is a tab or VCF column neither header declares. `--field-pairs` accepts named value pairs, and `--nearest-sites` settles `NEAREST` by transcription start site through the candidate's JSON cache. `--input` keys the VCF lines and JSON objects by input record through a line's first five columns, since Ensembl VEP writes a structural variant's empty or `0` QUAL, FILTER or INFO back as `.`, and a default or tab comparison reads the input records where a class reads more than a row carries: under the three structural-variant classes (`reference_skipped_record`, which names a record from them and refuses a run without them, `breakend_mate_local_read` and `bnd_own_chromosome_rows_lost`), under a class that reads a record's alleles, which finds the one line that carries a row's identifier, and under `splice_term_past_the_record_span` and `stop_codon_overlap_insertion_range` with their plugin value classes, `stop_codon_completed_by_polyadenylation`, `stop_codon_overlap_without_stop_codon`, `stop_codon_behind_last_residue_without_stop_codon`, `start_codon_overlap_without_start_codon`, `refseq_edited_model_equal_span_insertion_term` and `hgvsp_frameshift_count_at_the_first_residue`, which take the reference bases the allele of a row under an identifier stands in for, since such a row's first column spells no alleles; `--unordered` compares the shared entries without their order, for a flag set whose rows Ensembl VEP emits in hash order (`"unordered": true` in the corpus manifest). `--self-parity` holds one engine's four formats to each other column by column; `--warnings` compares two engines' `_warnings.txt` and skipped-variants files; a JSON report with a `verdict`, exit 1 on red. A Rust tool with its own Cargo root, which `feature_oracle.sh` starts; see its [README](../../tools/vep-formats-compare/README.md) |
| `../../tools/vep-stats-compare/` | Section-by-section comparator of two statistics reports, a Rust tool (the `<output>_summary.txt` that `--stats_text` writes; a corpus keeps Ensembl VEP's as `expected/<format>.summary.log`): `[VEP run statistics]`, `[Data version]` and `[General statistics]` row by row, the command line, the two times, the run time and the input and output rows compared by key alone and the paths of `Annotation sources` masked; every chart section as a set of `(key, count)` pairs, since Ensembl VEP sorts ties in hash order; `--tolerance N` lets a chart count differ by up to N; `--mask-row` sets one more row aside; the `Novel / existing variants` row of a run under a frequency check is held to its documented difference and named under `documented`, and so is a key of `[Consequences (all)]` that names no consequence type in a run under a `--custom` source keyed `Consequence`; the same sections in the same order on both sides; a JSON report with a `verdict`, exit 1 on red. `cargo test` there holds it to the Ensembl VEP report of the `GRCh38-statistics` corpus and to the goldens over its fixtures; the golden harness applies the same rules, those of the two documented differences excepted, to a manifest run's `stats_file`. See its [README](../../tools/vep-stats-compare/README.md) |
| `../../tools/vep-fastvep-normalize/` | Rewrites fastVEP's raw tab output so both engines' rows key alike before comparison: a Rust tool that streams the output across the host's cores. `run_clone_measurement.sh` runs the tree's release build on every fastVEP SNP/indel cell that scores, builds nothing, and stops before its first engine run where the build is missing (`cargo build --release --locked --manifest-path tools/vep-fastvep-normalize/Cargo.toml`) |
| `../../tools/vep-fastvep-plugin-tuples/` | Writes the supplementary columns of fastVEP's tab output (one per source loaded with `--sa-dir`) as plugin tuples (Location, Allele, Feature, Ensembl VEP field name, decoded value), the form in which a comparator's values are scored against Ensembl VEP's plugin and `--custom` fields, a Rust tool. Its data file `data/fastvep_field_map.tsv` names the Ensembl VEP field of each part it maps and the expected set a fixture source is scored against; its tests run the binary, over lines of fastVEP's own output among them |
| `../../tools/vep-file-census/` | The passes these launchers make over whole files, across the host's cores and in the bytes of the pipelines they stand for: the record and row counts of `generate_reference_output.sh` and of the Ensembl VEP cells of `run_clone_measurement.sh` (`grep -c`, and the annotated units of a VCF or JSON output), the contig-prefix pass of those cells, and the digests of inputs, outputs, the FASTA and every file of both caches (`sha256sum`) in `generate_reference_output.sh` and `run_concordance.sh`. A Rust tool: the first two launchers run the tree's release build, build nothing and stop before their first engine run where it is missing (`cargo build --release --locked --manifest-path tools/vep-file-census/Cargo.toml`); `run_concordance.sh` builds it through cargo |
| `count_warnings.py`          | Aggregates per-dataset ERROR/WARN/INFO counts from both engines' stderr logs into `reports/warning_counts.json` and a Markdown section (`--out-md`)                                                                                                                                    |
| `../../tools/vep-cache-parity-audit/` | Audits the Storable-to-JSON cache conversion for field coverage over a transcript sample, a Rust tool |
| `analyze_all_discordants.py` | Classifies the `discordant.tsv` of several datasets and builds a cross-dataset pattern map. Reads `<out>/suites/<id>/report(s)/discordant.tsv`, the file a `run_concordance.sh` run with `--work-dir <out>/suites/<id>` writes; see "Dataset severity" below.                              |
| `../../tools/vep-intermediates-compare/` | Joins a `trace_perl_intermediates.pl` TSV with vep-rs output on (location, allele, transcript) and reports which fields diverge; a Rust tool, run as `vep-intermediates-compare --perl-trace <tsv> --rust-output <output> [--discordant-only]` from its release build (`cargo build --release --locked --manifest-path tools/vep-intermediates-compare/Cargo.toml`) |
| `trace_perl_intermediates.pl` | Runs inside the Ensembl VEP container and dumps per-transcript intermediate values (CDS and translation coordinates, codons, peptides, consequence set) for a VCF of discordant variants                                                                                                                               |
| `trace_perl_variant.sh`      | Run Ensembl VEP on a single variant with injected `warn()` tracing in VariationEffect.pm, TranscriptVariationAllele.pm, and BaseTranscriptVariation.pm. Produces per-transcript predicate trace.                                                                                       |
| `generate_reference_output.sh` | Writes the Ensembl VEP ground truth `run_clone_measurement.sh` scores against (one `output.txt` per SNP/indel dataset and one `<vcf>.txt` per SV input, via Docker, with the harness's own invocation), a flag family's reference in any output format (`--flags`, `--format`, one container run per format), and a golden corpus's four expected outputs with its `provenance.json` (`--corpus`); the image digest is asserted before any run and `--stats` keeps Ensembl VEP's statistics report. See below. |
| `finalize_sv_gt_set.py` | Finalizes a per-VCF SV ground-truth set directory after a file is regenerated: refreshes that file's canonical-inputs manifest entry from the file, and writes the set's `provenance.json` (per-file digests of every reference output and canonical input, both assemblies), `PROVENANCE` and `STATUS` from the per-file sidecars `generate_reference_output.sh` writes. |
| `test_generate_reference_output.sh` | Tests for `generate_reference_output.sh`: a stub `docker` exercises the dry run and every fail-closed path                                                                                                                                                                       |
| `feature_oracle.sh`          | The comparison oracle: the harness that runs Ensembl VEP and vep-rs on the same inputs and compares every column of every row in every output format. Runs the comparison cases of a case table (input, cache, species, FASTA, flags, formats, Ensembl VEP image and mode, masks as classes or `<class>=<expected applied count>`, floors, field pairs, unordered, reference, refusal and tool per case; a `-` species is homo_sapiens, and the cache slice is read under `<species>/<version>_<assembly>`): Ensembl VEP in Docker once per case and format, offline or against Ensembl's database, then each candidate `vep` with the same flags (less `--genomes`, which selects Ensembl VEP's server) and `--species`, every format compared by `tools/vep-formats-compare` (the tree's release build, or the build `VEP_FORMATS_COMPARATOR` names as a path or a command on `PATH`; the runner builds nothing and stops before the first engine run where there is none), and a case naming `stats` has its two default runs' statistics reports compared by `tools/vep-stats-compare` (`VEP_STATS_COMPARATOR`, found the same ways); writes `<out>/<candidate_pin>/<case>/report.json` with a `verdict` and `<out>/status.tsv` (`case`, `candidate_pin`, `verdict`, `report_digest`, `run_utc`); a case's `reference` directory of stored outputs in the golden-corpus layout replaces its Ensembl VEP runs, as `--reference-dir` does for every case; a `refusal` directory of the two engines' expected error texts makes a case green when both refuse; a `fasta` cell marked `unindexed:` gives every engine run a copy of the file with no index beside it and holds the indexes the two engines write to each other byte for byte, and one marked `unindexed-readonly:` gives the copies in directories no engine can write, for a `refusal` case; `unordered` compares the rows as a multiset; a `tool` cell (`filter_vep:<args>`) compares Ensembl VEP's `filter_vep` with the candidate's `vep filter` over an annotated input in place of the engines; exit 1 on any red. A file `oracle_cases_on_demand.tsv` beside the case table (`case`, `reason`) lists population cases that run only when `--only` names them. `--discover` asks the comparator, in place of each masks cell, for the built-in classes (five that change what is compared only where named), the classes the case's manifest documents and the cell's own (a case its named-runs row keeps pooled, on a manifest of several runs over its input, is asked for none the manifest documents but the cell's own), and reports which applied, in a `discovery.tsv` per case and for the run, with a status of `discovered` or `incomplete` and no verdict (`--discover-classes` narrows the list). Every engine run and comparison is timed, and `report.json` and `run.json` carry a `sizing` block: rows compared, peak memory, the work disk's high-water mark and wall minutes. The header of the script documents the table's columns. |
| `test_feature_oracle.sh`     | Tests for `feature_oracle.sh`: a stub `docker` and a stub `vep` answer with the committed `tests/golden/116/GRCh38-release116` outputs, so the real format comparator decides each format's verdict, and the statistics checks compare through a stub; covers the table contract, the command shapes, every fail-closed path, the stored-reference, refusal, layout and species cases, the on-demand list, the discovery pass, the sizing block, and the verdicts                                             |
| `test_species_harness.py`    | The species as an input of the harness: `feature_oracle.sh` plans both engines under the case's species and slice, `generate_reference_output.sh` under a corpus manifest's `species` and `release`, both defaulting to homo_sapiens; every other species with a committed corpus under `tests/golden/116/` has one that covers the parity matrix row of its genome, found by its manifest, whose manifest, cache and reference record name its species, assembly and cache version, and whose reference generation is planned under them |
| `harness_genomes.sh`         | The harness's genomes as data: a file to source whose one function, `harness_genome_field <table> <species> <assembly> <column>`, prints one field of a genome's row of `harness_genomes.tsv` (species, assembly, the number Ensembl VEP's `--cache_version` takes, and the FASTA and converted cache as paths under the data root). The last row for a genome wins; `test_harness_genomes.sh` tests it, and [Genome checks](../../docs/genome-checks.md) states the table's columns. `generate_reference_output.sh` and `run_clone_measurement.sh` read each genome's FASTA, converted cache and cache version from the table, or from the one `VEP_HARNESS_GENOMES` names. |
| `extract_perl_reference.sh`  | Extracts five Ensembl VEP modules from the Docker image into a local directory for debugging and tracing.                                                                                                                                                                              |
| `../../tools/vep-prepare-benchmark-vcfs/` | Prepares the input VCFs of a `run_concordance.sh` run, a Rust tool the run builds through cargo: decompresses bgzipped files, normalizes chromosome names, keeps the canonical contigs, truncates for smoke mode, and writes `inputs/manifest.json` with each file's counts and digest |

## Comparing regulatory feature rows

Both field comparators accept `--keep-only regulatory_feature_rows` to compare only
`RegulatoryFeature` and `MotifFeature` consequence rows from **both** inputs:

```bash
tools/vep-fields-compare/target/release/vep-fields-compare \
  --perl perl.txt --rust rust.txt --out regulatory-fields.json \
  --keep-only regulatory_feature_rows

tools/vep-formats-compare/target/release/vep-formats-compare \
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

## Generating the Ensembl VEP ground truth

`run_clone_measurement.sh` only reads the ground-truth trees under `<data>/ground_truth/perl/`.
`generate_reference_output.sh` writes them with the same Docker image (`release_116.3`), cache
version (116) and `vep` flags the harness's Ensembl VEP engine uses, so the reference and the timed
runs cannot diverge; `--perl-image` and `--cache-version` move both together:

```bash
# Preflight everything and print the docker command for every output, running nothing.
scripts/concordance/generate_reference_output.sh --data-dir <data> --dry-run

# Generate both assemblies' SNP/indel and SV references (hours: ClinVar and gnomAD are large).
scripts/concordance/generate_reference_output.sh --data-dir <data>

# One assembly's SV set only, into an explicit location.
scripts/concordance/generate_reference_output.sh --data-dir <data> --assembly GRCh37 --datasets sv \
  --sv-gt-dir <data>/ground_truth/perl/sv_per_vcf
```

What it needs, all local (it fetches nothing and uses no cloud credentials):

- Docker and the `ensemblorg/ensembl-vep:release_116.3` image (`--perl-image` to override).
- The **full-genome** Ensembl VEP cache for each assembly, of the release the image pins (cache
  version 116, the `cache_version` of the genome's row of `harness_genomes.tsv`;
  `--cache-version` names one for every genome), extracted under one directory so it holds
  `homo_sapiens/116_GRCh37/` and `homo_sapiens/116_GRCh38/` (`--perl-cache-dir`,
  default `<data>/caches/perl/vep-cache`; another species' cache extracts beside them as
  `<species>/<version>_<assembly>/`, the layout Ensembl VEP resolves from `--species`,
  `--cache_version` and `--assembly`). Ensembl FTP:
  `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz`
  and `https://ftp.ensembl.org/pub/grch37/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh37.tar.gz`
  (about 24 GB each). A chromosome-subset cache is refused before anything runs: it does
  not fail a VEP run, it under-annotates and exits 0. A slice under homo_sapiens holds all
  25 primary contigs, and a slice of another species at least one chromosome directory.
- An indexed reference FASTA per genome: the `fasta` path of the genome's row of
  `harness_genomes.tsv`, under `<data>` (`<data>/reference/<asm>/genome.fa` for
  `homo_sapiens:GRCh37` and `homo_sapiens:GRCh38`), which
  `--fasta-for <species>:<assembly>=<path>`, repeatable, or `--fasta` with a single
  `--assembly` replaces. VEP receives `--fasta` on every call.
- The canonical inputs: `vcf/<asm>/all_variants_canonical/<base>.canonical.vcf.gz` for the
  six SNP/indel datasets and `<sv-gt-dir>/canonical_inputs/<asm>/*.vcf[.gz]` for the SV set.

Fork and buffer are part of the ground truth's identity: SNP/indel runs `--fork 16
--buffer_size 5000` and the SV set runs `--fork 4 --buffer_size 5000`, the settings of the
released references. VEP's transcript selection above `--max_sv_size` is
batch-dependent, so a different SV fork count changes the row count. `--fork` and `--sv-fork`
override them, but a reference generated at other settings is a different reference.

Every run is gated on annotation volume (rows >= input records, per dataset and per SV set,
the same floor the harness's `assert_perl_annotated` applies), refuses to overwrite an
existing reference without `--force`, and writes a `<name>.provenance.json` beside each
output recording the image digest, flags, input and output digests and counts. Before
anything runs, the pulled image's digest must match `--perl-image-digest` (default the
`release_116.3` tag's manifest-list digest, or `VEP_PERL_IMAGE_DIGEST`), so a moved tag
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
| `{plugins}`, `{generated}` | `{plugins}` is the plugin directory each engine reads: for Ensembl VEP the `--dir-plugins` directory, mounted read-only at `/work/plugins`, and for vep-rs the corpus's `plugin_data`. `{generated}` is the directory the manifest's `generated` block writes, which the golden tests fill before the run | Beside a `manifest.json`, the corpus's `plugin_data` and `generated` directories, read as `{corpus}/plugin_data` and `{corpus}/generated`, so a cell written as a run of the manifest writes it is that run; with no manifest beside the input the cell is refused |

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

- `summary.json`: machine-readable metrics (precision, recall, F1, tuple counts)
- `summary.md`: human-readable report
- `discordant.tsv`: every tuple that the output of Ensembl VEP or of vep-rs holds and the other lacks, one row each
- `discordant_open.tsv`: the rows of `discordant.tsv` that no excluding rule removes, in the same order and layout; its row count is the adjusted one-sided total
- `discordant.tsv.partial`: where the tuple comparison stops while it compares the output files, the rows written before the stop stand under this name in place of `discordant.tsv`, and the error line names it; where it stops on a key file that does not ascend, the directory holding the two key files is kept and named in that line as well
- `provenance.json`: exact versions, flags, plugin argument strings, and cache checksums used

## How It Works

1. **Input preparation**: `tools/vep-prepare-benchmark-vcfs`, which the run builds through cargo before any engine runs, normalizes and optionally truncates the input VCF
2. **Ensembl VEP**: Runs via Docker with matching flags (offline, cache, same plugins)
3. **vep-rs**: Runs the locally-built `target/release/vep` binary with equivalent flags
4. **Comparison**: the tuple comparator extracts semantic tuples from both outputs and computes set-based metrics. It is `tools/vep-outputs-compare`, which the run builds through cargo before any engine runs
5. **Reporting**: Results written to `<work-dir>/reports/`

## Dataset severity

`analyze_all_discordants.py` reads an optional `<out>/dashboard.json` of the form
`{"suites": [{"id": "<dataset-id>", "name": "<display name>", "severity": "<severity>"}]}`.
The severity decides how many representative discordant rows the pattern map keeps
per category and dataset:

| Severity   | Discordant samples retained |
| ---------- | --------------------------- |
| `blocking` | 10                          |
| `major`    | 5                           |
| `minor`    | 5                           |

A dataset with no severity set, or no `dashboard.json` at all, is treated as `minor`.
