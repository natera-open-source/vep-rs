# vep-cli

The command line of vep-rs: the `vep` binary, a drop-in for the Ensembl Variant Effect
Predictor (VEP) command line, and the `vep_cli` library behind it. The crate's two binaries
(`vep` and `bench_plugins`) call the library, `vep-cache-converter` lists it as a dependency,
and a program that embeds annotation holds its `Annotator` across many batches. It reads
variants in Variant Call Format (VCF) or VEP's other input formats, annotates them against a
JSON transcript cache, and writes VEP's default, tab, VCF, JSON or Parquet output through the
`vep-io` formatters.

## How it works

Follow one VCF record. `main.rs` parses the command line into `Args`; `Config::from_parsed`
resolves the defaults, the fields `--everything` and `--vcf` imply, and the automatic thread
counts; `Runner::run` refuses to overwrite an existing output without `--force_overwrite`,
builds a rayon thread pool of `--fork` workers, opens the indexed FASTA (the reference
genome's sequence file) when `--fasta` names one, loads each `--plugin` argument as a built-in
from `vep-builtins` or, when no built-in claims the name, as a C-ABI (application binary
interface) dynamic library through `vep-plugin`, enumerates the cache's `transcripts/<chr>/`
shards into a `LazyTranscriptIndexes` that parses a chromosome the first time a variant asks
for it, and opens the cache's variation files (`variations/<chr>/all_vars.gz`, one per
chromosome and read through its index, or `<start>-<end>_var.gz` region files parsed whole)
only when `--check_existing`, or a flag that implies it such as `--af` or `--pubmed`, needs
them. The run then proceeds as three stages joined by bounded channels. The reader
thread scans the input line by line: a `#` line is held for the VCF writer's header block, a
data line is appended to the current `LineBatch` with its line number, and the first record of
a new chromosome asks a helper thread to warm that chromosome's transcript index and the next
one in the input's tabix contig order. The coordinator, the calling thread, parses the batch's
lines on the pool: `parse_vcf_line_with_max_sv_size` routes the record through `vep-core`'s
structural-variant rules and splits a sequence record into one `InputVariant` per alternate
allele, alleles trimmed. Each variant then passes through `annotate_one`: `annotate_variant`
visits every transcript overlapping the variant's span and calls
`vep_effects::calculate_consequences` on each, a missense consequence gains its SIFT (Sorting
Intolerant From Tolerant) and PolyPhen (Polymorphism Phenotyping) predictions, a breakend
gains its mate's consequences, `variation_matcher` attaches the co-located known variants, and
a variant on a cached chromosome that overlaps nothing becomes `intergenic_variant`. Built-in
plugins see the whole batch in one `run_batch` call and dynamic-library plugins run per
consequence; `pick::apply_to_batch` applies the pick flags and summaries to each record's
alleles together; `render_record` hands each record's alleles to the chosen `vep-io` formatter,
or to `render_default_record` for VEP's default text format. Without plugins the coordinator
fuses annotation, filtering, rendering and freeing into one pass over the pool per batch. The
writer thread receives the rendered batches in sequence order and streams their bytes to the
output; a Parquet run writes a TSV (tab-separated values) intermediate and converts it with a
DuckDB subprocess once the pipeline has drained.

## Modules

| Module                      | Role                                                                                                                                                                                                                        |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/main.rs`               | The `vep` binary: parses `Args`, sets the log level from `--verbose` and `--quiet`, prints the version banner, builds `Config` and runs `Runner`; mimalloc is the global allocator under the default `mimalloc` feature     |
| `src/lib.rs`                | The library root: re-exports the embedding surface (`Annotator`, `Config`, `apply_to_record`, the transcript index types, `InputVariant`, `JsonOutputFormatter`) and declares the modules                                      |
| `src/args.rs`               | The clap definition of every flag, named as Perl VEP names them (underscores included); a flag the error table of `config.rs` governs takes its help text from that table                                                                     |
| `src/config.rs`             | `Config::from_parsed`: the error table of refused and no-effect flags, the fields the option sets imply, the `--fork` and `--decompression_threads` automatic caps, `--distance`, tilde expansion and the cache path the output header prints                       |
| `src/option_rules.rs`       | Ensembl VEP's option rules as data (`Config.pm:327-645`): the lists, option sets, valid values, requirements, incompatibilities and deprecations, and the flag state they are applied to and checked on |
| `src/cache_resolve.rs`      | `resolve`: the cache directory from `--json_cache`, else from `--dir_cache`, `--species`, `--refseq`/`--merged`, `--cache_version` and `--assembly` in Ensembl VEP's layout; requires `transcripts/` and `info.json` and an assembly that agrees with the cache, and picks up an indexed FASTA inside it unless `--fasta` or `--no_fasta` is given |
| `src/runner.rs`             | `Runner::run`: opens input and output, loads FASTA, plugins, transcripts and variations, picks the output formatter, drives the pipeline and finalizes Parquet; owns `annotate_one`, `annotate_batch` and `annotate_variant` |
| `src/pipeline.rs`           | The three-stage pipeline: reader thread, coordinator and writer thread over bounded channels, the resolution of each input chromosome name to the cache's through the synonym table, the chromosome prewarmer and the fused annotate-render pass |
| `src/warnings.rs`           | The warning `Sink`: every `warn!` event and skipped-variant report, written once per text to the warnings file (`--warning_file`, default `<output_file>_warnings.txt`) and echoed on standard error unless `--quiet`; the `--skipped_variants_file` written when the run ends |
| `src/render.rs`             | `render_record` and `render_default_record`: one input record's rows, through the chosen `vep-io` formatter or as VEP's default text format, appended to a render task's output bytes                                        |
| `src/output_compress.rs`    | `OutputCompression`: `--compress_output` gzip or bgzip encoding of the output in process, in place of the `gzip -c` or `bgzip -c` pipe Ensembl VEP opens; `CompressedOutput::finish` writes the stream's trailer once every line is written |
| `src/annotator.rs`          | The embedding entry point: `Annotator` over `AnnotationResources` and a scoped pool, the `RawRow` simple-allele input, the per-row `AnnotationError`, and the shared loaders for FASTA, plugins and transcripts            |
| `src/json_cache.rs`         | Deserializes the JSON cache: enumerates `transcripts/<chr>/<start>-<end>.json` shards and loads one chromosome's transcripts                                                                                                 |
| `src/variation_source.rs`   | `VariationSource`: the cache's variation files under `variations/<chr>/`, opened on first use and read through their index for the records a variant's window reaches, with a per-batch prefetch                           |
| `src/transcript_index.rs`   | Transcript interval indexes (UCSC binning, coitrees, sorted array) behind one `TranscriptIndex` surface, each chromosome's protein-coding transcription start sites beside it; `LazyTranscriptIndexes` loads each chromosome once; `MaterializedTranscriptIndexes` proves every chromosome loaded |
| `src/vcf_parser.rs`         | Parses one VCF data line into `InputVariant`s: the fixed fields, allele trimming and minimisation per alternate allele, structural records through `vep-core`'s routing, the `--max_sv_size` limit                          |
| `src/variation_matcher.rs`  | Matches input variants to the known variants of their window as Ensembl VEP does (by allele after trimming from either end, strand-aware, span-only for null alleles and under `--no_check_alleles`, `--failed` applied), and writes the Extra column fields (IMPACT, STRAND, the frequency and clinical fields, plugin keys) |
| `src/pick.rs`               | `FilterConfig`, `apply_to_record` and `apply_to_batch`: the record-level steps over an input record, the row filters of `src/row_filter.rs` first, then VEP's pick algorithm and summaries for `--pick`, `--per_gene`, `--pick_order`, `--most_severe`, `--summary` and the `--flag_pick` variants |
| `src/row_filter.rs`         | `RowFilters` and `affects_cds`: the `--no_intergenic` and `--coding_only` row filters, applied to each input record before the pick family as Ensembl VEP applies them |
| `src/filter_set.rs`         | `FilterSet`: Ensembl VEP's filter grammar (`field operator value` units joined by `and` and `or`, grouped by parentheses and inverted by `not`), compiled from `--transcript_filter` over the cache's transcript keys and evaluated on each overlapping transcript before its consequences are calculated |
| `src/dylib_plugins.rs`      | Loads the `--plugin` names no built-in claimed as C-ABI dynamic libraries and runs them over a batch, merging their JSON results into each consequence's `plugin_data`                                                        |
| `src/bin/bench_plugins.rs`  | The `bench_plugins` binary: times the phases of an annotation run with and without plugins and writes a JSON report                                                                                                          |

## Command line

`vep` has no subcommands: one invocation annotates one input. The flags are Perl VEP's,
underscores included, so an existing VEP command line parses; a flag the help marks "accepted
for VEP compatibility" parses and changes no annotation (the cache flags among them only name
the cache path the output header prints). The groups below name the flags that shape a run;
the full list, with every flag's parity note, is
[docs/cli-reference.md](../../docs/cli-reference.md).

| Group     | Flags                                                                                                                                                                                                                                                      |
| --------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Input     | `-i`/`--input_file` (default `STDIN`); `--format` vcf, ensembl, hgvs, region, or `guess`, which reads the file extension                                                                                                                                   |
| Output    | `-o`/`--output_file` (default `variant_effect_output.txt`, `STDOUT` for standard output); `--output_format` vep, vcf, json, tab or parquet, with `--vcf`, `--json` and `--tab` as shorthands; `--force_overwrite`; `--no_headers`; `--vcf_info_field` (default `CSQ`) |
| Cache     | `--dir_cache` (else `--dir`), `--species`, `--cache_version`, `--refseq` and `--merged` select the cache in Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout; `--json_cache <dir>` (alias `--full_cache_dir`) names the directory itself; `--show_cache_info` prints its source versions; `--no_fasta` ignores a FASTA inside it; `--synonyms <file>` replaces the cache's `chr_synonyms.txt`, through which an input chromosome name reaches the cache's |
| Reference | `--fasta`/`--fa`, an indexed FASTA; `--hgvs`, `--hgvsg` and `--shift_3prime 1`, which require it; `--assembly`, read by the assembly-sensitive plugins                                                                                                                              |
| Runtime   | `--fork`/`--threads` (0 = every logical CPU up to 32); `--decompression_threads` (0 = up to 10, for `.gz` input); `--buffer_size` (default 5000 variants per batch); `--max_sv_size` (-1 lifts the limit)                                                   |
| Fields    | `--everything` and the field flags it implies: `--symbol`, `--canonical`, `--mane`, `--hgvs`, `--sift`, `--polyphen`, `--numbers`, `--domains`, `--variant_class`, `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`, `--check_existing` and the rest |
| Picking   | `--pick`, `--pick_allele`, `--per_gene`, `--pick_allele_gene`, `--most_severe`, `--summary`, `--flag_pick`, `--flag_pick_allele`, `--flag_pick_allele_gene`                                                                                                 |
| Plugins   | `--plugin Name,key=value,...`, repeatable; `--dir_plugins`, where dynamic-library plugins are found                                                                                                                                                         |
| Parsing   | `--dont_skip` aborts on a malformed line instead of skipping it; `--allow_non_variant` passes reference-only VCF records through unchanged                                                                                                                  |
| Verbosity | `-v`/`--verbose` (debug logging); `-q`/`--quiet` (errors only, no progress line)                                                                                                                                                                           |

Three declared flags are refused with an error rather than ignored, because the feature they
belong to is not implemented: `--regulatory_gff`, `--extended_promoters` and
`--custom_suppress_filter`. `--regulatory` is accepted and warns that no regulatory consequence
is emitted.

```bash
vep -i input.vcf -o output.txt --json_cache /path/to/json_cache --offline --everything
```

`bench_plugins`, the crate's second binary, takes `--input`, `--json-cache`, repeated `--plugin`
arguments in the `vep --plugin` format, `--variants`, `--runs` and `--output-json`, and reports
per-phase timings as JSON.

## Public API

The re-exports in `lib.rs` are the embedding surface; the submodules stay `pub` for the in-tree
binaries and tests, and their contents may move between releases.

| Item                                                                                        | Use                                                                                                                                                                                                                 |
| ------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Config`, `Config::from_argv`                                                               | The resolved run configuration; `config::resolve_fork` applies the same automatic thread count to an embedder's own setting                                                                                        |
| `Annotator`                                                                                 | One loaded cache plus a scoped thread pool, reused across batches: `from_config` loads every chromosome up front, `annotate_batch` annotates `InputVariant`s in place, `try_annotate_rows` takes `RawRow`s         |
| `AnnotationResources`, `build_thread_pool`                                                  | The shared transcripts, variation data, effects configuration and plugins behind an `Annotator`, and the pool it runs on                                                                                           |
| `AnnotationStats`, `AnnotationError`, `RawRow`, `validate_raw_row`, `row_to_input_variant`  | The per-batch counts, the row-level error taxonomy, and the simple-allele input (chromosome, 1-based position, reference allele, alternate allele) with its validation and conversion                               |
| `apply_to_record`, `apply_to_batch`, `FilterConfig`, `RowFilters`, `PickCriterion` | The record-level steps over the annotated alleles of one input record: the `--coding_only` and `--no_intergenic` row filters, then VEP's pick algorithm and summaries |
| `LazyTranscriptIndexes`, `MaterializedTranscriptIndexes`, `TranscriptIndexImpl`, `TranscriptSet` | The per-chromosome transcript index, the guard that an `Annotator` requires before fanning out to the pool, the index implementation choice, and the transcript set (every transcript, or a GENCODE set) the indexes keep |
| `InputVariant`, `Consequence`, `Impact`, `TranscriptConsequence`                            | Re-exported from `vep-core`, so a caller constructs and inspects annotations without a direct dependency                                                                                                             |
| `JsonOutputFormatter`                                                                       | Re-exported from `vep-io` for an embedder that renders JSON                                                                                                                                                         |
| `runner::Runner`                                                                            | The file-to-file pipeline the `vep` binary drives: `Runner::new(config)` resolves the cache, then `.run()`                                                                                                           |

## Run it

```bash
cargo build -p vep-cli                     # target/debug/vep and target/debug/bench_plugins
cargo test -p vep-cli                      # unit tests plus the integration suites under tests/
cargo run -p vep-cli --bin vep -- --help   # the full flag list
```

`cargo build -p vep-cli --no-default-features --bin vep` builds `vep` on the system allocator
instead of mimalloc. Each integration suite under `tests/` covers one contract:
`integration_test.rs` runs the built binary as a subprocess; `golden.rs` and `format_parity.rs`
compare every output format on the golden corpora with the Ensembl VEP files committed beside
them; `pipeline_parity.rs` checks that the output bytes do not change with the batch size or the
thread count; `parity_matrix.rs` checks that every long flag has a row in
`tests/parity/matrix.tsv` at the repository root; `dylib_plugins.rs` exercises the
dynamic-library plugin path. `cargo bench -p vep-cli --bench cache_load_bench` times cache
loading and skips itself when `$VEP_JSON_CACHE_DIR` names no cache. The repository's standard
commands (the workspace test run, clippy, fmt, the release build) are in
[CONTRIBUTING.md](../../CONTRIBUTING.md).

## Related

- Depends on [vep-core](../vep-core/README.md) (variant, transcript and consequence types),
  [vep-effects](../vep-effects/README.md) (the consequence engine),
  [vep-io](../vep-io/README.md) (input parsers and output formatters),
  [vep-builtins](../vep-builtins/README.md) (the built-in plugins),
  [vep-plugin](../vep-plugin/README.md) (the C-ABI plugin loader) and
  [vep-fasta](../vep-fasta/README.md) (the indexed FASTA reader).
- Depended on by [vep-cache-converter](../vep-cache-converter/README.md).
- Behaviour: [docs/cli-reference.md](../../docs/cli-reference.md) (every flag),
  [docs/output-formats.md](../../docs/output-formats.md) (each format against Ensembl VEP's),
  [docs/plugins.md](../../docs/plugins.md) (plugin arguments and data files),
  [docs/cache-setup.md](../../docs/cache-setup.md) (building the `--json_cache` directory) and
  [docs/intended-divergences.md](../../docs/intended-divergences.md).
- The root [README](../../README.md) and [CONTRIBUTING.md](../../CONTRIBUTING.md).
