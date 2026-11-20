# vep-cli

The command line of vep-rs: the `vep` binary, a drop-in for the Ensembl Variant Effect
Predictor (VEP) command line, and the `vep_cli` library behind it. The crate's two binaries
(`vep`, and `bench_plugins` under the `bench` feature) call the library, `vep-cache-converter`
lists it as a dependency, and a program that embeds annotation holds its `Annotator` across
many batches. It reads variants in Variant Call Format (VCF) or VEP's other input formats,
annotates them against a JSON transcript cache, and writes VEP's default, tab, VCF, JSON or
Parquet output through the `vep-io` formatters.

```bash
cargo install vep-cli --locked
```

Documentation: https://natera-open-source.github.io/vep-rs/

## How it works

Follow one VCF record. `main.rs` parses the command line into `Args`; `Config::from_parsed`
resolves the defaults and the fields `--everything` and `--vcf` imply; `Runner::run` opens the
indexed FASTA (the reference genome's sequence file) when `--fasta` names one, `--fasta_dir`
holds one for the species and assembly or, without `--no_fasta`, the cache directory holds
one, enumerates the cache's `transcripts/<chr>/` shards into a `LazyTranscriptIndexes` that
parses and indexes a chromosome when it is first warmed (below), opens the input, resolving
the automatic `--decompression_threads` count against the input's shape and, under
`--format id`, listing the cache's variation files the identifiers are resolved through, and
refuses to overwrite an existing statistics report or output without `--force_overwrite`. It
then resolves the automatic `--fork` count against the same shape, builds a rayon thread pool
of `--fork` workers, loads each `--plugin` argument as a built-in
from `vep-builtins` or, when no built-in claims the name, as a C-ABI (application binary
interface) dynamic library through `vep-plugin`, and lists the cache's variation files
(`variations/<chr>/all_vars.gz`, one per chromosome and read through its index, or
`<start>-<end>_var.gz` region files parsed whole, each opened on first use) for the co-located
variants only when `--check_existing`, or a flag that implies it such as `--af` or `--pubmed`,
needs them. Every annotation source,
the regulatory features of `--regulatory` or `--regulatory_gff` among them, opens before the
statistics report and the output file are created, so a refused source leaves neither file
behind. The run then proceeds as a pipeline of stages that holds a bounded number of batches
in flight: a reader thread, the coordinator and a writer thread, with a prefetch stage between
the reader and the coordinator in the batched mode (below). The reader
thread scans the input line by line: a `#` line is held for the VCF writer's header block, a
data line is appended to the current `LineBatch` with its line number, and the first record of
a new chromosome asks a helper thread to warm that chromosome's transcript index and the next
one in the input's tabix contig order. The batch's lines are parsed on the pool, by the
streaming tasks of the coordinator (the calling thread) or as one pass of the batched mode's
prefetch stage (both below): `parse_data_line` reads the record once (its columns, its `INFO` view and its
tandem-repeat warnings), routes it through `vep-core`'s structural-variant rules (a
`<CNV:TR>` record expanded to explicit alleles when a FASTA is given) and splits a sequence
record into one `InputVariant` per alternate allele, alleles trimmed, a bi-allelic record's
from both ends under `--minimal`; the record then passes the checks of `validate`, its
deletion test reading the record's parsed `INFO` (the view, or a copy of the column for a
line that is not UTF-8), after which `split_minimal_record` trims each allele of a
multi-allelic record under `--minimal`. Under `--check_frequency` the prefetch stage
co-locates the batch first and `frequency_filter` drops each allele and record that the
population's frequencies rule out, as Ensembl VEP's variation source filters its buffer
before any transcript is read. Each remaining variant then passes through `annotate_one`:
`annotate_variant` visits every transcript overlapping the variant's span and calls
`vep_effects::calculate_consequences` on each, a single-residue substitution gains its SIFT (Sorting
Intolerant From Tolerant) and PolyPhen (Polymorphism Phenotyping) predictions, a breakend
gains its mate's consequences, `variation_matcher` attaches the co-located known variants,
`sv_source` the known structural variants of `--check_svs`, and a variant that overlaps
nothing becomes `intergenic_variant`. Under `--regulatory`,
`RegulatorySource::annotate_variant` then adds the variant's `RegulatoryFeature` and
`MotifFeature` rows after its transcript rows, keeping the `intergenic_variant` row of an
allele that overlaps no transcript. Built-in
plugins see the whole batch in one `BuiltinRegistry::annotate_batch` call and dynamic-library
plugins run per consequence whose feature type they declare (transcript rows unless the plugin
exports `vep_plugin_feature_types`), a plugin that declares `Intergenic` also on the intergenic
row of a variant that overlaps no transcript, whether or not the variant has regulatory or
motif rows; `pick::apply_to_batch` applies the pick flags and summaries to each record's
alleles together; `render_record` hands each record's alleles to the chosen `vep-io` formatter
(the tab format's byte path, `render_tab_record`, included), or to
`render_default_record` for VEP's default text format. Without plugins, `--custom` sources, a
known-variant source or the `--check_frequency` filter the coordinator streams: it warms the
chromosomes a batch names and hands
the batch to the pool as independent tasks that each annotate, filter, render and free their own
records (and parse them first, for a VCF input; the other formats are parsed on the reader
thread), then takes the next batch without waiting for them; with plugins or `--custom` sources,
with the cache's variation files or the `--check_svs` file open, or under the `--check_frequency`
filter, the prefetch stage parses each batch, warms its chromosomes and reads its sources ahead
while the coordinator annotates the batch before it, and each batch then crosses the pool
whole, one pass after another
(`--pipeline_mode batched` asks for that path on any run). The writer thread puts the rendered
chunks back in input order and writes them to the output; a Parquet run renders each chunk
into typed Arrow batches on the pool, and the writer thread hands them, in input order, to the
Parquet writer, whose own thread encodes each full block on the pool into the parts of the
output directory. The warnings a task's records draw while they are annotated and filtered
(`Transcript-assembly mismatch`, `Unable to assign consequence type`) travel to the writer with
its chunk and are reported once the batch is complete; the batched coordinator reports them
once the batch is annotated. In both modes the warnings file lists a record's warnings in input
order, after the skipped-record warnings of its batch.

## Modules

| Module                      | Role                                                                                                                                                                                                                        |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/main.rs`               | The `vep` binary: parses `Args`, sets the log level from `--verbose` and `--quiet`, builds `Config`, answers `--resolved_config` before any cache opens, prints the version banner, and runs `Runner`; mimalloc is the global allocator under the default `mimalloc` feature |
| `src/lib.rs`                | The library root: re-exports the embedding surface (`Annotator`, `Config`, `apply_to_record`, the transcript index types, `InputVariant`, `JsonOutputFormatter`) and declares the modules                                      |
| `src/args.rs`               | The clap definition of every flag, named as Ensembl VEP names them (underscores included); a flag the error table of `config.rs` governs takes its help text from that table                                                                  |
| `src/argv.rs`               | The command line before clap reads it: `prepass` lays the flags of the `--config` file, `<dir>/vep.ini` and the `VEP_<FLAG>` environment variables out around the command line in Ensembl VEP's priority (`Config.pm:693-698`), and resolves Ensembl VEP's option spellings against the clap definition |
| `src/config.rs`             | `Config::from_parsed`: the error table of refused and no-effect flags, the fields the option sets imply, the `--fork` and `--decompression_threads` automatic caps, `--distance`, tilde expansion and the cache root the cache directory is resolved under                       |
| `src/option_rules.rs`       | Ensembl VEP's option rules as data (`Config.pm:327-645`): the lists, option sets, valid values, requirements, incompatibilities and deprecations, and the flag state they are applied to and checked on |
| `src/cache_resolve.rs`      | `resolve`: the cache directory from `--json_cache`, else from `--dir_cache`, `--species`, `--refseq`/`--merged`, `--cache_version` and `--assembly` in Ensembl VEP's layout; requires `transcripts/` and `info.json` and an assembly that agrees with the cache, and picks up an indexed FASTA inside it unless `--fasta` or `--no_fasta` is given |
| `src/runner.rs`             | `Runner::run`: opens input and output, loads FASTA, plugins, transcripts and variations, picks the output formatter (the Parquet writer with its footer metadata and chromosome order included), and drives the pipeline; owns `annotate_one`, `annotate_batch` and `annotate_variant` |
| `src/input_shape.rs`        | `InputShape`: whether the input is an indexed BGZF file and, from four index-resolved windows of blocks, how many bytes its records inflate to and how many it holds, the signals `--decompression_threads auto`, `--fork 0` and `--resolved_config` read |
| `src/pipeline.rs`           | The staged pipeline: reader thread, coordinator (streaming by default; batched, behind a prefetch stage, on request or when plugins, `--custom` sources, a known-variant source or the `--check_frequency` filter are in play) and writer thread, a bounded number of batches in flight between them, the resolution of each input chromosome name to the cache's through the synonym table (a bare `M` renamed `MT` over the chromosomes of every annotation source, as Ensembl VEP renames it), the chromosome prewarmer, the fused annotate-render pass and the output sink behind `--writer_strategy` |
| `src/affinity.rs`           | The CPUs a thread may run on, as the kernel's scheduler affinity holds them: read for the run's own set, and set when `config::reserve_reader_core` gives the reader thread a core of its own, the reader onto its CPU and every other thread of the process onto the rest; Linux carries the calls, and on any other system nothing moves |
| `src/warnings.rs`           | The warning `Sink`: every `warn!` event and skipped-variant report, written once per text to the warnings file (`--warning_file`, default `<output_file>_warnings.txt`) and echoed on standard error unless `--quiet`; the `--skipped_variants_file` written when the run ends |
| `src/stderr.rs`             | Standard error as the `vep` run writes it: `write`, `write_or_drop` and `Stream`, which the banner, the status, progress and summary lines, the echo of a warning and a refusal's `ERROR:` line pass through; once the process calls `drop_refused_writes`, as the `vep` binary does at its start, a write the stream refuses is dropped, so an unwritable stream never ends a run |
| `src/stats.rs`, `src/stats/` | The run statistics Ensembl VEP keeps (`Stats.pm`): each pool thread counts the records it renders and their rows into its own `ThreadStats`, folded into the run's `RunStats` once the run is complete; a `Report` built from those is rendered as the text report of `--stats_text` and as the HTML report written by default (`html.rs`) |
| `src/render.rs`             | `render_record` and `render_default_record`: one input record's rows, through the chosen `vep-io` formatter (the tab byte path included) or as VEP's default text format, appended to a render task's output bytes; a row whose plugin wrote a key named as one of the thirteen columns takes every column from the `vep-io` field table, so the plugin's value prints in the column |
| `src/output_compress.rs`    | `OutputCompression`: `--compress_output` gzip or bgzip encoding of the output in process, in place of the `gzip -c` or `bgzip -c` pipe Ensembl VEP opens; `CompressedOutput::finish` writes the stream's trailer once every line is written |
| `src/annotator.rs`          | The embedding entry point: `Annotator` over `AnnotationResources` and a scoped pool, the `RawRow` simple-allele input, the per-row `AnnotationError`, and the shared loaders for plugins and transcripts            |
| `src/json_cache.rs`         | Deserializes the JSON cache: enumerates `transcripts/<chr>/<start>-<end>.json` shards and loads one chromosome's transcripts                                                                                                 |
| `src/regulatory/`           | `RegulatorySource` (`source.rs`): the cache's regulatory and motif features under `regulatory/<chr>/` with their binding matrices from `regulatory/binding_matrices.json` (a baked shard element's `source` and `ordinal` name the file its record came from and the record's place in that file; a converted cache's features are one file in shard order), or under `--regulatory_gff` the tabix-indexed GFF files of `gff.rs` (one source per key, the accepted types, the promoter widening of `--extended_promoters` with a bound that is not a whole number ignored and warned about, a record without `ID` or `Name` named by the file's own chromosome spelling with its coordinates and type, the warning on two distinct records of one file sharing an identifier, the motif matrix stubs); one chromosome loaded on first use into the start-sorted index of `index.rs`, which ranks each feature in the order it was met (its file, then its place in the file); `rows.rs` adds a record's `RegulatoryFeature` and `MotifFeature` rows with their terms and motif fields, a structural variant's rows carrying that met rank as their order, one regulatory row for a sequence variant inside two records of one file sharing an identifier, the later record's |
| `src/variation_source.rs`   | `VariationSource`: the cache's variation files under `variations/<chr>/`, opened on first use and read through their index for the records a variant's window reaches, with a per-batch prefetch                           |
| `src/variation_name_index.rs` | `resolve_names`: the `--format id` pass over the cache's variation files, every name of the list held in memory and every file read once, each name mapped to its variation features                                   |
| `src/transcript_lookup.rs`  | `TranscriptLookup`: the stable-id, translation-id and gene-name index over the run's transcript cache behind HGVS `c.`, `n.` and `p.` input, built on the first lookup from the reference sequences' shards and extended to the alternate loci, patches and LRG regions only when a stable id is found on none of them |
| `src/transcript_index.rs`   | The per-chromosome transcript interval index (UCSC binning), each chromosome's protein-coding transcription start sites beside it; `LazyTranscriptIndexes` loads each chromosome once; `MaterializedTranscriptIndexes` proves every chromosome loaded |
| `src/nearest.rs`            | `--nearest`: the protein-coding transcription start sites nearest a record, found as Ensembl VEP's `TranscriptTree::nearest` finds them over the chromosome's site table, and the identifiers reported for them                                                                                   |
| `src/vcf_parser.rs`         | Parses one VCF data line into `InputVariant`s: the fixed fields, allele trimming and minimisation per alternate allele, structural records through `vep-core`'s routing, the `--max_sv_size` limit                          |
| `src/delimiter.rs`          | `Delimiter`: the `--delimiter` value as Ensembl VEP's column-based parsers read it, a regular expression every data line is split on; the patterns format detection chooses and plain text are split directly, and any other value is compiled and split as a pattern |
| `src/validate.rs`           | Ensembl VEP's record validation, `validate_vf` (`Parser.pm:457-649`), over the alleles parsed from one input record: the `--chr` list, the coordinates, the `--max_sv_size` bound, the chromosome against the annotation sources (the cache's and every `--custom` or `--regulatory_gff` file's), the allele string, an insertion's coordinates, the reference allele's length, then `--check_ref` or `--lookup_ref`; a failing record is reported through the warning sink and dropped, or kept under `--dont_skip` |
| `src/individuals.rs`        | The per-individual readings of a VCF record: under `--individual` each named sample's genotype, read from the record's sample columns, becomes a copy with its own allele string and zygosity; under `--individual_zyg` the record stays one feature whose `ZYG` list carries every named sample's zygosity; and the reference-allele row of `--process_ref_homs` |
| `src/frequency_filter.rs`   | `FrequencyFilter`: the `--check_frequency` stage (`--freq_pop`, `--freq_freq`, `--freq_gt_lt`, `--freq_filter`, `--filter_common`), Ensembl VEP's population check, frequency lookup with the strand, remainder and matched-allele rules, the include and exclude verdicts, the `FREQS` entries and the rewritten allele string of a record that lost alleles |
| `src/variation_matcher.rs`  | Matches input variants to the known variants of their window as Ensembl VEP does (by allele after trimming from either end, strand-aware, span-only for null alleles and under `--no_check_alleles`, `--failed` applied), and writes the Extra column fields (IMPACT, STRAND, the frequency and clinical fields, plugin keys) |
| `src/sv_source.rs`          | `SvSource`: Ensembl's structural-variant VCF of `--check_svs <file>`, read through its index for the top-level records overlapping each variant, with a per-batch prefetch; the names fill `SV` and the JSON `colocated_structural_variants`                                                                                  |
| `src/pick.rs`               | `FilterConfig`, `apply_to_record` and `apply_to_batch`: the record-level steps over an input record, the row filters of `src/row_filter.rs` first, then VEP's pick algorithm and summaries for `--pick`, `--per_gene`, `--pick_order`, `--most_severe`, `--summary` and the `--flag_pick` variants |
| `src/row_filter.rs`         | `RowFilters` and `affects_cds`: the `--no_intergenic` and `--coding_only` row filters, applied to each input record before the pick family as Ensembl VEP applies them |
| `src/filter_set.rs`         | `FilterSet`: Ensembl VEP's filter grammar (`field operator value` units joined by `and` and `or`, grouped by parentheses and inverted by `not`), compiled from `--transcript_filter` over the cache's transcript keys and evaluated on each overlapping transcript before its consequences are calculated |
| `src/filter_regex.rs`       | `translate`: Perl's regular expression dialect, as Ensembl VEP reads a filter's pattern (`/$pattern/i`, `FilterSet.pm:714-715`), rewritten for the `regex` crate that `filter_set.rs` compiles it with; the pattern runs over the value's bytes under Perl's default rules, and over their Latin-1 reading once a `\p{}` or `\P{}` class, a `\x{}`, `\o{}` or octal escape above 255 or a `\N{U+}` escape puts it under Perl's Unicode rules, where a case-insensitive match compares full case folds, or a `u`, `a` or `aa` flag switches a group; a construct the crate lacks (a back-reference, look-around) is refused naming the pattern |
| `src/filter_cli.rs`         | `vep filter`: Ensembl VEP's `filter_vep` over a VCF or tab output of `vep` or of Ensembl VEP, filtered line by line through the grammar of `filter_set.rs`, with `--only_matched`, `--soft_filter`, `--count`, `--list`, `--start`, `--limit` and `--test` as `filter_vep` reads them; an option is read after one dash or two, in any case, abbreviated to any prefix that names one option |
| `src/filter_input.rs`       | The lines `vep filter` reads as `filter_vep` reads them: the header block of an annotated output, which names the fields a filter may read, the VCF or tab format of its data lines, and each line parsed into the record a filter is evaluated on |
| `src/recode.rs`             | `vep recode`: a variant translated between its identifiers and notations as Ensembl VEP's `variant_recoder` translates it, from the cache alone; the tool is an annotation run under the recoder's fixed configuration, its JSON objects regrouped into one object per input line keyed by allele, each carrying the fields `--fields` names (`id`, `hgvsg`, `hgvsc`, `hgvsp` and `spdi` by default) and the ones `--vcf_string`, `--var_synonyms`, `--mane_select` and `--ga4gh_vrs` add |
| `src/recode_cli.rs`         | The command line of `vep recode`: `variant_recoder`'s options less the ones that name an Ensembl database connection or a RefSeq alignment file, which are refused by name; an option is read after one dash or two, in any case, abbreviated to any prefix that names one option, and every error is one `ERROR: [<component>] <message>` line and exit 1 |
| `src/dylib_plugins.rs`      | Loads the `--plugin` names no built-in claimed as C-ABI dynamic libraries and runs them over a batch, merging their JSON results into each consequence's `plugin_data`, or the variant's on an intergenic row |
| `src/glob.rs`               | `expand`: the `--dir_plugins` value (else `<dir>/Plugins`) expanded as Ensembl VEP's `glob` expands its plugins directory (`File::Glob::bsd_glob` under the `glob` operator's default flags): a backslash quotes the byte after it, the first brace group splits the pattern into alternatives, and each alternative has its leading `~` replaced, is walked against the file system a segment at a time and has its matches sorted; the value is one pattern, spaces included |
| `src/bin/bench_plugins.rs`  | The `bench_plugins` binary, built under the `bench` feature: times the phases of an annotation run with and without plugins and writes a JSON report                                                                         |

## Command line

`vep` annotates one input per invocation; `vep filter` and `vep recode`, the companion tools
(Ensembl VEP's `filter_vep` and `variant_recoder`), read their own options
([vep filter](../../docs/cli/filter.md), [vep recode](../../docs/cli/recode.md)). The flags are
Ensembl VEP's, underscores included, so an existing VEP command line parses; a flag the help
marks "Accepted without effect" parses and changes no annotation. The groups below name the
flags that shape a run; the full list, with every flag's parity note, is
[docs/cli-reference.md](../../docs/cli-reference.md).

| Group     | Flags                                                                                                                                                                                                                                                      |
| --------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Input     | `-i`/`--input_file` (default `STDIN`, in any case; a gzip or BGZF file is decompressed whatever its name) or `--input_data`/`--id`, the records as a string; `--format` vcf, ensembl, hgvs, spdi, region, id (a list of variant identifiers, resolved from the cache's variation files), or `guess`, which reads the first data line as Ensembl VEP does and which standard input cannot use                                                                                                                                   |
| Output    | `-o`/`--output_file` (default `variant_effect_output.txt`, `STDOUT` for standard output); `--output_format` vep, vcf, json, tab or parquet, with `--vcf`, `--json` and `--tab` as shorthands; `--force_overwrite`; `--no_headers`; `--vcf_info_field` (default `CSQ`) |
| Cache     | `--dir_cache` (else `--dir`), `--species`, `--cache_version`, `--refseq` and `--merged` select the cache in Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout; `--json_cache <dir>` (alias `--full_cache_dir`) names the directory itself; `--show_cache_info` prints its source versions (and, under `--regulatory_gff`, the first of its `file`, `emars` and `motifs` files as `regulatory_gff`); `--no_fasta` ignores a FASTA inside it; `--synonyms <file>` replaces the cache's `chr_synonyms.txt`, through which an input chromosome name reaches the cache's |
| Reference | `--fasta`/`--fa`, an indexed FASTA; `--hgvs`, `--hgvsg` and `--shift_3prime 1`, which require it; `--assembly`, read by the assembly-sensitive plugins                                                                                                                              |
| Runtime   | `--fork`/`--threads` (0 = every logical CPU up to 32; on an inflate-bound indexed input, the logical CPUs less the decompression threads, up to 32 and never fewer than half the capped CPU count); `--decompression_threads` (0 = up to 10, for BGZF input; `auto`, the default on an aarch64 host, reads an indexed BGZF input's record width); `--buffer_size` (default 5000 variants per batch); `--max_sv_size` (-1 lifts the limit) |
| Fields    | `--everything` and the field flags it implies: `--symbol`, `--canonical`, `--mane`, `--hgvs`, `--sift`, `--polyphen`, `--numbers`, `--domains`, `--variant_class`, `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`, `--check_existing` and the rest |
| Picking   | `--pick`, `--pick_allele`, `--per_gene`, `--pick_allele_gene`, `--most_severe`, `--summary`, `--flag_pick`, `--flag_pick_allele`, `--flag_pick_allele_gene`                                                                                                 |
| Plugins   | `--plugin Name,key=value,...`, repeatable; `--dir_plugins`, where dynamic-library plugins are found; `--custom <file>,<name>,<format>,...`, repeatable, a tabix-indexed VCF, BED or GFF source, or a bigWig, read as Ensembl VEP reads it ([Custom annotation](../../docs/custom-annotation.md))                                                                                                                                                         |
| Parsing   | `--dont_skip` keeps a record that fails a check; `--allow_non_variant` keeps reference-only VCF records, which the VCF output carries as their line and the JSON output as their object while the default and tab outputs write no row for them; `--chr`, `--check_ref` and `--lookup_ref` filter and check records as Ensembl VEP does                                 |
| Verbosity | `-v`/`--verbose` (debug logging); `-q`/`--quiet` (no status or progress lines and no standard-error echo of a warning; a plugin's transcript-set warning under `--refseq` or `--merged`, and dbNSFP's column warning and per-row `ERROR:` lines for a malformed transcript list, are printed whatever `--quiet` says) |

`--regulatory` annotates the regulatory and motif features of a cache converted with them, and
warns once when the cache holds none. `--regulatory_gff file=<gff3>,motifs=<gff3>,emars=<gff>`
reads the Ensembl regulatory annotation files in place of the cache's regulatory data, and
`--extended_promoters` widens their promoters to their extended bounds
([Regulatory and motif features](../../docs/regulatory.md)).

```bash
vep -i input.vcf -o output.txt --json_cache /path/to/json_cache --offline --everything
```

`bench_plugins`, the crate's second binary, takes `--input`, `--json-cache`, repeated `--plugin`
arguments in the `vep --plugin` format, `--variants`, `--runs` and `--output-json`, and reports
per-phase timings as JSON. It is a development tool behind the default-off `bench` feature, so
`cargo install vep-cli` installs `vep` alone.

## Public API

The re-exports in `lib.rs` are the embedding surface; the `pub` submodules stay so for the
in-tree binaries and tests (the rest are `pub(crate)`), and their contents may move between
releases.

| Item                                                                                        | Use                                                                                                                                                                                                                 |
| ------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Config`, `Config::from_argv`                                                               | The resolved run configuration; `config::resolve_fork` applies the same automatic thread count to an embedder's own setting                                                                                        |
| `Annotator`                                                                                 | One loaded cache plus a scoped thread pool, reused across batches: `from_config` loads every chromosome up front (the regulatory tree whole under `--regulatory`, so a corrupt shard fails construction), `annotate_batch` annotates `InputVariant`s in place, `try_annotate_rows` takes `RawRow`s |
| `AnnotationResources`, `PredictionConfig`, `build_thread_pool`                              | The shared transcripts, variation data, regulatory source, effects configuration and plugins behind an `Annotator`, the SIFT and PolyPhen output modes and the HumDiv matrix choice (`--humdiv`) it reads from `PredictionConfig`, and the pool it runs on |
| `AnnotationStats`, `AnnotationError`, `RawRow`, `validate_raw_row`, `row_to_input_variant`  | The per-batch counts, the row-level error taxonomy, and the simple-allele input (chromosome, 1-based position, reference allele, alternate allele) with its validation and conversion                               |
| `apply_to_record`, `apply_to_batch`, `FilterConfig`, `RowFilters`, `PickCriterion` | The record-level steps over the annotated alleles of one input record: the `--coding_only` and `--no_intergenic` row filters, then VEP's pick algorithm and summaries |
| `LazyTranscriptIndexes`, `MaterializedTranscriptIndexes`, `TranscriptSet` | The per-chromosome transcript index, the guard that an `Annotator` requires before fanning out to the pool, and the transcript set (every transcript, or a GENCODE set) the indexes keep |
| `InputVariant`, `Consequence`, `Impact`, `TranscriptConsequence`                            | Re-exported from `vep-core`, so a caller constructs and inspects annotations without a direct dependency                                                                                                             |
| `JsonOutputFormatter`                                                                       | Re-exported from `vep-io` for an embedder that renders JSON                                                                                                                                                         |
| `runner::Runner`                                                                            | The file-to-file pipeline the `vep` binary drives: `Runner::new(config)` resolves the cache, then `.run()`                                                                                                           |

## Run it

```bash
cargo build -p vep-cli                     # target/debug/vep
cargo test -p vep-cli                      # unit tests plus the integration tests under tests/
cargo run -p vep-cli --bin vep -- --help   # the full flag list
cargo build -p vep-cli --features bench --bin bench_plugins   # target/debug/bench_plugins
```

`cargo build -p vep-cli --no-default-features --bin vep` builds `vep` on the system allocator
instead of mimalloc and without the Parquet writer, so `--output_format parquet` is refused;
adding `--features parquet` keeps the writer. Each integration test under `tests/` covers one contract; among them,
`integration_test.rs` runs the built binary as a subprocess; `golden.rs` and `format_parity.rs`
compare every output format on the golden corpora with the Ensembl VEP files committed beside
them; `pipeline_parity.rs` checks that the output bytes do not change with the batch size, the
thread count, the pipeline mode, depth and task size or the writer strategy, and that the
warnings and skipped-variants files list skipped records and the warnings a record draws in
input order in both modes; `parity_matrix.rs` checks that every long flag has a row in
`tests/parity/matrix.tsv` at the repository root; `dylib_plugins.rs` exercises the
dynamic-library plugin path; `regulatory.rs` checks the regulatory and motif rows, the picks
and row filters over them and the cache shapes the loader refuses, against Ensembl VEP's rows
on a slice of a converted cache, and the rows `--regulatory_gff` writes from the GFF files
under `tests/fixtures/regulatory-gff/`. `cargo bench -p vep-cli --bench cache_load_bench` times cache
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
