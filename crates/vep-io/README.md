# vep-io

Input parsers and output formatters for vep-rs, the Rust reimplementation of the
Ensembl Variant Effect Predictor (VEP). The crate turns one line of Ensembl
default, region, SPDI (Sequence Position Deletion Insertion), genomic,
transcript or protein Human Genome Variation Society (HGVS) or variant identifier
input into
`vep_core::variant::InputVariant` values, detects an input's format from its
first data line, and renders annotated variants as VEP's `--tab`, Variant Call
Format (VCF), JSON (JavaScript Object Notation) and Parquet outputs through one
shared field table. `vep-cli` is its only dependent in the workspace, and
`vep-core` is the only sibling crate it depends on.

## How it works

Take the Ensembl-format line `21 25000100 25000100 A/G/T 1 rs123` (one record,
two ALT alleles). `EnsemblParser::next_variant` splits the columns on the
context's delimiter, upper-cases the allele string, reads the strand (`-` or
`-1` is reverse, any other token or none forward) and yields two `InputVariant`
values in ALT order, `allele_index` 0 and 1, both marked `minimised` because the
site is multi-allelic, both carrying `rs123` as the record's name. A bi-allelic
indel written un-minimised (`21 100 101 AC/A`) is trimmed to `C/-` at 101 as
Ensembl VEP trims it, and an allele column without a `/` (`DUP`) is a
structural variant over the line's span. A VCF record follows the same path
through the `vep-cli` crate's line parser. Annotation happens outside this
crate, and the two variants return as one record slice.
`OutputFormatter::format_record` receives that slice: `fields::record_rows`
orders the rows transcript-major (by transcript id) and allele-minor (ALT order),
then adds one intergenic row per allele that overlapped nothing, and
`fields::write_field_value` answers every column name from the row as bytes
(`fields::field_value` collects them into a `String`), reading plugin data for a
name the table lacks. The tab formatter prints each column with `-` for an
absent value, as lines through `format_record` or as bytes through
`tab::render_tab_record`; the VCF formatter folds every row into a `|`-joined CSQ entry (the
INFO key VEP uses for its annotations), joins the entries with `,` and appends
`CSQ=...` to the record's own INFO column, a VCF input's line kept as the bytes
it was read from and any other input's line built as Ensembl VEP's
`to_VCF_record` builds it (`output::vcf_output::to_vcf_record`: the anchor base
before an indel from the run's reference or `N`, reverse-strand alleles
reverse-complemented, a structural class as `<DEL>` with `END=`); the JSON formatter buckets the rows into
`transcript_consequences` and its sibling arrays under VEP's lowercased, renamed
keys; the Parquet formatter writes one tab-separated row per (allele, consequence)
with `chrom`, `pos`, `end`, `ref` and `alt` ahead of the `--tab` columns, for the
`duckdb` subprocess the command-line runner spawns at end of run.

The `vep` binary reads VCF through its own line parser in `vep-cli`, which keeps
the input's header lines for VCF output passthrough; this crate has no VCF
parser of its own. The Ensembl, region, SPDI, HGVS and identifier parsers are
the ones the binary runs; the identifier parser yields the mappings its caller
resolved from the cache's variation files (`vep-cli`'s `variation_name_index`),
one record per mapping. VEP's default text format has no
formatter struct here: `vep-cli` renders its rows from this crate's field table
(`record_rows_into`, `extra_keys_with`, `write_field_value`, `write_terms`,
`write_versioned_id`) and takes its header from `default_format_header`.

## Modules

| Module               | Role                                                                                                                                                                                                                                   |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `error`              | `IoError`: `Io`, `VcfParse` (a parse failure of the Ensembl parser, and the JSON formatter's serialisation failure), `SkippedLine` (a line Ensembl VEP's parser warns about and reads past, the warning's text as the message), `Unsupported` (input refused by design, the reason as the whole message; the `:r.` HGVS notation) and `CacheUnreadable` (a transcript cache an input names a record in but a lookup cannot read, the cache's own refusal as the message; the run stops and prints it under the `cache` component). |
| `input`              | The `InputParser` trait (`next_variant` yields one `InputVariant` per ALT allele and `None` at end of input; `take_empty_lines` the empty lines read past; `take_warnings` the warnings the parser raised, in order; `record_line_number` the line number of the record last yielded, from a parser whose format counts its lines other than by records read), `Delimiter` (the column delimiter as `--delimiter` names it or detection chooses it), `ParserContext` (delimiter, the run's reference through the `ReferenceSequence` trait, synonyms, `--max_sv_size`, the annotation source's chromosomes and its transcripts through the `TranscriptSource` trait), `minimise_alleles` and `empty_lines_warning` (the text `skip_empty_lines` warns a run of empty lines with). |
| `input::detect`      | `detect_format`: Ensembl VEP's format detection from the first data line, returning the verdict, the delimiter it split the line on and the bytes read to reach it.                                                                 |
| `input::ensembl`     | `EnsemblParser`: VEP's four- to six-column default format split on the context's delimiter; strand tokens, class words, upper-casing and the minimisation of a bi-allelic indel as Ensembl VEP reads them; one variant per ALT of the allele string. |
| `input::region`      | `RegionParser`: `chr:start-end[:strand]/allele`, the allele after a `/` or a `:`, matched whole against Ensembl VEP's region expression; the strand is `:1`, `:+1` or `:-1`; the reference allele is `N` per base of the span for the run's reference lookup, `-` for an insertion; a word outside `ACGTN-` names a structural class; a line outside the expression or a class the table lacks is a `SkippedLine`. |
| `input::spdi`        | `SpdiParser`: `sequence:position:deletion:insertion`, the position 0-based, a deletion or insertion written as bases or as a length; the reference allele is checked against or read from the context's reference, the sequence name resolves to the source's chromosome through its synonyms, and a refused notation is a `SkippedLine` carrying Ensembl VEP's text, as is a record whose minimised span exceeds `--max_sv_size`. |
| `input::hgvs`        | `HgvsParser`: HGVS notation resolved against the context's reference sequence as Ensembl VEP resolves it through its database; `input::hgvs::genomic` reads the `g.` and `m.` shapes (substitution, `del`, `dup`, `delins`, `ins`, `inv`, `=`) with the genome's bases, the reference check and `--lookup_ref`; `input::hgvs::transcript` reads `c.` and `n.` on the transcript the context's `TranscriptSource` finds by stable id or, through Ensembl VEP's pick order, by gene name, placing the position through the transcript's exon mapper and keeping its strand; `input::hgvs::protein` reads `p.` on the transcript found by translation id, by transcript id or by gene name (every transcript of the gene tried in turn), turning the residue change into the one shortest nucleotide path to it (a substitution or `delins`) or into the deletion of its codons; `r.` fails with `Unsupported`, and a notation Ensembl VEP's resolver refuses or `validate_vf` skips is a `SkippedLine` carrying its text. |
| `input::id`          | `IdParser`: one variant identifier per line, whitespace removed, the list read whole; `names` for the caller to resolve and `supply_mappings` to hand back one `NameMapping` list per name; a name without a mapping is a `No variant found` warning and a run of empty lines `skip_empty_lines`'s, both through `take_warnings`, the empty lines numbered as Ensembl VEP counts the list's lines. |
| `output`             | The `OutputFormatter` trait: `write_header`, `format_variant`, `format_record` (every allele split from one input record together), `format_record_bytes` (the same lines as bytes, for a format that writes input bytes back) and `finish`. |
| `output::fields`     | The field table every text format shares: the thirteen default columns, the flag-to-field order, header descriptions, `record_rows` ordering, `write_field_value` and its `String` collector `field_value`, the `Extra` column keys, the default and tab headers and the VCF `##VEP=` meta line. |
| `output::frequency`  | The allele frequency fields of a row (`AF`, the 1000 Genomes and gnomAD population columns, `MAX_AF` and `MAX_AF_POPS`), filled from a record's co-located variants as Ensembl VEP fills them.                                      |
| `output::tab`        | `TabOutputFormatter` and `render_tab_record`: the default columns, then the active flag fields, then plugin fields, each as its own column, `-` where a value is absent.                                                                                        |
| `output::vcf_output` | `VcfOutputFormatter`: the input line, or the line `to_vcf_record` builds for a non-VCF input, with one CSQ entry per row, VEP's escaping, an earlier CSQ entry replaced, and the header block built from the input's own header lines. |
| `output::json`       | `JsonOutputFormatter`: one JSON Lines object per record, keys lowercased and renamed as VEP renames them, numeric-looking values emitted as numbers except for identifier keys.                                                         |
| `output::parquet`    | `ParquetOutputFormatter`: the tab-separated intermediate the command-line runner hands to `duckdb`; a `-` becomes an empty cell except in `ref`, `alt` and the allele columns `fields::dash_is_a_value` names, where a dash is a real allele. |

The field groups that `--symbol`, `--hgvs`, `--af` and the other boolean flags
switch on are listed once, in VEP's order, in `output::fields`; a group's fields
appear only when its flag is set, and the `user` group (`IMPACT`, `DISTANCE`,
`STRAND`, `FLAGS`) is always on. The format specifications, with headers and
examples, are in [docs/output-formats.md](../../docs/output-formats.md).

## Public API

Input:

| Item                                           | Phrase                                                                                                                                       |
| ---------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `InputParser::next_variant`                    | The one method every parser implements.                                                                                                      |
| `detect::detect_format(reader)`                | The format of the first data line and the bytes read to find it, for the caller to replay.                                                   |
| `EnsemblParser::new(reader, &context)`, `SpdiParser::new(reader, &context)`, `RegionParser::new(reader, &context)`, `HgvsParser::new(reader, &context, lookup_ref)` | Line parsers over any `BufRead`; blank lines are skipped, and so are `#` lines except by the SPDI parser, which reads them as Ensembl VEP's does. |
| `IdParser::read(reader)`, `IdParser::names`, `IdParser::supply_mappings`, `IdParser::attempted_line` | The identifier list read whole, its names, the resolved mappings handed back before the first `next_variant`, and the line of the name whose mapping was last built, for a skip warning to quote. |
| `Delimiter::from_flag`, `ParserContext`         | The `--delimiter` value as a pattern (a tab, a space, a run of spaces, `\s+`, literal text, or a regular expression split by the caller's compiled splitter); what the parsers share.         |

Output:

| Item                                                    | Phrase                                                                                                                                              |
| ------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| `OutputFormatter`                                       | `write_header`, `format_variant`, `format_record`, `format_record_bytes`, `finish`; the default `format_record` formats each allele on its own.      |
| `fields::FieldOptions`                                  | One `bool` per field group, mirroring the command-line flags, plus `no_escape` (HGVSp keeps `=` instead of `%3D`).                                  |
| `fields::ExtraFieldsPlan::new(options)`                 | The row-invariant part of the `Extra` column, computed once per run.                                                                                |
| `fields::flag_fields`, `fields::field_description`      | The active optional fields in VEP's order; VEP's description of a field for the header blocks.                                                      |
| `fields::record_rows`, `fields::record_rows_into`       | The rows one input record produces, in VEP's order; the `_into` form reuses a caller-owned vector.                                                   |
| `fields::field_value(field, variant, tc, options, sep)` | The value of one named field for one row, `sep` joining list values (`,` in text formats, `&` inside a CSQ entry).                                  |
| `fields::write_field_value(out, field, variant, tc, options, sep)` | The same value appended to an output vector, nothing when the field is absent on the row; `field_value` collects it. |
| `fields::write_terms(out, terms, vocabulary)`, `fields::write_versioned_id(out, id, version, flag)` | Consequence terms `,`-joined in the configured vocabulary; a stable id with `.<version>` appended when the flag asks for it. |
| `tab::render_tab_record(out, record, &formatter, rows)` | One record's tab lines appended as bytes, each column through `write_field_value`; `format_record` returns the same lines. |
| `fields::extra_keys`, `fields::extra_keys_with`         | The `Extra` keys of a row: active flag fields first, then overlap and plugin keys alphabetically.                                                    |
| `fields::default_format_header`, `fields::tab_header`   | The `##` header lines and the `#`-prefixed column line of the default and tab formats.                                                               |
| `fields::tab_columns`, `fields::vcf_vep_meta_line`      | The tab format's column list; the `##VEP=` meta line with version, time, cache directory and source versions.                                       |
| `TabOutputFormatter::new(options, plugin_fields, selected_fields, header_lines, no_headers)` | `selected_fields` is the `--fields` list, which replaces the columns; `columns()` returns the output order.                      |
| `VcfOutputFormatter::new(info_field_name, options, plugin_headers, selected_fields)` | `plugin_headers` pairs each plugin field with its header description; `with_run_info` sets the cache and command line the meta lines carry; `header_lines(input_headers)` builds the header block. |
| `JsonOutputFormatter::new(options, plugin_fields, assembly)` | `format_record_value` builds the `serde_json::Value` without a string round trip; `format_variant_value` treats one allele as a whole record.    |
| `ParquetOutputFormatter::new(options, plugin_fields)`   | `header_columns()` is `KEY_COLUMNS` then the tab columns; the Parquet step itself runs in the command-line runner.                                   |
| `IoError`                                               | The crate's error type, re-exported at the root with `InputParser` and `OutputFormatter`.                                                           |

The flags behind `FieldOptions` are documented in
[docs/cli/output-control.md](../../docs/cli/output-control.md); the input
formats and their selection in the CLI reference's [Input](../../docs/cli/input-output.md#input)
section.

## Run it

```bash
cargo build -p vep-io
cargo test -p vep-io
cargo doc --no-deps -p vep-io
cargo bench -p vep-io --bench json_format_bench
```

The unit tests live beside the code in each module. The benchmark compares
`format_variant_value` (the embedding path, a direct `Value`), `format_variant`
(the command-line path, one JSON line) and a synthetic serialise-then-parse
round trip; a healthy build keeps the round-trip lane materially slower than the
direct lane. The output parity tests against Ensembl VEP's golden corpora run from
`vep-cli` (`cargo test -p vep-cli --test golden --test format_parity`). The
repository's standard commands are in [CONTRIBUTING.md](../../CONTRIBUTING.md).

## Related

- [vep-core](../vep-core/README.md): `InputVariant`, `TranscriptConsequence`,
  coordinate parsing, structural-variant routing and record construction, the
  cache info the headers print.
- [vep-cli](../vep-cli/README.md): the `vep` binary, which selects the parser
  and formatter from its flags, renders the default format from this crate's
  field table and runs the Parquet `duckdb` step.
- [docs/output-formats.md](../../docs/output-formats.md): the five output formats,
  their headers, row order and parity with Ensembl VEP.
- [docs/cli-reference.md](../../docs/cli-reference.md): every flag, including the
  input and output selection.
- [docs/intended-divergences.md](../../docs/intended-divergences.md): where the
  output differs from Ensembl VEP on purpose.
- [Root README](../../README.md): installation, quick start and the published
  concordance figures.
