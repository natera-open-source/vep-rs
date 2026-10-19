# vep-io

Input parsers and output formatters for vep-rs, the Rust reimplementation of the
Ensembl Variant Effect Predictor (VEP). The crate turns one line of Ensembl
default, region, SPDI (Sequence Position Deletion Insertion) or genomic Human
Genome Variation Society (HGVS) input into `vep_core::variant::InputVariant`
values, detects an input's format from its
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
`fields::field_value` answers every column name from the row, reading plugin data
for a name the table lacks. The tab formatter prints each column with `-` for an
absent value; the VCF formatter folds every row into a `|`-joined CSQ entry (the
INFO key VEP uses for its annotations), joins the entries with `,` and appends
`CSQ=...` to the record's own INFO column; the JSON formatter buckets the rows into
`transcript_consequences` and its sibling arrays under VEP's lowercased, renamed
keys; the Parquet formatter writes one tab-separated row per (allele, consequence)
with `chrom`, `pos`, `end`, `ref` and `alt` ahead of the `--tab` columns, for the
`duckdb` subprocess the command-line runner spawns at end of run.

The `vep` binary reads VCF through its own line parser in `vep-cli`, which keeps
the input's header lines for VCF output passthrough; this crate has no VCF
parser of its own. The Ensembl, region, SPDI and HGVS parsers are the ones the
binary runs. VEP's default text format has no
formatter struct here: `vep-cli` renders its rows from this crate's field table
(`record_rows_into`, `extra_keys_with`, `field_value`) and takes its header from
`default_format_header`.

## Modules

| Module               | Role                                                                                                                                                                                                                                   |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `error`              | `IoError`: `Io`, `VcfParse` (a parse failure of the Ensembl, region or HGVS parser), `SkippedLine` (a line Ensembl VEP's parser warns about and reads past, the warning's text as the message), `NotImplemented` (transcript-level HGVS) and `Unsupported` (input refused by design). |
| `input`              | The `InputParser` trait (`next_variant` yields one `InputVariant` per ALT allele and `None` at end of input; `take_empty_lines` the empty lines read past), `Delimiter` (the column delimiter as `--delimiter` names it or detection chooses it), `ParserContext` (delimiter, reference sequence, synonyms, the annotation source's chromosomes) and `minimise_alleles`. |
| `input::detect`      | `detect_format`: Ensembl VEP's format detection from the first data line, returning the verdict, the delimiter it split the line on and the bytes read to reach it.                                                                 |
| `input::ensembl`     | `EnsemblParser`: VEP's four- to six-column default format split on the context's delimiter; strand tokens, class words, upper-casing and the minimisation of a bi-allelic indel as Ensembl VEP reads them; one variant per ALT of the allele string. |
| `input::region`      | `RegionParser`: `chr:start-end:strand/allele` or `chr:start-end/allele`; the reference allele is `-` and the strand defaults to forward.                                                                                              |
| `input::spdi`        | `SpdiParser`: `sequence:position:deletion:insertion`, the position 0-based, a deletion or insertion written as bases or as a length; the reference allele is checked against or read from the context's reference, the sequence name resolves to the source's chromosome through its synonyms, and a refused notation is a `SkippedLine` carrying Ensembl VEP's text. |
| `input::hgvs`        | `HgvsParser`: genomic `:g.` substitutions, deletions, insertions, delins and duplications (a duplication is read as an insertion of the duplicated bases); `:c.`, `:p.`, `:n.` and `:r.` expressions fail with `NotImplemented`.       |
| `output`             | The `OutputFormatter` trait: `write_header`, `format_variant`, `format_record` (every allele split from one input record together) and `finish`.                                                                                       |
| `output::fields`     | The field table every text format shares: the thirteen default columns, the flag-to-field order, header descriptions, `record_rows` ordering, `field_value`, the `Extra` column keys, the default and tab headers and the VCF `##VEP=` meta line. |
| `output::tab`        | `TabOutputFormatter`: the default columns, then the active flag fields, then plugin fields, each as its own column, `-` where a value is absent.                                                                                        |
| `output::vcf_output` | `VcfOutputFormatter`: the input line with one CSQ entry per row, VEP's escaping, an earlier CSQ entry replaced, and the header block built from the input's own header lines.                                                           |
| `output::json`       | `JsonOutputFormatter`: one JSON Lines object per record, keys lowercased and renamed as VEP renames them, numeric-looking values emitted as numbers except for identifier keys.                                                         |
| `output::parquet`    | `ParquetOutputFormatter`: the tab-separated intermediate the command-line runner hands to `duckdb`; a `-` becomes an empty cell except in `Allele`, `ref` and `alt`, where a dash is a real allele.                                     |

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
| `EnsemblParser::new(reader, &context)`, `SpdiParser::new(reader, &context)`, `RegionParser::new`, `HgvsParser::new` | Line parsers over any `BufRead`; blank lines are skipped, and so are `#` lines except by the SPDI parser, which reads them as Ensembl VEP's does. |
| `Delimiter::from_flag`, `ParserContext`         | The `--delimiter` value as a pattern (a tab, a space, a run of spaces, `\s+` or literal text; a regular expression is refused); what the parsers share.         |

Output:

| Item                                                    | Phrase                                                                                                                                              |
| ------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| `OutputFormatter`                                       | `write_header`, `format_variant`, `format_record`, `finish`; the default `format_record` formats each allele on its own.                            |
| `fields::FieldOptions`                                  | One `bool` per field group, mirroring the command-line flags, plus `no_escape` (HGVSp keeps `=` instead of `%3D`).                                  |
| `fields::ExtraFieldsPlan::new(options)`                 | The row-invariant part of the `Extra` column, computed once per run.                                                                                |
| `fields::flag_fields`, `fields::field_description`      | The active optional fields in VEP's order; VEP's description of a field for the header blocks.                                                      |
| `fields::record_rows`, `fields::record_rows_into`       | The rows one input record produces, in VEP's order; the `_into` form reuses a caller-owned vector.                                                   |
| `fields::field_value(field, variant, tc, options, sep)` | The value of one named field for one row, `sep` joining list values (`,` in text formats, `&` inside a CSQ entry).                                  |
| `fields::extra_keys`, `fields::extra_keys_with`         | The `Extra` keys of a row: active flag fields first, then overlap and plugin keys alphabetically.                                                    |
| `fields::default_format_header`, `fields::tab_header`   | The `##` header lines and the `#`-prefixed column line of the default and tab formats.                                                               |
| `fields::tab_columns`, `fields::vcf_vep_meta_line`      | The tab format's column list; the `##VEP=` meta line with version, time, cache directory and source versions.                                       |
| `TabOutputFormatter::new(options, plugin_fields, header_lines, no_headers)` | `columns()` returns the output order.                                                                                            |
| `VcfOutputFormatter::new(info_field_name, options, plugin_fields)` | `with_run_info` sets the cache and command line the meta lines carry; `header_lines(input_headers)` builds the header block.               |
| `JsonOutputFormatter::new(options, plugin_fields, assembly)` | `format_record_value` builds the `serde_json::Value` without a string round trip; `format_variant_value` treats one allele as a whole record.    |
| `ParquetOutputFormatter::new(options, plugin_fields)`   | `header_columns()` is `KEY_COLUMNS` then the tab columns; the Parquet step itself runs in the command-line runner.                                   |
| `IoError`                                               | The crate's error type, re-exported at the root with `InputParser` and `OutputFormatter`.                                                           |

The flags behind `FieldOptions` are documented in
[docs/cli-reference.md](../../docs/cli-reference.md#output-control); the input
formats and their selection in its [Input](../../docs/cli-reference.md#input)
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
