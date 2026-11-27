# vep-intermediates-compare

Joins a trace of Ensembl VEP's intermediate values with vep-rs's default-format output of the same
variants and reports, row by row, which values diverge. The trace is the tab-separated file
`scripts/concordance/trace_perl_intermediates.pl` prints inside the Ensembl VEP container: one row
per variant, allele and transcript with the coding and translation coordinates, the reference and
alternative codon and peptide, and the consequence terms. It is the tool for a small VCF of
discordant variants: it shows whether a disagreement starts at the coordinates, the codons, the
peptides or only the consequence.

```text
vep-intermediates-compare --perl-trace TRACE.tsv --rust-output OUTPUT.txt [--discordant-only] [--output REPORT]
```

`--discordant-only` prints only the rows with a diverging value. `--output` names the report file;
the default `-` is standard output.

The command line is read as CPython's `argparse` reads it. An option may be shortened to any prefix
that fits it alone (`--perl`, `--disc`), takes its value as the next argument or after `=`, and a
later value replaces an earlier one. A value that starts with `-` is taken only when it is `-`
alone, starts like a negative number in ASCII digits or holds a space; after `=` any value is taken.
`-h` and `--help` print the usage and the options.

## Contract

**Inputs.** Both files are read as text: `\n`, `\r\n` and a lone `\r` end a line, as universal
newlines read them, and the bytes are UTF-8 whatever the locale. A byte sequence that is not UTF-8
anywhere in either file, a comment line included, refuses the run; a gzip file is not gunzipped
and is refused the same way. Each input path goes through `str(pathlib.Path(arg))`; the report
path is used as given.

The trace is read as `csv.DictReader(file, delimiter="\t")` reads it. A field that starts with `"`
runs to its closing quote across tabs and line ends, `""` inside it is one quote, and text after
the closing quote is kept. The first record is the header and the columns are found by name, the
last one where a name stands twice; a blank line is skipped. The three key columns are `location`,
`allele` and `transcript_id`: a header without one of them refuses the run at its first row. The
value columns are `cds_start`, `cds_end`, `ref_codon`, `alt_codon`, `ref_peptide`, `alt_peptide`
and `consequence_set`: one the header lacks reads as the empty string. A field longer than 131,072
code points, `csv.field_size_limit()`, refuses the run.

Of the output, a line that starts with `#` and a row of fewer than 14 tab-separated columns are
skipped. Location, Allele and Feature key a row; Consequence, CDS_position, Amino_acids and Codons
are its values. No field is stripped.

**Join.** Under one key the last row of each side counts, and the keys are ordered as Python
orders tuples of `str`, by code point, a shorter string before any string it starts. `Perl rows`
and `Rust rows` are the keys of each side, `Unique keys` the keys of both together.

**Comparison**, for every key both sides hold:

- `consequence`: both texts are split on `,` and `&`, empty terms dropped, and compared as sets.
- `cds_start`: the trace's value against the text before the first `-` of CDS_position, when both
  are present and CDS_position is not `-`.
- `ref_codon` and `alt_codon`: against the two sides of Codons under `str.upper()`, the full
  Unicode mapping of the toolchain's tables, when the trace has a reference codon and Codons holds
  exactly one `/`.
- `ref_peptide` and `alt_peptide`: against the two sides of Amino_acids, exactly, when the trace
  has a reference peptide and Amino_acids holds exactly one `/`.

**Report.** Three counts and a blank line; one line per compared key in key order, `DISCORD` or
`match`, the key, the diverging fields in the order above and both sides' values; then the summary:
the four totals, each diverging field with its count and its share of the discordant rows as
`f"{x:.1f}"`, and the 20 most frequent consequence disagreements. Both lists keep
`Counter.most_common`'s order: the larger count first and, between equal counts, the one counted
first. The labels `Perl` and `Rust` in the report and in the option names stand for the trace and
the output.

**A trace row shorter than its header** holds `None` in the columns it does not reach, as
`DictReader` fills them. Such a value prints as `None`, and an absent `alt_peptide` diverges. Three
shapes end the report early with exit code 1 and one line on standard error, the report written so
far kept:

- two keys that differ first in a field one of them lacks have no order, as `sorted()` has none
  between `None` and a `str`: the report ends after its three counts;
- a compared row without `consequence_set`: the report ends before that row;
- a compared row with a reference codon and without `alt_codon`, where Codons holds one `/`: the
  same.

**Exit codes.** 0 for a complete report. 2 for a usage error, with the usage on standard error. 1
for any other failure, reported as `ERROR: [vep-intermediates-compare] <detail>` on standard
error: an input that cannot be read or is refused above, a report path that cannot be written, a
`sort` failure, or a report that ends early. An input failure writes nothing and leaves an
existing report file untouched, since both inputs are read to their ends before the report file
is opened; the report may therefore be written over one of the inputs.

**Memory.** Each side is read as a stream and reduced to its key and value columns, and both
sides are ordered at once by a child `sort -S 1G` each, under `LC_ALL=C` with `LC_COLLATE`
removed: a side larger than the buffer is sorted in runs that `sort` keeps under `TMPDIR` and
merges. The merge of the two sides holds one row of each. The row lines wait in a file until the
counts that open the report are known; that file and the sorted sides live in a directory under
`TMPDIR` that is removed on exit. The resident memory is therefore bounded by the two `sort`
buffers and one entry per distinct consequence disagreement, whatever the size of the inputs. A
host without `sort` orders each side in memory.

## Layout

- `src/cli.rs`: the command line as `argparse` reads it, the usage line and the help text.
- `src/input.rs`: lines under universal newlines, the UTF-8 check, the `csv` reader and the header
  lookup of the trace, the 14-column rows of the output, the `sort` child.
- `src/rowcodec.rs`: the row each side hands to `sort`, encoded so that bytewise order is the
  order of the keys as tuples of strings.
- `src/compare.rs`: the sorted merge, the last row under a key, the comparison and the report's
  lines.
- `src/main.rs`: the two sides in parallel, the temporary directory, the report and the exit code.
- `tests/golden.rs` and `tests/fixtures/`: synthetic inputs that reach each branch (`trace.tsv`
  and `output.txt` for the comparison, the `*_trace.tsv` files for headers in another order, short
  rows and refused traces) and the report each run must write, byte for byte, as
  `expected_<run>.txt`. `cargo test` runs the binary over the fixtures and compares standard
  output, the report file and the exit code.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-intermediates-compare/Cargo.toml
cargo test --locked --manifest-path tools/vep-intermediates-compare/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-intermediates-compare/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-intermediates-compare/Cargo.toml -- --check
```

A run after the tracer, from the repository root:

```bash
tools/vep-intermediates-compare/target/release/vep-intermediates-compare \
    --perl-trace perl_intermediates.tsv --rust-output rust_output.txt --discordant-only
```
