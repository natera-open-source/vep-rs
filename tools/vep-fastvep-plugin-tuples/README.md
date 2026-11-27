# vep-fastvep-plugin-tuples

The supplementary columns of fastVEP's tab output as plugin tuples. A fastVEP run with `--sa-dir`
adds one tab column per loaded source after the built-in columns. This tool writes those columns as
plugin tuples, the form in which a comparator's values are scored against Ensembl VEP's plugin and
`--custom` fields: tab-separated, five columns, no header: the row's Location, Allele and Feature,
the Ensembl VEP field name, the decoded value. It reads the tab output as fastVEP writes it or as
the harness normalizes it: the two carry the same supplementary columns, and the input is never
changed.

```text
vep-fastvep-plugin-tuples <fastvep.tab> <tuples.tsv> [--field-map <map.tsv>] [--fields <name,name>]
```

- `<fastvep.tab>`: the fastVEP tab output.
- `<tuples.tsv>`: the tuple file to write.
- `--field-map`: the field map (default: `data/fastvep_field_map.tsv` of this tool, at the path the
  tool was built from; a binary run away from its source tree names its map).
- `--fields`: the Ensembl VEP field names to keep, joined by commas, for a run whose database or
  reference holds fewer fields than the map.

The command line is read as CPython's `argparse` reads it. An option may be shortened to any prefix
that fits it alone (`--field-m`; `--field` fits both and is refused), takes its value as the next
argument or after `=`, and a later value replaces an earlier one. A value that starts with `-` is
taken only when it is `-` alone, starts like a negative number or holds a space. The two paths may
stand before, between or after the options, and every argument after the first `--` is a path.
`-h` and `--help` print the usage and the options.

## Contract

**Lines and the column header.** A line ends at its line feed, with the carriage return before it
where the file has one. Any other carriage return is part of its cell and is written `%0D`, in the
Location, Allele and Feature cells too. The column header is the last line before the first row
that starts with one `#`; a line that starts with `##` is never the header. Every line after the
header that does not start with `#` is a row, an empty line included. A row of fewer than five
cells writes nothing, and a row that ends before the last mapped column is counted as too short.

**Cells, entries and parts.** A cell of `-` and an empty cell hold nothing. Any other cell holds
one entry, or several joined with `,`. An entry is a `|`-joined string whose first part, the
alternate allele or the gene symbol, is dropped. The field map names the Ensembl VEP field of each
part after it, one row per field:

| Column           | Holds                                                                                |
| ---------------- | ------------------------------------------------------------------------------------ |
| `fastvep_column` | header name of the supplementary column                                              |
| `part`           | `N`, counted from 1 after the dropped part; `N-M` for parts N to M joined with `\|`; `N.word` or `N.score` for one half of a part written `word(score)` |
| `vep_field`      | Ensembl VEP field name the value is written under                                    |
| `source`         | plugin or `--custom` source the field belongs to                                     |
| `kind`           | `text` or `number`, the value rule a scorer applies                                  |
| `note`           | free text                                                                            |
| `list_separator` | `-`, or the separator Ensembl VEP writes where the column joins a list with `&`      |
| `reference_set`  | `-`, or the expected set of the source's committed corpus whose reference declares the field |

A `word(score)` part splits at its last `(` and must end with `)`. A range whose end lies past the
entry joins the parts that are there, and its end may lie far past any entry.

**Values.** A value is percent-decoded, except that a tab and every character that breaks a line
stay encoded, so a tuple stays one line. A number keeps the text fastVEP printed. An empty part
writes no tuple. A populated part without a row in the map, and the populated half of a
`word(score)` part without one, is dropped and counted. Equal values from several entries of one
cell are written once. `--fields` keeps the named fields alone.

**The summary line.** One line on stderr, `fastvep_plugin_tuples: wrote N tuples from R rows (S too
short); ...`, gives the tuples written, the rows read and the rows that end before the last mapped
column, the tuples of each kept field of the mapped columns (a field that wrote none reads 0), the
mapped columns, the dropped parts and halves by column (`COLUMN:3=2`, `COLUMN:1.word=5`, in the
order of column, part number and part), the values `--fields` left out, the parts a `word(score)`
row could not split, the cells with several entries and the `FV_` columns the map has no row for.
A cell repeated on the next row is counted again in every counter.

**Refusals.** A run ends with one line, `ERROR: [fastvep_plugin_tuples] <message>`, the name in
brackets a fixed part of the line, when `--fields` names no field, as a whole or beside a comma;
when the map cannot be read, is not UTF-8 text or is malformed (`<map>:<line>: <fault>`), or lacks
a field `--fields` names; when the input cannot be opened, has no column header line, names none of
the map's columns, names one of them more than once or lacks the column of a field `--fields`
names; and when the tuple file is the input or cannot be opened. The faults are looked for in that
order, and the tuple file is opened only after them, so a refused run leaves an earlier tuple file
as it was.

**Bytes.** Wherever a byte could be written in more than one way it follows a CPython rule, named
here by its CPython name.

- The input and the tuple file are UTF-8 whatever the locale, and a byte that is not UTF-8, raw or
  percent-encoded, is written as it is (`surrogateescape`).
- An escape is decoded as `urllib.parse.unquote_to_bytes` decodes it: `%` and two hexadecimal
  digits of either case are one byte, and any other `%` stays. The escapes are decoded among the
  value's own bytes, so the two read as one UTF-8 text.
- The characters that stay encoded are the tab and those `str.splitlines()` ends a line at (`\n`,
  `\v`, `\f`, `\r`, U+001C to U+001E, U+0085, U+2028, U+2029), each written as
  `urllib.parse.quote` writes it: `%` and two uppercase hexadecimal digits for each of its UTF-8
  bytes.
- The field map is strict UTF-8, and its lines end where `str.splitlines()` ends them. A part's
  digits are the decimal digits of any script, as `\d` of `re` and `int()` read them. A part or a
  kind in a refusal is written as `repr()` writes a `str`.
- The map is the file `pathlib.Path` reads its name as (`map.tsv/` is `map.tsv`, the empty name is
  the directory `.`), and a line names it as the command line wrote it. The input and the tuple
  file are opened by their names as written; whether they are one file is decided as
  `pathlib.Path.samefile` decides it.
- The names `--fields` lists that the map lacks are written once each in the order `sorted()`
  gives a list of `str`.
- stderr is UTF-8. A byte of an argument or of the column header that is not UTF-8 is written as
  `backslashreplace` writes the `str` that `surrogateescape` reads it as: `\udcXX`.
- The reason a file cannot be opened is `os.strerror` of the error number.

The printable code points and the decimal digits are those of Unicode 16.0.0, the database of
CPython 3.14.

**Exit codes.** 0 when the tuple file and the summary line are written. 1 for a refusal. 1, with
one line on stderr that is no refusal, for a failure after the files are opened: a read or a write
that fails ends with `OSError: [Errno N] <strerror>`, under the name of the class CPython raises
for the error number, and leaves a partial tuple file once rows are written; a part of more than
4300 digits ends with the `ValueError` line of `int()`; a thread bound that is no whole number of
one or more ends with `ERROR: [vep-fastvep-plugin-tuples] <detail>`. 2 on a command line that is
refused: the usage line, then `vep-fastvep-plugin-tuples: error: <message>`, on stderr.

**Memory and cores.** The input is read as a stream, in blocks of whole lines of about 256 KiB.
Each block's rows are written as tuples on a worker thread, and the tuples of a block are written
when every block before it is written, so the tuple file is in input order on any number of
threads. A worker takes no block further than one block a thread past the last one written, so the
blocks in memory are bounded by the threads; neither the input nor the tuple file is ever held
whole. The threads are the host's available cores, no more than a thread a block of a regular
file, or the whole number `VEP_FASTVEP_PLUGIN_TUPLES_THREADS` names; the thread count changes no
byte of the outputs.

## Layout

- `data/fastvep_field_map.tsv`: the field map of fastVEP 0.5.0's source columns.
- `src/main.rs`: the run, the order of the refusals, the exit code.
- `src/cli.rs`: the command line as `argparse` reads it, the usage line and the help text;
  `src/cli_argparse_tests.rs` holds it to `tests/fixtures/command_lines.tsv`, the table of what
  CPython's `argparse` reads from each of its command lines.
- `src/map.rs`: the field map as rules, the names `--fields` keeps, the refusals of a malformed map.
- `src/columns.rs`: the cells of a line, the header's mapped columns and what each part of one is
  read as.
- `src/cells.rs`: one cell: its entries and parts, the percent-decoding, the halves of a
  `word(score)` part, the parts no rule reads.
- `src/rows.rs`: the header lines, the blocks, the worker threads, the ordered writing, the
  summary line.
- `src/pyfmt.rs`: the CPython text rules named above; `src/unicode.rs`: the two Unicode tables.
- `src/*_tests.rs`: unit tests of the module each is named for.
- `tests/layouts.rs`, `tests/values.rs`, `tests/rows.rs`, `tests/header.rs`,
  `tests/command_line.rs` and `tests/field_map.rs`: the binary over fabricated tab output: one
  cell in the layout of each source column, lines as fastVEP 0.5.0 wrote them over the
  repository's synthetic fixtures, escapes and lists, empty parts, rows and carriage returns,
  headers, options, files that cannot be opened and malformed maps. `tests/field_map.rs` also holds
  the shipped map to the references and the input records of the repository's golden corpora under
  `tests/golden/116/`. `tests/common/` holds what the test files share.
- `tests/golden.rs` and `tests/fixtures/`: `cases.tsv` lists command lines over the fabricated
  `inputs/` and `maps/`, and `expected/` holds each case's stderr and its tuple file where it
  writes one. `cargo test` runs the binary over the cases, on one thread and on several, and
  compares bytes; a usage case is held to its exit code alone. A `.gitattributes` beside this page
  keeps the fixtures' line endings as committed.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-fastvep-plugin-tuples/Cargo.toml
cargo test --locked --manifest-path tools/vep-fastvep-plugin-tuples/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-fastvep-plugin-tuples/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-fastvep-plugin-tuples/Cargo.toml -- --check
```

The goldens under `tests/fixtures/expected/` are the contract in bytes, and `cargo test` holds the
binary to them. A change that is meant to move the contract rewrites the goldens it moves with the
binary's own output, so that the diff of the goldens shows the change. From this directory, which
makes a line that names a fixture hold its relative path as the goldens do:

```bash
target/release/vep-fastvep-plugin-tuples tests/fixtures/inputs/sources.tab \
    tests/fixtures/expected/sources.tuples.tsv 2> tests/fixtures/expected/sources.stderr.txt
```
