# vep-discordant-classes

Classes the rows of a discordant-tuple table. A concordance comparison lists the tuples that one
engine wrote and the other did not, a row a tuple. This tool pairs the reference-only rows of such a
table with its scored-only rows and names, for every row, why it has no equal on the other side: the
end of its Location alone differs, its consequence set differs, the other engine annotated the
variant under other transcripts, or the other engine wrote nothing for the variant. It reads the
table of the tuple comparator (`discordant.tsv` of `tools/vep-outputs-compare`) and the table of the
structural-variant comparison, writes the table again with the class in a last column,
`mismatch_type`, and prints a summary on stdout. Both are fixed to the byte.

```text
vep-discordant-classes <discordant.tsv> [-o <classified table>]
```

- The one argument is the table to class.
- `-o`, `--output`: the table to write. The default is `discordant_classified.tsv` in the input's
  directory.

The command line is read as CPython's `argparse` reads it. `--output` may be shortened to any prefix
(`--out`), takes its value as the next argument or after `=`, and `-o` takes it as the rest of its
own argument as well (`-oclassified.tsv`); a later value replaces an earlier one. A value that starts
with `-` is taken only when it is `-` alone, starts like a negative number or holds a space. `--`
ends the options. `-h` and `--help` print the usage and the options.

## Contract

**Columns.** The first record names the columns. The table of the structural-variant comparison
names the engine of a row `direction`, its file `source_file` and its consequence set `consequence`;
where a table lacks one of those names, the tuple comparator's is read in its place: `source`,
`file_name`, `consequence_set`. Beside these three a table needs `location`, `allele` and `feature`,
and every other column is carried into the written table unread. A row is reference-only where its
engine column holds `only_perl` or `missing_in_rust` and scored-only on any other value.

**Classes.** A row's key values are its Location, the start of its Location, its Allele, its Feature
and its consequence set, each compared as text. The start is the text between the first colon of the
Location and the first colon or `-` behind it, the whole Location where it holds no colon; the
chromosome is no part of it. Rows are taken in file order, and each row is paired once at most.

1. `coordinate_only`: a reference-only row and a scored-only row with one start, Allele, Feature and
   consequence set and two different Locations. A reference-only row looks through the scored-only
   rows of its key from where the reference-only row of that key before it stopped, and takes the
   first whose Location differs from its own. The rows it passes on the way, each at its own
   Location, are not looked at again, and neither is any row of the key once a reference-only row
   has found no partner.
2. `consequence_swap`: among the rows the first pairing left, a reference-only row and a scored-only
   row with one Location, Allele and Feature. Each reference-only row takes the first scored-only
   row of its key that is not taken; the two consequence sets are not compared.
3. `transcript_selection`: a row left unpaired whose start and Allele the other engine wrote, in any
   of its rows, under other Features alone.
4. `true_missing` for every other reference-only row and `true_extra` for every other scored-only
   row.

**The written table.** One header row, the names read and `mismatch_type` behind them, then every row
of the table in file order with its class in that last column. Rows end in CR LF. A row is written
with as many fields as the header names: one that ends early is filled with empty fields and the
fields behind the header's last column are dropped. An empty line is no row. A column whose name the
header holds more than once shows, in each of its places, the value of the last such column. A
column named `_start` shows the start of the row's Location. A table that already holds a
`mismatch_type` column is not classed again: each row keeps the value it holds there, in that column
and in the appended one, and only the rows the second pairing takes are written `consequence_swap`.

**The summary.** Stdout carries a title between two rules of 70 `=`, `Total discordants: <rows>`, and
three blocks, each under its own line and a blank line: the rows of each class with their share of
all rows in percent, the rows of each class and engine value, and the thirty largest counts of file
and class with their share. Every block runs from its largest count down, equal counts in the order
in which the table first shows them. A class, an engine value or a file that the summary would have
to print for a row that ends before that column cannot be printed: the summary ends with the line
before it and the run with exit code 1, behind a whole written table.

**Bytes.** Wherever a byte could be read or written in more than one way it follows a CPython rule,
named here by its CPython name.

- The input is text read with universal newlines: strict UTF-8 under any locale, and `\r\n` and a
  lone `\r` each read as `\n`. A byte order mark is not skipped, so the first column name of a table
  that starts with one is not the name it shows.
- Records are those of `csv.reader` in the `excel` dialect with a tab delimiter. A field that opens
  with `"` runs to its closing quote across tabs and line endings, `""` inside it is one quote, text
  behind the closing quote is appended to the field, and a quote left open ends with the input. A
  `"` anywhere else is a character of its field. No blank is stripped. A field holds at most 131,072
  code points (`csv.field_size_limit()`).
- Rows are those of `csv.DictReader`: the first record is the header, an empty line included, and
  empty lines behind it are skipped.
- The written table is `csv.writer`'s under `QUOTE_MINIMAL`: a field is written between quotes, its
  quotes doubled, only where it holds a tab, a quote or a line ending.
- The summary's lines are `f"  {class:25} {count:>8} ({share:5.1f}%)"`,
  `f"  {class:25} {engine:12} {count:>8}"` and `f"  {file:35} {class:25} {count:>8} ({share:5.1f}%)"`:
  a text is filled with blanks up to its width in code points and never cut, and a share is
  `100 * count / total` to one decimal of its exact binary value, a tie going to the even digit.
  The order of a block is `Counter.most_common()`'s.
- The default table is `str(pathlib.Path(input).parent / "discordant_classified.tsv")`.

**Order of work.** The input is read to its end, and every row classed, before the written table is
made or, where it exists, emptied. A run that refuses its input therefore leaves an existing table
as it was and prints nothing. A table that names the input itself, by its name or through a link,
is written over it with the classes of the rows it held; the default table of an input named
`discordant_classified.tsv` is that case.

**Exit codes.** 0 when the table and the summary are written. 2 on a command line that is refused:
the usage line, then `vep-discordant-classes: error: <message>`, on stderr. 1 on any other failure
but one, with the single line `ERROR: [vep-discordant-classes] <detail>` on stderr:

- `<input>: <the system's error>` for an input that cannot be opened or read, `<input>: is a
  directory`, and `<output>: <the system's error>` for a table that cannot be made or written;
- `<input>: the byte at offset <n> is not UTF-8`, the offset counted from 0; a gzipped table is this
  case;
- `<input>: no header row` for an input of no bytes, and `<input>: Unsupported discordant TSV schema;
  missing required columns: <names>`, the names in code point order;
- `<input>: row <n> ends before its location column`, the rows counted from 1 behind the header;
- `<input>: field larger than field limit (131072)`;
- `<input>: more than 4294967294 rows`, the bound of the tool;
- `<input>: the summary ends at a row that ends before its <column> column`;
- `<input>: the table changed while it was read`, where the second read does not find the rows of
  the first;
- a thread bound that is no count.

A summary that stdout does not take ends the run with exit code 120, the code CPython ends with
when it cannot flush stdout at exit, and the line `ERROR: [vep-discordant-classes] stdout: <the
system's error>`.

## Threads and memory

The class of a row depends on rows anywhere in the table, so the table is read twice and its text is
never held. The first read keeps four numbers a row, one for each of its Location, Allele, Feature
and consequence set, and each distinct value of those columns once; the pairing works on the
numbers; the second read writes each row with its class and counts the summary. Memory follows the
rows and the distinct values of the four columns, not the bytes of the table. Of the two halves of a
table, its reference-only rows and its scored-only rows, the pairings index the smaller: the rows of
the larger half are held only under the keys the smaller half holds.

Two inputs are held whole instead, since they can be read only once: an input that is no regular
file (a pipe, `/dev/stdin`) and an input that the written table replaces.

The work is spread over a few threads, each taking its share in row order, so every bound writes the
same bytes. In the first read the calling thread cuts the text into records and up to four threads
number the four columns. The two pairings run on one thread, in file order as their definition has
it, beside a thread that collects what each engine wrote; the classes are then filled in across the
bound. In the second read the calling thread lays out the rows, one thread writes them and one
counts. `VEP_DISCORDANT_CLASSES_THREADS=<n>` bounds the threads at work, `1` keeping the whole run on
the calling thread; unset or empty, the bound is the host's available cores. It is an environment
variable and not an option, since the command line names the table and the bound changes no byte
written, and a value that is not a whole number of 1 or more is refused.

## Layout

- `src/main.rs`: the run: the input as a file read twice or held whole, the two reads, the order in
  which the files are opened, the error lines and the exit codes.
- `src/cli.rs`: the command line as `argparse` reads it, the usage line and the help text;
  `src/cli_parse_tests.rs` holds its cases.
- `src/text.rs`: the input as text, in runs: strict UTF-8 and universal newlines.
- `src/csv.rs`: the reader of records and the writer of a field; `src/csv_tests.rs` holds both to
  what the `csv` module makes of the same text, at every run length.
- `src/table.rs`: the header, the two schemas, what each written column shows, and the start of a
  Location.
- `src/intern.rs`: the distinct values of a column as numbers, the batches of values on their way to
  a thread, and sets of rows.
- `src/scan.rs`: the first read: the header and the numbers of every row.
- `src/classes.rs`: the two pairings and the class of every row; `src/classes_tests.rs` holds them to
  tables worked by hand, at every bound.
- `src/report.rs`: the second read: the written table, the counts and the summary's text.
- `src/pyfmt.rs`: the CPython formatting rules named above, each with a unit test that holds its
  bytes.
- `src/threads.rs`: the thread bound.
- `tests/fixtures/`: made-up tables (`<name>.tsv`), and for each the goldens that hold the written
  table and the summary to the byte (`expected_<name>.tsv`, `expected_<name>.stdout`; a run that
  writes neither has neither) and its exit code in `cases.tsv`. A `.gitattributes` beside this page
  keeps their line endings as committed.
- `tests/golden.rs`: the binary over every fixture under five thread bounds, and the error line of
  each refused fixture.
- `tests/command_line.rs`: the binary's exit code and bytes for the default table, a table written
  over its input, a pipe, a device, files that cannot be read or written, refused command lines, the
  thread bound and help. `tests/common/` holds what the test files share.
- `tests/input_shapes.rs`: the binary over tables made at run time and held to bytes laid out row by
  row: a table of many runs under every bound, every kind of line ending, quoted fields across the
  runs, and fields at and past the limit.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-discordant-classes/Cargo.toml
cargo test --locked --manifest-path tools/vep-discordant-classes/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-discordant-classes/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-discordant-classes/Cargo.toml -- --check
```

The goldens are the contract in bytes, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the goldens it moves with the binary's own output, so that the
diff of the goldens shows the change:

```bash
cd tools/vep-discordant-classes
target/release/vep-discordant-classes tests/fixtures/tuple.tsv -o tests/fixtures/expected_tuple.tsv \
    > tests/fixtures/expected_tuple.stdout
```

## Harness

`scripts/concordance/analyze_all_discordants.py` runs this tree's release build of the tool,
`tools/vep-discordant-classes/target/release/vep-discordant-classes`, on every `discordant.tsv` it
finds, and reads the classified tables into its pattern map. It builds nothing: a build that is
missing or is not an executable stops the script before it classes a table, with one line that names
the build command above. A table the tool refuses is skipped with a warning that carries the tool's
error line.
