# vep-parquet-to-tab

Projects a vep-rs Parquet output back to the `--tab` format: the inverse of the writer's cell
contract (the module documentation of `crates/vep-io/src/output/parquet`), cell by cell. It reads one
part, a directory of parts or a `--parquet_partition chrom` directory across the host's cores and
writes the rows in file order. `--canonical` prints the canonical form of a Parquet or a `--tab`
output, in which two outputs that hold the same rows are equal byte for byte.

```text
vep-parquet-to-tab [--canonical] [--output OUT] PATH
vep-parquet-to-tab [--canonical] [--output OUT] --input PATH
vep-parquet-to-tab --canonical --types-of PARQUET [--output OUT] TAB
```

- `PATH`, or `--input PATH` (`-i`): a Parquet file or directory; with `--canonical`, a `--tab`
  output too, plain or gzip. `--input` wins over the positional when both are given and it is not
  empty.
- `--output OUT` (`-o`): the file to write; `-`, the default, is standard output.
- `--canonical`: the canonical form in place of the rows in file order.
- `--types-of PARQUET`: with `--canonical` on a `--tab` output, the Parquet output it is compared
  with; a column that output stores as strings stays as written. No other run reads it.

The command line reads as Python's `argparse` reads one: a long option under any prefix that names
it alone (`--canon`, `--out`), `--name=value`, `-ivalue` and `-i=value`, `--` before a path that
starts with a dash, the last value of a repeated option. An argument that looks like an option and
names none, a second path, a missing value or no path at all is a usage error.

## Contract

**Input.** A Parquet output is a regular file that opens with the Parquet magic bytes, or a
directory. A directory stands for every file under it whose name ends in `.parquet`, at any depth,
hidden files and directories included; a symbolic link to a directory is not followed. The files
are read in path order, the paths compared component by component, so `a/z.parquet` comes before
`a-c/x.parquet`. Every file leads with key columns (any of `contig_index`, `chrom`, `pos`, `end`,
`ref`, `alt`) and then holds the tab columns; the files of one output name the same tab columns.

**Rows.** The output is the header line, `#` and the tab column names joined by tabs, then one line
per row: the files in path order, each file's row groups in order. The `##` meta lines of a run are
not stored in Parquet and are not written.

**Column names.** The writer stores a tab column named as a key column as `tab_<name>`, so that no
two columns of a file share a name: the `ref` and `alt` columns the committed synthetic dbNSFP
fixture declares under `ALL` are stored as `tab_ref` and `tab_alt`, and a tab column whose name is
a key column's behind one or more `tab_` is stored behind one more (`tab_ref` as `tab_tab_ref`).
It stores a tab column whose name an earlier column carries in other letters behind `tab_` as well
(`tab_REF` for a custom source named `REF`), and states the `--tab` name of every tab column in
the footer, under the key `tab_columns`, joined by tabs. The header line names each column as the
footer states it, and the spelling of a BOOLEAN column and the canonical form follow that name.
For a file that states no such list, or a list of another length, a stored name that is a key
column's behind one or more `tab_` prints without the first, and every other name prints as
stored.

**Cells.** The types are read from the Parquet schema alone, so a string, a list or a number is one
type whatever Arrow type its writer declared for it.

- A NULL cell prints `-`. A string prints as stored: a stored `-`, an empty string, a tab or a line
  end inside it included.
- An integer of any width prints in decimal.
- A floating-point number prints as the shortest decimal that reads back to the same value, in
  plain notation (`0.095`, `0.000098`, `1`, `-0`, `1000000000000000000000`), `NaN`, `inf` or `-inf`;
  of two shortest decimals equally close to the value, the one of larger magnitude
  (`126430982921026.63`). A 32-bit or 16-bit value prints as the double it widens to.
- A boolean prints its column's spelling: `YES` for `CANONICAL`, `1` for `PICK` and
  `GENCODE_PRIMARY` (false prints `-` for these three), `Y` and `N` for `HIGH_INF_POS`, `1` and `0`
  for any other column and inside a list.
- A list prints its elements joined by `,`, a NULL element as `-`, an empty list as an empty cell.
- A struct is a positions cell: its integer members `start`, `end` and `total` print `start` when
  `start` and `end` agree, `start-end` otherwise with `?` for a NULL member, `-` when both are
  NULL, then `/total` when `total` is set. A member the struct lacks is NULL, and any other member
  is not read.
- A column stored under any other type (binary, a date or time, a decimal, a list of lists or of
  structs, a map, a struct whose `start`, `end` or `total` is no integer) is refused, by file and
  column, before any output.

**Columns by place.** A tab column is read from its own field: the field at its place after the
leading key columns, whatever its name. A tab column named as a key column the file leads with
(`ref` or `alt` among a plugin's own columns) prints the cells of its own field, and the row holds
a cell for every column of the header. That holds for a part that stores the column as
`tab_<name>`, as the writer does, and for one that stores it under the key column's own name.

**The canonical form** (`--canonical`) is the header line, then the data lines sorted in the order
of their bytes, each line with its DOUBLE and INT32 cells reprinted:

- A row of a Parquet input is its `--tab` line as above. A `--tab` input is read as strict UTF-8; a
  line ends at `\n` alone and loses every `\r` before it, so a cell that holds another separator
  character (U+2028, a form feed) stays in its row; a `##` line and an empty line are dropped; the
  one line that starts with a single `#` is the header, wherever it stands; every other line is a
  data line. A file that opens with the gzip magic bytes is read as gzip, any number of members,
  the zero bytes between and after them skipped.
- A line is split at its tabs and each cell is reprinted by the header's column at its position; a
  cell past the last column stays as written. A `-` cell stays `-`.
- DOUBLE (`MAX_AF`, `MOTIF_SCORE_CHANGE`) and `LIST<DOUBLE>` (`AF` and the 27 `<population>_AF`
  columns; each token between commas): a token inside the number grammar prints as the shortest
  decimal of its value, as a floating-point cell does; any other token stays as written. The
  grammar is the one a Python `re` pattern and `float()` read: an optional sign, then digits with
  an optional fraction or a fraction alone and an optional exponent, or `inf`, `infinity` or `nan`
  in any case. A digit is any Unicode decimal digit (general category Nd of Unicode 16.0), as `\d`
  matches and `float()` reads; one line end may close the token, as `$` allows. `inf`, `infinity`
  or `nan` spelt with U+0130 or U+0131 in place of `i` is inside the case-insensitive pattern and
  is no number: the run ends on it with exit code 2 after the header line is written.
- INT32 (`STRAND`, `DISTANCE`, `ALLELE_NUM`, `TSL`, `HGVS_OFFSET`, `OverlapBP`, `MOTIF_POS`): an
  optional sign and ASCII digits print as their integer (`007` is `7`, `+7` is `7`, `-0` is `0`); a
  cell that is no INT32, an empty one included, stays as written.
- A string is no number, whatever its column's name: a field that a plugin or a `--custom` source
  declares under one of these names is a string column of a Parquet output. A cell of a Parquet
  input that is printed from a string column stays as stored. A column of lists of strings is a
  string column too.
- A `--tab` input does not state the types of its columns, so each is reprinted by its name. With
  `--types-of PARQUET` a column stays as written where every file of that Parquet output stores
  every column of its name as strings.
- A `--tab` cell that an INT32 column of a Parquet output could not hold is NULL there and prints
  `-`, so the two forms differ on that cell.

**Memory.** Rows are written a wave of row groups at a time, one row group a thread and sixteen at
most, each read from its file in one piece: a wave is read across the threads while the wave before
it is written, so at most two waves are held as text. The canonical form is a sort of every line, so
the canonical lines of the one input are held until they are sorted and written, and a `--tab`
input is held whole while its lines are reprinted.

**Exit codes and the error line.** 0 on success, and when the reader of standard output closes it.
2 on a usage error (a `usage:` line, then the error) and on any refusal, reported as one line,
`ERROR: [parquet_to_vep_tab] <detail>`, on standard error. 1 on a gzip input whose deflate data is
corrupt. The details, `<path>` being the argument as `str(pathlib.Path(arg))` prints it:

```text
<path>: no such file or directory
<path>: neither a regular file nor a directory
<path>: not a Parquet file
<path>: no .parquet file under the directory
<file>: its tab columns differ from <first file>'s
no tab columns in the Parquet schema
<path>: no # header line; a tab output needs its header
<path>: two header lines
<file>: column <name> is stored as <type>, which the cell contract has no tab text for
<file>: row group <n> cannot be read
could not convert string to float: '<token>'
```

`not a Parquet file` is a regular file without the magic bytes, read without `--canonical` or
named by `--types-of`; the last line is the token of a DOUBLE cell described above. The other details follow the reader that
met the failure:

- A file that cannot be opened or written, as `str(OSError)` prints it, the name as `repr()` quotes
  it: `[Errno 13] Permission denied: '<path>'`.
- A `--tab` input that is not UTF-8, as the codec's `UnicodeDecodeError` reads, the position
  counted in the decoded data: `'utf-8' codec can't decode byte 0xff in position 0: invalid start byte`.
- A broken gzip stream, in the words of Python's `gzip` module:

  ```text
  Compressed file ended before the end-of-stream marker was reached
  CRC check failed 0x<stored> != 0x<computed>
  Incorrect length of data produced
  Not a gzipped file (b'<the two bytes read>')
  Unknown compression method
  ```

- A Parquet file without a whole footer, and a part that cannot be opened, in the words of the
  Arrow C++ reader; any other failure of the Parquet reader is `<file>: <the reader's message>`:

  ```text
  Parquet file size is 0 bytes
  Parquet file size is <n> bytes, smaller than the minimum file footer (8 bytes)
  Parquet magic bytes not found in footer. Either the file is corrupted or this is not a parquet file.
  Parquet file size is <n> bytes, smaller than the size reported by footer's (<m>bytes)
  [Errno <n>] Failed to open local file '<file>'. Detail: [errno <n>] <reason>
  ```

The file `--output` names is made, or emptied, before the input is read, so a refused run leaves it
empty. A failure after the header line leaves the header and, in file order, the row groups read
before it.

## Threads

The row groups of a Parquet input, the blocks of a `--tab` input and the sort run across the host's
available cores. `VEP_PARQUET_TO_TAB_THREADS=<n>` bounds the threads at work; it is an environment
variable and not an option, and a value that is not a whole number of 1 or more is refused. Every
thread count writes the same bytes.

## Layout

- `src/main.rs`: the two forms of output, the waves of row groups and their writer, the sort, the
  exit codes and the error line.
- `src/args.rs`: the command line.
- `src/parquet_input.rs`: the files of an output in path order, their footers, the tab columns by
  their place and the row groups.
- `src/cells.rs`: the `--tab` text of a stored cell, by type.
- `src/canonical.rs`: the columns of each canonical type, the number grammar, the INT32 grammar and
  the canonical line.
- `src/tab_input.rs`: a `--tab` input as blocks of whole lines, its header and its decode errors.
- `src/gunzip.rs`: gzip members and their failures.
- `src/pyfmt.rs`: the CPython formatting rules listed above, each with a unit test.
- `tests/fixtures/`: Parquet files under the writer's schema and under the other types a cell can
  be stored as, `--tab` files, and files that are refused; `cases.tsv` lists the runs over them, and
  `expected/` holds the standard output and the error line of each.
- `tests/golden.rs`: runs the binary over every case on one thread, on three and on the default,
  and compares bytes.
- `tests/cell_contract.rs`: the cell of every type family, the header, every input shape, the
  canonical form, `--output` and the refusals, each against values the test states.
- `tests/writer_round_trip.rs`: annotates three golden corpora, and the synthetic dbNSFP fixture
  corpus under the flag set whose plugin names columns `ref` and `alt`, with the `vep` binary into
  `--tab` and into Parquet and holds the two canonical forms equal.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-parquet-to-tab/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

The round trips of `tests/writer_round_trip.rs` run the `vep` binary: `VEP_BINARY`, else the release
or debug build under the repository's cargo target directory. Without one they pass without
running and print `skipped`; `VEP_ROUND_TRIPS_MUST_RUN=1` fails them instead, so run them after
`cargo test` of the workspace when the Parquet writer or its reader changes.

## Comparing a Parquet output with a `--tab` output

```bash
vep=tools/vep-parquet-to-tab/target/release/vep-parquet-to-tab
cmp <($vep --canonical out.parquet) \
    <($vep --canonical --types-of out.parquet out.tab)
```

The rows of a Parquet output are in the file's order, not the input's, and a frequency is stored
as a number, so `0.0010` in the `--tab` output is `0.001` in both canonical forms. A column the
Parquet output stores as strings holds each cell as the `--tab` output writes it, and `--types-of`
keeps it so in the canonical form of the `--tab` output.
