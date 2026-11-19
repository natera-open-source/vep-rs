# vep-formats-remainder

The front end of `scripts/concordance/compare_vep_formats.py`, the per-column format comparator of
the comparison oracle, for its `default` and `tab` formats. It reads a reference output and a
candidate output on every core, keeps back the record groups the two write alike, and writes the
rest, the remainder, for the Python comparator, with the counts of what it kept back. Every rule,
class, mask and report text lives in the Python comparator: this tool decides one thing, whether a
group is the same text on both sides, and hands a group it cannot decide over whole. It exists
because the Python parses and pairs every row of both outputs on one core, and on a population
dataset nearly every row is the same on both sides.

```text
vep-formats-remainder --expected <reference output> --actual <candidate output> --out-dir <dir>
    [--force-locations <file>] [--force-text <file>] [--keep-values <column> <file>] [--unordered]
```

## Contract

**Reads.** Two outputs as the comparator's `read_text` and `text_lines` read them: UTF-8, gunzipped
when the name ends in `.gz` (every member), a line ended by `\n`, `\r\n` or a lone `\r`, a last line
without its end kept; zero bytes after a gzip member are padding. A gzip member is read as Python's
`gzip.open` reads it: its header checksum is skipped without a check, a flag bit the format reserves
is read past, and the checksum and length of its text are held to its trailer. Its deflate stream is
inflated by `zlib-rs`, which refuses the streams zlib refuses (a copy from before the stream's
start, which flate2's default decoder takes). Text that is not UTF-8, a file that cannot be read, a
line of 4 GiB or more and an output of more than 4,294,967,295 lines are refused in one line,
`ERROR: [vep-formats-remainder] <detail>`, exit 1; a usage error exits 2.

**Groups.** A record group is every data line of one side that shares a Location, wherever in the
file its lines stand. A group is kept back when its lines are the same text in the same order on
both sides and nothing below hands it over. The comparator pairs entries by a key that starts with
the Location and numbers a repeated key in line order, so a kept-back group pairs each line with
itself.

**Handed over whole.** The comparator names, at run time, what marks a group it must see:

- `--force-locations`: one Location per line (an empty line is the empty Location). Every group of
  a listed Location is handed over; the comparator lists the Locations its manifest documents,
  whose keys it checks on equal values too.
- `--force-text`: one text per line. Every group holding a reference line that contains a listed
  text is handed over; the comparator lists the text its reference-side rewrite looks for.
- `--keep-values <column> <file>`: one value per line. Every group holding a reference line whose
  `<column>` is not a listed value is handed over; the comparator lists the Feature_type values of
  the rows it keeps when a class drops reference rows.

A group is also handed over when one of its lines, on either side, holds a byte outside ASCII, a run
of more than 640 digits or more than 64 opening brackets (`[` or `{`). The comparator converts some
values before it compares them, a whole number or an object's text, and its interpreter takes a
digit outside ASCII for a digit without converting it, refuses a number longer than it converts (a
length that can be set as low as 640 digits) and stops at an object nested past its call depth. A
value it does not convert ends the comparison with the comparator's error, and such a line reaches
the comparator on the whole outputs and on the remainder alike.

**Order.** The comparator reports the first place where the keys both sides share come in another
order. Without `--unordered`, a group is also handed over when the kept-back lines would not stand
in one order on both sides, and when it stands between two lines of a handed-over group that both
sides hold; the place the comparator finds among the remainder's shared keys is then its place
among all of them, moved by the kept-back lines before it. `--unordered` is for a comparison that
reads no order.

**Everything handed over.** Where a line's text does not settle its entry, both outputs are written
out whole and nothing is kept back: a header line after a data line, a data line before any column
line, a column named twice, column lines that differ between the outputs, a data line whose width
is not the column line's, a key column (Location, Allele, Feature, Feature_type) holding `|` or `#`,
a Location holding a space (the comparator labels a record by its Location, a space and its
Uploaded_variation), an `Extra` pair named as a key column, and a `--force-text` text in the
reference's column line or in a key column of a reference line, where a rewrite would rename a
column or move a row to another key. `counters.json` names the shape under `whole`.

**Writes**, under `--out-dir`:

- `expected.txt` and `actual.txt`: each side's `#` lines as they stand, then every handed-over data
  line in file order, each ended by `\n`.
- `counters.json`: `whole` (the shape, or `null`); `kept_entries`, the lines kept back on each
  side; `kept_populated`, per field name the kept-back lines that carry a value, counted as the
  comparator's `populated_counts` counts entries (every column, every `key=value` pair of `Extra`
  as a field of its own, `-` and the empty text no value); `kept_before`, the reference remainder's
  data lines in file order as runs of `[lines, kept-back lines before them]`; `groups` and
  `groups_handed_over`; and per side the data lines, the lines handed over and their bytes.

The three files are the same bytes at every worker count. `VEP_FORMATS_REMAINDER_THREADS` bounds
the workers; unset, the tool uses every core available to the process.

## Layout

- `src/input.rs`: reading, gunzipping, the UTF-8 check and the lines, cut on every worker.
- `src/layout.rs`: the column line, the key columns' places and the shapes that hand everything over.
- `src/scan.rs`: per line, the Location's hash, the comparator's lists and the per-line shapes.
- `src/groups.rs`: the groups, their comparison and the two order conditions.
- `src/counters.rs`: the populated counts of the kept-back lines.
- `src/output.rs`: the remainder files and `counters.json`.
- `tests/remainder.rs`: planted pairs of outputs, every row written by the test, run at one worker,
  two and the default and held to the same bytes and to the remainder each shape must leave.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-formats-remainder/Cargo.toml
cargo test --locked --manifest-path tools/vep-formats-remainder/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-formats-remainder/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-formats-remainder/Cargo.toml -- --check
```

The root `.cargo/config.toml` applies, so a release binary carries the same `target-cpu` floor as
`vep`. CI runs the test, clippy and fmt commands above in its Test, Clippy and Format jobs and
builds the package on the declared `rust-version` in MSRV. Comparator Tests builds the release
binary and runs `scripts/concordance/test_compare_vep_formats.py`, whose hand-over tests compare
planted outputs on the whole outputs and through this tool and hold the report, the standard output
and the exit code of the two to the same bytes; Harness Tests builds it for the one check of
`scripts/concordance/test_feature_oracle.sh` that runs a case through the runner on both sides of
the switch.

## Harness

`compare_vep_formats.py` hands a `default` or `tab` comparison over to this binary when
`VEP_FORMATS_REMAINDER` names an executable (a path, or a command on `PATH`), and compares the whole
outputs otherwise. A comparison naming a structural-variant class or `individual_allele_order`, one
under `--keep-only`, every other format and the other two modes compare the whole outputs with the
variable set as well, and the comparator says which it did in one line on the standard error. A
value that names no executable, a front end that exits non-zero and a list file the comparator
cannot write for it end the comparison with the comparator's error, exit 2, before the report file
is opened. `feature_oracle.sh` starts the comparator with its own environment, so
exporting the variable before the runner hands over every such comparison of a run. Nothing exports
the variable by default.
