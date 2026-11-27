# vep-fields-compare

The per-column field comparator of the concordance harness. It reads two VEP default-format outputs,
a reference and a scored one, pairs their rows on the consequence key (Location, Allele, Feature,
Feature_type and the consequence set) and reports, for every other column and for every `KEY=VALUE`
pair of the Extra column, on how many of the paired rows the two outputs agree. It writes one JSON
report, `fields.json`, and a summary on stdout; both are fixed to the byte.

```text
vep-fields-compare --perl <reference output> --rust <scored output> --out <fields.json> \
    [--examples 20] [--tmp-dir DIR] [--keep-only regulatory_feature_rows]
```

- `--perl`, `--rust`: the reference output and the scored output. A name that ends in `.gz` is
  gunzipped on read.
- `--out`: the report. Its parent directories are created.
- `--examples`: how many mismatch shapes a field keeps, the most frequent first (default 20). A count
  at or below zero keeps none.
- `--tmp-dir`: where the directory of sort intermediates is made (default: the system's temporary
  directory).
- `--keep-only regulatory_feature_rows`: drops every row of both inputs whose Feature_type is neither
  `RegulatoryFeature` nor `MotifFeature` before sorting, so the key totals, the duplicate counts and
  the field comparison cover those rows alone; the report then ends with the key
  `"keep_only": "regulatory_feature_rows"`. It is the one class the option names.

The command line is read as CPython's `argparse` reads it. An option may be shortened to any prefix
that fits it alone (`--ex 5`), takes its value as the next argument or after `=`, and a later value
replaces an earlier one. A value that starts with `-` is taken only when it is `-` alone, starts like
a negative number or holds a space. `--examples` takes what `int()` takes in ASCII digits: blanks
around the number, a sign, single underscores between digits, at any length; a digit outside ASCII,
which `int()` reads, is refused. `-h` and `--help` print the usage and the options.

## Contract

**Rows and keys.** A line that starts with `#` is skipped, and so is a row under seven columns. A row
of seven to thirteen columns reads its missing value columns as empty, and columns past the
fourteenth are not read. A row's key is its Location, Allele, Feature and Feature_type and its
consequence set: the Consequence column split on `,`, emptied of empty terms, deduplicated and sorted
by code point. A key that one side alone holds is counted under `reference_only` or `scored_only` and
compared no further. A key that a side repeats counts once, keeps the first of its rows in the
bytewise order of key and values, and adds the others to `reference_duplicate_keys` or
`scored_duplicate_keys`.

**Values.** On a key both sides hold, Uploaded_variation, Gene, cDNA_position, CDS_position,
Protein_position, Amino_acids, Codons and Existing_variation are compared verbatim: `-` and the empty
string are different values and nothing is normalised. The Extra column is split on `;` into
`KEY=VALUE` parts: a part without `=` is a key with the empty value, a repeated key keeps its last
value, and `-` and the empty column hold no key. Every key either side holds is compared under
`extra`, and a key one side lacks compares as `<absent>` there.

**The report.** `fields.json` holds, in this order, `totals` (the seven key counts), `fields` (the
eight columns), `extra` (the Extra keys in code point order), `inputs` (the two input names),
`key`, `comparison` and, under `--keep-only`, `keep_only`. A field's block holds `compared`,
`equal`, `differ`, `agreement` (`equal / compared` rounded to six places, `null` where nothing was
compared), `distinct_mismatch_shapes` and `top_mismatches`, each a `reference` value, a `scored`
value and a `count`. Stdout holds one line of key totals and one line a field, an Extra key as
`Extra:<key>`.

**Bytes.** Wherever a byte could be written in more than one way it follows a CPython rule, named
here by its CPython name. Records end at `\n`, `\r\n` or a lone `\r`, as universal newlines read
them, and are decoded as UTF-8 with `errors="replace"`. Fields are stripped with `str.strip()`'s
whitespace set (U+001C to U+001F included). Every path option is read through
`str(pathlib.Path(arg))`, and the `inputs` block holds those strings; a byte of an input's name
outside UTF-8 is written as `surrogateescape` reads it (`\udcXX`, in the `inputs` block and in the
not-a-file line). The JSON is `json.dumps(indent=2, ensure_ascii=True)` and a final newline. An
agreement is `round(x, 6)` written as `float.__repr__`, `None` on stdout where the report holds
`null`. The stdout counts are grouped as `f"{n:,}"`. `top_mismatches` keeps `Counter.most_common`'s
order: by count, equal counts in the order the merge first met them.

**Inputs.** An input whose name ends in `.gz` is read as the `gzip` module reads it: member after
member, the zero bytes that pad the file after a member skipped, a zero-byte file an empty input, and
any other byte where a member would start a corrupt stream. Both sides are sorted bytewise by a child
`sort -S 20%` under `LC_ALL=C` with `LC_COLLATE` removed; the two sorts run in parallel and inherit
`TMPDIR`, the rows reach them as UTF-8 under any locale, and a host without `sort` sorts each side
in memory to the same order. The sorted sides live in a directory `vep_fields_<8 characters>` under
`--tmp-dir`, removed when the run ends, and the merge reads them as two streams.

**Exit codes.** 0 when the report and the summary are written. 2 on a command line that is refused
(the usage line, then `vep-fields-compare: error: <message>`, on stderr) and when an input is not a
file (`ERROR: [compare_vep_fields] not a file: <p>` on stderr, the name in brackets a fixed part of
the line). 1 on any other failure: a `sort` that fails, an input that cannot be read or is a corrupt
gzip stream, a `--tmp-dir` that is no directory, an `--out` that cannot be written. Such a run writes
one line, `ERROR: [vep-fields-compare] <detail>`, on stderr and no summary.

## Layout

- `src/main.rs`: the run: the two sorts in parallel, the merge, the report and the summary.
- `src/cli.rs`: the command line as `argparse` reads it, the usage line and the help text.
- `src/gzip.rs`: the members of a gzip file and the zero padding after them.
- `src/input.rs`: records, lossy decoding, the `--keep-only` row class, the prepared row (five key
  columns, nine value columns) and the `sort` child.
- `src/merge.rs`: the sorted merge, duplicate-key collapse, per-field and per-Extra-key counts and
  mismatch shapes.
- `src/report.rs`: the JSON emitter and the stdout lines.
- `src/pyfmt.rs`: the CPython rules named above, each with a unit test that holds its bytes.
- `src/main_compare_tests.rs`, `src/input_consequence_set_tests.rs` and `src/cli_parse_tests.rs`: unit
  tests of the comparison over the `pair_*` fixture, of the consequence set in the key and of the
  command line.
- `tests/fixtures/`: made-up output pairs (`reference.txt` and `scored.txt`, then
  `<name>_reference.txt` and `<name>_scored.txt`) and the goldens `expected_<name>_fields.json` and
  `expected_<name>_stdout.txt` that hold the report and the summary of each to the byte. A
  `.gitattributes` beside this page keeps their line endings as committed.
- `tests/golden.rs`: the binary over `reference.txt` and `scored.txt` under two `--examples`
  counts, a pair with no shared key, the `pair_*` pair, gzipped inputs, input names outside UTF-8
  and the `regulatory_*` pair with `--keep-only` and without.
- `tests/input_shapes.rs`: the binary over pairs that each hold one shape of input (a zero-byte
  file, repeated keys, Extra keys on one side, ragged rows, line endings, blanks and bytes outside
  UTF-8, a byte order mark, tied mismatch shapes) and over a generated pair of 20,001 keys.
- `tests/command_line.rs`: the binary's exit code and bytes for abbreviated options, counts, refused
  command lines, a missing, unreadable or directory input, an unwritable `--out`, a host without
  `sort`, a `sort` that fails, and help. `tests/common/` holds what the test files share.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-fields-compare/Cargo.toml
cargo test --locked --manifest-path tools/vep-fields-compare/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-fields-compare/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-fields-compare/Cargo.toml -- --check
```

CI runs the test, clippy and fmt
commands above in its Test, Clippy and Format jobs, builds the package on the declared `rust-version`
in MSRV and audits its lockfile.

The goldens are the contract in bytes, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the goldens it moves with the binary's own output, so that the
diff of the goldens shows the change. From this directory, which makes the `inputs` block hold the
fixture paths as the goldens do:

```bash
target/release/vep-fields-compare --perl tests/fixtures/reference.txt --rust tests/fixtures/scored.txt \
    --out tests/fixtures/expected_fields.json > tests/fixtures/expected_stdout.txt
```

## Harness

`scripts/concordance/run_clone_measurement.sh` runs this tree's release build of the tool,
`tools/vep-fields-compare/target/release/vep-fields-compare`, on each SNP/indel cell of the vep-rs
engine, after the tuple comparison: the cell's two outputs as `--perl` and `--rust`,
`report/fields.json` as `--out` and the cell's directory as `--tmp-dir`. It reads no environment
variable for the comparator and builds nothing: a build that is missing or is not an executable stops
the script before its first engine run, with one line that names the build command above. The field
comparison is not asked for under `--walltime-only` or `--skip-fields`, nor for another engine. It
bears no F1: a comparison that fails is recorded as a `fields_compare_failed` row of the run's CSV
and the cell completes.
