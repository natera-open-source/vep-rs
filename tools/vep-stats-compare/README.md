# vep-stats-compare

Section-by-section comparison of two VEP statistics reports (the `<output>_summary.txt` that
`--stats_text` writes; a golden corpus keeps Ensembl VEP's as `expected/<format>.summary.log`): one
Ensembl VEP report against a candidate's, a JSON report with a `verdict`, and exit 0 on `green`, 1 on
`red`, 2 on an error.

```text
vep-stats-compare --expected <Ensembl VEP report> --actual <candidate report> [--tolerance N] [--mask-row LABEL ...] [--report <stats.json>]
```

## Contract

**The rules.** The grammar is `Stats.pm:806-831` at release/116.3: a `[Section]` line opens each
section and a row is `key<TAB>value`, a row without a tab being a key with no value. The first three
sections are tabular (`[VEP run statistics]`, `[Data version]`, `[General statistics]`); every other
section is a chart of `key<TAB>count` rows.

- A tabular section compares row by row: the same keys in the same order, every value verbatim. The
  rows no second run reproduces compare by key alone: `Command line options`, `Start time`, `End time`,
  `Run time`, `Input file`, `Input data`, `Output file`, and every `--mask-row LABEL`, wherever the
  label occurs. In `Annotation sources` every path (a token carrying `/`, up to a blank or `;`) reads
  `[PATH]` before the values compare.
- One row of `[General statistics]` is held to a documented difference,
  [`stats_existing_variants_check_frequency`](../../docs/divergences/stats-existing-variants-check-frequency.md):
  `Novel / existing variants`, where the expected report ran a frequency check and counts every
  variant as existing. The expected report's `Command line options` row names `--check_frequency` or
  `--filter_common` as a word of its own (`Config.pm:734-749`; a flag before the value `0` is off),
  and its value is `0 (0.0) / N (100.0)`, N its `Variants processed`. The actual value then agrees
  when it is the row `Stats.pm:727-734` writes for some count of existing variants among the N: `-`
  for none, else `a (p) / b (q)` with `a + b = N` and each percentage `100 * (count / N)` under
  `%.1f`. Any other pair of values compares verbatim, and a `--mask-row` of the label sets both
  values aside unread, as for any row.
- One kind of key of `[Consequences (all)]` is held to a second documented difference,
  [`stats_consequences_custom_key_named_consequence`](../../docs/divergences/stats-consequences-custom-key-named-consequence.md):
  a key the expected report alone holds that names no consequence type, where the expected
  report's `Command line options` row names a `--custom` source whose key is `Consequence`: the
  word `--custom` before a value whose last `short_name=` item is `Consequence`, or whose second
  item is where the value holds no `=` (`AnnotationSourceAdaptor.pm:351-363`). Such a key is the
  name of a record of that source, which Ensembl VEP counts once for each row of the allele the
  record names. A consequence type is a name of `data/consequence_terms.tsv` in any of its three
  vocabularies. A key that names one, a key the actual report holds too and a count that is no
  integer compare as in any chart, and so does a key of any other section.
- A chart section compares as a set of `(key, count)` pairs, since Ensembl VEP sorts ties in hash
  order. Under `--tolerance N` (N at least 1) the two counts of a key may differ by up to N, a key on
  one side alone counting as 0 on the other; at the default 0 the sets must be equal. A count is an
  integer where `int()` reads the value as one and the text as written where it does not; a text
  count agrees only with the same text.
- Both files must carry the same sections in the same order. A section on one side alone, the shared
  sections in another order, a section title or a row key printed twice (its first copy is the one
  compared) and a file without a section are each a reason for `red`.

**The outputs.** Stdout carries `stats: <verdict>`, then one line per reason and one per row set
aside under a documented difference (`[<section>]: <row> set aside as <class>`), each indented by two.
`--report` writes the JSON report, its directory made first: `mode` and `format` (both `stats`),
`verdict`, `reasons`, `tolerance`, `masked_rows` (sorted), `sections` (each file's titles in order and
`order_agrees`), `compared` (`sections` and `rows`), `documented` (the rows set aside under a
documented difference, each with its `class`, `section`, `row`, `expected` and `actual`, a chart
key with its count as written and no `actual`; an empty list where none was), `differences`
(section title to its list of
differences, each a `kind` of `value`, `count`, `missing`, `unexpected`, `order`, `duplicate`,
`section_missing` or `section_unexpected` with `row`, `side`, `expected` and `actual` where they
apply), then `expected` and `actual`, the two paths. A `red` verdict writes the report as a `green`
one does; an error writes none.

**The errors.** One line on stderr, `ERROR: [compare_vep_stats] <detail>`, and exit 2. The details, in
the order they are looked for:

1. `cannot read <path>: <reason>` for the expected file, then for the actual file. The reason is the
   text of the exception Python raises for the same file: `[Errno N] <strerror>: <repr of the path>`
   for a file that cannot be opened, `'utf-8' codec can't decode byte 0x.. in position N: <reason>`
   (or `bytes in position N-M`) for text that is not UTF-8, and for a `.gz` file
   `Not a gzipped file (<repr of its first bytes>)`, `Unknown compression method`,
   `Compressed file ended before the end-of-stream marker was reached`,
   `CRC check failed 0x.. != 0x..`, `Incorrect length of data produced` or
   `Error -3 while decompressing data: <zlib's name of the fault>`. A file that cannot be read to its
   end is reported before text that is not UTF-8.
2. `--tolerance N is negative`.
3. `line N precedes the first [section] header: <repr of the line>` for the expected file, then for
   the actual file: a non-blank line before the first section means the file is no statistics report.
4. `cannot write <path>: <strerror>` for the report.

A usage error (an unknown option, a missing `--expected` or `--actual`, a `--tolerance` that is no
integer) prints a usage line and exits 2; `-h` and `--help` print the options and exit 0. The command
line is read as Python's `argparse` reads it: an option by its name or any unambiguous prefix of it,
its value as the next argument or after `=`, a value that starts with `-` only where it opens like a
negative number or carries a space, the last value of a repeated option, every `--mask-row`.

**The bytes.** Every byte of stdout, of an error line and of the report follows the CPython rule that
defines it:

- A file whose name ends in `.gz` is read as Python's `gzip` module reads one: every member in turn,
  each held to its CRC and length, zero bytes between and after members skipped, an empty file an
  empty text.
- Text is strict UTF-8 read with universal newlines, split into lines as `str.splitlines()` splits it
  (`\n`, `\r`, `\r\n`, `\v`, `\f`, U+001C to U+001E, U+0085, U+2028, U+2029); a line is blank where
  `str.strip()` leaves nothing (U+001C to U+001F included). The blanks around a path in
  `Annotation sources` are those of `\s` in `re`, the same set.
- A count and `--tolerance` are read as `int()` reads a `str`: blanks at both ends (U+001C to U+001F
  excepted), one sign, decimal digits of any script with single underscores between them, an integer
  of any size up to 4300 digits.
- The report is `json.dumps(indent=1)` under `ensure_ascii=True`, followed by one newline. A path is
  written as `str(pathlib.Path(argument))`, and its `.gz` ending is read from that form. `masked_rows`
  is in the order `sorted()` gives a list of `str`. An argument that is not UTF-8 is written as the
  `str` Python decodes it to under `surrogateescape`: `\udcxx` for each undecodable byte, in the report
  and on stderr.
- A `repr` in an error line is `repr()` of a `str` or of `bytes`: the quote chosen by the quotes in
  the text, `\t`, `\n` and `\r` by name, and each code point `str.isprintable()` refuses as `\xNN`,
  `\uNNNN` or `\UNNNNNNNN`.

The printable code points and the decimal digits are those of Unicode 16.0.0, the database of
CPython 3.14.

**Memory and cores.** Each file is read as a stream on a thread of its own and is never held whole.
The parsed rows of both reports are held until the comparison ends, because a section of one file is
looked up by title in the other and a row by key.

## Layout

- `src/main.rs`: the run, the order of the errors, the report file, the exit code.
- `src/cli.rs`: the command line.
- `src/input.rs`: one file as a stream: gzip members, strict UTF-8, lines, sections.
- `src/report.rs`: sections, the comparison rules, the report and the stdout lines;
  `src/report_tests.rs` holds them to an Ensembl VEP 116.3 report with one change at a time, the
  rule of the documented row to two made-up reports of a run under `--check_frequency`, and the
  rule of the documented key to both engines' reports of a run under a source keyed `Consequence`.
- `src/terms.rs` and `data/consequence_terms.tsv`: the names a consequence type is counted under,
  one row for each entry of Ensembl VEP's consequence table with its three vocabularies.
- `src/json.rs`: the report's values and their text.
- `src/pyfmt.rs`: the CPython text rules listed above, each with a unit test pinning Python's output.
- `src/unicode.rs`: the two Unicode tables.
- `tests/golden.rs` and `tests/fixtures/`: `cases.tsv` lists the command lines over the fixture
  reports (`reference.summary.log` and `predictions.summary.log` are the Ensembl VEP 116.3 reports of
  the `GRCh38-statistics` corpus; the others are those reports with changes, small reports of counts
  and of text outside ASCII, the made-up `frequency_*` reports of a run under a frequency check, the
  `consequence_key_expected` and `consequence_key_actual` reports, which are both engines' reports
  of one record of the `GRCh38-custom` corpus under a source keyed `Consequence` with their paths,
  times and command lines rewritten, the two other `consequence_key_*` reports, which are those
  with one change each, and gzip files, whole and damaged), and `expected/` holds each case's
  stdout, stderr and report where it writes one. `cargo test` runs the binary over the cases and
  compares bytes; a usage case is held to its exit code alone.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-stats-compare/Cargo.toml
cargo test --locked --manifest-path tools/vep-stats-compare/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-stats-compare/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-stats-compare/Cargo.toml -- --check
```

CI runs the test, clippy and fmt
commands above in its Test, Clippy and Format jobs, builds the package on the declared `rust-version`
in MSRV and audits its lockfile in Advisory Reachability.

The goldens under `tests/fixtures/expected/` are the contract. A change that moves a byte of them
changes the contract, and the case's golden changes with it in the same commit.

The inflater is `flate2` on its `zlib-rs` backend, which names a damaged stream as zlib does. It
inflates into 256 bytes at a time: with 260 or more bytes of room that backend takes its fast loop,
which reports every fault it meets as `repeated call with bad state`, and three cases of `cases.tsv`
(`distance_too_far_back`, `invalid_distance_code`, `invalid_length_code`) hold the fault's own name.
