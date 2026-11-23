# vep-divergence-reevaluate

Re-runs every intended divergence on Ensembl VEP and on vep-rs and says whether each is still
present. Every entry of the category tables of `docs/intended-divergences.md` has a reproduction,
`scripts/divergences/reproductions/<entry>.toml`: the input, the flags, the cache or data the example
reads, the field or row where the engines differ, both engines' recorded values and the
`implementation` table that names the code, the tests and the comparison rules that implement the
entry. The schema of a record, the fixed arguments each engine is given and what each compare kind
reads are
[`scripts/divergences/reproductions/README.md`](../../scripts/divergences/reproductions/README.md).

```text
vep-divergence-reevaluate --image <Ensembl VEP image tag or digest> --vep <vep-rs binary>
    --perl-cache-dir <dir> [--perl-cache-dir <dir> ...] [--plugins-dir <dir>]
    [--tool-dir <dir> ...] [--extra-dir NAME=PATH ...] [--entry <name> ...] [--jobs <n>]
    [--work-dir <dir>] [--keep] [--timeout <seconds>] [--expect-all-present] [--json <file>]
    [--reproductions <dir>]
vep-divergence-reevaluate --check [--reference 116.3] [--page <page.md>] [--reproductions <dir>]
vep-divergence-reevaluate --list [--reproductions <dir>]
```

- A run executes each reproduction on both engines, Ensembl VEP first, reads the compared value
  from each output and prints one line an entry as it ends, `<entry>\t<verdict>[\t<detail>]`:
  `still-present` when the engines differ and each prints its recorded value, `closed` when the
  engines agree, `changed` when they differ in another way, `error` when a record, a run or a read
  failed. A record that holds further runs, the tables of its `also` list, is run once for each:
  the entry reads the verdict every run reads, `error` where a run failed and `changed` where its
  runs read differently, and its detail names a further run as `also[<place>]`, counted from 1. A
  blank line and a table of the verdict counts follow, then, where an entry is `closed` or
  `changed`, a blank line and a table of the code, tests and comparison rules its record names, and
  `outputs under <dir>` when the outputs are kept. `--json` writes the entries in name order, each
  with its verdict, both value lists and its detail.
- Ensembl VEP runs in Docker from `--image`; `--ensembl-vep <command>` runs a local installation
  instead, its `filter_vep` and `variant_recoder` beside its `vep`. Its cache is the
  `<species>[_merged|_refseq]/<version>_<assembly>` directory found under the `--perl-cache-dir` roots
  in order, the first that holds every chromosome the input names. `--plugins-dir` holds its plugin
  modules, a `--tool-dir` a program a record's `run.tools` names, and `--extra-dir NAME=PATH` a
  directory a record names as `{NAME}`.
- `--jobs` entries run at once (default 1), and one engine run is ended after `--timeout` seconds
  (default 1800). The engines write under `--work-dir`, `<entry>/stage` for what an entry prepares
  and `<entry>/out/<engine>` for an engine's output, its `command.txt`, `stdout.txt` and
  `stderr.txt`, and a further run's under `<entry>/also-<place>`; without `--work-dir` they write
  under a temporary directory that is removed unless
  `--keep` is given.
- Each engine runs in its entry's output directory. A relative `--vep`, `--perl-cache-dir`,
  `--plugins-dir` or `--work-dir` is read from the directory the tool is called in. So are the first
  and the last word of an `--ensembl-vep` command, its program and its script, where the word names a
  directory; any other path in that command has to be absolute. A `--vep` value written without a
  directory, `vep` and not `./vep`, is a command found on `PATH`.
- `--check` runs nothing. It prints `<n> entries have a reproduction verified against <reference>:
  <a> observed, <b> unobserved, <c> malformed` when every entry of the page's tables has a record
  that passes the schema, every entry's record holds an `implementation` table whose rows each
  name a file of the repository and an item that file defines, the record's `verified_against` and
  the page's `Verified against` cell equal `--reference`, no record lacks a row, every record's
  input class agrees with its entry's tier and section, and every `Tier and count` cell whose
  figure reads `<n> of <total> cache transcripts` names the same total. Otherwise each violation is
  one line on stderr; cells that
  name more than one total are one line, each total with its number of rows and the first entry
  that names it, in page order.
- `--list` prints each record's name, `verified_against` and `page`, tab-separated, and comes before
  `--check` when both are given.

The command line is read as CPython's `argparse` reads it: an option may be shortened to any prefix
that fits it alone, takes its value as the next argument or after `=`, and a later value replaces an
earlier one unless the option is repeatable. `--jobs` and `--timeout` take what `int()` takes in
ASCII digits. `-h` and `--help` print the usage and the options.

## Contract

Every line the tool prints and every file it writes is fixed to the byte. Where a byte could be
written in more than one way it follows a CPython rule, named here by its CPython name.

**Records.** A record is TOML 1.0, read with `toml_edit`, and a `\r\n` line end reads as `\n`. A
table's keys stand in the order the document first names them, which is the order violations are
listed in. An integer is held to 64 bits, and a file with a larger one is refused. A file that is not TOML is refused in one line,
`<path>: <what the reader refused> (at line <l>, column <c>)`, and bytes that are not UTF-8 in the
words of the `utf-8` codec's `UnicodeDecodeError`. The records of a directory are its `*.toml`
entries in name order.

**Values in messages.** A value a message quotes is written as `repr()` writes it (`'text'`, `5`,
`True`, `['a', 'b']`, `{'key': 1}`, `datetime.date(2020, 1, 1)`), and a value compared or passed as
an argument as `str()` writes it. A list of compared values is written as `repr()` of a `list` of
`str`, and a count of rows as `format(n, ",")`.

**Patterns.** `compare.pattern`, `compare.where.pattern` and each `compare.mask` are regular
expressions of CPython's `re` module for a `str` pattern, compiled with no flag: its syntax, its
inline flags, its Unicode categories and case folding, and the match `re.search` finds; a mask is
applied as `re.sub`. A pattern `re.compile` refuses is a schema violation that ends in the message
of `re.PatternError`, its position included. Three things differ from `re`. `\N{name}` is refused as an
undefined character name, since the tool carries no table of names. A conditional on a group from
inside that group reads the group as not matched. A repeat count of 4,294,967,295 or more is a
violation of the record, in the words `the repetition number is too large`. One behaviour of `re`
is kept on purpose: a search holds a match's first character to the set the pattern opens on, read
by the pattern's own flags, so a set under `(?a:` at the start of a pattern is searched by the
Unicode categories.

**Outputs.** An output file, a warnings file and a statistics file are read as text: a file that
starts with the gzip magic is decompressed whatever its name, the bytes are decoded as UTF-8 with
`errors="backslashreplace"`, and lines end as universal newlines end them. An engine's standard
streams are decoded the same way and split as `str.splitlines()` splits them. A line of a JSON
output is read as `json.loads` reads it (`NaN` and `Infinity` are numbers, an integer has any
length, a repeated key keeps its first place and its last value); a line that is not JSON is
skipped, a record that is no object is the entry's error, and half of a surrogate pair reads as
U+FFFD. A JSON value that is no string is compared as `json.dumps(value, sort_keys=True,
separators=(",", ":"))`.

**Comparison.** Each value has trailing whitespace removed (`str.rstrip()`) unless `strip = false`,
every absolute path replaced by `<path>` unless `mask_paths = false` and every `mask` match by
`<masked>`; the list is sorted by code point unless `sorted = false`. An entry's two outputs are
read whole, since its values are sorted and compared as two lists; an output is a record's few
lines.

**Commands.** `--ensembl-vep` is split as `shlex.split` splits it. A `prepare` step runs as `bash
-euo pipefail -c <step>` in the stage, and both of its streams must be UTF-8. An engine runs in its
output directory with the tool's own standard input unless the record pipes the input.
`command.txt` holds the command as `shlex.join` writes it and `exit <status>`, a signal as its
negative number. A container is named `reevaluate-<entry>-<process id>-<eight hex digits>` and is
removed when its run is ended at the time allowed. Where a `--vep`, `--perl-cache-dir`,
`--plugins-dir`, `--work-dir` or `--ensembl-vep` path is read from the calling directory, it reaches
the engines as `Path.absolute()` writes it: joined to that directory, no symbolic link resolved,
`..` kept.

**Paths.** A path option is read through `str(pathlib.Path(arg))`, and a record's path through
`Path.resolve()`: absolute, every symbolic link resolved, a component that does not exist kept. The
repository a record's relative path, `{repo}` and `REPO` name, and the default page and record
directory, is the nearest directory above the executable that holds
`scripts/divergences/reproductions`, else the checkout the tool was built in.

**Errors.** An entry that cannot run is an `error` whose detail is the last line of the failure: a
schema violation list joined by `; `, a refusal of the tool's own, a file error as `str(OSError)`
(`[Errno 2] No such file or directory: '<path>'`), a prepare step's last line of stderr, or
`Command '[...]' timed out after <n> seconds`. A failure outside an entry is one line on stderr,
`ERROR: [reevaluate] <detail>`, the name in brackets a fixed part of the line.

**Items.** A row of a record's `implementation` table is `<path>::<item>`: a path relative to the
repository and the chain of names from the file's top level to the item, joined by `::`
(`function`, `Type::method` for a method of any `impl` block of the type, `Trait::method`,
`tests::function`, a constant or a static by its name). A Rust file is read lexically, strings,
characters and comments passed over, and inside an attribute plain string literals alone; a row
means the first item of its name in the file, and a function written inside another function's
body has no name of its own. A braced constant expression inside a header's parentheses or
generic brackets is no body. An item spans the lines
from its first doc comment or attribute to its closing brace or semicolon, counted by line feeds.
Of a Python or a shell file `--check` asks only that it defines a function or a class of the
name. `scripts/ci/admission.py` reads items by the same rules, and one case table holds both
readers to one answer.

**Results.** The `--json` file is `json.dumps(results, indent=1)` and a final newline, each entry
an object of `entry`, `verdict`, `ensembl_vep`, `vep_rs` (each a list, `null` for an entry that did
not run) and `detail`. An entry found `closed` or `changed` whose record holds an `implementation`
table carries one key more, `implementation`: an object of the table's lists that hold a row, in
the order `code`, `tests`, `rules`. The run's second table lists the same rows, one a line, as
`| <entry> | <verdict> | <list> | <row> |` under the header `| entry | verdict | implementation |
item |`.

**Exit codes.** 0 when a run ends, whatever its verdicts, when `--check` finds no violation and when
`--list` prints the records. 1 when `--expect-all-present` is given and an entry is not
`still-present`, when `--check` finds a violation, and on any failure outside an entry: an
`--extra-dir` that is not `NAME=PATH` with a lowercase name and a directory, a `--tool-dir` that is
no directory, no `docker` on `PATH` for `--image`, a record that is missing or is not TOML, a page
that cannot be read, a work directory or a `--json` file that cannot be written. 2 on a command line
that is refused (the usage, then `vep-divergence-reevaluate: error: <message>`, on stderr), which
includes a run without `--vep` or without one of `--image` and `--ensembl-vep`.

## Layout

- `src/main.rs`: the three modes, the entries' runs across `--jobs` threads and the report.
- `src/cli.rs`: the command line as `argparse` reads it, the usage line and the help text.
- `src/record.rs`: the records of a directory, TOML into values with keys in document order.
- `src/schema.rs`: the keys and types of a record, what each compare kind needs, the form of an
  input class's basis and the shape of the `implementation` table.
- `src/item.rs`: the named items of a Rust source text and the lines each spans.
- `src/page.rs`: the index page's sections and entry rows and the rule that holds each entry's
  input class to its tier and section, read as `scripts/divergences/check_divergence_order.py` reads
  the page.
- `src/check.rs`: `--check`, the `implementation` tables held to the repository's files among it.
- `src/engine.rs`: the host paths and container mounts, the input file, the Ensembl VEP cache, the
  `prepare` steps, each engine's arguments and its run.
- `src/extract.rs`: what each compare kind reads, the normalising of values and the verdict.
- `src/evaluate.rs`: one entry from its record to its verdict.
- `src/pyre/`: regular expressions as `re` reads and runs them: the parser with its messages
  (`parse.rs`), the compiler (`compile.rs`) and the matcher (`vm.rs`).
- `src/pyfmt.rs`, `src/pyjson.rs`, `src/pypath.rs`, `src/value.rs`: the CPython rules named above,
  each with tests that hold its bytes.
- `src/gzip.rs`, `src/repo.rs`, `src/tempdir.rs`, `src/dump.rs`: a gzip file's members, the
  repository root, a temporary work directory, and the two readers' output as JSON.
- `src/*_tests.rs`, `src/pyre/tests.rs`, `src/testing.rs`: the unit tests and what they share.
- `examples/readers.rs`: prints what the page reader and the TOML reader read, as JSON, so another
  reader of the same page or records can be held to them.
- `tests/runs.rs`, `tests/modes.rs`, `tests/common/`: the binary over fixture reproductions through
  stand-in engines written as shell scripts, and the modes that run nothing.
- `tests/fixtures/`: `reproductions/` and `outputs/` (one entry still present, one closed),
  `kinds/` (one reproduction a compare kind, with the outputs its stand-in engines write),
  `implementation/` (the source text the two fixture records' tables name), `toml/`
  and `page.md` with their expected dumps, the case tables (`regex_cases.jsonl`,
  `validate_cases.jsonl`, `page_cases.jsonl`, `py_cases.jsonl` and `item_span_cases.jsonl`, which
  a test of `scripts/ci/test_admission.py` reads as well) and `golden/`, the bytes each
  scenario of the tests writes, with `{{ROOT}}`, `{{FIXTURES}}` and `{{REPO}}` for the paths that
  differ from checkout to checkout. A `.gitattributes` beside this page keeps their line endings as
  committed.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
m=tools/vep-divergence-reevaluate/Cargo.toml
cargo build --release --locked --manifest-path $m
cargo test --locked --manifest-path $m
cargo clippy --all-targets --locked --manifest-path $m -- -D warnings
cargo fmt --all --manifest-path $m -- --check
```

CI runs the test, clippy and fmt commands above in its Test, Clippy and Format jobs, builds the
package on the declared `rust-version` in MSRV and audits its lockfile. `cargo test` runs no engine
and no container: a stand-in writes the output a fixture holds. One of its tests holds the
repository's own records and index page to each other, as `--check` does.

The goldens and the case tables are the contract in bytes, and `cargo test` holds the tool to them.
A change that is meant to move the contract changes them in the same commit.

## Running it

```bash
cargo run --release --locked --manifest-path tools/vep-divergence-reevaluate/Cargo.toml -- --check
tools/vep-divergence-reevaluate/target/release/vep-divergence-reevaluate \
    --image ensemblorg/ensembl-vep@sha256:<digest> --vep target/release/vep \
    --perl-cache-dir <cache root> --plugins-dir <VEP_plugins checkout> --jobs 6 \
    --expect-all-present
```

A `prepare` step may call `bgzip`, `tabix`, `bigtools` and `cargo`, which the host then needs: a
step that converts a cut of an Ensembl VEP cache builds and runs `tools/vep-storable-to-json`. One
that hard-links a cache into the stage needs `--work-dir` on the cache root's filesystem.
