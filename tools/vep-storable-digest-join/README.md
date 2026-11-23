# vep-storable-digest-join

Joins two digest files of `tools/vep-storable-digest` (one line per transcript of an Ensembl VEP
Storable cache: `stable_id contig start end md5_full md5_content`) on the transcript id and reports
the ids in one cache only, the ids whose `Data::Dumper` text is identical, those identical in
content alone and those that differ in content.

```text
vep-storable-digest-join --version-a 113 --dump-a 113.tsv --version-b 115 --dump-b 115.tsv \
    --out-dir <dir> [--parity-csv <table.csv> --assembly GRCh37 [--sweep-id sweep-<UTC>]]
```

- `--version-a`, `--version-b`: the names of the two caches, used in file names, keys and messages.
- `--dump-a`, `--dump-b`: the two digest files.
- `--out-dir`: where the lists and the summary go. It is made with its parents.
- `--parity-csv` with `--assembly`: a table with an `assembly` column; the row of that assembly
  gains the counts. Each needs the other.
- `--sweep-id`: the id of the run that wrote the two digest files, `sweep-` and a UTC time
  (`sweep-<8 digits>T<6 digits>Z`); written to the table's row.

The command line is read as CPython's `argparse` reads it. An option may be shortened to any prefix
that fits it alone (`--out o`), takes its value as the next argument or after `=`, and a later value
replaces an earlier one. A value that starts with `-` is taken only when it is `-` alone, starts
like a negative number or holds a space. `-h` and `--help` print the usage and the options.

## Contract

Every output is fixed to the byte: the files, standard output, the lines on standard error and the
exit code. `tests/fixtures/expected/` holds the reference bytes of 40 runs made by CPython 3.14, and
`cargo test` holds the binary to them. Wherever a byte could be written in more than one way it
follows a CPython rule, named here by its CPython name.

**Checks before anything is read.** A `--sweep-id` that is not empty and does not match
`sweep-\d{8}T\d{6}Z` in full (`\d` as in a `str` pattern: a decimal digit of any script, Unicode
16.0.0) is `ERROR: [compare_storable_cache_versions] sweep id must be sweep-<UTC>: <id>` and exit
code 2. `--parity-csv` without a non-empty `--assembly`, or the reverse, is `ERROR:
[compare_storable_cache_versions] --parity-csv and --assembly go together` and exit code 2.

**A digest file** is read as UTF-8 text with universal newlines: `\r\n` and a lone `\r` end a line
as `\n` does. Every line must hold six tab-separated columns; any other line, an empty one
included, ends the run with `ERROR: [compare_storable_cache_versions] malformed dump line in <file>:
<the line as repr() writes it>` and exit code 1, the first file's before the second's. The file's
name is written as `str(pathlib.Path(...))` gives it. Of two lines with one id the later stands.

**The lists.** Four files under `--out-dir`, each `transcript_id\tcontig\tstart\tend` and one line
an id in code point order, with the coordinates of the first file where it holds the id:
`only_in_<a>.tsv`, `only_in_<b>.tsv`, `identical_after_strip_<a>_<b>.tsv` (ids in both whose
`md5_full` differs and whose `md5_content` is equal) and `differing_<a>_<b>.tsv` (ids in both whose
`md5_content` differs). With one version for both files the second `only_in` list is written over
the first.

**The summary** is `summary_<a>_vs_<b>.json`, `json.dumps(summary, indent=2)` and a newline, with
its keys in this order: `version_a`, `version_b`, `storable_<a>_transcripts`,
`storable_<b>_transcripts`, `in_both`, `only_in_<a>`, `only_in_<b>`, `identical_full_dump`,
`identical_content_only`, `differing_content`, `content_digest_normalises` (the five keys and the
integer rule of the content digest) and `content_identical` (no id on one side only and none
differing). A key written twice keeps its first place and its last value, as a `dict`'s does. A
character outside ASCII is written as `\uNNNN`.

**Standard output** is one line, its counts with a comma every three digits:

```text
[compare_storable_cache_versions] <a>: 195,232 transcripts; <b>: 195,232; both 195,232; only <a> 0; only <b> 0; identical 0; identical after strip 194,922; differing 310
```

**The parity table** is read after the lists, the summary and the line are written. It is read as
`csv.DictReader` reads a file opened with `newline=""` in the `excel` dialect: the first record
names the columns; a field in quotes runs across commas and line ends, `""` in it is one quote and
text behind its closing quote is appended; an empty line is passed over; a short record lacks its
last values and what a long record holds past the names is dropped; of two columns with one name
the last one's value is the name's. The columns `storable_<a>_transcripts`,
`storable_<b>_transcripts`, `storable_<a>_<b>_identical_transcripts` (identical in full or in
content), `storable_<a>_<b>_differing_transcripts`, `storable_<a>_only`, `storable_<b>_only` and,
with a sweep id, `storable_<a>_<b>_sweep_id` are added behind the last column where the table lacks
them, and the one row whose `assembly` is `--assembly` takes their values. No such row, or several,
is `ERROR: [compare_storable_cache_versions] <table> has <n> rows for <assembly>` and exit code 1,
with the table left as it stands. The table is written as `csv.DictWriter` writes it with
`lineterminator="\n"` (a field in quotes when it holds a comma, a quote, `\n` or `\r`, its quotes
doubled; a record of one empty field as `""`) to the file `Path.with_suffix(".csv.tmp")` names
beside it, which then takes the table's place.

**Cores and memory.** The two digest files are read on two threads. Both are held in memory, one
map from id to the five other columns each: the join needs every id of both.

**Outside the contract.** `argparse`'s own usage and error text is not reproduced beyond the exit
code 2 and the message after `error: `; the usage line names this tool. Where CPython ends on a
traceback the tool prints one line `ERROR: [vep-storable-digest-join] <detail>` and exits with 1
as CPython does: a digest file or a table that cannot be read or is not UTF-8, an output directory
that cannot be made, a file that cannot be written, a table without an `assembly` column, a field
of more than 131,072 characters, a table path without a name, a version or an assembly that is not
UTF-8.

## Layout

- `src/argparse.rs`: a command line of long options with one value each, as `argparse` reads it.
- `src/py.rs`: text as CPython writes and reads it (`pathlib`'s path text, universal newlines,
  `repr()`, a JSON string, a grouped integer).
- `src/csv.rs`: the `excel` dialect's reader and writer and `DictReader`'s rows.
- `src/join.rs`: the join; `src/main.rs`: the command line.
- `tests/golden.rs` with `tests/fixtures/`: the binary against the reference bytes of 40 runs.
- `tools/vep-storable-transcript-classes` compiles `src/argparse.rs` and `src/py.rs` by `#[path]`.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest.

```bash
cargo build --release --locked --manifest-path tools/vep-storable-digest-join/Cargo.toml
cargo test --locked --manifest-path tools/vep-storable-digest-join/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-storable-digest-join/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-storable-digest-join/Cargo.toml -- --check
```

The reference bytes are the contract, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the files it moves under `tests/fixtures/expected/<case>/` with
the binary's own output, so that the diff of the fixtures shows the change. An output whose name
holds a character outside printable ASCII, a double quote or a backslash is kept as
`<case>/files/<plain name>`, and `<case>/files.tsv` spells the path the run writes it to.
