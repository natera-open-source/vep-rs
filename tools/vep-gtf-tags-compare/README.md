# vep-gtf-tags-compare

Compares the transcript tags of a GTF with the transcript attributes of a JSON cache.
`vep-cache-builder --gtf` turns the `tag` values on each GTF `transcript` line into the four cache
attributes VEP keeps (`basic` or `gencode_basic` to `gencode_basic`, `gencode_primary`,
`cds_start_NF`, `cds_end_NF`); a cache converted from an Ensembl VEP cache holds the same four from
the Ensembl core database. This tool reads a GTF the way the builder does and the transcript shards
of a cache the way the runtime does, and reports for each of the four codes how many transcripts
agree on its presence and which disagree. It prints a summary on stdout and, with `--out`, writes a
JSON report; both are fixed to the byte.

```text
vep-gtf-tags-compare --gtf <GTF> --cache <JSON cache directory> [--out report.json] \
    [--max-mismatches 50]
```

- `--gtf`: the GTF. A name that ends in `.gz` is gunzipped on read.
- `--cache`: the cache directory, the one that holds `transcripts/<chromosome>/`.
- `--out`: the JSON report. Its directory has to exist, and an empty value writes no report.
- `--max-mismatches`: how many disagreeing transcripts a code keeps as examples in the report
  (default 50). A count at or below zero keeps none. Stdout shows the first five of them.

The command line is read as CPython's `argparse` reads it. An option may be shortened to any prefix
(`--g`, `--max`), takes its value as the next argument or after `=`, and a later value replaces an
earlier one. A value that starts with `-` is taken only when it is `-` alone, starts like a negative
number in ASCII digits or holds a space. `--max-mismatches` takes what `int()` takes in ASCII
digits: blanks around the number, a sign, single underscores between digits, up to 4,300 digits; a
digit outside ASCII, which `int()` reads, is refused. `-h` and `--help` print the usage and the
options.

## Contract

**The GTF side.** A line counts when it does not start with `#`, has nine tab-separated columns or
more and its third column is `transcript`. Its ninth column is split on `;`. Each part is stripped
and, where it holds a space, split at the first one into a key and a value; the value is stripped
and then freed of the double quotes at its ends. The last `transcript_id` names the line, cut at
its first `.`; a line without one is passed over, and an empty id is an id like any other. Each
`tag` whose value is one of the five kept tags adds its code. The lines of one id add up to one set
of codes.

**The cache side.** The shards are the `.json` and `.json.gz` files under `transcripts/` of the
cache directory, or under the cache directory itself where it holds no `transcripts` directory.
They are read in one order: top down, the files of a directory sorted by name, then its
subdirectories in the order the directory lists them; a symbolic link to a directory is not
followed. A shard is one JSON document. Its transcripts are the elements of a top-level array, or
of every array that is a member of a top-level object, in member order, a repeated member name
keeping its first place and its last value; any other document holds none. An element Python
counts false (`null`, `false`, zero, an empty string, array or object) is passed over; any other
element that is no object ends the run. A transcript's id is `str(value).split(".")[0]` of its
`stable_id`: a string as it stands, `null` as `None`, `true` as `True`, a number as `repr` writes
it. A transcript without an id, or with an empty one, is passed over. Its codes are the `code`
strings among the four that the objects of its `attributes` array carry; `attributes` that is
missing, false to Python, a string or an object carries none. A stable id stored in several shards
counts once, with the codes of its first copy in reading order. `attributes` that is `true` or a
number other than zero ends the run where the copy is the first of its id, and is read past on a
later copy.

**The comparison.** The transcripts both sides hold are compared code by code, in the order of
their ids by code point. A code counts the transcripts on which the two sides agree and disagree
and the transcripts that carry it on each side, and keeps its first disagreements as examples while
it holds fewer than `--max-mismatches`. A transcript only one side holds is counted and is no
disagreement.

**The report.** One object: `transcripts` (`gtf`, `cache`, `shared`, `gtf_only`, `cache_only`),
`per_code` (the four codes in the order above, each with `agree`, `disagree`, `present_in_gtf`,
`present_in_cache` and `examples`, an example a `transcript` id with `gtf` and `cache` as
booleans) and `exact_agreement`, true when no code disagrees. Stdout holds one line of transcript
counts, then one line a code and under it the code's first five examples as
`<id>: gtf=True cache=False`.

**Bytes.** Wherever a byte could be read or written in more than one way it follows a CPython rule,
named here by its CPython name. The GTF is text with universal newlines: `\n`, `\r\n` and a lone
`\r` each end a line. Both inputs are UTF-8, and a byte outside UTF-8 anywhere in them ends the
run. GTF parts are stripped with `str.strip()`'s whitespace set (U+001C to U+001F included). A
shard is read by the grammar of `json.load`: strict JSON, plus `NaN`, `Infinity` and `-Infinity`,
plus a `\u` escape of a lone surrogate, without a byte order mark, without a control character
inside a string, and without an integer of more than 4,300 digits. The report is
`json.dump(indent=1)` with `ensure_ascii` escapes and a final newline. The stdout counts are
grouped as `f"{n:,}"`, the code is padded as `{code:<16}` and its four counts as `{n:>9,}`,
`{n:>7,}`, `{n:>8,}` and `{n:>8,}`.

**Inputs.** An input whose name ends in `.gz` is read as the `gzip` module reads it: member after
member, the zero bytes that pad the file after a member skipped, a zero-byte file an empty input,
and any other byte where a member would start a corrupt stream. The GTF and every shard are read
as streams: a shard through a window that grows to the longest string the shard holds, each value
the comparison does not read skipped under the same grammar. The run holds the ids and code sets
of both sides and one window a thread, never a shard's document and never the cache.

**Exit codes.** 0 when every code agrees on every shared transcript. 1 when a code disagrees; the
summary and the report are written as on exit 0. 2 on a command line that is refused (the usage
line, then `vep-gtf-tags-compare: error: <message>`, on stderr) and when `--gtf` names nothing or
`--cache` no directory (`ERROR: [compare_gtf_tags_to_cache] --gtf must be a file and --cache a
directory` on stderr, the name in brackets a fixed part of the line). 1 as well on any other
failure: an input that cannot be read, is no UTF-8 or is a corrupt gzip stream, a shard outside
the grammar, an element or `attributes` value that ends the run. Such a run writes one line,
`ERROR: [vep-gtf-tags-compare] <file>: <detail>`, on stderr, no summary and no report. An `--out`
that cannot be written fails the same way after the summary is on stdout.

**Bounds.** A `stable_id` that is an array or an object ends the run, where Python makes an id of
its `repr`. Arrays and objects nest to any depth, where Python stops at its recursion limit. Text
is UTF-8 under any locale, where Python reads a `.gz` input and writes stdout in the locale's
encoding.

## Threads

The inputs are read across the host's available cores: the GTF is one task and every shard one
more, and the tasks share the threads. `VEP_GTF_TAGS_COMPARE_THREADS=<n>` bounds the threads at
work, `1` keeping the whole run on the calling thread; it is an environment variable and not an
option, so that the command line holds the four options alone, and a value that is not a whole
number of 1 or more is refused. Every thread count writes the same bytes: the shards are merged in
reading order once they are read, so the first copy of a stable id is the same copy on any count.

## Layout

- `src/main.rs`: the command line as `argparse` reads it, the reading tasks and the exit codes.
- `src/input.rs`: how a file is opened, the members of a gzip file and the zero padding after
  them, and the records of a text.
- `src/gtf.rs`: the `transcript` lines of a GTF, their ids and their kept tags.
- `src/json.rs`: the stream scanner of a JSON document under the grammar above.
- `src/cache.rs`: the shards in reading order, the transcripts of a shard and the merge that keeps
  the first copy of every stable id.
- `src/report.rs`: the comparison, the JSON report and the stdout lines.
- `src/pyfmt.rs`: the CPython rules named above, each with a unit test that holds its bytes.
- `src/threads.rs`: the thread bound and the indexed fan-out.
- `tests/fixtures/`: made-up GTFs and caches (`agree/`, `disagree/`, `shapes/` with every value
  shape the two readers meet, `many/` with counts past a thousand, `edge/`, and `faults/` with one
  refused input each), `cases.tsv` with one command line a row, and `expected/` with the stdout,
  stderr and report of each row to the byte. A `.gitattributes` there keeps the line endings as
  committed.
- `tests/golden.rs`: the binary over every row of `cases.tsv` on one thread, on two and on the
  host's default; a generated cache of 60 shards that store each transcript up to three times,
  over five thread counts; the thread bound; the failure line.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-gtf-tags-compare/Cargo.toml
cargo test --locked --manifest-path tools/vep-gtf-tags-compare/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-gtf-tags-compare/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-gtf-tags-compare/Cargo.toml -- --check
```

The root `.cargo/config.toml` applies, so a release binary carries the same `target-cpu` floor as
`vep` (`neoverse-v1` on aarch64 Linux, `x86-64-v3` on x86_64 Linux).

The goldens are the contract in bytes, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the goldens it moves with the binary's own output, so that the
diff of the goldens shows the change. From `tests/fixtures`, for a row of `cases.tsv` named
`agree`:

```bash
../../target/release/vep-gtf-tags-compare --gtf agree/genes.gtf --cache agree/cache \
    --out expected/agree.json > expected/agree.stdout
```
