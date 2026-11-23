# vep-cache-diff

The structural comparison of two JSON caches. It reads two cache directories, a baseline and a
builder's, pairs their transcript shards on the shard's path and their variation files on the
chromosome, and reports how many transcripts and variations the two caches share, how many of those
agree, and under which fields the others differ; the two `info.json` are compared key by key. It
writes one JSON report and a summary on stdout; both are fixed to the byte.

```text
vep-cache-diff --baseline-dir <cache> --builder-dir <cache> --output-report <report.json>
```

- `--baseline-dir`: the cache the other is held against, for instance one converted from an Ensembl
  VEP cache.
- `--builder-dir`: the cache under test, for instance one `vep-cache-builder` wrote.
- `--output-report`: the report. Its parent directories are created.

The command line is read as CPython's `argparse` reads it. An option may be shortened to any prefix
that fits it alone (`--ba`, `--bu`, `--o`; `--b` fits two and is refused), takes its value as the next
argument or after `=`, and a later value replaces an earlier one. A value that starts with `-` is
taken only when it is `-` alone, starts like a negative number or holds a space. `-h` and `--help`
print the usage and the options.

## Contract

**Transcripts.** A shard is an entry under `transcripts/` of a cache whose name ends in `.json` or
`.json.gz`, at any depth, as `Path.rglob` lists them: an entry of any type counts, a directory is
entered unless it is a symbolic link, and a name that starts with `.` is a name like any other. The
two sides are paired on the shard's path under `transcripts/` with a `.gz` suffix removed, so a
gzipped shard pairs with a plain one; where one side holds both under one name, the gzipped file is
read. A link that leads nowhere is a shard without records.

A shard is one JSON document. Its records are the objects among the items of a top-level array. A
top-level object with a member `transcripts` stands for that member: an array's objects, or an
object's members. Without such a member, and for the members of an object under `transcripts`, each
member gives the objects of its array, or is one record where it is an object; a repeated member
name keeps its first place and its last value. Any other document holds no record. A record's id
is its `stable_id`, or its `id` where `stable_id` is missing or false to Python (`null`, `false`,
zero, an empty string, array or object); a record without an id is passed over, and of two records
of one id in a shard the later one counts.

Under each shard name the ids of the two sides are compared. An id one side alone holds is counted.
Two records of one id are compared in this order, and each difference is counted under its name:

- the members `start`, `end`, `strand`, `biotype`, `source`, `gene_stable_id`, `gene_symbol`,
  `hgnc_id` and `canonical`, as Python's `!=` compares what `json.load` returns: a missing member
  is `null`, `1`, `1.0` and `true` are equal, an object equals an object of the same members in any
  order, and `NaN` differs from `NaN` as a member and equals it inside an array or object;
- `exon_count`, `intron_count` and `mapper_pair_count`: the lengths of `sorted_exons` and `introns`
  where they are arrays and of `mapper` where it is an array or an object, zero otherwise, each read
  from the record's feature cache or, where the cache's member is missing or false to Python, from
  the record itself. The feature cache is the member `variation_effect_feature_cache` where Python
  counts it true, else the member `vefc`, else an empty object;
- the members `introns`, `sorted_exons`, `five_prime_utr`, `three_prime_utr`, `codon_table`,
  `protein_features` and `seq_edits` of the feature cache, compared as the nine members are: under
  `vefc_<member>_count` where both sides hold an array or an object and the two lengths differ,
  under `vefc_<member>` otherwise.

**Variations.** A variation file is `variations/<chromosome>/all_vars.gz`, as `Path.glob` lists
them: `<chromosome>` is a directory or a link to one. Its lines are tab-separated cells in the
order of the `variation_cols` of the cache's `info.json`: an array of names, the characters of a
string or the keys of an object, as `list()` reads them, and no column where the member is missing
or false to Python or the cache has no `info.json`. A line keeps as many cells as the cache declares
columns; an empty cell and a cell that is `.` alone are no value, and of two columns of one name the
last cell with a value counts. A line without a `variation_name` is passed over, and of two lines of
one name in a file the later one counts.

Under each chromosome the names of the two sides are compared, whatever the `chr` column holds. A
name one side alone holds is counted. Two records of one name are compared on `failed`, `somatic`,
`start`, `end`, `allele_string`, `strand`, `clin_sig`, `phenotype_or_disease`, `clin_sig_allele`,
`clin_sig_ref_allele`, `clinical_impact`, `pubmed` and `var_synonyms`, each difference counted
under the column's name, and together, under `frequencies`, on the frequency columns: those the
baseline declares, other than the thirteen, `chr` and `variation_name`, that the builder declares
as well.

**info.json.** The two documents match where Python's `==` holds for them. Where it does not, each
key either holds, in the order of its code points, whose two values differ is listed with `str()`
of both values cut to 200 characters; a key one side lacks reads as `None` there. A cache without
the file is named, and two caches without it match.

**The report.** One object with, in this order, `baseline_dir` and `builder_dir` (the two
directories as `Path.resolve()` gives them), `transcripts`, `variations` and `info_json`. The two
comparison blocks open with `total_baseline`, `total_builder`, `shared`, `only_baseline`,
`only_builder`, `matched`, `differed` and `field_mismatches`, and close with `match_rate`
(`matched / shared` rounded to six places, over 1 where nothing is shared). `field_mismatches`
holds each name that was counted, the largest count first. Among equal counts the transcript names
stand in the order the comparison first met them, going through the shard names, the ids under
each and a record's differences in sorted order; the variation names stand in the order of the
first chromosome each was met in, and within one chromosome in the order listed above. The
`transcripts` block ends with `regions_with_diffs`: the first fifty shard names, in sorted order,
under which an id is held by one side alone or differs, each with its counts and the first ten ids
of each side alone. The `variations` block holds `frequency_columns_compared` before its rate.
`info_json` holds `matches`, `diffs` and, where neither cache has the file, a `note`. Stdout holds
the two directories, one line of counts for each comparison, `matches=True` or `matches=False`, and
the report's path.

**Bytes.** Wherever a byte could be read or written in more than one way it follows a CPython rule,
named here by its CPython name. Every path option goes through `Path(arg).resolve()` without
`strict`. Shard names and chromosomes sort as `sorted()` sorts `str`: by code point, a byte of a
file name outside UTF-8 standing as the lone surrogate `surrogateescape` gives it and written
`\udcXX` in the report. A shard and an `info.json` are read by the grammar of `json.load`: strict
JSON, plus `NaN`, `Infinity` and `-Infinity`, plus a `\u` escape of a lone surrogate, without a byte
order mark, without a control character inside a string, and without an integer of more than 4,300
digits; an integer keeps every digit. A variation file is text with universal newlines: `\n`,
`\r\n` and a lone `\r` each end a line. The report is `json.dump(indent=2)` with `ensure_ascii`
escapes and a final newline. A rate is `round(x, 6)` written as `float.__repr__` in the report and
as `f"{x:.4f}"` on stdout. A value of `info.json` is written as `str()` writes it: a string as it
stands, anything else as `repr()`, a string inside it in the quotes `repr()` chooses and with the
escapes `str.isprintable()` asks for under Unicode 16.0.0.

**Inputs.** A file whose name ends in `.gz` is read as the `gzip` module reads it: member after
member, the zero bytes that pad the file after a member skipped, a zero-byte file an empty input,
and any other byte where a member would start a corrupt stream. Neither cache is ever held in
memory. A shard is read whole, one pair of shards a thread, and of its records only the members
the comparison reads are kept; every other value is read past under the same grammar. A variation
file is read as a stream. Each line is reduced to its name and the cells the comparison reads (to
its name alone where the other cache has no record of the chromosome to hold it against) and
routed, by a hash of the name, into one of 256 partitions of its side; a chromosome that holds more
than 256 MiB of such records writes its partitions to a temporary file and joins them one at a
time, in several passes where a partition is itself past that bound, and a record that a later
line of its name replaces is dropped where a partition is read back. The temporary files stand in a
directory `vep_cache_diff_<process id>_<8 characters>` under the system's temporary directory
(`TMPDIR`), which the run removes; a run that stays under the bound makes none. A chromosome's
file holds its names and compared cells of both sides once and is removed when the chromosome is
joined, so the directory needs that much room for each chromosome a thread is at work on.

**Exit codes.** 0 when the report and the summary are written. 2 on a command line that is refused
(the usage line, then `vep-cache-diff: error: <message>`, on stderr). 1 when a cache directory is
missing (`ERROR: [diff_caches] baseline directory not found: <path>` or `... builder directory not
found: <path>` on stderr, the name in brackets a fixed part of the line, the baseline looked at
first) and on any other failure: a shard or an `info.json` outside the grammar or outside UTF-8, a
corrupt gzip stream, a file that cannot be read, an `info.json` that is no object, `variation_cols`
that is a number or `true`, a shared id whose feature cache is no object, an id that is an array or
an object, a report that cannot be written. Such a run writes one line,
`ERROR: [vep-cache-diff] <file>: <detail>`, on stderr and no report; stdout holds the lines of the
steps that ended before the failure, so a failure among the shards leaves `Comparing
transcripts...` as its last line.

**Bounds.** An id that is a number or `true` ends the run, where Python keeps it as a key. An item
of `variation_cols` that is no string ends the run, where Python compares a column named by a
number, a boolean or `null`. A value the comparison reads may nest 1,000 arrays and objects deep;
a value it reads past nests to any depth, where Python stops when its stack is used up. Text is
UTF-8 under any locale, where Python reads `info.json` and the variation files in the locale's
encoding, and a path outside UTF-8 is written on stdout as its bytes, which Python does under the
C, POSIX and C.UTF-8 locales alone. Among variation names of equal counts Python's order follows
the iteration of a `set` and can change from run to run; the order above is one of the orders
Python gives, and the only one where its order cannot change.

## Threads and memory

The comparison runs across the host's available cores: every pair of shards is one task and every
chromosome one more, the largest chromosomes first. `VEP_CACHE_DIFF_THREADS=<n>` bounds the threads
at work, `1` keeping the whole comparison on the calling thread.
`VEP_CACHE_DIFF_BUFFER_BYTES=<n>` sets the bytes of variation records one chromosome holds in
memory before it writes them to its temporary file (default 268435456). Both are environment
variables and not options, so that the command line holds the three options alone, and a value that
is not a whole number of 1 or more is refused before anything is read. Every thread count and every
bound writes the same bytes. What a run holds in memory grows with the threads at work: each
chromosome a thread compares holds up to its bound of records, and at most twice that while a
buffer grows, and each pair of shards its two texts.

## Layout

- `src/main.rs`: the run: the three comparisons in order, the stdout lines, the report's blocks and
  the exit codes.
- `src/cli.rs`: the command line as `argparse` reads it, the usage line and the help text.
- `src/paths.rs`: `Path.resolve()`, the shards and the variation files as `rglob` and `glob` list
  them, and the two sides paired on their names.
- `src/input.rs`: the members of a gzip file and the zero padding after them, and the records of a
  text.
- `src/json.rs`: the reader of a JSON document under the grammar above, the values it keeps, and
  Python's `==`, `!=` and truth for them.
- `src/transcripts.rs`: the records of a shard, what is kept of one, the comparison of two records
  and the counts of one shard name.
- `src/variations.rs`: the columns, the record of a line, the partitions, the temporary file and the
  join of one chromosome.
- `src/info.rs`: the comparison of the two `info.json`.
- `src/report.rs`: the JSON emitter and the stdout line of a comparison.
- `src/pyfmt.rs` and `src/unicode.rs`: the CPython rules named above, each with a unit test that
  holds its bytes, and the code points `str.isprintable()` refuses.
- `src/threads.rs`: the two environment variables and the indexed fan-out.
- `tests/fixtures/`: made-up cache pairs under `caches/` (`identical/`, `transcripts/` with every
  container, every mismatch name and more shard names than the report lists, `variations/` with
  columns in two orders, every mismatch name and every line end, caches with `info.json` on one
  side, `faults/` with one refused input each, `bounds/` with the inputs of the Bounds paragraph),
  `cases.tsv` with one pair a row, and `expected/` with the stdout, stderr and report of each row to
  the byte. A `.gitattributes` there keeps the bytes as committed.
- `tests/golden.rs`: the binary over every row of `cases.tsv` under four bounds of threads and
  memory, a generated pair of 60 shards and five chromosomes whose counts are known from how it is
  made, and the temporary directory.
- `tests/command_line.rs`: the binary's exit code and bytes for shortened options, refused command
  lines, help, the report's path, the two environment variables, names outside UTF-8, caches made of
  links and the inputs of the Bounds paragraph. `tests/common/` holds what the test files share.

## Build and test

Its own Cargo root outside the repository workspace; every command names the manifest:

```bash
cargo build --release --locked --manifest-path tools/vep-cache-diff/Cargo.toml
cargo test --locked --manifest-path tools/vep-cache-diff/Cargo.toml
cargo clippy --all-targets --locked --manifest-path tools/vep-cache-diff/Cargo.toml -- -D warnings
cargo fmt --all --manifest-path tools/vep-cache-diff/Cargo.toml -- --check
```

The root `.cargo/config.toml` applies, so a release binary carries the same `target-cpu` floor as
`vep` (`neoverse-v1` on aarch64 Linux, `x86-64-v3` on x86_64 Linux).

The goldens are the contract in bytes, and `cargo test` holds the binary to them. A change that is
meant to move the contract rewrites the goldens it moves with the binary's own output, so that the
diff of the goldens shows the change. From `tests/fixtures`, for a row of `cases.tsv` named
`transcripts`, with the directory's own path put back as `{{FIXTURES}}` and the report's as
`{{REPORT}}`:

```bash
../../target/release/vep-cache-diff --baseline-dir caches/transcripts/baseline \
    --builder-dir caches/transcripts/builder --output-report report.json \
    | sed -e "s#$PWD/report.json#{{REPORT}}#" -e "s#$PWD#{{FIXTURES}}#g" > expected/transcripts.stdout
sed -e "s#$PWD#{{FIXTURES}}#g" report.json > expected/transcripts.json && rm report.json
```

## Harness

`scripts/validation/compare_cache_outputs.sh` and `scripts/validation/validate_builder_cache.sh` run
the tool for their structural comparison and read its report. `compare_cache_outputs.sh` runs this
tree's release build, `tools/vep-cache-diff/target/release/vep-cache-diff`, and builds nothing: a
build that is missing or is not an executable stops the script before its first engine run, with
one line that names the build command above. `validate_builder_cache.sh` builds the tool through
cargo, as it builds the cache builder it runs, and a build that fails stops it before a cache is
generated. `scripts/validation/test_cache_comparison_callers.sh` holds both scripts to this, and
runs the comparison step of `compare_cache_outputs.sh` once on the tree's own release build where
there is one.
