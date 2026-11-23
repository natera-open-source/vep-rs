# vep-genome-table

The reader, writer and checker of the genome table under `tests/genomes`, and the digest its
rows carry. [docs/genome-checks.md](../../docs/genome-checks.md) describes the four files of a
table directory, their columns and the genome key.

```text
vep-genome-table check <dir> [--complete]
vep-genome-table row <dir> <genome key> [--field <column>]
vep-genome-table keys <dir> [--set <set>] [--outcome <prefix>]
vep-genome-table put <dir> <row file>...
vep-genome-table digest <output file>
```

Build it with `cargo build --release --locked --manifest-path tools/vep-genome-table/Cargo.toml`.
The library `vep_genome_table` carries the same operations for a Rust caller.

## Subcommands

- `check` holds each of the four files to its header and to the grammar of every field. No key
  stands twice, each row is in the file its set names, the rows stand in bytewise key order, and
  every `differs:` and `excluded:` value is declared in `outcomes.tsv`. An `exact` row carries a
  digest and records above zero. With `--complete`, every set of `sets.tsv` has exactly its
  `caches` rows and every row names a set of `sets.tsv`. A pass prints `<n> rows in <m> sets: ok`.
- `row` prints the header line and the row of one genome key, or one field of it.
- `keys` prints one key a line in key order. `--set` keeps the rows of one set; `--outcome` keeps
  the rows whose outcome starts with the prefix (`exact`, `differs`, `excluded:`).
- `put` reads row files, each the table's header line and rows, plain or gzip. A row takes the
  place of the row with its key, the last one given winning, and both row files are rewritten in
  key order. The compressed file is written with no name and a zero time, so equal rows give
  equal bytes. A row that fails the field rules refuses the whole command and nothing is written.
  It prints `<n> rows in <m> sets: <a> added, <r> replaced`.
- `digest` prints the SHA-256 of a default-format output, plain or gzip, with every line that
  starts with `##` left out.

## Exit codes and errors

| Code | Meaning                                                          |
| ---- | ---------------------------------------------------------------- |
| 0    | the command did what it says                                     |
| 1    | a refusal: a failed check, a key no row has, an unreadable file  |
| 2    | a command line the tool cannot read; the usage follows the error |

Every error is one line on standard error, `ERROR: [vep-genome-table] <what>`. A failure of a
table file reads `ERROR: [vep-genome-table] <file>:<line>: <what>`, one line a failure, in file
order and then line order. Lines count from 1 in the file's text, decompressed for
`bacteria.tsv.gz`; line 0 marks a file that cannot be read or written. `row` and `keys` refuse a
table whose row files fail `check`'s rules for them.
