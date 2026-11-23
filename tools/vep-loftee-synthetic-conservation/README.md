# vep-loftee-synthetic-conservation

A writer of made-up conservation files for the LoFTEE plugin on GRCh37. The plugin's `master`
branch reads a score for a base from one file and scores for an exon from another, and vep-rs
is compared with Ensembl VEP on runs of the plugin over the same two files. This tool reads one
cache in the JSON layout `tools/vep-storable-to-json` writes, through the cache reader of
[vep-seeded-variants](../vep-seeded-variants/README.md), and writes both files for every exon of
every transcript of the cache. Every value is made up, a function of a seed and the row's key,
and every file carries `synthetic` in its name.

```text
vep-loftee-synthetic-conservation --cache <dir> --seed <n> --out-dir <dir>
    [--threads <n>] [--per-sequence <file.tsv>] [--force]
    [--skip-transcripts-without-exons]
vep-loftee-synthetic-conservation --version
```

## What it writes

A run writes three files into the directory `--out-dir` names:

| File | What it holds |
| --- | --- |
| `GERP_scores.synthetic.txt.gz` | one row for each base of an exon: the sequence, the position and a score |
| `GERP_scores.synthetic.txt.gz.tbi` | the tabix index of those rows |
| `phylocsf_gerp.synthetic.sql` | a SQLite database of two tables, with a row of each for every exon |

The plugin's own source states what it reads of each file:

- `LoF.pm:138-141` keeps the command `tabix <gerp_file>`. `gerp_dist.pl:92-96` runs it on a
  region, splits each row it returns on tabs and adds up the third field, and
  `gerp_dist.pl:119-121` reads the third field of one position.
- `gerp_dist.pl:129` runs `SELECT * FROM gerp_exons where transcript_id = ? AND exon_num = ?`
  with the transcript's stable id and the exon's place in the transcript, counted from 1.
  `gerp_dist.pl:133` reads the column `gerp` of the row.
- `LoF.pm:570` runs `SELECT * FROM phylocsf_data WHERE transcript = ? AND exon_number = ?`
  with the stable id and the exon number of the output row. `LoF.pm:293-296` read the columns
  `corresponding_orf_score` and `max_score` of the row.

### The table

A row is `sequence<TAB>position<TAB>score`, the position counted from 1. The rows cover the
bases that lie in an exon of any transcript of the cache, whatever its biotype, each base once:
the sequences in the cache's order, and a sequence's bases in ascending position. The file is
BGZF-compressed. Its index is the one `tabix -s 1 -b 2 -e 2` writes: the generic format, the
sequence in column 1 and the position in column 2 as begin and end. An exon base above position
536,870,911, the last one a region read through a tabix index names, fails the run.

### The database

```sql
CREATE TABLE gerp_exons (transcript_id TEXT, exon_num INTEGER, gerp REAL);
CREATE TABLE phylocsf_data (transcript TEXT, exon_number INTEGER,
                            corresponding_orf_score REAL, max_score REAL);
CREATE INDEX gerp_exons_key ON gerp_exons (transcript_id, exon_num);
CREATE INDEX phylocsf_data_key ON phylocsf_data (transcript, exon_number);
```

The two tables are the ones `scripts/fixtures/loftee_conservation_database.py` writes for a
golden corpus. Each holds one row for every exon of every transcript record, the exons of a
record numbered from 1. Each table has an index on the two columns its query names, since a
query of a genome-wide table without one reads the whole table. An index changes no row a query
returns.

A transcript record is one stable id, start, end and `dbID`, as the cache reader tells records
apart. A record filed in several region files has its rows once. A stable id under two `dbID`s
is two records, and both have their rows: an exon number that both records have stands in two
rows under one key, and the two rows hold equal values.

## The values

Each value is drawn from the stream that the seed, the name of the value's table and the row's
key name, so it does not move when the cache gains or loses a transcript. A value has at most
three decimals. Ensembl VEP prints a number with up to fifteen significant digits, so both
engines print such a value as the same text.

| Value | Key | Rule |
| --- | --- | --- |
| the score of a base | sequence, position | a step of 0.001, 0.01, 0.1 or 1, then a multiple of it from -6 to 9 |
| the total of an exon, `gerp` | stable id, exon number | a multiple of 0.01 from -300 to 900 |
| the PhyloCSF scores of an exon | stable id, exon number | one of three outcomes, each as likely as the others |

- The scores of the bases average 1.5, so the bases of a run of about 120 add up to the 180
  that `LoF.pm:268` holds the walk's total against. A score is written as its shortest decimal:
  `-0.5`, `3`, `2.125`.
- The total of an exon is drawn from its own key and never from the scores of the exon's bases,
  so a total read from the other file shows as a difference.
- The first PhyloCSF outcome is an annotated score that is not negative, an exact 0 one time in
  eight, which `LoF.pm:295-296` leave without a flag. The second is a negative annotated score
  under a positive highest score, flagged `PHYLOCSF_UNLIKELY_ORF`. The third is a negative
  annotated score under a highest score that is negative or, one time in three, an exact 0,
  flagged `PHYLOCSF_WEAK`.
- Both PhyloCSF scores of an exon are multiples of one step. An annotated score lies from -300
  to 400, a highest score is at most 400, and the highest score is never below the annotated
  one.

## Options

- `--cache` names the converted cache directory, the one `--json_cache` takes.
- `--seed` is a whole number from 0 to 2^64 - 1.
- `--out-dir` names the directory of the three files, which is created where it is absent. A
  file is written under a temporary name beside its own and renamed when all three are whole,
  so a run that fails before that leaves none of them.
- `--threads` is the number of threads that read the cache and of threads that compress the
  table; the default is one a processor. The same cache and seed give the same table and the
  same index, byte for byte, on any number of threads.
- `--per-sequence` writes a table of the sequences: the header `sequence`, `transcripts`,
  `exons`, `exon_bases`, then one line a sequence in the order of the table's rows. The
  file is written after the three files are in place: a run that cannot write it fails and
  leaves them.
- `--force` overwrites a file that exists. Without it one of the three files, or the file
  `--per-sequence` names, that exists fails the run before a transcript is read. So does a
  directory under one of those names, with the option or without it.
- `--skip-transcripts-without-exons` leaves out the transcript records the cache reader gives
  no exon for: a record that places an exon outside its own span, and one stored across the
  origin of a circular sequence. Without the option such a record fails the run, with a line
  that names its sequence and stable id and counts the records of its kind.
- `--version` prints the tool's version.

## What a run prints

A run prints tab-separated lines on standard output, a key and its value, in this order:

| Line | Its value |
| --- | --- |
| `seed` | the seed |
| `sequences` | the sequences of the cache |
| `transcripts` | the transcript records the cache reader gives |
| `transcripts_without_exons` | the records left out; the line is printed under `--skip-transcripts-without-exons` alone |
| `exons` | the exons of the records that have one, which is the rows of each table of the database |
| `exon_bases` | the bases that lie in an exon, each counted once, which is the rows of the table |
| `GERP_scores.synthetic.txt.gz` | `rows=<n>` and `rows_sha256=<hex>` |
| `phylocsf_gerp.synthetic.sql` | `rows=<n>`, `gerp_exons=<n>`, `phylocsf_data=<n>` and `rows_sha256=<hex>` |

The `rows_sha256` of the table is the sha256 of its rows as text, before compression. The
bytes of a SQLite file depend on the library that writes it, so the `rows_sha256` of the
database is taken over its rows, read back from the written file. Its lines are
`gerp_exons<TAB>transcript_id<TAB>exon_num<TAB>gerp` in the order of those three columns, then
`phylocsf_data<TAB>transcript<TAB>exon_number<TAB>corresponding_orf_score<TAB>max_score` in the
order of those four. A score stands as the shortest decimal that reads back as the same double,
and `rows` counts the rows of both tables.

## Exit codes

A run that cannot read its command line prints
`ERROR: [vep-loftee-synthetic-conservation] <what>` and the usage and exits 2. A run that fails
prints one such line and exits 1: a cache that cannot be read, a file that exists, a transcript
record without an exon, or an exon base above the last position of a tabix index.

## Building and testing

```bash
cargo build --release --locked \
    --manifest-path tools/vep-loftee-synthetic-conservation/Cargo.toml
cargo test --locked \
    --manifest-path tools/vep-loftee-synthetic-conservation/Cargo.toml
cargo clippy --all-targets --locked \
    --manifest-path tools/vep-loftee-synthetic-conservation/Cargo.toml \
    -- -D warnings
cargo fmt --all \
    --manifest-path tools/vep-loftee-synthetic-conservation/Cargo.toml \
    -- --check
```

The tests write made-up caches into a scratch directory and read the three files of a run
back. One of them also reads regions of the table with the `tabix` program, where the machine
has one under `~/.local/bin` or on the search path.
