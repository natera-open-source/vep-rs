# Binary interval queries and verification

`vep-cache-converter convert-plugin-data` converts a tabix-indexed plugin
file into a `.vpd` data file and a `.vpdi` index. Select the columns used to
build the tabix index: `--chr-col`, `--start-col` and `--end-col` count from
zero. Without `--end-col`, or when it names the start column, every record
is a point. `--zero-based` converts
the start to one-based coordinates; interval ends stay inclusive at their
original value.

```bash
vep-cache-converter convert-plugin-data \
  --input spans.tsv.gz --output spans.vpd \
  --end-col 2 --verify-all
```

A record matches when its start is at or before the query end and its end
is at or after the query start. For the repository's
[`spans.tsv.gz` fixture](../../crates/vep-builtins/tests/data/spans.tsv.gz),
`21:250-255` returns both `span150_260` and `span180_320`.
`21:320-320` returns the latter alone, and `21:321-399` returns nothing.
Single queries and batch lookups preserve the same records, their order
and their multiplicity, including identical source rows.

The binary files retain layout version 1. Opening a store scans record
bounds and field lengths once, without decoding column text, and builds
an in-memory prefix of maximum interval ends. The prefix adds eight bytes
per indexed start. Queries search both starts and reach, then filter
exact overlaps. Point-only queries restrict candidates to starts within
the query. Invalid record bounds, index order, counts and offsets fail
at open.

## What verification checks

Conversion queries the original tabix file and compares coordinates,
ordered records and every column with the binary store. Probes include
record starts, interior midpoints, endpoints, complete spans, adjacent
bases and gaps between covered regions. Each probe runs alone and in an
isolated batch, then in batches of at most 64 probes, reversed batches
with a repeated query, and batches containing alternate probes.

By default, probes come from at most 1,000 evenly spaced source records,
plus a point after each chromosome's final interval. `--verify-all`
derives probes from every source record; it does not query every base.
`--no-verify` skips these comparisons. Verification requires the source's
tabix index. Failures identify the query mode, chromosome and range,
then the differing count, record bounds or column values.

Existing version 1 stores remain readable. Reconvert an older store if
same-start records with different ends need their original source order:
earlier converters sorted those rows by end. Newly converted stores keep
source order among equal starts.
