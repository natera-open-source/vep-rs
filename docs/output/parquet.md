# Parquet format

This page describes the Parquet output, a vep-rs extension with no Ensembl VEP counterpart: its layout, its schema, the two row shapes and the writer's choices. Selected with `--output_format parquet`. The run writes a TSV intermediate at `<output>.tsv.tmp`, one row per allele x transcript consequence in VEP's row order, and post-processes it to Parquet through the `duckdb` CLI. The intermediate is deleted on success unless the hidden `--keep_intermediate` flag is given.

The output is a directory, so `-o STDOUT` is refused (`ERROR: [output] --output_format parquet needs a file path: Parquet cannot be written to STDOUT`). `--fields` and `--compress_output` are refused with it too: the column set is fixed, and the files carry their own compression. `--no_headers` leaves the intermediate's column line in place, since DuckDB reads the schema from it.

## Requirements

The `duckdb` CLI, 1.3.2 or newer, must be on `PATH`. Install it with `brew install duckdb` on macOS, or from the CLI archive at <https://duckdb.org/docs/installation/> (on Linux the `duckdb-cli-linux-*.tar.gz` release asset; 1.x releases shipped `duckdb_cli-linux-*.zip`). vep-rs checks for it after the FASTA, the plugins and the variation cache are loaded and before the first variant is read. It fails there when `duckdb` is missing or too old, so no annotation is spent on a run it cannot finalize.

## Directory layout

A Hive-partitioned directory keyed by `chrom`, with the partition column also stored inside the files:

```
out.parquet/
├── chrom=1/data_0.parquet
├── chrom=2/data_0.parquet
└── chrom=X/data_0.parquet
```

Read it with `hive_partitioning=false` to keep `chrom` as the stored VARCHAR. Hive type inference would otherwise read `chrom=21` as a number.

```sql
read_parquet('out.parquet/**/*.parquet', hive_partitioning=false)
```

## Schema

Five key columns, then exactly the default `--tab` columns in `--tab` order:

| Column                                        | Type    | Notes                                     |
| --------------------------------------------- | ------- | ----------------------------------------- |
| `chrom`                                       | VARCHAR | input chromosome name                     |
| `pos`, `end`                                  | BIGINT  | the record's VEP start and end            |
| `ref`, `alt`                                  | VARCHAR | this row's reference and alternate allele |
| `Uploaded_variation` ... `Existing_variation` | VARCHAR | the 13 default columns                    |
| `IMPACT`, `FLAGS`, flag and plugin fields     | VARCHAR | `-` is stored as NULL                     |
| `DISTANCE`                                    | INTEGER |                                           |
| `STRAND`                                      | TINYINT |                                           |

In the allele columns (`ref`, `alt`, `Allele`, and `GIVEN_REF` and `USED_REF` under `--use_transcript_ref`) a `-` is a real allele, a deletion's ALT or an insertion's REF, and is kept.

## Row shape (`--parquet_shape`)

`nested`, the default, writes one row per (record, allele). `Uploaded_variation`, `Location`, `Allele` and `Existing_variation` stay scalar; every per-consequence column is a LIST in the flat row order, aligned index-wise across columns. Two input records at one site stay two rows when they differ in ID, END or another per-record column; identical ID-less duplicate lines merge into one row. Query a LIST column with `UNNEST` or `list_contains`:

```sql
SELECT chrom, pos, ref, alt
FROM read_parquet('out.parquet/**/*.parquet', hive_partitioning=false)
WHERE list_contains(Feature, 'ENST00000366667')
  AND list_contains(IMPACT, 'HIGH');
```

`flat` writes one row per allele x transcript consequence, every column scalar, so predicates push down directly:

```sql
SELECT chrom, pos, ref, alt, Consequence, IMPACT
FROM read_parquet('out.parquet/**/*.parquet', hive_partitioning=false)
WHERE Feature = 'ENST00000366667' AND IMPACT = 'HIGH';
```

## Writer options

Every option is chosen for the reader, not the writer:

| Option           | Value                                                                                             |
| ---------------- | ------------------------------------------------------------------------------------------------- |
| Page format      | Parquet v2                                                                                        |
| Compression      | ZSTD, level 9                                                                                     |
| Row group size   | 122,880 rows; the hidden `--parquet_row_group_size` overrides it                                 |
| Sort order       | `(chrom, pos, ref, alt)`, so the `pos` statistics of a row group are tight                        |
| Dictionary       | Every column; the dictionary size limit equals the row group size                                 |
| Bloom filters    | Every dictionary-encoded column, at a 0.1% false-positive rate                                    |
| Partition column | `chrom`, stored inside the files as well as in the directory names                                |
| Footer metadata  | `vep_rs_version`, `vep_api_version`, `assembly`, `cache`, `command_line`, every `source_<key>` of the cache |

DuckDB writes a Bloom filter only for a dictionary-encoded column, so the dictionary size limit is set to the row group size, which a flat column's distinct values can never exceed. The deprecated `DICTIONARY_COMPRESSION_RATIO_THRESHOLD` option is passed as well; DuckDB 1.2 and later ignore it. A column that is NULL throughout a row group has neither a dictionary nor a Bloom filter.

In the nested shape a LIST column carries a Bloom filter only when its elements are dictionary-encoded, which needs their distinct count in the row group to fit the limit (DuckDB 1.2 wrote no Bloom filter for nested types at all). So a column whose elements are mostly distinct, such as `Feature` or `HGVSc`, usually has none. The scalar columns, whose distinct count cannot exceed the row count, keep theirs.

`parquet_round_trips_to_the_tab_output` in `crates/vep-cli/tests/format_parity.rs` asserts, on every corpus, that both shapes round-trip to the `--tab` rows. On the flat shape it also asserts that every populated column of every row group carries a dictionary, a Bloom filter and statistics, and that the footer holds `vep_rs_version`, `assembly` and `command_line`. The page version, the compression, the row group size, the sort order and the false-positive rate are not asserted.

## Round-trip to `--tab`

The adapter script projects the tab columns back out: NULL to `-`, integers to text, nested rows unnested. Rows come out in Parquet order; compare against a `--tab` file as a multiset of lines.

```bash
scripts/adapters/parquet_to_vep_tab.py <dir|file> [--output tab.txt]
```

## Debug flags

Both flags are hidden from `--help`.

`--keep_intermediate` preserves the TSV intermediate after the DuckDB step for inspection.

`--parquet_row_group_size N` sets the rows per row group, rounded up to a multiple of 2,048, DuckDB's vector size. The dictionary size limit follows the rounded value; otherwise a row group with a distinct value in every row would lose its dictionary and its Bloom filter. The flag exists so a small input can be written as several row groups under test:

```
parquet_row_groups_prune_on_bloom_filters
```

This test writes a corpus at 2,048-row groups and checks with `parquet_bloom_probe` that an absent `Feature` inside the column's value range is excluded by every row group, and a single-row `Feature` by every row group but its own.
