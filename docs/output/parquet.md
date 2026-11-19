# Parquet format

The Parquet output is a vep-rs extension with no Ensembl VEP counterpart: a directory of typed parts to query with DuckDB, Polars, Spark or Athena.

Selected with `--output_format parquet`, the run writes the Parquet files itself, with no intermediate file and no external tool. This page covers:

- the directory the run writes;
- the typed schema and its inverse to the `--tab` format;
- the row order;
- the writer's choices;
- how the common readers query the result.

## Directory layout

The path given to `--output_file` becomes a directory holding the parts, numbered in the order they were written, and a manifest written last:

```
out.parquet/
├── part-00000.parquet
├── part-00001.parquet
├── part-00002.parquet
└── _vep_rs_manifest.json
```

A part closes once it holds 512 MiB of compressed data and, when the input order is [checked](#row-order-and-the-ordering-precondition), at every chromosome boundary. So under the check every part holds one chromosome and the parts follow the input's chromosome order. A chromosome is one run of records: the spellings and synonyms that resolve to one cache chromosome (`21` and `chr21`) share a part. A run with no rows writes one empty part.

`_vep_rs_manifest.json` lists the run that wrote the directory and its parts:

- `vep_rs_version`;
- `run_id`, the same value as every part's footer carries;
- the parts in order with their row counts;
- the total row count.

Its leading underscore keeps it out of a `*.parquet` glob and out of Spark's file listing.

The output is a directory, so `-o STDOUT` is refused (`ERROR: [output] --output_format parquet needs a file path: Parquet cannot be written to STDOUT`). `--fields` and `--compress_output` are refused with it too: the column set is fixed, and the parts carry their own compression. The text formats' header block has no counterpart here, so `--no_headers` changes nothing; the run's identity is in each part's [footer](#footer-metadata).

### `--parquet_partition chrom`

`--parquet_partition chrom` lays the parts out in one directory per chromosome as the input spells it, a part closing at every change of spelling:

```
out.parquet/
├── chrom=21/part-00000.parquet
├── chrom=22/part-00001.parquet
├── chrom=X/part-00002.parquet
└── _vep_rs_manifest.json
```

The spelling is escaped as Hive escapes a partition value (`chrom=HLA-A%2A01%3A01`). `chrom` is still a column of every part, so a reader need not parse the directory name. A reader that infers Hive partitions types `chrom=21` as a number while the column stores a string, so the [recipes](#reading-the-output) below switch the inference off. The flag takes `chrom` alone and applies to `--output_format parquet` alone (`ERROR: [output] --parquet_partition applies to --output_format parquet`).

## Schema

Six key columns, none nullable, then exactly the default `--tab` columns in `--tab` order, each nullable and typed as the next section states:

| Column         | Type   | Value                                                                                                        |
| -------------- | ------ | ------------------------------------------------------------------------------------------------------------ |
| `contig_index` | INT32  | the rank of the record's chromosome in the footer's `contig_order`                                           |
| `chrom`        | STRING | the chromosome as the input spells it                                                                        |
| `pos`, `end`   | INT64  | the record's span, 1-based and inclusive, the smaller coordinate first (an insertion's two flanking bases)   |
| `ref`, `alt`   | STRING | the reference allele and the row's alternate allele, `-` for an empty allele                                 |

Every row is one (allele, feature) pair, the rows of each input record in the row order of the `--tab` format. An allele overlapping no feature has its intergenic row. The tab columns are the ones the run's flags produce:

- `--symbol` adds `SYMBOL`;
- `--max_af` adds `MAX_AF`;
- a plugin adds its fields;
- a custom source adds its columns.

Each is a nullable column of the type below.

## The cell contract

Every tab column has a type with two functions, `parse(text)` from the `--tab` cell to the stored value and `to_tab_text(value)` back. For every cell the tab format writes, `to_tab_text(parse(t))` is `t`, except that a DOUBLE prints as its shortest round-trip decimal. `0.0950` is stored as 0.095 and prints `0.095`.

An INT32 `TSL` cell is stored as its integer. So a support level a GFF3 wrote with leading zeros (`007`) prints without them, and one past 2147483647 is NULL. The levels of an Ensembl gene set are one digit.

A reader that projects a row through `to_tab_text` column by column gets the `--tab` row.

| Columns                                                                                          | Type                                           | `--tab` cell                                                                                                                                                                                                                                       |
| ------------------------------------------------------------------------------------------------ | ---------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `STRAND`, `DISTANCE`, `ALLELE_NUM`, `TSL`, `HGVS_OFFSET`, `OverlapBP`, `MOTIF_POS`               | INT32                                          | the decimal integer                                                                                                                                                                                                                                |
| `MAX_AF`, `MOTIF_SCORE_CHANGE`                                                                   | DOUBLE                                         | the decimal number; it prints as the shortest decimal that reads back to the same value                                                                                                                                                            |
| `CANONICAL`                                                                                      | BOOLEAN                                        | `YES` is true; false prints `-`, which the tab format never writes                                                                                                                                                                                 |
| `PICK`, `GENCODE_PRIMARY`                                                                        | BOOLEAN                                        | `1` is true; false prints `-`                                                                                                                                                                                                                      |
| `HIGH_INF_POS`                                                                                   | BOOLEAN                                        | `Y` is true and `N` false                                                                                                                                                                                                                          |
| the 28 population frequency columns (`AF`, `AFR_AF`, `gnomADe_AF`, `gnomADg_NFE_AF`, ...)        | LIST\<DOUBLE\>                                 | the `,`-separated numbers of the cell; it prints them `,`-joined                                                                                                                                                                                   |
| `SOMATIC`, `PHENO`                                                                               | LIST\<BOOLEAN\>                                | the `,`-separated `0` and `1` flags                                                                                                                                                                                                                |
| `MANE`, `Consequence`, `CLIN_SIG`, `PUBMED`, `FLAGS`, `TRANSCRIPTION_FACTORS`                    | LIST\<STRING\>                                 | the cell split on `,`; it prints the elements `,`-joined                                                                                                                                                                                           |
| `EXON`, `INTRON`, `cDNA_position`, `CDS_position`, `Protein_position`                            | STRUCT\<start INT32, end INT32, total INT32\>  | `[span][/total]`: `span` is `start`, `start-end`, `start-?`, `?-end` or `-`, a `?` or `-` member NULL and a single position setting `start` and `end` alike; `total` is the `--total_length` suffix. It prints `start` when the two agree, `start-end` otherwise, `?` for a NULL member beside a set one, `-` for neither, then `/total` when set |
| every other column: the identifiers, HGVS notations, codons, `IMPACT`, `IND`, `ZYG`, plugin and custom fields | STRING                             | the cell as written, each invalid UTF-8 sequence replaced by U+FFFD                                                                                                                                                                                 |

A `-` cell is NULL, and NULL prints `-`, in every column but the allele columns `Allele`, `GIVEN_REF` and `USED_REF`. There `-` is an empty allele (a deletion's alternate allele, an insertion's reference allele) and is stored as the text `-`. A cell the tab format prints empty is stored as the empty string. A token outside its column's grammar, such as letters in an integer column, is stored as NULL.

`Existing_variation` is one STRING per row, the co-located identifiers `,`-joined as the tab format writes them (`rs1234,rs5678`). So an equality test finds a row only when its whole cell equals the value. A cell holding several identifiers needs `contains` or `LIKE`, which the Bloom filter below cannot serve.

## Reading the output

The recipes read both layouts, the flat directory and the `--parquet_partition chrom` one, with `chrom` as the stored string. Every example row set is the default run's; `SYMBOL` needs `--symbol`, `MAX_AF` needs `--max_af`, each a column once its flag is given.

### DuckDB

The `**` glob matches the parts in the directory and in its `chrom=` subdirectories; `hive_partitioning = false` keeps the stored `chrom`.

```sql
SELECT chrom, pos, ref, alt, Feature, Consequence, IMPACT
FROM read_parquet('out.parquet/**/*.parquet', hive_partitioning = false)
WHERE Feature = 'ENST00000366667' AND IMPACT = 'HIGH'
ORDER BY contig_index, pos;
```

A region is a range on `pos` within one chromosome, pruned by the row-group statistics. A list column is tested with `list_contains`, a frequency list with `list_max` or `list_min`:

```sql
SELECT chrom, pos, ref, alt, Consequence, MAX_AF
FROM read_parquet('out.parquet/**/*.parquet', hive_partitioning = false)
WHERE chrom = '21' AND pos BETWEEN 42000000 AND 43000000
  AND list_contains(Consequence, 'missense_variant')
  AND MAX_AF < 0.01;
```

The footer is read with `parquet_kv_metadata`, and `parquet_bloom_probe` shows which row groups a Bloom filter excludes for a value:

```sql
SELECT decode(key) AS key, decode(value) AS value
FROM parquet_kv_metadata('out.parquet/part-00000.parquet');

SELECT file_name, row_group_id, bloom_filter_excludes
FROM parquet_bloom_probe('out.parquet/**/*.parquet',
                         'Feature', 'ENST00000366667');
```

### Polars

```python
import polars as pl

rows = pl.scan_parquet("out.parquet/**/*.parquet", hive_partitioning=False)
hits = (
    rows.filter(
        pl.col("Consequence").list.contains("missense_variant")
        & (pl.col("IMPACT") == "MODERATE")
    )
    .select("chrom", "pos", "ref", "alt", "Feature", "Consequence")
    .sort("contig_index", "pos")
    .collect()
)
```

A frequency column is a list of floats: `pl.col("gnomADg_AF").list.max() < 0.01` keeps the rows whose every co-located variant is rare.

### Spark

`recursiveFileLookup` reads the parts of either layout and switches partition discovery off, so `chrom` is the stored string; Spark skips the underscore-prefixed manifest by itself.

```python
from pyspark.sql import functions as F

rows = spark.read.option("recursiveFileLookup", "true").parquet("out.parquet")
hits = (
    rows.filter(
        F.array_contains(rows.Consequence, "missense_variant")
        & (rows.Gene == "ENSG00000141510")
    )
    .select("chrom", "pos", "ref", "alt", "Feature", "Consequence")
    .orderBy("contig_index", "pos")
)
```

### Athena

A table over the directory names and types each column as the schema states:

| Schema type          | Athena type                           |
| -------------------- | ------------------------------------- |
| INT32                | `int`                                 |
| INT64                | `bigint`                              |
| STRING               | `string`                              |
| DOUBLE               | `double`                              |
| BOOLEAN              | `boolean`                             |
| a list               | `array<...>`                          |
| the positions struct | `struct<start:int,end:int,total:int>` |

`end` is a reserved word and is quoted; `LOCATION` is the bucket and prefix the directory was copied to, as an `s3` URL. The default run's columns:

```sql
CREATE EXTERNAL TABLE vep_output (
  contig_index int,
  chrom string,
  pos bigint,
  `end` bigint,
  ref string,
  alt string,
  Uploaded_variation string,
  Location string,
  Allele string,
  Gene string,
  Feature string,
  Feature_type string,
  Consequence array<string>,
  cDNA_position struct<start:int, `end`:int, total:int>,
  CDS_position struct<start:int, `end`:int, total:int>,
  Protein_position struct<start:int, `end`:int, total:int>,
  Amino_acids string,
  Codons string,
  Existing_variation string,
  IMPACT string,
  DISTANCE int,
  STRAND int,
  FLAGS array<string>
)
STORED AS PARQUET
LOCATION '<bucket URL>/out.parquet/';
```

Athena reads the parts under the location, the `chrom=` directories of a partitioned output included, and ignores the underscore-prefixed manifest. The table declares no partition column, since `chrom` is already a column of every part. A run with more flags adds its columns to the list in `--tab` order.

## Round-trip to `--tab`

`scripts/adapters/parquet_to_vep_tab.py` is the cell contract's inverse in Python, on `pyarrow`. It prints the `#` header line of the tab columns, then one line per row through `to_tab_text`. It reads one part, a directory of parts or a `--parquet_partition chrom` directory. The `##` meta lines of a live run are not stored in Parquet and are not reproduced, and the rows come out in file order. So compare against a `--tab` file as a multiset of lines.

```bash
scripts/adapters/parquet_to_vep_tab.py out.parquet [--output out.tab]
scripts/adapters/parquet_to_vep_tab.py --canonical out.parquet
scripts/adapters/parquet_to_vep_tab.py --canonical out.tab
```

`--canonical` prints the canonical form, the header then the data lines sorted with every DOUBLE cell as its shortest round-trip decimal and every INT32 cell as its integer. It takes a `--tab` output (plain or gzip) as well. Two outputs with equal canonical forms hold the same rows. `parquet_round_trips_to_the_tab_output` in `crates/vep-cli/tests/format_parity.rs` asserts that equality on every golden corpus, through the crate's own reader.

## Row order and the ordering precondition

Every row group is sorted by (`contig_index`, `pos`). `contig_index` ranks the chromosomes in the footer's `contig_order`. That list holds, in this order:

1. the chromosome list of the cache's own metadata, in its karyotype's order (for a human cache the autosomes by number, then X, Y and MT);
2. the cache's other chromosomes, by the same rule;
3. under `--fasta`, the reference's other sequences in its index's order.

A chromosome in none of those lists shares the rank after the last listed one, and its rows are told apart by their `chrom` spelling. `ORDER BY contig_index, pos` is genome order whatever the input's chromosome spellings.

The writer checks every row's record against the rows before it. A chromosome's records stay together, in any order of chromosomes. Within a chromosome the position the input wrote (a VCF record's POS, an Ensembl-format line's start) never decreases.

The first record that breaks this ends the run with exit 1. The run leaves the output path [as it was](#overwriting-an-output): no directory at a new path, the earlier output where one stood. The error names the record's line and the two positions: `ERROR: [parquet] input is not sorted by chromosome and position at line 2 (21:50 after 21:100); sort the input or pass --no_check_variants_order`.

A refused run's warnings file holds the warnings of every input batch (`--buffer_size` records) read through the one that holds the refused record, in input order and at any worker count. `--pipeline_mode batched` adds those of the later batches it had read.

The check reads the written position rather than the stored `pos`, which is the annotated span's start. A VCF deletion's start is one past its POS once the anchor base is trimmed. A file sorted by POS therefore steps back in `pos` between a deletion and a substitution at the same POS, and `--minimal` trims further.

Those steps are sorted away before the rows are encoded. So a sorted input yields parts wholly in (`contig_index`, `pos`) order, and the footer's `sort` key says `contig_index,pos`. A part whose row groups are each in order but not in order with one another says `contig_index,pos within each row group`.

### How the check differs from Ensembl VEP's

This check is a documented subset of Ensembl VEP's. Ensembl VEP checks the order of every run whatever the output format. It counts a record whose start lies more than 100 bp before the record before it on the same chromosome. It ends the run once that count passes 100 (`InputBuffer.pm`, `Constants.pm`). It checks no chromosome order, and exempts identifier, HGVS, SPDI and region input (`Config.pm` lines 640-645).

vep-rs checks the Parquet output alone, strictly per record, and exempts the same four input formats, whose rows are written unchecked. Every text format annotates an input in any order and ignores the flag.

### Skipping the check

`--no_check_variants_order` skips the check. Then:

- the rows are encoded in input order;
- no part closes at a chromosome boundary (a `--parquet_partition chrom` output still closes one at every change of spelling);
- no row group declares a sort;
- the footer's `sort` key is `none`.

The run writes one warning to standard error and the warnings file: `WARNING: input order not checked: Parquet row groups are not sorted and parts are not rolled at contig boundaries`.

## Writer options

Every option is chosen for the reader:

| Option          | Value                                                                                                                                                      |
| --------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Page format     | Parquet 2 data pages                                                                                                                                       |
| Compression     | ZSTD, level 9, on every column chunk                                                                                                                       |
| Row groups      | at most 128 MiB of encoded data or 1,048,576 rows                                                                                                          |
| Parts           | closed at 512 MiB compressed and, under the order check, at every chromosome boundary                                                                     |
| Sort            | `contig_index` then `pos`, ascending, declared in every row group's `sorting_columns` under the order check                                                |
| Statistics      | page-level statistics, the column index and the offset index on every column, so a reader skips pages by `pos`, `chrom` or any typed column's min and max |
| Bloom filters   | `Uploaded_variation`, `Existing_variation`, `Feature`, `Gene` and `SYMBOL`, at a 1 % false-positive rate, and no other column                              |
| Footer metadata | the keys of the next section, and the `ARROW:schema` key Arrow readers use to restore the list and struct field names                                      |

A Bloom filter answers an equality test on the whole cell. `WHERE Feature = 'ENST00000366667'` skips every row group whose filter rules the value out. `parquet_file_layout_follows_the_declared_writer_properties` in `crates/vep-cli/tests/format_parity.rs` asserts that a part's filter admits every `Feature` of the part and rules out at least 190 of 200 identifiers absent from it. A range, prefix or `contains` predicate is served by the statistics alone, and an identifier inside a multi-valued `Existing_variation` cell is not found by equality.

## Footer metadata

Each part's footer carries the run's identity in its key-value metadata, the same values in every part of a run:

| Key                                     | Value                                                                                                                                            |
| --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `vep_rs_version`, `vep_api_version`     | the vep-rs release and the Ensembl VEP API version it follows                                                                                    |
| `cache`, `assembly`, `source_<name>`    | the cache directory as the output headers print it, the assembly, one key per cache source with its version (`source_gencode`, `source_dbSNP`)   |
| `command_line`                          | the command as run                                                                                                                               |
| `coordinate_system`                     | `1-based inclusive`                                                                                                                              |
| `contig_order`                          | the `,`-joined chromosome list `contig_index` ranks                                                                                              |
| `layout`                                | `flat`: one row per (allele, feature) pair                                                                                                       |
| `sort`                                  | `contig_index,pos`, `contig_index,pos within each row group`, or `none` when the order was not checked                                           |
| `run_id`                                | an identifier unique to the run, equal in every part and in the manifest                                                                         |

## Overwriting an output

A run never writes into the output path. It writes its parts and manifest into a temporary directory beside it, `.<name>.vep-rs-tmp-<run_id>`. It renames that directory to the output path once the run has succeeded. A directory at the output path is therefore a completed run.

A run that fails or is refused removes its temporary directory, naming it in its error line when the removal fails. It leaves the output path as it was. A path with no final name to place that directory beside, `.` or `..`, is refused as an output that cannot be written.

An interrupted run can leave the temporary directory beside the output, where no later run reads or removes it. Its leading dot keeps it out of a shell glob over the output's directory and out of a `pyarrow` dataset opened there. DuckDB's glob matches it.

The temporary directory is created beside the output and renamed over it. So the output's parent directory must be writable and the output path cannot be a mount point. Under `--force_overwrite`, a run into a mount point ends at its last step with `ERROR: [output] Could not write to output file <path>`, the path as it was.

### Warnings, skipped-variants and statistics files

Name the warnings file, the skipped-variants file and the statistics report outside the output directory. A new path has no such directory until the run has succeeded. A file written into an earlier output is an entry vep-rs did not write, which refuses that replacement at its last step and every later one. An output path written with a trailing separator (`out.parquet/`) places their default names inside it, so write the path without one.

### Replacing an earlier output

An existing output refuses the run, as every output file does. `--force_overwrite` replaces a previous run's output only when every entry of the directory is one vep-rs writes. Such an entry is a part whose footer carries `vep_rs_version`, the manifest, or a `chrom=` directory of such parts. Anything else refuses the run before it starts, with exit 1 and every file still in place:

```
ERROR: [parquet] output directory contains files vep-rs did not write
```

The earlier output stays in place until the new one is complete, so a failed or refused replacement leaves it untouched. The statistics report beside it is a file of its own, emptied as a replacement starts. The warnings file is rewritten once the replacement draws a warning.

Once the run has succeeded it:

1. checks the directory again, since an entry can arrive while it runs;
2. renames the earlier output aside;
3. renames the new one in;
4. removes the earlier one.

The directory at the output path is then a new one with the earlier directory's permission bits. Its owner, group and access-control list are those a directory made beside the output gets.

If the earlier output cannot be removed at that point, the run ends with exit 1 and names the directory beside the output that holds what is left of it. The new output is in place and complete.

A regular file at the output path is replaced only when it is a part vep-rs wrote. Any other file is refused with `ERROR: [parquet] output path <path> is a file vep-rs did not write`. A symbolic link at the output path, or inside the earlier output, refuses the run by name: `ERROR: [parquet] output path <path> is a symbolic link`. What it points at is neither read nor changed.

## Removed flags and 0.5.0 output

`--parquet_shape`, `--keep_intermediate` and `--parquet_row_group_size` are refused as unknown arguments (exit 2, `error: unexpected argument '--parquet_shape' found`). The writer has one row shape, writes no intermediate file, and closes a row group by its encoded size.

A Parquet output written by 0.5.0 differs from the one this page describes, so a query written against it changes meaning here. 0.5.0:

- wrote one `chrom=<c>/data_0.parquet` file per chromosome through the `duckdb` CLI;
- had five key columns (no `contig_index`) and a nested row shape by default;
- stored every tab column but `DISTANCE` and `STRAND` as a VARCHAR, so `WHERE gnomADg_AF < '0.01'` compared strings;
- sorted the rows by the `chrom` string;
- put a Bloom filter on every dictionary-encoded column;
- wrote no `coordinate_system`, `contig_order`, `layout`, `sort` or `run_id` key in the footer.

Both outputs project to the same `--tab` rows.

## Where next

- [VEP default and tab formats](vep-and-tab.md): the columns and the row order of the `--tab` format.
- [Engine identity](../output-formats.md#engine-identity): where each format names the vep-rs release.
- [CLI reference](../cli-reference.md): every flag, output format and plugin.
