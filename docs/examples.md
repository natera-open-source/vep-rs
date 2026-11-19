# Examples

Copy a command for a common shape of a vep-rs run: offline annotation, VCF, JSON or Parquet output, clinical plugins, parallel annotation or filtered output.

In each command `/path/to/json_cache` stands for the cache directory, the one holding `info.json` and `transcripts/`. Step 2 of the [Quick start](getting-started/quick-start.md#2-get-a-transcript-cache) installs one. `input.vcf` and `large_input.vcf` stand for your own VCF.

## Basic offline annotation

```bash
vep -i input.vcf -o output.txt \
    --offline --cache \
    --json_cache /path/to/json_cache \
    --species homo_sapiens --assembly GRCh37
```

## VCF output with common fields

```bash
vep -i input.vcf -o annotated.vcf \
    --vcf --offline --cache \
    --json_cache /path/to/json_cache \
    --everything
```

## JSON output

```bash
vep -i input.vcf -o output.json \
    --json --offline --cache \
    --json_cache /path/to/json_cache \
    --everything
```

## Clinical plugins

The plugin data files come from each plugin's provider: the support table on [Plugins](plugins.md) names every file and its provider.

```bash
cadd_snv=/path/to/your/snv.tsv.gz
cadd_indels=/path/to/your/indels.tsv.gz
vep -i input.vcf -o output.txt \
    --offline --cache --json_cache /path/to/json_cache \
    --assembly GRCh38 \
    --plugin CADD,snv=$cadd_snv,indels=$cadd_indels \
    --plugin REVEL,file=/path/to/revel/revel_scores.tsv.gz \
    --plugin SpliceAI,snv=/path/to/your/spliceai_snv.vcf.gz \
    --everything --force
```

`--assembly` names the assembly the run is checked against. CADD and dbscSNV
use it to reject a data file built from the other genome; REVEL uses it to select
which of its two coordinate columns to read and rejects a copy whose tabix index
was built on the other column. Without the flag the cache's own assembly serves
those guards. Ensembl VEP 116.2 builds its plugins before it reads the cache, so
there REVEL, and dbscSNV without its own assembly parameter, fail to load without
`--assembly` and the run goes on without them.

SpliceAI also accepts `indel=`, `cutoff=` and `split_output=`. Its Illumina
files, available through BaseSpace (`SpliceAI.pm` lines 101-104), are for
academic use only under Illumina's terms, so the comparison with Ensembl VEP
runs on synthetic VCF files of the layout the module reads (`SpliceAI.pm`
lines 326-329; see [Plugins](plugins.md)).

## Parallel annotation

```bash
vep -i large_input.vcf -o output.txt \
    --offline --cache --json_cache /path/to/json_cache \
    --fork 8 --buffer_size 10000
```

## Filtered output (one consequence per allele per gene)

```bash
vep -i input.vcf -o output.txt \
    --offline --cache --json_cache /path/to/json_cache \
    --pick_allele_gene --canonical \
    --everything
```

## Parquet output

```bash
vep -i input.vcf -o out.parquet \
    --offline --cache --json_cache /path/to/json_cache \
    --output_format parquet
```

The output is a directory of typed Parquet parts, one row per (allele, feature)
pair, sorted by chromosome and position. The input must be sorted the same way,
or `--no_check_variants_order` writes it unchecked. Read it with
`read_parquet('out.parquet/**/*.parquet', hive_partitioning = false)` in DuckDB,
or as the [Parquet](output/parquet.md#reading-the-output) page describes for
Polars, Spark and Athena.

## Where next

- [CLI reference](cli-reference.md): every flag, documented.
- [Plugins](plugins.md): the support table and how a plugin is enabled.
- [Parquet](output/parquet.md): the columns of the Parquet output and how to read it.
- [Quick start](getting-started/quick-start.md): a first run, from the bundled example to your own VCF.
