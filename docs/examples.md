# Examples

Every flag is documented in the [CLI reference](cli-reference.md); these are the common shapes of a run.

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

```bash
cadd_snv=/path/to/cadd/whole_genome_SNVs.tsv.gz
cadd_indels=/path/to/cadd/gnomad.genomes.r4.0.indel.tsv.gz
vep -i input.vcf -o output.txt \
    --offline --cache --json_cache /path/to/json_cache \
    --assembly GRCh38 \
    --plugin CADD,snv=$cadd_snv,indels=$cadd_indels \
    --plugin REVEL,file=/path/to/revel/revel_scores.tsv.gz \
    --plugin SpliceAI,snv=/path/to/spliceai/spliceai_scores.raw.snv.vcf.gz \
    --everything --force
```

`--assembly` is not optional here. CADD and dbscSNV use it to reject a data file
built from the other genome; REVEL uses it to select which of its two coordinate
columns to read and rejects a copy whose tabix index was built on the other
column. Without it those guards are inert and REVEL falls back to the GRCh37
column with a warning, so a GRCh38-indexed file annotates nothing. SpliceAI also
accepts `indel=` and `cutoff=`; only the openly published GRCh38 SNV file is
validated, and the indel files, distributed through Illumina BaseSpace, are
unsupported data rather than a rejected parameter (see [Plugins](plugins.md)).

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
    --output_format parquet --parquet_shape flat
```

The `duckdb` CLI must be on `PATH`; the output is a directory, read with
`read_parquet('out.parquet/**/*.parquet', hive_partitioning=false)` as
the [Parquet](output/parquet.md#directory-layout) page describes.
