# CLI reference

The `vep` binary accepts the Ensembl VEP command-line interface. Most flags have the same name and behaviour. A flag a group page lists under "Accepted for compatibility" is parsed and has no effect on the annotation, so an existing Ensembl VEP command line is accepted; the cache flags among them only name the cache path the output header prints. The annotation cache is always `--json_cache`.

```bash
# Show all flags
vep --help

# Short form
vep -i input.vcf -o output.txt --offline --json_cache /path/to/json_cache
```

## Option groups

Every option has an entry on one of these pages, with a summary table per topic at the top of the page. An entry's anchor is `opt_` and the flag name, for example [`#opt_max_sv_size`](cli/cache-species-runtime.md#opt_max_sv_size) on the Cache, species and runtime page.

| Group | Covers |
| --- | --- |
| [Input and output](cli/input-output.md) | The input file and format, the output file and format, the header lines and the Parquet shape. |
| [Cache, species and runtime](cli/cache-species-runtime.md) | The JSON cache and FASTA, species and assembly, the batch size, thread counts and distances. |
| [Output control](cli/output-control.md) | The field flags, `--everything` among them, that choose which annotation columns appear. |
| [Allele frequencies and co-located variants](cli/allele-frequencies.md) | Population frequencies, PubMed identifiers and gene-phenotype associations. |
| [Filtering, picking and transcript sets](cli/filtering-picking.md) | The pick and flag-pick families, the allele columns and the transcript-set flags. |
| [Plugins](cli/plugins.md) | Loading built-in and dynamic-library plugins. |
| [Statistics, verbosity and miscellaneous](cli/statistics-verbosity-misc.md) | The statistics flags, logging, `--dont_skip`, `--allow_non_variant` and `--mirna`. |
| [Compatibility with Ensembl VEP](cli/compatibility.md) | What vep-rs does not implement: refused flags, terms never assigned, input formats not read, database mode. |

## Flags with a documented divergence

| Flag | Note |
| --- | --- |
| `--max_sv_size` | Ensembl VEP 116 drops an oversize record before annotation and writes nothing for it in any format; vep-rs annotates it: rows in the default, `--tab` and `parquet` outputs, its VCF line without consequences, no JSON object. See [class 5 on the Defect C page](divergences/defect-c.md#5-a-record-ensembl-vep-drops-before-annotation) of [the intended divergences](intended-divergences.md) and the [entry](cli/cache-species-runtime.md#opt_max_sv_size). |
| `--everything` | Also switches on `--allele_number` and `--hgvsg`, which Ensembl VEP's `--everything` does not, so the `ALLELE_NUM` and `HGVSg` header keys differ there ([entry](cli/output-control.md#opt_everything)). |
| `--uniprot` | The `UNIPARC` and `UNIPROT_ISOFORM` columns are added but never computed ([entry](cli/output-control.md#opt_uniprot)). |
| `--domains` | The `DOMAINS` column is added but never populated ([entry](cli/output-control.md#opt_domains)). |
| `--minimal` | Per-allele minimisation of multi-allelic records is not implemented ([entry](cli/filtering-picking.md#opt_minimal)). |
| `--pick_order` | The pick order is fixed to Ensembl VEP's default and the final `length` tie-break is not applied ([entry](cli/filtering-picking.md#opt_pick_order)). |
| `--dont_skip` | Aborts the run on a variant that fails parsing; Ensembl VEP keeps the variant ([entry](cli/statistics-verbosity-misc.md#opt_dont_skip)). |

## Examples

### Basic offline annotation

```bash
vep -i input.vcf -o output.txt \
    --offline --cache \
    --json_cache /path/to/json_cache \
    --species homo_sapiens --assembly GRCh37
```

### VCF output with common fields

```bash
vep -i input.vcf -o annotated.vcf \
    --vcf --offline --cache \
    --json_cache /path/to/json_cache \
    --everything
```

### JSON output

```bash
vep -i input.vcf -o output.json \
    --json --offline --cache \
    --json_cache /path/to/json_cache \
    --everything
```

### Clinical plugins

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

### Parallel annotation

```bash
vep -i large_input.vcf -o output.txt \
    --offline --cache --json_cache /path/to/json_cache \
    --fork 8 --buffer_size 10000
```

### Filtered output (one consequence per allele per gene)

```bash
vep -i input.vcf -o output.txt \
    --offline --cache --json_cache /path/to/json_cache \
    --pick_allele_gene --canonical \
    --everything
```
