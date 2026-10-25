# What vep-rs does not do

This page is for a reader who brings an Ensembl VEP command line to vep-rs and wants to know, in one place, which Ensembl VEP features vep-rs does not have and what a run does when asked for one. Every flag Ensembl VEP accepts is honoured, refused with an error that names it, or accepted without effect with its help text saying so. The [CLI reference](../cli-reference.md) and its [Compatibility with Ensembl VEP](../cli/compatibility.md) page carry the per-flag detail; each family below gets its flags, what a run does with them, and the release that adds the feature, if any.

## RefSeq and merged caches

`--refseq` and `--merged` select the cache, the species directory suffixed `_refseq` or `_merged`, as Ensembl VEP does. Under `--json_cache`, `--merged` changes nothing and `--refseq` only puts RefSeq transcripts first when an HGVS notation names a gene. `--exclude_predicted` adds Ensembl VEP's own filter, `not stable_id match ^X._`, which drops RefSeq's predicted `XM_` and `XR_` models. A converted RefSeq or merged cache annotates every transcript it carries, including those Ensembl VEP drops without `--all_refseq`, each as an Ensembl transcript. The `REFSEQ_MATCH`, `BAM_EDIT`, `REFSEQ_OFFSET` and `SOURCE` columns are not written, a BAM-edited transcript's edited sequence is not used, and the `ensembl` and `refseq` pick criteria tie; `--all_refseq` and `--bam` are refused as unknown arguments (exit 2).

Planned for the next release, 0.6.0.

## Regulatory and motif features

`--regulatory`, alone or through `--everything`, is accepted and warns that regulatory-feature annotation is not implemented: the motif columns stay empty, no `RegulatoryFeature` or `MotifFeature` row is written, and the six regulatory consequence terms are never assigned ([the divergence page](../divergences/regulatory-rows.md) states it). `--cell_type` is refused as an unknown argument. `--regulatory_gff` and `--extended_promoters` are each refused with an error naming the flag; the error for `--regulatory_gff` begins `ERROR: [cli] --regulatory_gff is not supported: it would annotate regulatory features and motifs from a GFF3 file in place of the cache`.

Planned for the next release, 0.6.0.

## Custom annotation sources

`--custom <file>,<name>,<format>,...` is accepted and the source is not read, so the run proceeds and no custom field appears; the flag's help text says so. `--gff`, `--gtf`, `--bigwig`, `--phyloP` and `--phastCons`, which Ensembl VEP reads as `--custom` sources, and `--custom_multi_allelic`, `--ucsc_assembly` and `--ucsc_data_root` are refused as unknown arguments (`error: unexpected argument '--gff' found`, exit 2). `--custom_suppress_filter` is refused with an error that begins `ERROR: [cli] --custom_suppress_filter is not supported: it would keep a --custom VCF source's FILTER column out of the output`.

Planned for the next release, 0.6.0.

## Plugins

`--plugin <Name>,<params>`, `--dir_plugins` and `--safe` work as in Ensembl VEP. A name resolves to one of the built-in plugins listed on [Plugins](../plugins.md), else to a shared library under `--dir_plugins`; a name that resolves to neither is skipped with a warning, or refused under `--safe` (an error that begins `ERROR: [plugin] Failed to load plugin <Name>`). The rest of Ensembl VEP's plugin catalogue has no vep-rs counterpart: those plugins are Perl modules, which vep-rs cannot load, so `--plugin` with one of their names is that warning or that error. The built-in plugins ship as [Plugins](../plugins.md) documents them; this release publishes no measurement of their agreement with the Ensembl VEP plugin modules, and 0.6.0 measures it.

Planned for the next release, 0.6.0.

## The statistics report

`--no_stats`, `--stats_file` (`--sf`) and `--stats_text` are accepted and change nothing: no statistics file is written, with or without them. `--stats_html` is refused as an unknown argument (`error: unexpected argument '--stats_html' found`, exit 2). Ensembl VEP writes the report, a summary of the run's records and consequences, as HTML beside the output by default and as text under `--stats_text`.

Planned for the next release, 0.6.0.

## Parquet columns

`--output_format parquet` writes a fixed column set, five key columns and then the default `--tab` columns ([Parquet](../output/parquet.md)); `--fields` beside it is refused with an error: `ERROR: [output] --fields applies to --tab and --vcf output; the Parquet column set is fixed`. Parquet is vep-rs's own format, with no Ensembl VEP counterpart; a column set chosen with `--fields` is what it lacks.

Planned for the next release, 0.6.0.

## The companion tools

Ensembl VEP ships `filter_vep`, which filters an annotated output by field expressions, `haplo`, which computes transcript haplotypes from phased genotypes, and `variant_recoder`, which translates a variant between its identifiers and notations. vep-rs's own programs are `vep`, `vep-cache-builder` and `vep-cache-converter`, so the release archives and the container image carry no program by those names.

Planned for the next release, 0.6.0.

## Progress and debug output

`--no_progress` and `--debug` are refused as unknown arguments (`error: unexpected argument '--no_progress' found`, exit 2). Ensembl VEP's `--no_progress` suppresses its progress bars and `--debug` prints plugin-loading detail; neither changes an annotation.

Planned for the next release, 0.6.0.

## Ensembl database mode

`--database`, `--host`, `--port`, `--user`, `--password`, `--registry`, `--db_version`, `--genomes` and `--is_multispecies` are each refused with an error naming the flag; the error for `--database` begins `ERROR: [cli] --database is not supported: vep-rs reads a local cache only and never connects to an Ensembl database`. `--lrg` and `--ambiguous_hgvs`, which need the database, are refused with them. Database mode is refused for good and is not planned: vep-rs reads a local cache only.

## 1000 Genomes pilot coordinates

`--gp` is refused with an error naming it, which begins `ERROR: [cli] --gp is not supported: it reads a GRCh37 chr:pos from the INFO/GP field that some 1000 Genomes VCFs carry`. In Ensembl VEP the flag reads a record's chromosome and position from the `GP` field of its INFO column, the convention of the 1000 Genomes pilot VCFs, and skips a record without that field (`Config.pm` line 120, `Parser/VCF.pm` lines 271-285). Not planned.

## Very long lines

An input, cache variation or FASTA index line longer than 64 MiB is refused under the `input`, `cache` or `fasta` component (exit 1), with `longer than the 67108864-byte limit for one line` in the message, where Ensembl VEP reads such a line whole. Such a line drops a plugin at load (the `Failed to instantiate plugin` warning, or the `plugin` error under `--safe`); met by a query, it is a `plugin execution error` warning, `--safe` or not, and the batch loses that plugin's values and every later plugin's. `vep-cache-builder` refuses a GFF3, GTF, protein FASTA or core-dump line longer than 16 MiB under the `cache-builder` component, with `16777216-byte` in the message, and `vep-cache-converter` a plugin data line over 64 MiB.

## Bioconda

There is no Bioconda package: vep-rs is installed from a pre-built binary on the [Releases](https://github.com/natera-open-source/vep-rs/releases) page, a build from source or the container image, as [Installation](installation.md) describes. The Bioconda recipe is maintained outside this repository.

Planned for the next release, 0.6.0.

To ask for a feature, or to report one this page leaves out, open an issue on the [issue tracker](https://github.com/natera-open-source/vep-rs/issues).
