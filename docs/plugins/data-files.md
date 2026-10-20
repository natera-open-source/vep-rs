# Preparing data files {#data-file-preparation}

This page states the file format every tabix plugin expects, the `.vpd` store
that replaces tabix, and the script that stages the files per assembly and tier.
The download and indexing steps for one plugin are on its own page:
[annotation-store plugins](annotation-store.md), [LoFTEE](loftee.md),
[flat-file plugins](flat-file.md) and [SpliceAI](spliceai.md).

## Tabix indexes and `.vpd` stores

All tabix plugins require bgzipped, tabix-indexed data files (`.tsv.gz` +
`.tsv.gz.tbi` or `.vcf.gz` + `.vcf.gz.tbi`). The one exception is a `.vpd`
binary store built from the data file by `vep-cache-converter
convert-plugin-data`, with its `.vpdi` index beside it. When a store named after
the data file with `.vpd` in place of `.gz` sits beside it
(`whole_genome_SNVs.tsv.vpd` beside `whole_genome_SNVs.tsv.gz`), the plugin
reads its records from the store instead of through tabix.

## The `setup_plugin_data.sh` script

Stage plugin data files with `scripts/data/setup_plugin_data.sh`. It does not
download from the upstream databases itself. `--source upstream`, the default,
prints, per plugin, the upstream URL and the preparation each file needs (header
fixes, per-assembly tabix indexes, decompression) for you to run. `--source s3`
syncs an already-prepared copy from an S3 bucket of your own, named by
`--bucket` or `$VEP_BENCHMARK_BUCKET`, that you seeded from those steps:

```bash
# print sources for chr21 subsets
scripts/data/setup_plugin_data.sh --assembly GRCh38 --tier quick \
    --source upstream

# sync genome-wide files from your bucket
scripts/data/setup_plugin_data.sh --assembly GRCh38 --tier full \
    --source s3 --bucket my-bucket
```

The genome-wide `full` tier is roughly 135 GB for GRCh38, 87 GB of it CADD.

`--assembly` is required: it selects which files are staged and which plugins
apply, and an unsupported combination is skipped with its reason rather than
substituted from the wrong build. Files land under
`~/.vep/plugin_data/<assembly>/<tier>/` by default, in per-plugin
subdirectories (`cadd/whole_genome_SNVs.tsv.gz`, `revel/revel.tsv.gz`, ...).
That is the layout `scripts/data/plugin_data_layout.sh` defines for the
data-preparation and concordance scripts alike.
