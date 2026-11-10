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
convert-plugin-data`, with its `.vpdi` index and `.vpd.provenance.json`
manifest beside it. The source file must remain present.

A store named `whole_genome_SNVs.tsv.vpd` beside
`whole_genome_SNVs.tsv.gz` is selected only after its manifest verifies.
Verification hashes the complete compressed source and both binary files with
SHA256 and checks the chromosome, start, end and zero-based conversion settings
against the plugin's settings. Reference-dependent end settings that the
converter cannot express also fall back to tabix.

A stale source, swapped data or index, mismatched settings, missing legacy
manifest or corrupt manifest produces a warning and falls back to the requested
tabix source. That fallback requires its tabix index. Regenerate older stores
with `vep-cache-converter convert-plugin-data` to enable automatic selection.
Normal startup logging names the selected backend and paths; a verified binary
selection also logs the source, data and index SHA256 identities.

Verification reads all three files at each store open, including the full
compressed source. It uses bounded hashing buffers but can add substantial
startup I/O for genome-wide datasets; timestamps and lengths never replace
content verification. Hashing touches every mapped data page, and the index
remains loaded in memory. The binary layout and query algorithms are unchanged.

Conversion hashes the compressed stream it parses and rehashes the source before
publication. It rejects a changed source or settings. Data and index are staged
beside the destination, and the completion manifest is published last, after
writing and any requested concordance check succeed. `--no-verify` skips
concordance only; content identity and publication checks still run.

A failed conversion before publication keeps an existing completed store.
A failure during publication leaves no manifest that certifies a partial pair,
so subsequent annotation uses tabix. Conversion needs temporary disk space for
the new pair alongside any existing pair. A destination lock prevents concurrent
converters from mixing their files. After a killed converter, remove its
`.vpd.conversion-lock` directory only after confirming no converter still uses it.
Keep sources and published files immutable during annotation; replace stores
through the converter rather than editing mapped files in place.

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
