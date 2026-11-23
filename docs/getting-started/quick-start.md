# Quick start

Start with a freshly built `vep`: run the bundled example, install a transcript cache, then annotate your own VCF, with plugins where you need them.

## 1. Try the bundled example

This needs no download. It annotates `examples/homo_sapiens_GRCh37.vcf`, 175 records cut from the GRCh37 golden corpus, against the pruned cache that ships beside that corpus in `tests/golden/`. The repository's `examples/README.md` describes both example files:

```bash
target/release/vep \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --assembly GRCh37 \
  -i examples/homo_sapiens_GRCh37.vcf -o output.txt
```

## 2. Get a transcript cache

For real data you need a JSON transcript cache. One command installs a release's cache and its indexed genome FASTA under `$HOME/.vep` in the layout `vep` resolves. The command shows GRCh38. `GRCh37` in its place installs that assembly without the 1000 Genomes frequencies, PubMed identifiers and variant synonyms, whose files the GRCh37 tree lacks, and takes the SIFT and PolyPhen matrices with `--include-predictions --mysql-port 3337`. [Cache setup](../cache-setup.md) explains each input:

```bash
# Ensembl release 116, GRCh38
cargo build --release -p vep-cache-builder
target/release/vep-cache-builder --auto homo_sapiens GRCh38 --release 116 \
  --dir_cache $HOME/.vep
```

The cache lands at `$HOME/.vep/homo_sapiens/116_GRCh38` with the release's known variants, their 1000 Genomes frequencies, PubMed identifiers, synonyms and failed flag, SIFT and PolyPhen matrices, and regulatory and motif features with their cell-type activity and binding matrices. It is about 32 GB of downloads checked against Ensembl's `CHECKSUMS` files, with no `wget`, `gunzip` or `samtools` step; `--chromosomes` narrows the variation VCFs and not the table dumps. A species with one variation VCF also fetches the `variation_feature` table for a sequence of its cache that file holds no record on. The installation's last lines name each part it left out and the gnomAD frequencies, which no release file carries ([Cache setup](../cache-setup.md#build-a-cache)).

The inputs can also be named one flag at a time. A GTF alone, with the genome, builds the cache Ensembl VEP's `--gtf` would annotate from instead ([A GTF or GFF3 as the transcript model](../cache-from-gtf-or-gff3.md)). Or convert an existing Ensembl VEP cache with `tools/vep-storable-to-json` ([scripts/README.md](../../scripts/README.md)).

## 3. Annotate

Run the same command you run with Ensembl VEP. Name the cache root with `--dir_cache`, as Ensembl VEP does. Or name the cache directory itself, the one holding `info.json` and `transcripts/`, with `--json_cache`:

```bash
target/release/vep \
  --dir_cache $HOME/.vep --cache_version 116 --assembly GRCh38 \
  -i input.vcf -o output.txt
```

## Usage

```bash
# With annotation plugins: REVEL on data staged by
# scripts/data/setup_plugin_data.sh, CADD on your own file
plugin_data=$HOME/.vep/plugin_data/grch38/full
cadd_snv=/path/to/your/snv.tsv.gz
target/release/vep \
  --json_cache $HOME/.vep/homo_sapiens/116_GRCh38 \
  --assembly GRCh38 \
  --plugin CADD,snv=$cadd_snv \
  --plugin REVEL,file=$plugin_data/revel/revel.tsv.gz \
  -i input.vcf -o output.txt --force_overwrite

# JSON output with HGVS notation
target/release/vep \
  --json_cache $HOME/.vep/homo_sapiens/116_GRCh38 \
  --assembly GRCh38 --json \
  --hgvs --everything \
  -i input.vcf -o output.json
```

**Plugins.** `--plugin <Name>,<params>` resolves a built-in plugin first, and otherwise loads a C-ABI shared library from `--dir_plugins`. Support is per reference genome: `--assembly` selects the right column and refuses a data file whose own assembly contradicts it. The support table, the data files and how to write a plugin are in [Plugins](../plugins.md).

**HGVS.** The FASTA `--auto` places in the cache directory is found without a flag. `--fasta` takes a `samtools faidx`-indexed genome, plain or compressed with `bgzip` (Ensembl's `dna_index` files), and gives `--hgvs` the reference context for the 3' shifting of indels. `--hgvs` without a FASTA is refused, as Ensembl VEP refuses it offline.

**Output formats.** The first four match Ensembl VEP's layouts, and the fifth is vep-rs's own:

- VEP default: the tab-delimited text format with the `Extra` column.
- `--vcf`: annotations in the `CSQ` INFO field.
- `--json`: one object per variant.
- `--tab`.
- `--output_format parquet`: a directory of typed Parquet parts written natively, with no external tool.

Field-level specs are in [Output formats](../output-formats.md).

**All flags.** The [CLI reference](../cli-reference.md) is the complete flag reference, including the flags accepted for Ensembl VEP compatibility but not implemented.
