# Quick start

This page takes a freshly built `vep` from its first run on the bundled corpus to annotating your own VCF with a transcript cache and plugins.

## 1. Try the bundled corpus

This needs no download. It runs the bundled 1,008-record excerpt of the benchmark datasets against the pruned GRCh37 cache that ships in `tests/golden/`:

```bash
target/release/vep \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --assembly GRCh37 \
  -i tests/golden/116/GRCh37/variants.vcf -o output.txt
```

## 2. Get a transcript cache

For real data you need a JSON transcript cache. Build one natively with `vep-cache-builder` (GRCh38 shown; [Cache setup](../cache-setup.md) explains each input and gives the GRCh37 command):

```bash
# Ensembl release 116, GRCh38
wget https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/\
Homo_sapiens.GRCh38.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/\
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

cargo build --release -p vep-cache-builder
target/release/vep-cache-builder \
  --species homo_sapiens --assembly GRCh38 --release 116 \
  --gtf Homo_sapiens.GRCh38.116.gtf.gz \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/json_cache/homo_sapiens/116_GRCh38
```

Or convert an existing Ensembl VEP cache with `scripts/data/storable_to_json.pl` ([scripts/README.md](../../scripts/README.md)).

## 3. Annotate

Run the same VEP command you run, with `--json_cache` naming the cache directory. That is the directory holding `info.json` and `transcripts/`, so `116_GRCh38` itself, not its parent:

```bash
target/release/vep \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/116_GRCh38 \
  --assembly GRCh38 \
  -i input.vcf -o output.txt
```

## Usage

```bash
# With annotation plugins, on data staged by scripts/data/setup_plugin_data.sh
plugin_data=$HOME/.vep/plugin_data/grch38/full
target/release/vep \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/116_GRCh38 \
  --assembly GRCh38 \
  --plugin CADD,snv=$plugin_data/cadd/whole_genome_SNVs.tsv.gz \
  --plugin REVEL,file=$plugin_data/revel/revel.tsv.gz \
  -i input.vcf -o output.txt --force_overwrite

# JSON output with HGVS notation
target/release/vep \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/116_GRCh38 \
  --assembly GRCh38 --json \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --hgvs --everything \
  -i input.vcf -o output.json
```

**Output formats.** VEP default (the tab-delimited text format with the `Extra` column), `--vcf` (annotations in the `CSQ` INFO field), `--json` (one object per variant) and `--tab` match Ensembl VEP's layouts. `--output_format parquet` is vep-rs's own, written through the `duckdb` CLI (1.3.2 or newer on PATH). Field-level specs are in [Output formats](../output-formats.md).

**Plugins.** `--plugin <Name>,<params>` resolves one of the eleven built-in plugins first (ten with supported data), and otherwise loads a C-ABI dynamic library from `--dir_plugins`. Support is per reference genome: `--assembly` selects the right column and refuses a data file whose own assembly contradicts it. The support matrix, the data files and how to write a plugin are in [Plugins](../plugins.md).

**HGVS.** `--fasta` takes a `samtools faidx`-indexed genome, plain or compressed with `bgzip` (Ensembl's `dna_index` files), and gives `--hgvs` the reference context for the 3' shifting of indels; `--hgvs` without a FASTA is refused, as Ensembl VEP refuses it offline.

**All flags.** The [CLI reference](../cli-reference.md) is the complete flag reference, including the flags accepted for Ensembl VEP compatibility but not implemented.
