# Cache setup

This page explains the two ways to get the JSON transcript cache that `--json_cache` reads: building one with `vep-cache-builder`, or converting an existing Ensembl VEP cache.

vep-rs annotates from a JSON transcript cache passed as a directory with `--json_cache`. That directory is the one holding `info.json` and `transcripts/`: the release directory itself, not its parent. For example:

```
--json_cache $HOME/.vep/json_cache/homo_sapiens/116_GRCh38
```

## Build a cache

`vep-cache-builder` builds a cache natively, with no Perl dependency. It downloads the release's Ensembl GFF3 and protein FASTA, or reads local copies passed with `--gff3` and `--fasta`, and writes the cache directory. For GRCh37 the download fetches the release-87 gene set, the one Ensembl publishes for GRCh37 from release 88 onward. The GFF3 it retrieves is therefore `Homo_sapiens.GRCh37.87.gff3.gz` whatever `--release` says.

Two more inputs decide what the cache can express:

- `--gtf`, the Ensembl GTF of the same release, is the only source of the `cds_start_NF` and `cds_end_NF` transcript attributes (the `FLAGS` output field), which the GFF3 does not carry. Without it the cache has no FLAGS. The GTF header also supplies the `assembly` and `genebuild` source versions the output headers print.
- `--genome-fasta`, the primary-assembly genome (decompressed, with a `samtools faidx` index beside it), supplies the CDS sequence behind every codon-level consequence: `missense_variant`, `synonymous_variant`, `stop_gained`, HGVSp and the rest. Without it, coding variants fall back to `coding_sequence_variant`.

Add `--include-predictions` to embed the SIFT and PolyPhen matrices that `--sift` and `--polyphen` read; see [SIFT and PolyPhen predictions](#sift-and-polyphen-predictions).

The GRCh38 command is in the [Quick start](getting-started/quick-start.md#2-get-a-transcript-cache). For GRCh37 the same download works. Ensembl ships the GRCh37 gene set frozen at release 87, as `Homo_sapiens.GRCh37.87.gff3.gz` and `Homo_sapiens.GRCh37.87.gtf.gz` under `pub/grch37/release-116/`, and the builder names that release-87 file for any `--release` from 88 onward. Local copies can be passed instead, as below.

```bash
target/release/vep-cache-builder \
  --species homo_sapiens --assembly GRCh37 --release 116 \
  --gff3 Homo_sapiens.GRCh37.87.gff3.gz \
  --gtf Homo_sapiens.GRCh37.87.gtf.gz \
  --genome-fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/json_cache/homo_sapiens/116_GRCh37
```

Variation data, which `--check_existing` and the `--af*` flags read, is optional (`--include-variations`). The builder downloads Ensembl's per-chromosome variation VCFs itself, one file per chromosome from this directory, under `pub/grch37/release-116/` for GRCh37:

```
https://ftp.ensembl.org/pub/release-116/variation/vcf/homo_sapiens/
homo_sapiens-chr<N>.vcf.gz
```

It fetches chromosomes 1-22, X, Y and MT unless `--chromosomes` narrows the list. Local copies are read instead from the directory passed with `--variation-vcf-dir`.

`scripts/validation/compare_gtf_tags_to_cache.py` compares the four attributes a GTF yields with those a cache carries and lists every disagreement.

## Convert an existing Ensembl VEP cache

`scripts/data/storable_to_json.pl` reads an Ensembl VEP Storable cache directory, for example `~/.vep/homo_sapiens/116_GRCh38`, and writes the same JSON layout for chromosomes 1-22, X, Y and MT. Its `info.json` is derived from the cache's `info.txt`. The script runs inside the `ensemblorg/ensembl-vep` container ([scripts/README.md](../scripts/README.md) gives the invocation), and its output goes straight to `--json_cache`.

The published concordance figures were measured on such a conversion ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)), which holds the transcript set identical between the two engines. The converted caches carry 195,232 transcripts (GRCh37) and 505,231 (GRCh38). The GRCh37 conversion is of the cache-version-113 Storable that Ensembl VEP release 115 ships for that assembly. It differs from the cache-version-115 build in 310 of those transcripts (the `storable_113_*` columns of `manuscript/data/transcript_set_parity.csv`). A converted cache carries no variation data and no SIFT or PolyPhen matrices.

The conversion carries each translation's sequence edits (`seq_edits`: the initial-methionine, selenocysteine and amino-acid-substitution edits Ensembl VEP applies to the reference peptide at annotation time), and vep-rs applies them the same way. A cache converted without that key, the published conversions among them, still loads. On such a cache, a variant at codon 1 of a transcript with an alternative initiation codon is called from the cached peptide rather than from the edit. Converting the cache again restores Ensembl VEP's call. A natively built cache carries no edit list either: its peptide is the protein FASTA's sequence, with the edits already applied, and none of the builder's inputs carries the list.

## SIFT and PolyPhen predictions

vep-rs supports Ensembl VEP's `--sift` and `--polyphen` flags with their `p` (prediction), `s` (score) and `b` (both) modes, from matrices embedded in the cache. `vep-cache-builder --include-predictions` fetches per-transcript prediction matrix blobs from Ensembl's public MySQL server and stores them base64-encoded and gzip-compressed beside each transcript. The server is `ensembldb.ensembl.org`, anonymous read-only; port 3306 serves the GRCh38 databases and `--mysql-port 3337` the GRCh37 ones. The database is `homo_sapiens_variation_<release>_<assembly>`, the tables `translation_md5`, `protein_function_predictions` and `attrib`, keyed on the MD5 of the full peptide.

At runtime a pure-Rust decoder (`crates/vep-core/src/prediction.rs`) looks up `(protein_position, alt_amino_acid)` for each missense variant and prints the prediction as Ensembl VEP prints it: `deleterious(0.01)`, `probably_damaging(0.998)`. The score is `prob/1000`, so `0.01` rather than `0.010`. PolyPhen uses the HumVar matrix, as Ensembl VEP does by default (`--humdiv` is not accepted). A prediction needs a protein position and a single alternate amino acid, so in practice a `missense_variant`, which itself needs the CDS sequence `--genome-fasta` supplies.

```bash
vep-cache-builder --release 116 --assembly GRCh38 --output-dir /path/to/cache \
  --gff3 ... --gtf ... --include-predictions \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
# GRCh37: --assembly GRCh37 --mysql-port 3337

vep -i input.vcf --sift b --polyphen b --json_cache /path/to/cache --offline
```

A cache converted from an Ensembl VEP Storable cache carries no matrices, because the converter reads a `predictions_data` field that only the native builder writes. On a converted cache `--sift` and `--polyphen` print nothing. The `dbNSFP` plugin is a separate, unsupported route to third-party SIFT and PolyPhen scores ([Plugins](plugins.md)).
