# Cache setup: building or converting the JSON transcript cache

vep-rs annotates from a JSON transcript cache passed as a directory with `--json_cache`: the directory that holds `info.json` and `transcripts/` (for example `$HOME/.vep/json_cache/homo_sapiens/115_GRCh38` itself, not its parent). There are two ways to get one.

**Build one natively with `vep-cache-builder`**, with no Perl dependency. It downloads the
release's Ensembl GFF3 and protein FASTA, or reads local copies passed with `--gff3` and
`--fasta`, and writes the cache directory. For GRCh37 the download fetches the release-87 gene set, the one Ensembl publishes for GRCh37 from release 88 onward, so the GFF3 it retrieves is `Homo_sapiens.GRCh37.87.gff3.gz` whatever `--release` says. Two more inputs decide what the cache can express:

- `--gtf`, the Ensembl GTF of the same release, is the only source of the `cds_start_NF`
  and `cds_end_NF` transcript attributes (the `FLAGS` output field), which the GFF3 does
  not carry; without it the cache has no FLAGS. The GTF header also supplies the
  `assembly` and `genebuild` source versions the output headers print.
- `--genome-fasta`, the primary-assembly genome (decompressed, with a `samtools faidx`
  index beside it), supplies the CDS sequence behind every codon-level consequence
  (`missense_variant`, `synonymous_variant`, `stop_gained`, HGVSp and the rest); without
  it, coding variants fall back to `coding_sequence_variant`.

Add `--include-predictions` to embed the SIFT and PolyPhen matrices that `--sift` and
`--polyphen` read; see [SIFT and PolyPhen predictions](#sift-and-polyphen-predictions).

The GRCh38 command is in the [README's quick start](../README.md#quick-start). For GRCh37
pass the GFF3 explicitly with `--gff3`: Ensembl ships the GRCh37 gene set frozen at
release 87, as `Homo_sapiens.GRCh37.87.gff3.gz` and `Homo_sapiens.GRCh37.87.gtf.gz`
under `pub/grch37/release-115/`, and the builder's own GRCh37 GFF3 download looks for
a release-115 file name that is not there.

```bash
target/release/vep-cache-builder \
  --species homo_sapiens --assembly GRCh37 --release 115 \
  --gff3 Homo_sapiens.GRCh37.87.gff3.gz \
  --gtf Homo_sapiens.GRCh37.87.gtf.gz \
  --genome-fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/json_cache/homo_sapiens/115_GRCh37
```

Variation data, which `--check_existing` and the `--af*` flags read, is optional
(`--include-variations`). The builder downloads Ensembl's per-chromosome variation VCFs
itself (`https://ftp.ensembl.org/pub/release-115/variation/vcf/homo_sapiens/homo_sapiens-chr<N>.vcf.gz`,
under `pub/grch37/release-115/` for GRCh37; chromosomes 1-22, X, Y and MT unless
`--chromosomes` narrows the list), or reads local copies from the directory passed with
`--variation-vcf-dir`. `scripts/validation/compare_gtf_tags_to_cache.py`
compares the four attributes a GTF yields with those a cache carries and lists every
disagreement.

**Convert an existing Perl VEP cache.** `scripts/data/storable_to_json.pl` reads a Perl VEP
Storable cache directory (for example `~/.vep/homo_sapiens/115_GRCh38`) and writes the same
JSON layout for chromosomes 1-22, X, Y and MT, with an `info.json` derived from the cache's
`info.txt`; it runs inside the `ensemblorg/ensembl-vep` container
([scripts/README.md](../scripts/README.md) gives the invocation), and its output goes straight
to `--json_cache`. The published concordance figures were measured on that conversion ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)), which holds the transcript set identical between the two engines; the converted caches carry 195,232 transcripts (GRCh37) and 505,231
(GRCh38); the GRCh37 conversion is of the cache-version-113 Storable that VEP release 115 ships
for that assembly, which differs from the cache-version-115 build in 310 of those transcripts
(the `storable_113_*` columns of `manuscript/data/transcript_set_parity.csv`). A converted cache carries no variation data and no SIFT
or PolyPhen matrices.
`vep-cache-converter --json-dir <converted> --output-dir <dest>` copies a converted cache
into the runtime layout and sanity-checks it on the way.

## SIFT and PolyPhen predictions

vep-rs supports Perl VEP's `--sift` and `--polyphen` flags with their `p` (prediction), `s` (score) and `b` (both) modes from matrices embedded in the cache. `vep-cache-builder --include-predictions` fetches per-transcript prediction matrix blobs from Ensembl's public MySQL server (`ensembldb.ensembl.org`, anonymous read-only; port 3306 serves the GRCh38 databases and `--mysql-port 3337` the GRCh37 ones; database `homo_sapiens_variation_<release>_<assembly>`, tables `translation_md5`, `protein_function_predictions` and `attrib`, keyed on the MD5 of the full peptide) and stores them base64-encoded and gzip-compressed beside each transcript. At runtime a pure-Rust decoder (`crates/vep-core/src/prediction.rs`) looks up `(protein_position, alt_amino_acid)` for each missense variant and prints the prediction as Perl prints it (`deleterious(0.01)`, `probably_damaging(0.998)`; the score is `prob/1000`, so `0.01` rather than `0.010`). PolyPhen uses the HumVar matrix, as Perl VEP does by default (`--humdiv` is not accepted). A prediction needs a protein position and a single alternate amino acid, so in practice a `missense_variant`, which itself needs the CDS sequence `--genome-fasta` supplies.

```bash
vep-cache-builder --release 115 --assembly GRCh38 --output-dir /path/to/cache \
  --gff3 ... --gtf ... --include-predictions \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
# GRCh37: --assembly GRCh37 --mysql-port 3337

vep -i input.vcf --sift b --polyphen b --json_cache /path/to/cache --offline
```

A cache converted from a Perl VEP Storable cache carries no matrices: the converter reads a `predictions_data` field that only the native builder writes, so `--sift` and `--polyphen` print nothing on a converted cache. The `dbNSFP` plugin is a separate, unsupported route to third-party SIFT and PolyPhen scores ([plugins.md](plugins.md#not-supported-dbnsfp)).
