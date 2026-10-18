# Cache setup: building or converting the JSON transcript cache

vep-rs annotates from a JSON transcript cache: a directory that holds `info.json` and `transcripts/` (for example `$HOME/.vep/homo_sapiens/116_GRCh38`). There are two ways to name it and two ways to get one.

## Naming the cache

Ensembl VEP's cache flags select the directory the way Ensembl VEP does, so a command line
written for it works unchanged. The cache root is `--dir_cache`, else a directory given to
`--cache`, else `--dir` (default `$HOME/.vep`); under it the species directory is `--species`
(default `homo_sapiens`, the Latin name only), suffixed `_refseq` or `_merged` under those
flags; under that, the directory named `<cache_version>_<assembly>`, where `--cache_version`
defaults to the VEP release. When one directory of that version exists, it is the cache and its
assembly is read from its name; when several exist, `--assembly` chooses; an `--assembly` that
disagrees with the directory name, or a directory name or `--assembly` that disagrees with the
cache's own `info.json`, is an error. Every
failure to resolve a cache exits 1 before the input is opened with one `ERROR: [cache] ...`
line, in Ensembl VEP's words where it has them (`Cache directory <dir> not found`, `No cache
found for homo_sapiens, version 116`, `Multiple assemblies found for cache version 116 ...`,
the assemblies listed in name order; the metadata mismatch names `info.json`, the file a
JSON cache carries where a Storable cache has `info.txt`).

`--json_cache <dir>` (alias `--full_cache_dir`) names the cache directory itself instead,
whatever `--dir_cache`, `--species` and `--cache_version` say; only `--assembly` is checked
against its `info.json`. Either way the directory must hold `transcripts/` and `info.json`, so a
Storable cache or a mistyped path fails rather than annotating nothing.

`--assembly` is compared exactly and case-sensitively, as Ensembl VEP compares it, with one
addition: `hg19` and `hg38` are read as `GRCh37` and `GRCh38`. When `--fasta` is absent, the
first `.fa` by name inside the cache directory that has its `.fai` index beside it is used as
the reference, unless `--no_fasta` is given; a `.fa.gz` or an unindexed `.fa` there is not read
(Ensembl VEP reads a `.fa.gz` through Bio::DB::HTS and builds a missing index itself).
`--show_cache_info` resolves the cache by the same flags, refuses a named input file that does
not exist as Ensembl VEP does, prints the cache's source versions as `<key>\t<value>` lines
sorted by key, and exits without reading the input.

```bash
# Ensembl VEP's layout: $HOME/.vep/homo_sapiens/116_GRCh38
vep -i input.vcf -o output.txt --offline --cache --dir_cache $HOME/.vep --cache_version 116 --assembly GRCh38

# The directory itself
vep -i input.vcf -o output.txt --offline --json_cache /path/to/116_GRCh38
```

## Getting a cache

**Build one natively with `vep-cache-builder`**, with no Perl dependency. It downloads the
release's Ensembl GFF3 and protein FASTA, or reads local copies passed with `--gff3` and
`--fasta`, and writes the cache directory. For GRCh37 the download fetches the release-87 gene set, the one Ensembl publishes for GRCh37 from release 88 onward, so for any such `--release` the GFF3 it retrieves is `Homo_sapiens.GRCh37.87.gff3.gz`. The GFF3 alone gives each transcript its
exons, CDS and ENSP identifier, its biotype, CCDS identifier and transcript support level,
its GENCODE basic and primary, canonical and MANE memberships, and the gene's symbol and
version; Ensembl VEP's own dumper drops `artifact` and readthrough transcripts, and so
does the builder. Three more inputs decide what else the cache can express:

- `--core-dumps`, a directory of the release's Ensembl core and variation MySQL table
  dumps (`https://ftp.ensembl.org/pub/release-<N>/mysql/homo_sapiens_core_<N>_38/` and
  `homo_sapiens_variation_<N>_38/` for release `<N>`, under `pub/grch37/release-<N>/`
  with `_37` for GRCh37; any table missing from the directory is downloaded into it),
  supplies what the GFF3 does not carry: the MANE Select and MANE Plus Clinical accessions, APPRIS, the
  `cds_start_NF` and `cds_end_NF` flags, the mature miRNA arms and the `ncRNA`
  structure, the RefSeq, UniProt, TrEMBL, UniParc and UniProt isoform identifiers, the
  protein features behind `DOMAINS`, the gene symbol's source and HGNC identifier, the
  codon table of each sequence region, the `source_versions` of the output headers and
  the gene phenotype count behind `GENE_PHENO`. Without it the `SYMBOL_SOURCE`, `HGNC_ID`,
  `MANE_SELECT`, `MANE_PLUS_CLINICAL`, `APPRIS`, `miRNA`, `RefSeq`, `SWISSPROT`,
  `TREMBL`, `UNIPARC`, `UNIPROT_ISOFORM`, `DOMAINS` and `GENE_PHENO` columns are empty on
  a native cache, and on GRCh37, whose release-87 GFF3 carries no
  `transcript_support_level` and no `Ensembl_canonical` tag, so are `TSL` and
  `CANONICAL`. A native cache lists a translation's protein features by feature
  identifier where the Ensembl cache keeps the database's row order, so `DOMAINS` may
  name the same features in another order.
- `--gtf`, the Ensembl GTF of the same release, supplies, when no core dumps are read,
  the `cds_start_NF` and `cds_end_NF` transcript attributes (the `FLAGS` output field)
  and, from its header, the `assembly` and `genebuild` source versions.
- `--genome-fasta`, the primary-assembly genome (decompressed, with a `samtools faidx`
  index beside it), supplies the CDS sequence behind every codon-level consequence
  (`missense_variant`, `synonymous_variant`, `stop_gained`, HGVSp and the rest); without
  it, coding variants fall back to `coding_sequence_variant`.

Add `--include-predictions` to embed the SIFT and PolyPhen matrices that `--sift` and
`--polyphen` read; see [SIFT and PolyPhen predictions](#sift-and-polyphen-predictions).

The GRCh38 command is in the [README's quick start](../README.md#quick-start). For GRCh37
the same download works: Ensembl ships the GRCh37 gene set frozen at release 87, as
`Homo_sapiens.GRCh37.87.gff3.gz` and `Homo_sapiens.GRCh37.87.gtf.gz` under
`pub/grch37/release-116/`, and the builder names that release-87 file for any `--release`
from 88 onward. Local copies can be passed instead, as below.

```bash
target/release/vep-cache-builder \
  --species homo_sapiens --assembly GRCh37 --release 116 \
  --gff3 Homo_sapiens.GRCh37.87.gff3.gz \
  --core-dumps core_dumps_grch37 \
  --genome-fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/json_cache/homo_sapiens/116_GRCh37
```

A table missing from the `--core-dumps` directory is downloaded into it: about 520 MB for
GRCh38 and 1.6 GB for GRCh37, whose `phenotype_feature.txt` Ensembl publishes
uncompressed (1.5 GB). The tables are read once per build, in about half a minute and
3 GB of memory for GRCh38.

Variation data, which `--check_existing` and the `--af*` flags read, is optional
(`--include-variations`). The builder downloads Ensembl's per-chromosome variation VCFs
itself (`https://ftp.ensembl.org/pub/release-116/variation/vcf/homo_sapiens/homo_sapiens-chr<N>.vcf.gz`,
under `pub/grch37/release-116/` for GRCh37; chromosomes 1-22, X, Y and MT unless
`--chromosomes` narrows the list), or reads local copies from the directory passed with
`--variation-vcf-dir`. `scripts/validation/compare_gtf_tags_to_cache.py`
compares the four attributes a GTF yields with those a cache carries and lists every
disagreement.

**Convert an existing Perl VEP cache.** `scripts/data/storable_to_json.pl` reads a Perl VEP
Storable cache directory (for example `~/.vep/homo_sapiens/116_GRCh38`) and writes the same
JSON layout for every chromosome directory the cache holds, the assembled chromosomes and the
patch, alternate-locus, unplaced-scaffold and `LRG_` directories alike, as Ensembl VEP itself
lists them; beside the transcripts it writes an `info.json` derived from the cache's `info.txt`
(with `contig_order`, the converted contigs as autosomes by number, then X, Y and MT, then every
other contig by name) and copies the cache's `chr_synonyms.txt`. It runs inside the
`ensemblorg/ensembl-vep` container ([scripts/README.md](../scripts/README.md) gives the
invocation), and its output goes straight to `--json_cache`. Every field the cache dumper keeps
on a transcript is carried, flattened or dropped under a rule the parity audit names
(`scripts/concordance/cache_parity_audit.pl`): carried are the gene, transcript and protein
identifiers (UniProt isoforms included, and on a merged cache the source), the attributes, the
exons and introns, the translation with its sequence edits, the translateable sequence and
peptide, the codon table, the exon coordinate mapper's pairs and the SIFT and PolyPhen prediction
matrices; flattened are the UTR objects, to their sequence strings, and each
protein feature's analysis, to the label Ensembl VEP's `DOMAINS` column prints; dropped are the
Perl object references (the slice, the gene object, the exon coordinate mapper's trees), the
per-exon and per-intron strand, which is the transcript's, a sequence edit's description and the
dump-time flags. The variation files are copied as they are into `variations/<chr>/`: each
chromosome's `all_vars.gz` with its tabix index (the bgzip file the Ensembl VEP dumper writes,
with the chromosome first, `.` for an empty cell and the columns of the cache's
`variation_cols` line) or, for a cache dumped without the tabix conversion, its
`<start>-<end>_var.gz` region files. vep-rs opens each chromosome's `all_vars.gz` on first use
and reads it through its index, or parses a region file whole the first time a variant falls in
it, so `--check_existing` and the frequency flags read the converted cache as Ensembl VEP reads
the original; a cache whose `variations/` directory holds no such file warns once that the
co-location columns stay empty. A conversion made by an earlier release carries no
`variations/<chr>/all_vars.gz`: convert the cache again, or copy each chromosome's `all_vars.gz`
and its `.csi` from the Ensembl VEP cache into `variations/<chr>/`. A cache written with the Sereal serialiser (`serialiser_type sereal` in its
`info.txt`, or `.sereal` shards) is refused with an error naming the serialiser; only Storable
caches convert, and a run in which any shard fails to read or write stops with an error before
`info.json` is written.
The published concordance figures were measured on such a conversion ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)), which holds the transcript set identical between the two engines; the caches those figures were measured on carry 195,232 transcripts (GRCh37) and 505,231
(GRCh38), the canonical chromosomes only, and the GRCh37 conversion is of the cache-version-113 Storable that VEP release 115 ships
for that assembly, which differs from the cache-version-115 build in 310 of those transcripts
(the `storable_113_*` columns of `manuscript/data/transcript_set_parity.csv`).

The conversion carries each translation's sequence edits (`seq_edits`: the initial-methionine,
selenocysteine and amino-acid-substitution edits VEP applies to the reference peptide at
annotation time), and vep-rs applies them the same way. A cache converted without that key, the
published conversions among them, still loads, but on it a variant at codon 1 of a transcript
with an alternative initiation codon is called from the cached peptide rather than from the edit;
convert such a cache again to get VEP's call. A natively built cache carries no edit list either:
its peptide is the protein FASTA's sequence, with the edits already applied, and the builder
writes none.

## SIFT and PolyPhen predictions

vep-rs supports Perl VEP's `--sift` and `--polyphen` flags with their `p` (prediction), `s` (score) and `b` (both) modes from matrices embedded in the cache. `vep-cache-builder --include-predictions` fetches per-transcript prediction matrix blobs from Ensembl's public MySQL server (`ensembldb.ensembl.org`, anonymous read-only; port 3306 serves the GRCh38 databases and `--mysql-port 3337` the GRCh37 ones; database `homo_sapiens_variation_<release>_<assembly>`, tables `translation_md5`, `protein_function_predictions` and `attrib`, keyed on the MD5 of the full peptide) and stores them base64-encoded and gzip-compressed beside each transcript. At runtime a pure-Rust decoder (`crates/vep-core/src/prediction.rs`) looks up `(protein_position, alt_amino_acid)` for each missense variant and prints the prediction as Perl prints it (`deleterious(0.01)`, `probably_damaging(0.998)`; the score is `prob/1000`, so `0.01` rather than `0.010`). PolyPhen uses the HumVar matrix, as Perl VEP does by default (`--humdiv` is not accepted). A prediction needs a protein position and a single alternate amino acid, so in practice a `missense_variant`, which itself needs the CDS sequence `--genome-fasta` supplies.

```bash
vep-cache-builder --release 116 --assembly GRCh38 --output-dir /path/to/cache \
  --gff3 ... --gtf ... --include-predictions \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
# GRCh37: --assembly GRCh37 --mysql-port 3337

vep -i input.vcf --sift b --polyphen b --json_cache /path/to/cache --offline
```

A cache converted from a Perl VEP Storable cache carries the matrices the Storable holds, the SIFT, PolyPhen HumVar and PolyPhen HumDiv blobs stored for a translated transcript, base64-encoded under the same `predictions_data` key and in the same binary layout the native builder writes, so `--sift` and `--polyphen` read either kind of cache alike; `info.json` keeps the `sift` and `polyphen` markers only when the conversion wrote a matrix of that tool. The `dbNSFP` plugin is a separate, unsupported route to third-party SIFT and PolyPhen scores ([plugins.md](plugins.md#not-supported-dbnsfp)).
