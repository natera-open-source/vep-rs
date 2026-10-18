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
first `.fa` or `.fa.gz` by name inside the cache directory that has its indexes beside it is
used as the reference, unless `--no_fasta` is given; an unindexed file there is not read
(Ensembl VEP builds a missing index itself). `--show_cache_info` resolves the cache by the same
flags, refuses a named input file that does not exist as Ensembl VEP does, prints the cache's
source versions as `<key>\t<value>` lines sorted by key, and exits without reading the input.

## Reference FASTA

`--fasta` takes a plain FASTA with the `.fai` that `samtools faidx` writes, or a FASTA
compressed with `bgzip` with the `.fai` and `.gzi` that `samtools faidx` writes for it: Ensembl's
`dna_index/` files are of that shape. A gzip file that is not BGZF is refused with
`ERROR: [fasta] <file> is gzip but not BGZF: recompress with bgzip and index with samtools faidx`,
because its `.fai` offsets cannot be reached without inflating the whole file; no index is built
(Ensembl VEP builds one). The `.gzi` is checked against the file's block layout when it is opened,
so a file recompressed without being indexed again is refused rather than read. A directory given
to `--fasta` is read as Ensembl VEP's Bio::DB::Fasta reads one, every indexed `.fa` in it (and the
`.fa.gz` vep-rs adds, which Bio::DB::Fasta does not read), their sequences together; two files naming
the same sequence are refused (Ensembl VEP with Bio::DB::HTS installed takes a file only). Under
`VEP_FASTA_STORAGE=memory` a bgzipped FASTA is inflated whole into memory.

`--fasta_dir <dir>` chooses, as Ensembl VEP does, the `.fa` or `.fa.gz` of the directory whose
name carries the species and the assembly (case-insensitively), a file named `toplevel` before
any other and a compressed file before a plain one; the choice stands in for `--fasta` and for
the cache directory's FASTA. A directory holding no such file leaves those in force and warns
`WARNING: No compatible FASTA files found in <dir>`, written to the run's warnings file and echoed
on standard error unless `--quiet`; one that cannot be opened is passed over without a message, as
Ensembl VEP passes it over.

`--hgvs` or `--hgvsg` with no FASTA in reach exits 1 with
`ERROR: [fasta] Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file (see --fasta)`,
Ensembl VEP's offline refusal. `--everything` without a FASTA instead prints
`INFO: Disabling --hgvs; using --offline and no FASTA file found` on standard error (withheld
under `--quiet`) and drops the HGVS notations it switched on, as Ensembl VEP does; an explicit
`--hgvsg` beside it is still refused.

```bash
# Ensembl VEP's layout: $HOME/.vep/homo_sapiens/116_GRCh38
vep -i input.vcf -o output.txt --offline --cache --dir_cache $HOME/.vep --cache_version 116 --assembly GRCh38

# The directory itself
vep -i input.vcf -o output.txt --offline --json_cache /path/to/116_GRCh38
```

## Chromosome names

An input may name a chromosome as the cache does (`21`), with a `chr` prefix (`chr21`), or by
an accession (`NC_000021.9`, `CM000683.2`). Each name is resolved to the cache's own before any
lookup, as Ensembl VEP's `get_source_chr_name` resolves it: the name itself when the cache holds
it, else a synonym of it the cache holds, else the name with a leading `chr` removed or added
when the cache holds that, else a bare `M` as `MT`; a name that resolves to nothing skips the
record with Ensembl VEP's warning (`Chromosome NC_000022.11 not found in annotation sources or
synonyms`). The output keeps the input's spelling in `Location`, the VCF `CHROM` column and the
JSON `seq_region_name`; Ensembl VEP alone renames a bare `M` to `MT` in `Location` and
`seq_region_name`. The synonyms come from the cache's `chr_synonyms.txt`, the file the Ensembl
VEP cache ships and the converter copies (the native builder writes it from the core
`seq_region_synonym` table with `--core-dumps`), or from `--synonyms <file>`, which replaces the
cache's file rather than adding to it. The file holds, per line, whitespace-separated names of
one chromosome, the first the reference name and every other a synonym of it, in both
directions; a file that cannot be read is `ERROR: [config] Could not read synonyms file <file>:
<reason>` for a `--synonyms` file and `ERROR: [cache] ...` for the cache's own. A cache without
the file (a native build without `--core-dumps`, or a conversion made before the converter
carried it) resolves no synonym, so there `chrM` and an accession are skipped as Ensembl VEP
skips them on such a cache; convert the cache again, or pass `--synonyms`. Two cases differ from
Ensembl VEP, which checks a record against the cache's names widened by their synonyms and then
looks the unwidened name up: a bare `M` on a cache whose file links `chrM` but not `M`, and a
`chr`-prefixed accession (`chrNC_000021.9`), pass its check and annotate as intergenic; here
`M` resolves to `MT` and the prefixed accession is skipped. A reference FASTA whose records are
named by accession is read under the cache's names through the same table, so `--fasta` may
name a RefSeq or GenBank assembly file beside an Ensembl cache.

`info.json` carries `species` and `assembly` in every cache the converter and the builder write
(the resolver compares `--assembly` with it), the `sift` and `polyphen` markers only for a tool
whose matrices the cache holds, `source_versions` (on a native cache the core `meta` table's
`assembly`, `gencode` and `genebuild`, or the GTF header's; the dbSNP, ClinVar, COSMIC and gnomAD
versions an Ensembl cache's `info.txt` lists come from the variation database, which the builder
does not read, so a native cache's headers print none), `variation_cols` and `var_type tabix`
only when the cache holds variation files (`variation_cols` names their columns), and
`contig_order`, the cache's chromosome directories as autosomes by number, then X, Y and MT, then
every other contig by name.

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
  --output-dir $HOME/.vep/homo_sapiens/116_GRCh37
```

A table missing from the `--core-dumps` directory is downloaded into it: about 520 MB for
GRCh38 and 1.6 GB for GRCh37, whose `phenotype_feature.txt` Ensembl publishes
uncompressed (1.5 GB). The tables are read once per build, in about half a minute and
3 GB of memory for GRCh38.

Variation data, which `--check_existing` and the `--af*` flags read, is optional
(`--include-variations`). The builder downloads Ensembl's variation VCFs itself from
`https://ftp.ensembl.org/pub/release-116/variation/vcf/homo_sapiens/` (under
`pub/grch37/release-116/` for GRCh37), or reads local copies under the same names from the
directory passed with `--variation-vcf-dir`: the per-chromosome germline files
`homo_sapiens-chr<N>.vcf.gz` (the download covers chromosomes 1-22, X, Y and MT unless
`--chromosomes` narrows the list; from a local directory every `*-chr*.vcf.gz` file is read) and,
each with its `.csi` index, the genome-wide `homo_sapiens_somatic.vcf.gz` (every
record of it is flagged `somatic`), `homo_sapiens_phenotype_associated.vcf.gz` and
`homo_sapiens_clinically_associated.vcf.gz` (whose records the per-chromosome files already carry
at release 116; a record found only there is added) and `1000GENOMES-phase_3.vcf.gz`, the 1000
Genomes frequencies. A genome-wide file the directory or the site lacks is left out with a
warning naming what the cache then lacks. The records are written per chromosome as
`variations/<chr>/all_vars.gz` with a `.csi` on the start column, the file Ensembl VEP's own
cache carries, so `--check_existing` and the frequency flags read a native cache exactly as they
read a converted one. Each record takes the form the variation database holds: an indel's
alleles with the shared anchor base removed and `-` for an emptied allele; an insertion the
export marks `TSA=insertion` two bases after the base its VCF line names (`start` one past
`end`), where the Ensembl cache places it; a COSMIC variant without alleles as `COSMIC_MUTATION`
on its span and every HGMD-PUBLIC variant, alleles or not, as `HGMD_MUTATION`, the string the
variation database holds for that source; `clin_sig` from the record's `CLIN_<term>`
flags (a term the export declares no flag for, `likely_risk_allele`, is absent) and
`phenotype_or_disease` from `E_Phenotype_or_Disease`. The Ensembl cache keeps dbSNP's untrimmed
alleles for about seven indels in ten (`CCCCC/CCCC` at 5030241 where the export gives `C/-` at
5030242; 1,058,020 of chromosome 21's 1,500,839 indel records, 7 of every 100 records), and the
export's trimmed form is what a native cache holds: an input co-locates with such a record when
its own spelling trims to the same alleles and position, where the Ensembl cache answers the
untrimmed spelling instead. A native cache carries no `failed` flag (the export leaves the failed
variants out, the `dbSNP_novariation` records among them, and the few it keeps are written
unflagged), no `pubmed` or `var_synonyms` (the export carries none) and no gnomAD
frequency: Ensembl publishes no gnomAD VCF for release 116, so gnomAD frequencies come from the
Ensembl variation cache, by conversion. A native cache built by an earlier release holds JSON
variation files that nothing reads; build it again. The 1000 Genomes file exists for GRCh38
only, so a GRCh37 native cache has no frequency column; where it is read, the `AFR`, `AMR`,
`EAS`, `EUR` and `SAS` columns hold each allele's continental frequency as the file gives it and
`AF` the whole panel's: the file states no panel-wide frequency, so `AF` is the allele's count
over the panel's 5,008 chromosomes, each population's count recovered from its frequency and
sample count with two chromosomes per sample (an approximation on X and Y; a population the file
gives no value for drops out of both the count and the total). An allele whose frequency is 0 in
every population gets no entry: the release file lists no frequency for an allele absent from
its call set, and the export does not tell such an allele from one the release lists with a
count of 0, for which the Ensembl cache prints `0` (3,408 of the 1,106,899 chromosome 21 records
the Ensembl cache holds frequencies for). The release's own panel-wide frequency, which the
Ensembl cache prints, differs from the recovered count's for 443 of those records, 315 of them
by one unit in the last significant digit. `scripts/validation/compare_gtf_tags_to_cache.py`
compares the four attributes a GTF yields with those a cache carries and lists every
disagreement; `scripts/validation/diff_caches.py` compares two caches' transcripts, variation
files and `info.json`.

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
