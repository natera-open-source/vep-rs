# Cache setup

This page explains the two ways to get the JSON transcript cache vep-rs annotates from: building one with `vep-cache-builder`, or converting an existing Ensembl VEP cache. It then describes how a run names the cache, the reference FASTA, chromosome names, what `info.json` carries, and the SIFT and PolyPhen matrices.

The cache is a directory holding `info.json` and `transcripts/`: the release directory itself, not its parent, for example `$HOME/.vep/homo_sapiens/116_GRCh38`. A run names it with `--json_cache`, or with Ensembl VEP's own cache flags when the directory sits in Ensembl VEP's layout ([Naming the cache](#naming-the-cache)). For example:

```
--json_cache $HOME/.vep/homo_sapiens/116_GRCh38
```

## Build a cache

`vep-cache-builder` builds a cache natively, with no Perl dependency. It downloads the release's Ensembl GFF3 and protein FASTA, or reads local copies passed with `--gff3` and `--fasta`, and writes the cache directory. For GRCh37 the download fetches the release-87 gene set, the one Ensembl publishes for GRCh37 from release 88 onward. The GFF3 it retrieves is therefore `Homo_sapiens.GRCh37.87.gff3.gz` for any such `--release`.

The GFF3 alone gives each transcript its exons, CDS and ENSP identifier, its biotype, CCDS identifier and transcript support level, its GENCODE basic and primary, canonical and MANE memberships, and the gene's symbol and version. Ensembl VEP's own dumper drops `artifact` and readthrough transcripts, and so does the builder. Three more inputs decide what else the cache can express:

- `--core-dumps`, a directory of the release's Ensembl core and variation MySQL table dumps, supplies what the GFF3 does not carry. Ensembl publishes the dumps under `https://ftp.ensembl.org/pub/release-<N>/mysql/homo_sapiens_core_<N>_38/` and `homo_sapiens_variation_<N>_38/` for release `<N>`, under `pub/grch37/release-<N>/` with `_37` for GRCh37; any table missing from the directory is downloaded into it.

  The tables hold the MANE Select and MANE Plus Clinical accessions, APPRIS, the `cds_start_NF` and `cds_end_NF` flags, the mature miRNA arms and the `ncRNA` structure, and the RefSeq, UniProt, TrEMBL, UniParc and UniProt isoform identifiers. They also hold the protein features behind `DOMAINS`, the gene symbol's source and HGNC identifier, the codon table of each sequence region, the `source_versions` of the output headers and the gene phenotype count behind `GENE_PHENO`.

  Without it the `SYMBOL_SOURCE`, `HGNC_ID`, `MANE_SELECT`, `MANE_PLUS_CLINICAL`, `APPRIS`, `miRNA`, `RefSeq`, `SWISSPROT`, `TREMBL`, `UNIPARC`, `UNIPROT_ISOFORM`, `DOMAINS` and `GENE_PHENO` columns are empty on a native cache, and on GRCh37, whose release-87 GFF3 carries no `transcript_support_level` and no `Ensembl_canonical` tag, so are `TSL` and `CANONICAL`. A native cache lists a translation's protein features by feature identifier where the Ensembl cache lists them by analysis in Perl hash order, so `DOMAINS` may name the same features in another order. In the JSON output the `domains` array can also hold different objects: an entry that does not match Ensembl VEP's `(\w+):(\w+)` pattern, such as `Low_complexity_(Seg):seg`, repeats the `db` and `name` of the last matched entry, or is dropped when none precedes it.

  A native cache lists `cds_start_NF` before `cds_end_NF` on every transcript carrying both, whatever order the core dumps or the GTF carried, so `FLAGS` prints `cds_start_NF,cds_end_NF`. The Ensembl cache keeps the order its database returned, `cds_end_NF` first on 7 of the 2,946 GRCh38 transcripts that carry both at release 116 and on none of GRCh37's, and Ensembl VEP prints that order.
- `--gtf`, the Ensembl GTF of the same release, supplies the `cds_start_NF` and `cds_end_NF` transcript attributes (the `FLAGS` output field) when no core dumps are read; the GFF3 does not carry them. Its header also supplies the `assembly` and `genebuild` source versions the output headers print.
- `--genome-fasta`, the primary-assembly genome (decompressed, with a `samtools faidx` index beside it), supplies the CDS sequence behind every codon-level consequence: `missense_variant`, `synonymous_variant`, `stop_gained`, HGVSp and the rest. Without it, coding variants fall back to `coding_sequence_variant`.

Add `--include-predictions` to embed the SIFT and PolyPhen matrices that `--sift` and `--polyphen` read; see [SIFT and PolyPhen predictions](#sift-and-polyphen-predictions).

The GRCh38 command is in the [Quick start](getting-started/quick-start.md#2-get-a-transcript-cache). For GRCh37 the same download works: the release-87 files are `Homo_sapiens.GRCh37.87.gff3.gz` and `Homo_sapiens.GRCh37.87.gtf.gz` under `pub/grch37/release-116/`. Local copies can be passed instead, as below.

```bash
target/release/vep-cache-builder \
  --species homo_sapiens --assembly GRCh37 --release 116 \
  --gff3 Homo_sapiens.GRCh37.87.gff3.gz \
  --core-dumps core_dumps_grch37 \
  --genome-fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/homo_sapiens/116_GRCh37
```

A table missing from the `--core-dumps` directory is downloaded into it: about 520 MB for GRCh38 and 1.6 GB for GRCh37, whose `phenotype_feature.txt` Ensembl publishes uncompressed (1.5 GB). The tables are read once per build, in about half a minute and 3 GB of memory for GRCh38.

### Variation data

Variation data, which `--check_existing` and the `--af*` flags read, is optional (`--include-variations`). The builder downloads Ensembl's variation VCFs itself from this directory, under `pub/grch37/release-116/` for GRCh37, or reads local copies under the same names from the directory passed with `--variation-vcf-dir`:

```
https://ftp.ensembl.org/pub/release-116/variation/vcf/homo_sapiens/
```

The files are the per-chromosome germline files `homo_sapiens-chr<N>.vcf.gz` and four genome-wide files, each with its `.csi` index. The download covers chromosomes 1-22, X, Y and MT unless `--chromosomes` narrows the list; from a local directory every `*-chr*.vcf.gz` file is read. The genome-wide files are `homo_sapiens_somatic.vcf.gz` (every record of it is flagged `somatic`), `homo_sapiens_phenotype_associated.vcf.gz` and `homo_sapiens_clinically_associated.vcf.gz` (whose records the per-chromosome files already carry at release 116; a record found only there is added) and `1000GENOMES-phase_3.vcf.gz`, the 1000 Genomes frequencies. A genome-wide file the directory or the site lacks is left out with a warning naming what the cache then lacks.

The records are written per chromosome as `variations/<chr>/all_vars.gz` with a `.csi` on the start column, the file Ensembl VEP's own cache carries. `--check_existing` and the frequency flags therefore read a native cache exactly as they read a converted one.

Each record takes the form the variation database holds, as far as the export preserves it. A deletion's alleles have the anchor base the export prepended removed, with `-` for the emptied allele. An insertion the export marks `TSA=insertion` sits two bases after the base its VCF line names (`start` one past `end`), where the Ensembl cache places it. An indel (`TSA=indel`, every allele non-empty) keeps the alleles the export wrote, at the line's `POS`, since the export prepends an anchor only to alleles that do not share a first base. So an `A/AA` record stays `A/AA`, while a `CC/CTTCT` line, whose alleles differ once the shared base goes, is the database's `C/TTCT` one base on.

A COSMIC variant without alleles is written as `COSMIC_MUTATION` on its span, and every HGMD-PUBLIC variant, alleles or not, as `HGMD_MUTATION`, the string the variation database holds for that source. `clin_sig` comes from the record's `CLIN_<term>` flags (a term the export declares no flag for, `likely_risk_allele`, is absent) and `phenotype_or_disease` from `E_Phenotype_or_Disease`.

The Ensembl cache keeps dbSNP's alleles over the whole run an indel can shift along (`CCCCC/CCCC` at 5030241-5030245), which the export trims to `CC/C` at 5030241. Such records are 751,114 of chromosome 21's `indel` records, 5 of every 100 records. A native cache holds the trimmed form, which minimises from the right end to the same alleles and position as the Ensembl cache's record. So an input spelled at the left end of the run (the VCF convention, `AC/A` at 5030240) co-locates with it exactly as Ensembl VEP co-locates it.

An input anchored on the run's first base (`CC/C` at 5030241) minimises to the record's second base. The trimmed form also reaches that base and the Ensembl cache's record does not, so the input co-locates with the native record alone. An input spelled at the run's right end co-locates with neither, its start lying outside the one-base window the start-indexed files answer.

A native cache carries no `failed` flag: the export leaves the failed variants out, the `dbSNP_novariation` records among them, and the few it keeps are written unflagged. It carries no `pubmed` or `var_synonyms` (the export carries none) and no gnomAD frequency: Ensembl publishes no gnomAD VCF for release 116, so gnomAD frequencies come from the Ensembl variation cache, by conversion. A native cache built by an earlier vep-rs release holds JSON variation files that nothing reads; build it again.

The 1000 Genomes file exists for GRCh38 only, so a GRCh37 native cache has no frequency column. Where it is read, the `AFR`, `AMR`, `EAS`, `EUR` and `SAS` columns hold each allele's continental frequency as the file gives it, and `AF` the whole panel's. The file states no panel-wide frequency, so `AF` is the allele's count over the panel's 5,008 chromosomes. Each population's count is recovered from its frequency and sample count with two chromosomes per sample (an approximation on X and Y); a population the file gives no value for drops out of both the count and the total.

An allele whose frequency is 0 in every population gets no entry. The release file lists no frequency for an allele absent from its call set, and the export does not tell such an allele from one the release lists with a count of 0. For the latter the Ensembl cache prints `0` (3,408 of the 1,106,899 chromosome 21 records the Ensembl cache holds frequencies for). The release's own panel-wide frequency, which the Ensembl cache prints, differs from the recovered count's for 443 of those records, 315 of them by one unit in the last significant digit.

### Checking a cache

`scripts/validation/compare_gtf_tags_to_cache.py` compares the four attributes a GTF yields with those a cache carries and lists every disagreement. `scripts/validation/diff_caches.py` compares two caches' transcripts, variation files and `info.json`.

## Convert an existing Ensembl VEP cache

`scripts/data/storable_to_json.pl` reads an Ensembl VEP Storable cache directory, for example `~/.vep/homo_sapiens/116_GRCh38`, and writes the same JSON layout for every chromosome directory the cache holds: the assembled chromosomes and the patch, alternate-locus, unplaced-scaffold and `LRG_` directories alike, as Ensembl VEP itself lists them. Beside the transcripts it writes an `info.json` derived from the cache's `info.txt`, with `contig_order` (the converted contigs as autosomes by number, then X, Y and MT, then every other contig by name), and copies the cache's `chr_synonyms.txt`. The script runs inside the `ensemblorg/ensembl-vep` container ([scripts/README.md](../scripts/README.md) gives the invocation), and its output goes straight to `--json_cache`.

Every field the cache dumper keeps on a transcript is carried, flattened or dropped under a rule the parity audit names (`scripts/concordance/cache_parity_audit.pl`):

- Carried: the gene, transcript and protein identifiers (UniProt isoforms included, and on a merged cache the source), the attributes, the exons and introns, and the translation with its sequence edits. Also carried: the translateable sequence and peptide, the codon table, the exon coordinate mapper's pairs, and the SIFT and PolyPhen prediction matrices.
- Flattened: the UTR objects, to their sequence strings, and each protein feature's analysis, to the label Ensembl VEP's `DOMAINS` column prints.
- Dropped: the Perl object references (the slice, the gene object, the exon coordinate mapper's trees), the per-exon and per-intron strand, which is the transcript's, a sequence edit's description and the dump-time flags.

The variation files are copied as they are into `variations/<chr>/`: each chromosome's `all_vars.gz` with its tabix index or, for a cache dumped without the tabix conversion, its `<start>-<end>_var.gz` region files. `all_vars.gz` is the bgzip file the Ensembl VEP dumper writes, with the chromosome first, `.` for an empty cell and the columns of the cache's `variation_cols` line.

vep-rs opens each chromosome's `all_vars.gz` on first use and reads it through its index, or parses a region file whole the first time a variant falls in it. `--check_existing` and the frequency flags therefore read the converted cache as Ensembl VEP reads the original. A cache whose `variations/` directory holds no such file warns once that the co-location columns stay empty. A conversion made with an earlier version of the converter carries no `variations/<chr>/all_vars.gz`: convert the cache again, or copy each chromosome's `all_vars.gz` and its `.csi` from the Ensembl VEP cache into `variations/<chr>/`.

A cache written with the Sereal serialiser (`serialiser_type sereal` in its `info.txt`, or `.sereal` shards) is refused with an error naming the serialiser; only Storable caches convert. A run in which any shard fails to read or write stops with an error before `info.json` is written.

The published concordance figures were measured on such a conversion ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)), which holds the transcript set identical between the two engines. The caches those figures were measured on carry 195,232 transcripts (GRCh37) and 505,231 (GRCh38), the canonical chromosomes only. The GRCh37 conversion is of the cache-version-113 Storable that Ensembl VEP release 115 ships for that assembly. It differs from the cache-version-115 build in 310 of those transcripts (the `storable_113_*` columns of `manuscript/data/transcript_set_parity.csv`).

The conversion carries each translation's sequence edits (`seq_edits`: the initial-methionine, selenocysteine and amino-acid-substitution edits Ensembl VEP applies to the reference peptide at annotation time), and vep-rs applies them the same way. A cache converted without that key, the published conversions among them, still loads. On such a cache, a variant at codon 1 of a transcript with an alternative initiation codon is called from the cached peptide rather than from the edit. Converting the cache again restores Ensembl VEP's call. A natively built cache carries no edit list either: its peptide is the protein FASTA's sequence, with the edits already applied, and the builder writes none.

## Naming the cache

Ensembl VEP's cache flags select the directory the way Ensembl VEP does, so a command line written for it works unchanged ([Cache, species and runtime](cli/cache-species-runtime.md) lists the flags). The cache root is `--dir_cache`, else a directory given to `--cache`, else `--dir` (default `$HOME/.vep`). Under it the species directory is `--species` (default `homo_sapiens`, the Latin name only), suffixed `_refseq` or `_merged` under those flags. Under that sits the directory named `<cache_version>_<assembly>`, where `--cache_version` defaults to the VEP release.

When one directory of that version exists, it is the cache and its assembly is read from its name; when several exist, `--assembly` chooses. An `--assembly` that disagrees with the directory name, or a directory name or `--assembly` that disagrees with the cache's own `info.json`, is an error.

Every failure to resolve a cache exits 1 before the input is opened, with one `ERROR: [cache] ...` line in Ensembl VEP's words where it has them. Those words are `Cache directory <dir> not found`, `No cache found for homo_sapiens, version 116` and `Multiple assemblies found for cache version 116 ...`, the assemblies listed in name order. The metadata mismatch names `info.json`, the file a JSON cache carries where a Storable cache has `info.txt`.

`--json_cache <dir>` (alias `--full_cache_dir`) names the cache directory itself instead, whatever `--dir_cache`, `--species` and `--cache_version` say; only `--assembly` is checked against its `info.json`. Either way the directory must hold `transcripts/` and `info.json`, so a Storable cache or a mistyped path fails rather than annotating nothing.

`--assembly` is compared exactly and case-sensitively, as Ensembl VEP compares it, with one addition: `hg19` and `hg38` are read as `GRCh37` and `GRCh38`. When `--fasta` is absent, the first `.fa` or `.fa.gz` by name inside the cache directory that has its indexes beside it is used as the reference, unless `--no_fasta` is given. An unindexed file there is not read (Ensembl VEP builds a missing index itself). `--show_cache_info` resolves the cache by the same flags, refuses a named input file that does not exist as Ensembl VEP does, prints the cache's source versions as `<key>\t<value>` lines sorted by key, and exits without reading the input.

```bash
# Ensembl VEP's layout: $HOME/.vep/homo_sapiens/116_GRCh38
vep -i input.vcf -o output.txt --offline --cache \
  --dir_cache $HOME/.vep --cache_version 116 --assembly GRCh38

# The directory itself
vep -i input.vcf -o output.txt --offline --json_cache /path/to/116_GRCh38
```

## Reference FASTA

`--fasta` takes a plain FASTA with the `.fai` that `samtools faidx` writes, or a FASTA compressed with `bgzip` with the `.fai` and `.gzi` that `samtools faidx` writes for it. Ensembl's `dna_index/` files are of that shape. A gzip file that is not BGZF is refused with `ERROR: [fasta] <file> is gzip but not BGZF: recompress with bgzip and index with samtools faidx`, because its `.fai` offsets cannot be reached without inflating the whole file. No index is built (Ensembl VEP builds one).

The `.gzi` is checked against the file's block layout when it is opened, so a file recompressed without being indexed again is refused rather than read. A directory given to `--fasta` is read as Ensembl VEP's Bio::DB::Fasta reads one: every indexed `.fa` in it (and the `.fa.gz` vep-rs adds, which Bio::DB::Fasta does not read), their sequences together. Two files naming the same sequence are refused (Ensembl VEP with Bio::DB::HTS installed takes a file only). Under `VEP_FASTA_STORAGE=memory` a bgzipped FASTA is inflated whole into memory.

`--fasta_dir <dir>` chooses, as Ensembl VEP does, the `.fa` or `.fa.gz` of the directory whose name carries the species and the assembly (case-insensitively), a file named `toplevel` before any other and a compressed file before a plain one. The choice stands in for `--fasta` and for the cache directory's FASTA. A directory holding no such file leaves those in force and warns `WARNING: No compatible FASTA files found in <dir>`, written to the run's warnings file and echoed on standard error unless `--quiet`. A directory that cannot be opened is passed over without a message, as Ensembl VEP passes it over.

`--hgvs` or `--hgvsg` with no FASTA in reach exits 1 with `ERROR: [fasta] Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file (see --fasta)`, Ensembl VEP's offline refusal. `--everything` without a FASTA instead prints `INFO: Disabling --hgvs; using --offline and no FASTA file found` on standard error (withheld under `--quiet`) and drops the HGVS notations it switched on, as Ensembl VEP does. An explicit `--hgvsg` beside it is still refused. `--use_transcript_ref` without a FASTA exits 1 with `ERROR: [fasta] Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref`, since a non-coding transcript's cDNA is read from the FASTA.

## Chromosome names

An input may name a chromosome as the cache does (`21`), with a `chr` prefix (`chr21`), or by an accession (`NC_000021.9`, `CM000683.2`). Each name is resolved to the cache's own before any lookup, as Ensembl VEP's `get_source_chr_name` resolves it, in this order:

1. the name itself, when the cache holds it;
2. else a synonym of it the cache holds;
3. else the name with a leading `chr` removed or added, when the cache holds that;
4. else a bare `M` as the cache's `MT`, under whatever name the cache holds it.

A name that resolves to nothing skips the record with Ensembl VEP's warning (`Chromosome NC_000022.11 not found in annotation sources or synonyms`).

The output keeps the input's spelling in `Location`, in the JSON `seq_region_name` and in a VCF input's `CHROM` column, which is written as read. A bare `M` the cache knows as `MT` but not as `M` or `chrM` prints `MT` in `Location` and `seq_region_name`, as Ensembl VEP renames it; [A bare `M` under a cache whose synonyms name `chrM`](divergences/bare-m-chromosome.md) describes the cases where the two engines then differ.

The synonyms come from the cache's `chr_synonyms.txt`, the file the Ensembl VEP cache ships and the converter copies; the native builder writes it from the core `seq_region_synonym` table with `--core-dumps`. `--synonyms <file>` replaces the cache's file rather than adding to it. The file holds, per line, whitespace-separated names of one chromosome, the first the reference name and every other a synonym of it, in both directions. A file that cannot be read is `ERROR: [config] Could not read synonyms file <file>: <reason>` for a `--synonyms` file and `ERROR: [cache] ...` for the cache's own.

A cache without the file (a native build without `--core-dumps`, or a conversion made with an earlier version of the converter) resolves no synonym, so there `chrM` and an accession are skipped, as Ensembl VEP skips them on such a cache. Convert the cache again, or pass `--synonyms`.

Two cases differ from Ensembl VEP, which checks a record against the cache's names widened by their synonyms and then looks the unwidened name up. The first is a bare `M` on a cache whose file links `chrM` but not `M`. The second is a name that reaches the cache only through a synonym of its `chr`-stripped or `chr`-prefixed spelling (`chrNC_000021.9`, or `Q` under a line `21 chrQ`). Both pass its check and annotate as intergenic; in vep-rs `M` resolves to `MT` and the other name is skipped with the not-found warning.

A reference FASTA whose records are named by accession is read under the cache's names through the same table, so `--fasta` may name a RefSeq or GenBank assembly file beside an Ensembl cache.

## What `info.json` carries

`info.json` carries `species` and `assembly` in every cache the converter and the builder write; the resolver compares `--assembly` with `assembly`. Its other keys:

- `sift` and `polyphen`, markers written only for a tool whose matrices the cache holds.
- `source_versions`: on a native cache the core `meta` table's `assembly`, `gencode` and `genebuild`, or the GTF header's. The dbSNP, ClinVar, COSMIC and gnomAD versions an Ensembl cache's `info.txt` lists come from the variation database, which the builder does not read, so a native cache's headers print none.
- `variation_cols` and `var_type tabix`, only when the cache holds variation files; `variation_cols` names their columns.
- `contig_order`: the cache's chromosome directories as autosomes by number, then X, Y and MT, then every other contig by name.

## SIFT and PolyPhen predictions

vep-rs supports Ensembl VEP's `--sift` and `--polyphen` flags with their `p` (prediction), `s` (score) and `b` (both) values, the only three accepted, from matrices embedded in the cache. `vep-cache-builder --include-predictions` fetches per-transcript prediction matrix blobs from Ensembl's public MySQL server and stores them base64-encoded and gzip-compressed beside each transcript. The server is `ensembldb.ensembl.org`, anonymous read-only; port 3306 serves the GRCh38 databases and `--mysql-port 3337` the GRCh37 ones. The database is `homo_sapiens_variation_<release>_<assembly>`, the tables `translation_md5`, `protein_function_predictions` and `attrib`, keyed on the MD5 of the full peptide.

At runtime a pure-Rust decoder (`crates/vep-core/src/prediction.rs`) looks up `(protein_position, alt_amino_acid)` and prints the prediction as Ensembl VEP prints it: `deleterious(0.01)`, `probably_damaging(0.998)`. The score is `prob/1000`, so `0.01` rather than `0.010`; under `s` the bare score is printed, under `p` the bare prediction.

A prediction is looked up for every single-residue substitution inside the coding sequence, the shape Ensembl VEP predicts on: a `missense_variant`, a `start_lost` whose codon still translates (`M/T`), a same-length substitution of several bases within one codon. It is never looked up for a synonymous change (whose peptide allele string is the one residue, `C`, not `C/C`), a stop gain or loss, an unknown alternate residue, an insertion or a deletion. The codon-level consequences themselves need the CDS sequence `--genome-fasta` supplies.

PolyPhen reads the HumVar matrix, as Ensembl VEP does by default; `--humdiv` reads the HumDiv matrix instead, and a cache without that matrix prints an empty column, with no error in either engine.

vep-rs checks the tools against the cache's metadata before any record is read, as Ensembl VEP does. A cache whose `info.json` (or Ensembl VEP's `info.txt`) names no `source_versions` entry for the tool stops a run that asks for it, with Ensembl VEP's text: `ERROR: [cache] SIFT not available` or `ERROR: [cache] PolyPhen not available`. Under `--everything` or `--safe` the run instead drops the tool with `INFO: disabling SIFT` on standard error (withheld under `--quiet`) and writes no column for it.

A natively built cache names no tool version; its `sift` and `polyphen` markers, written only when the build embedded a matrix of that tool, count as availability instead. A native cache built without `--include-predictions` is therefore the one case the flags refuse that Ensembl VEP's own caches never present.

```bash
vep-cache-builder --release 116 --assembly GRCh38 --output-dir /path/to/cache \
  --gff3 ... --gtf ... --include-predictions \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
# GRCh37: --assembly GRCh37 --mysql-port 3337

vep -i input.vcf --sift b --polyphen b --json_cache /path/to/cache --offline
```

A cache converted from an Ensembl VEP Storable cache carries the matrices the Storable holds: the SIFT, PolyPhen HumVar and PolyPhen HumDiv blobs stored for a translated transcript. They are base64-encoded under the same `predictions_data` key and in the same binary layout the native builder writes, so `--sift` and `--polyphen` read either kind of cache alike. `info.json` keeps the `sift` and `polyphen` markers only when the conversion wrote a matrix of that tool. The `dbNSFP` plugin is a separate, unsupported route to third-party SIFT and PolyPhen scores ([Plugins](plugins.md)).
