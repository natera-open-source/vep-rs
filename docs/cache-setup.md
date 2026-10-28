# Cache setup

The JSON transcript cache vep-rs annotates from is built with `vep-cache-builder` or converted from an existing Ensembl VEP cache. [Cache checksum manifest](cache-checksums.md) covers integrity checks, and [SIFT and PolyPhen predictions](cache-predictions.md) covers the prediction matrices.

The cache is a directory holding `info.json` and `transcripts/`: the release directory itself, not its parent, for example `$HOME/.vep/homo_sapiens/116_GRCh38`. A run names it with `--json_cache`, or with Ensembl VEP's own cache flags when the directory sits in Ensembl VEP's layout ([Naming the cache](#naming-the-cache)). For example:

```
--json_cache $HOME/.vep/homo_sapiens/116_GRCh38
```

## Build a cache

`vep-cache-builder` builds a cache natively, with no Perl dependency. It downloads the release's Ensembl GFF3 and protein FASTA, or reads local copies passed with `--gff3` and `--fasta`, and writes the cache directory. For GRCh37 the download fetches the release-87 gene set, the one Ensembl publishes for GRCh37 from release 88 onward.

The GFF3 alone gives each transcript its exons, CDS and ENSP identifier, its biotype, CCDS identifier and transcript support level, its GENCODE basic and primary, canonical and MANE memberships, and the gene's symbol and version. Ensembl VEP's own dumper drops `artifact` and readthrough transcripts, and so does the builder. Three more inputs decide what else the cache can express:

- `--core-dumps`, a directory of the release's Ensembl core and variation MySQL table dumps, supplies what the GFF3 does not carry. Ensembl publishes the dumps under `https://ftp.ensembl.org/pub/release-<N>/mysql/homo_sapiens_core_<N>_38/` and `homo_sapiens_variation_<N>_38/` for release `<N>`, under `pub/grch37/release-<N>/` with `_37` for GRCh37; any table missing from the directory is downloaded into it.

  The tables hold the MANE Select and MANE Plus Clinical accessions, APPRIS, the `cds_start_NF` and `cds_end_NF` flags, the mature miRNA arms and the `ncRNA` structure, and the RefSeq, UniProt, TrEMBL, UniParc and UniProt isoform identifiers. They also hold the protein features behind `DOMAINS`, the gene symbol's source and HGNC identifier, the codon table of each sequence region, the `source_versions` of the output headers and the gene phenotype count behind `GENE_PHENO`.

  Without it the `SYMBOL_SOURCE`, `HGNC_ID`, `MANE_SELECT`, `MANE_PLUS_CLINICAL`, `APPRIS`, `miRNA`, `RefSeq`, `SWISSPROT`, `TREMBL`, `UNIPARC`, `UNIPROT_ISOFORM`, `DOMAINS` and `GENE_PHENO` columns are empty on a native cache, and on GRCh37, whose release-87 GFF3 carries no `transcript_support_level` and no `Ensembl_canonical` tag, so are `TSL` and `CANONICAL`. A native cache lists a translation's protein features by feature identifier where the Ensembl cache lists them by analysis in Perl hash order, so `DOMAINS` may name the same features in another order. The JSON `domains` array can also differ: an entry that does not match Ensembl VEP's `(\w+):(\w+)` pattern, such as `Low_complexity_(Seg):seg`, repeats the `db` and `name` of the last matched entry, or is dropped when none precedes it.

  A native cache lists `cds_start_NF` before `cds_end_NF` on every transcript carrying both, whatever order the core dumps or the GTF carried. The Ensembl cache keeps the order its database returned, `cds_end_NF` first on 7 of the 2,946 GRCh38 transcripts that carry both at release 116, and Ensembl VEP prints that order.
- `--gtf`, the Ensembl GTF of the same release, supplies the `cds_start_NF` and `cds_end_NF` transcript attributes (the `FLAGS` output field) when no core dumps are read; the GFF3 does not carry them. Its header also supplies the `assembly` and `genebuild` source versions the output headers print.
- `--genome-fasta`, the primary-assembly genome (decompressed, with a `samtools faidx` index beside it), supplies the CDS sequence behind every codon-level consequence: `missense_variant`, `synonymous_variant`, `stop_gained`, HGVSp and the rest. Without it, coding variants fall back to `coding_sequence_variant`.

`--include-predictions` embeds the SIFT and PolyPhen matrices ([SIFT and PolyPhen predictions](cache-predictions.md)).

A failed build prints one `ERROR: [cache-builder] <message>` line and exits 1, writing nothing outside the output directory. GFF3 lines the builder does not take as written:

| GFF3 line | Outcome |
| --- | --- |
| A seqid that cannot name a directory of the cache (empty, `.`, `..`, or holding a path separator) | refused: `invalid chromosome name '<seqid>' at line <n>` |
| A feature whose start is 0 or whose end precedes its start | refused: `GFF3 coordinates <start>-<end> at line <n> are not 1-based and ordered` |
| A transcript whose exons mix ranked and unranked lines | its exons are ordered by start |

The GRCh38 command is in the [Quick start](getting-started/quick-start.md#2-get-a-transcript-cache). For GRCh37 the same download works: the release-87 files are `Homo_sapiens.GRCh37.87.gff3.gz` and `Homo_sapiens.GRCh37.87.gtf.gz` under `pub/grch37/release-116/`. Local copies can be passed instead:

```bash
target/release/vep-cache-builder \
  --species homo_sapiens --assembly GRCh37 --release 116 \
  --gff3 Homo_sapiens.GRCh37.87.gff3.gz \
  --core-dumps core_dumps_grch37 \
  --genome-fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/homo_sapiens/116_GRCh37
```

The downloads are about 520 MB for GRCh38 and 1.6 GB for GRCh37, whose `phenotype_feature.txt` Ensembl publishes uncompressed (1.5 GB). The tables are read once per build, in about half a minute and 3 GB of memory for GRCh38.

### Variation data

Variation data, which `--check_existing` and the `--af*` flags read, is optional (`--include-variations`). The builder downloads Ensembl's variation VCFs itself from this directory, under `pub/grch37/release-116/` for GRCh37, or reads local copies under the same names from the directory passed with `--variation-vcf-dir`:

```
https://ftp.ensembl.org/pub/release-116/variation/vcf/homo_sapiens/
```

The files are the per-chromosome germline files `homo_sapiens-chr<N>.vcf.gz` and four genome-wide files, each with its `.csi` index. The download covers chromosomes 1-22, X, Y and MT unless `--chromosomes` narrows the list; from a local directory every `*-chr*.vcf.gz` file is read. The genome-wide files are `homo_sapiens_somatic.vcf.gz` (every record of it is flagged `somatic`), `homo_sapiens_phenotype_associated.vcf.gz` and `homo_sapiens_clinically_associated.vcf.gz` (whose records the per-chromosome files already carry at release 116; a record found only there is added) and `1000GENOMES-phase_3.vcf.gz`, the 1000 Genomes frequencies. A genome-wide file the directory or the site lacks is left out with a warning naming what the cache then lacks.

The records are written per chromosome as `variations/<chr>/all_vars.gz` with a `.csi` on the start column, the file Ensembl VEP's own cache carries. `--check_existing` and the frequency flags therefore read a native cache exactly as they read a converted one.

Each record takes the form the variation database holds, as far as the export preserves it. A deletion's alleles have the anchor base the export prepended removed, with `-` for the emptied allele. An insertion the export marks `TSA=insertion` sits two bases after the base its VCF line names (`start` one past `end`), where the Ensembl cache places it. An indel (`TSA=indel`, every allele non-empty) keeps the alleles the export wrote, at the line's `POS`, since the export prepends an anchor only to alleles that do not share a first base. So an `A/AA` record stays `A/AA`, while a `CC/CTTCT` line, whose alleles differ once the shared base goes, is the database's `C/TTCT` one base on.

A COSMIC variant without alleles is written as `COSMIC_MUTATION` on its span, and every HGMD-PUBLIC variant, alleles or not, as `HGMD_MUTATION`, the string the variation database holds for that source. `clin_sig` comes from the record's `CLIN_<term>` flags (a term the export declares no flag for, `likely_risk_allele`, is absent) and `phenotype_or_disease` from `E_Phenotype_or_Disease`.

The Ensembl cache keeps dbSNP's alleles over the whole run an indel can shift along (`CCCCC/CCCC` at 5030241-5030245), which the export trims to `CC/C` at 5030241. Such records are 5 of every 100 of chromosome 21's `indel` records. A native cache holds the trimmed form, which minimises from the right end to the same alleles and position as the Ensembl cache's record. So an input spelled at the left end of the run (the VCF convention, `AC/A` at 5030240) co-locates with it exactly as Ensembl VEP co-locates it.

An input anchored on the run's first base (`CC/C` at 5030241) minimises to the record's second base. The trimmed form also reaches that base and the Ensembl cache's record does not, so the input co-locates with the native record alone. An input spelled at the run's right end co-locates with neither, its start lying outside the one-base window the start-indexed files answer.

A native cache carries no `failed` flag: the export leaves the failed variants out, the `dbSNP_novariation` records among them, and the few it keeps are written unflagged. It carries no `pubmed` or `var_synonyms` (the export carries none) and no gnomAD frequency: Ensembl publishes no gnomAD VCF for release 116, so gnomAD frequencies come from the Ensembl variation cache, by conversion. A native cache built by an earlier vep-rs release holds JSON variation files that nothing reads; build it again.

The 1000 Genomes file exists for GRCh38 only, so a GRCh37 native cache has no frequency column. Where it is read, the `AFR`, `AMR`, `EAS`, `EUR` and `SAS` columns hold each allele's continental frequency as the file gives it, and `AF` the whole panel's. The file states no panel-wide frequency, so `AF` is the allele's count over the panel's 5,008 chromosomes. Each population's count is recovered from its frequency and sample count with two chromosomes per sample (an approximation on X and Y); a population the file gives no value for drops out of both the count and the total.

An allele whose frequency is 0 in every population gets no entry. The release file lists no frequency for an allele absent from its call set, and the export does not tell such an allele from one the release lists with a count of 0. For the latter the Ensembl cache prints `0` (among the 1,106,899 chromosome 21 records the Ensembl cache holds frequencies for). The release's own panel-wide frequency, which the Ensembl cache prints, differs from the recovered count's for 443 of those records, 315 of them by one unit in the last significant digit.

### Checking a cache

`scripts/validation/compare_gtf_tags_to_cache.py` compares the four attributes a GTF yields with those a cache carries and lists every disagreement. `scripts/validation/diff_caches.py` compares two caches' transcripts, variation files and `info.json`.

## Convert an existing Ensembl VEP cache

`scripts/data/storable_to_json.pl` reads an Ensembl VEP Storable cache directory and writes the same JSON layout for every chromosome directory the cache holds: the assembled chromosomes and the patch, alternate-locus, unplaced-scaffold and `LRG_` directories alike, as Ensembl VEP itself lists them. It writes the regulatory build as `regulatory/<chr>/<start>-<end>.json` shards ([Regulatory and motif features](regulatory.md#the-cache)), an `info.json` derived from the cache's `info.txt`, and a copy of the cache's `chr_synonyms.txt`. The script runs inside the `ensemblorg/ensembl-vep` container ([scripts/README.md](../scripts/README.md) gives the invocation), and its output goes straight to `--json_cache`.

Every field the cache dumper keeps on a transcript is carried, flattened or dropped under a rule the parity audit names (`scripts/concordance/cache_parity_audit.pl`):

- Carried: the gene, transcript and protein identifiers (UniProt isoforms included, and on a merged cache the source), the attributes, the exons and introns, and the translation with its sequence edits. Also carried: the translateable sequence and peptide, the codon table, the exon coordinate mapper's pairs, and the SIFT and PolyPhen prediction matrices.
- Flattened: the UTR objects, to their sequence strings, and each protein feature's analysis, to the label Ensembl VEP's `DOMAINS` column prints.
- Dropped: the Perl object references (the slice, the gene object, the exon coordinate mapper's trees), the per-exon and per-intron strand, which is the transcript's, a sequence edit's description and the dump-time flags.

The variation files are copied as they are into `variations/<chr>/`: each chromosome's `all_vars.gz` with its tabix index or, for a cache dumped without the tabix conversion, its `<start>-<end>_var.gz` region files. `all_vars.gz` is the bgzip file the Ensembl VEP dumper writes, with the chromosome first, `.` for an empty cell and the columns of the cache's `variation_cols` line.

vep-rs opens each chromosome's `all_vars.gz` on first use and reads it through its index, or parses a region file whole the first time a variant falls in it. `--check_existing` and the frequency flags therefore read the converted cache as Ensembl VEP reads the original. A cache whose `variations/` directory holds no such file warns once that the co-location columns stay empty. A conversion made with an earlier version of the converter carries no `variations/<chr>/all_vars.gz`: convert the cache again, or copy each chromosome's `all_vars.gz` and its `.csi` from the Ensembl VEP cache into `variations/<chr>/`.

A region file writes an empty cell for a value the record lacks. Ensembl VEP reads an empty `clin_sig_allele` or `clinical_impact` cell inside a line as a value: `CLIN_SIG` then lists no term for the record under the default `--clin_sig_allele 1`, and `--clinvar_somatic_classification` writes its key with an empty value. vep-rs reads the two cells the same way.

A cache written with the Sereal serialiser (`serialiser_type sereal` in its `info.txt`, or `.sereal` shards) is refused with an error naming the serialiser; only Storable caches convert. A run in which any shard fails to read or write stops with an error before `info.json` is written.

The published concordance figures were measured on such a conversion ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)), which holds the transcript set identical between the two engines. The caches those figures were measured on carry 195,232 transcripts (GRCh37) and 505,231 (GRCh38), the canonical chromosomes only. The GRCh37 conversion is of the cache-version-113 Storable that Ensembl VEP release 115 ships for that assembly. It differs from the cache-version-115 build in 310 of those transcripts (the `storable_113_*` columns of `manuscript/data/transcript_set_parity.csv`).

The conversion carries each translation's sequence edits (`seq_edits`: the initial-methionine, selenocysteine and amino-acid-substitution edits Ensembl VEP applies to the reference peptide at annotation time), and vep-rs applies them the same way. A cache converted without that key, the published conversions among them, still loads. On such a cache, a variant at codon 1 of a transcript with an alternative initiation codon is called from the cached peptide rather than from the edit. Converting the cache again restores Ensembl VEP's call. A natively built cache carries no edit list either: its peptide is the protein FASTA's sequence, with the edits already applied, and the builder writes none.

## Naming the cache

Ensembl VEP's cache flags select the directory as Ensembl VEP does ([Cache, species and runtime](cli/cache-species-runtime.md)). The cache root is `--dir_cache`, else a directory given to `--cache`, else `--dir` (default `$HOME/.vep`). Under it the species directory is `--species` (default `homo_sapiens`, the Latin name only), suffixed `_refseq` or `_merged` under those flags. Under that sits the directory named `<cache_version>_<assembly>`, where `--cache_version` defaults to the VEP release.

When one directory of that version exists, it is the cache and its assembly is read from its name; when several exist, `--assembly` chooses. An `--assembly` that disagrees with the directory name, or a directory name or `--assembly` that disagrees with the cache's own `info.json`, is an error.

Every failure to resolve a cache exits 1 before the input is opened with one `ERROR: [cache] ...` line, in Ensembl VEP's words where it has them. Those words are `Cache directory <dir> not found`, `No cache found for homo_sapiens, version 116` and `Multiple assemblies found for cache version 116 ...`, the assemblies listed in name order. The metadata mismatch names `info.json`, the file a JSON cache carries where a Storable cache has `info.txt`.

`--json_cache <dir>` (alias `--full_cache_dir`) names the cache directory itself instead, whatever `--dir_cache`, `--species` and `--cache_version` say; only `--assembly` is checked against its `info.json`. Either way the directory must hold `transcripts/` and `info.json`, so a Storable cache or a mistyped path fails rather than annotating nothing. `transcripts/` must also hold a shard under some chromosome directory: a cache with none is refused, where Ensembl VEP writes every record as `intergenic_variant` ([A cache holding no transcript shard](divergences/damaged-cache-shard.md#a-cache-holding-no-transcript-shard)).

`--assembly` is compared exactly and case-sensitively, as Ensembl VEP compares it, with one addition: `hg19` and `hg38` are read as `GRCh37` and `GRCh38`. When `--fasta` is absent, the first `.fa` or `.fa.gz` by name inside the cache directory that has its indexes beside it is used as the reference, unless `--no_fasta` is given. An unindexed file there is not read (Ensembl VEP builds a missing index itself). `--show_cache_info` resolves the cache by the same flags, refuses a named input file that does not exist as Ensembl VEP does, prints the cache's source versions as `<key>\t<value>` lines sorted by key, and exits without reading the input.

```bash
# Ensembl VEP's layout: $HOME/.vep/homo_sapiens/116_GRCh38
vep -i input.vcf -o output.txt --offline --cache \
  --dir_cache $HOME/.vep --cache_version 116 --assembly GRCh38

# The directory itself
vep -i input.vcf -o output.txt --offline --json_cache /path/to/116_GRCh38
```

## Reference FASTA

`--fasta` takes a plain FASTA with its `samtools faidx` `.fai`, or a `bgzip`-compressed FASTA with the `.fai` and `.gzi` `samtools faidx` writes for it, the shape of Ensembl's `dna_index/` files. No index is built (Ensembl VEP builds one). Three files are refused at open: a gzip file that is not BGZF, whose `.fai` offsets need the whole file inflated; a `.gzi` that does not match the file's blocks, recompressed and not re-indexed; an `.fai` naming more bases than the file holds. The table adds five flag conditions.

| Condition | Line printed |
| --- | --- |
| A gzip file that is not BGZF | `ERROR: [fasta] <file> is gzip but not BGZF: recompress with bgzip and index with samtools faidx` |
| A `.gzi` that does not describe the file's blocks | `ERROR: [fasta] BGZF index <file>.gzi does not describe <file> (<the first disagreement>): index the file again with samtools faidx` |
| An `.fai` naming bases past the file's end | `ERROR: [fasta] FASTA index names <n> bases beyond the end of <file>` |
| `--fasta 0` or `--fasta ''`, a value Ensembl VEP reads as the flag unset ([compatibility page](cli/compatibility.md#where-configpm-is-not-followed)) | `ERROR: [fasta] Specified FASTA file/directory <value> not found` |
| `--hgvs` or `--hgvsg` with no FASTA in reach | `ERROR: [fasta] Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file (see --fasta)` |
| `--use_transcript_ref` without a FASTA | `ERROR: [fasta] Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref` |
| `--everything` without a FASTA | `INFO: Disabling --hgvs; using --offline and no FASTA file found` |
| `--fasta_dir` holding no compatible file | `WARNING: No compatible FASTA files found in <dir>` |

A `--fasta` directory is read as Bio::DB::Fasta reads one: every indexed `.fa` in it, plus the `.fa.gz` vep-rs adds, as one set of sequences; two files naming one sequence are refused (Ensembl VEP with Bio::DB::HTS takes a file only).

On GRCh37 and GRCh38 a read of `Y` inside a pseudoautosomal region is served from `X` at the region's offset, as Ensembl VEP serves it: Ensembl's primary-assembly FASTA carries `N` there. The redirect follows the FASTA's own name for the sequence, so a file naming it `chrY` is read as written, as Ensembl VEP reads it.

The redirect is armed by the assembly the cache resolves to, its `info.json` first, then its directory name, then `--assembly`; Ensembl VEP arms it from its assembly parameter (`BaseVEP.pm` line 525, `Utils/FastaSequence.pm` line 254). That parameter is set by `--assembly`, `--grch37` (`Config.pm` lines 399-403) or the layout scan (`CacheDir.pm` lines 303-356), which runs for `--dir_cache` and `--dir`, not under `--full_cache_dir`. Under `--offline --full_cache_dir` (`--json_cache` here) without `--assembly` or `--grch37` it stays unset, so Ensembl VEP reads a `Y` span inside the region as written and vep-rs reads `X`.

`--fasta_dir <dir>` chooses, as Ensembl VEP does, the `.fa` or `.fa.gz` whose name carries the species and assembly (case-insensitively), `toplevel` before any other name and a compressed file before a plain one, in place of `--fasta` and of the cache directory's FASTA. A directory holding no such file leaves those in force and warns (the table's line); one that cannot be opened is passed over silently, as Ensembl VEP does.

Two files tied on both preferences are ordered by name, the first taken, where Ensembl VEP's sort on the two preferences alone keeps the directory listing's order (`BaseVEP.pm` lines 500-509), so on such a directory the two engines can read different files. Here the filter reads the same assembly as the redirect, where Ensembl VEP reads its assembly parameter alone (line 497). With that parameter unset, Ensembl VEP's filter passes every file naming the species (line 508), and its preferences and the listing's order choose among every assembly's files.

The table's `ERROR:` lines exit 1; `--use_transcript_ref` needs the FASTA for a non-coding transcript's cDNA. `--everything` without a FASTA instead prints its INFO line on standard error (withheld under `--quiet`) and drops the HGVS notations it switched on, as Ensembl VEP does; an explicit `--hgvsg` beside it is still refused.

Under `--lookup_ref` a contig the FASTA lacks skips the record at every position with the could-not-fetch warning, where Ensembl VEP's one-base placeholder slice (`BaseVEP.pm` lines 424-434) gives a one-base record at position 1 the reference `N`. Under `--check_ref` the same contig skips the record at every position with the sub-slice warning and then the mismatch warning; an insertion, naming no reference base, passes under either flag in both engines (`Parser.pm` lines 600-602 and 636). So Ensembl VEP annotates a one-base record there whose reference is written `N`, and refuses one written otherwise with the mismatch warning alone (`Parser.pm` lines 614-631).

Under `--hgvs` the same contig gives these notations; neither engine warns:

| Record on the contig under `--hgvs` | Ensembl VEP, whose transcript slice reads `N` per base there (`Utils/FastaSequence.pm` lines 597-601) | vep-rs, which reads no reference base there |
| --- | --- | --- |
| A substitution | `c.764N>T`, the `N` standing for the genome's `G` | no HGVSc |
| A deletion | `c.632del`, as over the genome | no HGVSc |
| An insertion | `c.307_308insC` | `c.307_308insC` |
| HGVSp | as over the genome | as over the genome |

### A `Y` read on a FASTA without `X`

With the pseudoautosomal redirect in force and no `X` in the FASTA, a `Y` read's bases inside the region are `N`, as Ensembl VEP reads them (`Utils/FastaSequence.pm` lines 597-601). The bases outside the region come from `Y`, so the read has the length asked for. Ensembl VEP returns the region's `N`s alone for a read that straddles the region's edge, two bases for a four-base request. Neither engine warns on the read: Ensembl VEP's `FASTA does not contain sequence for X` (line 599) is raised at warning level 2000, above the verbosity of 1999 a run sets (`BaseRunner.pm` line 491).

## Chromosome names

An input may name a chromosome as the cache does (`21`), with a `chr` prefix (`chr21`), or by an accession (`NC_000021.9`, `CM000683.2`). Each name is resolved to the cache's own before any lookup, as Ensembl VEP's `get_source_chr_name` resolves it, in this order:

1. the name itself, when the cache holds it;
2. else a synonym of it the cache holds;
3. else the name with a leading `chr` removed or added, when the cache holds that;
4. else a bare `M` as the cache's `MT`, under whatever name the cache holds it.

A name that resolves to nothing skips the record with Ensembl VEP's warning (`Chromosome NC_000022.11 not found in annotation sources or synonyms`).

The output keeps the input's spelling in `Location`, in the JSON `seq_region_name` and in a VCF input's `CHROM` column, which is written as read. A bare `M` the cache knows as `MT` but not as `M` or `chrM` prints `MT` in `Location` and `seq_region_name`, as Ensembl VEP renames it; [A bare `M` under a cache whose synonyms name `chrM`](divergences/bare-m-chromosome.md) describes where the two engines then differ.

The synonyms come from the cache's `chr_synonyms.txt`, the file the Ensembl VEP cache ships and the converter copies; the native builder writes it from the core `seq_region_synonym` table with `--core-dumps`. `--synonyms <file>` replaces the cache's file instead of extending it. Each line holds the whitespace-separated names of one chromosome, the first the reference name and the rest its synonyms, in both directions.

A file that cannot be read is `ERROR: [config] Could not read synonyms file <file>: <reason>` for a `--synonyms` file and `ERROR: [cache] ...` for the cache's own (a directory included, which Ensembl VEP opens as an empty table, `BaseVEP.pm` line 644). A file holding a byte that is not UTF-8 is refused the same way, the reason `stream did not contain valid UTF-8`, where Ensembl VEP reads the file as bytes.

The `--synonyms` value is a file name as given (`-`, `0`, an empty value, surrounding whitespace, a leading `>` or `|`), where Ensembl VEP's two-argument `open` reads `-` as standard input, trims the whitespace and takes a false value as unset ([compatibility page](cli/compatibility.md#where-configpm-is-not-followed)).

A cache without the file (a native build without `--core-dumps`, or a conversion by an earlier converter) resolves no synonym, so there `chrM` and an accession are skipped, as Ensembl VEP does on such a cache. Convert the cache again, or pass `--synonyms`.

Three cases differ from Ensembl VEP, which checks a record against the cache's names widened by their synonyms, a set it widens again in place at every miss (`Parser.pm` line 782), and then looks the unwidened name up one hop. Each passes its check there and annotates as intergenic; in vep-rs a bare `M` resolves to `MT` and every other such name is skipped as not found, whatever record came before it.

The first is a bare `M` on a cache whose file links `chrM` but not `M`. The second is a name that reaches the cache only through a synonym of its `chr`-stripped or `chr`-prefixed spelling (`chrNC_000021.9`). The third is a name two hops from a cache chromosome under a `--synonyms` table whose lines chain or open with whitespace (`21 ALPHA` then `ALPHA BETA`: `BETA`), admitted once an earlier record missed.

A reference FASTA whose records are named by accession is read under the cache's names through the same table, so `--fasta` may name a RefSeq or GenBank assembly file beside an Ensembl cache.

Beyond the table, a record is reached by its name with the prefix `chr` added or removed and, for a record named `MT`, `M`, `chrM` or `chrMT`, by any of the four names, matched in any letter case: aliases vep-rs builds in. Ensembl VEP tries the exact name, a synonym and the `chr` variants alone (`BaseVEP.pm` lines 413 and 689-738, `Utils/FastaSequence.pm` lines 380-389), so a cache's `MT` over a FASTA record named `chrM`, with no synonym linking the two, reaches no sequence there. That chromosome gets the one-base placeholder slice described above (`BaseVEP.pm` lines 425-426), so under `--check_ref` a record such as `MT 3308 T G` is skipped with the sub-slice and mismatch warnings, where vep-rs annotates it.

## What `info.json` carries

`info.json` carries `species` and `assembly` in every cache the converter and the builder write; the resolver compares `--assembly` with `assembly`. Its other keys:

- `sift` and `polyphen`, markers written only for a tool whose matrices the cache holds; a converted cache keeps a marker only beside the tool's `source_versions` entry ([SIFT and PolyPhen predictions](cache-predictions.md)).
- `source_versions`: on a native cache the core `meta` table's `assembly`, `gencode` and `genebuild`, or the GTF header's. The dbSNP, ClinVar, COSMIC and gnomAD versions an Ensembl cache's `info.txt` lists come from the variation database, which the builder does not read, so a native cache's headers print none.
- `variation_cols` and `var_type tabix`, only when the cache holds variation files; `variation_cols` names their columns.
- `regulatory`, written only when the cache holds the `regulatory/<chr>/` tree.
- `contig_order`: the cache's chromosome directories as autosomes by number, then X, Y and MT, then every other contig by name.
