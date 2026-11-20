# Cache setup

The JSON transcript cache vep-rs annotates from is built with `vep-cache-builder` or converted from an Ensembl VEP cache ([Cache conversion](cache-conversion.md)).

[Cache checksum manifest](cache-checksums.md) covers integrity checks, [SIFT and PolyPhen predictions](cache-predictions.md) the prediction matrices.

The cache is a directory holding `info.json` and `transcripts/`, the release directory itself (`$HOME/.vep/homo_sapiens/116_GRCh38`). `--json_cache` names it, or Ensembl VEP's cache flags find it ([Naming the cache](#naming-the-cache)). [Other species](cache-species.md) covers the genomes beside human, the four compared with Ensembl VEP among them.

## Build a cache

`vep-cache-builder` builds a cache natively, with no Perl dependency, and one command installs a release. `--auto` fetches every input from the Ensembl site, checks each against its directory's `CHECKSUMS` file, builds the cache as `<dir_cache>/<species>/<release>_<assembly>` and places the indexed genome FASTA inside it:

```bash
target/release/vep-cache-builder --auto homo_sapiens GRCh38 --release 116 \
  --dir_cache $HOME/.vep
```

The inputs are the release's GFF3, GTF, protein FASTA, core and variation table dumps and genome FASTA: the `dna_index/` file with its published `.fai` and `.gzi`, else the `dna/` file, recompressed to bgzip and indexed by the builder. A file whose BSD `sum` differs from its `CHECKSUMS` line stops the installation, naming the file. The downloads land under `<dir_cache>/tmp/` and go once the cache is in place; an existing release directory is never written into.

`--include-predictions` adds the SIFT and PolyPhen matrices ([SIFT and PolyPhen predictions](cache-predictions.md)), `--include-variations` the variation data below. `--site` names another server or a local directory holding the site's `pub/release-<N>/` tree. For GRCh37 the same command reads `pub/grch37/release-<N>/`, whose gene set is the release-87 one Ensembl publishes from release 88 onward. The cache is built from the release files, never from Ensembl VEP's prebuilt cache archive, whose Storable format `vep` does not read.

### A build from local copies

Local copies build the same cache into `--output-dir`, one input per flag. The model is the release's Ensembl GFF3, downloaded or passed with `--gff3`, or any local GFF3 or GTF ([A GTF or GFF3 as the transcript model](cache-from-gtf-or-gff3.md)). The peptides come from the release's protein FASTA, downloaded or passed with `--fasta`.

A GRCh37 build from local copies:

```bash
target/release/vep-cache-builder \
  --species homo_sapiens --assembly GRCh37 --release 116 \
  --gff3 Homo_sapiens.GRCh37.87.gff3.gz \
  --fasta Homo_sapiens.GRCh37.pep.all.fa.gz \
  --core-dumps core_dumps_grch37 \
  --genome-fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/homo_sapiens/116_GRCh37
```

The downloads are about 1.6 GB for GRCh38 and 2.6 GB for GRCh37 (its `phenotype_feature.txt` is published uncompressed). The tables are read in about half a minute and 3 GB of memory.

### What each input supplies

The Ensembl GFF3 alone gives each transcript:

- its exons
- its CDS and ENSP identifier
- its biotype
- its CCDS identifier and transcript support level
- its GENCODE basic and primary, canonical and MANE memberships
- the gene's symbol and version

With the release's protein FASTA or core dumps, the builder drops `artifact` and readthrough transcripts, as Ensembl VEP's dumper does. Three more inputs decide the rest:

- `--core-dumps` is a directory of the release's core and variation MySQL table dumps (a missing table is downloaded into it). It supplies what the GFF3 lacks: the columns below, each sequence region's codon table ([below](#the-codon-table-of-each-sequence-region)), the headers' `source_versions` and the translation's sequence edits.

  Without it the `SYMBOL_SOURCE`, `HGNC_ID`, `MANE_SELECT`, `MANE_PLUS_CLINICAL`, `APPRIS`, `miRNA`, `RefSeq`, `SWISSPROT`, `TREMBL`, `UNIPARC`, `UNIPROT_ISOFORM`, `DOMAINS` and `GENE_PHENO` columns are empty on a native cache, and on GRCh37, whose release-87 GFF3 carries no `transcript_support_level` and no `Ensembl_canonical` tag, so are `TSL` and `CANONICAL`. `DOMAINS` may list a translation's protein features in another order. A native cache orders them by feature identifier, Ensembl VEP by [the region-file copy its run reads](divergences/domains-order-region-file-copy.md). In the JSON `domains` array an entry outside Ensembl VEP's `(\w+):(\w+)` pattern (`Low_complexity_(Seg):seg`) repeats the last matched entry's `db` and `name`. It is dropped when none precedes it.

  A native cache lists `cds_start_NF` before `cds_end_NF` on every transcript carrying both, whatever order the core dumps or the GTF carried. The Ensembl cache keeps the order its database returned, `cds_end_NF` first on 7 of the 2,946 GRCh38 transcripts that carry both at release 116, and Ensembl VEP prints that order.
- `--gtf`, beside the Ensembl GFF3, is the Ensembl GTF of the same release. It supplies the `cds_start_NF` and `cds_end_NF` transcript attributes (the `FLAGS` output field) when no core dumps are read; the GFF3 does not carry them. Its header also supplies the `assembly` and `genebuild` source versions the output headers print. Without `--gff3`, the GTF is the transcript model itself ([A GTF or GFF3 as the transcript model](cache-from-gtf-or-gff3.md)).
- `--genome-fasta` is the genome with its `samtools faidx` indexes, plain or bgzip-compressed. It supplies the CDS sequence behind every codon-level consequence: `missense_variant`, `synonymous_variant`, `stop_gained`, HGVSp and the rest. Without it, coding variants fall back to `coding_sequence_variant`.

  Without `--fasta` and `--include-predictions` it also supplies the peptides. Each is the translation of the transcript's coding sequence under the region's codon table, as Ensembl VEP translates a transcript read from a GFF3 or GTF. The stop codon is dropped and a start codon of the table is read as `M`. The protein FASTA's peptide carries the database's sequence edits instead, so a selenocysteine reads `U` there and `*` in the translation. Pass `--fasta` to keep those residues.

### The codon table of each sequence region

Each transcript's `codon_table` is the genetic code its consequences and peptide are read under. It is its sequence region's `codon_table` attribute in the core dumps, as Ensembl VEP's dumper reads it from the slice. A region without the attribute is table 1, the nuclear code, as the dumper writes it.

Without `--core-dumps` nothing states an organelle genome's table, so `--codon-table <region>=<table>` names it, repeatable or comma-joined. The region is spelled as the GFF3 or GTF spells it (`chrM` and `chrMT` read as `MT`). The table is one of the genetic codes `Bio::Tools::CodonTable` ships:

```bash
--codon-table MT=2                   # a vertebrate's mitochondrion, no dumps
--codon-table Mito=3                 # yeast
--codon-table mitochondrion_genome=5 # the fly
--codon-table Pt=11,Mt=1             # the Arabidopsis plastid and mitochondrion
```

A region named like an organelle genome, in any letter case, whose table neither the dumps nor the flag states stops the build before any shard is written. The line is `ERROR: [cache-builder] codon table of sequence region <name> unknown: pass --core-dumps or --codon-table <name>=<n>`. The names are:

- `MT`
- `chrM`
- `chrMT`
- `Mito`
- `MtDNA`
- `mitochondrion_genome`
- `Mt`
- `Pt`
- `ChrC`
- `ChrM`
- `chloroplast`
- `plastid`

A region the dumps and the flag give two different tables stops the build too. So does a table the module lacks (7, 8, 17 to 20 and anything past 24).

The release 116 core dumps carry the attribute for the human and the fly mitochondrion (`MT` 2, `mitochondrion_genome` 5). So a build with `--core-dumps` needs no flag for them. Where the dumps carry no attribute for an organelle-named region, the build stops the same way and the line names the flag alone.

### GFF3 or GTF lines the builder does not take as written

A failed build prints one `ERROR: [cache-builder] <message>` line and exits 1, writing nothing outside the output directory. GFF3 or GTF lines the builder does not take as written:

| GFF3 or GTF line | Outcome |
| --- | --- |
| A file from which no transcript is read | refused before the build writes anything: `No transcript was read from <file>` |
| A file whose every transcript is `artifact` or readthrough, on a build that reads a protein FASTA or core dumps | refused before the build writes anything: `No transcript was kept from <file>: <n> read, <a> set aside as artifact and <r> as readthrough, as Ensembl VEP's cache holds neither`; a build from the genome alone keeps both kinds |
| A line that is not UTF-8 | refused: `<file> line <n>: GFF3 line is not UTF-8` (`GTF` for a GTF; a protein FASTA's `FASTA line is not UTF-8`) |
| A line that is no comment and ends before its fourth column, the start | refused: `<file> line <n>: GFF3 line ends at column <k>, before its start` (`GTF` for a GTF), as the index build fails such a line; a blank line is skipped, a GFF3's `##FASTA` directive is a comment, and nothing is read from the first line that begins with `>` and has a sequence behind it, the next line that is neither blank nor a comment having fewer than four columns, or a start or end the index build fails (a line longer than 16777216 bytes, only where no tab stands in it, its 16777217th byte aside), to the file's end, a line the log names; a `>` line before a line of four or more columns whose start and end the index build parses, or with nothing behind it, is read as any other line |
| A record inside a GFF3's sequence section, nine or more columns whose fourth and fifth are whole numbers, neither a comment nor led by `>` | refused: `<file> line <n>: GFF3 annotation after the sequence section, which opens at line <m>`; the section ends the annotation, so the record would go unread |
| A byte-order mark at the start of a line | taken off before the line is read, in a GFF3 or GTF model, a GTF read for its tags, a regulatory file and the protein FASTA: a file with marks builds the cache of the same file without them |
| A line that opens with a tab | read at the columns the file writes, its seqid empty, as the index build reads it: refused for the type that stands where a start is read, `<file> line <n>: GFF3 start "<type>" is not a number`, or, on a line the model reads, for the empty seqid; a line of a skipped type whose own start and end the index build parses is skipped; white space and a tab before a record refuse it for the type alike, the record's type kept or skipped |
| A seqid that cannot name a directory of the cache (empty, `.`, `..`, or holding a path separator), on a line the model reads | refused: `<file> line <n>: invalid chromosome name '<seqid>'` |
| A feature whose start or end is not a number, whose start is 0 or whose end precedes its start, on a line the model reads; a start or end no number leads (`abc`, a lone sign, an empty end before a score of `.`) or an end below zero (`-9`), on any line of the model's file | refused: `<file> line <n>: GFF3 start "<value>" is not a number`, `<file> line <n>: GFF3 coordinates <start>-<end> are not 1-based and ordered` (`GTF` for a GTF); a line of a type the model skips is read no further than its type, its identifier, whether a number leads its start and end and whether its end is below zero, as Ensembl VEP reads only the start and the type of one, and a line of four to eight columns no further than its start and end; an empty or blank column is led by the number of the column after it, as `tabix` reads one; Ensembl VEP builds its model without a kept line its `tabix` could not place, or dies, a documented intended divergence, [`gxf_malformed_coordinate_line_refused`](divergences/gxf-malformed-coordinate-line-refused.md) |
| A gene, transcript or CDS line whose `version` (a transcript line's `transcript_version` and, on a build of the cache's content, a GTF CDS line's `protein_version` too) is not a whole number, or is one above 4294967295, the largest the cache holds | refused: `<file> line <n>: version "<text>" is not a whole number`, or `<file> line <n>: version <digits> is above 4294967295, the largest version the cache holds`; Ensembl VEP prints the text as written, [`gff_version_not_a_number`](divergences/gff-version-not-a-number.md); a protein FASTA header whose identifier ends in `.` and digits above the limit is refused the same way, and one ending in `.` followed by anything but digits (`.cds1`, `.p1`) states no version, as does one whose identifier is the CDS line's `protein_id` whole |
| A GTF line repeating an attribute key read as one value (`gene_id`, `transcript_id`, `gene_name`, `transcript_biotype`, a version, `protein_id` and the rest; `tag` is a list) | refused: `<file> line <n>: GTF attribute <key> is repeated; Ensembl VEP reads a repeated key as a list`; Ensembl VEP prints a symbol, biotype or protein identifier as the values joined by commas, a version (and a protein identifier under `--protein_version`) as the list's address, no `TSL` for a repeated support level, and loses a record whose identifier it was, [`gtf_repeated_attribute_key`](divergences/gtf-repeated-attribute-key.md) |
| A GTF value whose `%XX` escapes spell bytes that are not UTF-8 (`"%C3"`) | refused: `<file> line <n>: GTF attribute <key> "<value>" decodes to text that is not UTF-8`; Ensembl VEP carries the bytes |
| A GFF3 value whose `%XX` escapes spell bytes that are not UTF-8 (`Name=caf%FF`), on a gene, transcript, exon or CDS line the model reads | refused: `<file> line <n>: GFF3 attribute <key> "<value>" decodes to text that is not UTF-8`; a line of a type the model skips refuses nothing, as Ensembl VEP reads no attribute of one; Ensembl VEP carries the bytes |
| A GFF3 source column whose `%XX` escapes spell bytes that are not UTF-8, on a gene or transcript line | refused: `<file> line <n>: GFF3 source "<value>" decodes to text that is not UTF-8`; a type column that does so matches no type the model reads, so its line is skipped, as it is on Ensembl VEP |
| A GTF value holding `%XX` escapes or quotes inside the outer pair (`"A%2CB"`, `"a"b"`) | read as Ensembl VEP reads it: every quote removed, every escape decoded, a tab or line feed included (`A,B`, `ab`); `"%30"` is the absent `0` |
| A GFF3 value holding a `%XX` escape (`Name=A%2CB`) | decoded (`A,B`); Ensembl VEP's `--gff` source prints the escape, [`gff_attribute_escapes_read_raw`](divergences/gff-attribute-escapes-read-raw.md). A value with a raw second `=` ends before it on both engines |
| A GFF3 pair written after a space (`; Name=x`) or without `=` (`Name`) | its key is ` Name`, so `Name` is unset, as it is to Ensembl VEP, which trims nothing; a pair without `=` unsets its key on both engines |
| A translation that starts or ends in no exon of its transcript, on either strand: its start is the 5' end of the first CDS line in transcript order and its end the 3' end of the last, a GTF's `stop_codon` lines counted with its CDS lines | refused before the build writes anything: `<file> line <n>: CDS <start>-<end> of transcript <id> starts the translation at <base>, in no exon of the transcript` or `... ends the translation at <base>, in no exon of the transcript`; the line and the span are those of the CDS line that sets the base; a CDS line that leaves its exon between the two ends is read as written; Ensembl VEP builds without the line or on its span as written, a documented intended divergence, [`gxf_malformed_coordinate_line_refused`](divergences/gxf-cds-line-outside-its-exons.md) |
| A transcript whose exon `rank` values run against its strand, its coding sequence opening before the exon ranked first | refused before the build writes anything: `<file> line <n>: the coding sequence of transcript <id>, <start>-<end>, starts before <start>-<end>, the first exon in transcript order that it overlaps` (`ends after` on the reverse strand); the line is the CDS line's that opens the coding sequence and the first span the whole coding sequence's |
| A transcript whose translation's end would be counted outside its exon: exon `rank` values that run against position at the translation's end, or an exon nested in the last coding exon and coming after it | refused before the build writes anything: `<file> line <n>: CDS <start>-<end> of transcript <id> ends the translation at base <offset> of the <length>-base exon <start>-<end>, the last in transcript order that the coding sequence overlaps`; the line and the span are those of the CDS line that sets the translation's end |
| A transcript whose exons mix ranked and unranked lines | its exons are ordered by start |
| A line of a type outside the [transcript-model page](cache-from-gtf-or-gff3.md)'s list (`five_prime_UTR`, `biological_region`, `start_codon`) | skipped, unless no number leads its start or end or its end is below zero (the row of a start or end that is not a number); the skipped types are logged once with their counts |
| A record whose parent list (`Parent`, or a GTF's `gene_id` or `transcript_id`) opens with a parent in no line that was read | placed as its other parents place it, else dropped with its children; a warning names the parents |
| A record whose parent list opens with `0` and names nothing set (`Parent=0,0`; a GTF transcript's `gene_id "0"`) | dropped, without a warning |
| A record whose parent list opens with an empty field (`Parent=,g1`) | kept as a record of its own, or as the transcript of a gene a later field names ([`gff_parent_list_with_an_empty_field`](divergences/gff-parent-list-with-an-empty-field.md)) |
| A transcript whose biotype no attribute or type rule names | dropped, with a warning naming it |
| A `protein_coding` transcript without an exon line | dropped, with a warning naming it |
| A transcript with no exon or CDS line | dropped |

### A GTF or another GFF3 as the transcript model

`vep --gff <file>` and `vep --gtf <file>`, which Ensembl VEP reads as an annotation source in place of its cache, have no in-process counterpart. vep-rs reads transcripts from its cache alone. The same annotation is a `vep-cache-builder` run on the file and the genome, then a `vep` run against that cache. [A GTF or GFF3 as the transcript model](cache-from-gtf-or-gff3.md) gives the commands and the cache's contents.

### Variation data

Variation data, which `--check_existing` and the `--af*` flags read, is optional (`--include-variations`). The builder downloads Ensembl's variation VCFs from this directory (`pub/grch37/release-116/` for GRCh37), or reads local copies under the same names from `--variation-vcf-dir`:

```
https://ftp.ensembl.org/pub/release-116/variation/vcf/homo_sapiens/
```

The files are the per-chromosome germline files `homo_sapiens-chr<N>.vcf.gz` and four genome-wide files, each with its `.csi`. The download covers chromosomes 1-22, X, Y and MT under `homo_sapiens` unless `--chromosomes` narrows it. For another species it covers its GFF3's `##sequence-region` names. A GTF declares none, so a GTF-only species names its chromosomes with `--chromosomes` or the build stops. From a local directory every `*-chr*.vcf.gz` is read; a `##fileformat` outside major version 4 stops the build.

The four genome-wide files:

| Genome-wide file | Contents |
| --- | --- |
| `homo_sapiens_somatic.vcf.gz` | every record flagged `somatic` |
| `homo_sapiens_phenotype_associated.vcf.gz` | records the per-chromosome files already carry at release 116; a record found only there is added |
| `homo_sapiens_clinically_associated.vcf.gz` | records the per-chromosome files already carry at release 116; a record found only there is added |
| `1000GENOMES-phase_3.vcf.gz` | the 1000 Genomes frequencies |

A genome-wide file the directory or the site lacks is left out with a warning naming what the cache then lacks.

The records are written per chromosome as `variations/<chr>/all_vars.gz` with a `.csi` on the start column, the file Ensembl VEP's own cache carries. So `--check_existing` and the frequency flags read a native cache exactly as a converted one.

#### The form of a record

Each record takes the form the variation database holds, as far as the export (Ensembl's variation VCFs) preserves it. A deletion's alleles have the anchor base the export prepended removed, with `-` for the emptied allele. An insertion the export marks `TSA=insertion` sits two bases after the base its VCF line names (`start` one past `end`), where the Ensembl cache places it.

An indel (`TSA=indel`, every allele non-empty) keeps the alleles the export wrote, at the line's `POS`. The export prepends an anchor only to alleles that do not share a first base. So an `A/AA` record stays `A/AA`, while a `CC/CTTCT` line, whose alleles differ once the shared base goes, is the database's `C/TTCT` one base on.

A COSMIC variant without alleles is written as `COSMIC_MUTATION` on its span. Every HGMD-PUBLIC variant, alleles or not, is written as `HGMD_MUTATION`, the string the variation database holds for that source. `clin_sig` comes from the record's `CLIN_<term>` flags (a term the export declares no flag for, `likely_risk_allele`, is absent) and `phenotype_or_disease` from `E_Phenotype_or_Disease`.

#### An indel that can shift along a run

The Ensembl cache keeps dbSNP's alleles over the whole run an indel can shift along (`CCCCC/CCCC` at 5030241-5030245), which the export trims to `CC/C` at 5030241. Such records are 5 of every 100 of chromosome 21's `indel` records. A native cache holds the trimmed form, which minimises from the right end to the same alleles and position as the Ensembl cache's record. So an input spelled at the left end of the run (the VCF convention, `AC/A` at 5030240) co-locates with it exactly as Ensembl VEP co-locates it.

An input anchored on the run's first base (`CC/C` at 5030241) minimises to the record's second base. The trimmed form also reaches that base and the Ensembl cache's record does not, so the input co-locates with the native record alone. An input spelled at the run's right end co-locates with neither, its start lying outside the one-base window the start-indexed files answer.

#### What a native cache lacks

A native cache carries no `failed` flag. The export leaves the failed variants out, the `dbSNP_novariation` records among them, and the few it keeps are written unflagged. The cache carries no `pubmed` or `var_synonyms` (the export carries none) and no gnomAD frequency. Ensembl publishes no gnomAD VCF for release 116, so gnomAD frequencies come from the Ensembl variation cache, by conversion. A native cache built by a vep-rs release before 0.5.0 holds JSON variation files that nothing reads; build it again.

#### The 1000 Genomes frequencies

The 1000 Genomes file exists for GRCh38 only, so a GRCh37 native cache has no frequency column. Where it is read, the `AFR`, `AMR`, `EAS`, `EUR` and `SAS` columns hold each allele's continental frequency as the file gives it, and `AF` the whole panel's. The file states no panel-wide frequency, so `AF` is the allele's count over the panel's 5,008 chromosomes. Each population's count is recovered from its frequency and sample count with two chromosomes per sample (an approximation on X and Y). A population the file gives no value for drops out of both the count and the total.

An allele whose frequency is 0 in every population gets no entry. The release file lists no frequency for an allele absent from its call set, and the export does not tell such an allele from one the release lists with a count of 0. For the latter the Ensembl cache prints `0` (among the 1,106,899 chromosome 21 records the Ensembl cache holds frequencies for). The release's own panel-wide frequency, which the Ensembl cache prints, differs from the recovered count's for 443 of those records, 315 of them by one unit in the last significant digit.

### Regulatory data

`--regulatory-gff`, `--emars-gff` and `--motif-gff` bake the Ensembl regulatory annotation files into the cache's regulatory tree ([Regulatory features from GFF files](regulatory.md#regulatory-features-from-gff-files)); `--extended-promoters` widens the promoters:

```bash
target/release/vep-cache-builder --release 116 --assembly GRCh38 \
  --output-dir cache --gff3 Homo_sapiens.GRCh38.116.gff3.gz \
  --genome-fasta genome.fa.gz --codon-table MT=2 \
  --regulatory-gff Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz \
  --motif-gff Homo_sapiens.GRCh38.motif_features.v116.gff3.gz \
  --emars-gff Homo_sapiens.GRCh38.EMARs.v116.gff.gz
```

### Checking a cache

`scripts/validation/compare_gtf_tags_to_cache.py` compares the four attributes a GTF yields with a cache's and lists every disagreement. `scripts/validation/diff_caches.py` compares two caches' transcripts, variation files and `info.json`.

## Convert an existing Ensembl VEP cache

[Cache conversion](cache-conversion.md) describes the converter, `scripts/data/storable_to_json.pl`, which converts an Ensembl VEP Storable cache to this layout, what it carries and the caches the published figures were measured on. The corpora under `tests/golden/116/` compare both engines on converted mouse, fruit fly, Arabidopsis and yeast caches ([Other species](cache-species.md)).

### RefSeq and merged caches

Ensembl also ships a RefSeq and a merged cache. Each converts with that converter into `homo_sapiens_refseq/116_GRCh38` or `homo_sapiens_merged/116_GRCh38` under the cache root, where `--refseq` or `--merged` finds it. The converter reads the set from the species directory's suffix, so it is given the directory inside that layout:

```bash
tar -xzf homo_sapiens_refseq_vep_116_GRCh38.tar.gz
docker run --rm -v $(pwd):/input:ro -v $HOME/.vep:/output \
  -v $(pwd)/scripts/data:/scripts:ro \
  ensemblorg/ensembl-vep:release_116.2 \
  perl /scripts/storable_to_json.pl \
  /input/homo_sapiens_refseq/116_GRCh38 /output/homo_sapiens_refseq/116_GRCh38
vep -i input.vcf -o output.txt --offline --cache --refseq \
  --dir_cache $HOME/.vep --cache_version 116 --assembly GRCh38 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

`--refseq` annotates the RefSeq accessions alone, `--merged` every Ensembl model plus the RefSeq accessions among the RefSeq-sourced ones, and `--all_refseq` every model ([Filtering, picking and transcript sets](cli/filtering-picking.md#transcript-sets)). Both are BAM-edited, so a run needs `--fasta` without `--use_given_ref` ([Cache conversion](cache-conversion.md#refseq-and-merged-caches)). vep-rs builds no RefSeq model natively and produces no merged cache.

## Naming the cache

The cache flags select the directory as Ensembl VEP does ([Cache, species and runtime](cli/cache-species-runtime.md)):

1. The cache root is `--dir_cache`, else a directory given to `--cache`, else `--dir` (default `$HOME/.vep`).
2. Under it the species directory is `--species` (default `homo_sapiens`, the Latin name only), suffixed `_refseq` or `_merged` under those flags.
3. Under that sits `<cache_version>_<assembly>`, where `--cache_version` defaults to the VEP release.

When one directory of that version exists, it is the cache and its assembly is read from its name; when several exist, `--assembly` chooses. An `--assembly` that disagrees with the directory name, or either disagreeing with the cache's `info.json`, is an error.

Every failure to resolve a cache exits 1 before the input is opened with one `ERROR: [cache] ...` line, in Ensembl VEP's words where it has them. Those words are `Cache directory <dir> not found`, `No cache found for homo_sapiens, version 116` and `Multiple assemblies found for cache version 116 ...`, the assemblies in name order. The metadata mismatch names `info.json`, a JSON cache's counterpart of `info.txt`.

`--json_cache <dir>` (alias `--full_cache_dir`) names the cache directory itself instead, whatever `--dir_cache`, `--species` and `--cache_version` say. Only `--assembly` is checked against its `info.json`. Either way the directory must hold `transcripts/` and `info.json`, so a Storable cache or a mistyped path fails instead of annotating nothing. `transcripts/` must also hold a shard under some chromosome directory: a cache with none is refused, where Ensembl VEP writes every record as `intergenic_variant` ([A cache holding no transcript shard](divergences/damaged-cache-shard.md#a-cache-holding-no-transcript-shard)).

`--assembly` is compared exactly and case-sensitively, as Ensembl VEP compares it, with one addition: `hg19` and `hg38` are read as `GRCh37` and `GRCh38`. Without `--fasta`, the first indexed `.fa` or `.fa.gz` by name inside the cache directory is the reference, unless `--no_fasta` is given. An unindexed file there is not read (Ensembl VEP builds a missing index itself).

`--show_cache_info` resolves the cache by the same flags and refuses a named input file that does not exist, as Ensembl VEP does. It prints the cache's source versions as `<key>\t<value>` lines sorted by key, and exits without reading the input.

```bash
# Ensembl VEP's layout: $HOME/.vep/homo_sapiens/116_GRCh38
vep -i input.vcf -o output.txt --offline --cache \
  --dir_cache $HOME/.vep --cache_version 116 --assembly GRCh38

# The directory itself
vep -i input.vcf -o output.txt --offline --json_cache /path/to/116_GRCh38
```

## Reference FASTA

`--fasta` takes a plain FASTA with its `samtools faidx` `.fai`, or a `bgzip`-compressed FASTA with the `.fai` and `.gzi` `samtools faidx` writes for it, the shape of Ensembl's `dna_index/` files and of the file `vep-cache-builder --auto` places in the cache directory. `vep` builds no index (Ensembl VEP builds one); the builder writes both for a `dna/` file it installs. Three files are refused at open: a gzip file that is not BGZF, whose `.fai` offsets need the whole file inflated; a `.gzi` that does not match the file's blocks, recompressed and not re-indexed; an `.fai` naming more bases than the file holds. So is a BGZF file holding a malformed block.

| Condition | Line printed |
| --- | --- |
| A gzip file that is not BGZF | `ERROR: [fasta] <file> is gzip but not BGZF: recompress with bgzip and index with samtools faidx` |
| A `.gzi` that does not describe the file's blocks | `ERROR: [fasta] BGZF index <file>.gzi does not describe <file> (<the first disagreement>): index the file again with samtools faidx` |
| An `.fai` naming bases past the file's end | `ERROR: [fasta] FASTA index names <n> bases beyond the end of <file>` |
| A BGZF block declaring more than 65,536 uncompressed bytes, the most a block holds | `ERROR: [fasta] malformed BGZF block at offset <n> in <file>`, `<n>` the block's byte offset in the file |
| `--fasta 0` or `--fasta ''`, a value Ensembl VEP reads as the flag unset ([compatibility page](cli/compatibility.md#where-configpm-is-not-followed)) | `ERROR: [fasta] Specified FASTA file/directory <value> not found` |
| `--hgvs` or `--hgvsg` with no FASTA in reach | `ERROR: [fasta] Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file (see --fasta)` |
| `--check_ref` without a FASTA | `ERROR: [fasta] Cannot check reference sequences (--check_ref) without a FASTA file (see --fasta)` |
| `--lookup_ref` without a FASTA | `ERROR: [fasta] Cannot look up reference sequences (--lookup_ref) without a FASTA file (see --fasta)` |
| `--use_transcript_ref` without a FASTA, which it reads a non-coding transcript's cDNA from | `ERROR: [fasta] Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref` |
| `--shift_3prime 1` without a FASTA, since the shift reads the reference | `ERROR: [fasta] Cannot shift variants (--shift_3prime 1) without a FASTA file (see --fasta)` |
| `--everything` without a FASTA | `INFO: Disabling --hgvs; using --offline and no FASTA file found` |
| `--fasta_dir` holding no compatible file | `WARNING: No compatible FASTA files found in <dir>` |

The table's `ERROR:` lines exit 1. `--everything` without a FASTA prints its INFO line on standard error (withheld under `--quiet`) and drops the HGVS notations it enabled, as in Ensembl VEP. An explicit `--hgvsg` beside it is refused after that line.

### A `--fasta` directory

A `--fasta` directory is read as Bio::DB::Fasta reads one: every indexed `.fa` in it, plus the `.fa.gz` vep-rs adds, as one set of sequences. Two files naming one sequence are refused (Ensembl VEP with Bio::DB::HTS takes one file).

### A `Y` read inside a pseudoautosomal region

On GRCh37 and GRCh38 a read of `Y` inside a pseudoautosomal region is served from `X` at the region's offset, as in Ensembl VEP. Ensembl's primary-assembly FASTA carries `N` there. The redirect follows the FASTA's own name for the sequence, so a file naming it `chrY` is read as written in both engines.

The redirect is armed by the assembly the cache resolves to, its `info.json` first, then its directory name, then `--assembly`; Ensembl VEP arms it from its assembly parameter (`BaseVEP.pm` line 525, `Utils/FastaSequence.pm` line 254). That parameter is set by `--assembly`, `--grch37` (`Config.pm` lines 399-403) or the layout scan (`CacheDir.pm` lines 303-356), which runs for `--dir_cache` and `--dir`, not under `--full_cache_dir`. Under `--offline --full_cache_dir` (`--json_cache` here) without `--assembly` or `--grch37` it stays unset, so Ensembl VEP reads a `Y` span inside the region as written and vep-rs reads `X`.

### The file `--fasta_dir` chooses

`--fasta_dir <dir>` chooses, as Ensembl VEP does, the `.fa` or `.fa.gz` whose name carries the species and assembly (case-insensitively), `toplevel` before any other name and compressed before plain. That file replaces `--fasta` and the cache directory's FASTA. Ensembl VEP's test (`BaseVEP.pm` line 506, `/\.fa.gz$/`) also opens a `.fa<any character>gz` file by content; vep-rs passes it over. A directory without one leaves those in force and warns (the table's line); an unopenable one is passed over silently, as in Ensembl VEP.

Two files tied on both preferences are ordered by name, the first taken, where Ensembl VEP's sort on the two preferences alone keeps the directory listing's order (`BaseVEP.pm` lines 500-509), so on such a directory the two engines can read different files. Here the filter reads the same assembly as the redirect, where Ensembl VEP reads its assembly parameter alone (line 497). With that parameter unset, Ensembl VEP's filter passes every file naming the species (line 508), and its preferences and the listing's order choose among every assembly's files.

### A contig the FASTA lacks

Under `--lookup_ref` a contig the FASTA lacks skips the record at every position with the could-not-fetch warning, where Ensembl VEP's one-base placeholder slice (`BaseVEP.pm` lines 424-434) gives a one-base record at position 1 the reference `N`. Under `--check_ref` the same contig skips the record at every position with the sub-slice warning and then the mismatch warning; an insertion, naming no reference base, passes under either flag in both engines (`Parser.pm` lines 600-602 and 636). So Ensembl VEP annotates a one-base record there whose reference is written `N`, and refuses one written otherwise with the mismatch warning alone (`Parser.pm` lines 614-631).

Under `--hgvs` the same contig gives these notations; neither engine warns:

| Record on the contig under `--hgvs` | Ensembl VEP, whose transcript slice reads `N` per base there (`Utils/FastaSequence.pm` lines 597-601) | vep-rs, which reads no reference base there |
| --- | --- | --- |
| A substitution | `c.764N>T`, the `N` standing for the genome's `G` | no HGVSc |
| A deletion | `c.632del`, as over the genome | no HGVSc |
| An insertion | `c.307_308insC` | `c.307_308insC` |
| HGVSp | as over the genome | as over the genome |

### A `Y` read on a FASTA without `X`

With the pseudoautosomal redirect in force and no `X` in the FASTA, a `Y` read's bases inside the region are `N`, as Ensembl VEP reads them (`Utils/FastaSequence.pm` lines 597-601). The bases outside the region come from `Y`, so the read keeps its length. Ensembl VEP returns the region's `N`s alone for a read that straddles the region's edge, two bases for a four-base request. Neither engine warns on the read: Ensembl VEP's `FASTA does not contain sequence for X` (line 599) is raised at warning level 2000, above the verbosity of 1999 a run sets (`BaseRunner.pm` line 491).

## Chromosome names

An input may name a chromosome as the cache does (`21`), with a `chr` prefix (`chr21`), or by an accession (`NC_000021.9`, `CM000683.2`). Each name is resolved to the cache's own before any lookup, as Ensembl VEP's `get_source_chr_name` resolves it, in this order:

1. the name itself, when the cache holds it;
2. else a synonym of it the cache holds;
3. else the name with a leading `chr` removed or added, when the cache holds that;
4. else a bare `M` as the cache's `MT`, under whatever name the cache holds it.

A name that resolves to nothing skips the record with Ensembl VEP's warning (`Chromosome NC_000022.11 not found in annotation sources or synonyms`). Not so when a `--custom` or `--regulatory_gff` file [lists it](custom-annotation.md#chromosomes).

The output keeps the input's spelling in `Location`, the JSON `seq_region_name` and a VCF input's `CHROM` column. A bare `M` the sources know as `MT` but not as `M` or `chrM` prints `MT` in `Location` and `seq_region_name`, as Ensembl VEP renames it. [A bare `M` under a cache whose synonyms name `chrM`](divergences/bare-m-chromosome.md) describes where the engines then differ.

### The synonyms file

The synonyms come from the cache's `chr_synonyms.txt`, copied from the Ensembl VEP cache by the converter. The native builder writes it from the core `seq_region_synonym` table with `--core-dumps`. `--synonyms <file>` replaces the cache's file rather than extending it. Each line holds one chromosome's whitespace-separated names, the reference name first, then its synonyms, in both directions.

An unreadable file is `ERROR: [config] Could not read synonyms file <file>: <reason>` for a `--synonyms` file and `ERROR: [cache] ...` for the cache's own (a directory included, which Ensembl VEP opens as an empty table, `BaseVEP.pm` line 644). A file holding a byte that is not UTF-8 is refused the same way, the reason `stream did not contain valid UTF-8`, where Ensembl VEP reads the file as bytes.

The `--synonyms` value is a file name as given (`-`, `0`, an empty value, surrounding whitespace, a leading `>` or `|`), where Ensembl VEP's two-argument `open` reads `-` as standard input, trims the whitespace and takes a false value as unset ([compatibility page](cli/compatibility.md#where-configpm-is-not-followed)).

A cache without the file resolves no synonym: a native build without `--core-dumps`, or the output of a converter older than 0.5.0. So `chrM` and an accession are skipped there, as in Ensembl VEP. Convert the cache again, or pass `--synonyms`.

### Three cases that differ from Ensembl VEP

Three cases differ from Ensembl VEP, which checks a record against the cache's names widened by their synonyms, a set it widens again in place at every miss (`Parser.pm` line 782), and then looks the unwidened name up one hop. Each passes its check there and annotates as intergenic; in vep-rs a bare `M` resolves to `MT` and every other such name is skipped as not found, whatever record came before it.

The first is a bare `M` on a cache whose file links `chrM` but not `M`. The second is a name that reaches the cache only through a synonym of its `chr`-stripped or `chr`-prefixed spelling (`chrNC_000021.9`). The third is a name two hops from a cache chromosome under a `--synonyms` table whose lines chain or open with whitespace (`21 ALPHA` then `ALPHA BETA`: `BETA`), admitted once an earlier record missed.

### A FASTA whose names differ from the cache's

A reference FASTA named by accession is read under the cache's names through the same table. So `--fasta` may name a RefSeq or GenBank assembly beside an Ensembl cache.

Beyond the table, a record is reached by its name with the prefix `chr` added or removed and, for a record named `MT`, `M`, `chrM` or `chrMT`, by any of the four names, matched in any letter case: aliases vep-rs builds in. Ensembl VEP tries the exact name, a synonym and the `chr` variants alone (`BaseVEP.pm` lines 413 and 689-738, `Utils/FastaSequence.pm` lines 380-389), so a cache's `MT` over a FASTA record named `chrM`, with no synonym linking the two, reaches no sequence there. That chromosome gets the one-base placeholder slice described above (`BaseVEP.pm` lines 425-426), so under `--check_ref` a record such as `MT 3308 T G` is skipped with the sub-slice and mismatch warnings, where vep-rs annotates it.

## What `info.json` carries

`info.json` carries `species` and `assembly` in every cache the converter and the builder write; the resolver compares `--assembly` with `assembly`. Its other keys:

- `sift` and `polyphen`, markers written only for a tool whose matrices the cache holds. A converted cache keeps a marker only beside the tool's `source_versions` entry ([SIFT and PolyPhen predictions](cache-predictions.md)).
- `source_versions`: on a native cache the core `meta` table's `assembly`, `gencode` and `genebuild`, or the GTF header's. The dbSNP, ClinVar, COSMIC and gnomAD versions an Ensembl cache's `info.txt` lists come from the variation database. The builder takes no version from it, so a native cache's headers print none.
- `variation_cols` and `var_type tabix`, only when the cache holds variation files; `variation_cols` names their columns.
- `contig_order`: the chromosome directories: autosomes by number, X, Y, MT, then every other contig by name; a `--core-dumps` build lists the dumps' karyotype first, by rank.
- `regulatory`, `regulatory_schema` and `cell_types`, only with a regulatory tree: shards under `regulatory/<chr>/` beside the table `regulatory/binding_matrices.json` ([Regulatory and motif features](regulatory.md)).
