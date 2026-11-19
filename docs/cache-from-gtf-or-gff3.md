# A GTF or GFF3 as the transcript model

A GTF or GFF3 serves as the transcript model in two steps: `vep-cache-builder` builds a cache from the file and the genome, then `vep` annotates against it.

```bash
vep-cache-builder --species homo_sapiens --assembly GRCh38 --release 116 \
  --gtf annotation.gtf.gz --genome-fasta genome.fa --codon-table MT=2 \
  --output-dir gtf_cache
vep -i input.vcf -o output.txt --offline \
  --json_cache gtf_cache --fasta genome.fa
```

In the commands, `genome.fa` is the assembly's genome FASTA with its `samtools faidx` index. `annotation.gtf.gz` and `input.vcf` are the user's own files.

`--gff3 <file>` takes a GFF3 the same way, plain, gzipped or bgzipped, so the indexed file Ensembl VEP needs serves both engines. Nothing is downloaded: `--release` labels the cache and `--assembly` names the assembly a run is checked against. `--codon-table MT=2` names the mitochondrion's genetic code, which no GTF or GFF3 carries; a file without an organelle region needs no such flag.

`vep --gff <file>` and `vep --gtf <file>`, which Ensembl VEP reads as an annotation source in place of its cache, have no in-process counterpart. vep-rs reads transcripts from its cache alone. Ensembl VEP's file source prints a `SOURCE` field naming the file on every transcript row, which a run on the built cache does not, a documented intended divergence, [`file_cache_writes_no_source_field`](divergences/file-cache-writes-no-source-field.md).

An Ensembl GTF built this way carries what Ensembl VEP's own `--gtf` reading carries, not what its cache carries. For the cache's content, build from the Ensembl GFF3 with the GTF beside it, as the [Quick start](getting-started/quick-start.md#2-get-a-transcript-cache) does, and `--core-dumps` for the rest.

## How the file is read

The file is read as Ensembl VEP's `BaseGXF`, `GFF.pm` and `GTF.pm` read it:

### Lines read and lines skipped

- Lines. The lines read, by file:

  | File | Lines read |
  | --- | --- |
  | A GTF | `gene`, `transcript`, `exon`, `CDS` and `stop_codon` |
  | A GFF3 | `gene`, `pseudogene`, `ncRNA_gene`, `mRNA`, `transcript`, `ncRNA`, `lnc_RNA`, `miRNA`, `primary_transcript`, `tRNA`, `rRNA`, `snRNA`, `snoRNA`, the gene segments, `exon`, `CDS` and the rest of Ensembl VEP's list |

  A GTF's stop codon is read as one more CDS segment, so the coding sequence runs through it as an Ensembl GFF3's does. Every other line is skipped and the skipped types are logged once with their counts. Ensembl VEP 116.2's `--gff` source skips an `ncRNA_gene` line and drops every transcript whose `Parent` names one, a documented intended divergence, [`file_cache_reads_ncrna_gene_lines`](divergences/file-cache-reads-ncrna-gene-lines.md).
- Ensembl lines. An Ensembl GFF3 types each line by its biotype. So a TEC transcript is an `unconfirmed_transcript` line and a scaRNA one a `scRNA` line, neither in the list. From the genome alone those lines are skipped as Ensembl VEP's `--gff` source skips them, and their transcripts are absent from the cache.

  A build of the cache's content is one reading the release's protein FASTA or core dumps (the Transcripts bullet below). It reads every `gene:` and `transcript:` line, the prefix in any letter case, whatever its type, since the Ensembl cache carries those transcripts.

### How records are linked, dropped and merged

- The tree. A GFF3 record is named by its `ID` (else `Name`, `id` or `name`) and belongs to each record of its `Parent` list. A GTF's exons belong to their `transcript_id` and a transcript to its `gene_id`.

  A record with no parent is a gene when its type names one, else a transcript with no gene line. A build of the cache's content reads a `gene:` line as a gene whatever its type. The types that name a gene are:
  - `gene`
  - `pseudogene`
  - `ncRNA_gene`
  - a gene segment
  - `processed_transcript`
  - `RNA`
- Absent values. An attribute whose value is empty or `0` is absent, as it is to Ensembl VEP. So each rule of these lists takes its next alternative. `ID=0;Name=g1` names the record `g1`, and a GFF3 `Parent=0` alone is no parent. A GTF transcript whose source column is `0` and whose attributes name no biotype has none.
- Dropped records. A record whose parent list opens with a parent in no line, or in a skipped line, is placed as its other parents place it, else dropped with its children. The parents are named in a warning, each case apart. A record under a transcript (a `miRNA` under a `primary_transcript`) is dropped without one.

  A record whose parent list opens with `0` and names nothing set (`Parent=0,0`; a GTF transcript's `gene_id "0"`) is dropped without a warning. Ensembl VEP drops it too (`BaseGXF.pm` lines 389 and 430-435). A `protein_coding` transcript without an exon line, and a transcript whose biotype no rule names, are dropped and named.
- An empty first parent field. A record whose list opens with an empty field is a record of its own, or the transcript of a gene a later field names. Ensembl VEP's `--gff` source gives such a transcript the gene's identifiers or none by hash order, a documented intended divergence, [`gff_parent_list_with_an_empty_field`](divergences/gff-parent-list-with-an-empty-field.md).
- CDS lines sharing a start. Two CDS lines of one transcript at one start are both read: every base either covers is coding, and the translation is named by the earlier line of the file on the forward strand and by the later on the reverse. Ensembl VEP keeps one of the two lines alone and loses what the other stated, a documented intended divergence, [`gff_cds_lines_sharing_a_start`](divergences/gff-cds-lines-sharing-a-start.md).

### Biotype and identifiers

- Biotype. A GFF3 transcript's is its `biotype`, `transcript_type` or `transcript_biotype`. Without one, an `mRNA` is `protein_coding`, an `ncRNA` takes its `ncrna_class` and a `V_gene_segment` is `IG_V_gene`. A transcript of a gene described `microRNA ...` is `miRNA`; otherwise `gbkey` is the biotype. A GTF transcript's is its `transcript_biotype`, `transcript_type` or `biotype`, else its source column.
- Identifiers. The transcript identifier is `transcript_id`, else the record's identifier, either without a `gene:` or `transcript:` prefix. The gene identifier is the transcript line's `gene_id`, else the `GeneID` of its `Dbxref` (`7157` for a RefSeq TP53 transcript), else the gene line's. The symbol is the transcript line's `gene_name` or `gene`, else the gene line's `name` or `Name`.

  The transcript version is `version`, else `transcript_version`, and the gene version is the gene line's `version`; a GTF's `gene_version` is not read. A `version` that is not a whole number, or is above 4294967295, the largest the cache holds, refuses the build, naming the file and line, as does a protein FASTA header's version above it; Ensembl VEP prints the text, [`gff_version_not_a_number`](divergences/gff-version-not-a-number.md).

### The translation and attribute values

- The translation. Its identifier and version are the `protein_id` and `version` of the CDS line that opens the coding sequence, the line Ensembl VEP reads them from. A line without `protein_id` gives the translation none, and its rows no HGVSp, unless a build of the cache's content reads the core dumps, which name it; Ensembl VEP writes an HGVSp on the empty identifier, a documented intended divergence, [`gff_cds_without_protein_id_hgvsp`](divergences/gff-cds-without-protein-id-hgvsp.md). The version is that line's `version`, else 1. Without one, a build of the cache's content takes a GTF line's `protein_version`, else the core dumps' where they hold it, else the header's, else none; a `version` that is not a whole number refuses the build, where Ensembl VEP prints its text, [`gff_version_not_a_number`](divergences/gff-version-not-a-number.md).
- Attribute values. A GTF value is read as `GTF.pm` reads it: the text from the key to the next `;`, every `"` removed wherever it stands, the `%XX` escapes decoded, so `"A%2CB"` is `A,B` and `"%30"` the absent `0`. A line repeating a key read as one value is refused, naming the key; Ensembl VEP reads the list: [`gtf_repeated_attribute_key`](divergences/gtf-repeated-attribute-key.md). A GFF3 value is the text between its pair's first and second `=`, as `GFF3.pm` cuts it, its `%XX` escapes decoded; Ensembl VEP's `--gff` source prints them as written, a documented intended divergence, [`gff_attribute_escapes_read_raw`](divergences/gff-attribute-escapes-read-raw.md). A line's end and padding blanks are not read; Ensembl VEP keeps the blanks of an unterminated last value ([section 2 there](divergences/gtf-repeated-attribute-key.md#2-a-last-value-padded-with-blanks-at-the-end-of-its-line)).

### Attributes, miRNA arms, kept transcripts and sequence

- Attributes. These are read:
  - `Ensembl_canonical`
  - `MANE_Select`
  - `MANE_Plus_Clinical`
  - `gencode_basic` (spelled `basic` in the GRCh37 release-87 files and the GENCODE GTF)
  - `gencode_primary`
  - `transcript_support_level`

  So is the CCDS identifier, an Ensembl GFF3's `ccdsid` or a GTF's `ccds_id` as the line carries it, which Ensembl VEP's `--gff` and `--gtf` sources never print or rank under `--pick`, a documented intended divergence, [`file_cache_reads_ccdsid`](divergences/file-cache-reads-ccdsid.md).

  These are not read, so their columns stay empty:
  - `cds_start_NF` and `cds_end_NF`
  - a RefSeq GFF3's CCDS
  - APPRIS
  - the symbol's source
  - the HGNC identifier
  - the xref lists
  - protein features
  - the gene phenotype count
  - the sequence edits
- Mature miRNA arms. No GTF or GFF3 carries them, so `mature_miRNA_variant` is never assigned on a cache built without `--core-dumps`. On such a cache a variant inside an arm reads `non_coding_transcript_exon_variant` or the intron term. Ensembl VEP reads the arms from the `miRNA` transcript attribute of its database (`VariationEffect.pm` line 540), which `--core-dumps` supplies ([Cache setup](cache-setup.md#build-a-cache)).
- Transcripts. Every transcript is kept, readthrough and `artifact` ones included, as Ensembl VEP's `--gff` and `--gtf` sources keep them (`BaseGXF.pm` `_create_transcript`). A build that reads the release's protein FASTA or core dumps is the cache's content and leaves both out, as Ensembl VEP's cache does. Such a build passes `--fasta`, `--core-dumps`, or `--include-predictions`, which reads the protein FASTA for its peptide keys ([Cache setup](cache-setup.md#build-a-cache)).
- Sequence. The coding sequence, UTRs and peptide come from the genome ([`--genome-fasta`](cache-setup.md#build-a-cache)), the peptide under the region's codon table. That table is 1 unless `--codon-table` names another. A region named like an organelle genome stops the build until the flag names its table ([The codon table of each sequence region](cache-setup.md#the-codon-table-of-each-sequence-region)). Such names are `MT`, `Mito`, `mitochondrion_genome` and the rest.

  Ensembl VEP's `--gff` and `--gtf` sources take table 2 for a region named `MT` and 1 for every other (`BaseGXF.pm` lines 765 to 767). So `--codon-table MT=2` gives a human cache the table Ensembl VEP annotates that source under. That table reaches every codon Ensembl VEP translates for a row; the stored reference peptide stays under the standard code ([below](#the-reference-peptide-of-a-mitochondrial-transcript-read-from-a-gtf-or-gff3), tier 2).

## Transcript filters on the built cache

`--transcript_filter` and `--exclude_predicted` run on the built cache as on any other. Ensembl VEP applies the same filter set to its `--gtf` and `--gff` sources. So a filter on a key both carry keeps the same transcripts on both engines.

`--gencode_basic` and `--gencode_primary` filter the built cache as they filter every cache, where Ensembl VEP's `--gtf` and `--gff` sources accept both flags and annotate every transcript, a documented intended divergence, [`file_cache_honours_gencode_basic`](divergences/file-cache-honours-gencode-basic.md). A transcript Ensembl VEP builds from the file holds fewer keys than its cache's, and the built cache holds the cache's keys. [The transcript filter page](divergences/transcript-filter-fields.md#a-gtf-or-gff3-as-the-transcript-source) lists the keys that read differently.


## The reference peptide of a mitochondrial transcript read from a GTF or GFF3

One behaviour of Ensembl VEP 116.2's `--gtf` and `--gff` sources is not reproduced by a cache built from the same file, a documented class. This section names the rows a comparison sets aside, and no comparator applies the class, since no committed run reaches its shape, as [Intended divergences from Ensembl VEP](intended-divergences.md) lists it; its [How to read the counts](intended-divergences.md#how-to-read-the-counts) defines the terms.

The class carries its tier under that page's rule: tier 1 is a shape the chromosome 21 population cuts carry with no flag, tier 2 one that needs a flag or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts are over the cuts' 15,646,481 records, the 111 golden corpora and their 392 flag sets.

**The defect in one sentence.** Ensembl VEP's `--gtf` and `--gff` sources translate every transcript's reference peptide under the standard genetic code, the `MT` transcripts included, while stamping those transcripts with the mitochondrial table it then uses for the altered codon.

**Why this is a defect.** The source builds each transcript's peptide once, and to do so it detaches the transcript from its sequence region, because the translation would otherwise look the region's codon table up in a database it does not have. Detached, the translation falls back to the standard code. Two lines later the same routine records the table as 2 for a region named `MT`, and every later step honours that stamp: the altered codon, the consequence terms, the residue a substitution's HGVSp names.

The stored peptide therefore disagrees with the transcript's own stamp, and a codon the two codes read differently gets a reference residue the mitochondrial ribosome does not make wherever the notation reads the stored peptide.

**What vep-rs does instead.** vep-rs translates the peptide of a cache built from a GTF or GFF3 under the region's table, the one `--codon-table` names, so the reference residue and the altered residue of an `MT` row come from the same code.

### The source lines, the class and how often it occurs

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `ensembl` core and `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/BaseGXF.pm`, lines 761-768, `_add_translation`: line 763 deletes the transcript's slice, line 764 stores `$translation->seq`, lines 765-768 stamp `codon_table` 1, or 2 for `MT`.
- `Bio/EnsEMBL/Transcript.pm`, lines 2103-2112, `translate`: the codon table is read from the slice's attributes when the slice is defined, else 1.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 692-736, `peptide`: a row's reference residue is translated from its codon under the stamped table (line 712 the table, line 736 the translation), so a substitution reads the same on both engines. The stored peptide, read through `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm` line 1117, reaches the HGVSp paths of the same `TranscriptVariationAllele.pm`: an insertion's flanking residues (lines 2368-2374, `_get_surrounding_peptides`), a frameshift's and a stop loss's reference residues (lines 2320-2336 and 2477-2494).

**What correcting it would change.** The `MT` transcripts of such a source would carry the table-2 peptide, as the same transcripts carry it in Ensembl VEP's own cache, and their HGVSp would agree with vep-rs's.

**The class it produces.** `gxf_organelle_peptide_standard_table`, documented only: the `HGVSp` of a row on an `MT` transcript of a cache built from a GTF or GFF3 whose notation reads a reference residue from the stored peptide, where Ensembl VEP's residue is the standard-table reading of the codon. Those notations are an insertion's flanking residues and a frameshift's or a stop loss's reference residue. `Consequence`, `Amino_acids` and a substitution's HGVSp agree, since they come from the codon under table 2. The `GRCh38-gtf-builder` corpus's GTF holds chromosome 21 alone and no comparison case annotates an `MT` record against a cache built from a GTF or GFF3, so no manifest names the class, no comparator applies it and no row is set aside.

**How often the shape occurs.** Tier 2: the shape is an organelle indel on a GTF- or GFF-built source under `--hgvs`; the flanking residue's codon decides whether the stored peptide reads differently. No chromosome 21 cut carries a mitochondrial record; 247 of the 487 mitochondrial records of the committed corpora (507,187 per million) are indels, in 7 corpora, and no committed flag set runs `--hgvs` on a cache built from a GTF or GFF3.

### An in-frame insertion after codon 31 of MT-ND1

**In short.** Codon 31 of MT-ND1 (`ENST00000361390`, `MT:3307-4262`, forward strand) is `ATA`, methionine under the vertebrate mitochondrial code and isoleucine under the standard one. The insertion `MT 3399 . A AGCT`, c.93_94insGCT, adds an alanine after it. On a cache built from the Ensembl GTF with `--codon-table MT=2`, vep-rs writes `inframe_insertion`, `Amino_acids` `-/A` and `HGVSp=ENSP00000354687.1:p.Met31_Gln32insAla`. Ensembl VEP 116.2 reading the same GTF through `--gtf` holds the standard-table peptide, whose residue 31 is `I`, and writes `p.Ile31_Gln32insAla` beside the same terms and `Amino_acids`; from its own cache, whose peptide was translated under table 2, it writes `p.Met31_Gln32insAla`.

The substitution `MT 3399 . A G`, which makes the codon `ATG`, is `synonymous_variant`, `M`, `p.Met31%3D` on both engines and both routes, since a substitution's residue is translated from the codon under table 2.

## Where next

- [Cache setup](cache-setup.md#build-a-cache): the builder's other inputs and the GFF3 or GTF lines it refuses.
- [Quick start](getting-started/quick-start.md#2-get-a-transcript-cache): a build of the cache's content from the Ensembl GFF3 with the GTF beside it.
- [The transcript filter page](divergences/transcript-filter-fields.md#a-gtf-or-gff3-as-the-transcript-source): the filter keys that read differently on a GTF or GFF3 source.
- [Intended divergences](intended-divergences.md): where the two engines differ, and why.
