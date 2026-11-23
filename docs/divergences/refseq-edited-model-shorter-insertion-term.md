# The term of an insertion three bases short of the bases an edited RefSeq model holds and the genome lacks

**The defect in one sentence.** Where a RefSeq model edited to the genome holds bases the genome lacks, inside its coding sequence, an insertion of those bases with three left out stands on all of them and takes three bases out of the model, in frame. Ensembl VEP reads the allele's length class from the genomic span and asks `inframe_deletion` only of a deletion. On a record whose reference allele holds bases, no predicate then accepts the row, and Ensembl VEP writes its default `intergenic_variant` inside a coding exon. vep-rs writes `inframe_deletion` and `IMPACT=MODERATE` on that row, beside the positions Ensembl VEP prints.

**Why this is a defect.** The row names a transcript, a gene and a strand and prints a cDNA, a CDS and a protein position, and `intergenic_variant` is defined as a variant between genes. Both engines map the insertion onto the model's added bases and build the same alternate coding sequence, with the allele in their place. With three bases fewer in than out, the reading frame is kept and one residue goes, which is what the description of `inframe_deletion` in Ensembl VEP's own table of terms states. Written as a plain insertion, the same bases between the same two genomic positions get another term from Ensembl VEP, `protein_altering_variant`, and neither term names the deletion.

**What vep-rs does instead.** vep-rs maps the insertion onto the same span and builds the same alternate coding sequence. It asks the in-frame deletion test as well of an allele that is longer than its genomic span and shorter than the CDS span it stands on, and it asks on the model's own codons. Where the test holds, the row carries `inframe_deletion` and `IMPACT=MODERATE` beside the positions Ensembl VEP prints, whether the record is a plain insertion or keeps both bases beside the point. Every other term is decided as on any row.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116 and `ensembl-vep` at release/116.3, the modules the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, line 446: the reference length is the genomic span, `($vf_end - $vf_start) + 1`, which is 0 for a plain insertion and 2 for a record that keeps both bases beside the point. Lines 369-387 set the allele's length against it: a longer allele is an `insertion` (line 384), whatever its span on the model.
- `Bio/EnsEMBL/Variation/Utils/Constants.pm`, lines 866-869: the `include` rule of `inframe_deletion` asks for `deletion => 1`. `_skip_oc` (`BaseVariationFeatureOverlapAllele.pm`, lines 317-333) skips a predicate whose rule the allele does not meet, so the predicate is never run on such an allele. Line 861 holds the term's description: "An inframe non synonymous variant that deletes bases from the coding sequence".
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1161-1185, `inframe_deletion`: the test holds when the alternate codons are shorter than the reference codons and are their first or last bases (line 1178), or what is left of them once whole codons are cut from inside (lines 1181-1185).
- The same file, the predicates that are asked: `inframe_insertion` (line 1111) needs alternate codons longer than the reference codons. `frameshift` (lines 1556-1564) takes the allele's length less the CDS span's, modulo three. `protein_altering_variant` (lines 375-393) returns false where the in-frame deletion test holds (line 388) and true at line 392 where no line before it refused. `coding_unknown` (lines 1626-1636) needs an undefined peptide or an `X`.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, line 285: with no predicate true, `$DEFAULT_OVERLAP_CONSEQUENCE` is assigned. `Bio/EnsEMBL/Variation/Utils/Constants.pm`, lines 636-650: the default is `intergenic_variant`, described as "A sequence variant located in the intergenic region, between genes".
- `Bio/EnsEMBL/Variation/TranscriptVariation.pm`, lines 441-468, `_get_ref_allele`: under `--use_transcript_ref`, which a cache of corrected models switches on (`Bio/EnsEMBL/VEP/CacheDir.pm` of `ensembl-vep`, lines 410-415), the reference allele of a record that holds bases is cut from the model's sequence over the row's cDNA span (lines 456-463). An insertion's is `-` (line 451).
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 872-874 and 877: on a model carrying an `_rna_edit` attribute, the reference allele is written over the CDS span before the reference codons are read. For the empty allele of a plain insertion this is the write of [`refseq_edited_model_inserted_span_codons`](refseq-edited-model-inserted-span-codons.md), which cuts the span out.
- `Bio/EnsEMBL/VEP/Parser/VCF.pm` of `ensembl-vep`, lines 301-308: on a record with more than one alternate allele, the first base is dropped only when the reference allele and every alternate allele start with the same base.

**What correcting it would change.** An `include` rule could let `inframe_deletion` be asked of an allele shorter than its span on the model, whatever its genomic span. The row of a record that keeps both bases beside the point would then carry `inframe_deletion` and `IMPACT=MODERATE`: the test holds on Ensembl VEP's own two codon strings there. The plain insertion's row needs the correction of [`refseq_edited_model_inserted_span_codons`](refseq-edited-model-inserted-span-codons.md) as well. Its reference codons are read after the span is cut out, and the test fails on them; on the model's own codons it holds. No position column changes: both engines print the same three.

**The class it produces.** `refseq_edited_model_shorter_insertion_term`: a `Transcript` row of such an allele on a record whose reference allele holds bases. Its `Consequence` is `intergenic_variant` in Ensembl VEP and `inframe_deletion` in vep-rs, with `IMPACT` `MODIFIER` against `MODERATE` and the same positions, amino acids and codons in both. The plain insertion's row, `protein_altering_variant` against `inframe_deletion`, is a row of [`refseq_edited_model_inserted_span_codons`](refseq-edited-model-inserted-span-codons.md): the write that entry describes decides its reference codons and, through them, its term. No comparison tool sets a row of the shape aside: a disagreement of this shape is reported as a difference and charged to vep-rs. The concordance figures are measured on the Ensembl cache, whose models hold no such bases, so they hold no row of the shape.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an insertion written where a RefSeq model holds four or more bases the genome lacks, inside its coding sequence, whose bases are the model's own there with three left out. The Ensembl caches of release 116 hold no such point, 0 in the 238,823 transcripts of GRCh37 and 0 in the 669,092 of GRCh38, so the divergence needs a cache that holds RefSeq models. Of the committed corpora, 7 hold models with added bases at some point. None of the 18,476 records in the `variants.vcf` files of all the corpora lies at a point of added bases: the nearest is 14 bases from one.

A set of 348 records written for the purpose on the RefSeq GRCh38 cache of release 116 gives 10,346 rows in Ensembl VEP under `--refseq`. Six of them have the shape, at two points: four between 19:14090084 and 19:14090085 on `NM_138352.3`, which holds 318 bases there, and two between 4:3074935 and 4:3074936 on `NM_002111.8`, which holds 6. Three are rows of records that keep both bases beside the point, with `intergenic_variant`: the rows of this class. Three are rows of plain insertions, with `protein_altering_variant`.

vep-rs writes `inframe_deletion` on all six, with the positions Ensembl VEP prints, and a row for each of the set's other 10,340. At a third point of the set, 1:16607258 in a 5' UTR, where two models hold 4 bases, both engines write `5_prime_UTR_variant` on all 12 rows of such records.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Three of the six bases `NM_002111.8` holds and the genome lacks, written two ways

**In short.** Two records written for this example insert the same three bases between 4:3074935 and 4:3074936 on GRCh38, one as a plain insertion and one with both neighbouring bases in its reference allele. The forward-strand RefSeq model `NM_002111.8` (HTT) holds six bases at that point that the genome lacks, its coding positions 111 to 116, and the three inserted bases are the first three of them. Under `--refseq` Ensembl VEP 116.3 writes `protein_altering_variant` on the model's row of the first record and `intergenic_variant` on that of the second. vep-rs writes `inframe_deletion` on both, with the same positions: the model's sequence loses three bases and one residue.

### 1.1 What a reader needs to know first

A RefSeq cache carries each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). A model the alignment corrected prints `BAM_EDIT=OK` on its rows, and an `_rna_edit` attribute names each difference as `<start> <end> <bases>` in cDNA coordinates. An insertion edit has a start one past its end: `256 255` with six bases says the model holds six bases after cDNA position 255 that the genome lacks. A run on such a cache switches on `--use_transcript_ref`, which reads a row's reference allele from the model in place of the input.

An insertion is written between two neighbouring genomic bases. Its span on a transcript runs from one past the base before it to one before the base after it, which on an ordinary model is empty. Where the model holds added bases between the two, the span is those bases. The row prints them as its `cDNA_position` and `CDS_position`, and the alternate coding sequence both engines build carries the allele in their place. An allele shorter than the span leaves the model shorter by the difference.

A VCF record can write one insertion in more than one way. With the base before the point in both alleles, both engines drop the shared base and read a plain insertion, whose reference allele is empty. A record with two alternate alleles is read whole unless every allele starts with the same base. The second record of this example holds both bases beside the point in its reference allele, and a second alternate allele that changes the first of them. Nothing is dropped, and the span on the model is the added bases and the two bases beside them.

### 1.2 The records

```
4	3074935	three_of_six	A	AGCA	.	.	.
4	3074935	three_of_six_neighbours_kept	AA	AGCAA,CA	.	.	.
```

The first record inserts `GCA` between 4:3074935 and 4:3074936. The second writes the two genomic bases there, `AA`, as its reference allele. Its first alternate allele is those two bases with the same three between them, and its second replaces the first base alone.

The RefSeq GRCh38 cache of release 116 holds, for `NM_002111.8`, the edit `256 255` with six bases: cDNA 256-261, coding positions 111-116, between 4:3074935 at cDNA 255 and 4:3074936 at cDNA 262. Section 1.2 of [`refseq_edited_model_inserted_span_codons`](refseq-edited-model-inserted-span-codons.md) lists the model's mapper pairs. The six bases are `GCAGCA`, the same three twice, so the inserted three are the model's six with three left out. Both engines run `--refseq` with `--fasta`, which a run under `--use_transcript_ref` requires.

### 1.3 What each engine writes

The model's rows for the inserted bases, in the default format. Ensembl VEP 116.3:

```
three_of_six	4:3074935-3074936	GCA	3064	NM_002111.8	Transcript	protein_altering_variant	256-261	111-116	37-39	QQP/QQ	caacagccg/caGCAa	-	IMPACT=MODERATE;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK
three_of_six_neighbours_kept	4:3074935-3074936	AGCAA	3064	NM_002111.8	Transcript	intergenic_variant	255-262	110-117	37-39	QQQ/QQ	cAGCAGCAA/cAGCAA	-	IMPACT=MODIFIER;STRAND=1;GIVEN_REF=AA;USED_REF=AGCAGCAA;BAM_EDIT=OK
```

vep-rs:

```
three_of_six	4:3074935-3074936	GCA	3064	NM_002111.8	Transcript	inframe_deletion	256-261	111-116	37-39	QQQ/QQ	cagcagcaa/caGCAa	-	IMPACT=MODERATE;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK
three_of_six_neighbours_kept	4:3074935-3074936	AGCAA	3064	NM_002111.8	Transcript	inframe_deletion	255-262	110-117	37-39	QQQ/QQ	cAGCAGCAA/cAGCAA	-	IMPACT=MODERATE;STRAND=1;GIVEN_REF=AA;USED_REF=AGCAGCAA;BAM_EDIT=OK
```

On the second record's row, the row of this class, the engines differ in `Consequence` and in the `IMPACT` that follows it. They agree on the three position columns, on `Amino_acids` and on `Codons`. On the first record's row they differ in `Consequence` and in the reference halves of `Amino_acids` and `Codons`, `QQP` against `QQQ`, the difference of [`refseq_edited_model_inserted_span_codons`](refseq-edited-model-inserted-span-codons.md).

The second record's other allele gets `protein_altering_variant` on the model in both engines. Both records' rows on `NM_001388492.1`, a model of the same gene without the six bases, are equal as well. For the inserted bases they carry `inframe_insertion`, `CDS_position` 110-111, `Protein_position` 37 and `Q/QQ`.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP chooses the predicates to ask from two lengths: the reference length, which line 446 takes from the genomic span, and the allele's length. The first record's span holds no base and its allele three. The second's holds two bases and its allele five. Both alleles are the longer side, so both are `insertion` alleles (line 384). The `include` rule of `inframe_deletion` asks for a `deletion`, and the predicate is skipped on both rows.

On the second record's row the reference allele is the model's own eight bases over cDNA 255-262, which the row prints as `USED_REF`. The reference codons are the model's nine bases over residues 37 to 39. The alternate codons hold six: coding position 109, then the five-base allele in place of positions 110 to 117. `inframe_insertion` needs alternate codons longer than the reference codons. `frameshift` subtracts the span from the allele's length: 5 - 8 = -3, a multiple of three.

`protein_altering_variant` asks the in-frame deletion test and returns false where it holds. It holds on this row: the alternate codons are the last six bases of the reference codons (line 1178). `coding_unknown` needs an undefined peptide or an `X`, and the row has three residues against two. With no term left, line 285 assigns the default, `intergenic_variant`.

The first record's reference allele is empty (line 451). Line 873 writes it over coding positions 111-116, which cuts the six bases out, and line 877 reads nine bases from position 109 of the shortened sequence: positions 109 and 110, then 117 to 123. The reference half of `Codons` is that string, and the reference half of `Amino_acids` its translation, `QQP`.

The alternate codons are positions 109 and 110, the three inserted bases and position 117, the same six as on the second record's row. They are neither the first nor the last six bases of that reference string, nor the string with three bases cut from inside. The in-frame deletion test fails, no line before it in `protein_altering_variant` refuses, and line 392 returns true.

### 1.5 Why vep-rs's row is the correct one

The records' flanking bases are 4:3074935 and 4:3074936, cDNA 255 and 262 on the model. Between them lie cDNA 256 to 261, coding positions 111 to 116: the six bases the genome has no position for. Coding positions 109 to 111 are the codon of residue 37, since 3 × 37 = 111, and 115 to 117 the codon of residue 39: the rows' `Protein_position` 37-39 in both engines.

The model's nine bases there are three glutamine codons, the `QQQ` Ensembl VEP prints as the reference on the second record's row. The alternate window keeps positions 109 and 110, holds the three inserted bases in place of 111 to 116, and keeps 117: 9 - 6 + 3 = 6 bases. They are the model's own codons of residues 38 and 39, two glutamines where the model has three.

Three bases stand on six: the same frame, three bases and one residue fewer, and no base the model does not hold. That is an `inframe_deletion` by the description Ensembl VEP's table gives the term. `intergenic_variant` says the variant lies between genes, and the same row places it on coding positions 110 to 117 of the model. `protein_altering_variant` follows from a reference string that leaves out the six bases the row's `CDS_position` names. On the model's own codons Ensembl VEP's test holds, and its predicate declines the term.

The model without the six bases shows the same reading from the other side. On `NM_001388492.1`, which holds no added bases there, the three bases are an addition, and both engines write `inframe_insertion`. The term follows the sequence of the model the row is about: three bases more on one model, three fewer on the other.

### 1.6 Try it yourself

With the records of 1.2 saved as `records.vcf` under the two VCF header lines, the command as vep-rs takes it:

```
vep -i records.vcf -o out.txt --offline --refseq \
  --json_cache /path/to/json_cache/homo_sapiens_refseq/116_GRCh38 \
  --fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --species homo_sapiens --assembly GRCh38
grep NM_002111.8 out.txt
```

The cache is the conversion of Ensembl's RefSeq GRCh38 cache of release 116 ([RefSeq and merged caches](../cache-setup.md#refseq-and-merged-caches)). vep-rs writes `inframe_deletion` on the model's two rows of the inserted bases. Ensembl VEP 116.3 takes `--refseq` and `--fasta` with the cache options of [Reproducing these examples](reproducing.md), its cache root holding `homo_sapiens_refseq/116_GRCh38`, and writes `protein_altering_variant` and `intergenic_variant`.

No committed GRCh38 corpus holds a model with such bases inside its coding sequence. The reproduction `scripts/divergences/reproductions/refseq_edited_model_shorter_insertion_term.toml` therefore copies the one region file that holds the model, `4/3000001-4000000`, from the RefSeq cache it is given and converts it with `tools/vep-storable-to-json` before it runs the two engines.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`tools/vep-formats-compare`) has no class of this shape, and no committed corpus holds a row of it. A disagreement of the shape is reported as a difference in `Consequence` and `IMPACT` and charged to vep-rs.

The plain insertion's row is reported as well, in `Consequence`. Its `Amino_acids` and `Codons` are the difference of [`refseq_edited_model_inserted_span_codons`](refseq-edited-model-inserted-span-codons.md), whose class sets both aside where a comparison names it and gives the cache vep-rs read.
