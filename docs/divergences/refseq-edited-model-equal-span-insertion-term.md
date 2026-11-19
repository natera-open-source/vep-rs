# The term of an insertion as long as the bases an edited RefSeq model holds and the genome lacks

**The defect in one sentence.** Where a RefSeq model edited to the genome holds bases the genome lacks, an insertion of as many bases written at that point stands on them one for one, changing residues and keeping the length. Ensembl VEP reads the allele's length class from the genomic span, never asks `missense_variant` of an insertion, and writes its default `intergenic_variant` on a row inside a coding exon. vep-rs writes `missense_variant` with `IMPACT=MODERATE` beside the same positions.

**Why this is a defect.** The row names a transcript, a gene and a strand and prints a cDNA, a CDS and a protein position, and `intergenic_variant` is defined as a variant between genes. Both engines map the insertion onto the model's added bases and build the same alternate coding sequence, with the allele in their place. With as many bases in as out, the reading frame is kept and the residues the span touches change, which is what the description of `missense_variant` in Ensembl VEP's own table of terms states. One or two bases inserted at the same point are a `frameshift_variant` in both engines.

**What vep-rs does instead.** vep-rs maps the insertion onto the same span and builds the same alternate coding sequence. It asks `missense_variant` as well of an allele that is longer than its genomic span and as long as its CDS span, and writes the term with `IMPACT=MODERATE` beside the positions Ensembl VEP prints. Every other term is decided as on any row: an allele equal to the model's own bases is a `synonymous_variant`, one that brings a stop codon a `stop_gained`, and one whose length differs by other than a multiple of three a `frameshift_variant`.

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `ensembl` at release/116, the modules the Ensembl VEP 116.2 image holds):

- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, line 446: the reference length is the genomic span, `($vf_end - $vf_start) + 1`, which is 0 for an insertion. Lines 371-386 set the allele's length against it: a longer allele is an `insertion` and `increase_length` (lines 384-386), whatever its span on the model.
- `Bio/EnsEMBL/Variation/Utils/Constants.pm`, lines 886-890: the `include` rule of `missense_variant` asks for `increase_length => 0`, so the predicate is never run on such an allele. Line 881 holds the term's description: "changes one or more bases, resulting in a different amino acid sequence but where the length is preserved".
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1085-1098, `missense_variant`: line 1097 is true when the two peptides differ and have one length, so the predicate, once asked, accepts the row.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, the predicates that are asked: `frameshift` (lines 1556-1564) takes the allele's length less the CDS span's, modulo three. `inframe_insertion` (line 1111) needs alternate codons longer than the reference codons. `protein_altering_variant` (line 384) leaves peptides of equal length to `synonymous_variant` and `missense_variant`. `coding_unknown` (lines 1625-1638) needs a peptide that is undefined or holds an `X`.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, line 285: with no predicate true, `$DEFAULT_OVERLAP_CONSEQUENCE` is assigned. `Bio/EnsEMBL/Variation/Utils/Constants.pm`, lines 636-649: the default is `intergenic_variant`.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 263-264, and `Bio/EnsEMBL/Mapper.pm`, lines 485-558, `map_insert`: an insertion's CDS span runs from one past the base before the insertion point to one before the base after it.

**What correcting it would change.** The `include` rule would let `missense_variant` be asked of an allele as long as its span on the model, whatever its genomic span. Line 1097 then accepts a row whose residues change and whose length is kept, as it does for a substitution, and the row carries `missense_variant` and `IMPACT=MODERATE`. No position column changes: both engines print the same three.

**The class it produces.** `refseq_edited_model_equal_span_insertion_term`: a `Transcript` row of an insertion whose `Consequence` is `intergenic_variant` in Ensembl VEP and `missense_variant` in vep-rs, with `IMPACT` `MODIFIER` against `MODERATE` and the same cDNA, CDS and protein positions in both. No comparison tool sets a row of the shape aside: a disagreement of this shape is reported as a difference and charged to vep-rs. The concordance figures are measured on the Ensembl cache, whose models hold no such bases, so they hold no row of the shape.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an insertion written where a RefSeq model holds bases the genome lacks, inside its coding sequence, with as many inserted bases as the model holds there. A set of 51 hand-written insertions at and beside 8 such places of the RefSeq GRCh38 cache gives 884 rows in each engine under `--refseq`. Of its 36 rows of the shape, 31 differ this way, each a three-base insertion on a three-base span. The other 5 do not: 3 one-base rows are `synonymous_variant` in both engines, and 2 on mitochondrial stop codons are `coding_sequence_variant` in Ensembl VEP and `stop_lost` in vep-rs, another difference.

Among the ClinVar GRCh38 records, 32 insertions reach bases a model holds and the genome lacks, in any region of the model, 3 as written (2 of them inside a coding sequence) and 29 once shifted 3'. vep-rs writes 146 rows for them, and one row differs this way: record 551063 on `NM_015120.4`. No committed corpus holds an insertion on such bases.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Three bases inserted on the three bases `NM_015120.4` holds and the genome lacks

**In short.** ClinVar record 551063 inserts three bases between 2:73385942 and 2:73385943 on GRCh38. The forward-strand RefSeq model `NM_015120.4` holds three bases at that point that the genome lacks, its coding positions 75 to 77. Both engines print `CDS_position` 75-77 and `Protein_position` 25-26 on the model's row, with the same alternate codons. Ensembl VEP 116.2 writes `intergenic_variant` and `IMPACT=MODIFIER` there; vep-rs writes `missense_variant` and `IMPACT=MODERATE`.

### 1.1 What a reader needs to know first

A RefSeq cache carries each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). A model the alignment corrected prints `BAM_EDIT=OK` on its rows, and an `_rna_edit` attribute names each difference as `<start> <end> <bases>` in cDNA coordinates. An insertion edit has a start one past its end: `186 185` with three bases says the model holds three bases after cDNA position 185 that the genome lacks. The model's coordinate mapper numbers them, so the genomic bases on either side, neighbours on the genome, lie four positions apart on the model.

An insertion is written between two neighbouring genomic bases. Its span on a transcript runs from one past the base before it to one before the base after it, which on an ordinary model is empty. Where the model holds added bases between the two, the span is those bases. The row prints them as its `cDNA_position` and `CDS_position`, and the alternate coding sequence both engines build carries the allele in their place.

### 1.2 The record

```
2	73385942	551063	A	AGGC	.	.	.
```

The record is ClinVar variation 551063: the bases `GGC` inserted between 2:73385942 and 2:73385943. The RefSeq GRCh38 cache of release 116 holds, for `NM_015120.4`, the edit `186 185` with three bases, a coding sequence that starts at cDNA 112, and two mapper pairs that meet at the insertion point:

| cDNA | Genome |
| --- | --- |
| 1-185 | 2:73385758-73385942 |
| 189-438 | 2:73385943-73386192 |

Both engines run `--refseq --use_given_ref`. The second flag lets a run read an edited cache without a reference FASTA; with a FASTA in its place each engine writes the same term.

### 1.3 What each engine writes

The model's row in the default format. Ensembl VEP 116.2:

```
551063	2:73385942-73385943	GGC	7840	NM_015120.4	Transcript	intergenic_variant	186-188	75-77	25-26	EE/EA	gaagag/gaGGCa	-	IMPACT=MODIFIER;STRAND=1;BAM_EDIT=OK
```

vep-rs:

```
551063	2:73385942-73385943	GGC	7840	NM_015120.4	Transcript	missense_variant	186-188	75-77	25-26	EE/EA	gaggaa/gaGGCa	-	IMPACT=MODERATE;STRAND=1;BAM_EDIT=OK
```

The rows agree on the three position columns, on `Amino_acids` and on the alternate half of `Codons`. They differ in `Consequence` and in the `IMPACT` that follows it. The reference half of `Codons` differs too, `gaagag` against `gaggaa`: that is a difference of its own, in the bases Ensembl VEP reads as the reference codons, and no part of this entry.

ClinVar record 1167157 inserts two bases at the same point, `2 73385942 1167157 A AGG`. Its row on the model is `frameshift_variant` with the same three positions in both engines.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP chooses the predicates to ask from two lengths: the reference length, which line 446 takes from the genomic span, and the allele's length. An insertion's genomic span holds no base, so the reference length is 0 and the three-base allele is an `increase_length` allele (lines 384-386). The `include` rule of `missense_variant` asks for `increase_length => 0`, and the predicate is skipped.

Every predicate that is asked refuses the row. `frameshift` subtracts the CDS span from the allele's length: 3 - (77 - 75 + 1) = 0, a multiple of three. `inframe_insertion` needs alternate codons longer than the reference codons, and both strings hold six bases.

`protein_altering_variant` returns false for two peptides of one length, leaving them to `synonymous_variant` and `missense_variant`. `coding_unknown` needs an undefined peptide or an `X`, and the row has two residues on each side. With no term left, line 285 assigns the default, `intergenic_variant`.

### 1.5 Why vep-rs's row is the correct one

The record's flanking bases are 2:73385942 and 2:73385943. The first mapper pair ends at 73385942, cDNA 185, and the second starts at 73385943, cDNA 189. Counted from the coding start at cDNA 112, those are coding positions 185 - 111 = 74 and 189 - 111 = 78. Between them lie cDNA 186 to 188, coding positions 75 to 77: the three bases the genome has no position for, and the row's `CDS_position` in both engines.

Coding positions 73 to 75 are the codon of residue 25, since 3 × 25 = 75, and 76 to 78 the codon of residue 26: the row's `Protein_position` 25-26. The model's six bases there are `gaggaa`, the reference codons vep-rs prints, two glutamates. The alternate window keeps positions 73 and 74, holds the three inserted bases in place of 75 to 77, and keeps 78: `gaGGCa`, six bases, a glutamate and an alanine. The first two inserted bases equal the model's own 75 and 76 and the third differs from 77, so one base of codon 26 changes.

Three bases stand on three: 6 - 3 + 3 = 6 bases, the same length and the same frame, and residue 26 changed. That is a `missense_variant` by the description Ensembl VEP's table gives the term. `intergenic_variant` says the variant lies between genes, and the same row places it on coding positions 75 to 77 of the model.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache <cache>/homo_sapiens_refseq/116_GRCh38 \
  --species homo_sapiens --assembly GRCh38 \
  --refseq --use_given_ref
grep NM_015120.4 out.txt
```

`<cache>` is the conversion of Ensembl VEP's RefSeq GRCh38 cache of release 116 ([Cache conversion](../cache-conversion.md)). vep-rs writes `missense_variant` on the model's row. Ensembl VEP 116.2 takes the same flags with the cache options of [Reproducing these examples](reproducing.md) and its RefSeq cache, and writes `intergenic_variant`. With `AGG` as the alternate allele, record 1167157, both write `frameshift_variant`.

No committed corpus holds a model with such bases inside its coding sequence. The entry's reproduction therefore converts the one region file of Ensembl VEP's cache that holds the model, `2/73000001-74000000`, before it runs the two engines.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`scripts/concordance/compare_vep_formats.py`) has no class of this shape, and no committed corpus holds a row of it. A disagreement of the shape is reported as a difference in `Consequence` and `IMPACT` and charged to vep-rs.

Two more columns of such a row differ for causes of their own and stay differences as well: the reference half of `Codons`, as 1.3 shows, and, under `--hgvs`, the `HGVSc`.
