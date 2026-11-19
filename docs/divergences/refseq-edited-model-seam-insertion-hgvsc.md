# The HGVSc of an insertion standing on the bases an edited RefSeq model holds and the genome lacks

**The defect in one sentence.** Under `--hgvs`, for an insertion on bases that a RefSeq model holds and the genome lacks, Ensembl VEP builds the notation on the genome and writes an `ins` of the allele that keeps every one of those bases. On the model's own sequence the allele stands in their place, so the `ins` describes a sequence the record does not give. vep-rs sets the allele against those bases, clips the ends the two share and writes the change that is left.

**Why this is a defect.** An HGVSc is a statement on the sequence of the transcript it names. The two genomic bases around the insertion point both have a position on the model, and the model's added bases lie between them. Whatever sequence the record gives on the genome, on the model it is the model's own with those bases exchanged for the allele. Ensembl VEP's `ins` keeps the added bases and puts the whole allele beside them. That sequence is longer than the record's by as many bases as the model adds at the point: 3 in section 1, 11, 1 and 1 in section 2, and 9 in section 3.

The page shows the insertion in three forms:

| Form of the insertion | Section | What its row prints |
| --- | --- | --- |
| Written on such bases inside a coding sequence | 1 | the span of the bases as `cDNA_position` and `CDS_position`; both engines build its alternate codons with the allele in place of the bases |
| Written on such bases in an untranslated region | 2 | the span as `cDNA_position` alone; the region's term |
| Moved onto such bases by the 3' shift | 3 | the positions of the point the record writes, and the length of the move as `HGVS_OFFSET` |

**What vep-rs does instead.** vep-rs reads the model's bases at the point from the model's cached sequence and sets the allele against them. It clips the ends the two share, as for any allele that replaces reference bases, and names what is left: a substitution, a deletion, a duplication, an insertion between two neighbours, an inversion or a `delins`. While the 3' shift is on, the default, it places that change 3'-most on the model's sequence. Under `--shift_hgvs 0` the change stays where the clip leaves it. vep-rs writes no HGVSc where the allele equals the model's own bases.

On a row that the 3' shift moves onto such bases, vep-rs makes the notation at the shifted point: it reads the allele, rotated by the shift, in place of the model's bases there. The row's term, `Codons`, `Amino_acids` and position columns are read at the point the record writes, where the allele is added whole, and both engines print them alike.

Two neighbouring entries start from the same bases. The entry [`refseq_edited_model_seam_shift_and_duplication_hgvs`](refseq-edited-model-seam-shift-and-duplication-hgvs.md) covers an insertion that the 3' shift carries off them, and a duplication across them. The entry [`refseq_edited_model_unminimised_allele_hgvsc`](refseq-edited-model-unminimised-allele-hgvsc.md) covers an allele that keeps a base it shares with its reference, and says where vep-rs writes an insertion that the 3' placement carries to an exon's last base.

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `ensembl` at release/116, the modules the Ensembl VEP 116.2 image holds; `TranscriptVariationAllele.pm` stands for `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`):

- `TranscriptVariationAllele.pm`, lines 1474-1483, `hgvs_transcript`: the call of `hgvs_variant_notation` builds the notation on the transcript's genomic slice, from `_slice_start` and `_slice_end`. On that slice the insertion's two positions are neighbours.
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 591-594: an insertion that copies no bases beside it on the slice gets the type `ins`, between those two positions.
- `TranscriptVariationAllele.pm`, lines 1508 and 1517: for an `NM_` or `XM_` model with `_rna_edit` attributes, `get_misalignment_offset` (lines 1594-1624) adds up edit lengths. Line 1610 leaves out an edit that starts after the row's cDNA start, and line 1617 adds the length of an insertion edit.
- `TranscriptVariationAllele.pm`, lines 1524-1525: both slice positions are moved by that sum and numbered by `_get_cDNA_position` (line 2755). Lines 2766-2794 turn each moved position into a genomic base and count it from the start of its exon, as if the exon held no added base. Lines 2862-2867 count from the coding start, and lines 2848-2852 write a position past the stop codon as `*` and its distance.
- `TranscriptVariationAllele.pm`, line 1441: the insertion is first moved 3' along the genome, unless the shift is off. On a model with an `_rna_edit`, lines 1510-1513 then write the record's bases back as the record spells them.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 147-148 and 263-264, and `Bio/EnsEMBL/Mapper.pm`, `map_insert`, lines 504-558: the row's cDNA and CDS span of an insertion runs from one past the base before the insertion point to one before the base after it. Where the model numbers the two bases more than one apart, the span is the positions between them.
- The sum and the exon count are the routine that the entry [`refseq_edited_model_hgvsc_offset`](refseq-edited-model-hgvsc.md) describes, on rows where the two engines' notations differ in their coordinates alone. On this page the kind of change differs as well.

**What correcting it would change.** The notation would be built on the model's sequence: at the row's own span, or at the shifted point for an insertion that the 3' shift moves. The `ins` that keeps the added bases would give way to the change the allele makes in their place. The correction concerns the `HGVSc` alone.

**The class it produces.** `refseq_edited_model_seam_insertion_hgvsc`: the `HGVSc` of a `Transcript` row of an insertion whose `CDS_position` is a span of bases the model holds and the genome lacks, where Ensembl VEP's value is an `ins` of the allele and vep-rs's the change the allele makes to those bases. Of the three forms, the class holds the first alone, the rows of section 1.

The format comparator applies the class's test where a comparison names the class, and no comparison case or corpus manifest does. A row of the shape is therefore reported as a difference in `HGVSc` and charged to vep-rs. Section 1.7 states what the test sets aside and what stays a difference. Sections 2.7 and 3.7 show that it accepts no row of the other two forms. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an insertion on bases that a RefSeq model holds and the genome lacks, written there or moved there by the 3' shift, in a run under `--hgvs`. The tier counts models: the committed caches hold 33 models with such bases, in 3 corpora, none of them with the bases inside its coding sequence. One committed flag set runs a RefSeq or merged cache under `--hgvs`, and no record of a committed corpus stands on such bases.

A set of 51 hand-written insertions, at and beside 8 such places inside coding sequences of the RefSeq GRCh38 cache, gives 884 rows in each engine under `--refseq --hgvs`. vep-rs writes this entry's notation on 114 of them: a `delins` over the span on 95, an inversion on 15 and a one-base substitution on 4. Of the 114, 108 are on nuclear models, where Ensembl VEP writes an `ins` on every one. The other 6 are on mitochondrial models, whose rows print no `BAM_EDIT`.

Among the ClinVar GRCh38 records, 32 insertions reach such bases, 3 as written and 29 once shifted 3'. vep-rs writes this entry's notation for 13 of them, 2 written there and 11 shifted there, all on `NM_015120.4` at 2:73385942: 5 insertions between neighbours, 4 duplications, 2 substitutions, 1 `delins` and 1 deletion. The other 19 (8 on `NM_015120.4`, 10 on `NM_002111.8` and 1 on `NR_163972.1`) are insertions that the genome shows as a duplication of the bases before the place, and they are no part of this entry. The 32 are not split by region here.

No row of the 114 or of the 13 has an allele equal to the model's own bases or a model without a cached sequence.

In an untranslated region, three records of the ClinVar GRCh37 file reach such bases as written, under `--refseq --hgvs`: 1049039, 212026 and 335960, the records of section 2. How many insertions of that file the 3' shift carries onto such bases is not counted.

The RefSeq caches of release 116 hold 603 such spans on GRCh37 and 1,218 on GRCh38, and 303 and 694 of them lie in an untranslated region. The concordance datasets are annotated on the Ensembl caches, which hold no edited model, so no row of them has the shape in any of the three forms.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. One base inserted on the three bases `NM_015120.4` holds and the genome lacks

**In short.** `2 73385942 seam_one_base A AC` on GRCh38 inserts one base between 2:73385942 and 2:73385943. At that point the forward-strand RefSeq model `NM_015120.4` holds three bases that the genome lacks, its coding positions 75 to 77. Both engines print `CDS_position` 75-77 and `frameshift_variant` on the model's row. Under `--hgvs`, Ensembl VEP 116.2 writes `c.77_78insC`, an insertion after the third of those bases that keeps all three. vep-rs writes `c.75_77delinsC`, the three bases replaced by the one. The class's test accepts the row where a comparison names the class.

### 1.1 What a reader needs to know first

A RefSeq cache carries each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). Where the model holds bases that the genome lacks, the cache records them as an `_rna_edit` attribute, `<start> <end> <bases>` in cDNA coordinates. Such a model is the edited model of this page, and the rows of `NM_015120.4` print `BAM_EDIT=OK`. An insertion edit has a start one past its end: `186 185` with three bases says that the model holds three bases after cDNA position 185. The model's mapper numbers those bases. The genomic bases before and after them, neighbours on the genome, are therefore four positions apart on the model.

An insertion is written between two neighbouring genomic bases. Its span on a transcript runs from one past the base before the point to one before the base after it. Where the model holds added bases between the two, the span is those bases. The row prints it as its `cDNA_position` and, inside the coding sequence, as its `CDS_position`. Where the model holds no base between the two, the row prints the two neighbours' positions, as `125-126` does in section 3. `REFSEQ_OFFSET` prints the sum of edit lengths that Ensembl VEP adds to the row's notation, as 1.4 describes.

HGVS describes a change on the transcript's own sequence, at its 3'-most position. An `ins` names two neighbouring positions and the bases put between them, and it keeps every reference base. A `delins` names a span and the bases that replace it. Under `--hgvs` both engines move an insertion 3' along the genome before they write its notation, unless `--shift_hgvs 0` is given. A row whose insertion moved prints the distance as `HGVS_OFFSET`.

### 1.2 The record

```
2	73385942	seam_one_base	A	AC	.	.	.
2	73385942	551063	A	AGGC	.	.	.
2	73385942	1167157	A	AGG	.	.	.
```

The first record is written for this example: one base inserted between 2:73385942 and 2:73385943. The other two are ClinVar records 551063 and 1167157, which insert three bases and two at the same point. For `NM_015120.4` the RefSeq GRCh38 cache of release 116 holds the edit `186 185` with three bases, a coding sequence that starts at cDNA 112, and two mapper pairs that meet at the point:

| cDNA | Genome |
| --- | --- |
| 1-185 | 2:73385758-73385942 |
| 189-438 | 2:73385943-73386192 |

Both engines run `--refseq --hgvs` with a reference FASTA. No row of the three records prints `HGVS_OFFSET`.

### 1.3 What each engine writes

The first record's row on the model, in the default format. Ensembl VEP 116.2:

```
seam_one_base	2:73385942-73385943	C	7840	NM_015120.4	Transcript	frameshift_variant	186-188	75-77	25-26	EE/DX	gaagag/gaCa	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=3;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK;HGVSc=NM_015120.4:c.77_78insC;HGVSp=NP_055935.4:p.Glu25AspfsTer102
```

vep-rs:

```
seam_one_base	2:73385942-73385943	C	7840	NM_015120.4	Transcript	frameshift_variant	186-188	75-77	25-26	EE/DX	gaggaa/gaCa	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=3;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK;HGVSc=NM_015120.4:c.75_77delinsC;HGVSp=NP_055935.4:p.Glu25AspfsTer102
```

The two rows agree in every column but the `HGVSc` and the reference half of `Codons`. This entry is the `HGVSc`. The codons are the subject of the entry [`refseq_edited_model_inserted_span_codons`](refseq-edited-model-inserted-span-codons.md). The `HGVSc` of the three records on the model:

| Record | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- |
| `seam_one_base`, one base | `c.77_78insC` | `c.75_77delinsC` |
| 551063, three bases | `c.77_78insGGC` | `c.77A>C` |
| 1167157, two bases | `c.77_78insGG` | `c.78del` |

All six rows print `CDS_position` 75-77. The two rows of record 551063 differ in their term as well: Ensembl VEP writes `intergenic_variant` and vep-rs `missense_variant`. That is the subject of the entry [`refseq_edited_model_equal_span_insertion_term`](refseq-edited-model-equal-span-insertion-term.md).

### 1.4 Why Ensembl VEP writes that

The steps on the first record's row, by line of `TranscriptVariationAllele.pm`:

| Step | Result |
| --- | --- |
| Lines 1474-1483 build the notation on the genomic slice, which starts at the model's first base, 2:73385758 | the insertion's two slice positions are 73385942 - 73385758 + 1 = 185, and 186 |
| `hgvs_variant_notation` types the allele between those two neighbours | an `ins` |
| The edit starts at cDNA 186, the row's cDNA start: line 1610 does not leave it out, and line 1617 adds its length | 3, the row's `REFSEQ_OFFSET` |
| Lines 1524-1525 move both positions by that sum | 188 and 189 |
| The genomic bases at those slice positions | 73385758 + 188 - 1 = 73385945, and 73385946 |
| Both lie in the model's first exon, which starts at cDNA 1: counted from the exon's start | cDNA 1 + (73385945 - 73385758) = 188, and 189 |
| Counted from the coding start, cDNA 112 | 188 - 111 = 77 and 189 - 111 = 78: `c.77_78insC` |

On the model's own sequence those two genomic bases are not 77 and 78. The mapper pair that starts at 73385943 with cDNA 189 places 73385945 at cDNA 189 + 2 = 191, coding position 191 - 111 = 80, and 73385946 at 81. Coding position 77 is the last of the three added bases, which no genomic base holds. The notation pairs two genomic bases, three positions past the insertion point, with numbers that belong in part to the added bases.

### 1.5 Why vep-rs's row is the correct one

The record's two flanking bases are 2:73385942 and 2:73385943. The first mapper pair ends at 73385942 with cDNA 185, and the second starts at 73385943 with cDNA 189: coding positions 185 - 111 = 74 and 189 - 111 = 78. Between them lie cDNA 186 to 188, coding positions 75 to 77, three bases for which the genome has no position. Both engines print that span as the row's `CDS_position`.

On the model, then, the record's sequence is the model's own with coding positions 75 to 77 exchanged for the allele. The row's other columns say so in both engines. The alternate half of `Codons`, `gaCa`, holds coding positions 73 and 74, the inserted base and position 78. That is four bases where the row's two codons hold six, 6 - 3 + 1 = 4: the three added bases are gone and the inserted base stands in their place.

`c.75_77delinsC` states that exchange. The inserted base is neither the first nor the last of the model's three, so the clip takes nothing and the change spans all three. `c.77_78insC` states a sequence that keeps all three added bases, three bases longer than the one the row's other columns describe.

The two ClinVar records show the clip and the 3' placement. Record 551063 inserts three bases. The first two equal the model's own at 75 and 76 and are clipped, which leaves the third against the base at 77: `c.77A>C`. Record 1167157 inserts two bases, the model's own at 75 and 76, so the clip leaves the deletion of the base at 77. The base at 78 is the same as the base at 77, and the 3' rule on the model writes `c.78del`. Under `--shift_hgvs 0` the deletion stays at 77. Ensembl VEP writes `c.77_78ins` and the whole allele for both.

### 1.6 Try it yourself

With the records of 1.2 saved as `record.vcf`, vep-rs runs them as:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache <cache>/homo_sapiens_refseq/116_GRCh38 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --species homo_sapiens --assembly GRCh38 --refseq --hgvs
grep NM_015120.4 out.txt | grep -o 'HGVSc=[^;]*'
```

`<cache>` is Ensembl VEP's RefSeq GRCh38 cache of release 116, converted for vep-rs ([Cache conversion](../cache-conversion.md)). The command prints the three values of the vep-rs column of 1.3. Ensembl VEP 116.2 takes the same flags with the cache options of [Reproducing these examples](reproducing.md) and its own RefSeq cache, and prints `c.77_78ins` with each allele.

The entry's reproduction is one run, of these three records. No committed cache holds a model with such bases inside its coding sequence, and no committed FASTA holds chromosome 2. The reproduction therefore converts the one region file of Ensembl VEP's cache that holds the model, `2/73000001-74000000`. It also writes a FASTA that holds the 435 genomic bases of the model's first exon and `N` elsewhere.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) holds the test of this class, the function `edited_model_seam_insertion_replacement`. It applies the test where a comparison names the class, which no comparison case or corpus manifest does. A row of the shape is therefore reported as a difference in `HGVSc`.

Where a comparison names the class, the test sets a differing `HGVSc` aside when the row's own fields show the shape. It has four conditions:

| Condition on | What must hold |
| --- | --- |
| 1. The row | it is an insertion's and carries `BAM_EDIT`, `OK` or `FAILED`, in both engines |
| 2. Ensembl VEP's value | it is an `ins` of the row's allele, on the transcript's strand, in `c.` numbering, at any pair of coordinates |
| 3. vep-rs's position and codons | its `CDS_position` names a span inside the codons that its `Protein_position` and `Codons` give, and the alternate half of its `Codons` is the reference half with the span's bases replaced by the allele |
| 4. vep-rs's notation | under the same reference name, it is a change that turns the reference half of its `Codons` into the alternate half, or it is absent where the allele is the model's own bases |

The third condition tells such a span from an ordinary insertion, whose alternate half is longer by the whole allele. The fourth judges the notation by what it says, not by its text, because the 3' placement reads bases of the model that the row does not carry. `c.75_77delinsC` on `gaggaa` gives `gaCa`. `c.78del`, one position past the span and inside the row's codons, turns `gaggaa` into that row's `gaGGa`.

Everything else stays a difference. The comparator's own tests (`scripts/concordance/test_compare_vep_formats.py`) give the class's test each of these and hold it refused:

- on a row of the shape, a vep-rs notation that gives another sequence, Ensembl VEP's own `ins` among them, or an `ins` between two positions that are not neighbours;
- any notation where the allele is the model's own bases, and a notation under another reference name;
- a row whose Ensembl VEP value is not an `ins` of the allele, a row without `BAM_EDIT`, and any row off such a span;
- a notation that names a position outside the row's codons, except an `ins` between two neighbours of which one lies one base outside them;
- a row in an untranslated region or on a non-coding model, which prints no `CDS_position`, `Protein_position` or `Codons` for the test to read.

On the 884 paired rows of the set of 51 insertions, `HGVSc` differs on 202. The test accepts 108, every row on a nuclear model's added bases, and refuses 94: the 6 mitochondrial rows, which carry no `BAM_EDIT`, and 88 rows beside such a place or shifted onto one. The test accepts both ClinVar rows of 1.3. A row beside such a place, or moved onto one by the shift, is off such a span: the moved row prints the written point's positions, as section 3 shows. How many of the 88 carry vep-rs's notation of section 3 is not counted here.

## 2. Three ClinVar insertions on the bases a model holds in an untranslated region

**In short.** ClinVar records 1049039, 212026 and 335960 on GRCh37 are insertions written exactly where a RefSeq model holds, in an untranslated region, bases the genome lacks. Both engines print the span of those bases as the row's `cDNA_position`, and the row has no `CDS_position`. Under `--refseq --hgvs`, Ensembl VEP 116.2 writes an `ins` of every inserted base, between two positions numbered as if the exon held no added base. vep-rs writes the change the inserted bases make to the model's own, a deletion, no notation and a duplication; `HGVSc` alone differs, and the comparison sets nothing aside for these rows.

### 2.1 What a reader needs to know first

Section 1.1 holds here. A row in an untranslated region prints no `CDS_position`, `Codons` or `Amino_acids`, and its term is the region's, so the span of the model's added bases shows in `cDNA_position` alone. An HGVSc position there is counted from the coding sequence: `c.*97` is 97 bases past the stop codon, and `c.-97` is 97 bases before the start codon.

### 2.2 The records

Three records of the ClinVar GRCh37 file, each adding bases between two neighbouring genomic bases:

| Record | Between | Model, strand | Inserted, region |
| --- | --- | --- | --- |
| 1049039 | 6:57512788 and 6:57512789 | `NM_000947.5`, forward | 2 bases, 3' untranslated |
| 212026 | 1:145511796 and 1:145511797 | `NM_005105.5`, forward | 1 base, 3' untranslated |
| 335960 | 2:38303017 and 2:38303018 | `NM_000104.4`, reverse | 2 bases, 5' untranslated |

Each model prints `BAM_EDIT=OK` and holds at that point the bases of one insertion edit: 11, 1 and 1. Both engines run `--refseq --hgvs`, and no row prints `HGVS_OFFSET`. Sections 2.4 and 2.5 work the first record through. The RefSeq GRCh37 cache of release 116 holds for its model, `NM_000947.5`:

| Of the model | Value |
| --- | --- |
| Genomic span | 6:57182428-57513376 |
| Coding sequence | cDNA 82-1611 |
| Insertion edit | `1698 1697`, 11 bases |
| First base of the last exon | 6:57512472, cDNA 1381 |
| Mapper pair before the point | cDNA 1381-1697 on 6:57512472-57512788 |
| Mapper pair after the point | cDNA 1709-1713 on 6:57512789-57512793 |

### 2.3 What each engine writes

Of the 18 columns and fields compared on each row, `HGVSc` alone differs. Both engines print the same term, the region's, which is `3_prime_UTR_variant` on the first row. Both print the span as `cDNA_position`, no `CDS_position` and the same `REFSEQ_OFFSET`: 11, 1 and 1. The `HGVSc` of each row, by its kind and positions:

| Record | `cDNA_position` | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- | --- |
| 1049039 | 1698-1708 | `ins` of 2 bases between `c.*97` and `c.*98` | `del` of `c.*87` to `c.*95` |
| 212026 | 3169 | `ins` of 1 base between `c.*2586` and `c.*2587` | no notation |
| 335960 | 277 | `ins` of 2 bases between `c.-97` and `c.-96` | `dup` of `c.-97` |

### 2.4 Why Ensembl VEP writes that

The cause is that of section 1, and the steps are those of 1.4. On the first record:

| Step | Result |
| --- | --- |
| The point's two slice positions | 57512788 - 57182428 + 1 = 330361, and 330362 |
| The edit starts at the row's cDNA start, 1698: line 1610 keeps it, line 1617 adds its length | 11, the row's `REFSEQ_OFFSET` |
| Lines 1524-1525 move both positions | 330372 and 330373 |
| The genomic bases at those positions | 6:57512799 and 6:57512800 |
| Counted from the first base of the last exon | 1381 + (57512799 - 57512472) = 1708, and 1709 |
| Lines 2848-2852, past the stop codon at cDNA 1611 | `c.*97` and `c.*98` |

On the model's own sequence `c.*97` is cDNA 1708, the last of the 11 added bases, and `c.*98` is cDNA 1709, the record's neighbour after the point. The two genomic bases the numbers were counted from, 6:57512799 and 6:57512800, lie past both mapper pairs of 2.2. The model's mapper numbers them cDNA 1716 and 1717, `c.*105` and `c.*106`, in pairs that 2.2 does not list. On the third record lines 2865 and 2867 count both positions before the start codon as negative numbers.

On all three rows the `ins` stands between exactly the two numbers that the slide by the row's `REFSEQ_OFFSET` gives. Read on the model, it gives another sequence than the record's.

### 2.5 Why vep-rs's row is the correct one

The first record's flanking bases are cDNA 1697 and 1709. Between them lie cDNA 1698 to 1708, the 11 bases the genome has no position for, 87 to 97 bases past the stop codon.

Both flanking bases have a position on the model, with the model's added bases between them. On either side of the point the model carries the genome's sequence. The two mapper pairs of 2.2 hold 322 bases: the model's sequence is the genome's on 313 of them, and the other 9 stand at substitution edits of the model. On the other two models the mapper pairs on either side of the point hold 4,371 and 371 bases, with 1 and 0 at such edits.

So whatever sequence the record gives on the genome, on the model it is the model's own with the bases between the two flanking bases exchanged for the inserted ones. An HGVSc is a statement on the model's sequence, and the row's notation is the description of that exchange.

The 2 inserted bases share no leading base with the model's 11 and equal its last 2. The clip leaves the deletion of the first 9, cDNA 1698 to 1706, `c.*87` to `c.*95`, which the 3' rule does not move. The lengths tell the two notations apart:

| Sequence on the model | Bases |
| --- | --- |
| The model's cDNA | 2,303 |
| The record's, which the deletion of 9 gives | 2303 - 11 + 2 = 2,294 |
| Ensembl VEP's `ins` | 2303 + 2 = 2,305 |

On the second record the inserted base is the model's own. The record's sequence on the model is the model's itself, and vep-rs writes no notation where Ensembl VEP names an insertion. On the third the model holds one added base at the point, and the 2 inserted bases are that base twice. vep-rs writes the duplication of the model's own base, a `dup` of `c.-97`. Each value is the smallest change, placed 3'-most, that the inserted bases make in place of the model's own.

### 2.6 Try it yourself

The records are taken from the ClinVar GRCh37 file by their numbers. The `vep` command of 1.6 runs them with `--assembly GRCh37`, the conversion of Ensembl VEP's RefSeq GRCh37 cache of release 116 and a GRCh37 FASTA. The rows of `NM_000947.5`, `NM_005105.5` and `NM_000104.4` then carry, in each engine, the `HGVSc` that 2.3 describes. This form has no reproduction of its own: the entry's reproduction is one run, on GRCh38, of the records of section 1, and it does not run these three.

### 2.7 How the comparator treats this class

Nothing is set aside for these rows: they are outside the class. As 1.7 says, the test refuses a row in an untranslated region: it reads vep-rs's `CDS_position` and `Protein_position` as numbers, and both are `-` on these rows. The rows meet the first two conditions of 1.7. Each is an insertion's row with `BAM_EDIT` in both engines, and Ensembl VEP's value is an `ins` of the row's allele in `c.` numbering.

The format comparator was run on the three rows of each engine three times: applying no class, applying this class, and applying seven edited-model classes together. The seven are this class and the six `refseq_edited_model_` classes that 3.7 lists. The comparator's `--masks` option names the classes a run applies. Each run reports `HGVSc` as a difference on 3 of 3 rows, and no class sets a row aside.

## 3. A ClinVar insertion that the 3' shift moves onto the bases five models hold

**In short.** The 3' shift moves ClinVar record 768677, 18 bases inserted on GRCh37, onto 9 bases that five RefSeq models hold and the genome lacks. Ensembl VEP 116.2 writes an `ins` of the 18, vep-rs a `dup` of the 9, and nothing is set aside.

### 3.1 What a reader needs to know first

Section 1.1 holds here. A 3' shift rotates an insertion's bases by the distance it moves them. The entry [`refseq_edited_model_shifted_insertion_bases`](refseq-edited-model-shifted-insertion-bases.md) describes the rotation, and that Ensembl VEP writes the bases unrotated on a model with an `_rna_edit`.

### 3.2 The record

| Of the record | Value |
| --- | --- |
| Record | 768677 of the ClinVar GRCh37 file |
| Inserted | 18 bases, one string of 9 twice |
| Between | 14:94582130 and 14:94582131 |
| The genome lets it move 3' | by 3 bases and no further |
| Rows whose insertion the shift moves onto such bases | 5 of the record's 6, one on each model below |

| Of the five models | Value |
| --- | --- |
| The models | `NM_001130080.3`, `NM_001288952.2`, `NM_001288956.2`, `NM_001366993.1`, `NM_001366994.1` |
| Strand and `BAM_EDIT` | forward, `OK` |
| Insertion edit of each, its one `_rna_edit` | 9 bases after 14:94582133: coding 129 to 137 |
| Coding start of `NM_001130080.3` | cDNA 121 |
| Its mapper pair before the edit | cDNA 242-248 on 14:94582127-94582133 |
| Its mapper pair after the edit | cDNA 258-403 on 14:94582134-94582279 |

Both engines run `--refseq --hgvs`.

### 3.3 What each engine writes

| On each of the 5 rows | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- |
| `CDS_position`, the two coding positions the written point lies between | 125-126 | 125-126 |
| `Protein_position` | 42 | 42 |
| `Codons` | 3 bases against 21 | 3 bases against 21 |
| Term | `inframe_insertion` | `inframe_insertion` |
| `REFSEQ_OFFSET` | none | none |
| `HGVS_OFFSET`, the bases the 3' shift moved the insertion | 3 | 3 |
| `HGVSc` | `ins` of 18 bases between `c.128` and `c.129` | `dup` of `c.129` to `c.137` |

Their `HGVSp` differs too: Ensembl VEP writes a `delins` of residues 44 to 46 and vep-rs a `dup` of them. This entry and this section judge the `HGVSc` alone.

With `--shift_hgvs 0`, vep-rs writes on these rows an `ins` of 18 bases between `c.125` and `c.126`, the written point, and no `HGVS_OFFSET`. This page does not show Ensembl VEP's rows under that flag.

### 3.4 Why Ensembl VEP writes that

| Step, by line of `TranscriptVariationAllele.pm` | Result |
| --- | --- |
| As in section 1, the notation is typed on the genome and numbered as if the exon held no added base | here at the shifted point |
| Line 1441 moves the insertion 3' along the genome, by 3 | between 14:94582133 and 14:94582134, where each model holds its 9 bases |
| Lines 1474-1483 type it on the genomic slice, where `Sequence.pm`, lines 576-588, finds no copy of the 18 bases beside the point and lines 591-594 set the type | an `ins` between two neighbours |
| Lines 1510-1513 write the record's bases back as the record spells them, not rotated by the 3 bases of the shift | all 18 bases |
| Line 1610: the row's cDNA start is the written point's, before the edit's start, so the edit is left out | nothing is slid, and no `REFSEQ_OFFSET` |
| Lines 2766-2794 count each position from the start of its exon, as if the exon held no added base, and lines 2862-2867 from the coding start. On `NM_001130080.3`, 14:94582133 is cDNA 248, coding 248 - 120 = 128, and 14:94582134 is counted as the next number | `c.128` and `c.129` |
| On the model's own sequence | `c.129` is the first of the 9 added bases, and 14:94582134 is cDNA 258, `c.138` |

### 3.5 Why vep-rs's row is the correct one

| On each of the five models | Value |
| --- | --- |
| The model's sequence | the genome's exons with the 9 added bases once |
| The record's 18 bases, rotated by the 3 bases of the shift | the 9 added bases twice |
| The record's sequence over the exons | the genome's exons with those 9 bases twice |
| Smallest change from the model's to the record's, placed 3'-most | a duplication of `c.129` to `c.137` |
| vep-rs, at the shifted point | finds the model's 9 bases between the two shifted neighbours and names their replacement by the 18 rotated bases, clipped and placed 3'-most |
| vep-rs's `dup` | the model's sequence and 9 bases: the record's, on 5 of 5 |
| Ensembl VEP's `ins` | the model's sequence and 18 bases: another, on 5 of 5 |

### 3.6 Try it yourself

The record is taken by its number and run as 2.6 says; the five models' rows then carry the `HGVSc` of 3.3. This form has no reproduction of its own: the entry's reproduction is one run, on GRCh38, and does not run this record.

### 3.7 How the comparator treats this class

Nothing is set aside for these rows: their `CDS_position`, 125-126, is the written point's, and 1.7 keeps any row off such a span a difference.

| Of the comparison and the reach | Value |
| --- | --- |
| The class | does not hold these rows: their `CDS_position` is no span of bases the model holds and the genome lacks |
| The comparator applying this class alone | `HGVSc` a difference on 5 of 5 rows |
| Applying each of the seven other classes named below alone, and all eight together | `HGVSc` a difference on 5 of 5 rows |
| The rows' `HGVSp` in those runs | a difference on 5 of 5 rows |
| Rows of the concordance datasets with this form | 0: the caches they are annotated on, one for each assembly, hold no model with such bases |
| Rows shown with it, under `--refseq --hgvs` on the RefSeq GRCh37 cache | these 5, of one record of the ClinVar GRCh37 file |
| Other insertions of the ClinVar GRCh37 file that the shift alone moves onto such bases | not counted |

The seven other classes are `start_cooccurrence_swap` and six edited-model classes: `refseq_edited_model_shifted_insertion_bases`, `refseq_edited_model_hgvsc_offset`, `refseq_edited_model_hgvs_offset`, `refseq_edited_model_tail_term`, `refseq_edited_model_deleted_span_term` and `refseq_edited_model_deleted_span_codons`.
