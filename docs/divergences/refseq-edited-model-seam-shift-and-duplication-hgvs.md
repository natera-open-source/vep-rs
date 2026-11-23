# An insertion the 3' shift carries off an edited RefSeq model's added bases, and a duplication across them

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP builds the notation of an insertion or a duplication on the genome, after the genome's own 3' shift, and never reads the bases an edited RefSeq model holds at that place and the genome lacks. The row's position columns are read on the model, so the notation and the columns of one row describe two different sequences. Ensembl VEP's notation keeps those bases; vep-rs's replaces them by the allele.

**Why this is a defect.** An HGVSc is a statement on the sequence of the transcript it names, and the HGVSp follows from it. Where a model holds added bases between two genomic neighbours, both engines print a span that is those bases, or holds them, as the row's `cDNA_position` and `CDS_position`. The row's columns read the allele in place of that span. Ensembl VEP's notation keeps every added base and adds the allele beside them, or names a duplicated range that runs into them. It states a sequence longer than the row's own by the length of the added bases.

**What vep-rs does instead.** vep-rs reads such a row where its record is written and describes the change on the model's own sequence. Three forms of record meet the added bases, and the three paragraphs headed in bold that follow take them in turn: the first and third are the class's two shapes, the second a row beside them.

**An insertion the 3' shift moves off those bases.** The genome holds no base for the model's added ones, so an allele that starts with the base after them shifts past them along the genome. Ensembl VEP 116.3 writes the notation at the shifted point: `HGVS_OFFSET`, an `ins` or a `dup` that keeps the added bases, and an `HGVSp` read there, beside a `CDS_position` that is those bases. vep-rs reads such a row where its record is written. Its `HGVSc` names the replacement of the added bases by the allele, its `HGVSp` is read from the same span, and it prints no `HGVS_OFFSET`.

Under `--shift_3prime 1`, the flag that moves an insertion or deletion 3' before its terms and positions are computed, vep-rs annotates such a row at the shifted point and writes that point's notation.

**A duplication across those bases.** An insertion one base past the added bases can copy, on the genome, the base before them and the base after. The genome shows a duplication of two neighbours; on the model those two bases lie the added bases' length apart. vep-rs numbers the duplication on the model. Where the model's bases beside the point equal the allele, the range is those bases, as many positions as bases copied: `c.77_78dup` on `NM_015120.4`, the range Ensembl VEP 116.3 writes on a forward-strand model. Where they do not, the allele is no duplication on the model and is an `ins` at the point; no row counted on this page is of that kind.

**A duplication spelt as a doubled reference across those bases.** A record can spell that copy another way. Its reference is the two or three genome bases that lie on both sides of the added bases, its first alternate is that reference twice, and a second alternate keeps the reference whole. Ensembl VEP 116.3 writes a `dup` of as many positions as the genome copies, `c.74_75dup` on `NM_015120.4`. vep-rs names the change the alternate makes to the model's span, after the ends the two share: `c.75_76delinsA`.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3 and `ensembl-variation` at release/116, the modules the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1746, and `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 263-264: the row's `CDS_position` is written first, from the span the mapper gives the record. Line 1760 of `OutputFactory.pm` makes the notation afterwards, by `hgvs_transcript`.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1441, 142 and 422-423: `hgvs_transcript` calls `_return_3prime`, which calls `_genomic_shift`, and the shift reads the genome's slice. Line 388, in `create_shift_hash`, keeps the places moved as `_hgvs_offset`.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1446 and 1474-1478: `hgvs_variant_notation` types the allele on the genomic slice, at the slice positions plus that offset.
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 584-585 and 592-594: there an insertion is a `dup` of the slice's bases before it, or an `ins` between the two slice positions. Lines 526, 535-536 and 602-606: an alternate that is the slice's reference twice is typed `dup`, with the two ends of the record's span as its positions.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 123 and 1492: no 3' shift moves a record that is neither an insertion nor a deletion, and a `dup` skips the clip of shared ends.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1524-1525, 1608, 1610 and 2790: both positions are moved by the offset of the model's edits and numbered by `_get_cDNA_position`, each counted along its exon from its genomic distance. An edit that starts after the row's first cDNA position adds nothing to that offset. No line between 1441 and 1525 reads the bases the model holds at the place, and none compares the shifted point with the span of the row's columns.
- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 306 and 308, and `Bio/EnsEMBL/VEP/Parser.pm`, lines 881 and 978: a record's first base is chopped only when every allele starts with it, and only a record of one alternate is minimised to an insertion at a point. A second alternate that starts with another base therefore keeps a doubled reference whole.
- `Bio/EnsEMBL/Variation/TranscriptVariation.pm`, lines 457 and 463: the reference a row is analysed against is cut from the model's own sequence over the record's span on the model, the added bases included.

**What correcting it would change.** The notation would be built on the model's sequence at the row's own span. A row of either shape would then carry an `HGVSc` and an `HGVSp` that state the exchange its `CDS_position` states, and no `HGVS_OFFSET` where the notation is read at the written place. No other column changes.

**The class it produces.** `refseq_edited_model_seam_shift_and_duplication_hgvs`: a `Transcript` row of one of two shapes on bases the model holds and the genome lacks. The first is an insertion written on those bases that the genome's 3' shift moves: the engines differ in `HGVSc`, `HGVSp` and `HGVS_OFFSET`. The second is a reference doubled across those bases: the engines differ in `HGVSc`. No comparison tool sets a row of either shape aside: a disagreement of these shapes is reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the fields.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a record written where a RefSeq model holds bases the genome lacks, in a run under `--hgvs`. The committed caches hold 49 models with such bases, in 7 corpora, and seven committed flag sets run a RefSeq or merged cache under `--hgvs`. No committed corpus and no comparison case holds a row of either shape. Eight committed corpora hold a RefSeq model with an edit of any kind, and none of their 791 records starts within 60 bases of a place where a model holds added bases.

A set of 73 hand-written records stands at four such places of the RefSeq GRCh38 cache, which 18 models hold: forward-strand `NM_015120.4` at 2:73385942 and `NM_002111.8` at 4:3074935, reverse-strand `NM_138352.3` at 19:14090084 and 15 more at 1:16567174. Those 15 hold a second such place, 1:16567177, three bases apart.

Among the 73, insertions of one, two, three and four bases that start with the genome's base after the added bases give 72 rows on those models under `--refseq --hgvs`. Both engines print the span of the added bases in the three position columns on all 72. Ensembl VEP prints `HGVS_OFFSET` on all 72, beside an `ins` on 54 and a `dup` of one base on 18. vep-rs writes a `delins` on 69 and a deletion on 3, and no `HGVS_OFFSET` on any. Its `HGVSp` is a frameshift on 54, a substitution of one residue on 16 and a `delins` on 2.

On the 8 rows of the two forward-strand models, Ensembl VEP's positions are the model's own, after the added bases. On the 64 rows of the 16 reverse-strand models they carry intron offsets, the numbering of the entry [`refseq_edited_model_hgvsc_offset`](refseq-edited-model-hgvsc.md).

The same set holds a copy of the genome's two bases beside the place, and at 1:16567174 alone a copy of three, written on the added bases: 33 rows, 18 and 15. The genome's shift moves them one place as well, so Ensembl VEP prints `HGVS_OFFSET` and the shifted point's notation on them. vep-rs writes a `delins` and no `HGVS_OFFSET` on all 33. The two-base copy written one base on gives 18 rows. vep-rs writes a `dup` of two positions on all 18, equal to Ensembl VEP's value on the 2 forward-strand rows: `c.77_78dup` on `NM_015120.4` and `c.116_117dup` on `NM_002111.8`. On the 16 reverse-strand rows Ensembl VEP's range carries intron offsets, as on the 64 rows above.

A second set of 60 records at the same four places holds 16 references of two and of three genome bases doubled across the added bases: 72 rows on the 18 models. Ensembl VEP writes a `dup` of as many positions as bases copied on all 72. On the 12 rows of the three models that hold a place alone, the row prints no `REFSEQ_OFFSET`, and Ensembl VEP's positions are plain and start at the row's first coding position. On the 60 rows of the 15 models at 1:16567174, the row prints a `REFSEQ_OFFSET` and the positions carry intron offsets. vep-rs writes a `delins` on 70 and a change of one base on 2.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Two bases inserted on the three bases `NM_015120.4` holds and the genome lacks, moved one place by the genome's shift

**In short.** `2 73385942 shifted_two_bases A AAC` on GRCh38 inserts two bases between 2:73385942 and 2:73385943. The forward-strand RefSeq model `NM_015120.4` holds three bases at that point that the genome lacks, its coding positions 75 to 77. The first inserted base equals the genome's base at 2:73385943, so the genome's 3' shift moves the insertion one place. Both engines print `CDS_position` 75-77 and `frameshift_variant`. Under `--hgvs`, Ensembl VEP 116.3 writes `HGVS_OFFSET=1`, `c.78_79insAC` and a frameshift named at residue 27. vep-rs writes `c.75_77delinsAC`, a frameshift named at residue 26 and no `HGVS_OFFSET`.

### 1.1 What a reader needs to know first

A RefSeq cache carries each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). Where the alignment found bases in the model that the genome lacks, the model holds them between two bases that are neighbours on the genome. The entry [`refseq_edited_model_seam_insertion_hgvsc`](refseq-edited-model-seam-insertion-hgvsc.md) describes the edit attribute that names them. An insertion written between those two neighbours has the added bases as its span on the model, and the row prints them as its `cDNA_position` and `CDS_position`. `REFSEQ_OFFSET` prints the summed length change of the model's edits that start at or before the row's first cDNA position.

Under `--hgvs` an insertion is shifted 3' as far as its bases repeat before its notation is written, and `HGVS_OFFSET` prints the places moved. An inserted allele whose first base equals the genomic base after the insertion point can be moved one place: the allele's remaining bases followed by that base, inserted one base on, give the same genome.

HGVS describes a change on the transcript's own sequence. An `ins` names two neighbouring positions and the bases put between them, and a `dup` names a range that is copied; both keep every reference base. A `delins` names a span and the bases that replace it.

### 1.2 The record

```
2	73385942	shifted_two_bases	A	AAC	.	.	.
2	73385942	copy_on_added_bases	A	AAA	.	.	.
2	73385943	copy_one_base_on	A	AAA	.	.	.
```

All three records are written for this example. The first inserts `AC` between 2:73385942 and 2:73385943. The genome's bases at 2:73385943 and 2:73385944 are `A` and `G`: the first inserted base equals the one and the second differs from the other, so the genome's shift moves the insertion one place and stops. The second record inserts `AA` at the same point, a copy of the genome's bases at 2:73385942 and 2:73385943, which the shift moves one place as well. The third writes that copy one base on, between 2:73385943 and 2:73385944.

The RefSeq GRCh38 cache of release 116 holds, for `NM_015120.4`, a coding sequence that starts at cDNA 112 and two mapper pairs that meet at the insertion point of the first two records:

| cDNA | Genome |
| --- | --- |
| 1-185 | 2:73385758-73385942 |
| 189-438 | 2:73385943-73386192 |

Both engines run `--refseq --hgvs` with a reference FASTA.

### 1.3 What each engine writes

The first record's row on the model, in the default format. Ensembl VEP 116.3:

```
shifted_two_bases	2:73385942-73385943	AC	7840	NM_015120.4	Transcript	frameshift_variant	186-188	75-77	25-26	EE/EX	gaagag/gaACa	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=3;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK;HGVSc=NM_015120.4:c.78_79insAC;HGVSp=NP_055935.4:p.Glu27GlnfsTer16;HGVS_OFFSET=1
```

vep-rs:

```
shifted_two_bases	2:73385942-73385943	AC	7840	NM_015120.4	Transcript	frameshift_variant	186-188	75-77	25-26	EE/EX	gaggaa/gaACa	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=3;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK;HGVSc=NM_015120.4:c.75_77delinsAC;HGVSp=NP_055935.4:p.Glu26GlnfsTer16
```

Every column but the HGVSc, the HGVSp, `HGVS_OFFSET` and the reference half of `Codons` agrees. The reference codons differ as the entry [`refseq_edited_model_inserted_span_codons`](refseq-edited-model-inserted-span-codons.md) describes, no part of this entry. Ensembl VEP's `ins` spells the two bases in the record's order, the subject of the entry [`refseq_edited_model_shifted_insertion_bases`](refseq-edited-model-shifted-insertion-bases.md). The three records on the model:

| Record | Field | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `shifted_two_bases` | `HGVSc` | `c.78_79insAC` | `c.75_77delinsAC` |
| `shifted_two_bases` | `HGVSp` | `p.Glu27GlnfsTer16` | `p.Glu26GlnfsTer16` |
| `shifted_two_bases` | `HGVS_OFFSET` | 1 | none |
| `copy_on_added_bases` | `HGVSc` | `c.77_78dup` | `c.75_76delinsA` |
| `copy_on_added_bases` | `HGVSp` | `p.Glu27LysfsTer16` | `p.Glu26LysfsTer16` |
| `copy_on_added_bases` | `HGVS_OFFSET` | 1 | none |
| `copy_one_base_on` | `HGVSc` | `c.77_78dup` | `c.77_78dup` |
| `copy_one_base_on` | `HGVSp` | `p.Glu27LysfsTer16` | `p.Glu27LysfsTer16` |
| `copy_one_base_on` | `HGVS_OFFSET` | none | none |

The first two records print `CDS_position` 75-77 in both engines, and the third prints 78-79. On the third record's row every column agrees.

### 1.4 Why Ensembl VEP writes that

The row's columns are written first. The record's two neighbours, 2:73385942 and 2:73385943, lie in two mapper pairs, at cDNA 185 and 189. The insertion's span runs from one past the first to one before the second, cDNA 186 to 188. From the coding start that is 186 - 111 = 75 to 188 - 111 = 77, the row's `CDS_position`.

The notation is made afterwards, and the insertion is shifted first. The genome has no base for cDNA 186 to 188, so on it the insertion point is followed directly by 2:73385943. That base equals the first inserted base, and the shift moves the point one place, between 2:73385943 and 2:73385944: `HGVS_OFFSET=1`.

`hgvs_variant_notation` then types the allele on the genomic slice at the shifted point. It copies no neighbouring bases there, so it is an `ins` between two slice positions. The slice starts at the model's first base, 2:73385758, so those positions are 73385943 - 73385758 + 1 = 186 and 187.

Lines 1524-1525 move both by the row's `REFSEQ_OFFSET`, 3, to 189 and 190, and count them along the exon as if it held no added base, the routine section 1.4 of the entry [`refseq_edited_model_seam_insertion_hgvsc`](refseq-edited-model-seam-insertion-hgvsc.md) describes. From the coding start that is 189 - 111 = 78 and 190 - 111 = 79: `c.78_79insAC`. The offset equals the three bases the exon count leaves out, so the two numbers are the model's own coding positions for 2:73385943 and 2:73385944. The HGVSp is read at the same point, after the codon that ends at coding position 78, residue 26, so the frameshift is named at residue 27.

### 1.5 Why vep-rs's row is the correct one

The record's flanking bases, 2:73385942 and 2:73385943, are coding positions 185 - 111 = 74 and 189 - 111 = 78. Between them lie coding positions 75 to 77, the three bases the genome has no position for, and both engines print that span as the row's `CDS_position`. Both also print the same alternate codons, `gaACa`: coding positions 73 and 74, the two inserted bases, and position 78. The model's six bases at 73 to 78 are `gaggaa`, vep-rs's reference codons, so the alternate holds 6 - 3 + 2 = 5 bases. The three added bases are gone and the two inserted bases stand in their place.

`c.75_77delinsAC` states that exchange. Neither inserted base is an end the two share: the first, `A`, differs from the model's base at 75, `g`, and the second, `C`, from its base at 77, `a`. `c.78_79insAC` states another sequence: all six bases kept and two added, 6 + 2 = 8 bases over the same codons, where the row's `Codons` column has five.

The HGVSp follows the same reading. In `gaACa` the first codon, `gaA`, codes the same residue as the model's `gag`, residue 25. The next codon starts with the inserted `C` where the model's codon `gaa` starts with `g`, so the first changed residue is 26: `p.Glu26GlnfsTer16`.

The genome's shift does not carry the row's change to an equal one on the model. On the genome, an insertion written at the point and its spelling one base on give one sequence. On the model they give two: at the written point the allele meets the three added bases, and one base on it meets none. The second and third records are such a pair. vep-rs reads each row where its record is written, as every other column of the row is read in both engines: `c.75_76delinsA` beside `CDS_position` 75-77, and `c.77_78dup` beside 78-79. Ensembl VEP writes `c.77_78dup` on both rows.

### 1.6 Try it yourself

With the records of 1.2 saved as `record.vcf`, the command as vep-rs takes it:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache <cache>/homo_sapiens_refseq/116_GRCh38 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --species homo_sapiens --assembly GRCh38 --refseq --hgvs
grep NM_015120.4 out.txt | grep -o 'HGVS[^;]*'
```

`<cache>` is the conversion of Ensembl VEP's RefSeq GRCh38 cache of release 116 ([Cache conversion](../cache-conversion.md)). vep-rs prints the values of its column in 1.3. Ensembl VEP 116.3 takes the same flags with the cache options of [Reproducing these examples](reproducing.md) and its RefSeq cache, and prints the values of its column.

The entry's reproduction prepares its inputs as the reproduction of `refseq_edited_model_seam_insertion_hgvsc` does. It converts the one region file of Ensembl VEP's cache that holds the model, `2/73000001-74000000`, and writes a FASTA of the 435 genomic bases of the model's first exon with `N` elsewhere. Both engines write the same rows on it as on the primary assembly.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`tools/vep-formats-compare`) has no class of this shape, and no committed corpus holds a row of it. A disagreement of the shape is reported as a difference in `HGVSc`, `HGVSp` and `HGVS_OFFSET` and charged to vep-rs.

The neighbouring class `refseq_edited_model_seam_insertion_hgvsc` has a test in the format comparator, the function `edited_model_seam_insertion_replacement`. No comparison case names that class, so no comparison applies the test, to these rows or to any other. The test is no rule for this shape either. It reads the `HGVSc` alone and asks that Ensembl VEP's value be an `ins` of the row's allele, which it is on 54 of the 72 rows of the first set.

Rows beside the shape are no part of this entry. The same allele written one base on, where its point lies between two neighbours of the model, gets one notation from both engines on a forward-strand model.

## 2. Two genome bases doubled across the three bases `NM_015120.4` holds and the genome lacks

**In short.** `2 73385942 doubled_two_bases AA AAAA,CA` on GRCh38 has as its reference the two genome bases on either side of the three bases `NM_015120.4` holds and the genome lacks. Its first alternate is that reference twice: on the genome, a duplication of two bases. For that alternate both engines print `CDS_position` 74-78, `USED_REF=AGGAA` and `frameshift_variant`. Under `--hgvs`, Ensembl VEP 116.3 writes `c.74_75dup`, two positions that start at the row's first and run into the added bases. vep-rs writes `c.75_76delinsA`, the change the alternate's four bases make to the row's five.

### 2.1 What a reader needs to know first

A VCF record whose alleles all start with one base loses that base, and a record of one alternate whose alleles differ in length is cut down to its smallest form, an insertion at a point for a doubled reference. A record with a second alternate that starts with another base is left as written. Its span on a transcript is then the span of its whole reference, from the first base to the last.

Where the reference lies on both sides of a model's added bases, that span on the model holds the reference's bases and the added ones. The row prints the genome's reference as `GIVEN_REF` and the model's bases over the span as `USED_REF`, and its `Codons` set the alternate against the model's bases.

### 2.2 The record

```
2	73385942	doubled_two_bases	AA	AAAA,CA	.	.	.
2	73385942	doubled_three_bases	AAG	AAGAAG,CAG	.	.	.
```

Both records are written for this example. Each has its reference twice as its first alternate, and a second alternate that starts with another base and keeps the reference whole. The first record's reference is the genome's bases at 2:73385942 and 2:73385943, and the second's the three bases 2:73385942 to 2:73385944. The model and its mapper pairs are those of 1.2: cDNA 185 and 189 are the first two genome bases, and cDNA 186 to 188 are the added ones.

### 2.3 What each engine writes

The row of the first record's first alternate on the model, in the default format. Ensembl VEP 116.3:

```
doubled_two_bases	2:73385942-73385943	AAAA	7840	NM_015120.4	Transcript	frameshift_variant	185-189	74-78	25-26	EE/EX	gAGGAA/gAAAA	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=AA;USED_REF=AGGAA;BAM_EDIT=OK;HGVSc=NM_015120.4:c.74_75dup;HGVSp=NP_055935.4:p.Glu26LysfsTer16
```

vep-rs:

```
doubled_two_bases	2:73385942-73385943	AAAA	7840	NM_015120.4	Transcript	frameshift_variant	185-189	74-78	25-26	EE/EX	gAGGAA/gAAAA	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=AA;USED_REF=AGGAA;BAM_EDIT=OK;HGVSc=NM_015120.4:c.75_76delinsA;HGVSp=NP_055935.4:p.Glu26LysfsTer16
```

Every column but the HGVSc agrees. The first alternates of the two records on the model:

| Record | `CDS_position` | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `doubled_two_bases` | 74-78 | `c.74_75dup` | `c.75_76delinsA` |
| `doubled_three_bases` | 74-79 | `c.74_76dup` | `c.75G>A` |

The second record's row prints `USED_REF=AGGAAG`, `synonymous_variant` and an HGVSp that ends in `=`, which the default format writes as `%3D`, in both engines. No row of the two records prints `HGVS_OFFSET`, and the rows of the second alternates agree in every column.

### 2.4 Why Ensembl VEP writes that

The record is neither an insertion nor a deletion, so no 3' shift moves it. `hgvs_variant_notation` reads the type on the genome's slice over the record's genomic span, 2:73385942 to 2:73385943. The alternate is the slice's reference twice, so the type is `dup` and its positions are the span's two ends. A `dup` skips the clip of shared ends and prints its positions and the word alone.

The row's first cDNA position is 185, before the edit that adds the three bases, which starts at 186. Line 1610 leaves such an edit out of the offset, so this row prints no `REFSEQ_OFFSET` and its positions are not moved. The first two rows of section 1 start at 186, on the edit, and are moved by 3.

Each end is numbered along its exon from its genomic distance alone, so the two ends stay as far apart as on the genome: one position. The first is the row's first coding position, 74. The second is 75, which on the model is the first of the three added bases. The second copied base is coding position 78. The result is `c.74_75dup`.

### 2.5 Why vep-rs's row is the correct one

The record's span on the genome is the two copied bases. The mapper pairs number its first and last base three positions further apart than the genome does, cDNA 185 and 189. The span on the model is therefore five bases, which the row's reference is cut from: `CDS_position` 74-78 and `USED_REF=AGGAA` in both engines, beside `GIVEN_REF=AA`.

The row therefore states that the alternate's four bases stand in place of those five. Both engines print that in `Codons`, `gAGGAA/gAAAA`, and both name the frameshift it makes, one base fewer. The two alleles share their first base and their last two. What is left is the model's bases at 75 and 76 replaced by one base: `c.75_76delinsA`.

`c.74_75dup` states another sequence: two bases more than the model's, and three more than the row's own. A duplication of two bases is a frameshift too, so the term does not tell the two readings apart. The HGVSp does. Both engines print `p.Glu26LysfsTer16`, the new stop 16 residues on, as the replaced span gives. `c.74_75dup` applied to the model's sequence gives 17.

The three-base record shows the same on a row whose term tells the readings apart. Its span on the model is six bases, and its alternate is six bases that differ from them at one position: `c.75G>A`, a change of one base that leaves the protein as it is. Both engines print `synonymous_variant` on the row. Ensembl VEP's `c.74_76dup` names a duplication of three coding bases, which is no synonymous change.

The same copy written as an insertion at the point is the second record of 1.2. vep-rs writes one notation for both spellings, `c.75_76delinsA`, beside one HGVSp.

### 2.6 Try it yourself

The command of 1.6, with the records of 2.2 saved as `record.vcf`, prints the values of 2.3 on the model's rows. The entry's reproduction holds the five records of 1.2 and 2.2.

### 2.7 How the comparator treats this class

Nothing is set aside. The format comparator has no class of this shape, and no committed corpus holds a row of it. A disagreement of the shape is reported as a difference in `HGVSc` and charged to vep-rs.

Rows beside the shape are no part of this entry. A doubled reference of two bases with no added bases between or beside them is a `dup` of two positions in both engines, on 18 of 18 rows of the second set. This section's two-base record without its second alternate is cut down by both engines to an insertion at a point, the form of the third record of 1.2. The second set holds it in two spellings at each place, one behind a base the parse chops, and vep-rs writes a `dup` of two positions on all 36 rows.
