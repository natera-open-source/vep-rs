# The HGVSc of an allele that keeps a base it shares with its reference on an edited RefSeq model

**The defect in one sentence.** Under `--hgvs`, on a RefSeq model that carries an `_rna_edit` attribute, Ensembl VEP clips the bases an allele shares with its reference out of the notation and then writes both alleles back whole. The HGVSc keeps the positions and the kind the clip left, and names a base that the record does not change. vep-rs writes the alleles the clip leaves, an `ins` of the inserted bases alone or a substitution of one base by one.

**Why this is a defect.** The clip leaves a notation of the bases that change: its positions and its kind describe those bases alone. The alleles written back are the row's whole ones, so the notation's bases do not fit its positions. An `ins` carries the base its allele shares with the reference and inserts one base more than the record does: the shared base twice. A substitution names two bases a side at a position that holds one. On a model without an `_rna_edit` attribute both engines write, for the same inserted allele, an `ins` of the inserted bases alone.

**What vep-rs does instead.** vep-rs spells the alleles the clip leaves: an `ins` of the inserted bases alone, a substitution of one base by one. Section 1 shows the substitution on two records of a committed corpus, and section 2 the insertion.

Three further forms arise where the record's bases stand at a place where the model differs from the genome. In sections 3 and 4 the genome's alleles, which the clip reads, and the model's reference allele, which is written back, differ as well. Section 3 shows a kept base that the model spells differently, section 4 bases that the model holds and the genome lacks, and section 5 an insertion that the 3' placement carries to an exon's last base.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116; `TranscriptVariationAllele.pm`, `TranscriptVariation.pm` and `Utils/Sequence.pm` are under `Bio/EnsEMBL/Variation/`):

- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines 306-309: in a record with more than one alternate, one of them of another length than the reference, the first base of the alleles is dropped only when every allele starts with it. `Bio/EnsEMBL/VEP/Parser.pm`, lines 860-881, minimises a record with one alternate whose length differs from its reference's, and `minimise_alleles` (lines 964-1019) leaves a record with a second alternate alone at lines 978-984.
- `TranscriptVariationAllele.pm`, lines 1474-1483, `hgvs_transcript`: the notation is first built against the sequence of the transcript's genomic slice. Line 1492 then calls `_clip_alleles` unless the type is `dup`.
- `_clip_alleles`, lines 2217-2221 and 2234-2237: a base the two alleles share at the 5' end, then at the 3' end, is taken off both, and the start or the end follows. Lines 2261-2266 make one base left a side a substitution, and line 2283 makes inserted bases an `ins`. Lines 2245-2248 write the clipped pair and its positions back.
- `TranscriptVariationAllele.pm`, lines 1508-1514: on a model with `_rna_edit` attributes both alleles are then written over from `feature_seq`, the reference allele of the row and the record's alternate, each whole. The positions and the type stay the clip's.
- `Bio/EnsEMBL/Variation/VariationFeatureOverlapAllele.pm`, `feature_seq`, lines 246-262: the allele as stored, on the transcript's strand.
- `TranscriptVariation.pm`, `_get_ref_allele`, line 463: the row's reference allele is cut from the model's sequence over the record's span, and the row prints it as `USED_REF`. `Bio/EnsEMBL/VEP/CacheDir.pm`, lines 412-414, turns that reading on for a cache built with a BAM file.
- `Utils/Sequence.pm`, `format_hgvs_string`, lines 654-656: a substitution prints the start alone and the two alleles as they stand. Lines 666-667 print a `delins` as its coordinates and the alternate.
- `TranscriptVariationAllele.pm`, line 1441, `_return_3prime` (line 123) and `Utils/Sequence.pm`, line 389: an allele that keeps a base has the class `indel`, which the 3' shift does not move.
- `TranscriptVariationAllele.pm`, lines 1516-1517 and 1524-1525: on such a model the two positions are moved by the summed length of the model's insertion edits, and numbered along the genome's exon, as if it held no added base. Line 1610 leaves out of the sum an edit that starts after the row's cDNA start. The entry [`refseq_edited_model_seam_insertion_hgvsc`](refseq-edited-model-seam-insertion-hgvsc.md) describes that routine.

**What correcting it would change.** The alleles written back at lines 1508-1514 would be the row's alleles without the bases the clip took. The notation's alleles would then fit the positions and the kind it keeps, and the rows of section 1 would carry vep-rs's notation.

**The class it produces.** `refseq_edited_model_unminimised_allele_hgvsc`: a differing `HGVSc` on a row that prints `BAM_EDIT=OK`, where the two values are two substitutions, two `ins` or two `delins` at the same coordinates, Ensembl VEP's spelling the row's whole alleles and vep-rs's the same alleles without the bases the clip takes. The format comparator sets a row of the class aside in a comparison that names the class, and counts it once; 1.7 says which rows.

In a comparison that does not name the class, a row of the shape is reported as a difference in `HGVSc` and charged to vep-rs. So is a row of the further forms of sections 3 to 5, which are outside the class. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is an allele that keeps a base it shares with its reference, on a model with an `_rna_edit` attribute, in a run under `--hgvs`. The committed caches hold 160 such models, in 10 corpora, and seven committed flag sets write `HGVSc` on such a cache.

The concordance datasets are annotated on the Ensembl caches, which hold no `_rna_edit` attribute, so no row of them has the shape. No ClinVar record has a second alternate, 0 of 4,398,127 on GRCh38 and 0 of 4,388,172 on GRCh37, so the inserted allele of section 2 comes from no ClinVar record. The rows of sections 2 to 5 are of constructed records on the RefSeq GRCh38 cache of release 116, run under `--refseq --hgvs`. The entry's reproduction runs the two records of section 1 alone.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. One base changed and the base beside it kept, on `NM_004449.4`

**In short.** `21 38575676 kept_base_before CT CC` and `21 38575677 kept_base_after TG CG` on GRCh38 both change the one base 21:38575677, and each keeps a neighbouring base that its reference and its alternate share. The changed base is coding base 24 of the reverse-strand RefSeq model `NM_004449.4`, which carries `_rna_edit` attributes. Under `--refseq --hgvs`, Ensembl VEP 116.3 writes `c.24AG>GG` for the first record and `c.24CA>CG` for the second, each naming two bases a side. vep-rs writes `c.24A>G` for both. The format comparator sets the difference aside in a comparison that names the entry's class.

### 1.1 What a reader needs to know first

A RefSeq cache carries each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). A model that the alignment edited carries `_rna_edit` attributes, and its rows print `BAM_EDIT=OK`. On such a cache a row is read against the model's own sequence: `GIVEN_REF` prints the record's reference allele, and `USED_REF` the model's bases over the record's span.

A record's alleles are trimmed as the record is read. Both engines minimise a record with one alternate allele whose length differs from its reference's: the base the two alleles share is dropped. In a record with two or more alternates, one of another length than the reference, Ensembl VEP's VCF parser drops the first base only when every allele starts with it. An allele therefore keeps a base it shares with its reference in two kinds of record: one whose reference and alternate have the same length and share an end base, and one with two or more alternates.

In HGVS a substitution replaces one nucleotide by one other, and a change of two or more neighbouring nucleotides is a `delins` over their range. An `ins` names two neighbouring positions and the bases put between them.

HGVS describes a change at its 3'-most position on the transcript's sequence, and names inserted bases that repeat the bases before them as a duplication: the 3' rule of this page. The 3' placement is vep-rs's use of that rule on the model's sequence. The 3' shift is the engines' move of an insertion or deletion along the genome before the notation is written, and a moved row prints the distance as `HGVS_OFFSET`.

### 1.2 The record

```
21	38575676	kept_base_before	CT	CC	.	.	.
21	38575677	kept_base_after	TG	CG	.	.	.
```

Both records change the one base 21:38575677. The first is written from 21:38575676 and keeps that base, before the changed one. The second is written from 21:38575677 and keeps the base after it. `NM_004449.4` is a reverse-strand model of the committed `GRCh38-refseq` corpus. It carries `_rna_edit` attributes and prints `BAM_EDIT=OK`. The changed base is its coding base 24, and the model spells both bases of each record as the genome does. Both engines run `--refseq --hgvs` with the corpus's reference FASTA.

### 1.3 What each engine writes

The two records' rows on the model, in the default format. Ensembl VEP 116.3:

```
kept_base_before	21:38575676-38575677	CC	2078	NM_004449.4	Transcript	synonymous_variant	296-297	24-25	8-9	PA	ccAGca/ccGGca	-	IMPACT=LOW;STRAND=-1;GIVEN_REF=CT;USED_REF=CT;BAM_EDIT=OK;HGVSc=NM_004449.4:c.24AG>GG;HGVSp=NP_004440.1:p.ProAla8%3D
kept_base_after	21:38575677-38575678	CG	2078	NM_004449.4	Transcript	synonymous_variant	295-296	23-24	8	P	cCA/cCG	-	IMPACT=LOW;STRAND=-1;GIVEN_REF=TG;USED_REF=TG;BAM_EDIT=OK;HGVSc=NM_004449.4:c.24CA>CG;HGVSp=NP_004440.1:p.Pro8%3D
```

vep-rs:

```
kept_base_before	21:38575676-38575677	CC	2078	NM_004449.4	Transcript	synonymous_variant	296-297	24-25	8-9	PA	ccAGca/ccGGca	-	IMPACT=LOW;STRAND=-1;GIVEN_REF=CT;USED_REF=CT;BAM_EDIT=OK;HGVSc=NM_004449.4:c.24A>G;HGVSp=NP_004440.1:p.ProAla8%3D
kept_base_after	21:38575677-38575678	CG	2078	NM_004449.4	Transcript	synonymous_variant	295-296	23-24	8	P	cCA/cCG	-	IMPACT=LOW;STRAND=-1;GIVEN_REF=TG;USED_REF=TG;BAM_EDIT=OK;HGVSc=NM_004449.4:c.24A>G;HGVSp=NP_004440.1:p.Pro8%3D
```

The rows of each record agree in every column but the `HGVSc`:

| Record | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `kept_base_before` | `c.24AG>GG` | `c.24A>G` |
| `kept_base_after` | `c.24CA>CG` | `c.24A>G` |

`GIVEN_REF` and `USED_REF` are equal on all four rows.

### 1.4 Why Ensembl VEP writes that

The steps on the two rows, by line of `TranscriptVariationAllele.pm`:

| Step | `kept_base_before` | `kept_base_after` |
| --- | --- | --- |
| The row's alleles on the transcript's strand, and their coding positions | `AG` against `GG`, 24-25 | `CA` against `CG`, 23-24 |
| Line 1492, `_clip_alleles`, takes the base the two share off both, and the positions follow | the last base, at 25 (lines 2234-2237) | the first base, at 23 (lines 2217-2221) |
| What the clip leaves | `A` against `G` at 24 | `A` against `G` at 24 |
| Lines 2261-2266: one base is left a side | a substitution | a substitution |
| Lines 1508-1514: the model carries `_rna_edit` attributes, so both alleles are written over with the row's whole alleles | `AG` and `GG` | `CA` and `CG` |
| `Utils/Sequence.pm`, lines 654-656, print the start alone and the two alleles as they stand | `c.24AG>GG` | `c.24CA>CG` |

On a model without an `_rna_edit` attribute lines 1508-1514 write nothing over, and the notation keeps the alleles the clip left.

### 1.5 Why vep-rs's row is the correct one

Each record names two neighbouring genomic bases and changes one of them, 21:38575677. On the model that base is coding position 24, and the record changes it from `A` to `G` on the transcript's strand. The kept base, coding position 25 in the first record and 23 in the second, is the same in the reference and in the alternate. The model spells it as the genome does: `GIVEN_REF` and `USED_REF` are equal on every row.

One nucleotide is replaced by one other, which HGVS writes as a substitution of that nucleotide: `c.24A>G`. vep-rs writes it for both records, which make the same change. Ensembl VEP's `c.24AG>GG` and `c.24CA>CG` each name two reference bases at a position that holds one, and they give one change two notations.

### 1.6 Try it yourself

```
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > record.vcf
printf '21\t38575676\tkept_base_before\tCT\tCC\t.\t.\t.\n' >> record.vcf
printf '21\t38575677\tkept_base_after\tTG\tCG\t.\t.\t.\n' >> record.vcf
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38-refseq/cache_root/homo_sapiens_refseq/116_GRCh38 \
  --fasta tests/golden/116/GRCh38-refseq/reference.fa.gz \
  --species homo_sapiens --assembly GRCh38 --refseq --hgvs
grep NM_004449.4 out.txt | grep -o 'HGVSc=[^;]*'
```

vep-rs prints `HGVSc=NM_004449.4:c.24A>G` twice. Ensembl VEP 116.3 takes the same flags with the cache options of [Reproducing these examples](reproducing.md) and its own RefSeq GRCh38 cache, and prints the two values of its column in 1.3.

The entry's reproduction, `scripts/divergences/reproductions/refseq_edited_model_unminimised_allele_hgvsc.toml`, holds the two records and both engines' values. `tools/vep-divergence-reevaluate` runs it on both engines and holds each value to the recorded one.

### 1.7 How the comparator treats these rows

The format comparator (`tools/vep-formats-compare`) applies `refseq_edited_model_unminimised_allele_hgvsc` when a comparison case names it in `masks`. Its test, the function `edited_model_clipped_alleles`, accepts a differing `HGVSc` on a row that prints `BAM_EDIT=OK` and whose `GIVEN_REF` and `USED_REF` are each the same in both engines. The two values must be two substitutions, two `ins` or two `delins`, under one sequence name and numbering, at the same coordinates. Ensembl VEP's must spell the row's whole alleles on the transcript's strand, and vep-rs's the same alleles without the bases the clip takes. Each accepted row is counted once under the class.

The class accepts both rows of 1.3. The same pair on a row of a model whose edit failed, a notation that vep-rs leaves out and a notation of another allele stay differences. So do two values that differ in kind or in coordinates, which is why the class holds no row of sections 3 to 5. Ensembl VEP's substitution on an allele that runs past the model's last base, whose reference names a base outside the transcript that the clip does not take, is a row of [`hgvs_clipped_reference_against_whole_alternate`](hgvs-clipped-reference-against-whole-alternate.md), section 8 of its page.

## 2. An inserted allele that keeps its reference base

**In short.** A record whose first alternate is its reference base followed by inserted bases, and whose second alternate does not start with that base, keeps the reference base on the inserted allele. On a model with an `_rna_edit` attribute, Ensembl VEP 116.3 writes an `ins` of the inserted bases and the kept base, one base more than the record inserts. Away from bases the model adds, vep-rs writes the `ins` of the inserted bases alone, as both engines do on a model without the attribute. The entry's class asks for two `ins` at the same coordinates; how many of these rows it holds is not counted here.

### 2.1 What each engine writes

A set of 18 constructed records with a second alternate gives 228 rows on which Ensembl VEP writes an `HGVSc` for such an allele. No row of the 228 prints `HGVS_OFFSET`: an allele that keeps a base is not moved by the 3' shift.

| Rows | Ensembl VEP 116.3 |
| --- | --- |
| 214 on models with an `_rna_edit` attribute: 147 in an intron, away from any bases the model adds, and 67 on a model that adds bases at the record's point | an `ins` of the record's inserted length plus one, on 214 of 214 |
| 14 on models without one: 12 that print no `BAM_EDIT` and 2 that print `BAM_EDIT=FAILED` | an `ins` of the record's inserted length, on 14 of 14 |

On `NM_001405680.2`, which prints `BAM_EDIT=OK`, the intron row of a record that inserts three bases carries, in Ensembl VEP, an `ins` of four bases between `c.2922+2` and `c.2922+3`. On `NM_001040152.2`, a model without the attribute, a record that inserts one base carries an `ins` of one base between `c.956` and `c.957`, in both engines.

vep-rs spells the allele the clip leaves, on a model with the attribute as on one without. Away from bases the model adds, it writes the `ins` of the record's inserted bases alone. Where the notation stands on the 67 rows is not shown in this section.

### 2.2 Why vep-rs's row is the correct one

The record's first alternate is its reference base followed by the inserted bases. The reference base stands in both alleles and is not inserted. Among the 67 rows, 59 print `Codons`, and on each of the 59 the codons add exactly the record's inserted length, in each engine. The same allele on a model without an `_rna_edit` attribute is written as an insertion of that length by both engines. A notation that inserts one base more describes a sequence the record does not give: the kept base twice.

## 3. A kept base that the model spells differently from the genome

**In short.** A record's two alleles share an end base, and the model spells that base differently from the genome. Ensembl VEP's clip takes the base off the genome's alleles, and the model's reference allele, written back whole, does not share it with the alternate. Ensembl VEP 116.3 writes a substitution of the model's whole reference allele by the record's whole alternate, two or three bases a side, at the one position the clip leaves. vep-rs writes a `delins` of two bases, or an inversion where those two bases are the reverse complement of the model's.

vep-rs clips a shared base only where the model's own base equals it. It writes a `delins` of two bases, from the base the model spells differently to the base the record changes, or an inversion where those two bases are the reverse complement of the model's. Nothing is set aside.

### 3.1 What each engine writes

One row of a set of 53 constructed records. The set holds records at five places where an edited model differs from the genome by one base and, for section 4, at three places where a model holds bases the genome lacks. The row:

| Of the row | Value |
| --- | --- |
| The record | two reference bases, 1:16057173-16057174, and one alternate of two bases that repeats the first and changes the second |
| The model | `NM_001165945.2`, forward strand, `BAM_EDIT=OK`; `3_prime_UTR_variant`, `cDNA_position` 2021-2022 |
| The model against the genome | another base at cDNA 2021, `c.*257`; the same base at cDNA 2022, `c.*258` |
| `GIVEN_REF` and `USED_REF` | two bases each, and they differ |
| Ensembl VEP 116.3 | a substitution at `c.*258` with two bases a side: the model's bases at `c.*257` and `c.*258` against the record's two |
| vep-rs | the record's two bases in place of the model's two, over `c.*257_*258` |

The set gives 798 rows whose clip takes a base the model spells differently. Ensembl VEP writes the substitution of the whole alleles on 798 of 798, with two bases a side on 266 and three on 532. vep-rs writes a `delins` of two bases on 699 of them and an inversion of two bases on 99. On the 532 rows the clip takes two bases; vep-rs drops a clipped base from the pair only where the model has it, and two bases are left on every row.

Under `--minimal`, where the repeated base is dropped as the record is read, both engines write a substitution of one base at `c.*258` for the record of the table.

### 3.2 Why vep-rs's row is the correct one

The record names two neighbouring genomic bases and changes one. The base it repeats is the genome's, and the model holds another base there. With the record's two bases in place of the model's two, both bases of the model change. That is a change of two neighbouring nucleotides, written over their range, `c.*257_*258`, the span the row's own `cDNA_position` covers.

Ensembl VEP's substitution names two reference bases at a position that holds one. Applied to the model at its printed position, it gives neither that sequence nor the one in which the changed base alone is exchanged. That holds on the 378 of its 798 rows whose position is exonic; the other 420 carry a position with an intron offset, which the entry [`refseq_edited_model_hgvsc_offset`](refseq-edited-model-hgvsc.md) describes.

## 4. A kept base beside bases the model holds and the genome lacks

**In short.** A record's reference covers both genomic neighbours of a place where the model holds bases the genome lacks, and its alleles share an end base. The model's reference allele, `USED_REF`, then holds the model's bases as well. Ensembl VEP 116.3 keeps the position or the span that its clip left on the genome's alleles and writes the whole alleles there. vep-rs writes the change from the first base the genome lacks to the last changed base. Nothing is set aside. On 60 `delins` rows and 12 substitution rows the smallest change on the model is shorter than that span, and which of the two the notation should name is not settled.

### 4.1 What each engine writes

Three rows show the three kinds the clip leaves. The first two are on `NM_015120.4`, of records whose reference is the two bases 2:73385942-73385943. The model numbers those bases cDNA 185 and 189 and holds three bases between them, `c.75` to `c.77`. Both rows print `cDNA_position` 185-189 and `CDS_position` 74-78, a `GIVEN_REF` of two bases and a `USED_REF` of five.

| The record's alternate | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| Two bases: the first kept, the second changed | a substitution at `c.75`, five reference bases against two | a `delins` of one base over `c.75_78` |
| Three bases, the first kept, in a record with a second alternate | a `delins` of three bases, the whole alternate, at the one position `c.75` | a `delins` of two bases, the ones the clip leaves, over `c.75_78` |

The third row is on `NM_002111.8`, of a record with a second alternate whose reference is the three bases 4:3074935-3074937 and whose first alternate holds two: the clip takes one base at each end. The model numbers the three reference bases cDNA 255, 262 and 263 and holds six bases of its own at cDNA 256 to 261. Both engines print `cDNA_position` 255-263 and a `USED_REF` of nine bases. Ensembl VEP writes `c.111del`, the first of the six bases. vep-rs writes a `del` of the seven positions `c.111` to `c.117`, cDNA 256 to 262.

The counts are of two sets of constructed records: the 53 records of section 3, and 24 records with two alternates each. By the kind the clip leaves:

| Kind | Rows | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| Substitution | 16 of the first set and 24 of the second | on the 16, a substitution at one position that names the model's whole reference allele; its values on the 24 are not shown here | on the 40, a `delins` of one base from the first base the genome lacks |
| `delins` | 75 of the second set | a `delins` of the record's whole alternate, the kept base included, over as many positions as the clip left | a `delins` from the first base the genome lacks to the last changed base, with the alternate bases the clip leaves |
| `del` | 1 of the second set | a `del` of one position | a `del` of seven positions |

### 4.2 Why vep-rs's row is the correct one

On each of the 75 `delins` rows the record's reference bases are the genome's, and the model holds bases of its own between the kept base and the first changed base. The model spells the kept base as the record does, and the row's `USED_REF` and `cDNA_position` cover all of these bases. Put the record's bases in place of the model's over that span and clip the kept base. What is left is a `delins` from the first base the genome lacks to the last changed base, with the alternate bases the clip leaves.

vep-rs's notation gives that sequence on 75 of 75 rows. Exchanging the changed bases alone gives it on none. Ensembl VEP's notation, applied to the model, gives neither sequence on 71 rows and has no exonic position on 4. On the 16 substitution rows Ensembl VEP's notation gives neither sequence either. On the `del` row the seven positions are the smallest change on the model at its 3'-most place, and `c.111del` gives neither sequence.

On 60 of the 75 rows, all at 1:16567177, the alternate's base at one end of vep-rs's span is the model's own there. The smallest change on the model is then one base shorter than that span: a `delins` on 45 rows and a `del` on 15. vep-rs names the longer span on all 75, and which of the two the notation should name on those 60 is not settled. On 12 of the 40 substitution rows the base the record writes is likewise the model's own at one end of the span, where the smallest change is a deletion, and the same question is open. The 15 models at 1:16567177 hold such bases at 1:16567174 too, three bases apart.

## 5. An insertion that the 3' placement carries to an exon's last base

**In short.** An allele that keeps a base stands at a place where the model holds bases the genome lacks, and the bases it adds are a copy that the 3' rule can carry to the exon's last base. Ensembl VEP 116.3 does not move such an allele. It writes an `ins` of the whole allele, the kept base included, between two consecutive positions at the record's point. vep-rs places the change 3'-most on the model's sequence, and carries an insertion past the exon's last base: onto the bases the model holds directly after that base where it holds any, or else onto the genome's intron. Nothing is set aside.

### 5.1 What each engine writes

One row of a set of 348 constructed records:

| Of the row | Value |
| --- | --- |
| The record | one reference base at 2:73385942; a first alternate of 255 bases, which are the kept base, the model's three own bases, the 250 exon bases after them and the intron's first base; a second alternate of one base |
| The model | `NM_015120.4`, forward strand, `BAM_EDIT=OK`; cDNA 1-185 on 2:73385758-73385942 and cDNA 189-438 on 2:73385943-73386192; cDNA 438 is the exon's last base |
| `cDNA_position` and `CDS_position`, both engines | 185 and 74 |
| Ensembl VEP 116.3 | an `ins` of 255 bases between `c.74` and `c.75`, and no `HGVS_OFFSET` |
| vep-rs | a `dup` of 251 bases, `c.78` to `c.327+1`: cDNA 189 to 438 and the intron's first base |

On 240 rows of the set, of an allele that keeps a base on a model that holds such bases at the record's point, Ensembl VEP writes an `ins` on 240 of 240. Its two positions are consecutive and stand 5 to 618 positions before the exon's last base. Its bases are the row's whole allele, and no row prints `HGVS_OFFSET`. Every other compared column of the 240 rows is equal in the two engines.

On 243 rows of an allele that keeps a base, vep-rs writes a notation that reaches past the exon's last base:

| Rows | vep-rs |
| --- | --- |
| 42, at five places | a duplication that ends 1 to 5 bases into the intron |
| 180, at two places where the model holds bases of its own directly after the exon's last base; this page does not show how these lie to the bases at the record's point | on 90 a duplication carried one base past the exon's last base, and on 90 an insertion between the first and the second base the model numbers after it |
| 21, at the five places of the 42 | an insertion between the exon's last base and the intron's first position |

The same changes written without the kept base take the same notation in vep-rs, on 14, 60 and 7 rows. So 28 rows carry an insertion between the exon's last base and the intron's first position: the 21, and the 7 of the same seven changes, each of which is written three ways that keep a base. Under `--shift_hgvs 0` vep-rs does not place the change, which stays where the clip leaves it.

Ensembl VEP's value is shown for the 240 rows: the 42, the 180, 3 of the 21, and 15 further rows on which vep-rs writes a duplication by the 3' rule. This page does not show it for the other 18 of the 21, or for the rows written without the kept base.

### 5.2 Why vep-rs's row is the correct one

On the row of 5.1 the record's point lies between 2:73385942 and 2:73385943. The model numbers those two bases cDNA 185 and 189 and holds three bases of its own between them. Three things follow for Ensembl VEP's `ins` between `c.74` and `c.75`, each from those numbers alone:

- `c.75`, cDNA 186, is not the record's second neighbour. It is the first of the bases the genome lacks, and two consecutive numbers cannot be the two neighbours of a point that the model numbers four apart.
- 255 bases are one more than the record adds. The kept base stands in the reference and in the alternate and is no inserted base.
- After the model's own three bases the record adds the exon's 250 bases and the intron's first base. A copy of those 251 bases gives the same sequence, and the 3' rule names the copy and not the written point.

vep-rs's `dup` of `c.78` to `c.327+1` is that copy. On each of the 240 rows the model numbers the record's two neighbours more than one apart, Ensembl VEP's two positions are consecutive, and its bases are the whole allele.

### 5.3 How the comparator treats these rows

The format comparator was run on the 240 rows of each engine. Applying the class `refseq_edited_model_seam_insertion_hgvsc`, it reports `HGVSc` as a difference on 240 of 240 rows and applies the class to none. Applying seven edited-model classes together, it applies none to any row: that class, `refseq_edited_model_shifted_insertion_bases`, `refseq_edited_model_hgvsc_offset`, `refseq_edited_model_hgvs_offset`, `refseq_edited_model_tail_term`, `refseq_edited_model_deleted_span_term` and `refseq_edited_model_deleted_span_codons`. Those runs did not apply this entry's class. By its test the class holds none of the rows of sections 3 to 5, whose two values differ in kind or in coordinates.
