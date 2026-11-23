# The HGVSc of an insertion whose placements run past the 3' shift's bound

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 shifts an insertion 3' at most 1001 positions less its count of bases, or 1,000 for an allele over 1,000 bases (over 1,001 on a reverse-strand transcript), so one repeating the reference further is described where that walk stops. vep-rs describes it at its most 3' placement, as a duplication or an insertion, reads the HGVSp at that placement too, and writes no HGVSc where it lies past the transcript's last base.

**Why this is a defect.** An insertion that repeats the bases beside it gives one changed sequence at every point of a run of placements, and a description names one of them. The HGVS nomenclature names the last: "3'rule: for all descriptions, the most 3' position possible of the reference sequence is arbitrarily assigned to have been changed" (`hgvs-nomenclature.org/stable/recommendations/general/`). Ensembl VEP follows that rule wherever the run ends at or before its bound: it walks to that end and describes the insertion there.

The bound is a count of the walk and no property of the sequence. Where the run is longer, the walk ends on a base that matches. The notation written there is an `ins` of the inserted bases turned by the steps walked, or a `dup` where the bases after or before the stop repeat them, beside an `HGVS_OFFSET` that is the bound. Either way the placement is short of the run's end, and where the change is a duplication at that end, an `ins` does not name it one. The nomenclature's Duplication page asks otherwise: "when a variant can be described as a duplication, it must be described as a duplication and not as, e.g., an insertion" (`recommendations/DNA/duplication/`).

The protein notation is the predicted consequence of the change the HGVSc describes, so it belongs to the same placement. Ensembl VEP reads it at the span where the walk stopped. Where that span lies in the coding sequence and the run ends outside it, its row carries a notation that the description at the last placement does not have. Where the record writes the insertion in an intron and the run ends between two coding bases, its row has none.

**What vep-rs does instead.** Under `--hgvs`, vep-rs walks to the end of the run and describes the insertion there: a duplication where the bases before that point repeat the inserted bases, an insertion otherwise. `HGVS_OFFSET` is the count of positions from the written point to it, and the HGVSIntronOffset plugin's two fields are the intron parts of that notation's coordinates. Where the last placement lies past the transcript's last base, the change has no description on that transcript, and vep-rs writes no HGVSc and neither plugin field. Under `--shift_3prime 1` the whole row is annotated at the span the bounded walk reaches, and vep-rs describes the insertion at that span.

vep-rs reads the HGVSp at the placement its HGVSc names. Between two coding bases it is the protein notation of the insertion there, also where the record writes the insertion in an intron inside the coding region and the codons the transcript holds give the notation, the reading of [`hgvsp_shifted_indel_written_span_gate`](hgvsp-shifted-indel-written-span-gate.md). A record written in a UTR keeps no HGVSp there, as in Ensembl VEP (line 1731). In an intron, in a UTR or past the transcript's last base the row has no HGVSp, as an insertion described there has none in either engine. Without `--shift_3prime 1`, every other column of the row, the Consequence among them, is annotated where the record writes the insertion, in both engines.

**Where it is in Ensembl VEP's source** (`ensembl-variation`, `ensembl-vep` and `VEP_plugins` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 427-437, `_genomic_shift`: the shift takes 1,000 bases of the chromosome on either side of the record. Line 427 sets the figure, lines 431-432 take the slice and lines 436-437 cut the two flanks.
- The same file, lines 303-306, `perform_shift`: the loop's limit is the flank's length less the allele's length, so the loop compares at most 1001 less that length times. For an allele longer than the flank, line 304 sets the negative limit to the flank's length, 1,000, and comparison 1,001 finds no base. Going 5', on a reverse-strand transcript, line 303 adds 1 to the limit and the loop counts from 1: an allele of 1,001 bases has the limit 0, which line 304 leaves, and is not walked.
- The same file, lines 323-348: a base that matches moves the insertion one position and turns its bases (lines 326 and 338-343), and the first base that differs ends the loop (line 347). The loop also ends at its limit, with no base differing. Line 388 stores the count of steps.
- The same file, `hgvs_transcript`: line 1441 runs the shift, line 1444 takes the turned bases and line 1446 the count of steps. Lines 1474-1483 ask for the notation at the insertion's coordinates plus that count (lines 1477-1478). Line 1472 ends the routine without a notation where the shifted end lies past the transcript's last base.
- The same file, `hgvs_protein`: line 1711 runs the same shift and lines 1715-1721 take its count of steps. Lines 1730-1738 return no notation unless the allele is coding where the record writes it (line 1731) and the span moved by that count has a translation start and end (lines 1732-1733).
- `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, lines 570-597: an insertion is a `dup` where as many bases after its point, or before it, equal the inserted bases (lines 577-587), and an `ins` otherwise (lines 592-594).
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1761-1764: the count of steps is printed as `HGVS_OFFSET`, times the strand, beside the notation.
- `HGVSIntronOffset.pm`, lines 62-72: the plugin writes the intron parts of the notation's coordinates, which `TranscriptVariationAllele.pm` stores at lines 1543 and 1545.

**What correcting it would change.** Ensembl VEP would go on along the reference while the bases match and describe the insertion at the end of its run. For the record of section 1 that is a duplication beside `HGVS_OFFSET=603`, where it writes an insertion beside `HGVS_OFFSET=399`. Lines 1732-1733 would find no translation start at that placement, which lies in an intron, and the row would carry no HGVSp.

**The class it produces.** `hgvs_insertion_left_where_the_shift_stops`: a differing HGVSc, HGVSp, `HGVS_OFFSET` and pair of HGVSIntronOffset fields on an insertion's row. Ensembl VEP's row carries an `ins` or a `dup` of the row's transcript beside an `HGVS_OFFSET` that is, in size, exactly the walk's bound for the allele's length: 1001 less the length, or 1,000 for an allele of more than 1,000 bases. The walk ended on its bound and not on a base that differs. vep-rs's row has one of two shapes.

Further 3': a `dup` or an `ins` in the same numbering, beside an `HGVS_OFFSET` of the same sign and a larger size, negative where Ensembl VEP's row carries none. A `dup` whose two coordinates are exonic names as many bases as the allele holds, or fewer around an intron where the comparator is given the cache. An `ins` writes the allele's bases as the transcript's strand holds them, turned by the steps its `HGVS_OFFSET` states. The HGVSc, the `HGVS_OFFSET` and the plugin's two fields are the class.

Past the last base: no HGVSc, no `HGVS_OFFSET` and neither plugin field. The HGVSc, the `HGVS_OFFSET` and the two fields are the class.

On a row of either kind the HGVSp is the class's in two shapes. In the first Ensembl VEP writes a notation and vep-rs none, and vep-rs's placement lies outside the coding sequence: it writes no HGVSc, or its `dup` ends, or its `ins` stands, at a coordinate with an intron part or in a UTR. An Ensembl VEP notation that names the stop and nothing else, on a coding sequence the row's flags show incomplete, is left to its own entry.

In the second shape Ensembl VEP writes none and vep-rs a frameshift, an insertion, a duplication or a deletion-insertion named at a residue other than the stop. The row then carries no CDS position, and vep-rs's `dup` ends, or its `ins` stands, on coding coordinates. Two notations that differ stay a difference.

One allele has the bound 0: 1,001 bases on a reverse-strand transcript. Ensembl VEP's row there carries no `HGVS_OFFSET`, and the class takes it the same way where the row is an insertion's own. Its Location names the two bases around the insertion's point, and the record's alleles, where the row spells them, hold no reference base. An allele that keeps bases of its span carries no `HGVS_OFFSET` either, and its row stays a difference.

A `dup` names a run of as many bases of the chromosome as the allele holds, after its point or before it, and each end of the run is mapped to a coordinate alone. Two exonic coordinates therefore count the allele's bases less those of every intron between them.

Given vep-rs's cache through `--json-cache`, the format comparator places the two runs from the row's Location, the `HGVS_OFFSET` beside the notation and the transcript's exon spans. On the row without an `HGVS_OFFSET`, Ensembl VEP's `dup` is the class's only where its two coordinates are a run's first and last base. A row whose transcript the cache does not hold on the row's strand stays a difference there. So does a `dup` with one coordinate in the 3' UTR and the other before it, whose distance the notation does not carry. A vep-rs `dup` of fewer bases than the allele holds is the class's where its coordinates are a run's.

Without the cache the comparator has no intron to count by. On the row without an `HGVS_OFFSET` it takes an Ensembl VEP `dup` of two exonic coordinates from 2 bases up to as many as the allele holds, whatever intron would explain the count. It leaves one count out: two fewer than the allele holds, which an allele that keeps the two bases of its span adds. A vep-rs `dup` of two exonic coordinates names as many bases as the allele holds, and one around an intron stays a difference.

Where the rows carry the plugin's fields, each row's two are the intron parts of its own notation's coordinates, in their printed order. A deletion's row stays a difference: both engines describe a deletion in a long run where the walk stops. So does an Ensembl VEP `HGVS_OFFSET` of any other size, and a vep-rs `HGVS_OFFSET` that is no larger than the bound or has the other sign. So do an `ins` whose bases are not the allele's turned by the steps walked, and plugin fields that are not their notation's own.

The class is built into the format comparator and applied once per row, and the golden corpus builder and the test harness name it. The concordance figures score consequence sets and never see these fields.

**How often the shape occurs.** Tier 2, `--hgvs` being needed. The shape is an insertion whose run of placements is longer than the walk. Two of the 15,646,481 records of [the index's five chromosome 21 sets](../intended-divergences.md#differences-outside-the-classes) carry it, 0.13 per million: ClinVar variation 4718230, an insertion of 713 bases, at 21:46924434 on GRCh37 and at 21:45504520 on GRCh38. The three population sets hold 993,227 records that insert bases, and on none of them is the run longer than the walk. No record of a committed corpus that holds a reference file carries the shape.

On the genome-wide ClinVar GRCh37 file Ensembl VEP 116.3 writes an `HGVS_OFFSET` equal to the bound on 236 rows of an insertion, from 61 records and 171 transcripts. vep-rs writes a duplication further 3' on 204 of them, an insertion further 3' on 10 and no HGVSc on 21. On the one row left the run ends exactly where the walk does, and the engines agree. On the GRCh38 file the walk reaches its bound on 1,289 rows, from 72 records and 864 transcripts: 1,168, 43, 77 and 1 in the same order.

On some of those rows the HGVSp differs as well:

| HGVSp of the row | GRCh37 | GRCh38 |
|---|---|---|
| Ensembl VEP a notation, vep-rs none: the last placement in an intron or in the 3' UTR | 6 | 31 |
| The same past the last base: the row of section 2, which the class of another entry takes | 1 | 0 |
| Ensembl VEP none, vep-rs a frameshift or an insertion: a record written in an intron | 14 | 75 |
| Two notations that differ | 0 | 0 |

The inserted alleles are hundreds to thousands of bases long, and every row is a forward-strand transcript's. Each of these records writes its insertion at the first point of its run, which on a reverse-strand transcript is the most 3' point.

## 1. A copy of 602 bases described 204 positions before the end of its run

**In short.** ClinVar variation 1050537, on GRCh37, inserts 602 bases between 13:32906922 and 13:32906923, a copy of the 602 bases after that point. The copy can stand at any of 604 points, the last between 32907525 and 32907526, 603 positions on. Ensembl VEP 116.3 walks 399 positions and writes an insertion between c.1706 and c.1707 of ENST00000380152 beside `HGVS_OFFSET=399`. vep-rs writes `c.1309_1909+1dup` beside `HGVS_OFFSET=603`, and an end offset of 1 where Ensembl VEP writes 0. Ensembl VEP reads a frameshift where its walk stops; vep-rs writes no HGVSp, the last placement lying in intron 10.

### 1.1 What a reader needs to know first

An insertion that repeats the bases after it can be written at more than one point. Moved one base 3', with its bases turned by one, it gives the same changed sequence. Both engines walk an insertion 3' along the chromosome, on the transcript's strand, while each base passed equals the next of the inserted bases, taken in a cycle. `HGVS_OFFSET` is the count of positions walked.

Ensembl VEP takes 1,000 bases of the reference on either side of the record for the walk, and at most 1001 steps less the count of inserted bases: 399 steps for 602 bases.

ENST00000380152 (ENSG00000139618) is a coding transcript of 27 exons on the forward strand, 32889611 to 32973347. Its exon 10 is 32906409 to 32907524, and intron 10 starts at 32907525. The record's point lies in exon 10 between c.1307 and c.1308, so the exon's last base, 32907524, is c.1909.

### 1.2 The record

```
13	32906922	1050537	A	<603 bases>	.	.	.
```

The alternate allele is the reference base followed by the 602 inserted bases, given here by its count. The whole line is an input line of the entry's reproduction, `scripts/divergences/reproductions/hgvs_insertion_left_where_the_shift_stops.toml`. Both engines run `--hgvs --plugin HGVSIntronOffset` on GRCh37 with a reference FASTA.

### 1.3 What each engine writes

Each row is the engine's own, with every run of 12 or more letters in the Allele, Amino_acids and Codons columns, and the inserted bases of the HGVSc, given by its count in angle brackets. Ensembl VEP 116.3, default format:

```
1050537	13:32906922-32906923	<602 bases>	ENSG00000139618	ENST00000380152	Transcript	frameshift_variant	1540-1541	1307-1308	436	K/<202 residues>	aag/<605 bases>	-	IMPACT=HIGH;STRAND=1;HGVSc=ENST00000380152.3:c.1706_1707ins<602 bases>;HGVSp=ENSP00000369497.3:p.Leu638LysfsTer23;HGVS_OFFSET=399;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

vep-rs:

```
1050537	13:32906922-32906923	<602 bases>	ENSG00000139618	ENST00000380152	Transcript	frameshift_variant	1540-1541	1307-1308	436	K/<202 residues>	aag/<605 bases>	-	IMPACT=HIGH;STRAND=1;HGVSc=ENST00000380152.3:c.1309_1909+1dup;HGVS_OFFSET=603;HGVS_IntronEndOffset=1;HGVS_IntronStartOffset=0
```

The HGVSc, the HGVSp, `HGVS_OFFSET` and `HGVS_IntronEndOffset` differ, and every other column and field is equal. The 602 bases of Ensembl VEP's `ins` are the record's, turned by 399. ENST00000544455, a second transcript of the gene with the same exon 10, carries the same notations in each engine.

### 1.4 Why Ensembl VEP writes that

Lines 427-437 take the 1,000 bases after 32906922, and the loop of lines 303-306 runs 1001 - 602 = 399 times. Each of the 399 bases it compares matches, so the walk ends on its limit, between 32907321 and 32907322, with 204 matching bases ahead of it. Line 388 stores 399.

`hgvs_transcript` adds 399 to the insertion's coordinates (lines 1477-1478), which gives the point between c.1706 and c.1707. The duplication test of `Utils/Sequence.pm` fails there. The 602 bases after the point are not the turned bases, since the run ends 204 positions on. The 602 bases before it reach 203 bases back past the record's own point, where the reference does not repeat the inserted bases. The notation is an `ins`, and both of its coordinates being exonic, the plugin writes 0 and 0.

`hgvs_protein` runs the same walk (line 1711) and reads the span 399 positions on. Between c.1706 and c.1707 that span has a translation start and end, so lines 1730-1738 let the notation stand. Read there, the inserted bases shift the frame: the first residue that differs is 638, and a stop stands 22 residues after it, `p.Leu638LysfsTer23`.

### 1.5 Why vep-rs's row is the correct one

The 602 bases after the record's point, 32906923 to 32907524, are the inserted bases. The base at 32907525 equals the first of them, and the base at 32907526 differs from the second. The run of placements therefore ends between 32907525 and 32907526, 603 positions from the record's point. By the nomenclature's rule that is the place the description takes: "the most 3' position possible of the reference sequence is arbitrarily assigned to have been changed" (`hgvs-nomenclature.org/stable/recommendations/general/`).

The 602 bases before that point, 32906924 to 32907525, are the inserted bases turned by one: the inserted bases are a copy standing directly 3' of its original. The nomenclature defines a duplication as a change where "a copy of one or more nucleotides is inserted directly 3' of the original copy of that sequence" (`recommendations/DNA/duplication/`), and asks that it be described as one. 32906924 is c.1309, and 32907525 is the first base of intron 10, c.1909+1: `c.1309_1909+1dup`. The intron parts of those two coordinates are 0 and 1.

The last placement lies between the first and the second base of intron 10. An insertion described there changes no coding base, and its row has no protein notation in either engine: lines 1732-1733 find no translation start for a span in an intron. vep-rs writes no HGVSp. Ensembl VEP's frameshift is read 204 positions before the run's end, which its walk reaches on any run shorter than the bound.

### 1.6 Try it yourself

Save the `text` lines of the reproduction named in 1.2 as `record.vcf` and annotate it against a GRCh37 cache with the reference FASTA. vep-rs's command is below; Ensembl VEP takes the same record with its own cache and the FASTA under `--offline --cache`:

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh37 --hgvs --plugin HGVSIntronOffset
```

vep-rs's row for ENST00000380152 carries `HGVSc=ENST00000380152.3:c.1309_1909+1dup`, `HGVS_OFFSET=603` and no HGVSp. Ensembl VEP 116.3's carries an `ins` between c.1706 and c.1707, `HGVS_OFFSET=399` and `HGVSp=ENSP00000369497.3:p.Leu638LysfsTer23`.

### 1.7 How the comparator treats this class

The format comparator's rule (`tools/vep-formats-compare`, `insertion_left_where_the_shift_stops`) accepts the HGVSc, the HGVSp, the `HGVS_OFFSET` and `HGVS_IntronEndOffset` of this row. Ensembl VEP's `HGVS_OFFSET` is 399, the bound for 602 bases, beside an `ins`. vep-rs's is 603, larger and of the same sign, beside a `dup`. Each row's two plugin fields are the intron parts of its own notation's coordinates. `HGVS_IntronStartOffset` is 0 on both rows and does not differ. vep-rs's `dup` ends at a coordinate with an intron part, so Ensembl VEP's HGVSp against none is the class's.

## 2. The same insertion on a transcript that ends inside the run

**In short.** ENST00000530893, a transcript of the same gene, ends at 32907428. The walk of section 1 stops 107 bases before that last base, and Ensembl VEP 116.3 writes an insertion between c.1337 and c.1338 beside `HGVS_OFFSET=399`, with the plugin's fields 0 and 0. The run's last placement lies 97 bases past the last base. vep-rs writes no HGVSc, no HGVSp, no `HGVS_OFFSET` and neither plugin field.

### 2.1 What a reader needs to know first

ENST00000530893 (ENSG00000139618) is a forward-strand transcript of 10 exons and 2,011 spliced bases, 32889642 to 32907428. It is flagged `cds_end_NF`: its coding sequence is incomplete at the 3' end. Its exon 10 starts at 32906409, as in ENST00000380152, and ends 96 bases earlier, at 32907428, the transcript's last base. The record's point lies between c.938 and c.939.

### 2.2 The record

The record is the one of section 1.2, with the same flags.

### 2.3 What each engine writes

The rows follow the rule of section 1.3. Ensembl VEP 116.3, default format:

```
1050537	13:32906922-32906923	<602 bases>	ENSG00000139618	ENST00000530893	Transcript	frameshift_variant	1505-1506	938-939	313	K/<202 residues>	aag/<605 bases>	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000530893.2:c.1337_1338ins<602 bases>;HGVSp=ENSP00000435699.2:p.Ter482IleextTer55;HGVS_OFFSET=399;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

vep-rs:

```
1050537	13:32906922-32906923	<602 bases>	ENSG00000139618	ENST00000530893	Transcript	frameshift_variant	1505-1506	938-939	313	K/<202 residues>	aag/<605 bases>	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF
```

The HGVSc, the HGVSp, the `HGVS_OFFSET` and the plugin's two fields differ. Ensembl VEP reads its HGVSp where the walk stops. vep-rs writes none: the last placement lies past the transcript's last base. It writes none on this transcript wherever the insertion is described, its coding sequence being incomplete, the entry [`hgvsp_stop_named_on_incomplete_coding_sequence`](hgvsp-stop-named-on-incomplete-coding-sequence.md).

### 2.4 Why Ensembl VEP writes that

The walk is the one of section 1.4: 399 steps, to the point between 32907321 and 32907322, which is between c.1337 and c.1338 of this transcript. That point lies 107 bases before the transcript's last base, so line 1472 lets the notation stand, and the duplication test fails as it does in section 1.4.

The bound alone keeps the notation inside the transcript. A walk to the end of the run adds 603 to the insertion's coordinates and carries its end 97 bases past the transcript's last base, where line 1472 ends the routine without a notation.

### 2.5 Why vep-rs's row is the correct one

The run's last placement is between 32907525 and 32907526, and the 602 bases it copies are 32906924 to 32907525. The transcript's last base is 32907428, so the placement and the last 97 of the copied bases lie outside the transcript. The nomenclature does not let a transcript's description reach them: "it is not allowed to describe variants in nucleotides beyond the boundaries of a transcript reference sequence, using that transcript reference sequence" (`hgvs-nomenclature.org/stable/background/numbering/`).

The change has no description on this transcript, so vep-rs writes no HGVSc. Without a notation there is no count of positions it was moved by and no intron part of its coordinates.

### 2.6 Try it yourself

The command of section 1.6 writes this row too. vep-rs's row for ENST00000530893 ends `FLAGS=cds_end_NF`. Ensembl VEP 116.3's goes on with an HGVSc between c.1337 and c.1338, an HGVSp, `HGVS_OFFSET=399` and the two plugin fields.

### 2.7 How the comparator treats this row

The rule accepts the HGVSc, the `HGVS_OFFSET` and both plugin fields of the row. Ensembl VEP's `HGVS_OFFSET` is the bound beside an `ins` whose two plugin fields are its own, and vep-rs writes none of the four. The HGVSp is left to its own class: under this class alone the row stays a difference in that field.

## 3. An insertion written in an intron whose run ends in the coding sequence

**In short.** ClinVar variation 1735214, on GRCh37, inserts 1,000 bases between 3:37047542 and 3:37047543, in intron 4 of ENST00000231790, a copy of the 1,000 bases after that point. The copy's last placement is between 37048542 and 37048543, 1,000 positions on, between c.441 and c.442 in exon 5. Ensembl VEP 116.3 walks one position and writes an insertion in the intron beside `HGVS_OFFSET=1` and no HGVSp. vep-rs writes `c.381-939_441dup` beside `HGVS_OFFSET=1000` and `p.Thr148AspfsTer54`.

### 3.1 What a reader needs to know first

ENST00000231790 (ENSG00000076242) is a coding transcript of 19 exons on the forward strand, 37034823 to 37092409, whose protein has 756 residues. Intron 4 is 37045966 to 37048481, 2,516 bases, and exon 5 is 37048482 to 37048554, c.381 to c.453. The record's point lies after the 1,577th base of intron 4. For 1,000 inserted bases the walk's bound is one step, 1001 less 1,000.

### 3.2 The record

```
3	37047542	1735214	A	<1001 bases>	.	.	.
```

The alternate allele is given by its count, as in section 1.2, and the whole line is an input line of the entry's reproduction.

### 3.3 What each engine writes

The rows follow the rule of section 1.3. Ensembl VEP 116.3, default format:

```
1735214	3:37047542-37047543	<1000 bases>	ENSG00000076242	ENST00000231790	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000231790.2:c.381-939_381-938ins<1000 bases>;HGVS_OFFSET=1;HGVS_IntronEndOffset=-938;HGVS_IntronStartOffset=-939
```

vep-rs:

```
1735214	3:37047542-37047543	<1000 bases>	ENSG00000076242	ENST00000231790	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000231790.2:c.381-939_441dup;HGVSp=ENSP00000231790.2:p.Thr148AspfsTer54;HGVS_OFFSET=1000;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=-939
```

The HGVSc, the HGVSp, `HGVS_OFFSET` and `HGVS_IntronEndOffset` differ, and every other column and field is equal.

### 3.4 Why Ensembl VEP writes that

The loop of lines 303-306 runs 1001 - 1000 = 1 time. The one base it compares matches, so the walk ends on its limit, between 37047543 and 37047544, with 999 matching bases ahead of it. That point lies in the intron, between the 939th and the 938th base before exon 5, and the duplication test of `Utils/Sequence.pm` fails there.

`hgvs_protein` returns no notation at line 1731: the allele lies in an intron where the record writes it, so it is not coding. The span one position on has no translation start either.

### 3.5 Why vep-rs's row is the correct one

The 1,000 bases after the record's point, 37047543 to 37048542, are the inserted bases, and the base at 37048543 differs from the first of them. The run of placements ends between 37048542 and 37048543, and the 1,000 bases before that point are the inserted bases: a copy standing directly 3' of its original. 37047543 is c.381-939 and 37048542 is c.441: `c.381-939_441dup`, with the intron parts -939 and 0.

c.441 is the last base of the codon of residue 147, so the copy stands between two codons. Its 1,000 bases are no whole number of codons, and the frame shifts after residue 147. Read from the reference and the transcript's coding sequence, the first residue that differs is 148, a threonine that becomes an aspartate, and a stop stands 53 residues after it: `p.Thr148AspfsTer54`. vep-rs writes the protein notation of an insertion described between two coding bases, as it does for a run shorter than the bound.

### 3.6 Try it yourself

The command of section 1.6 annotates this record too. vep-rs's row for ENST00000231790 carries `HGVSc=ENST00000231790.2:c.381-939_441dup`, `HGVS_OFFSET=1000` and `HGVSp=ENSP00000231790.2:p.Thr148AspfsTer54`. Ensembl VEP 116.3's carries an `ins` in intron 4 beside `HGVS_OFFSET=1` and no HGVSp.

### 3.7 How the comparator treats this row

The rule accepts the HGVSc, the HGVSp, the `HGVS_OFFSET` and `HGVS_IntronEndOffset` of the row. Ensembl VEP's `HGVS_OFFSET` is 1, the bound for 1,000 bases, and the row is the second shape of the class's HGVSp.
