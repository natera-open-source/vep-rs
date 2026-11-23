# The stop count of a frameshift named at the first residue

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 gives no count to a frameshift whose first changed residue is residue 1: it writes the unknown count, `fsTer?`, whether or not the new reading frame meets a stop codon. On a transcript at whose coding start no start codon stands, the count is owed, and vep-rs writes it. Such a transcript carries the `cds_start_NF` flag (coding sequence start not found), or its cache record shows no start codon.

**Why this is a defect.** The HGVS nomenclature's Frameshift page counts the position of the stop in the new reading frame from the first changed amino acid, "ending at the first stop codon encountered". It keeps the unknown count for a frame that "does not encounter a new translation termination (stop) codon". On the rows of this entry the new frame meets a stop codon, so the unknown count says what the sequence does not.

Ensembl VEP withholds the count because of a number, not because of anything in the sequence. Its count routine is handed the number of the residue before the first changed one, and it returns no count unless that number is true. Before residue 1 the number is 0, which Perl reads as false. Where the first changed residue has any later number, the same routine counts to the same stop codon. Section 1 shows one deletion on five transcripts, counted on four of them.

Residue 1 of a flagged transcript is no initiator. The flag says that the coding sequence goes on before the transcript's first base, so the first residue lies inside a protein and is numbered 1 because the transcript's record begins there. Ensembl VEP treats the row the same way: it writes no `start_lost` on a flagged transcript and describes the change as a frameshift. The frame after the change is the transcript's own, and the count is the place of its first stop codon.

**What vep-rs does instead.** vep-rs asks the transcript's cache record whether a start codon stands at the coding start, by the rule of [`start_codon_overlap_without_start_codon`](start-codon-overlap-without-start-codon.md). By that rule a start codon stands only on a record without the flag whose first codon is whole. It stands there where the first codon is a start codon of the genetic code, the recorded protein opens with a methionine, or the record carries a consensus coding sequence (CCDS) id.

Where no start codon stands, vep-rs counts a frameshift named at residue 1 as it counts any other: from the first changed residue, counted as 1, to the first stop codon of the new frame. The stop is sought through the alternate coding sequence and the transcript's 3' UTR, the bases Ensembl VEP's count reads, and no farther. Where the new frame meets no stop codon there, vep-rs keeps the unknown count, as Ensembl VEP writes it.

Where a start codon stands, Ensembl VEP's test of that number withholds the count as well, and there no count is owed. vep-rs keeps the unknown count, as Ensembl VEP writes it. Both engines read the count from a translation of the cached coding sequence with the allele written into it, begun at coding position 1. On such a transcript that translation is not shown to be the variant's protein.

The edited transcript can still hold the start codon and every codon after it. On the last row of section 2 it does, behind a 5' UTR that is one base longer. A deletion that runs from the 5' UTR into the start codon can take the frameshift form too, and the HGVS nomenclature gives a changed start codon's consequence no length. This entry charges no defect to Ensembl VEP on such a row and changes nothing in it.

A transcript without the flag whose record shows no start codon is counted as a flagged one is. For a frameshift in its first codon Ensembl VEP reads a lost start, and that row belongs to another entry. Section 3 says which rows of such a transcript are this entry's. This entry changes nothing in a row but the count.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the module the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2477-2526, `_stop_loss_extra_AA`: line 2483, its first statement, returns no count unless the position it is handed is true. Past it, line 2492 translates the alternate coding sequence and 3' UTR, line 2503 finds the first stop, and line 2508 takes its distance from the position.
- Lines 1898-2044, `_get_hgvs_protein_format`: line 2003 asks for the count of a frameshift and hands over the first changed residue's number less one, which is 0 for residue 1. Line 2006 writes `?` where no count came back, and line 2013 writes `fsTer` and the count.
- Lines 2320-2365, `_get_fs_peptides`: line 2340 starts at the first residue of the row's protein position, and lines 2348-2362 step on while the two translations agree. The residue they stop at is the first changed one.
- Lines 2161-2164: where `start_lost` holds, the description is rewritten as the lost start's, and the frameshift form is not reached.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, line 975: the start codon test returns 0 for a transcript that carries `cds_start_NF`, so no row of such a transcript holds `start_lost`.

**What correcting it would change.** A guard that asked whether the position was given, and not whether it is a true number, would let the routine count from residue 1. The rows of this entry would carry the count in place of `?`, and no other column would change.

**The class it produces.** `hgvsp_frameshift_count_at_the_first_residue`: Ensembl VEP's description is a frameshift at residue 1 with the unknown count, and vep-rs's is the same frameshift, with the same two residues, with a count. The transcript's cache record carries `cds_start_NF` or shows no start codon, and neither row holds `start_lost`. Section 1.7 says how a comparison treats the pair.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data carries the shape, and the divergence needs `--hgvs` or a flag that switches it on, such as `--everything`. That page is the index of the documented differences, and it counts each shape over 15,646,481 chromosome 21 records. The shape the record alone shows is an allele whose length differs from the reference allele's by a number that is no multiple of three. Of the records, 1,692,699 hold one, 108,184 per million. The rate is an upper bound: the transcript decides whether the first changed residue is residue 1 and whether a start codon stands.

The Ensembl caches of release 116 flag 14,140 coding transcripts `cds_start_NF` on GRCh38 and 13,798 on GRCh37. The coding sequence opens on a whole codon on 5,349 of the 14,140 and on 5,624 of the 13,798. On the others the first codon is incomplete and reads as the unknown residue, so the first changed residue is numbered 2 or more.

With `--hgvs` on the genome-wide ClinVar dataset, vep-rs writes a count for a frameshift named at residue 1 on 34 rows of flagged transcripts on GRCh38 and on 26 on GRCh37. Of the 60 rows, 52 carry `frameshift_variant`. The other 8 are rows of a deletion that begins before the transcript's first base, as the second record of section 2 does. One further row on GRCh38 stands on a transcript without the flag that shows no start codon. Ensembl VEP reads a lost start there, and the row belongs to the entry section 3.1 names.

The committed corpus `tests/golden/116/GRCh37-hgvs` holds four records of the class on ENST00000415023: three insertions between coding bases 1 and 2 and a deletion of five bases.

## 1. Two bases deleted from the first codon of a flagged transcript, on ENST00000415023

**In short.** `21 35247663 . GGA G` on GRCh37, a made-up record, deletes the two bases 21:35247664-35247665. They open the coding sequence of ENST00000415023, a transcript of the gene ITSN1 that is flagged `cds_start_NF`. Ensembl VEP 116.3 writes `p.Glu1LysfsTer?` and vep-rs writes `p.Glu1LysfsTer4`, each after the protein's identifier. The new frame reads three residues and then a stop codon, at coding bases 12 to 14. On four other transcripts of the gene the same residue has another number, and both engines count to that stop codon and write `fsTer4`.

### 1.1 What a reader needs to know first

`--hgvs` adds two notations to a row, `HGVSc` for the transcript and `HGVSp` for the protein. For a frameshift the protein notation names the first residue that changes, its number, the residue that replaces it, and `fsTer` with a count. `p.Glu1394LysfsTer4` says that residue 1394, a glutamate, becomes a lysine, and that the new reading frame ends on a stop codon four residues on, the changed residue counted as 1. `fsTer?` says that the place of the stop is unknown.

Both engines find the first changed residue by translating the coding sequence with the change made in it and comparing the result with the transcript's protein, residue by residue from the row's protein position. The first residue that differs is the first changed one.

A transcript flagged `cds_start_NF` is one whose annotation found no start for the coding sequence. Its record begins inside a coding sequence, and its first residue is whatever codon the record begins on. The rows of such a transcript print the flag in `FLAGS`.

### 1.2 The record

A made-up deletion of two bases, against the Ensembl GRCh37 cache of release 116:

```
21	35247663	first_codon_two_bases_deleted	GGA	G	.	.	.
```

Both engines drop the shared first base, so the record deletes 21:35247664-35247665. ENST00000415023 is a forward-strand transcript of five exons. Its first exon is 21:35247664-35247825, and the transcript has no 5' UTR and no 3' UTR: the coding sequence is its whole sequence, 519 bases, 173 codons. Coding base 1 is 21:35247664, and the transcript's cache record carries both `cds_start_NF` and `cds_end_NF`.

| Part | Coding bases | Genomic span |
| --- | --- | --- |
| Codon 1, a glutamate codon | 1 to 3 | 21:35247664-35247666 |
| The deleted bases | 1 and 2 | 21:35247664-35247665 |
| Codon 2, an asparagine codon | 4 to 6 | 21:35247667-35247669 |
| The new frame's stop codon | 12 to 14 | 21:35247675-35247677 |

The deletion can be written at no other place: the base before the two differs from the second of them, and the base after them from the first. Four other transcripts of the gene hold the two bases in their coding sequence. On ENST00000381318, the gene's canonical transcript, they are coding bases 4180 and 4181, the first two bases of codon 1394. On ENST00000420666, flagged `cds_start_NF`, they are the first two bases of codon 2, behind an incomplete first codon.

### 1.3 What each engine writes

The row on ENST00000415023 in the default format, under `--hgvs`. Ensembl VEP 116.3:

```
first_codon_two_bases_deleted	21:35247664-35247665	-	ENSG00000205726	ENST00000415023	Transcript	frameshift_variant	1-2	1-2	1	E/X	GAa/a	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_start_NF,cds_end_NF;HGVSc=ENST00000415023.1:c.1_2del;HGVSp=ENSP00000409800.1:p.Glu1LysfsTer?
```

vep-rs:

```
first_codon_two_bases_deleted	21:35247664-35247665	-	ENSG00000205726	ENST00000415023	Transcript	frameshift_variant	1-2	1-2	1	E/X	GAa/a	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_start_NF,cds_end_NF;HGVSc=ENST00000415023.1:c.1_2del;HGVSp=ENSP00000409800.1:p.Glu1LysfsTer4
```

The rows differ in the count of the `HGVSp` alone. Both carry `frameshift_variant`, the coding and protein positions and the same two residues. The record writes 8 rows in each engine, and the other 7 are equal in every column. Five of the 8 carry `frameshift_variant`, one for each transcript whose coding sequence holds the two bases. The table gives three of them. On ENST00000381285 and ENST00000399367 both engines write `fsTer4` as well, at residues 1394 and 1389.

| Transcript | First changed residue | Ensembl VEP | vep-rs |
| --- | --- | --- | --- |
| ENST00000381318 | 1394 | `p.Glu1394LysfsTer4` | `p.Glu1394LysfsTer4` |
| ENST00000420666 | 2 | `p.Glu2LysfsTer4` | `p.Glu2LysfsTer4` |
| ENST00000415023 | 1 | `p.Glu1LysfsTer?` | `p.Glu1LysfsTer4` |

### 1.4 Why Ensembl VEP writes that

The transcript carries `cds_start_NF`, so the start codon test returns 0 (`VariationEffect.pm` line 975) and the row holds no `start_lost`. The allele shifts the frame, and the protein notation takes the frameshift form.

`_get_fs_peptides` starts at the row's protein position, 1 (line 2340). The first codon of the edited sequence names a lysine and the transcript's first residue is a glutamate, so the two translations differ at once and the first changed residue is residue 1.

The frameshift form then asks for the count and hands over 1 - 1 = 0 (line 2003). `_stop_loss_extra_AA` returns at its first statement, because 0 is not true (line 2483). No sequence is built and no codon is read. Line 2006 writes `?` in the count's place.

On ENST00000381318 the first changed residue is 1394 and the routine is handed 1393. It translates the edited sequence and finds the first stop at residue 1397 (line 2503). Line 2508 gives 1397 - 1393 = 4. On ENST00000420666 it is handed 1, finds the stop at residue 5 and gives 5 - 1 = 4.

### 1.5 Why vep-rs's row is the correct one

The record removes coding bases 1 and 2 of ENST00000415023. The edited coding sequence begins with coding base 3, and its codons are coding bases 3 to 5, 6 to 8, 9 to 11 and 12 to 14. They name a lysine, a histidine and a proline, and the fourth is a stop codon of the standard genetic code. The first residue changes from a glutamate to a lysine, and the stop is the fourth codon counted from it: `p.Glu1LysfsTer4`.

The stop codon stands at 21:35247675-35247677, inside the first exon. On ENST00000381318 the same three bases are coding bases 4191 to 4193, and Ensembl VEP's own count reaches them: residue 1397, four counted from 1394. The frame, the stop codon and the count are the same on the three transcripts of the table. Only the number of the changed residue differs, and with it what line 2483 does.

`fsTer?` would be right for a frame that meets no stop codon. This one meets a stop codon after three residues, inside the transcript's own bases. The first residue is no initiator either: the flag says that the coding sequence goes on before 21:35247664, as it does on ENST00000381318, where 1,393 codons precede the same glutamate.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline --hgvs \
  --json_cache tests/golden/116/GRCh37-hgvs/json_cache \
  --fasta tests/golden/116/GRCh37-hgvs/reference.fa.gz \
  --species homo_sapiens --assembly GRCh37
grep ENST00000415023 out.txt
```

The committed GRCh37 cache holds the transcript and the committed reference holds chromosome 21. vep-rs writes `p.Glu1LysfsTer4` on the row. Ensembl VEP 116.3 takes the record with the cache options of [Reproducing these examples](reproducing.md), `--hgvs` and a GRCh37 reference, and writes `p.Glu1LysfsTer?`.

The entry's reproduction, `scripts/divergences/reproductions/hgvsp_frameshift_count_at_the_first_residue.toml`, runs this record and the three of section 2 on both engines. It holds nine rows to both recorded values: each record's row on ENST00000381318, the rows of `first_codon_two_bases_deleted` and `second_codon_two_bases_deleted` on ENST00000415023 and ENST00000420666, and the row of `deletion_from_before_the_first_base` on ENST00000419241. It reads one cache and so leaves out the record of section 3.

### 1.7 How the comparison treats this class

The format comparator (`tools/vep-formats-compare`) compares one Ensembl VEP output with one vep-rs output, every field of every row. Where a comparison names the class `hgvsp_frameshift_count_at_the_first_residue`, it sets aside an `HGVSp` whose two values are the pair of this entry. They are one frameshift at residue 1 of one protein, Ensembl VEP's with `?` and vep-rs's with a count, on rows that hold no `start_lost`. The class reads no other term of the rows.

The two rows do not show what stands at the coding start of an unflagged transcript, so the class reads the cache vep-rs annotated with, which the comparator takes as `--json-cache`. The transcript's record there must carry `cds_start_NF` or show no start codon. The comparator derives both descriptions itself: from the record's coding sequence and 3' UTR, from the row's allele, coding position and `HGVS_OFFSET`, and from the length of the input record's reference allele. The pair is set aside only where the two values are the derived ones.

A default or tab row takes that length from the records of `--input`. A deletion whose coding position is open at its start is read where its `HGVS_OFFSET` moves it, as on the second record of section 2.

A comparison that is given no cache sets nothing aside under this class, and each such row stays a difference. So does a count on a transcript whose start codon stands, a count the sequence does not give, a frameshift named at a later residue, and a default or tab row of a comparison that is given no input.

## 2. The other records of the reproduction

**In short.** Three more made-up records on transcripts of the same gene show what the count turns on. The same frame named at residue 2 is counted by both engines. A flagged transcript's row whose deletion begins before the transcript differs in the count as section 1's does. On a start codon that stands, both engines write the unknown count. Sections 2 and 3 number their parts as section 1 does and leave out the parts they do not need.

### 2.2 The records

```
21	35247666	second_codon_two_bases_deleted	AAA	A	.	.	.
21	35169856	deletion_from_before_the_first_base	TCCAT	T	.	.	.
21	35091133	base_inserted_before_a_start_codon	C	CA	.	.	.
```

### 2.3 What each engine writes

The `HGVSp` of each record's row on one transcript, after the protein's identifier:

| Record | Transcript | Ensembl VEP | vep-rs |
| --- | --- | --- | --- |
| `second_codon_two_bases_deleted` | ENST00000415023 | `p.Asn2HisfsTer3` | `p.Asn2HisfsTer3` |
| `deletion_from_before_the_first_base` | ENST00000419241 | `p.His1AsnfsTer?` | `p.His1AsnfsTer30` |
| `base_inserted_before_a_start_codon` | ENST00000381318 | `p.Met1AsnfsTer?` | `p.Met1AsnfsTer?` |

**The same frame named at residue 2.** The first record deletes 21:35247667-35247668, coding bases 4 and 5 of ENST00000415023. The first codon is untouched and the second changes from an asparagine to a histidine codon. From coding base 6 on, the edited sequence is read in the frame of section 1. Its stop codon is the same three bases, coding bases 12 to 14, three residues counted from residue 2. Ensembl VEP's routine is handed 1 and counts.

**A deletion that begins before the transcript's first base.** ENST00000419241, another transcript of the gene flagged `cds_start_NF`, begins at 21:35169858. The second record deletes 21:35169857-35169860: the base before the transcript and coding bases 1 to 3. Both engines write `coding_sequence_variant,5_prime_UTR_variant,NMD_transcript_variant` for a deletion that runs off the transcript's first base. They print `?-3` as its coding position and write no frameshift term, since the span has no coding position at that end. The base before the transcript is the same base as coding base 4, so the same deletion can be written at coding bases 1 to 4. Both engines print `HGVS_OFFSET=1` and make the protein notation on the deletion as written there.

Written there, the deletion changes residue 1 from a histidine to an asparagine, and the new frame's first stop is its codon 30, coding bases 92 to 94 of 234. vep-rs writes `p.His1AsnfsTer30`. On ENST00000381318, where the same histidine is residue 710, both engines write `p.His710AsnfsTer30`. The `HGVSc` of the row on ENST00000419241 differs too, as [`hgvs_shifted_deletion_clamped_span`](hgvs-shifted-deletion-clamped-span.md) documents.

**A start codon that stands.** On ENST00000381318, 21:35091134 is coding base 1, the first base of the start codon, and 21:35091133 is the last base of the 5' UTR. The third record inserts a copy of that first base between the two. Both engines write `5_prime_UTR_variant` and `HGVS_OFFSET=1`, and describe the insertion after coding base 1 as a frameshift at residue 1 with the unknown count.

The edited transcript holds the start codon whole, behind a 5' UTR that is one base longer, and every codon after it. The translation both engines read begins at the old coding position 1 and meets a stop codon at its codon 26. A count of 26 would give that length to a protein the record leaves as it was. vep-rs writes no such count, and the two engines' rows are equal.

## 3. A transcript without the flag that shows no start codon, on ENSMUST00000165454

**In short.** `JH584299.1 716294 . TTC T` on the mouse assembly GRCm39, a made-up record, deletes the base before ENSMUST00000165454 and its first coding base. The transcript carries no flag and shows no start codon. Ensembl VEP 116.3 writes `p.Leu1IlefsTer?` and vep-rs writes `p.Leu1IlefsTer16`. Neither engine writes `start_lost` or `start_retained_variant` on the row.

### 3.1 Which rows of such a transcript belong here

Ensembl VEP's start codon test reads the first three coding positions of every unflagged transcript as the start codon. For a frameshift in such a first codon it writes `start_lost` and the lost start's description: a residue, the number 1 and `?`. vep-rs writes `frameshift_variant` and the frameshift with its count. That pair of rows belongs to [`start_codon_overlap_without_start_codon`](start-codon-overlap-without-start-codon.md). Its page says how a comparison reads the two notations beside the terms.

A row of such a transcript is this entry's where Ensembl VEP reads no lost start and both engines describe the change as a frameshift at residue 1. The two protein notations then differ in the count alone. One such row is a deletion that begins before the transcript's first base. The start codon test asks both cDNA coordinates of the span (`VariationEffect.pm` line 981), and that span has none at one end.

### 3.2 The record

A made-up deletion of two bases, against the Ensembl GRCm39 cache of release 116:

```
JH584299.1	716294	deletion_from_before_an_unflagged_transcript	TTC	T	.	.	.
```

Both engines drop the shared first base, so the record deletes JH584299.1:716295-716296. ENSMUST00000165454 is a forward-strand transcript of four exons that begins at 716296. Its coding sequence opens on the transcript's first base and has 1,500 bases: 499 codons and a stop codon. The transcript's cache record carries no `cds_start_NF` and no CCDS id, which would name a complete coding sequence. Its first codon names a leucine and is no start codon of the genetic code, and its protein opens with that leucine. No start codon is shown.

The deletion removes the base before the transcript and coding base 1. The base before the transcript is the same base as coding base 2, so the same deletion can be written at coding bases 1 and 2.

### 3.3 What each engine writes

The row on ENSMUST00000165454 in the default format, under `--hgvs`. Ensembl VEP 116.3:

```
deletion_from_before_an_unflagged_transcript	JH584299.1:716295-716296	-	ENSMUSG00000091585	ENSMUST00000165454	Transcript	coding_sequence_variant,5_prime_UTR_variant	?-1	?-1	?-1	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000165454.3:c.2del;HGVSp=ENSMUSP00000126681.3:p.Leu1IlefsTer?;HGVS_OFFSET=1
```

vep-rs:

```
deletion_from_before_an_unflagged_transcript	JH584299.1:716295-716296	-	ENSMUSG00000091585	ENSMUST00000165454	Transcript	coding_sequence_variant,5_prime_UTR_variant	?-1	?-1	?-1	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENSMUST00000165454.3:c.1_2del;HGVSp=ENSMUSP00000126681.3:p.Leu1IlefsTer16;HGVS_OFFSET=1
```

Both rows carry the same terms, neither `start_lost` nor `start_retained_variant` among them, and `HGVS_OFFSET=1`. They differ in the count of the `HGVSp`. They differ in the `HGVSc` too, which is another entry's difference: [`hgvs_shifted_deletion_clamped_span`](hgvs-shifted-deletion-clamped-span.md). It is the one row the record writes in either engine.

### 3.5 Why vep-rs's row is the correct one

Written one base on, the deletion removes coding bases 1 and 2. The edited coding sequence begins with coding base 3. Its first codon names an isoleucine where the protein opens with a leucine, and its sixteenth codon, coding bases 48 to 50, is a stop codon: `p.Leu1IlefsTer16`. The transcript's record shows no start codon at 716296, so residue 1 is not shown to be an initiator here either, and the frame after the change is the transcript's own.

### 3.6 Try it yourself

With the record of 3.2 saved as `record.vcf`, from the repository's root:

```
vep -i record.vcf -o out.txt --offline --hgvs \
  --json_cache tests/golden/116/GRCm39-mouse/json_cache \
  --fasta tests/golden/116/GRCm39-mouse/reference.fa.gz \
  --species mus_musculus --assembly GRCm39
grep ENSMUST00000165454 out.txt
```

The committed mouse cache holds the transcript and the committed reference holds the scaffold. vep-rs writes `p.Leu1IlefsTer16` on the row, and Ensembl VEP 116.3, with its GRCm39 cache, writes `p.Leu1IlefsTer?`.
