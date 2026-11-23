# Stop codon terms at a coding end where no stop codon stands

**The defect in one sentence.** Ensembl VEP 116.3 reads the last three coding positions of every transcript without `cds_end_NF`, the cache's flag for a coding sequence whose 3' end was not found, as its stop codon, whatever those positions hold. For an insertion or deletion there it writes `stop_lost` or `stop_retained_variant` where the bases the change meets are no base of a stop codon. vep-rs writes a stop codon's terms only where a stop codon stands, and gives these rows `frameshift_variant`, an in-frame term (`inframe_insertion`, `inframe_deletion`) or `coding_sequence_variant`.

**Why this is a defect.** `stop_lost` says that the stop codon is gone after the change, and `stop_retained_variant` that it is still there. The Sequence Ontology's definitions of both terms speak of the terminator codon (SO:0001578, SO:0001567). A coding sequence whose last codon names a residue has no terminator codon at its end, and neither has one that ends inside a codon that the transcript's cache record does not complete. On such a transcript the two terms describe a codon that the cache does not hold.

The module asks the transcript one question before it reads those positions as a stop codon: whether it carries the flag. It reads no base that the reference holds there. On a transcript that carries the flag, its two stop codon overlap tests return 0, and Ensembl VEP writes the row that vep-rs writes here. Section 2 shows it on transcripts of one gene that share a coding end, one of them with the flag.

**What vep-rs does instead.** vep-rs reads the coding sequence that the cache records, in its reading frame from its first base. A stop codon stands at the coding end where that sequence is a whole number of codons and its last codon is a stop codon of the transcript's codon table. Some coding sequences end one or two bases into their stop codon on the transcript's last base, and the adenines added to the mRNA, its polyadenylation, complete the codon. Where the transcript's cache record shows such a model, those one or two bases are the stop codon's; [`stop_codon_completed_by_polyadenylation`](stop-codon-completed-by-polyadenylation.md) says what shows it.

This entry holds five kinds of coding end on a transcript without the flag. At each, a change at the last three coding positions can meet bases that are no base of a stop codon:

1. Whole codons whose last codon names a residue, with or without a 3' UTR after them.
2. A coding sequence that ends one or two bases into a codon and is followed by a 3' UTR.
3. A coding sequence that ends one or two bases into a codon on the transcript's last base, where those bases followed by adenines are no stop codon.
4. A coding sequence that ends the same way where those bases followed by adenines would be a stop codon, and the transcript's cache record does not show a model completed that way.
5. A model whose stop codon is completed by polyadenylation. Its stop codon is the one or two bases after the last whole codon, so the rest of the last three coding positions are bases of that whole codon.

Where it meets no base of a stop codon, an insertion or deletion takes the terms that ask for no stop codon. A change of length that is no multiple of three inside the coding sequence is `frameshift_variant`, and a change of whole codons is `inframe_insertion` or `inframe_deletion`. A span that leaves the coding sequence, into the 3' UTR, into an intron or past the transcript's end, is `coding_sequence_variant` beside the terms of the bases it reaches.

Two readings of the fourth kind were defensible: that polyadenylation completes a stop codon there, or that none is shown. vep-rs's maintainers decided the second. Every source that asserts a stop codon completed by polyadenylation marks it on the record: the mitochondrial reference records, and one consensus coding sequence id. No source marks one on these models. In the Ensembl caches of GRCh38 and GRCh37 all but one are `non_stop_decay` models, and that biotype's definition speaks of a coding sequence with no stop codon before its poly(A) tail.

At an end of the fourth kind, an allele inside the last whole codon and a span that removes the last bases with the transcript's end take the rows of the other kinds. An allele that starts in those last bases keeps `incomplete_terminal_codon_variant`, which both engines write.

Where a recorded coding sequence stops just short of the stop codon of its own frame, the codon's bases are recorded as the first of the 3' UTR and are not read as the stop codon. A change in them stays `3_prime_UTR_variant` in both engines, as the maintainers decided.

**Where it is in Ensembl VEP's source.** This page calls the last three positions of the coding sequence the old end; after an edit, the bases at the old end are the ones that then stand at those positions. The lines are those of `ensembl-variation` at release/116, the module the Ensembl VEP 116.3 image holds, all in `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`:

- Lines 1363-1385, `_overlaps_stop_codon`: line 1373 returns 0 for a transcript that carries `cds_end_NF`. Lines 1378-1381 overlap the span's cDNA coordinates with the last three cDNA positions of the coding sequence. No line reads what the reference holds at those positions.
- Lines 1433-1484, `_ins_del_stop_altered`, which line 1263 of `stop_lost` asks where the allele has no usable alternate peptide: line 1466 makes the edit in the coding sequence joined to the 3' UTR. Line 1469 answers that the stop codon is altered where the result is shorter than the coding sequence. Lines 1473-1480 answer the same where the three bases at the old end of the edited sequence are no stop codon.
- Lines 1389-1431, `_overlaps_stop_codon_cil`, the same overlap on the genome: the flag at line 1400, an insertion's start moved by the inserted length at lines 1408-1415, and the last three genomic coding positions at lines 1418-1427.
- Line 1318, `stop_retained`: an allele that changes the length and passes the genomic overlap is `stop_retained_variant` where `_ins_del_stop_altered_cil` answers 0. That routine, lines 1488-1539, makes the edit and the test of `_ins_del_stop_altered` behind the genomic overlap. It answers 0 where the edit leaves a stop codon at the old end, and at line 1508 where the span lacks a cDNA coordinate or a CDS start.
- Line 1549, `frameshift`: the routine returns 0 where `stop_retained` holds. A change inside the coding sequence whose length would give `frameshift_variant` then carries `inframe_insertion`, `inframe_deletion` or `protein_altering_variant` in its place.

**What correcting it would change.** A test of the bases beside the test of the flag, at lines 1373 and 1400, would end the overlap where the positions the allele meets are no base of a stop codon. `stop_lost` and `stop_retained_variant` would leave these rows, and `frameshift` would not be withheld behind `stop_retained`. The rows would be the ones Ensembl VEP writes on a transcript with the flag.

**The class it produces.** `stop_codon_overlap_without_stop_codon`: a `Transcript` row of an insertion or deletion on which Ensembl VEP writes `stop_lost` or `stop_retained_variant` and vep-rs writes neither. The transcript carries no `cds_end_NF`, and its cache record shows no base of a stop codon at the positions the change meets. 1.7 says how a comparison reads that.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape the input record alone shows is an insertion or deletion, and no command-line option is needed. That page is the index of the documented differences, and it counts each shape over 15,646,481 chromosome 21 records. Of them 2,154,724 have an allele of another length than the reference, 137,713 per million. The rate is an upper bound on how often the difference occurs: the transcript decides whether such a record meets a coding end of these kinds.

## 1. A deletion in a last codon that names a residue, on ENST00000243562

**In short.** On GRCh37, the record of 1.2 deletes 19:41125378, the middle base of the last codon of ENST00000243562, a forward-strand transcript of LTBP4. Its coding sequence is 186 whole codons, the last one names valine, and a 3' UTR follows. Ensembl VEP 116.3 writes `frameshift_variant,stop_lost`, and vep-rs writes `frameshift_variant`. Ensembl VEP reads the last three coding positions as a stop codon because the transcript carries no `cds_end_NF`; the transcript holds the codon of a residue there, so no stop codon is lost.

### 1.1 What a reader needs to know first

A coding sequence is read in codons, three bases at a time from its first base, and that reading is its frame. A complete one ends with a stop codon, three bases that end translation.

Some transcripts are recorded without a stop codon. Where the 3' end of a coding sequence was not found, the cache marks the transcript with the attribute `cds_end_NF`, and a row prints it as `FLAGS=cds_end_NF`. Other transcripts carry no mark and still end in one of the first four ways listed above, and the fifth kind holds only the first base or bases of its stop codon.

Where an allele has reference and alternate amino acids, Ensembl VEP reads the two terms from them. Where the alternate is missing or is the unknown residue `X`, as after a change of length that breaks a codon, a stop codon test decides. The test edits the coding sequence and the 3' UTR as the record says and reads the three bases at the old end.

A cDNA coordinate is a base's position in the spliced transcript, and a CDS coordinate its position in the coding sequence. A base outside the transcript's exons has neither, and a span that ends on one prints `?` there, as in `349-?`.

### 1.2 The record

A deletion of one base, written for this example, against the Ensembl GRCh37 cache of release 116:

```
19	41125377	deletion_in_the_last_codon_19_41125378	GT	G	.	.	.
```

Both engines drop the shared first base, so the record deletes 19:41125378. ENST00000243562 is a forward-strand transcript of 13 exons of the gene LTBP4. Its parts on the spliced transcript:

| Part | cDNA positions | Bases |
| --- | --- | --- |
| 5' UTR | 1 to 268 | 268 |
| Coding sequence | 269 to 826 | 558 |
| 3' UTR | 827 to 2569 | 1,743 |

The coding sequence is 558 = 3 × 186 bases, 186 whole codons. Its last codon is coding bases 556 to 558, at 19:41125377-41125379, and it names valine. The 3' UTR begins with the next base, 19:41125380, in the same exon. The transcript's translation has 186 residues and ends with that valine. Read on in the same frame, the 3' UTR holds its first stop codon at its bases 1,474 to 1,476.

The transcript carries no `cds_end_NF`. The record deletes coding base 557, the codon's middle base, at cDNA 825. vep-rs reads the transcript from the committed GRCh37 cache, a cut of the Ensembl cache that Ensembl VEP reads whole. Both engines write the rows of 1.3 with no option beyond their cache options.

### 1.3 What each engine writes

The row on ENST00000243562, in the columns `Location`, `Feature`, `Consequence`, `cDNA_position`, `CDS_position`, `Protein_position` and `IMPACT`. Ensembl VEP 116.3:

```
19:41125378 | ENST00000243562 | frameshift_variant,stop_lost | 825 | 557 | 186 | HIGH
```

vep-rs:

```
19:41125378 | ENST00000243562 | frameshift_variant | 825 | 557 | 186 | HIGH
```

The rows differ in the `Consequence` alone. Every other column is equal, the residue and codon cells among them. Ensembl VEP writes 33 rows for the record and vep-rs 30: the committed cache lacks three transcripts, on which Ensembl VEP writes `upstream_gene_variant`. The other 29 shared rows are equal in every column.

### 1.4 Why Ensembl VEP writes that

The deletion breaks a codon, and the alternate amino acid is the unknown residue: the row's `Amino_acids` cell ends in `X`. `stop_lost` therefore asks `_ins_del_stop_altered` (lines 1259 and 1263).

That routine asks `_overlaps_stop_codon` first (line 1444). The transcript carries no `cds_end_NF`, so line 1373 does not return. The span is cDNA 825 to 825, and lines 1378-1381 overlap it with the last three cDNA positions of the coding sequence, 824 to 826. The overlap holds.

The routine then joins the coding sequence to the 3' UTR and removes the deleted base (line 1466). The result is longer than the coding sequence, so line 1469 does not answer. Lines 1473-1480 read the three bases now at the old end: the codon's first base, its last base and the first base of the 3' UTR. They spell a codon of valine and no stop codon, and the routine answers that the stop codon is altered.

`stop_lost` is written. `frameshift_variant` stands beside it: the length changes by one base (line 1564), and `stop_retained`, which would withhold the term at line 1549, returns 0 once `stop_lost` holds (line 1297).

### 1.5 Why vep-rs's row is the correct one

The record removes one of the 558 coding bases. The reference holds the codon of valine at coding bases 556 to 558: the coding sequence is whole codons, and the last of them names a residue.

`stop_lost` needs a stop codon that the change removes or alters. None stands at the old end before the change, so none is lost. What the record shows is one base deleted inside the coding sequence. The length changes by a number that is no multiple of three, which is `frameshift_variant`.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --species homo_sapiens --assembly GRCh37
grep ENST00000243562 out.txt
```

vep-rs writes `frameshift_variant` on the row. Ensembl VEP 116.3 takes the record with the cache options of [Reproducing these examples](reproducing.md) and its GRCh37 cache, and writes `frameshift_variant,stop_lost`.

The entry's reproduction is `scripts/divergences/reproductions/stop_codon_overlap_without_stop_codon.toml`, which `tools/vep-divergence-reevaluate` runs on both engines. vep-rs's cache there is a copy of the committed GRCh37 cache, which holds the transcripts of sections 1, 3, 4 and 5. One transcript file of `tests/golden/116/GRCh37-plugin-downstream-own-span/json_cache` is added for those of section 2.

### 1.7 How the comparison treats this class

Where a comparison names the class `stop_codon_overlap_without_stop_codon`, it sets aside a `Transcript` row of an insertion or deletion on which Ensembl VEP writes `stop_lost` or `stop_retained_variant` and vep-rs writes neither, its other terms being the ones this page gives for the change. The two rows do not show what the coding end holds, so the class reads the cache vep-rs annotated with, the candidate's cache: the transcript must carry no `cds_end_NF` there and hold no base of a stop codon where the change meets it. A comparison that is given no cache sets nothing aside under this class, and each such row stays a difference.

Under `--hgvs` the class also sets aside an `HGVSp` pair. Where the transcript's last base cuts a frameshift's first changed codon short, Ensembl VEP writes the deleted residues with an extension of unknown length and vep-rs the first residue with `fs`.

```text
Ensembl VEP 116.3   p.<First><n>delextTer?   p.<First><n>_<Last><m>delextTer?
vep-rs              p.<First><n>fs           p.<First><n>fs
```

Rows beside the shape are no part of this entry:

- A substitution at the old end. Both engines read its terms from its codons, and the rows are equal.
- An allele that starts in the one or two bases after the last whole codon, at an end of the second, third or fourth kind. Both engines write `incomplete_terminal_codon_variant` there.
- A row on a transcript whose stop codon stands. [`stop_codon_removed_without_stop_codon_test`](stop-codon-removed-without-stop-codon-test.md), [`stop_retained_without_stop_codon_test`](stop-retained-without-stop-codon-test.md), [`stop_retained_three_base_stop_codon_test`](stop-retained-three-base-stop-codon-test.md) and [`stop_codon_overlap_insertion_range`](stop-codon-overlap-insertion-range.md) describe such rows. The third has the pair of terms that an edit which leaves a stop codon at the old end takes here, `stop_retained_variant` against `frameshift_variant`.
- An allele that touches the base or bases of a stop codon completed by polyadenylation: [`stop_codon_completed_by_polyadenylation`](stop-codon-completed-by-polyadenylation.md).
- A structural deletion over such a coding end: [`structural_deletion_stop_lost_without_stop_codon`](structural-deletion-stop-lost-without-stop-codon.md).

Sections 2 to 4 give one further kind of coding end each, under the numbers of this section's parts: .2 the record, .3 the rows, .5 the reason. Section 6 gives the fourth kind.

## 2. A coding sequence that ends inside a codon before a 3' UTR, on ENST00000419804

**In short.** On GRCh37, the record of 2.2 deletes X:152938095, the last base of the last whole codon of ENST00000419804, a reverse-strand transcript of PNCK. Its coding sequence is 42 codons and two more bases, and a 3' UTR follows. Ensembl VEP writes `frameshift_variant,stop_lost`, and vep-rs writes `frameshift_variant`. On ENST00000423545, which has the same coding end and carries `cds_end_NF`, both engines write `frameshift_variant`.

### 2.2 The record

```
X	152938094	deletion_in_the_last_whole_codon_X_152938095	GG	G	.	.	.
```

The record deletes X:152938095. On this strand a transcript runs from high to low coordinates.

| Part of ENST00000419804 | Where | Bases |
| --- | --- | --- |
| Coding sequence | cDNA 141 to 268 | 128 |
| Last whole codon, of alanine | X:152938097-152938095 | 3 |
| The two bases after it | X:152938094 and X:152938093 | 2 |
| 3' UTR | cDNA 269 to 561, from X:152938092 | 293 |

The coding sequence is 128 = 3 × 42 + 2 bases. The last whole codon is coding bases 124 to 126, and the record deletes base 126, at cDNA 266. Bases 127 and 128 begin a codon that the first base of the 3' UTR completes. That codon names leucine. The old end, coding bases 126 to 128, is the last base of one codon and the first two of the next: no codon of the frame.

ENST00000425526 and ENST00000458354 hold the same coding sequence before a 3' UTR of their own. ENST00000423545 holds the same 128 coding bases at the same genomic positions. It ends at its last coding base, X:152938093, and carries `cds_end_NF`. The committed cache of the corpus `GRCh37-plugin-downstream-own-span` holds all four.

### 2.3 What each engine writes

In the columns of 1.3 with `FLAGS` added, Ensembl VEP 116.3:

```
X:152938095 | ENST00000419804 | frameshift_variant,stop_lost | 266 | 126 | 42 | HIGH | -
X:152938095 | ENST00000425526 | frameshift_variant,stop_lost | 623 | 126 | 42 | HIGH | -
X:152938095 | ENST00000458354 | frameshift_variant,stop_lost | 324 | 126 | 42 | HIGH | -
X:152938095 | ENST00000423545 | frameshift_variant | 446 | 126 | 42 | HIGH | cds_end_NF
```

vep-rs:

```
X:152938095 | ENST00000419804 | frameshift_variant | 266 | 126 | 42 | HIGH | -
X:152938095 | ENST00000425526 | frameshift_variant | 623 | 126 | 42 | HIGH | -
X:152938095 | ENST00000458354 | frameshift_variant | 324 | 126 | 42 | HIGH | -
X:152938095 | ENST00000423545 | frameshift_variant | 446 | 126 | 42 | HIGH | cds_end_NF
```

Three rows differ, each in the `Consequence` alone. The fourth, on the transcript with the flag, is equal.

### 2.5 Why the rows differ, and why vep-rs's are the correct ones

On the three transcripts without the flag, line 1373 does not return, and the span, the cDNA position of coding base 126, overlaps the last three cDNA positions of the coding sequence. The edit of line 1466 leaves coding bases 127 and 128 and the first base of the 3' UTR at the old end. They spell the codon of leucine and no stop codon, so lines 1473-1480 answer that the stop codon is altered. On ENST00000423545 line 1373 returns 0, and `stop_lost` is not written.

The four transcripts hold the same bases at the same coding positions, and the record removes the same base of the same codon of alanine. The flag is the only thing that differs. The coding sequence of the three ends two bases into a codon, and the cache gives the base that completes the codon as untranslated. No stop codon stands at the old end, and vep-rs writes one row on all four.

## 3. A coding sequence that ends inside a codon on the transcript's last base, on ENST00000586663

**In short.** On GRCh37, ENST00000586663, a reverse-strand transcript of PEX12, ends two bases into a codon on its last base, and adenines complete those bases to the codon of lysine. A deletion of the last base of its last whole codon is `frameshift_variant,stop_lost` in Ensembl VEP and `frameshift_variant` in vep-rs.

### 3.2 The record

```
17	33902806	deletion_in_the_last_whole_codon_17_33902807	TC	T	.	.	.
```

| Part of ENST00000586663 | Where | Bases |
| --- | --- | --- |
| Coding sequence | cDNA 532 to 1607 | 1,076 |
| Last whole codon, of glutamate | 17:33902809-33902807 | 3 |
| The two bases after it | 17:33902806 and 17:33902805 | 2 |
| 3' UTR | none | 0 |

The coding sequence is 1,076 = 3 × 358 + 2 bases, and the transcript ends with its last coding base, 17:33902805. Its biotype is `non_stop_decay`. The two bases after the last whole codon, with an adenine after them, spell the codon of lysine: polyadenylation completes no stop codon here. The transcript carries no `cds_end_NF`.

The record deletes 17:33902807, coding base 1074, the last base of the last whole codon.

### 3.3 What each engine writes

In the columns of 1.3, Ensembl VEP 116.3:

```
17:33902807 | ENST00000586663 | frameshift_variant,stop_lost | 1605 | 1074 | 358 | HIGH
```

vep-rs:

```
17:33902807 | ENST00000586663 | frameshift_variant | 1605 | 1074 | 358 | HIGH
```

### 3.5 Why the rows differ, and why vep-rs's are the correct ones

The overlap of cDNA 1605 with the last three cDNA positions, 1605 to 1607, holds. The transcript has no 3' UTR, so the edited sequence is shorter than the coding sequence. Line 1469 answers that the stop codon is altered, and no base is read.

The old end is the last base of the codon of glutamate and the two bases after it: no codon of the frame, and no stop codon. The deletion loses none. It changes the length by one base inside the coding sequence, and vep-rs writes `frameshift_variant`.

## 4. The codon before a stop codon completed by polyadenylation, on ENST00000361453

**In short.** On GRCh37, the coding sequence of ENST00000361453, the mitochondrial gene MT-ND2, ends one base into its stop codon, which the adenines of the poly(A) tail complete. The last three coding positions are that base and two bases of the codon of leucine before it. A deletion or an insertion inside that codon is `frameshift_variant,stop_lost` in Ensembl VEP and `frameshift_variant` in vep-rs.

### 4.2 The records

```
MT	5509	deletion_in_the_last_whole_codon_MT_5510	TA	T	.	.	.
MT	5509	insertion_in_the_last_whole_codon_MT_5509	T	TC	.	.	.
```

| Part of ENST00000361453 | Where | Bases |
| --- | --- | --- |
| Coding sequence, the whole transcript | MT:4470-5511 | 1,042 |
| Last whole codon, of leucine | MT:5508-5510 | 3 |
| The base after it | MT:5511 | 1 |

The transcript is one exon on the forward strand with no UTR at either end. Its coding sequence is 1,042 = 3 × 347 + 1 bases, translated under the vertebrate mitochondrial code, codon table 2. The base at MT:5511 followed by two adenines is a stop codon of that table. The reference record of the human mitochondrial genome, `NC_012920.1`, annotates the stop codon of this coding sequence as completed by the adenines added to the mRNA.

The stop codon of this model is therefore one genomic base, MT:5511, coding base 1042. The first record deletes MT:5510, coding base 1041, the last base of the codon of leucine. The second inserts one base between MT:5509 and MT:5510, inside that codon. Neither touches MT:5511.

### 4.3 What each engine writes

In the columns of 1.3, the deletion and then the insertion. Ensembl VEP 116.3:

```
MT:5510 | ENST00000361453 | frameshift_variant,stop_lost | 1041 | 1041 | 347 | HIGH
MT:5509-5510 | ENST00000361453 | frameshift_variant,stop_lost | 1040-1041 | 1040-1041 | 347 | HIGH
```

vep-rs:

```
MT:5510 | ENST00000361453 | frameshift_variant | 1041 | 1041 | 347 | HIGH
MT:5509-5510 | ENST00000361453 | frameshift_variant | 1040-1041 | 1040-1041 | 347 | HIGH
```

### 4.5 Why the rows differ, and why vep-rs's are the correct ones

Ensembl VEP overlaps both spans with the last three cDNA positions, 1040 to 1042. The deletion lies on one of them and the insertion between two of them, 1040 and 1041, so both overlaps hold.

The transcript has no 3' UTR. The deletion leaves a sequence shorter than the coding sequence, and line 1469 answers that the stop codon is altered. After the insertion, lines 1473-1480 read the three bases at the old end, coding base 1040, the inserted base and coding base 1041. They spell a codon of serine, and the answer is the same.

Two of the three positions Ensembl VEP reads as the stop codon are bases of the codon of leucine. The stop codon begins after them, at MT:5511, and neither record touches that base or the adenines that follow it. Each record changes the length by one base inside the codon before the stop codon. vep-rs writes `frameshift_variant` and no stop term.

## 5. The other shapes the reproduction holds

Six more records of the reproduction follow. Five stand at an end of the first kind, and the last at the end of section 3. ENST00000287968 of NBPF1, on the reverse strand, is 164 whole codons before a 3' UTR, and its exon ends one base after the coding end. ENST00000593412 and ENST00000595468, on the reverse strand of chromosome 21, are 117 and 12 whole codons with no 3' UTR. Their last codons name glycine and glutamate, and the last codon of ENST00000287968 names proline.

| Record | What it changes | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `deletion_into_the_3_prime_utr_19_41125379` | ENST00000243562: the last coding base and the first 3' UTR base | `stop_lost,3_prime_UTR_variant`, `HIGH` | `coding_sequence_variant,3_prime_UTR_variant`, `MODIFIER` |
| `deletion_into_the_intron_1_16900981` | ENST00000287968: the last codon, the 3' UTR base after it and the first intron base | `splice_donor_variant,stop_retained_variant,3_prime_UTR_variant`, `HIGH` | `splice_donor_variant,coding_sequence_variant,3_prime_UTR_variant`, `HIGH` |
| `deletion_in_the_last_codon_21_47347680` | ENST00000593412: the middle base of the last codon | `frameshift_variant,stop_lost`, `HIGH` | `frameshift_variant`, `HIGH` |
| `deletion_past_the_transcript_end_21_47347678` | ENST00000593412: the last codon and the base past the transcript's end | `stop_retained_variant,3_prime_UTR_variant`, `LOW` | `coding_sequence_variant,3_prime_UTR_variant`, `MODIFIER` |
| `insertion_before_the_last_codon_21_34947433` | ENST00000595468: one base inserted before the last codon, leaving a stop codon at the old end | `inframe_insertion,stop_retained_variant`, `MODERATE` | `frameshift_variant`, `HIGH` |
| `insertion_before_the_last_three_coding_bases_17_33902807` | ENST00000586663: one base inserted before the last base of the last whole codon, leaving a stop codon at the old end | `protein_altering_variant,stop_retained_variant`, `MODERATE` | `frameshift_variant`, `HIGH` |

Each pair of rows differs in the `Consequence` and, where the table shows it, in the `IMPACT` that follows the terms. Every other column is equal.

- The deletion into the 3' UTR is `stop_lost` by the path of 1.4. The three bases at the old end after the edit spell a codon of valine.
- The two deletions that leave the transcript have no cDNA coordinate at their far end. `_overlaps_stop_codon` returns 0 at line 1376, so `stop_lost` is not written. The genomic overlap holds, `_ins_del_stop_altered_cil` returns 0 at line 1508 before it edits anything, and line 1318 writes `stop_retained_variant`. Each record removes the codon of a residue, proline on ENST00000287968 and glycine on ENST00000593412. No stop codon is shown before the change or after it, and vep-rs writes `coding_sequence_variant`.
- Both engines write `3_prime_UTR_variant` on ENST00000593412, which has no 3' UTR, for a span that crosses the transcript's last base. The term is no part of the difference.
- The insertion on ENST00000595468 stands before cDNA 34, the first of the last three cDNA positions, and not between two of them. The cDNA overlap fails, and `stop_lost` is not written. An insertion's genomic span starts at the base after it, here 21:34947434. On the reverse strand, line 1413 moves that start back by the inserted length to 21:34947433, the genomic base of coding base 34, and the genomic overlap of lines 1418-1427 holds.
- The edit leaves the inserted base and the first two bases of the last codon at the old end, a stop codon. `_ins_del_stop_altered_cil` answers 0, line 1318 writes `stop_retained_variant`, and line 1549 withholds `frameshift_variant`. That stop codon stood nowhere before the change, so none is kept, and vep-rs writes `frameshift_variant`. The insertion on ENST00000586663 is read the same way.

## 6. A last base that adenines would complete, on ENST00000605986

**In short.** ENST00000605986, a reverse-strand `non_stop_decay` transcript of AK9 on chromosome 6, ends one base into a codon on its last base, and that base followed by two adenines would be a stop codon. Nothing in its cache record shows a stop codon completed that way, so the end is of the fourth kind. The gnomAD v2.1.1 and v4.1 datasets each hold one deletion of 17 bases there, 6:109965696-109965712 on GRCh37 and 6:109644493-109644509 on GRCh38. It removes the last five coding bases and twelve bases past the transcript's end.

Ensembl VEP 116.3 writes `stop_retained_variant,3_prime_UTR_variant`, as for section 5's deletion past the transcript's end, and vep-rs writes `coding_sequence_variant,3_prime_UTR_variant`. No committed cache of a human assembly holds a transcript of this kind, so the reproduction runs no record on one.
