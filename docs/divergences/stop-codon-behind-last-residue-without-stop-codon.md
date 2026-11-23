# A stop codon written as retained where no stop codon stood

**The defect in one sentence.** Where a change keeps every residue of a protein and puts a stop codon directly behind the last one, Ensembl VEP 116.3 writes `stop_retained_variant` without asking whether a stop codon stood there before. On a transcript whose record shows no stop codon at its coding end, with the `cds_end_NF` flag or without it, its row is `inframe_insertion,stop_retained_variant`. vep-rs writes `stop_gained,inframe_insertion` on that row.

**Why this is a defect.** `stop_retained_variant` says that a change touches the stop codon and a stop codon still stands in its place. The Sequence Ontology defines the term on the terminator codon, which remains (SO:0001567). A coding sequence whose last codon names a residue has no terminator codon at its end, and neither has one that ends inside a codon the record does not complete or one the cache marks as cut off at its 3' end. On such a transcript a stop codon behind the last residue is new: none stood there to be retained.

The row also loses the term that describes the change. Ensembl VEP writes `stop_gained` where the alternate amino acids hold a stop and the reference amino acids hold none, and these rows meet that test. The module withholds the term because it has answered `stop_retained_variant` first. The row's impact falls with it, from `HIGH` to `MODERATE`.

**What vep-rs does instead.** vep-rs reads the coding sequence that the cache records, in its reading frame from its first base. A stop codon stands at the coding end where that sequence is a whole number of codons and its last codon is a stop codon of the transcript's codon table. Where the sequence ends one or two bases into a codon on the transcript's last base, and the record shows a model whose stop codon polyadenylation completes, those one or two bases are the stop codon's ([`stop_codon_completed_by_polyadenylation`](stop-codon-completed-by-polyadenylation.md)). Where a stop codon stands, vep-rs makes Ensembl VEP's comparison and writes `stop_retained_variant` as Ensembl VEP does.

This entry holds the coding ends where none stands:

1. Whole codons whose last codon names a residue, with or without a 3' UTR after them.
2. A coding sequence that ends one or two bases into a codon and is followed by a 3' UTR.
3. A coding sequence that ends one or two bases into a codon on the transcript's last base, where the record shows no such model.
4. Any coding end of a transcript that carries `cds_end_NF`.

Bases that the record gives as 3' UTR are not read as the stop codon, even where they are the stop codon of the coding sequence's own frame. Such a transcript is of the first or second kind, and its row is the one this page derives.

There vep-rs does not read a stop behind an unchanged protein as a stop codon kept. `stop_retained_variant` is not written, nothing withholds `stop_gained`, and the row carries it beside the in-frame term.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the module the Ensembl VEP 116.3 image holds), all in `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`:

- Lines 1290-1324, `stop_retained`: line 1313 takes an allele whose alternate amino acids are defined, not empty and without the unknown residue `X`, and line 1315 hands it to `ref_eq_alt_sequence`. The comment at line 1314 names the case the branch was written for, an in-frame insertion of a stop just before the stop.
- Lines 1326-1361, `ref_eq_alt_sequence`: line 1330 reads the transcript's cached peptide. Line 1344 answers 1 where the allele starts past the end of that peptide and its alternate amino acids open with a stop. Lines 1346-1347 write the alternate amino acids into a copy of the peptide at the allele's protein positions. The first clause of line 1357 answers 1 where the copy equals the peptide up to the peptide's length and a stop comes directly after.
- The second clause of that test, line 1358, holds only where the reference amino acids hold a stop, and is no part of this entry.
- Neither answer reads the coding end. No line of the routine reads `cds_end_NF`, a base or the length of the coding sequence.
- Lines 1211-1232, `stop_gained`: line 1221 returns 0 where `stop_retained` holds. Line 1228 is the term's own test, that the alternate amino acids hold a stop and the reference amino acids none, and it is not reached then.

**What correcting it would change.** A test of the coding end beside lines 1344 and 1357, asking whether the transcript's coding sequence ends on a stop codon, would end both answers on these transcripts. `stop_retained_variant` would leave the rows, line 1221 would not return, and line 1228 would write `stop_gained`. On a transcript whose stop codon stands nothing would change.

**The class it produces.** `stop_codon_behind_last_residue_without_stop_codon`: a `Transcript` row on which Ensembl VEP writes `stop_retained_variant` from its comparison of the peptides and vep-rs writes `stop_gained` in its place. The transcript's record shows no stop codon at its coding end, with `cds_end_NF` or without it. 1.7 says how a comparison reads that.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data does not carry the shape. That page is the index of the documented differences, and it counts each shape over 15,646,481 chromosome 21 records and over the committed corpora. The shape is a record that keeps every residue and puts a stop codon directly behind the last one, on a transcript that shows no stop codon at its coding end. No row of this kind is seen in the measured datasets or in the committed corpora, and the records of this page were written for the entry.

## 1. Three bases that bring a stop codon behind a last codon of leucine, on ENST00000545697

**In short.** On GRCh37, the record of 1.2 inserts three bases inside the last codon of ENST00000545697, a forward-strand transcript of LTBP4 whose coding sequence is 899 whole codons before a 3' UTR. The codon still names leucine after the change, and a stop codon stands directly behind it. Ensembl VEP 116.3 writes `inframe_insertion,stop_retained_variant`, and vep-rs writes `stop_gained,inframe_insertion`. Ensembl VEP finds the protein unchanged with a stop behind it and calls a stop codon retained; the record holds none there before the change, so one is gained.

### 1.1 What a reader needs to know first

A coding sequence is read in codons, three bases at a time from its first base, and that reading is its frame. A complete one ends with a stop codon, three bases that end translation. `stop_gained` says that a change makes a stop codon where the reference has none. `stop_retained_variant` says that a change touches the stop codon and a stop codon still stands in its place.

Some transcripts are recorded without a stop codon. Where the 3' end of a coding sequence was not found, the cache marks the transcript with the attribute `cds_end_NF`, and a row prints it as `FLAGS=cds_end_NF`. Other transcripts carry no mark and still end in one of the first three ways listed above.

Ensembl VEP translates the codons a change touches, before the change and after it. The row prints the two results in its `Amino_acids` cell, as in `L/L*`, where `L` is leucine and `*` a stop. This page calls them the reference and the alternate amino acids. The cache also holds each transcript's whole translation, its peptide. A peptide never holds the stop that ends it: where a stop codon closes the coding sequence, translation drops it.

A cDNA coordinate is the position of a base in the spliced transcript, and a CDS coordinate is its position in the coding sequence. An insertion lies between two bases, and the row prints both, as in `2696-2697`.

### 1.2 The record

Three bases inserted, written for this example, against the Ensembl GRCh37 cache of release 116:

```
19	41133637	stop_codon_inserted_in_the_last_codon_19_41133637	T	TGTA	.	.	.
```

Both engines drop the shared first base, so the record inserts three bases between 19:41133637 and 19:41133638. ENST00000545697 is a forward-strand transcript of 31 exons of the gene LTBP4. Its parts on the spliced transcript:

| Part | cDNA positions | Bases |
| --- | --- | --- |
| 5' UTR | 1 to 1530 | 1,530 |
| Coding sequence | 1531 to 4227 | 2,697 |
| 3' UTR | 4228 to 4776 | 549 |

The coding sequence is 2,697 = 3 × 899 bases, 899 whole codons. Its last codon is coding bases 2695 to 2697, at 19:41133636-41133638, and it names leucine. The 3' UTR begins with the next base, 19:41133639, in the same exon.

The transcript's peptide has 899 residues and ends with that leucine. Read on in the same frame, the 3' UTR holds its first stop codon at its bases 277 to 279. The transcript carries no `cds_end_NF`.

The record inserts its three bases between coding bases 2696 and 2697, the codon's second and third base, at cDNA 4226 and 4227. The first inserted base repeats the codon's third base, so with the codon's first two bases it spells the codon of leucine again. The other two inserted bases and the codon's old third base, which follows them, spell a stop codon. The edited coding sequence is 2,700 bases: the 899 codons of the reference and a stop codon behind them. The committed GRCh37 cache holds the transcript, and both engines write the rows of 1.3 with no flag beyond their cache options.

### 1.3 What each engine writes

The row on ENST00000545697, in the columns `Location`, `Feature`, `Consequence`, `cDNA_position`, `CDS_position`, `Protein_position`, `Amino_acids` and `IMPACT`. Ensembl VEP 116.3:

```
19:41133637-41133638 | ENST00000545697 | inframe_insertion,stop_retained_variant | 4226-4227 | 2696-2697 | 899 | L/L* | MODERATE
```

vep-rs:

```
19:41133637-41133638 | ENST00000545697 | stop_gained,inframe_insertion | 4226-4227 | 2696-2697 | 899 | L/L* | HIGH
```

The rows differ in the `Consequence` and in the `IMPACT` that follows it: Ensembl VEP's `stop_retained_variant` stands against vep-rs's `stop_gained`. Every other column is equal. Both rows print the same amino acids at residue 899, leucine against leucine and a stop.

Ensembl VEP writes 24 rows for the record, five of them on transcripts the committed cache lacks. vep-rs writes 19, one for each transcript that cache holds there. The other 18 rows the outputs share are equal in every column.

### 1.4 Why Ensembl VEP writes that

The allele has amino acids on both sides, `L` and `L*`, and the alternate holds no unknown residue. `stop_lost` therefore reads the two (line 1259). The reference amino acids hold no stop, so no stop is lost (line 1260).

`stop_retained` runs next. The alternate amino acids are defined, not empty and without `X`, so line 1313 holds and line 1315 asks `ref_eq_alt_sequence`. The overlap with the last three coding positions is the other branch of that test, line 1318, and is not asked for this allele.

That routine reads the transcript's peptide, 899 residues (line 1330). The allele's protein positions begin and end at residue 899 (lines 1332-1333). Line 1344 asks whether the allele starts past the peptide's end, and it does not. Lines 1346-1347 copy the peptide and write `L*` in place of residue 899. The copy has 900 characters: the 899 residues and a stop.

Line 1349 takes the copy's first 899 characters, which equal the peptide. Line 1351 takes what is left, the stop. The first clause of line 1357 holds: the protein is unchanged up to its old length, and a stop comes directly after it. The routine answers 1, and `stop_retained_variant` is written.

`stop_gained` returns 0 at line 1221, because `stop_retained` holds, and its own test at line 1228 is not reached. `inframe_insertion` is written from the amino acids: the alternate begins with the reference (line 1127).

No line of this path reads the coding end. A peptide of 899 residues looks the same whether a stop codon follows the 899th codon or not, since translation drops a closing stop. The comparison finds a stop behind an unchanged protein, and it cannot tell a stop codon that was kept from one that the record makes.

### 1.5 Why vep-rs's row is the correct one

The reference holds the codon of leucine at coding bases 2695 to 2697: the coding sequence is whole codons, and the last of them names a residue. The transcript's peptide ends with that leucine, at residue 899. The base after the coding end is the first base of the 3' UTR, and the cache records it and the 548 bases after it as untranslated.

`stop_retained_variant` needs a stop codon that stands before the change and after it. None stands at the coding end before the change. The stop codon behind residue 899 is made by the record: two of its three inserted bases, and the codon's old third base after them.

What the record shows is three bases inserted inside the coding sequence that keep every residue and add a stop codon. The length changes by one whole codon, which is `inframe_insertion`. The alternate amino acids hold a stop and the reference amino acids hold none, which is `stop_gained` by Ensembl VEP's own test at line 1228.

Read on through the 3' UTR, the reference frame meets its first stop codon at bases 277 to 279 of the 3' UTR. The cache records those bases as untranslated, and the record's stop codon stands before them.

Ensembl VEP's own row agrees where its comparison is not reached. Section 3 writes one edit at two places on one transcript. At the second place an earlier test answers before the peptides are compared, and Ensembl VEP writes `stop_gained`.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --species homo_sapiens --assembly GRCh37
grep ENST00000545697 out.txt
```

vep-rs writes `stop_gained,inframe_insertion` on the row. Ensembl VEP 116.3 takes the record with the cache options of [Reproducing these examples](reproducing.md) and its GRCh37 cache, and writes `inframe_insertion,stop_retained_variant`.

The entry's reproduction, `scripts/divergences/reproductions/stop_codon_behind_last_residue_without_stop_codon.toml`, runs every record of this page on both engines and holds each row to both recorded values. It also holds a second record on ENST00000545697, which replaces the codon's last base by four bases as the record of section 2 does. Ensembl VEP writes `inframe_insertion,stop_retained_variant` for it and vep-rs `stop_gained,inframe_insertion`. Its cache is a copy of the committed GRCh37 cache. One transcript file of the committed cache `tests/golden/116/GRCh37-plugin-downstream-own-span/json_cache` is added to it, which holds the transcripts of section 2.

### 1.7 How the comparison treats this class

A comparison that names the class `stop_codon_behind_last_residue_without_stop_codon` sets aside a `Transcript` row on which Ensembl VEP writes `stop_retained_variant` from its comparison of the peptides and vep-rs does not, where vep-rs's terms are the ones this page derives. The two rows do not show what the coding end holds. The class therefore reads the cache vep-rs annotated with, and takes a row only where that cache shows no stop codon at the transcript's coding end, with `cds_end_NF` or without it. A comparison that is given no cache sets nothing aside under this class, and each such row stays a difference.

Rows beside the shape are no part of this entry:

- A row on a transcript whose stop codon stands. There a stop codon directly behind the last residue takes the place of the one that stood there, and vep-rs writes `stop_retained_variant` as Ensembl VEP does.
- An insertion between the last whole codon and the one or two bases after it. Both engines read it in an incomplete terminal codon and write the fourth row of section 3.
- An insertion or deletion at such a coding end on which Ensembl VEP writes `stop_lost` or `stop_retained_variant` from its overlap with the last three coding positions, which it asks only where the alternate amino acids are missing, empty or hold the unknown residue: [`stop_codon_overlap_without_stop_codon`](stop-codon-overlap-without-stop-codon.md).
- A substitution of one base in the last whole codon. Both engines read its terms from its codons, and the rows are equal.

## 2. A coding sequence that ends inside a codon before a 3' UTR, with the flag and without it, on ENST00000419804

**In short.** On GRCh37, the record of 2.2 replaces X:152938095 by four bases. That base ends the last whole codon of ENST00000419804, a reverse-strand transcript of PNCK whose coding sequence is 42 codons and two more bases before a 3' UTR. The four bases keep the codon's alanine and put a stop codon behind it. Ensembl VEP writes `inframe_insertion,stop_retained_variant` and vep-rs `stop_gained,inframe_insertion`, each on ENST00000423545 as well, which has the same coding end and carries `cds_end_NF`.

### 2.2 The record

```
X	152938095	last_whole_codon_base_replaced_by_four_bases_X_152938095	G	TTAT	.	.	.
```

The record replaces one base by four, a substitution of several bases. The reference base equals neither the first nor the last of the four, so neither engine drops a shared base, and the row's position is the one base replaced. On this strand a transcript runs from high to low coordinates.

| Part of ENST00000419804 | Where | Bases |
| --- | --- | --- |
| Coding sequence | cDNA 141 to 268 | 128 |
| Last whole codon, of alanine | X:152938097-152938095 | 3 |
| The two bases after it | X:152938094 and X:152938093 | 2 |
| 3' UTR | cDNA 269 to 561, from X:152938092 | 293 |

The coding sequence is 128 = 3 × 42 + 2 bases. The last whole codon is coding bases 124 to 126, and the record replaces base 126, at cDNA 266. In the transcript's direction the first of the four bases completes the codon of alanine as the replaced base did, and the other three are a stop codon. Coding bases 127 and 128 follow them. In the reference they begin a codon that the first base of the 3' UTR completes, and that codon names leucine.

ENST00000425526 and ENST00000458354 hold the same coding sequence before a 3' UTR of their own. ENST00000423545 holds the same 128 coding bases at the same genomic positions. It ends at its last coding base, X:152938093, and carries `cds_end_NF`. The committed cache of the corpus `GRCh37-plugin-downstream-own-span` holds all four.

### 2.3 What each engine writes

In the columns of 1.3 with `FLAGS` added, Ensembl VEP 116.3:

```
X:152938095 | ENST00000419804 | inframe_insertion,stop_retained_variant | 266 | 126 | 42 | A/A* | MODERATE | -
X:152938095 | ENST00000425526 | inframe_insertion,stop_retained_variant | 623 | 126 | 42 | A/A* | MODERATE | -
X:152938095 | ENST00000458354 | inframe_insertion,stop_retained_variant | 324 | 126 | 42 | A/A* | MODERATE | -
X:152938095 | ENST00000423545 | inframe_insertion,stop_retained_variant | 446 | 126 | 42 | A/A* | MODERATE | cds_end_NF
```

vep-rs:

```
X:152938095 | ENST00000419804 | stop_gained,inframe_insertion | 266 | 126 | 42 | A/A* | HIGH | -
X:152938095 | ENST00000425526 | stop_gained,inframe_insertion | 623 | 126 | 42 | A/A* | HIGH | -
X:152938095 | ENST00000458354 | stop_gained,inframe_insertion | 324 | 126 | 42 | A/A* | HIGH | -
X:152938095 | ENST00000423545 | stop_gained,inframe_insertion | 446 | 126 | 42 | A/A* | HIGH | cds_end_NF
```

All four rows differ, each in the `Consequence` and the `IMPACT`. The flag changes neither engine's row.

### 2.5 Why the rows differ, and why vep-rs's are the correct ones

The peptide of each transcript is the 42 residues of its whole codons: translation leaves out the two bases after them. The alternate amino acids are alanine and a stop at residue 42. Written into the peptide, they leave its 42 residues as they were with a stop after them, and the first clause of line 1357 holds.

`ref_eq_alt_sequence` reads no flag, so it answers the same on ENST00000423545. Ensembl VEP's two overlaps with the stop codon do read the flag. Lines 1373 and 1400 return 0 for a transcript that carries `cds_end_NF`, and no stop codon term is written from them there.

The four transcripts hold the same bases at the same coding positions. Their coding sequence ends two bases into a codon. On three of them the record gives the base that completes the codon as untranslated, and the fourth is marked as cut off. No stop codon stands at any of the four coding ends, so the stop codon behind alanine is gained on all four.

The in-frame term is `inframe_insertion` in both engines. Four bases take the place of one, and the alternate amino acids begin with the reference (line 1127).

## 3. A coding sequence that ends on the transcript's last base, on ENST00000600312, ENST00000593412 and ENST00000586663

**In short.** On GRCh37, three transcripts without `cds_end_NF` end with their last coding base. ENST00000600312 and ENST00000593412, on chromosome 21, are 38 and 117 whole codons, and ENST00000586663 of PEX12 is 358 codons and two more bases, which adenines complete to the codon of lysine. Three bases inserted in the last whole codon of each keep its residue and put a stop codon behind it. Ensembl VEP writes `inframe_insertion,stop_retained_variant` on each row, and vep-rs writes `stop_gained,inframe_insertion`.

### 3.2 The records

```
21	37417817	stop_codon_inserted_in_the_last_codon_21_37417817	G	GGTA	.	.	.
21	47347679	stop_codon_inserted_in_the_last_codon_21_47347679	T	TTAT	.	.	.
17	33902807	stop_codon_inserted_in_the_last_whole_codon_17_33902807	C	CTAC	.	.	.
17	33902806	stop_codon_inserted_after_the_last_whole_codon_17_33902806	T	TCTA	.	.	.
```

| Transcript and strand | Coding sequence | Last whole codon | After it |
| --- | --- | --- | --- |
| ENST00000600312, forward | 114 bases, 38 codons | 21:37417816-37417818, of arginine | the transcript ends |
| ENST00000593412, reverse | 351 bases, 117 codons | 21:47347681-47347679, of glycine | the transcript ends |
| ENST00000586663, reverse | 1,076 bases, 358 codons and two bases | 17:33902809-33902807, of glutamate | 17:33902806 and 17:33902805, where the transcript ends |

Each of the first three records inserts three bases between the second and third base of the last whole codon, as the record of 1.2 does. On ENST00000586663 that is between coding bases 1073 and 1074. GENCODE gives that transcript the biotype `non_stop_decay`, for a transcript whose poly(A) tail is attached to its coding sequence with no stop codon before it.

The fourth record writes the edit of the third one base further on, between coding bases 1074 and 1075. The three bases it inserts there are the stop codon that the third record makes, and the two records leave the same edited sequence.

### 3.3 What each engine writes

In the columns of 1.3, Ensembl VEP 116.3:

```
21:37417817-37417818 | ENST00000600312 | inframe_insertion,stop_retained_variant | 113-114 | 113-114 | 38 | R/R* | MODERATE
21:47347679-47347680 | ENST00000593412 | inframe_insertion,stop_retained_variant | 350-351 | 350-351 | 117 | G/G* | MODERATE
17:33902807-33902808 | ENST00000586663 | inframe_insertion,stop_retained_variant | 1604-1605 | 1073-1074 | 358 | E/E* | MODERATE
17:33902806-33902807 | ENST00000586663 | stop_gained,inframe_insertion,incomplete_terminal_codon_variant | 1605-1606 | 1074-1075 | 358-359 | -/* | HIGH
```

vep-rs:

```
21:37417817-37417818 | ENST00000600312 | stop_gained,inframe_insertion | 113-114 | 113-114 | 38 | R/R* | HIGH
21:47347679-47347680 | ENST00000593412 | stop_gained,inframe_insertion | 350-351 | 350-351 | 117 | G/G* | HIGH
17:33902807-33902808 | ENST00000586663 | stop_gained,inframe_insertion | 1604-1605 | 1073-1074 | 358 | E/E* | HIGH
17:33902806-33902807 | ENST00000586663 | stop_gained,inframe_insertion,incomplete_terminal_codon_variant | 1605-1606 | 1074-1075 | 358-359 | -/* | HIGH
```

The first three rows differ in the `Consequence` and the `IMPACT`. The fourth is equal in both engines.

### 3.5 Why the rows differ, and why vep-rs's are the correct ones

The first three rows take the path of 1.4. Each peptide ends with the last whole codon's residue, the alternate amino acids are that residue and a stop, and the first clause of line 1357 holds.

On the first two transcripts the last codon names a residue, and the transcript ends with it. On ENST00000586663 the two bases after the last whole codon begin a codon that the transcript does not complete, and adenines would complete it to the codon of lysine. No stop codon stands at any of the three coding ends, so each record's stop codon is gained.

The fourth row is Ensembl VEP's own reading of the third record's edit. Its insertion lies after the last whole codon, between residues 358 and 359, where the coding sequence holds two bases of a codon. `stop_retained` returns 0 for an allele in such an incomplete codon at line 1296, before it compares any peptide, and `stop_lost` does the same at line 1240. Nothing withholds `stop_gained`, and line 1228 writes it.

So Ensembl VEP gives one edited sequence two rows, by the place the record writes the edit at: a stop codon retained, and a stop codon gained. Only the second describes a coding sequence that held no stop codon.

## 4. Transcripts that carry `cds_end_NF`, on ENST00000451489 and ENST00000405964

**In short.** On GRCh37, ENST00000451489 of MAP3K7CL and ENST00000405964 of C21orf91 carry `cds_end_NF`: the cache records each as cut off at its 3' end. Three bases inserted in the last whole codon of each keep its residue and put a stop codon behind it. Ensembl VEP writes `inframe_insertion,stop_retained_variant` beside `FLAGS=cds_end_NF`, and vep-rs writes `stop_gained,inframe_insertion`.

### 4.2 The records

```
21	30547065	stop_codon_inserted_in_the_last_codon_21_30547065	A	AATA	.	.	.
21	19167526	stop_codon_inserted_in_the_last_whole_codon_21_19167526	C	CTAC	.	.	.
```

| Transcript and strand | Coding sequence | Last whole codon | After it |
| --- | --- | --- | --- |
| ENST00000451489, forward | cDNA 275 to 556, 282 bases, 94 codons | 21:30547064-30547066, of glutamate | the transcript ends |
| ENST00000405964, reverse | cDNA 170 to 896, 727 bases, 242 codons and one base | 21:19167528-19167526, of glutamine | 21:19167525, where the transcript ends |

Each record inserts three bases between the second and third base of the last whole codon, as the record of 1.2 does: between coding bases 281 and 282 of the first transcript, and between 725 and 726 of the second.

### 4.3 What each engine writes

In the columns of 1.3 with `FLAGS` added, Ensembl VEP 116.3:

```
21:30547065-30547066 | ENST00000451489 | inframe_insertion,stop_retained_variant | 555-556 | 281-282 | 94 | E/E* | MODERATE | cds_end_NF
21:19167526-19167527 | ENST00000405964 | inframe_insertion,stop_retained_variant | 894-895 | 725-726 | 242 | Q/Q* | MODERATE | cds_end_NF
```

vep-rs:

```
21:30547065-30547066 | ENST00000451489 | stop_gained,inframe_insertion | 555-556 | 281-282 | 94 | E/E* | HIGH | cds_end_NF
21:19167526-19167527 | ENST00000405964 | stop_gained,inframe_insertion | 894-895 | 725-726 | 242 | Q/Q* | HIGH | cds_end_NF
```

### 4.5 Why the rows differ, and why vep-rs's are the correct ones

The flag marks a coding sequence whose 3' end the annotation did not find, so the record holds no stop codon for the transcript. Ensembl VEP's two overlaps with the stop codon return 0 for the flag, at lines 1373 and 1400. `ref_eq_alt_sequence` is reached through line 1313 and reads no flag. It finds the peptide unchanged with a stop behind it and answers that a stop codon is retained, on a transcript that the same module reads elsewhere as having none.

The stop codon of each record stands behind the last residue the record holds, where the reference has no stop codon. vep-rs writes `stop_gained` beside `inframe_insertion`, the row the same change takes at every other coding end where no stop codon stands.
