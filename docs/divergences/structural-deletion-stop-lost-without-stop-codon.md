# `stop_lost` for a structural deletion over a coding end with no stop codon

**The defect in one sentence.** For a structural deletion, Ensembl VEP 116.3 writes `stop_lost` wherever the span holds one of the three genomic positions at the 3' end of a transcript's coding region. It tests no flag and reads no base, so the term also stands where the span takes no base of a stop codon. vep-rs writes the row's other terms there, without `stop_lost`.

Most often no stop codon stands at that end: the transcript carries `cds_end_NF`, or its coding sequence ends on a codon for an amino acid, or inside a codon that its record does not show completed to a stop codon. Otherwise the last one or two coding bases begin a stop codon that the poly(A) tail completes, and the span ends before them.

**Why this is a defect.** The term says that a stop codon is lost, and on these rows the record removes none the transcript shows. `cds_end_NF` says that the transcript's coding end was not found. On the others the cache's own coding sequence shows a codon for an amino acid, or part of a codon, where a stop codon would stand.

The same routine does more for an allele given with its bases. Its stop codon tests return 0 on a transcript that carries `cds_end_NF`, and where the allele has reference and alternate amino acids it looks for the stop among them. So a deletion over the end of such a flagged transcript is `stop_lost` when its record writes `<DEL>` and is not when its record spells the bases.

**What vep-rs does instead.** vep-rs asks first whether a stop codon stands at the transcript's coding end. It reads the answer from the cache record: the flag, the codon table and the coding sequence, in the frame its first base sets. A stop codon stands where the transcript carries no `cds_end_NF` and one of two holds:

1. The coding sequence is a whole number of codons and its last codon is a stop codon of the table. The stop codon is the last three coding positions.
2. The coding sequence ends one or two bases into a codon on the transcript's last base, those bases followed by adenines are a stop codon of the table, and the record shows the model complete. It shows that by a mitochondrial genetic code (table 2, 5, 9 or 13) or by a consensus coding sequence (CCDS) id. The stop codon is those one or two positions, and the poly(A) tail completes it, as on six of the human mitochondrial transcripts.

A transcript that carries `cds_end_NF` shows no stop codon even where its cached coding sequence ends on one. The record says that its coding end was not found, and Ensembl VEP's own tests for an allele given with its bases answer by the flag there too. Nine rows of the two 1000 Genomes datasets stand on such transcripts, and 1.6 gives a record over one.

A structural deletion is `stop_lost` where its span takes a genomic base of that stop codon. Where no stop codon stands, no span loses one. A span that ends in the last whole codon before a stop codon of the second kind takes no base of it either, though it holds one of Ensembl VEP's three positions. The row keeps the terms its other tests give: `feature_truncation`, `coding_sequence_variant` or a frameshift or in-frame term, and the UTR, intron and splice terms.

Two readings are vep-rs's own. A stop codon is read in the coding sequence the cache records: three bases that the cache gives as 3' UTR are not read as one. The coding sequence is read in the frame its first base sets, so last bases that spell a stop codon across two codons of that frame are not read as one; section 4 shows such an end.

vep-rs keeps the two tests Ensembl VEP makes before its overlap. The term's filter asks it only of a span that reaches the coding sequence, and line 1240 excludes a span that starts in an incomplete last codon. vep-rs does not make the second test on a transcript of the second kind, where the bases after the last whole codon are the stop codon's own.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the module the Ensembl VEP 116.3 image holds), in `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm` unless another file is named:

- Lines 1234-1288, `stop_lost`. Line 1240 returns 0 for an allele that starts in an incomplete last codon (`partial_codon`, lines 1588-1614). It is the one test the routine makes before its two parts.
- Lines 1268-1280, the structural variant part. Line 1269 returns 0 unless the variant is a deletion (`deletion`, lines 192-216, which also holds for a copy number loss). Lines 1271-1272 read the coding region's genomic start and end.
- Line 1275, for the forward strand, asks whether the span overlaps `coding_region_end - 2` to `coding_region_end`. Line 1278, for the reverse strand, asks the same of `coding_region_start` to `coding_region_start + 2`. No line of the part reads an attribute of the transcript or a base of its sequence.
- Lines 1373 and 1400, `_overlaps_stop_codon` and `_overlaps_stop_codon_cil`: the stop codon tests of an allele given with its bases return 0 on a transcript that carries `cds_end_NF`.
- `Bio/EnsEMBL/Variation/Utils/Constants.pm`, lines 749-766: the term's own description is a change of "at least one base of the terminator codon (stop)". Lines 757-759 are the term's filter: it is asked only of a variant that reaches the coding sequence.

**What correcting it would change.** A structural part that made the flag test of lines 1373 and 1400, and then held the span against the bases of the stop codon that the coding sequence shows, would take `stop_lost` off these rows. Their other terms would stay. In the two structural-variant sets that is 2,117 rows on GRCh38 and 596 on GRCh37.

**The class it produces.** `structural_deletion_stop_lost_without_stop_codon`: a `Transcript` row of a structural deletion whose Consequence holds `stop_lost` in Ensembl VEP and the same terms without it in vep-rs, where the span takes no base of a stop codon at the transcript's coding end. Either none stands there, or the span ends before the one or two bases that begin it. The rows differ in no other column: `feature_truncation` keeps `IMPACT` at `HIGH` in both. 1.7 says how a comparison recognises the row.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data carries the shape, and the divergence needs no flag or plugin. The shape is a structural deletion, or a copy number loss, whose span holds one of the last three coding positions of a transcript and no base of a stop codon. A span that covers the whole transcript is not of the shape: its row carries `transcript_ablation` and no coding term.

That page is the index of the documented differences, and it counts the part of a shape that the record alone shows over 15,646,481 chromosome 21 records. Of them 440, all of 1000 Genomes Phase 3, are structural records with an allele that names a deletion or a copy number of zero, 28.1 per million. The rate is an upper bound: the transcript decides whether the divergence then occurs.

Counts are rows of the shape that vep-rs writes without `stop_lost`; Ensembl VEP 116.3 writes the term on every one of them but one. The two structural-variant sets, 16 files for each assembly, hold 2,713: 2,117 on GRCh38 and 596 on GRCh37, rows of 737 records on 90 transcripts. Every one stands on a transcript that carries `cds_end_NF`. The sets' gnomAD-SV and 1000 Genomes files hold 106 of them, rows of 59 records, and their synthetic files the rest.

The two genome-wide 1000 Genomes datasets, whose records include symbolic deletions, hold 757 more: 507 of Phase 3 on GRCh37 and 250 of the high-coverage dataset on GRCh38. Of them 739 stand on transcripts that carry `cds_end_NF` and 18 on transcripts without the flag. On one of the 250 Ensembl VEP writes no row: it drops the record as longer than `--max_sv_size` ([Defect C](defect-c.md)).

The committed corpora hold 29 rows of the shape in their default-format outputs, all on flagged transcripts: 19 in `GRCh37`, 4 in `GRCh38`, 3 in `GRCh38-synonyms`, 2 in `GRCh38-input-region` and 1 in `GRCh38-pick-summary`.

Set aside and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion over the end of a KCNJ15 fragment that carries `cds_end_NF`

**In short.** `21 39671500 <DEL>` with `END=39671700`, on GRCh37, deletes 200 bases inside the coding exon of KCNJ15. The fragment ENST00000419868 carries `cds_end_NF` and ends inside the span, at 21:39671600. Ensembl VEP 116.3 writes `stop_lost,feature_truncation,coding_sequence_variant,3_prime_UTR_variant`, and vep-rs writes the same row without `stop_lost`. The fragment's last codon names Cys, and the gene's stop codon starts 609 bases after the span, on the complete transcript, where neither engine writes `stop_lost`.

### 1.1 What a reader needs to know first

A coding sequence ends with a stop codon, three bases that end translation. `stop_lost` says that a variant changes that codon, so that translation runs on past it.

Some transcripts of a cache are fragments: their coding sequence breaks off before the stop codon. Ensembl marks a transcript whose coding end was not found with the attribute `cds_end_NF`, and both engines print it in the row as `FLAGS=cds_end_NF`. Other transcripts carry no flag, yet their cached coding sequence ends on a codon for an amino acid, or one or two bases into a codon. Sections 2 to 4 take those.

A structural variant is a record whose allele is a class name, such as `<DEL>` or `<CN0>`, with its span's end in `END`. Both engines start the span one base after `POS`. The record gives no bases, so each consequence term is read from where the span lies on the transcript. `feature_truncation` says that the deletion removes part of the transcript.

### 1.2 The record

The record is hand-written, against the Ensembl GRCh37 cache of release 116:

```
21	39671500	kcnj15_fragment_end	N	<DEL>	.	.	SVTYPE=DEL;END=39671700
```

It deletes 21:39671501-39671700, 200 bases. Two forward-strand transcripts of KCNJ15 matter here:

| Transcript | Coding sequence | Its last codon |
| --- | --- | --- |
| ENST00000328656, complete | 375 codons and a stop codon | the stop codon, 21:39672309-39672311 |
| ENST00000419868, `cds_end_NF` | 139 codons | codon 139, 21:39671598-39671600, which names Cys |

Both coding sequences start at the same base, 21:39671184, in the same exon, so codon 139 is the same three bases on both. The fragment ends with that codon, at 39671600, and has no 3' UTR. The complete transcript reads on through 236 more codons to its stop codon: 375 - 139 = 236.

The span holds the fragment's last 100 bases, 39671501-39671600, and the 100 bases after them. The stop codon starts 609 bases after the span's last base: 39672309 - 39671700 = 609.

### 1.3 What each engine writes

The fragment's row in the default format. Ensembl VEP 116.3:

```
kcnj15_fragment_end	21:39671501-39671700	deletion	ENSG00000157551	ENST00000419868	Transcript	stop_lost,feature_truncation,coding_sequence_variant,3_prime_UTR_variant	422-?	318-?	106-?	-	-	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;OverlapBP=100;OverlapPC=3.63
```

vep-rs:

```
kcnj15_fragment_end	21:39671501-39671700	deletion	ENSG00000157551	ENST00000419868	Transcript	feature_truncation,coding_sequence_variant,3_prime_UTR_variant	422-?	318-?	106-?	-	-	-	IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;OverlapBP=100;OverlapPC=3.63
```

The rows differ in one column, the `Consequence`, by `stop_lost` alone. Both carry `IMPACT=HIGH` and `FLAGS=cds_end_NF`. Both carry `3_prime_UTR_variant`, which each engine writes for a span that runs past the last coding base; that term is no part of this entry.

On the complete transcript, ENST00000328656, both engines write `frameshift_variant,feature_truncation` over CDS positions 318-517. The record writes 19 rows in each engine, and three differ, each as the row above does. They are the rows of three fragments that carry `cds_end_NF` and end inside the span: ENST00000419868, and ENST00000398927 and ENST00000398928, which end at 39671653. The other 16 rows are equal in every column.

### 1.4 Why Ensembl VEP writes that

The allele is a structural variant's, so `stop_lost` takes its structural part (line 1268). The record is a deletion (line 1269), and the fragment has a coding region (lines 1271-1272). The fragment is on the forward strand, so line 1275 compares the span with the last three positions of its coding region, 39671598-39671600. The span, 39671501-39671700, holds all three, and the routine answers 1.

Nothing else is asked. The routine does not read the transcript's attributes, where `cds_end_NF` stands, and it does not read the three bases. Line 1240 does not stop it either: the span starts in codon 106 of 139 whole codons.

On the complete transcript the same line compares the span with 39672309-39672311. They do not overlap, and the routine answers 0.

### 1.5 Why vep-rs's row is the correct one

The record removes 39671501-39671700. The three bases that Ensembl VEP's test finds inside it, 39671598-39671600, are codon 139 of the gene's coding sequence. They name Cys on the fragment and on the complete transcript alike. The gene's stop codon is at 39672309-39672311, and 39672309 is above 39671700: no stop codon lies in the span.

The fragment's cache record agrees. Its coding sequence has 139 codons, its translation has 139 amino acids and ends with Cys, and `cds_end_NF` says that its coding end was not found. vep-rs reads that no stop codon stands at this coding end and writes no `stop_lost`.

What the deletion does is in the other terms: `feature_truncation` and `coding_sequence_variant` on the fragment's row, and `frameshift_variant` on the complete transcript's, in both engines.

### 1.6 Try it yourself

From the repository's root, the command as vep-rs takes it:

```
printf '##fileformat=VCFv4.2\n' > record.vcf
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >> record.vcf
printf '21\t39671500\tkcnj15_fragment_end\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=39671700\n' >> record.vcf
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --species homo_sapiens --assembly GRCh37
grep -E 'ENST00000(419868|328656)' out.txt
```

The committed GRCh37 cache holds both transcripts, and vep-rs prints the two rows without `stop_lost`. Ensembl VEP 116.3 takes the same record with the cache options of [Reproducing these examples](reproducing.md) and its GRCh37 cache, and prints `stop_lost` on the fragment's row.

A copy number loss reads the same way, here over the end of a fragment on the reverse strand:

```
21	44485650	cbs_fragment_end	N	<CN0>	.	.	SVTYPE=DEL;END=44485760
```

ENST00000441030, a fragment of CBS, carries `cds_end_NF` and ends at 21:44485726 on a codon that names Lys. Ensembl VEP 116.3 writes `stop_lost,feature_truncation,coding_sequence_variant,3_prime_UTR_variant` on its row, by lines 1269 and 1278. vep-rs writes the same row without `stop_lost`.

The entry's reproduction, `scripts/divergences/reproductions/structural_deletion_stop_lost_without_stop_codon.toml`, runs these two records, the records of sections 2 to 4 and two more on both engines. It holds ten rows to their recorded values. The two more:

```
MT	5490	nd2_before_the_tail	N	<DEL>	.	.	SVTYPE=DEL;END=5510
1	235331830	arid4b_flagged_stop	N	<DEL>	.	.	SVTYPE=DEL;END=235331850
```

The first is the record of 3.1. The second deletes the stop codon, 1:235331840-235331842, on which two reverse-strand transcripts of ARID4B end. ENST00000349213 carries `cds_end_NF`: Ensembl VEP 116.3 writes `stop_lost` among its row's terms, and vep-rs does not. ENST00000264183 carries no flag, and both engines write `stop_lost` on its row.

### 1.7 How the comparator treats this class

A comparison that names the class `structural_deletion_stop_lost_without_stop_codon` sets aside a `Transcript` row of a structural deletion whose two sets of terms are equal but for Ensembl VEP's `stop_lost`. It reads vep-rs's cache to see that the span takes no base of a stop codon at the coding end of the row's transcript. Either none stands there, by the flag or by the coding sequence, or the span ends before the one or two bases that begin it. A comparison that is given no cache sets nothing aside under this class, and such a row stays a difference.

Rows beside the shape are no part of this entry. On a transcript where a stop codon stands at the coding end, both engines write `stop_lost` for a structural deletion that takes a base of it. This class sets aside neither a row that differs in another term as well nor a row of an allele given with its bases. A start codon's terms on a structural variant are no part of this entry.

## 2. A coding sequence of whole codons whose last codon names an amino acid

**In short.** `19 41125370 <DEL>` with `END=41125390`, on GRCh37, deletes 20 bases of an LTBP4 exon. On ENST00000243562, which carries no flag, they are the last three of 186 codons and the first 11 bases of the 3' UTR. The last codon names Val. Ensembl VEP 116.3 writes `frameshift_variant,stop_lost,feature_truncation,3_prime_UTR_variant`, and vep-rs writes the same row without `stop_lost`.

### 2.1 What a reader needs to know first

A cache can hold a transcript whose coding sequence is a whole number of codons, ends on a codon for an amino acid and carries no `cds_end_NF`. A 3' UTR may follow the coding sequence, as here, or the transcript may end with it, as on the second record of 2.2.

### 2.2 The records

Both records are hand-written, against the Ensembl GRCh37 cache of release 116:

```
19	41125370	ltbp4_coding_end	N	<DEL>	.	.	SVTYPE=DEL;END=41125390
21	34947400	ap000304_end	N	<DEL>	.	.	SVTYPE=DEL;END=34947440
```

ENST00000243562 is a forward-strand transcript of LTBP4 with 13 exons. Its coding sequence has 558 bases, 186 codons: 3 × 186 = 558. Codon 186 is 19:41125377-41125379 and names Val. The 3' UTR starts at the next base, 41125380, inside the same exon, and has 1,743 bases. The first record deletes 41125371-41125390: codons 184 to 186 and 11 bases of the 3' UTR.

ENST00000595468 is a reverse-strand transcript with two exons and no UTR. Its coding sequence has 36 bases, 12 codons. On this strand the transcript runs from high to low coordinates, so codon 12 is 21:34947433-34947431; it names Glu, and the transcript ends with it. The second record deletes 34947401-34947440: the transcript's last 10 bases and the 30 bases after them.

### 2.3 What each engine writes

The two rows in the default format. Ensembl VEP 116.3:

```
ltbp4_coding_end	19:41125371-41125390	deletion	ENSG00000090006	ENST00000243562	Transcript	frameshift_variant,stop_lost,feature_truncation,3_prime_UTR_variant	818-837	550-?	184-?	-	-	-	IMPACT=HIGH;STRAND=1;OverlapBP=20;OverlapPC=0.12
ap000304_end	21:34947401-34947440	deletion	ENSG00000268265	ENST00000595468	Transcript	stop_lost,feature_truncation,coding_sequence_variant,3_prime_UTR_variant	27-?	27-?	9-?	-	-	-	IMPACT=HIGH;STRAND=-1;OverlapBP=10;OverlapPC=4.65
```

vep-rs:

```
ltbp4_coding_end	19:41125371-41125390	deletion	ENSG00000090006	ENST00000243562	Transcript	frameshift_variant,feature_truncation,3_prime_UTR_variant	818-837	550-?	184-?	-	-	-	IMPACT=HIGH;STRAND=1;OverlapBP=20;OverlapPC=0.12
ap000304_end	21:34947401-34947440	deletion	ENSG00000268265	ENST00000595468	Transcript	feature_truncation,coding_sequence_variant,3_prime_UTR_variant	27-?	27-?	9-?	-	-	-	IMPACT=HIGH;STRAND=-1;OverlapBP=10;OverlapPC=4.65
```

Each pair of rows differs in the `Consequence` alone, by `stop_lost`. Neither transcript carries a flag. The first pair carries `frameshift_variant` in both engines, a term that is no part of this entry. On LTBP4's complete transcript ENST00000308370 the first record lies inside the coding sequence, at CDS positions 3388-3407, and both engines write `frameshift_variant,feature_truncation`.

### 2.4 Why Ensembl VEP writes that

On ENST00000243562 line 1275 compares the span, 41125371-41125390, with the last three coding positions, 41125377-41125379. On ENST00000595468, a reverse-strand transcript, line 1278 compares the span, 34947401-34947440, with the first three genomic positions of the coding region, 34947431-34947433. Each span holds its three positions, and the routine answers 1 without reading what those bases are.

### 2.5 Why vep-rs's row is the correct one

On ENST00000243562 the three positions are codon 186, which names Val, and the cache's translation has 186 amino acids and ends with Val. The cache gives every base after 41125379 as 3' UTR. The first record removes three codons for amino acids and 11 untranslated bases, and no stop codon is among them. On the complete transcript the same 20 bases lie in codons 1130 to 1136 of 1,623. Its stop codon is 19:41135453-41135455 and starts 10,063 bases after the span's last base.

On ENST00000595468 the three positions are codon 12, which names Glu, and the translation has 12 amino acids and ends with Glu. The transcript ends there and has no base after the codon.

vep-rs reads on both records whole codons whose last is no stop codon, so neither row carries `stop_lost`.

### 2.6 Try it yourself

Save the two records of 2.2 after the two header lines of 1.6 and run the command of 1.6 on the file. `grep -E 'ENST00000(243562|595468)' out.txt` selects the two rows. The committed cache holds both transcripts, and the reproduction holds both rows.

## 3. A coding sequence that ends inside a codon on the transcript's last base

**In short.** `19 42470700 <DEL>` with `END=42470760`, on GRCh37, deletes the last 27 bases of the ATP1A3 transcript ENST00000441343 and the 33 bases after them. The transcript carries no flag, and its coding sequence ends two bases into a codon on its last base. Ensembl VEP 116.3 writes `stop_lost,feature_truncation,coding_sequence_variant,3_prime_UTR_variant`, and vep-rs writes the same row without `stop_lost`.

### 3.1 What a reader needs to know first

A coding sequence whose length is no multiple of three ends one or two bases into a codon. vep-rs reads a stop codon there only where those bases are the transcript's last, adenines would complete them to a stop codon and the record shows the model complete. On this transcript adenines complete the two bases to a codon for Lys, and the record shows no such model: its code is the standard one and it carries no CCDS id.

Where the record does show it, the class takes a deletion that ends before those bases. On the mitochondrial transcript ENST00000361453 the one last base, MT:5511, begins the stop codon. The record `nd2_before_the_tail` of 1.6 deletes MT:5491-5510 on GRCh37: Ensembl VEP 116.3 writes `frameshift_variant,stop_lost,feature_truncation` and vep-rs `frameshift_variant,feature_truncation`.

ENST00000441343 has the biotype `non_stop_decay`, GENCODE's class for a transcript whose poly(A) tail follows the coding sequence with no stop codon before it.

### 3.2 The record

The record is hand-written, against the Ensembl GRCh37 cache of release 116:

```
19	42470700	atp1a3_end	N	<DEL>	.	.	SVTYPE=DEL;END=42470760
```

ENST00000441343 is a reverse-strand transcript with 22 exons and no 3' UTR. Its coding sequence has 3,680 bases: 1,226 codons and two bases more, since 3 × 1,226 = 3,678. Codon 1,226 is 19:42470738-42470736 and names Glu. The two bases after it are 42470735 and 42470734, the transcript's last base.

The record deletes 42470701-42470760. In the transcript's direction the span starts at 42470760, coding position 3654, in codon 1,218. It runs through the transcript's last base and 33 bases beyond it.

### 3.3 What each engine writes

The row in the default format. Ensembl VEP 116.3:

```
atp1a3_end	19:42470701-42470760	deletion	ENSG00000105409	ENST00000441343	Transcript	stop_lost,feature_truncation,coding_sequence_variant,3_prime_UTR_variant	3792-?	3654-?	1218-?	-	-	-	IMPACT=HIGH;STRAND=-1;OverlapBP=27;OverlapPC=0.10
```

vep-rs:

```
atp1a3_end	19:42470701-42470760	deletion	ENSG00000105409	ENST00000441343	Transcript	feature_truncation,coding_sequence_variant,3_prime_UTR_variant	3792-?	3654-?	1218-?	-	-	-	IMPACT=HIGH;STRAND=-1;OverlapBP=27;OverlapPC=0.10
```

The rows differ in the `Consequence` alone, by `stop_lost`. The record writes five rows in each engine, and the other four are equal in every column.

### 3.4 Why Ensembl VEP writes that

Line 1240 asks whether the span starts in an incomplete last codon. It starts in codon 1,218, a whole codon, so the routine goes on. Line 1278 compares the span with the first three genomic positions of the coding region, 42470734-42470736. Those are the two bases after the last whole codon and the third base of codon 1,226. The span holds them, and the routine answers 1.

### 3.5 Why vep-rs's row is the correct one

The three positions are no codon of the transcript: one is the last base of the codon for Glu, and two start a codon the transcript does not finish. The cache's translation has 1,226 amino acids and ends with Glu.

No stop codon is begun there either. The transcript ends on those two bases, and what a poly(A) tail adds after them is adenines. Followed by an adenine the two bases name Lys, so they begin no stop codon. vep-rs reads that no stop codon stands at this coding end and writes no `stop_lost`.

### 3.6 Try it yourself

Save the record of 3.2 after the two header lines of 1.6 and run the command of 1.6 on the file. `grep ENST00000441343 out.txt` selects the row. The committed cache holds the transcript, and the reproduction holds the row.

## 4. A coding sequence that ends inside a codon before a 3' UTR

**In short.** `2 109276100 <DEL>` with `END=109276140`, on GRCh37, deletes 40 bases of an exon of LIMS1. On ENST00000434274, which carries no flag, the coding sequence ends inside the span, two bases into its codon 21, and a 3' UTR follows. Ensembl VEP 116.3 writes `stop_lost` among the row's terms, and vep-rs writes the same row without it.

### 4.1 What a reader needs to know first

A transcript's frame is set by its first coding base: bases 1 to 3 are codon 1, bases 4 to 6 codon 2, and so on. Three neighbouring bases that are not one codon of that frame are never read together in translation, whatever they spell.

### 4.2 The record

The record is hand-written, against the Ensembl GRCh37 cache of release 116:

```
2	109276100	lims1_coding_end	N	<DEL>	.	.	SVTYPE=DEL;END=109276140
```

ENST00000434274 is a forward-strand transcript with three exons and the biotype `nonsense_mediated_decay`. Its coding sequence has 62 bases: 20 codons and two bases more, since 3 × 20 = 60. Codon 20 is 2:109276122-109276124 and names Ala. The two bases after it are 109276125 and 109276126. The 3' UTR starts at 109276127, inside the same exon, and has 666 bases.

The record deletes 109276101-109276140: coding positions 37 to 62 and the first 14 bases of the 3' UTR. Coding position 37 is the first base of codon 13: 3 × 12 + 1 = 37.

### 4.3 What each engine writes

The row in the default format. Ensembl VEP 116.3:

```
lims1_coding_end	2:109276101-109276140	deletion	ENSG00000169756	ENST00000434274	Transcript	frameshift_variant,stop_lost,feature_truncation,3_prime_UTR_variant,NMD_transcript_variant	88-127	37-?	13-?	-	-	-	IMPACT=HIGH;STRAND=1;OverlapBP=40;OverlapPC=0.82
```

vep-rs:

```
lims1_coding_end	2:109276101-109276140	deletion	ENSG00000169756	ENST00000434274	Transcript	frameshift_variant,feature_truncation,3_prime_UTR_variant,NMD_transcript_variant	88-127	37-?	13-?	-	-	-	IMPACT=HIGH;STRAND=1;OverlapBP=40;OverlapPC=0.82
```

The rows differ in the `Consequence` alone, by `stop_lost`. The record writes 16 rows in each engine, and the other 15 are equal in every column.

### 4.4 Why Ensembl VEP writes that

The span starts in codon 13, a whole codon, so line 1240 lets the routine go on. Line 1275 compares the span with the last three coding positions, 109276124-109276126, which it holds, and the routine answers 1.

### 4.5 Why vep-rs's row is the correct one

The three positions are the third base of codon 20 and the two bases after it. Read together they spell a stop codon, but the transcript's frame never reads them together. Codon 20 names Ala, and the cache's translation has 20 amino acids and ends with Ala. Both engines number the row in that frame: CDS position 37 is protein position 13.

The two bases after codon 20 start a codon that the transcript's next base, the first of the 3' UTR, would complete to a codon for Glu. The cache gives that base as untranslated. The transcript's frame holds no stop codon at 109276124-109276126, so the deletion loses none there. vep-rs writes the row's other terms without `stop_lost`.

### 4.6 Try it yourself

Save the record of 4.2 after the two header lines of 1.6 and run the command of 1.6 on the file. `grep ENST00000434274 out.txt` selects the row. The committed cache holds the transcript, and the reproduction holds the row.
