# A splice term written from positions past the end of the record's span

**The defect in one sentence.** On a record whose alternate allele is longer than its reference allele, Ensembl VEP 116.3 compares the two alleles index by index to the longer allele's length. It places every index past the reference allele's end on the genome after the record's span. A splice donor, splice acceptor, donor fifth base, donor region, polypyrimidine tract or intron term is then written from bases the record does not replace. vep-rs does not write such a term where the record's span holds none of the bases the term names.

**Why this is a defect.** A record names the bases it replaces. Its position and the length of its reference allele give its span, and the alternate allele stands in the place of that span. The bases after the span stay as the reference has them, however long the alternate allele is.

Each of the six terms names bases of an intron. Five name bases at one of its ends: the two at its 5' end, the two at its 3' end, its fifth base, its third to sixth base, or the tract before its 3' end. The intron term names any base of the intron. This page calls those bases the term's window. A term that rests only on positions after the span states a change at bases the record leaves as they were.

**What vep-rs does instead.** vep-rs compares the alleles the same way and then holds each of the six terms to the record. A term is withdrawn where Ensembl VEP's test for it is met only past the end of the record's span and the span holds no base of a window of that term on the transcript. Where the span holds such a base, the term stays, as Ensembl VEP writes it, even where the two alleles hold the same base there. Where no differing index meets the test, a shared base earns no term in either engine; [Defect A](defect-a.md) counts such rows.

For `intron_variant` the window is the whole intron. The term stays where the span holds a base of the intron, even one of the two at either end. It is withdrawn where the span holds no base of the intron and Ensembl VEP's test is met only past the span.

The terms that follow from the span itself are written. Where no donor, acceptor, fifth base or donor region term is left on the row and a base of the span lies in the window of `splice_region_variant`, that term is written. Where the fifth base term is withdrawn and a base of the span lies in the donor region's window, `splice_donor_region_variant` is written. No other term of the row changes.

An allele that opens with its whole reference allele changes no base of its span. vep-rs reads it as the insertion it writes, between the span's last base and the base after it, and writes that insertion's splice and intron terms in place of the span's. Section 1.5 works both lengths.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the modules the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/Variation/VariationFeatureOverlapAllele.pm`, lines 405-460, `_get_differing_regions`: where both alleles hold more than one base (line 418), line 428 combines the two allele strings byte by byte, and line 429 lists the indices at which they differ. The shorter string counts as padded, so every index past its end differs.
- The same routine, lines 434-450: neighbouring indices are joined into regions, counted from index 0. No line holds a region's end to the reference allele's length. Lines 452-454 give every other record one region, from index 0 to the reference length minus one.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariationAllele.pm`, lines 117-118, `_intron_effects`: each region is placed on the genome at the record's start plus its two indices, with no bound at the record's end.
- The same routine, lines 143-150: a region that overlaps an intron between its third base and its third-last base sets the flag the intron term reads.
- The same routine, lines 176 and 180: a region that overlaps the first two or the last two bases of an intron sets the flag the donor and acceptor terms read. Lines 184 and 196 set the fifth base's flags, lines 188 and 200 the donor region's, and lines 155, 159, 192 and 204 the tract's.
- The same routine, lines 215-216: the flag of `splice_region_variant` is not set once a flag of line 176 or 180 is.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 560-614: the predicates of the five splice terms each pick a flag by the transcript's strand. Line 603 returns 0 for the donor region beside the fifth base, and lines 620-624 return 0 for `splice_region_variant` beside a donor, acceptor, donor region or fifth base term. Lines 629-634 hold the intron term's predicate, which reads the flag of line 149.

**What correcting it would change.** With the six terms held to the record's span in this way, Ensembl VEP would not write them on the rows where the two engines differ. On the genome-wide ClinVar dataset that is 2,856 terms on 1,787 rows on GRCh38 and 625 terms on 386 rows on GRCh37. They are 701 and 162 donor terms, 700 and 142 acceptor terms, 626 and 138 fifth base terms, 101 and 25 donor region terms, 133 and 38 tract terms, and 595 and 120 intron terms.

`splice_region_variant` would be written on 1,061 and 240 of those rows, and `splice_donor_region_variant` on 268 and 50. `IMPACT`, which is read from a row's most severe term, would follow on 1,077 and 238 rows.

**The class it produces.** `splice_term_past_the_record_span`. Every comparison program that names the class states it in the same sentence:

| Class | Its definition |
| --- | --- |
| `splice_term_past_the_record_span` | On a record whose alternate allele is longer than its reference allele, Ensembl VEP writes a splice donor, splice acceptor, donor fifth base, donor region, polypyrimidine tract or intron term from positions past the end of the record's span; where no base of the span lies in a window of that term, vep-rs does not write it, unless it is a term of an insertion written after the span. |

A row of the class is a `Transcript` row on which Ensembl VEP writes one of the six terms in this way and vep-rs does not, or the rows differ by a term of an insertion written after the span. Section 1.7 says which rows a comparison sets aside and which it leaves.

**How often the shape occurs.** Tier 1 on [Intended divergences from Ensembl VEP](../intended-divergences.md), the page that lists every documented difference and orders them by how often their input occurs: ordinary data carries the shape, and the divergence needs no flag or plugin. The record shows two shapes:

| Shape, as both engines read the record | Counted below |
| --- | --- |
| an alternate allele longer than a reference allele of two or more bases, an insertion written after the span among them | yes |
| in a record of two or more alternate alleles, an alternate allele that opens with a reference allele of one base | no: 796 further records, all of 1000 Genomes |

Both engines drop a first base that every allele of the record shares and, on a record of one alternate allele, the bases its two alleles share at either end.

That page counts each entry's shape over 15,646,481 chromosome 21 records of gnomAD, 1000 Genomes and ClinVar. Of them 442 have the first shape, 28.2 per million. That rate is an upper bound: the transcript decides whether a position past the end of the span lies in a window that the span does not reach. The terms change on 12 rows of those chromosome 21 records. On the genome-wide ClinVar dataset they change on 1,787 rows on GRCh38 and 386 on GRCh37.

Those 1,787 and 386 rows are rows of 90 and 80 records. A comparison that is given the input records sets all of them aside under the class. One that is not sets aside 558 and 128 of them and leaves 1,229 and 258 as differences. Section 1.7 says what the comparison reads.

On 63 and 20 rows Ensembl VEP's test for a term is met only past the end of the span, while the span holds a base of the term's window at which the two alleles hold the same base. The record's span reaches the window there, and both engines write the term. On 288 and 59 rows the same holds of the intron term: its test is met only past the span, the span itself holds one or two bases of the intron's end, and both engines write it.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A donor term from the intron's first base, ClinVar 817939 on ENST00000646124

**In short.** ClinVar 817939 on GRCh38 replaces two bases, 22:20996964-20996965, with four. On the forward-strand transcript ENST00000646124 of LZTR1 they are the third-last and second-last base of exon 20. Ensembl VEP 116.3 lays the four alternate bases along the genome from 20996964, and the fourth falls on 20996967, the intron's first base. It writes `splice_donor_variant,frameshift_variant`. The record replaces no base of the intron, and vep-rs writes `frameshift_variant,splice_region_variant`. A second record, ClinVar 267720, replaces the two donor bases of an intron of BRCA2. Both engines write `splice_donor_variant` for it, and Ensembl VEP alone adds `splice_donor_region_variant`, from the base after the span.

### 1.1 What a reader needs to know first

An intron is cut out of a transcript between two exons. In the transcript's orientation its 5' end follows one exon and its 3' end precedes the next. Seven Sequence Ontology terms name a change in an intron or near its ends:

- `splice_donor_variant` (SO:0001575): "A splice variant that changes the 2 base pair region at the 5' end of an intron."
- `splice_acceptor_variant` (SO:0001574): "A splice variant that changes the 2 base region at the 3' end of an intron."
- `splice_donor_5th_base_variant` (SO:0001787): "A sequence variant that causes a change at the 5th base pair after the start of the intron in the orientation of the transcript."
- `splice_donor_region_variant` (SO:0002170): "A sequence variant that falls in the region between the 3rd and 6th base after splice junction (5' end of intron)."
- `splice_polypyrimidine_tract_variant` (SO:0002169): "A sequence variant that falls in the polypyrimidine tract at 3' end of intron between 17 and 3 bases from the end (acceptor -3 to acceptor -17)."
- `intron_variant` (SO:0001627): "A transcript variant occurring within an intron."
- `splice_region_variant` (SO:0001630): "A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron."

This page calls the bases a term names its window. For `intron_variant` that is every base of the intron. Ensembl VEP writes the narrowest splice term: the donor region is not written beside the fifth base, and `splice_region_variant` is not written beside a donor, acceptor, donor region or fifth base term. It writes `intron_variant` beside any of them.

A VCF record gives a position, a reference allele and an alternate allele. The reference allele's bases lie on the genome from the position on, one position each. Those positions are the record's span, and the record replaces them with the alternate allele. A pure insertion is a record whose reference allele is empty once both engines drop the bases its two alleles share at either end. It replaces no base.

Index 0 is each allele's first base, index 1 its second, and so on. In a record whose alternate allele is the longer one, the alternate allele's last indices have no reference base of their own index. This page is about such a record whose reference allele still holds two or more bases once those shared bases are dropped. Line 418 of `VariationFeatureOverlapAllele.pm` sends a pair of alleles to the comparison of its lines 428-429 only where both hold more than one base. A one-base reference allele that opens a longer allele, the second shape, gets one region, itself (lines 452-454).

### 1.2 The record

Two records of ClinVar, against the Ensembl GRCh38 cache of release 116:

```
22	20996964	817939	AA	TTGG	.	.	.
13	32380145	267720	GGT	ACAG	.	.	.
```

Neither record is a pure insertion: the two alleles of each differ at index 0, and at every other index of the span. Both transcripts are on the forward strand, so each intron's 5' end is its lowest position.

| Part | ClinVar 817939 | ClinVar 267720 |
| --- | --- | --- |
| Transcript | ENST00000646124, LZTR1 | ENST00000380152, BRCA2 |
| Reference allele | 2 bases | 3 bases |
| Alternate allele | 4 bases | 4 bases |
| The record's span | 22:20996964-20996965 | 13:32380145-32380147 |
| Exon before the intron | exon 20, 22:20996886-20996966 | exon 24, 13:32380007-32380145 |
| Intron | 22:20996967-20997231 | 13:32380146-32394688 |
| Donor's two bases | 22:20996967-20996968 | 13:32380146-32380147 |
| Third to sixth base | 22:20996969-20996972 | 13:32380148-32380151 |
| Fifth base | 22:20996971 | 13:32380150 |
| Third base to third-last base | 22:20996969-20997229 | 13:32380148-32394686 |

The first record's span is two bases of exon 20: 20996964 is the exon's third-last base and 20996965 its second-last. The exon's last base, 20996966, and every base of the intron lie after the span.

The second record's span is three bases: 32380145, the last base of exon 24, and 32380146 and 32380147, the intron's first two bases. The intron's third base, 32380148, is the first position after the span.

### 1.3 What each engine writes

The row of the first record on ENST00000646124 in the default format. Ensembl VEP 116.3:

```
817939	22:20996964-20996965	TTGG	ENSG00000099949	ENST00000646124	Transcript	splice_donor_variant,frameshift_variant	2479-2480	2404-2405	802	K/LX	AAg/TTGGg	-	IMPACT=HIGH;STRAND=1
```

vep-rs:

```
817939	22:20996964-20996965	TTGG	ENSG00000099949	ENST00000646124	Transcript	frameshift_variant,splice_region_variant	2479-2480	2404-2405	802	K/LX	AAg/TTGGg	-	IMPACT=HIGH;STRAND=1
```

The rows differ in the `Consequence` column alone. Ensembl VEP writes `splice_donor_variant` where vep-rs writes `splice_region_variant`, each beside `frameshift_variant`.

The row of the second record on ENST00000380152. Ensembl VEP 116.3:

```
267720	13:32380145-32380147	ACAG	ENSG00000139618	ENST00000380152	Transcript	splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant	9455-?	9256-?	3086-?	-	-	-	IMPACT=HIGH;STRAND=1
```

vep-rs:

```
267720	13:32380145-32380147	ACAG	ENSG00000139618	ENST00000380152	Transcript	splice_donor_variant,coding_sequence_variant,intron_variant	9455-?	9256-?	3086-?	-	-	-	IMPACT=HIGH;STRAND=1
```

Both rows carry `splice_donor_variant`. Ensembl VEP's row also carries `splice_donor_region_variant`.

The repository holds a cut of the GRCh38 cache, `tests/golden/116/GRCh38/json_cache`, with the transcripts of both records. With that cache the first record writes 28 rows and the second 14. The first record's rows, by kind of transcript:

| Rows | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| 7, coding | `splice_donor_variant,frameshift_variant` | `frameshift_variant,splice_region_variant` |
| 8, 3' UTR of a decay transcript | `splice_donor_variant,3_prime_UTR_variant,NMD_transcript_variant` | `splice_region_variant,3_prime_UTR_variant,NMD_transcript_variant` |
| 6, non-coding exon | `splice_donor_variant,non_coding_transcript_exon_variant` | `splice_region_variant,non_coding_transcript_exon_variant` |
| 7, others | equal | equal |

On the 14 rows of the second and third kind, `IMPACT` is `HIGH` in Ensembl VEP's row and `LOW` in vep-rs's, the impact of each row's most severe term. The seven other rows are six `downstream_gene_variant` rows and one `intron_variant` row, of a transcript that holds the record inside a longer intron.

All 14 rows of the second record carry `splice_donor_variant` in both engines and `splice_donor_region_variant` in Ensembl VEP's alone. Their other terms and their `IMPACT`, `HIGH`, are equal.

### 1.4 Why Ensembl VEP writes that

Both alleles of the first record hold more than one base, so line 418 of `VariationFeatureOverlapAllele.pm` sends them to the comparison of lines 428-429. The reference allele has indices 0 and 1, and the alternate allele indices 0 to 3. Indices 0 and 1 hold different bases. Indices 2 and 3 lie past the reference allele's end and count as differing. Lines 434-450 join the four into one region, index 0 to 3.

Lines 117-118 of `BaseTranscriptVariationAllele.pm` place that region from the record's start: 20996964 + 0 to 20996964 + 3, which is 20996967. The span ends at 20996965. The placed region runs two positions past it, onto 20996966, the exon's last base, and 20996967, the intron's first.

Line 176 tests the region against the intron's first two bases, 20996967-20996968. The region's end, 20996967, lies on the first of them, and the flag of the intron's start is set. On strand 1 the donor's predicate reads that flag (`VariationEffect.pm` lines 567-568), and the row carries `splice_donor_variant`. The flag also keeps line 215 from setting the flag of `splice_region_variant`, and line 620 of `VariationEffect.pm` returns 0 for that term.

Line 144 tests the region against the intron from its third base, 20996969, to its third-last. The region ends two positions before 20996969, the intron term's flag is not set, and the row carries no `intron_variant` in either engine.

The second record's alleles differ at indices 0 to 2, and index 3 lies past the reference allele's end. The one region, index 0 to 3, is placed on 32380145 to 32380148. Line 176 finds the intron's first two bases, 32380146-32380147, in it. They are positions of the span, indices 1 and 2.

Line 188 tests the region against the intron's third to sixth base, 32380148-32380151. The region's end, 32380148, lies on the first of them, and the donor region's flag is set from index 3 alone. Line 184 tests the fifth base, 32380150, which the region does not reach. Line 603 of `VariationEffect.pm` therefore lets the donor region stand, and the row carries `splice_donor_region_variant` beside `splice_donor_variant`.

Line 144 tests the region against the intron from its third base, 32380148, to its third-last. The region's end lies on 32380148. The intron term's flag is set from index 3 alone as well, and the row carries `intron_variant`.

### 1.5 Why vep-rs's row is the correct one

The first record replaces positions 20996964 and 20996965 and no other. In the edited sequence its four alternate bases stand between 20996963 and 20996966:

| Positions | In the edited sequence | On ENST00000646124 |
| --- | --- | --- |
| up to 20996963 | the reference's bases | exon 20 |
| 20996964-20996965 | the four alternate bases | exon 20, third-last and second-last base |
| 20996966 | the reference's base | exon 20, last base |
| 20996967-20996968 | the reference's bases | the intron's first two bases |

`splice_donor_variant` is the term of a variant that "changes the 2 base pair region at the 5' end of an intron". The intron's first two bases are 20996967 and 20996968, and the record changes neither. No base of the span lies in the donor's window, so the term does not follow from the record.

The span does hold two of the exon's last three bases, 20996964 and 20996965. A change "within 1-3 bases of the exon" is what `splice_region_variant` states, and vep-rs writes it. The exon gains two bases, which `frameshift_variant` states in both engines' rows.

The second record replaces 32380145 to 32380147. The span holds both bases of the donor's window, 32380146 and 32380147, so `splice_donor_variant` follows from the record, and both engines write it. The donor region's window is 32380148-32380151. Every one of those positions lies after the span and holds the reference's base in the edited sequence. No base of the span lies in that window, and vep-rs does not write `splice_donor_region_variant`.

`intron_variant` is the term of a variant "occurring within an intron", and its window is the whole intron. The span holds two bases of the intron, 32380146 and 32380147. The term therefore follows from the record, and vep-rs writes it as Ensembl VEP does. Ensembl VEP's own test for it is met only at 32380148, past the span.

**What stays alike.** Three kinds of record keep their terms in both engines. The first is a record whose span itself holds a base of a window of the term, as the second record's span holds the donor's two bases, which are also bases of the intron. Both engines write the term, and nothing differs.

The second is a pure insertion. It replaces no base, and its terms rest on the two bases beside the inserted ones. Where the inserted bases lie between an exon's edge and a donor or acceptor site, or between such a site and the rest of the intron, Ensembl VEP writes `splice_region_variant`. Its source names that case in a comment (`BaseTranscriptVariationAllele.pm` lines 212-214) and tests it at `VariationEffect.pm` lines 97-100.

For the tract's term Ensembl VEP tests an insertion at the two bases beside it, and either one inside the tract's window is enough (`BaseTranscriptVariationAllele.pm` lines 152-161). The donor, acceptor, fifth base and donor region terms need both bases inside their window. That is Ensembl VEP's stated handling of insertions. vep-rs writes the same terms, and a pure insertion is no part of this difference.

The third is a record whose two alleles have the same length. No index of it lies past the reference allele's end, every region lies inside the span, and such a record is no part of this difference either. Inside the span both engines compare the two alleles index by index: a base of a window at an index where they hold the same base draws no term, and two changes with a shared base between them are read apart.

**An insertion written after the span.** One such record of a committed corpus, then the terms that can differ:

| Part | Line 166 of `tests/golden/116/TAIR10-arabidopsis/variants.vcf` |
| --- | --- |
| Record | 4:10523299, a reference allele of three bases and two alternate alleles |
| Span | 4:10523300-10523301, once both engines drop the first base, which every allele shares |
| Second alternate allele | the span's two bases, then two further bases |
| Transcript | AT4G19230.2, forward strand, with the exon 4:10523181-10523301 and the intron 4:10523302-10523411 |
| The two further bases | between 4:10523301, the exon's last base, and 4:10523302, the intron's first |
| Ensembl VEP | `splice_donor_variant,frameshift_variant`: it places the two bases on 4:10523302-10523303, the donor's two bases |
| vep-rs | `frameshift_variant,splice_region_variant` |
| The same insertion as a record of one alternate allele | `frameshift_variant,splice_region_variant` in both engines |

| Part | The fourth record of the entry's reproduction |
| --- | --- |
| Record | 22:20996967, a reference allele of two bases and two alternate alleles |
| Span | 22:20996968, once both engines drop the first base, which every allele shares |
| First alternate allele | the span's one base, then three further bases |
| Second alternate allele | none of the span: it deletes the base at 22:20996968, and both engines write `splice_donor_variant` for it |
| Transcript | ENST00000646124, forward strand, with the intron 22:20996967-20997231 |
| The three further bases | between 22:20996968, the donor's second base, and 22:20996969, the intron's third |
| Ensembl VEP | `splice_donor_variant`, with `IMPACT` `HIGH`: it does not compare the two alleles and reads the base at 22:20996968 as changed, and that base lies before the intron's third base, where the intron term's window starts |
| vep-rs | `splice_region_variant,intron_variant`, with `IMPACT` `LOW`: the insertion stands before the intron's third base, where Ensembl VEP's reading of an insertion sets both terms |
| The same insertion as a record of one alternate allele | `splice_region_variant,intron_variant` in both engines |

| vep-rs's row beside Ensembl VEP's | Terms |
| --- | --- |
| Can hold, from the two bases beside the insertion | `splice_region_variant`, `splice_donor_region_variant`, `splice_polypyrimidine_tract_variant`, `intron_variant`; where the kept base is an intron of one base, `splice_donor_variant` or `splice_acceptor_variant`, the term both engines write for the same insertion as a record of one alternate allele |
| Can lack, where only the bases after the span give it | each of the six terms |
| Can lack, where only a one-base span that the allele keeps gives it | each of the six terms and `splice_region_variant` |

### 1.6 Try it yourself

With the two records of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38/json_cache \
  --species homo_sapiens --assembly GRCh38
grep -E 'ENST00000646124|ENST00000380152' out.txt
```

The committed GRCh38 cache holds the transcripts of both records, and no flag beyond the cache is needed. vep-rs writes `frameshift_variant,splice_region_variant` on the first row and `splice_donor_variant,coding_sequence_variant,intron_variant` on the second. Ensembl VEP 116.3 takes the same records with the cache options of [Reproducing these examples](reproducing.md) and its GRCh38 cache. It writes `splice_donor_variant,frameshift_variant` on the first row and adds `splice_donor_region_variant` to the second.

The entry's reproduction, `scripts/divergences/reproductions/splice_term_past_the_record_span.toml`, runs four records on both engines and holds fourteen rows to their recorded values:

| Record | Rows held |
| --- | --- |
| ClinVar 817939 | one coding, one decay and one non-coding row, and its `intron_variant` row, which is equal in both engines |
| ClinVar 267720 | two rows |
| A record of two alternate alleles at 22:20997227, the first an insertion written after the span 22:20997228-20997229 | that allele's four rows: its three further bases stand between the tract's last base and the acceptor's first, Ensembl VEP writes `splice_acceptor_variant`, and vep-rs writes `splice_region_variant`, `splice_polypyrimidine_tract_variant` and `intron_variant` |
| A record of two alternate alleles at 22:20996967, the first an insertion written after the one base 22:20996968 | that allele's four rows: on three Ensembl VEP writes `splice_donor_variant` and vep-rs writes `splice_region_variant` and `intron_variant`, and the fourth is `intron_variant` in both engines |

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) compares one Ensembl VEP output with one vep-rs output, every field of every row. A comparison applies a class when its configuration lists the class's name. Where `splice_term_past_the_record_span` is listed, the comparator sets aside the rows described below and counts each once under that name. The class's definition is the sentence in the table at the top of this page. The SNP/indel comparator (`tools/vep-outputs-compare`) and the golden corpus builder (`tools/vep-golden-corpus`) hold the class under the same name and definition.

The class's predicate, `splice_term_past_the_record_span_pair`, reads the two engines' rows and the record's two alleles. A row's first column gives the alleles where it spells them, and the VCF and JSON output formats hold them in fields of their own.

Where the first column holds an identifier, as in the rows of 1.3, the comparison reads the alleles from the input records. The SNP/indel comparator takes those as `--input-vcf` or finds them under `--input-dir`, and the format comparator takes them as `--input`. The structural-variant comparator finds them under `--input-dir`, the rows of each output read from the file it was written from. The golden corpus builder reads its corpus's own records, and for this class no comparison reads a cache.

With the alleles, the predicate takes a `Transcript` row of a record that is no pure insertion. The record's reference allele holds a base for each position of the row's `Location`, two or more unless the allele opens with it, and its alternate allele is longer. Ensembl VEP's `Consequence` holds one or more of the six terms that vep-rs's does not hold. vep-rs's may hold `splice_region_variant` or `splice_donor_region_variant` where Ensembl VEP's does not. The two sets are otherwise equal, and so are the two rows' `cDNA_position` and `CDS_position`.

| On a row of an insertion written after the span | The predicate reads |
| --- | --- |
| The record | an allele that opens with its whole reference allele, of one base or more, which the record's alleles alone show: without them no row is read this way |
| vep-rs's `Consequence` | the terms of the insertion between the span's last base and the base after it: it may also hold `splice_polypyrimidine_tract_variant` and `intron_variant` where Ensembl VEP's does not, and after a reference allele of one base `splice_donor_variant` or `splice_acceptor_variant`, where Ensembl VEP's holds none of the six and no `splice_region_variant` |
| Ensembl VEP's `Consequence` | may hold any of the six terms where vep-rs's does not, and need hold none: with the record's alleles, a pair that an earlier rule only counts is put to this predicate too |
| Ensembl VEP's `Consequence` after a reference allele of one base | may also hold `splice_region_variant` where vep-rs's does not: Ensembl VEP reads that base as changed |
| Either `Consequence` printed under `--terms display` | the three names the seven terms print there, `ESSENTIAL_SPLICE_SITE`, `SPLICE_SITE` and `INTRONIC`, read with the record's alleles alone; where no corpus manifest documents the row, the `IMPACT` beside it stays a difference |

The comparison is not given where a transcript's introns lie. It places an intron beside the record's span from what the two rows and the record's alleles show, at each place the rows' `cDNA_position` allows. It takes the row where, at one such place, the lines of Ensembl VEP's source cited above give Ensembl VEP's terms and the rule of this page gives vep-rs's.

That has a limit. Where more than one such place gives Ensembl VEP's row, the comparison cannot tell them apart. On some rows where vep-rs keeps one of the six terms, a row from which that kept term were missing would therefore be set aside as well. On the genome-wide ClinVar dataset that is so for 1,280 of the terms of the six that vep-rs keeps on GRCh38, and for 278 on GRCh37. Those kept terms are counted on the rows of records with a longer alternate allele that are no pure insertion: on each of those rows Ensembl VEP or vep-rs writes one of the five splice terms or `splice_region_variant`.

Given the input records, the comparison sets aside every row on which the terms change: 1,787 rows on GRCh38 and 386 on GRCh37 of the genome-wide ClinVar dataset. It leaves none of them as a difference. Of the two records of this page, it sets aside the 21 rows of the first that differ and the 14 of the second. On a row it sets aside, the format comparator also sets aside the `IMPACT` that follows the two sets of terms, by the rule `term_pair_following_cell`.

Where no record can be read, the first column holding an identifier and no input records given, the predicate takes only what the two rows decide alone. The `Location` then holds three or more positions, `cDNA_position` gives a number for each end of the span, as far apart as the span's two ends, and vep-rs's row holds none of the six terms. A comparison run that way sets aside 558 rows on GRCh38 and 128 on GRCh37, and it leaves 1,229 and 258 as differences.

Of the rows left that way, 173 and 31 have a `Location` of two neighbouring positions, which the row of an insertion between two bases shares. On 134 and 27 the span crosses an exon's end, and on 922 and 200 it lies in an intron. The rows of this page's two records are among them: the first record's `Location` is two neighbouring positions, and the second's span crosses an exon's end.

With or without the record, the comparison leaves every other row and reports a difference on it. It leaves a row of a pure insertion, a row of a record of equal lengths, and a row of a record whose reference allele is the longer one. It leaves a row whose term has a window base inside the record's span, within the limit stated above. It also leaves a row that differs in one more term or one more column than the class names.

One exception stands. A row that also differs by the terms of a stop codon class, such as [`stop_codon_removed_without_stop_codon_test`](stop-codon-removed-without-stop-codon-test.md), is set aside where that class's rule takes the rest of the difference. One row of the genome-wide ClinVar dataset on GRCh38 is of that kind.
