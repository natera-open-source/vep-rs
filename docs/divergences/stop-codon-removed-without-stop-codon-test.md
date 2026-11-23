# `stop_retained_variant` for a span that holds the whole stop codon

**The defect in one sentence.** For an insertion or deletion whose span holds the whole stop codon and lacks a cDNA coordinate, at one end or at both, or a CDS start, Ensembl VEP 116.3 makes no stop codon test. It writes `stop_retained_variant`, which says that a stop codon stays in place, without reading a base. vep-rs reads the sequence the record gives and writes `stop_lost` unless that sequence shows the stop codon kept.

**Why this is a defect.** The term states the result of a test that was not made. Two routines of the module return 0 before they build any sequence, because the span has no coordinate to make the edit at. Their callers read the first 0 as a stop codon that is not lost and the second as a stop codon that is kept.

The same module answers the other way elsewhere when it has no sequence to read. A stop codon test whose edited sequence is shorter than the coding sequence answers that the codon is altered. A deletion given as a structural variant is `stop_lost` where its genomic span overlaps the stop codon, and no coordinate is asked for.

**What vep-rs does instead.** vep-rs reads the sequence the record gives: the transcript's own exon bases with the record's edit made in them. The edit takes out the bases the record removes and puts in its alternate bases, where it has any. The stop codon is lost unless that sequence shows it kept, and nothing is read of bases that no record gives.

It asks three things in turn:

1. A deletion whose removed bases repeat the bases beside them can be written at more than one place, and each place gives the same edited genome. Can the deletion be written at a place where it takes no coding base and no stop base? Then the codon is kept, and the row keeps `stop_retained_variant`.
2. Does an exon base stand at or beside each end of the edit, at one of the places it can be written at? Then the transcript's sequence across the edit is given, and vep-rs reads it at the old stop position, the positions of the stop codon.
3. Otherwise the edit removes a splice site or the transcript's end together with exon bases. Only the exon bases before the edit are given. Where an exon base stands directly before the edit, the alternate bases are given after it as well.

A given sequence that ends before the old stop position, or holds no stop codon there, gives `stop_lost`. One that holds a stop codon there keeps `stop_retained_variant` where every residue before it, each amino acid of the protein, is the one the unedited transcript gives. One that holds a stop codon there behind a changed residue gets neither stop term. Its row carries `coding_sequence_variant` instead, or no term in the stop term's place beside `start_lost`.

vep-rs writes no `frameshift_variant` beside the `stop_lost` of these rows. Ensembl VEP's frameshift test asks for both CDS coordinates (lines 1541-1564, the condition at line 1551), and no row of the class has both. Three bounds are vep-rs's own. It looks for the other places a deletion can be written at up to 10,000 positions each way along the reference sequence. Without a reference file it reads the record's own place alone. A base of that sequence that is not one of the four ends the search.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the module the Ensembl VEP 116.3 image holds), all in `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`:

- Lines 1234-1288, `stop_lost`: the routine reads the peptide alleles where the allele has some. Otherwise line 1263 takes the answer of `_ins_del_stop_altered`.
- Lines 1375-1376, `_overlaps_stop_codon`: the routine reads the span's cDNA start and end and returns 0 unless both are set. `_ins_del_stop_altered` returns 0 on that answer at line 1444. Its line 1453 returns 0 unless the cDNA start, the cDNA end and the CDS start are all set.
- Line 1318, `stop_retained`: an allele that changes the length and overlaps the stop codon is `stop_retained_variant` where `_ins_del_stop_altered_cil` answers 0.
- Lines 1488-1539, `_ins_del_stop_altered_cil`: line 1507 reads the same three coordinates, and line 1508 returns 0 unless all are set. The edit of line 1521 and the read of the three bases at the old stop position, lines 1528-1535, come after that return.
- Line 1469: an edited sequence shorter than the coding sequence is answered as altered, with no base read.
- Lines 1269, 1275 and 1278, the structural variant part of `stop_lost`: a deletion is `stop_lost` where its genomic span overlaps the stop codon's three genomic positions.

**What correcting it would change.** A test that told "not made" from "made, and the stop codon kept" would take `stop_retained_variant` off these rows. That alone leaves `coding_sequence_variant` in its place on 5,080 rows on GRCh38 and 1,000 on GRCh37, and no term in its place beside `start_lost` on 45 and 6. A reading of the sequence the record gives, as vep-rs makes it, writes `stop_lost` on 5,068 and 1,002 of them. The splice, UTR, intron and NMD terms of the rows are the same either way.

**The class it produces.** `stop_codon_removed_without_stop_codon_test`: a `Transcript` row of an insertion or deletion on which Ensembl VEP writes `stop_retained_variant` and vep-rs `stop_lost`, every other term being equal. Both rows lack a cDNA position, whole or at one end, or a CDS start. The golden corpora record each such row under the class's name. The format comparator sets a row of the class aside in a comparison that names the class; 1.7 says which rows and fields.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data carries the shape, and the divergence needs no flag or plugin.

Counts are rows of the genome-wide ClinVar dataset, on GRCh38 first and on GRCh37 second. Ensembl VEP 116.3 writes `stop_retained_variant` with no test made on 5,125 and 1,006 rows whose record's span holds the whole stop codon. They are rows of 931 and 451 records on 3,444 and 645 transcripts: 4,532 and 883 are deletions, and 593 and 123 replace the removed bases by other bases. vep-rs writes `stop_lost` on 5,068 and 1,002 of them, the rows of this class. It keeps `stop_retained_variant` on 55 and 4 and writes neither stop term on 2 and 0.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion of the stop codon and the two intron bases before it, on ENST00000635895

**In short.** `19 13208045 1993671 ACTACT A` on GRCh38, ClinVar 1993671, deletes the five bases 19:13208046-13208050. On the reverse-strand transcript ENST00000635895 they are the last two bases of an intron and the stop codon, which opens the last exon. Ensembl VEP 116.3 writes `splice_acceptor_variant,stop_retained_variant`, and vep-rs writes `splice_acceptor_variant,stop_lost`. The span begins in the intron and has no cDNA coordinate there, so Ensembl VEP makes no stop codon test and reads the missing answer as a stop codon kept. The record removes the codon with the splice site before it, and no record gives what the transcript reads next: no stop codon is shown at the old stop position.

### 1.1 What a reader needs to know first

A coding transcript ends its coding sequence with a stop codon, three bases that end translation, where its record shows one there. For a change at that codon Ensembl VEP writes one of two terms. `stop_lost` says that no stop codon stands at the old stop position after the change, and `stop_retained_variant` says that one does. This page calls the positions of that codon the old stop position.

On most transcripts the old stop position is the last three positions of the coding sequence. On some the genome carries only one or two bases of the stop codon. The rest are adenines of the poly(A) tail, or bases that a RefSeq model holds and the genome lacks. There a span that starts before the genome's bases of the codon and holds them holds the whole stop codon. [`stop_codon_completed_by_polyadenylation`](stop-codon-completed-by-polyadenylation.md) says which transcripts the poly(A) tail completes and covers a span that starts inside the codon.

Where an allele has reference and alternate amino acids, Ensembl VEP reads the two terms from them. Where it has none, as on the row of this page, a stop codon test decides. The test edits the coding sequence and the 3' UTR as the record says and reads the three bases at the old stop position.

The edit needs coordinates. A cDNA coordinate is the position of a base in the spliced transcript, and a CDS coordinate is its position in the coding sequence. A base of an intron has neither. A span that begins or ends on such a base has no cDNA coordinate at that end, and the row prints `?` there, as in `?-7057`.

An intron ends at its acceptor site, its last two bases, which mark where the next exon starts. Both engines write `splice_acceptor_variant` for a change there.

### 1.2 The record

ClinVar 1993671, a deletion of five bases, against the Ensembl GRCh38 cache of release 116:

```
19	13208045	1993671	ACTACT	A	.	.	.
```

Both engines drop the shared first base, so the record deletes 19:13208046-13208050. Its reference allele gives the bases: `A` at 13208045 and `CTACT` at 13208046-13208050. ENST00000635895 is a reverse-strand transcript of 47 exons of the gene CACNA1A. Three of its parts matter here:

| Part | Genomic span | Bases |
| --- | --- | --- |
| Intron 46 | 19:13208049-13208755 | 707 |
| Exon 47, the last | 19:13206446-13208048 | 1,603 |
| Stop codon | 19:13208046-13208048 | 3 |

On this strand the transcript runs from high to low coordinates and carries the complement of the genome's bases. Intron 46 ends with 13208050 and 13208049, which the transcript carries as `A` and `G`: the acceptor site. Exon 47 opens with 13208048, 13208047 and 13208046, which the transcript carries as `T`, `A` and `G`: the stop codon. The 1,603 - 3 = 1,600 bases after it are the 3' UTR.

The coding sequence runs from cDNA 272 to cDNA 7057, 7057 - 272 + 1 = 6,786 bases. That is 2,261 codons and the stop codon: 3 × 2,261 + 3 = 6,786. The stop codon is coding bases 6784 to 6786, at cDNA 7055 to 7057. Exon 46 ends one base before it, at cDNA 7054 and coding base 6783.

In the transcript's direction the deleted span begins at 13208050, an intron base without a cDNA coordinate. It ends at 13208046, cDNA 7057 and coding base 6786. The span holds the acceptor site and the whole stop codon, and no other base. The deletion can be written at no other place. The committed GRCh38 cache holds the transcript, and both engines write the rows of 1.3 with no flag beyond their cache options.

### 1.3 What each engine writes

The row on ENST00000635895 in the default format. Ensembl VEP 116.3:

```
1993671	19:13208046-13208050	-	ENSG00000141837	ENST00000635895	Transcript	splice_acceptor_variant,stop_retained_variant	?-7057	?-6786	?-2262	-	-	-	IMPACT=HIGH;STRAND=-1
```

vep-rs:

```
1993671	19:13208046-13208050	-	ENSG00000141837	ENST00000635895	Transcript	splice_acceptor_variant,stop_lost	?-7057	?-6786	?-2262	-	-	-	IMPACT=HIGH;STRAND=-1
```

The rows differ in one column, the `Consequence`: Ensembl VEP's `stop_retained_variant` stands against vep-rs's `stop_lost`. Both rows carry `splice_acceptor_variant` and `IMPACT=HIGH`. Both print the span's missing start as `?`: `?-7057` for the cDNA position, `?-6786` for the CDS position and `?-2262` for the protein position, since 6786 / 3 = 2262.

The record writes 30 rows in each engine. Eight differ, each in the `Consequence` alone and each as this row does. They are the rows of eight transcripts whose stop codon lies at 13208046-13208048, directly behind an intron that ends at 13208049. The other 22 rows are equal in every column. On 11 of them the five bases lie inside a coding exon, the row prints a position at both ends, and both engines write `frameshift_variant`.

### 1.4 Why Ensembl VEP writes that

The allele has no amino acids, the row's `Amino_acids` column being `-`. `stop_lost` therefore asks `_ins_del_stop_altered` (line 1263).

That routine asks whether the span overlaps the stop codon on the cDNA (line 1444). `_overlaps_stop_codon` reads the span's cDNA start and end (line 1375). The span begins at 13208050, an intron base, so the cDNA start is not set, and line 1376 returns 0. `_ins_del_stop_altered` returns that 0, and line 1263 stores it: the stop codon is not lost. Past the overlap, line 1453 would have returned the same 0 for want of the cDNA start.

`stop_retained` runs next. The allele changes the length, and its genomic span, 13208046-13208050, overlaps the stop codon at 13208046-13208048. Line 1318 then asks `_ins_del_stop_altered_cil` and writes the term where the answer is 0.

That routine reads the cDNA start, the cDNA end and the CDS start (line 1507). The cDNA start is not set, and line 1508 returns 0. The edit of line 1521 and the read of the three bases at the old stop position, lines 1528-1535, are never reached.

The routine answers 0 in two cases. Either it built the edited sequence and found a stop codon at the old stop position, or it returned before it built anything. Line 1318 cannot tell the two apart and writes `stop_retained_variant` for both. On this row no sequence was built and no base was read.

### 1.5 Why vep-rs's row is the correct one

The record removes 13208046-13208050, five bases. Three of them, 13208046 to 13208048, are the stop codon, `TAG` on the transcript, at coding bases 6784 to 6786. The other two, 13208049 and 13208050, are the acceptor site of intron 46. After the deletion the intron has 707 - 2 = 705 bases and exon 47 has 1,603 - 3 = 1,600, every one a base of the 3' UTR. No base of the stop codon is left, and the splice site that joined exon 46 to it is gone.

vep-rs asks what sequence the record gives. The deletion can be written at no other place, so it cannot be moved off the codon. In the transcript's direction the span's first base, 13208050, is an intron base, and so is 13208051 before it. No exon base stands at or beside that end, so the transcript's sequence across the edit is not given. What the transcript reads after exon 46 depends on a splice that no record shows.

What is given is the exon bases before the edit. They end with exon 46, at coding base 6783, the last base of codon 2,261: 3 × 2,261 = 6783. The old stop position is coding bases 6784 to 6786. The given sequence ends before it, so no stop codon is shown there, and the codon that stood there is deleted. vep-rs writes `stop_lost`.

`stop_retained_variant` would need a stop codon at the old stop position. Ensembl VEP's row asserts one and its test read no base. The module's own handling of a structural deletion agrees with vep-rs: it calls `stop_lost` from the genomic overlap alone (lines 1269, 1275 and 1278), and the stop codon at 13208046-13208048 lies inside 13208046-13208050.

Both rows keep `splice_acceptor_variant`, which the two intron bases give, and `HIGH`, the impact of that term and of `stop_lost` alike.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38/json_cache \
  --species homo_sapiens --assembly GRCh38
grep ENST00000635895 out.txt
```

The committed GRCh38 cache holds the transcript. vep-rs writes `splice_acceptor_variant,stop_lost` on the row. Ensembl VEP 116.3 takes the record with the cache options of [Reproducing these examples](reproducing.md) and its GRCh38 cache, and writes `splice_acceptor_variant,stop_retained_variant`. The command passes no reference file, so vep-rs reads the deletion at the record's own place alone, the one place this deletion has.

The entry's reproduction, `scripts/divergences/reproductions/stop_codon_removed_without_stop_codon_test.toml`, runs the record on both engines and holds the row to both recorded values.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) compares one Ensembl VEP output with one vep-rs output, every field of every row. A comparison case is one such comparison: an input, the flags and the classes it applies. Where a comparison case names the class `stop_codon_removed_without_stop_codon_test`, the comparator sets a row of the class aside and counts it once under that name. The rows of the class are the ones its predicate, `stop_codon_removed_without_stop_codon_test_pair`, accepts. It accepts a `Transcript` row of an insertion or deletion whose two term sets are equal but for Ensembl VEP's `stop_retained_variant` against vep-rs's `stop_lost`. Both rows must lack a cDNA position, whole or at one end, or a CDS start.

A comparison that is given the cache vep-rs annotated with, the candidate's cache, reads the row's transcript there too. The predicate then accepts a row only where that record shows a stop codon at the coding end, a whole one or one that polyadenylation completes. It accepts none where no stop codon stands or where the cache gives no record for the row. A comparison that is given no cache reads the two rows alone.

Of Ensembl VEP 116.3's rows against vep-rs's for the genome-wide ClinVar dataset, the class sets aside 5,067 rows on GRCh38 and 1,002 on GRCh37.

On a row the class sets aside, the comparator also sets aside the fields that follow the two terms, by the rule `term_pair_following_cell`. The rule takes a field when its two values are the ones the two terms give, so that the terms explain the field's difference. The fields are counted in outputs of the same dataset that both engines write with the plugins verified on real data loaded and the HGVS notations on. The class has 5,068 and 1,002 rows there: one more on GRCh38, a row that also differs in a splice term and that the comparison above does not count under this class. On them the rule sets aside one field:

- `IMPACT` is each engine's most severe term's, `HIGH` for `stop_lost` against `LOW`: it is set aside on 182 and 128 rows. `IMPACT` differs on no other row of the class. Of the class's rows, 44 and 6 carry `start_lost`, a `HIGH` term, in both engines' rows.

The rule takes no field that the two terms do not decide. Every other field of a row the class sets aside is compared, and a difference in it is reported. `HGVSc` is one field that the two terms do not decide. In the same outputs it differs on 548 and 153 rows of the class, and the comparator reports each of those differences.

The golden corpora record each row of the class under its name. The corpus builder (`tools/vep-golden-corpus`) writes the name as the row's `expected_divergence` in the corpus manifest. The manifests of 19 committed corpora hold 593 such entries: 490 in the lists of the corpora's own runs and 103 in the lists that their flag sets carry.

Rows beside the shape are no part of this entry. A row with a cDNA position at both ends and a CDS start, a row with another differing term, a row with the two terms the other way round and a substitution's row stay differences. On a transcript whose stop codon stands, rows of the same two returns on which vep-rs writes neither `stop_lost` nor `stop_retained_variant` belong to [`stop_retained_without_stop_codon_test`](stop-retained-without-stop-codon-test.md). Their span does not hold the whole stop codon, or the sequence given shows a stop codon behind a changed residue. That entry's class sets aside 272 and 62 rows in all. The 36 and 8 on which vep-rs writes `stop_gained` are rows of [`stop_codon_settled_before_unknown_bases`](stop-codon-settled-before-unknown-bases.md).

On a transcript where no stop codon stands at the coding end, vep-rs writes neither stop term for a span of this shape, and the row belongs to [`stop_codon_overlap_without_stop_codon`](stop-codon-overlap-without-stop-codon.md).

The SNP/indel comparator (`tools/vep-outputs-compare`) compares tuples: a row's Location, Allele, Feature, Feature_type and set of `Consequence` terms. Two tuples, one from each engine, that differ in the set of terms alone are a pair. The comparator's 116.3 registry holds a pair rule under the class's name, with the same predicate, `stop_codon_removed_without_stop_codon_test_pair`. The rule sets aside both tuples of a pair whose sets differ by Ensembl VEP's `stop_retained_variant` against vep-rs's `stop_lost` alone.

A tuple carries no cDNA or CDS position, so the rule reads the rows of the two outputs that carry the pair's two tuples. The Allele, the record's alleles in the first column and the Location show whether the pair is an insertion's or a deletion's. `cDNA_position` and `CDS_position` show a row's cDNA positions and its CDS start. The pair is set aside when it is an insertion's or a deletion's and every such row lacks a cDNA position, whole or at one end, or a CDS start. Given the candidate's cache, the rule sets the pair aside only where the transcript's record shows a stop codon at the coding end.

The structural-variant comparator (`tools/vep-sv-compare`) applies it to the pairs of the structural-variant sets.
