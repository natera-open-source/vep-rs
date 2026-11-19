# `stop_retained_variant` on an insertion or deletion whose stop codon test is not made

**The defect in one sentence.** For an insertion or deletion, Ensembl VEP 116.2's test of the edited stop codon returns 0 before it tests anything when the allele holds a letter that names no base, or when the span has no cDNA start, no cDNA end or no CDS start. The caller reads that 0 as a stop codon kept, and the row is written `stop_retained_variant`. On the rows this page describes, vep-rs writes `coding_sequence_variant`, `frameshift_variant` or nothing in that term's place; a row on which it writes `stop_lost` belongs to `stop_codon_removed_without_stop_codon_test`.

**Why this is a defect.** The test answers 0 in two cases that mean different things. In one it built the edited sequence and found a stop codon at the old stop position. In the other it returned before it built anything. The caller cannot tell the two apart and writes the term for both, so in the second case the row states the result of a test that was not made. The routine's twin behind `stop_lost` has the same two returns. Its caller reads the 0 of a test that was not made as a stop codon that is not lost, so it does not write `stop_lost`.

**What vep-rs does instead.** vep-rs makes the test wherever the span has both cDNA coordinates and a CDS start, on the three bases the edited sequence holds at the old stop position. A letter that names no base matters only where it is one of those three bases, which then are no stop codon. It writes `stop_retained_variant` only where they are a stop codon and every codon before them reads the reference's residue. A span the test cannot be made on keeps nothing: the row carries what the other tests give, `coding_sequence_variant`, `frameshift_variant` or no coding term beside `start_lost`.

A span that holds all three bases of the stop codon without a coordinate, where vep-rs writes `stop_lost`, is the subject of [`stop_codon_removed_without_stop_codon_test`](stop-codon-removed-without-stop-codon-test.md).

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the modules the Ensembl VEP 116.2 image holds):

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, line 1318, `stop_retained`: the arm taken when the alternate peptide is undefined, empty or holds an `X` (line 1313). It sets the term for an allele that changes the length when `_overlaps_stop_codon_cil` holds and `_ins_del_stop_altered_cil` returns 0.
- The same file, lines 1418-1427, `_overlaps_stop_codon_cil`: the record's genomic span is overlapped with the last three genomic positions of the coding region.
- The same file, lines 1488-1539, `_ins_del_stop_altered_cil`: line 1495 sets the answer to 0. Line 1498 returns it unless the allele is unambiguous DNA. Lines 1507-1508 read the span's cDNA start, cDNA end and CDS start and return it unless all three are set.
- The same routine, lines 1521-1535, reached only after both returns: line 1521 writes the allele into the coding sequence and 3' UTR. Line 1524 answers "altered" where the edited sequence is shorter than the coding sequence, and lines 1528-1535 read the three bases at the old stop position.
- `Bio/EnsEMBL/Variation/VariationFeatureOverlapAllele.pm`, lines 293-302, `seq_is_unambiguous_dna`: the test of line 1498, against the pattern of line 74, which admits `A`, `C`, `G`, `T` and `-` alone.
- `VariationEffect.pm`, lines 1443 and 1453, `_ins_del_stop_altered`: the twin has the same two returns. Line 1263 stores its answer as `stop_lost`, so a test that was not made reads as a stop codon not lost.
- The same file, lines 1549 and 1635: `frameshift` and `coding_unknown`, the tests behind `frameshift_variant` and `coding_sequence_variant`, read `stop_retained`.

**What correcting it would change.** A routine that made the test wherever the record allows it, and reported where it could not, would write `stop_retained_variant` only on a row it had tested. A row it could not test would carry what the other tests give: `coding_sequence_variant` on most, and `frameshift_variant` on an insertion whose length is no multiple of three. `IMPACT` would follow the row's most severe term, and a plugin that runs beside `frameshift_variant` would write its fields.

**The class it produces.** `stop_retained_without_stop_codon_test`: a `Transcript` row of an insertion or deletion whose Consequence holds `stop_retained_variant` in Ensembl VEP's row and not in vep-rs's, with `coding_sequence_variant`, `frameshift_variant` or nothing in its place. Both rows show that the test is not made: an allele with a letter other than `A`, `C`, `G` and `T`, a cDNA position missing whole or at one end, or no CDS start. A row on which vep-rs writes `stop_lost` belongs to [`stop_codon_removed_without_stop_codon_test`](stop-codon-removed-without-stop-codon-test.md).

The golden corpora record each row of the class under its name. The format comparator sets a row of the class aside in a comparison that names the class; 1.7 says which rows and fields.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape the record alone shows is an insertion or deletion, and no flag is needed. That page is the index of the documented differences, and it counts each shape over 15,646,481 chromosome 21 records. Of them 2,154,724 have an allele of another length than the reference (an insertion, a deletion, or a replacement of the removed bases by another number of bases), 137,713 per million. The transcript decides whether the divergence then occurs.

The rows are counted on the genome-wide ClinVar dataset, GRCh38 first and GRCh37 second. Both engines annotate the whole file of each assembly, and every row on which Ensembl VEP 116.2 writes `stop_retained_variant` for an insertion or deletion is compared with vep-rs's row of the same record and transcript. The engines differ in this entry's way on 308 rows of GRCh38 and 70 of GRCh37. vep-rs writes `coding_sequence_variant` in the term's place on 217 and 50 of them and `frameshift_variant` on 90 and 20. On 1 row of GRCh38, which carries `start_lost`, it writes nothing in the term's place.

The rows are of three kinds. On 142 and 30 the allele holds a letter that names no base. On 164 and 40 the span lacks a cDNA coordinate or a CDS start and holds one or two of the stop codon's three bases.

On 2 rows of GRCh38 the span lacks a coordinate and holds all three. For such a span vep-rs reads the sequence the record gives: the transcript's own exon bases with the record's edit made in them. The page of `stop_codon_removed_without_stop_codon_test` describes that reading. On these 2 rows the sequence shows a stop codon at the old stop position behind a changed residue, and vep-rs writes neither stop term. The three kinds sum to the class: 142 + 164 + 2 = 308 and 30 + 40 + 0 = 70.

The rows of the first kind are 134 and 28 insertions and 8 and 2 replacements, records that delete bases and insert others. The replacements lack a cDNA coordinate as well, so the test cannot be made on them. Read on the three bases the edited sequence holds at the old stop position, 124 and 24 of the insertions hold three bases that form no stop codon. On 8 and 4 one of the three is a letter that names no base, so they are no stop codon either. The other 2, on GRCh38, hold a stop codon there; vep-rs also asks the codons before it and does not write the term.

These rows are part of the 5,431 and 1,076 on which Ensembl VEP 116.2 writes the term through one of the two returns. On 5,068 and 1,002 of the others vep-rs writes `stop_lost`, the rows of `stop_codon_removed_without_stop_codon_test`. On the last 55 and 4 both engines write `stop_retained_variant`: 308 + 5,068 + 55 = 5,431 and 70 + 1,002 + 4 = 1,076.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion of the stop codon's first base, four coding bases and an intron's last base, on ENST00000987744

**In short.** ClinVar record 1498654 deletes six bases on GRCh38, 6:57193930-57193935, on the reverse-strand transcript ENST00000987744. On the transcript they are the last base of an intron, the first four bases of the next exon and the first base of the stop codon. The span opens in the intron, so it has no cDNA start, and Ensembl VEP 116.2 returns from its stop codon test before it reads a base. It writes `splice_acceptor_variant,stop_retained_variant,NMD_transcript_variant`. vep-rs writes `splice_acceptor_variant,coding_sequence_variant,NMD_transcript_variant`, and every other column is equal. With the five coding bases gone, the old stop position reads `ATT`, which is no stop codon.

### 1.1 What a reader needs to know first

A coding sequence ends with a stop codon, three bases that end translation. A change that touches those bases can remove the stop codon or leave one in place. Ensembl VEP has a term for each: `stop_lost` where the stop codon is lost, and `stop_retained_variant` where one stands at the old stop position after the change.

Where a change has an alternate peptide without an `X`, the letter a peptide carries for a codon that gives no single residue, Ensembl VEP decides from the reference and alternate peptides. For an insertion or deletion whose alternate peptide is undefined, empty or holds an `X`, it asks two routines. `_overlaps_stop_codon_cil` asks whether the record's genomic span reaches the stop codon. `_ins_del_stop_altered_cil` writes the allele into the transcript's sequence and reads the three bases at the old stop position. This page calls the second routine the stop codon test.

A record's span is the stretch of reference bases it removes or replaces, and the test needs the span's place on the transcript. A transcript's bases are numbered along its exons, the cDNA positions, and its coding bases from the first coding base, the CDS positions. A base of an intron has neither. A row whose span opens or closes outside the exons prints `?` for that end, as in `?-1043`.

An intron ends at its acceptor site, whose last two bases mark where the next exon starts. An allele can also hold a letter other than `A`, `C`, `G` and `T`, one that names no single base. `coding_sequence_variant` is the general coding term: both engines write it for a change in the coding sequence that no more specific coding term describes.

### 1.2 The record

ClinVar record 1498654, against the Ensembl GRCh38 cache of release 116:

```
6	57193929	1498654	TAAAAAC	T	.	.	.
```

Both engines drop the shared first base, so the record deletes six bases, 6:57193930-57193935. ENST00000987744 is a reverse-strand transcript of RAB23 with six exons, and its biotype is `nonsense_mediated_decay`. Four of its parts matter here:

| Part | Genomic span | Bases |
| --- | --- | --- |
| Exon 4 | 6:57196450-57196606 | 157 |
| Intron 4 | 6:57193935-57196449 | 2,515 |
| Exon 5 | 6:57193842-57193934 | 93 |
| Stop codon | 6:57193928-57193930 | 3 |

On this strand the transcript runs from high to low coordinates, so intron 4 ends at 57193935 and exon 5 opens at 57193934. Exon 5 opens at cDNA 1039 and the coding sequence at cDNA 641, so the exon's first base is coding base 1039 - 641 + 1 = 399. The coding sequence has 405 bases, 134 codons and the stop codon: 3 × 134 + 3 = 405. The stop codon is coding bases 403 to 405, at 57193930, 57193929 and 57193928.

Residue n is translated from coding bases 3n - 2 to 3n. Codon 133 is bases 397 to 399 and reads `AAT`. Codon 134, the last before the stop codon, is bases 400 to 402 and reads `TTT`. The stop codon reads `TAA`. The six deleted bases, in the transcript's order:

| Genomic position | On the transcript | What it is |
| --- | --- | --- |
| 57193935 | `G` | last base of intron 4, in its acceptor site |
| 57193934 to 57193931 | `T`, `T`, `T`, `T` | coding bases 399 to 402 |
| 57193930 | `T` | coding base 403, first of the stop codon |

The deletion takes five coding bases, 399 to 403, and one intron base. Neither engine takes a flag beyond its cache for this record.

### 1.3 What each engine writes

The row on ENST00000987744 in the default format. Ensembl VEP 116.2:

```
1498654	6:57193930-57193935	-	ENSG00000112210	ENST00000987744	Transcript	splice_acceptor_variant,stop_retained_variant,NMD_transcript_variant	?-1043	?-403	?-135	-	-	-	IMPACT=HIGH;STRAND=-1
```

vep-rs:

```
1498654	6:57193930-57193935	-	ENSG00000112210	ENST00000987744	Transcript	splice_acceptor_variant,coding_sequence_variant,NMD_transcript_variant	?-1043	?-403	?-135	-	-	-	IMPACT=HIGH;STRAND=-1
```

The rows differ in one term of the Consequence column: `stop_retained_variant` in Ensembl VEP's row, `coding_sequence_variant` in vep-rs's. The position columns are equal, `?-1043`, `?-403` and `?-135`, and the `?` is the end of the span that lies in the intron. `IMPACT` is `HIGH` in both rows.

The record writes 27 rows in each engine, and the other 26 are equal in every column. Twenty of them are rows of protein-coding transcripts of the same gene whose exon opens at the same base, 57193934. There the same six bases are the acceptor site's last base and five coding bases with no stop codon near: bases 482 to 486 of 714 on 19 transcripts, and 446 to 450 of 678 on one. On all 20 both engines write `splice_acceptor_variant,coding_sequence_variant`.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP's row prints `-` for its codons and amino acids: the span has no start on the coding sequence, so there is no alternate peptide to compare. `stop_retained` therefore takes its second arm, line 1318. That arm holds when the allele changes the length, the span overlaps the stop codon, and the stop codon test returns 0.

The first two hold. The record deletes six bases. The last three genomic positions of the coding region are 57193928 to 57193930, and the span reaches 57193930, so `_overlaps_stop_codon_cil` holds (lines 1418-1427).

The stop codon test then runs. Line 1495 sets its answer to 0. Line 1498 passes: a deletion's allele is `-`, which the pattern admits. Lines 1507-1508 ask for the span's cDNA start, cDNA end and CDS start. The span opens at 57193935, an intron base, so it has no cDNA start and no CDS start: the row prints `?-1043` and `?-403`. The routine returns its 0 there. Lines 1521-1535 do not run, so no sequence is edited and no codon is read.

Back on line 1318 the 0 reads as "not altered", `stop_retained` is set, and the row carries `stop_retained_variant`. `coding_unknown` reads `stop_retained` (line 1635), so `coding_sequence_variant` is left off the row. The twin behind `stop_lost` returns 0 for the same span (line 1453), and line 1263 reads that as a stop codon not lost. One span that neither routine tested is "kept" to one caller and "not lost" to the other.

The other return works the same way. An allele that holds a letter outside the pattern stops at line 1498, before any base is read, and line 1318 writes the term.

### 1.5 Why vep-rs's row is the correct one

The span runs from 57193930 to 57193935: 57193935 - 57193930 + 1 = 6 bases. Exon 5 covers 57193842 to 57193934, so 57193935 lies outside it, in intron 4, and has no cDNA position. The other five bases are cDNA 1039 to 1043, since exon 5 opens at cDNA 1039 and 1039 + 4 = 1043. In coding positions they are 399 to 403: 1039 - 641 + 1 = 399 and 1043 - 641 + 1 = 403.

Coding base 399 is the last base of codon 133, since 3 × 133 = 399. Bases 400 to 402 are codon 134, the last before the stop codon. Base 403 is the first base of the stop codon, since 3 × 135 - 2 = 403. The record removes one base of codon 133, the whole of codon 134 and one base of the stop codon.

`stop_retained_variant` says that a stop codon stands at the old stop position after the change. Read with the five coding bases gone and nothing else changed, the sequence closes up. Bases 397 and 398, `A` and `A`, are followed by bases 404 and 405, `A` and `A`, and then by the 3' UTR, which opens `G`, `T`, `A`, `T`, `T`. Codon 133 reads `AAA` where it read `AAT`, and codon 134 reads `AGT` where it read `TTT`. The old stop position, bases 403 to 405, holds `ATT`, the third to fifth bases of the 3' UTR. `ATT` is no stop codon, so on that reading two residues change and translation runs on.

That reading takes the transcript to be spliced as before. The record also removes the last base of the acceptor site, which both engines report as `splice_acceptor_variant`, and it does not show how the exon is joined without that base. What the record does show is that a base of the stop codon and the codon before it are gone. Nothing in it shows a stop codon kept.

Ensembl VEP did not find otherwise, since its test returned on line 1508 without reading a base. On this record its own rows agree with vep-rs's wherever the span reaches no stop codon. On the 20 protein-coding transcripts of 1.3 the same deletion takes the same acceptor base and five coding bases, and Ensembl VEP writes `splice_acceptor_variant,coding_sequence_variant`. vep-rs writes those two terms on all 21 transcripts, beside `NMD_transcript_variant` on this one.

vep-rs cannot make the test on this span either, for want of the cDNA start. It reads that as keeping nothing, so the general coding term stands in the row. It does not write `stop_lost`: for a span without its coordinates it writes that term only where the span holds all three bases of the stop codon, and this span holds one.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38/json_cache \
  --species homo_sapiens --assembly GRCh38
grep ENST00000987744 out.txt
```

The committed GRCh38 cache holds the transcript. vep-rs writes `splice_acceptor_variant,coding_sequence_variant,NMD_transcript_variant` on the row. Ensembl VEP 116.2 takes the same record with the cache options of [Reproducing these examples](reproducing.md) and its GRCh38 cache, and writes `splice_acceptor_variant,stop_retained_variant,NMD_transcript_variant`.

The entry's reproduction, `scripts/divergences/reproductions/stop_retained_without_stop_codon_test.toml`, runs the record on both engines. It holds the row of ENST00000987744 to its recorded value in seven columns: the location, the transcript, the Consequence, the three position columns and `IMPACT`.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) compares one Ensembl VEP output with one vep-rs output, every field of every row. A comparison case is one such comparison: an input, the flags and the classes it applies. Where a comparison case names the class `stop_retained_without_stop_codon_test`, the comparator sets a row of the class aside and counts it once under that name. It recognises a row of the class by the predicate `stop_retained_without_stop_codon_test_pair`.

The golden corpora are the committed inputs under `tests/golden/116/`, each with Ensembl VEP 116.2's outputs and a manifest beside it. The corpus builder classes a corpus row of the shape under the entry's name. The manifests of 12 corpora name 161 rows under `stop_retained_without_stop_codon_test`, each with the consequence set vep-rs writes. They are rows of those corpora, counted apart from the 308 and 70 rows of the ClinVar dataset. The golden harness checks vep-rs's row of each against that recorded set.

The class takes the rows of the shape that "The class it produces" gives above. A row with both cDNA positions, a CDS start and an allele of bases is not a row of the class, and neither is a row where vep-rs keeps `stop_retained_variant` or writes a term other than `coding_sequence_variant` or `frameshift_variant` in its place, such as `stop_lost`. The comparator compares such a row like any other and reports its differences.

Of Ensembl VEP 116.2's rows against vep-rs's on the genome-wide ClinVar dataset, the class sets aside the 308 rows of GRCh38 and the 70 of GRCh37. On a row the class sets aside, the comparator also sets aside the fields that follow the two sets of terms, by the rule `term_pair_following_cell`. The rule takes a field when its two values are the ones the two engines' sets of terms give, so that the terms explain the field's difference. The fields it takes on those rows are listed below.

The fields are counted in the runs that measure the plugins verified on real data on the dataset. In them both engines annotate it with `--hgvs` and those plugins switched on, Downstream and NMD among them. LOFTEE is run on GRCh38 alone. Each count gives GRCh38 first and GRCh37 second.

- `IMPACT` is set aside on 142 and 34 rows. On these rows each engine prints the impact of its own row's most severe term, and the two impacts differ. On the class's other rows the two impacts are equal.
- Downstream's `DownstreamProtein` is set aside on 90 and 20 rows. It stands in vep-rs's row alone, where the row carries `frameshift_variant`.
- Downstream's `ProteinLengthChange` is set aside on the same rows, where it stands in vep-rs's row alone: 90 and 20.
- `NMD` is set aside on 76 and 19 rows. It stands in vep-rs's row alone where `frameshift_variant` lets the plugin run.
- LOFTEE's `LoF` is set aside on 76 rows of GRCh38. On those rows `frameshift_variant` lets LOFTEE run on vep-rs's row, and the field stands in that row alone.
- LOFTEE's `LoF_info` is set aside on the same rows, where it stands in vep-rs's row alone: 76.
- LOFTEE's `LoF_filter` is set aside on 1 of them, where it stands in vep-rs's row alone.
- `HGVSp` is set aside on 2 rows of GRCh38 and 0 of GRCh37. On both rows Ensembl VEP writes a notation of no change. vep-rs writes an extension's notation, the form a frameshift's notation takes where the frameshift begins at the stop codon. On both rows the record deletes an intron's last base and the stop codon's first base. The notation writes the same deletion one base further along the transcript, on the stop codon's first two bases, where it lies inside one exon. Where a deletion as the notation writes it runs across exons, vep-rs writes no frameshift's notation, and its `HGVSp` is the one Ensembl VEP writes.

The rule takes nothing else: every other field of the row is compared, and a difference in it is reported. On rows of the class that holds for `NMD` where Ensembl VEP writes it and vep-rs does not: 33 and 5 rows whose two term sets both hold a splice site term. On 1 row of GRCh38 it holds for `HGVSc`, `HGVS_OFFSET`, `HGVS_IntronStartOffset` and `HGVS_IntronEndOffset`, which vep-rs alone writes. `LoF_flags` differs on no row.

The SNP/indel comparator (`scripts/concordance/compare_vep_outputs.py`) compares tuples: a row's Location, Allele, Feature, Feature_type and set of `Consequence` terms. Two tuples, one from each engine, that differ in the set of terms alone are a pair. The comparator's 116.2 registry holds a pair rule under the class's name, with the same predicate, `stop_retained_without_stop_codon_test_pair`. The rule sets aside both tuples of a pair in which Ensembl VEP's set holds `stop_retained_variant` and vep-rs's is that set without it, with nothing, `coding_sequence_variant` or `frameshift_variant` in its place.

A tuple carries no cDNA or CDS position, so the rule reads the rows of the two outputs that carry the pair's two tuples. The Allele, the record's alleles in the first column and the Location show whether the pair is an insertion's or a deletion's. `cDNA_position` and `CDS_position` show a row's cDNA positions and its CDS start. The pair is set aside when it is an insertion's or a deletion's and every such row shows that the test is not made. A row shows that by an allele with a letter other than `A`, `C`, `G` and `T`, a cDNA position missing whole or at one end, or no CDS start.

`tools/vep-outputs-compare`, the Rust build of this comparator, applies the same rule to the same rows. The structural-variant comparator (`scripts/validation/compare_sv_concordance.py`) applies it to the pairs of the structural-variant sets.
