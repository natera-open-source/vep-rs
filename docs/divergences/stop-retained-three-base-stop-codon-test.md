# A stop codon written as retained although the protein before it changes

**The defect in one sentence.** For an insertion or deletion at a stop codon, Ensembl VEP 116.2 writes the allele into the coding sequence and the 3' UTR after it and reads three bases alone: the ones that then stand where the stop codon stood, the old stop position. Where they spell a stop codon, as when a deletion brings bases of the 3' UTR there, the row is written `stop_retained_variant`, although a residue before it changes or an earlier stop codon ends the protein. On such a row vep-rs writes `coding_sequence_variant`, or `frameshift_variant`, with `stop_gained` where the changed peptide gains a stop.

**Why this is a defect.** `stop_retained_variant` says that the stop codon stays, and Ensembl VEP keeps `frameshift_variant`, `stop_gained` and `coding_sequence_variant` off a row that carries it. On these rows the edited sequence does not give the reference protein. All but 2 are rows of deletions and of replacements, records that delete bases and insert another number in their place: 101 rows of the genome-wide ClinVar dataset on GRCh38 and 27 on GRCh37. On 99 of the 101 and on all 27, the record's span holds the whole stop codon and ends in the 3' UTR. There the stop codon the row calls retained is deleted, and the one the test finds is made of 3' UTR bases.

**What vep-rs does instead.** vep-rs makes the same test and asks one thing more: whether every codon before the old stop position reads the reference's residue in the edited sequence. It writes `stop_retained_variant` only where both hold. Elsewhere the row carries the terms that Ensembl VEP's other tests give once the stop codon is not read as retained: `coding_sequence_variant`, or `frameshift_variant`, with `stop_gained` where the changed peptide gains a stop.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the module the Ensembl VEP 116.2 image holds):

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1313-1318, `stop_retained`: where the alternate peptide is undefined, empty or holds an `X`, the term is written for an allele that changes the length, overlaps the stop codon and passes `!_ins_del_stop_altered_cil(@_)`.
- The same file, lines 1488-1539, `_ins_del_stop_altered_cil`: line 1508 makes the test only for a span with both cDNA coordinates and a CDS start. Line 1521 writes the allele into the coding sequence followed by the 3' UTR.
- The same routine, line 1524: an edited sequence shorter than the coding sequence is altered. Lines 1528-1534 translate the three bases at the old stop position, and line 1535 answers "altered" unless they are a stop codon. Nothing asks whether the bases before that position are the reference's.
- The same file, lines 1549, 1221 and 1635: `frameshift`, `stop_gained` and `coding_unknown` each read `stop_retained`, so their terms leave a row whose stop codon is read as retained.
- The same file, lines 1408-1415, in `_overlaps_stop_codon_cil`: an insertion's start is its end plus one, so its span holds no base. The routine moves the genomic start by the inserted length: it adds the length on the forward strand (line 1410) and subtracts it on the reverse strand (line 1413). On the forward strand the start moves further past the end, the case that `stop_codon_overlap_insertion_range` describes.
- The same lines, on the reverse strand: the move gives a span as long as the insertion, which lies after the insertion site in the transcript's direction. A reverse-strand insertion inside the stop codon, or before it in the transcript's direction by fewer bases than it inserts, therefore overlaps the codon and reaches the test. The 2 insertion rows counted under **How often the shape occurs** are reached in this way.
- The same file, lines 1065 and 1071, `_ins_del_start_altered`, for contrast: the start codon's test compares everything the edit could have changed. Its edited sequence is the 5' UTR followed by the coding sequence, with the allele written in. Line 1065 answers "not altered" early only where the 5' UTR and the start codon are unchanged. Otherwise line 1071 answers "altered" where the coding sequence differs from the last bases of the edited sequence, as many as the coding sequence has.

**What correcting it would change.** A test that also asked whether the bases before the old stop position are the reference's would answer "altered" on these rows. `stop_retained_variant` would leave them, and the row's other terms would follow from the module's own tests as they stand. That is `coding_sequence_variant` (lines 1616-1638) where the record runs out of the coding sequence, and `frameshift_variant` where it has both CDS coordinates. `stop_lost` would stay unwritten: lines 1473-1480 ask whether a stop codon stands at the old stop position, and one does. `IMPACT` and the plugin fields that are read from the terms would follow; 1.7 lists them.

**The class it produces.** `stop_retained_three_base_stop_codon_test`: a `Transcript` row of an insertion or deletion on which the test is made. Its columns show that: `cDNA_position` gives a start and an end, `CDS_position` a start, and the allele is `-` or letters of the four bases alone. Its two consequence sets differ in one of two ways. Ensembl VEP's holds `stop_retained_variant` where vep-rs's holds `coding_sequence_variant`. Or Ensembl VEP's holds `stop_retained_variant`, alone or beside an in-frame or protein-altering term, where vep-rs's holds `frameshift_variant` with or without `stop_gained`. The golden corpora record each such row under the class's name. The format comparator sets a row of the class aside in a comparison that names the class; 1.7 says which rows and fields.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data carries the shape, and the divergence needs no flag or plugin. The entry's row on that page, the index of the intended divergences, gives a rate, 137,713 per million. It is the share of records with an allele of another length than the reference (an insertion, a deletion, or a replacement of the removed bases by another number of bases) among the chromosome 21 records the index counts shapes on. It is an upper bound, since the transcript decides whether the divergence occurs.

In the genome-wide ClinVar dataset the two engines' rows differ in this way on 103 rows on GRCh38 and on 27 on GRCh37. On GRCh38, 101 rows belong to 20 records, deletions on 78 rows and replacements on 23, and 2 rows belong to two insertions on reverse-strand transcripts. On GRCh37 the 27 rows belong to 12 records, deletions on 23 rows and replacements on 4. On the 101 and 27 rows, a residue changes in a protein of the same length on 77 and 24, and the protein is shorter on 24 and 3.

The rows by the terms that differ:

| Ensembl VEP 116.2 | vep-rs | GRCh38 | GRCh37 |
| --- | --- | --- | --- |
| `stop_retained_variant` | `coding_sequence_variant` | 99 | 27 |
| `stop_retained_variant` | `frameshift_variant` | 1 | 0 |
| `inframe_deletion,stop_retained_variant` | `frameshift_variant` | 1 | 0 |
| `inframe_insertion,stop_retained_variant` | `stop_gained,frameshift_variant` | 1 | 0 |
| `protein_altering_variant,stop_retained_variant` | `stop_gained,frameshift_variant` | 1 | 0 |

The committed golden corpora hold 37 rows of the class, in seven corpora. The worked record below is one of them.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion of the last sense codon and the stop codon's first base, on ENST00001021815

**In short.** `21 26841549 rs761640407 CATTT C` on GRCh38 deletes the four bases 21:26841550-26841553 from the reverse-strand transcript ENST00001021815. They are the last sense codon, c.1228-1230, and the first base of the stop codon, c.1231. In the edited sequence residue 410 reads glutamate where the reference reads lysine, and the three bases after it are 3' UTR bases that spell a stop codon. Ensembl VEP 116.2 reads those three bases alone and writes `inframe_deletion,stop_retained_variant,NMD_transcript_variant` with `IMPACT=MODERATE`. vep-rs writes `frameshift_variant,NMD_transcript_variant` with `IMPACT=HIGH`. Every other column of the row is equal.

### 1.1 What a reader needs to know first

A transcript's coding sequence ends with its stop codon, the three bases at which translation ends. The bases after it are the 3' UTR. A coding base is counted from the first base of the coding sequence, and residue n is translated from coding bases 3n - 2 to 3n. This page calls the last three coding positions of the reference the old stop position: after an insertion or deletion, other bases can stand there.

`stop_retained_variant` says that the stop codon stays: the variant touches the codon, and the edited sequence gives the reference protein residue for residue.

Where Ensembl VEP cannot read a changed peptide for an insertion or deletion, it decides the term by a test on the sequence. This page calls it the stop codon test. The test writes the allele into the coding sequence followed by the 3' UTR and looks at the old stop position. A stop codon there is read as the stop codon retained. The answer reaches beyond the one term: with it, Ensembl VEP writes no `frameshift_variant`, `stop_gained` or `coding_sequence_variant` on the row.

Both engines write `IMPACT` from a row's most severe term: `MODERATE` for `inframe_deletion`, `HIGH` for `frameshift_variant`.

On a reverse-strand transcript the coding sequence runs from high to low genomic coordinates, and its bases are the complements of the genome's.

### 1.2 The record

A record of the committed GRCh38 golden corpus, `tests/golden/116/GRCh38/variants.vcf`, against the Ensembl GRCh38 cache of release 116:

```
21	26841549	rs761640407	CATTT	C	.	.	.
```

Both engines drop the shared first base, so the record deletes the four bases 21:26841550-26841553. ENST00001021815.1 is a reverse-strand transcript of ADAMTS1 with nine exons. Its fourth exon, 21:26840998-26841570, holds the end of the coding sequence:

| Part | Genomic span | On the transcript |
| --- | --- | --- |
| Codon 410, the last sense codon | 21:26841551-26841553 | c.1228-1230, `AAA` |
| Stop codon | 21:26841548-26841550 | c.1231-1233, `TGA` |
| First four bases of the 3' UTR | 21:26841544-26841547 | `ATGA` |
| The deletion | 21:26841550-26841553 | c.1228-1231, `AAAT` |

The coding sequence runs from cDNA 456 to cDNA 1688: 1688 - 456 + 1 = 1,233 bases, 410 codons and the stop codon, 3 × 410 + 3 = 1,233. The stop codon is cDNA 1686-1688, and 1686 - 456 + 1 = 1231. Codon 410 is coding bases 3 × 410 - 2 = 1228 to 1230 and reads lysine. The 3' UTR follows to the transcript's last base, cDNA 5056: 5,056 - 1,688 = 3,368 bases.

The record's reference allele agrees with the table. Read on the transcript's strand, from 26841553 down to 26841549, `CATTT` is `AAATG`: codon 410 and the first two bases of the stop codon. The deleted bases are cDNA 1683-1686, coding bases 1683 - 456 + 1 = 1228 to 1686 - 456 + 1 = 1231.

The record needs no flag beyond the cache.

### 1.3 What each engine writes

The row on ENST00001021815 in the default format. Ensembl VEP 116.2:

```
rs761640407	21:26841550-26841553	-	ENSG00000154734	ENST00001021815	Transcript	inframe_deletion,stop_retained_variant,NMD_transcript_variant	1683-1686	1228-1231	410-411	K*/X	AAATga/ga	-	IMPACT=MODERATE;STRAND=-1
```

vep-rs:

```
rs761640407	21:26841550-26841553	-	ENSG00000154734	ENST00001021815	Transcript	frameshift_variant,NMD_transcript_variant	1683-1686	1228-1231	410-411	K*/X	AAATga/ga	-	IMPACT=HIGH;STRAND=-1
```

The rows differ in the `Consequence` column and in the `IMPACT` of the last column. Both engines place the deletion alike: cDNA 1683-1686, CDS 1228-1231 and protein 410-411. Both write `K*/X` for the amino acids, a lysine and the stop against a changed peptide that holds an `X`, and `AAATga/ga` for the codons. Both write `NMD_transcript_variant`.

The record writes 18 rows in each engine, and the other 17 are equal in every column.

### 1.4 Why Ensembl VEP writes that

The alternate peptide of the row holds an `X`, which line 1313 asks about, so `stop_retained` takes the arm of line 1318. The arm asks three things: that the allele changes the length, that it overlaps the stop codon, and that the stop codon test does not answer "altered".

The overlap holds. `_overlaps_stop_codon_cil` takes the last three genomic positions of the coding region as the stop codon (lines 1418-1422), here 21:26841548-26841550. The deletion's span, 26841550-26841553, holds 26841550.

The test is made. The row has both cDNA coordinates, 1683 and 1686, and a CDS start, 1228, which is what line 1508 asks for.

Line 1521 writes the allele into the coding sequence followed by the 3' UTR, 1,233 + 3,368 = 4,601 bases. The edit starts at the CDS start, base 1228, and replaces 1686 - 1683 + 1 = 4 bases with nothing. The edited sequence has 4,601 - 4 = 4,597 bases. It is not shorter than the coding sequence, so line 1524 does not answer "altered".

Lines 1528-1534 take three bases at the offset 1,233 - 3 = 1230, which are bases 1231 to 1233 of the edited sequence. Four bases were removed before them, so they are the reference's bases 1235 to 1237. Those are the second, third and fourth base of the 3' UTR, `TGA`. They translate to a stop, and line 1535 answers "not altered". Nothing in the routine asks what stands before those three bases.

Line 1318 then sets `stop_retained`, and the row carries `stop_retained_variant`. `frameshift` returns 0 at line 1549 for such an allele, so `frameshift_variant` is not written. `_get_codon_alleles` returns nothing while `frameshift` holds (line 815). Here it does not hold, so the two codon alleles are formed, the reference's and the changed one: the two sides of the codons in 1.3. The in-frame deletion test (lines 1173 and 1178) asks whether the changed allele is shorter than the reference's and is its beginning or its end. Here it is its end. The row carries `inframe_deletion` for four deleted bases, and that term's impact, `MODERATE`.

### 1.5 Why vep-rs's row is the correct one

The deletion removes coding bases 1228 to 1231: 1231 - 1228 + 1 = 4 bases. In the edited sequence, coding base 1227 is followed by the stop codon's second and third base, `G` and `A`, and then by the 3' UTR, which opens `A`, `T`, `G`, `A`.

Codon 410 therefore reads `GAA`, glutamate, where the reference reads `AAA`, lysine. The three bases after it, `TGA`, are a stop codon. They are the second, third and fourth base of the 3' UTR.

The edited protein has 410 residues, as the reference protein has, and its last residue is another amino acid. The transcript's own stop codon does not stay. Its first base is deleted, and its other two bases open codon 410. The stop codon that ends the edited protein exists only because of the four bases the 3' UTR opens with.

Four is no multiple of three: 4 = 3 + 1. Every base after the deletion is read one frame off, and the shifted frame happens to end after one codon. A row that reads `inframe_deletion,stop_retained_variant` describes a deletion that keeps the reading frame and the stop codon. This record changes residue 410 and removes the stop codon's first base.

vep-rs makes the test of 1.4 and asks the further question. The three bases at the old stop position are a stop codon, but codon 410, a codon before them, does not read the reference's residue. `stop_retained_variant` is not written. The row's terms then follow from the tests Ensembl VEP itself applies to a row without that term. The record has both CDS coordinates, and four deleted bases are no multiple of three, so the row carries `frameshift_variant`. Its impact is that term's, `HIGH`.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38/json_cache \
  --species homo_sapiens --assembly GRCh38
grep ENST00001021815 out.txt
```

The committed GRCh38 cache holds the transcript. vep-rs writes `frameshift_variant,NMD_transcript_variant` and `IMPACT=HIGH`. Ensembl VEP 116.2 takes the record with the cache options of [Reproducing these examples](reproducing.md) and its GRCh38 cache, with no further flag. It writes `inframe_deletion,stop_retained_variant,NMD_transcript_variant` and `IMPACT=MODERATE`.

The entry's reproduction, `scripts/divergences/reproductions/stop_retained_three_base_stop_codon_test.toml`, runs the record on both engines. It holds the row's location, terms, three positions and impact to their recorded values.

### 1.7 How the comparator treats this class

The golden corpora record the class. The golden corpus builder (`scripts/golden/build_golden_corpus.py`) classes a corpus row of the shape under the entry's name. The manifests of seven golden corpora list 37 such rows in their `divergences`, each with the terms vep-rs writes. The golden harness and the format comparator check vep-rs's row against a listed row in place of Ensembl VEP's, as [How to read the counts](../intended-divergences.md#how-to-read-the-counts) describes.

The format comparator (`scripts/concordance/compare_vep_formats.py`) compares one Ensembl VEP output with one vep-rs output, every field of every row. A comparison case is one such comparison: an input, the flags and the classes it applies. Where a comparison case names the class `stop_retained_three_base_stop_codon_test`, the comparator sets a row of the class aside and counts it once under that name. The class's predicate, `stop_retained_three_base_stop_codon_test_pair`, recognises the row's shape. Of the columns that differ, the predicate reads the `Consequence` alone. Of Ensembl VEP 116.2's rows against vep-rs's for the genome-wide ClinVar dataset, the class sets aside 103 rows on GRCh38 and 27 on GRCh37.

On a row the class sets aside, the comparator also sets aside the fields that follow the two sets of terms, by the rule `term_pair_following_cell`. The rule takes a field when its two values are the ones the two sets of terms give, so that the terms explain the field's difference.

The counts below are the rows on which the field is set aside, GRCh38 first and GRCh37 second. They are counted in outputs that both engines write with the plugins verified on real data loaded and the HGVS notations on. LOFTEE is run on GRCh38 alone.

The fields the rule names:

- `IMPACT` is the impact of a row's most severe term.
- [Downstream](../plugins/transcript-model.md#downstream) is a plugin that writes on a row with `frameshift_variant` and no splice-site term. `DownstreamProtein` is the changed translation up to its first stop, and `ProteinLengthChange` is that translation's length minus the reference protein's.
- `NMD` is the field of the NMD plugin, not the term `NMD_transcript_variant`. The plugin writes `NMD_escaping_variant` on stop-gained, frameshift and splice-site rows that escape nonsense-mediated decay.
- LOFTEE is the plugin that the [support table](../plugins.md#support-table) of the plugins page lists as [LoFTEE](../plugins/loftee.md). It annotates `stop_gained`, `frameshift_variant`, `splice_donor_variant` and `splice_acceptor_variant`, and its fields are `LoF`, `LoF_filter`, `LoF_flags` and `LoF_info`. The rule names all four.
- `HGVSp` is the protein notation that `--hgvs` adds.

A plugin field stands in vep-rs's row alone where vep-rs's terms let the plugin write and Ensembl VEP's do not.

The fields set aside on the 103 and 27 rows of `stop_retained_three_base_stop_codon_test`:

- `IMPACT` is set aside on 92 and 21 rows: the rows on which the two engines' most severe terms have different impacts.
- Downstream's `DownstreamProtein` is set aside on 4 rows of GRCh38, where it stands in vep-rs's row alone and that row carries `frameshift_variant`.
- Downstream's `ProteinLengthChange` is set aside on the same 4 rows, where it stands in vep-rs's row alone.
- `NMD` is set aside on 2 rows of GRCh38, where it stands in vep-rs's row alone.
- LOFTEE's `LoF` and `LoF_info` are set aside on 2 rows and `LoF_filter` on 1, where they stand in vep-rs's row alone. `LoF_flags` is set aside on no row.
- `HGVSp` is set aside on 4 rows of GRCh38, where it is a frameshift's notation in vep-rs's row and another notation in Ensembl VEP's: a deletion-insertion's on 3 and an insertion's on 1.

On the 92 and 21 rows the impacts are `LOW` in Ensembl VEP's row against `MODIFIER` in vep-rs's on 88 and 21. On GRCh38 they are also `LOW` against `HIGH` on 1 row and `MODERATE` against `HIGH` on 3. On the other 11 and 6 rows of the class the two impacts are equal. The rule names no other field. Every other field of the row is compared, and a difference in one is reported.

Rows beside the shape are no part of this entry. A row whose stop codon test is not made is the shape of [`stop_retained_without_stop_codon_test`](stop-retained-without-stop-codon-test.md). A row whose span holds the whole stop codon and lacks a cDNA position or a CDS start is the shape of [`stop_codon_removed_without_stop_codon_test`](stop-codon-removed-without-stop-codon-test.md). A substitution's row, a row on which vep-rs keeps the term, and a row with any other differing term are reported as differences.

The SNP/indel comparator (`scripts/concordance/compare_vep_outputs.py`) compares tuples: a row's Location, Allele, Feature, Feature_type and set of `Consequence` terms. Two tuples, one from each engine, that differ in the set of terms alone are a pair. The comparator's 116.2 registry holds a pair rule under the class's name, with the same predicate, `stop_retained_three_base_stop_codon_test_pair`. The rule sets aside both tuples of a pair whose sets differ in one of two ways. Ensembl VEP's holds `stop_retained_variant` where vep-rs's holds `coding_sequence_variant`, or Ensembl VEP's holds it, alone or beside an in-frame or protein-altering term, where vep-rs's holds `frameshift_variant` with or without `stop_gained`.

A tuple carries no cDNA or CDS position, so the rule reads the rows of the two outputs that carry the pair's two tuples. The Allele, the record's alleles in the first column and the Location show whether the pair is an insertion's or a deletion's. `cDNA_position` and `CDS_position` show a row's cDNA positions and its CDS start. The pair is set aside when it is an insertion's or a deletion's and every such row shows that the test is made: both cDNA positions, a CDS start and an allele of `-` or the four bases alone.

`tools/vep-outputs-compare`, the Rust build of this comparator, applies the same rule to the same rows. The structural-variant comparator (`scripts/validation/compare_sv_concordance.py`) applies it to the pairs of the structural-variant sets.
