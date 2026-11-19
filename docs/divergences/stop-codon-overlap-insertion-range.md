# An insertion that keeps the stop codon on a forward-strand transcript

**The defect in one sentence.** Before Ensembl VEP 116.2 tests whether an insertion overlaps the stop codon, it moves the insertion's start by the inserted length, forward on a forward-strand transcript and back on a reverse-strand one. On the reverse strand the test then holds for an insertion inside the codon, or before it by fewer bases than it inserts, whatever the inserted bases are. On the forward strand it holds only for one base inserted between the codon's first and second base. An insertion that keeps the stop codon is therefore `stop_retained_variant` on one strand and `coding_sequence_variant` or `frameshift_variant` at the same coding position on the other. vep-rs writes `stop_retained_variant` for it on either strand.

**Why this is a defect.** A transcript's strand does not change what an insertion does to its coding sequence. For one coding position the two moves answer differently. On the reverse strand the moved span is the inserted length long and ends at the insertion site, which is what the comment above the routine announces: "consider insertion length". On the forward strand the move takes the start further past the end, and the span meets the codon in one case alone. A forward-strand insertion that leaves every residue before the stop codon as it was, and a stop codon at the old stop position, then misses the term that the same insertion carries on a reverse-strand transcript.

**What vep-rs does instead.** vep-rs reads an insertion on the coding sequence, the same way on either strand. It takes the coding position that follows the insertion, f, and the coding length, n. The insertion reaches the stop codon where it lies before or inside the codon, f ≤ n, and every codon before the old stop position reads the reference's residue in the edited sequence. The inserted bases need not repeat the coding bases before the codon. Where they differ, the insertion reaches the codon as long as every codon before the old stop position keeps its residue. It does not reach the codon where a residue before it changes.

vep-rs then makes the stop codon test: the three bases at the old stop position of the edited sequence are read, and a stop codon there gives `stop_retained_variant`. An insertion that the record places after the codon has no coding position and no coding term.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, the module the Ensembl VEP 116.2 image holds):

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1408-1415, `_overlaps_stop_codon_cil`: for an insertion of the letters `A`, `C`, `G`, `T` or `N`, the start is moved by the inserted length. Line 1410 adds it on strand 1, and line 1413 subtracts it on strand -1. The comment of line 1388 says the routine is to "consider insertion length".
- The same routine, line 1403: the span is the variation feature's start and end. An insertion's start is its end plus one.
- The same routine, lines 1424-1427: the moved span is tested for overlap with the last three bases of the genomic coding region.
- The same file, line 1318, `stop_retained`: for an allele that changes the length, the term needs this overlap and a stop codon test that finds the codon unaltered. Lines 1313-1315 decide `stop_retained_variant` for an allele whose alternate peptide is whole, with every codon complete. They do not test the overlap.
- The same file, lines 1521 and 1528-1535, `_ins_del_stop_altered_cil`: the stop codon test writes the allele into the coding sequence, followed by the 3' UTR, at its CDS start. It then reads the three bases at the old stop position.
- The same file, line 1549, `frameshift`: the routine returns 0 where `stop_retained` holds.
- The same file, lines 1616-1638, `coding_unknown`: the general term `coding_sequence_variant`, written where neither stop codon term, `stop_retained_variant` or `stop_lost`, holds.

**What correcting it would change.** An overlap read at the insertion's coding position, the same on both strands, would give the forward-strand rows of this class `stop_retained_variant`, the term the reverse-strand rows of the same shape carry. The cells that are read from the terms would follow. They are `IMPACT`, the fields that the Downstream, NMD and LOFTEE plugins write on a row whose terms let them run, and the form of the `HGVSp`.

**The class it produces.** `stop_codon_overlap_insertion_range`: a `Transcript` row of an insertion with both cDNA positions and a CDS start in both engines. Its two `Consequence` sets are equal but for one of two pairs. The first is Ensembl VEP's `coding_sequence_variant` against vep-rs's `stop_retained_variant`, on a row without `splice_region_variant`. The second is Ensembl VEP's `frameshift_variant`, with or without `stop_gained`, against vep-rs's `stop_retained_variant`, with or without `inframe_insertion`. The golden corpora record each such row under the class's name. The format comparator sets a row of the class aside in a comparison that names the class. Section 1.7 says which rows and fields.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data carries the shape, and the divergence needs no flag or plugin. The shape the record alone shows is an insertion. Of the 15,646,481 chromosome 21 records over which the index, the page linked above, counts each entry's shape, 996,643 are insertions, 63,698 per million. That rate is an upper bound: the transcript decides whether an insertion lies before or inside a stop codon and keeps it.

On the genome-wide ClinVar dataset, reading an insertion on the coding sequence in place of the moved span changes the `Consequence` on 1,166 of 125,572,915 rows on GRCh38 and on 182 of 33,483,969 on GRCh37. Every one gains `stop_retained_variant`. Of these, 1,159 and 180 are rows of forward-strand transcripts. The other 7 and 2 are rows of reverse-strand transcripts whose stop codon an intron splits. Of Ensembl VEP 116.2's rows against vep-rs's, the class sets aside 1,158 rows on GRCh38 and 180 on GRCh37. The other 8 and 2 carry `splice_region_variant` and have the shape of the class `stop_codon_overlap_genomic_span`.

On 1,077 of the forward-strand rows on GRCh38 and 169 on GRCh37, the inserted bases repeat the transcript's bases from the insertion through the codon's last base. The same change can then also be written as an insertion after the stop codon. On 82 and 11 rows it cannot. Section 1.5 says what that means for the worked record. Nine golden corpora record 270 entries under the class.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Four bases inserted inside a stop codon, ClinVar 1025066 on ENST00000381136

**In short.** ClinVar 1025066 on GRCh38 inserts four bases between 9:12709180 and 9:12709181, the first and second base of the stop codon of the forward-strand transcript ENST00000381136. The first two inserted bases repeat the codon's second and third base. The edited coding sequence therefore holds the same stop codon at the same place, and the protein keeps its 247 residues. Ensembl VEP 116.2 writes `inframe_insertion,coding_sequence_variant`, and vep-rs writes `inframe_insertion,stop_retained_variant`. Ensembl VEP's overlap test moves the insertion's start four bases forward, past the codon, and fails. A second record inserts four bases at the same place of the stop codon of the reverse-strand transcript ENST00000340345. Both engines write `inframe_insertion,stop_retained_variant` for it.

### 1.1 What a reader needs to know first

A coding sequence ends with a stop codon, three bases at which translation ends. Coding positions are counted from the first base of the coding sequence, and a row's `CDS_position` gives them. In a coding sequence of n bases the stop codon is positions n - 2 to n. The codons before it give the protein's residues.

`stop_retained_variant` is the term for a change that touches the stop codon and leaves a stop codon in its place. For an insertion or deletion Ensembl VEP decides it in two steps. The first asks whether the allele overlaps the stop codon, and this page calls it the overlap test. The second, the stop codon test, writes the allele into the coding sequence and reads the three bases at the old stop position. The term is written where the overlap holds and those three bases are a stop codon.

Where the overlap fails, the stop codon test is not asked. The row then carries the terms of Ensembl VEP's other coding tests. One of them is the general term `coding_sequence_variant`, written where no specific term holds.

A VCF record writes an insertion with the base before it. Both engines drop that shared base and give the row a `Location` of two neighbouring positions, the bases beside the insertion. Ensembl VEP stores such an insertion with a start one above its end.

A forward-strand transcript runs toward higher genomic coordinates, and a reverse-strand transcript toward lower ones. The coding position that follows an insertion is therefore the higher of its two neighbours on the forward strand and the lower on the reverse strand.

### 1.2 The record

Two records, against the Ensembl GRCh38 cache of release 116. The first is ClinVar 1025066. The second is the record `stop_ins4_twostops_rev` of the golden corpus `tests/golden/116/GRCh38-release116`:

```
9	12709180	1025066	T	TAACA	.	.	.
21	30282169	stop_ins4_twostops_rev	C	CTATT	.	.	.
```

Both engines drop the shared first base, so each record is four bases inserted between two neighbouring positions. Each insertion lies between the first and second base of a stop codon:

| Part | ClinVar 1025066 | Second record |
| --- | --- | --- |
| Transcript | ENST00000381136 | ENST00000340345 |
| Strand | forward | reverse |
| Stop codon | 9:12709180-12709182 | 21:30282170-30282168 |
| Coding positions of the codon | 742-744 | 763-765 |
| Bases beside the insertion | 9:12709180, 9:12709181 | 21:30282170, 21:30282169 |
| Coding position after it | 743 | 764 |

ENST00000381136 is a transcript of TYRP1 with five exons. Exons 1 to 4 hold 118 + 168 + 180 + 147 = 613 bases, and exon 5 is 9:12708977-12710597. Position 12709180 is base 12709180 - 12708977 + 1 = 204 of exon 5, cDNA 613 + 204 = 817. The coding sequence opens at cDNA 76, so that base is coding position 817 - 76 + 1 = 742. The coding sequence has 744 bases, 247 codons and the stop codon: 3 × 247 + 3 = 744. Positions 742 to 744, 9:12709180-12709182, are the stop codon, `TAA`.

The first record inserts `AACA` between 9:12709180 and 9:12709181, coding positions 742 and 743. The coding position that follows it is f = 743, and n = 744.

ENST00000340345 is a transcript of KRTAP24-1 with one exon, 21:30281309-30282950, on the reverse strand. Its cDNA 1 is 30282950, and its coding sequence is cDNA 19 to 783, 783 - 19 + 1 = 765 bases: 3 × 254 + 3 = 765. The stop codon is cDNA 781 to 783. Its first base is 30282950 - 781 + 1 = 30282170, coding position 781 - 19 + 1 = 763, and its last base is 30282168. On the transcript the codon reads `TGA`.

The second record inserts `TATT` between 21:30282169 and 21:30282170, which the transcript carries as `AATA`. On this strand 30282170 is coding position 763 and 30282169 is 764, so f = 764 and n = 765.

Three other transcripts of TYRP1 hold the same stop codon at coding positions 1612 to 1614: ENST00000388918, ENST00001078217 and ENST00001138527. The first record lies between their coding positions 1612 and 1613.

### 1.3 What each engine writes

The row of the first record on ENST00000381136 in the default format. Ensembl VEP 116.2:

```
1025066	9:12709180-12709181	AACA	ENSG00000107165	ENST00000381136	Transcript	inframe_insertion,coding_sequence_variant	817-818	742-743	248	*/*QX	taa/tAACAaa	-	IMPACT=MODERATE;STRAND=1
```

vep-rs:

```
1025066	9:12709180-12709181	AACA	ENSG00000107165	ENST00000381136	Transcript	inframe_insertion,stop_retained_variant	817-818	742-743	248	*/*QX	taa/tAACAaa	-	IMPACT=MODERATE;STRAND=1
```

The rows differ in the `Consequence` column alone. Ensembl VEP writes `coding_sequence_variant` where vep-rs writes `stop_retained_variant`, each beside `inframe_insertion`. The position columns, `Amino_acids` and `Codons` are equal. `Amino_acids` opens with the stop, `*`, on both sides of its slash.

The row of the second record on ENST00000340345, which both engines write alike:

```
stop_ins4_twostops_rev	21:30282169-30282170	TATT	ENSG00000188694	ENST00000340345	Transcript	inframe_insertion,stop_retained_variant	781-782	763-764	255	*/**X	tga/tAATAga	-	IMPACT=MODERATE;STRAND=-1
```

The first record writes 13 rows in each engine and the second writes one. Five are rows of coding transcripts whose stop codon holds the insertion:

| Transcript | Strand | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- | --- |
| ENST00000381136 | forward | `coding_sequence_variant` | `stop_retained_variant` |
| ENST00000388918 | forward | `coding_sequence_variant` | `stop_retained_variant` |
| ENST00001078217 | forward | `coding_sequence_variant` | `stop_retained_variant` |
| ENST00001138527 | forward | `coding_sequence_variant` | `stop_retained_variant` |
| ENST00000340345 | reverse | `stop_retained_variant` | `stop_retained_variant` |

Every one of the five carries `inframe_insertion` as well, and an `IMPACT` of `MODERATE`. The other nine rows of the first record are equal in every column. They are rows of two transcripts without a coding sequence and of seven transcripts of another gene that hold the record in an intron.

### 1.4 Why Ensembl VEP writes that

The first record changes the length of the allele. The alternate peptide of its row, the right side of `Amino_acids`, ends in `X`, a codon that is not whole. The term is therefore decided by line 1318, which asks the overlap test first.

The overlap test takes the span of the variation feature (line 1403). The insertion lies between 12709180 and 12709181, so its end is 12709180 and its start 12709181. Its inserted bases are among the letters for which lines 1408-1415 move the start. ENST00000381136 is on strand 1, so line 1410 adds the inserted length: 12709181 + 4 = 12709185.

Lines 1424-1427 test the moved span against the last three bases of the genomic coding region, 12709180 to 12709182. The moved start, 12709185, lies above 12709182, the codon's last base. The overlap fails. Line 1318 is false, and the stop codon test is not asked. The row's coding term is the general one of lines 1616-1638, `coding_sequence_variant`.

On the forward strand a one-base insertion at this place would pass, since 12709181 + 1 = 12709182 is the codon's last base. Every longer insertion moves the start past the codon.

The second record lies between 30282169 and 30282170, so its end is 30282169 and its start 30282170. ENST00000340345 is on strand -1, and line 1413 subtracts the inserted length: 30282170 - 4 = 30282166. The moved span, 30282166 to 30282169, is four bases long and ends at the insertion site. The last three bases of the genomic coding region are 30282168 to 30282170. The span holds two of them, 30282168 and 30282169, and the overlap holds.

The stop codon test follows. It writes the four bases into the coding sequence at the CDS start (line 1521) and reads the three bases at the old stop position (lines 1528-1535). They are a stop codon, so line 1318 is true and the row carries `stop_retained_variant`.

Moved back as on the reverse strand, the first record's start would be 12709181 - 4 = 12709177. The span 12709177 to 12709180 would hold 12709180, the codon's first base, and the overlap would hold. The two records lie at the same place of their stop codons and insert the same number of bases. The direction of the move alone separates their answers.

### 1.5 Why vep-rs's row is the correct one

The coding sequence of ENST00000381136 has 744 bases, and its stop codon, `TAA`, is positions 742 to 744. The record places `AACA` between positions 742 and 743. The edited sequence is the reference's positions 1 to 742, the four inserted bases, and the reference's positions from 743 on:

| Edited positions | Bases | Read as |
| --- | --- | --- |
| 1-741 | the reference's 1-741 | 247 codons, unchanged |
| 742-744 | the reference's 742, `T`, then the inserted `AA` | `TAA`, a stop codon |
| 745-748 | the inserted `CA`, then the reference's 743-744, `AA` | after the stop codon |

The first two inserted bases, `AA`, repeat the reference's positions 743 and 744. Positions 742 to 744 of the edited sequence therefore read `TAA`, the stop codon the reference has there. The 741 bases before position 742, 3 × 247, are untouched. The edited sequence gives the same 247 residues and ends at the same stop position. The edited sequence is four bases longer, 744 + 4 = 748, and its positions 745 to 748 lie after the stop codon.

A stop codon that the change touches and leaves in place is what `stop_retained_variant` states. `coding_sequence_variant` states only that the coding sequence changed. Ensembl VEP's own row shows the stop, `*`, first on both sides of its `Amino_acids` column.

Ensembl VEP's stop codon test would give the same answer. It reads the three bases at the old stop position, here `TAA`. The test is not reached, because the overlap test fails before it. On the second record, where the overlap holds, the same test gives the term in both engines.

The second record repeats the arithmetic on the reverse strand. The coding sequence of ENST00000340345 has 765 bases, and its stop codon, `TGA`, is positions 763 to 765. The record places `AATA` between positions 763 and 764. Positions 763 to 765 of the edited sequence hold the codon's first base and the first two inserted bases, `TAA`, a stop codon. The 762 bases before it, 3 × 254, are untouched, and the protein keeps its 254 residues.

For one coding position, between the first and second base of the stop codon, the two records make the same kind of change with the same effect on the protein. vep-rs reads both on the coding sequence. For the first record f = 743 and n = 744, and for the second f = 764 and n = 765. In both, f ≤ n and no coding base lies between the insertion and the codon. Both then pass the stop codon test.

One more fact belongs to the first record. Its first two inserted bases repeat the codon's last two, so the same change can also be written two bases on, between 9:12709182 and 9:12709183, after the codon. A record written there is an insertion after the coding sequence, with no coding position and no coding term. ClinVar 1025066 names the place inside the codon, and both engines read it there: `Location` and `CDS_position`, 742-743, are equal in the two rows. Under `--shift_3prime 1` both engines move an insertion of this kind to its last place. The 1,077 and 169 such rows of the genome-wide ClinVar dataset then carry no coding position in either engine.

### 1.6 Try it yourself

With the two records of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38/json_cache \
  --species homo_sapiens --assembly GRCh38
grep -E 'ENST00000381136|ENST00000340345' out.txt
```

The committed GRCh38 cache holds the transcripts of both records, and no flag beyond the cache is needed. vep-rs writes `inframe_insertion,stop_retained_variant` on both rows. Ensembl VEP 116.2 takes the same records with the cache options of [Reproducing these examples](reproducing.md) and its GRCh38 cache. It writes `inframe_insertion,coding_sequence_variant` on the first row and `inframe_insertion,stop_retained_variant` on the second.

The entry's reproduction, `scripts/divergences/reproductions/stop_codon_overlap_insertion_range.toml`, runs both records on both engines. It holds five rows to their recorded values: the four forward-strand rows that differ and the reverse-strand row that does not.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) compares one Ensembl VEP output with one vep-rs output, every field of every row. A comparison case is one such comparison: an input, the flags and the classes it applies. Where a comparison case names the class `stop_codon_overlap_insertion_range`, the comparator sets a row of the class aside and counts it once under that name. The predicate that recognises such a row is `stop_codon_overlap_insertion_range_pair`.

Of the columns that differ, the predicate reads the `Consequence` alone. The class takes a `Transcript` row of an insertion of the letters `A`, `C`, `G`, `T` and `N` with both cDNA positions and a CDS start on both sides, whose two sets of terms differ by one of the two pairs of the class. Of Ensembl VEP 116.2's rows against vep-rs's for the genome-wide ClinVar dataset, the class sets aside 1,158 rows on GRCh38 and 180 on GRCh37.

On a row the class sets aside, the comparator also sets aside the cells that follow the two terms, by the rule `term_pair_following_cell`. An `IMPACT` follows where each side's value is the impact of that side's most severe term. A field of the Downstream plugin follows where one side alone writes it, the side whose terms hold `frameshift_variant` and no splice site term. A field of the NMD or LOFTEE plugin follows where one side alone writes it, the side whose terms hold `stop_gained`, `frameshift_variant` or a splice site term.

An `HGVSp` follows where one side alone writes a frameshift's notation, and that side is the one without `stop_retained_variant`. The rule names no other cell. Every other cell of the row is compared, and a difference in it is reported.

On a row of the second pair the two engines reach the `HGVSp` by different routes. Ensembl VEP's terms hold `frameshift_variant`, so it walks the two translations residue by residue (`Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2348-2357). Where both reach a stop at one residue before any residue differs, it writes the stop at its residue number as unchanged (line 1911).

vep-rs compares the alternate protein with the reference protein, as [The HGVSp of a change at the stop codon](hgvsp-stop-codon-protein-comparison.md) describes. It does so for an insertion inside the stop codon, for one before the last residue's codon that brings a stop codon of its own to the stop's residue, and for one inside the last residue's codon whose length is no whole number of codons. Where it finds the two equal it writes the same value. The cell is then equal in both rows, and the rule is not asked.

Inside the last residue's codon, vep-rs writes the stop at its residue number as unchanged for an insertion whose length is no whole number of codons and that leaves the protein as it was. For one of a whole number of codons that leaves the protein as it was it writes no `HGVSp`.

vep-rs writes an insertion inside that codon that changes a residue or moves the stop as that change. On a made-up transcript whose stop is residue 12, it writes `p.Val11AlafsTer2` for five bases there that change residue 11 before a stop codon of their own, and `p.Ter12LysextTer3` for five that keep residue 11 and put a residue in the stop's place. On a reverse-strand transcript the two engines' terms are equal, `inframe_insertion,stop_retained_variant`, for five or eight bases there that leave the protein as it was, and the row is none of this class. Ensembl VEP writes no `HGVSp` on it and vep-rs writes the stop as unchanged.

The golden corpus builder (`scripts/golden/build_golden_corpus.py`) classes a corpus row of the shape under the class's name, in the `divergences` list of the corpus manifest, with the terms vep-rs writes. Nine corpora hold 270 such entries. On a corpus, the golden harness checks vep-rs's row against the recorded terms in place of Ensembl VEP's.

Rows beside the shape belong to other entries. The first pair beside `splice_region_variant` is the shape of `stop_codon_overlap_genomic_span`. Of the rows that vep-rs's reading of an insertion changes, the class of that name sets aside 8 on GRCh38 and 2 on GRCh37. One of them is a row of this mechanism, ClinVar 2839751 on ENST00001038946. Its stop codon, 5:138930833-138930835, is not split and opens its exon. The opposite pair of terms, Ensembl VEP's `stop_retained_variant` against vep-rs's `coding_sequence_variant` or `frameshift_variant` on a reverse-strand insertion before the codon, is the shape of [`stop_retained_three_base_stop_codon_test`](stop-retained-three-base-stop-codon-test.md).

Under `--shift_3prime 1` vep-rs's reading of an insertion changes 87 rows on GRCh38 and 13 on GRCh37 of the genome-wide ClinVar dataset. The comparator sets 80 and 11 of them aside under the class [`stop_codon_unshifted_span`](stop-codon-unshifted-span.md). That class is Ensembl VEP's `coding_sequence_variant` against vep-rs's `stop_retained_variant` on a row whose allele the flag moved. The other 7 and 2 have the shape of `stop_codon_overlap_genomic_span`, so none is left to this class. Of the 82 forward-strand rows of GRCh38 that cannot be written after the codon, the flag moves 2 into an intron, where both engines write `intron_variant` and no coding position.

The SNP/indel comparator (`scripts/concordance/compare_vep_outputs.py`) compares tuples: a row's Location, Allele, Feature, Feature_type and set of `Consequence` terms. Two tuples, one from each engine, that differ in the set of terms alone are a pair. The comparator's 116.2 registry holds a pair rule under the class's name, with the same predicate, `stop_codon_overlap_insertion_range_pair`. The rule sets aside both tuples of a pair whose sets differ by one of the two pairs of terms of the class.

The first is Ensembl VEP's `coding_sequence_variant` against vep-rs's `stop_retained_variant`, where neither set holds `splice_region_variant`. The second is Ensembl VEP's `frameshift_variant`, with or without `stop_gained`, against vep-rs's `stop_retained_variant`, with or without `inframe_insertion`.

A tuple carries no cDNA or CDS position, so the rule reads the rows of the two outputs that carry the pair's two tuples. The Allele, the record's alleles in the first column and the Location show whether the pair is an insertion's. `cDNA_position` and `CDS_position` show a row's cDNA positions and its CDS start. The pair is set aside when it is an insertion's of the letters `A`, `C`, `G`, `T` and `N` and every such row has both cDNA positions and a CDS start.

Under an identifier the rule reads an insertion from a Location of two neighbouring positions and an Allele of another length than two, so two bases replaced by a sequence of another length are set aside with the insertions.

`tools/vep-outputs-compare`, the Rust build of this comparator, applies the same rule to the same rows. The structural-variant comparator (`scripts/validation/compare_sv_concordance.py`) applies it to the pairs of the structural-variant sets.
