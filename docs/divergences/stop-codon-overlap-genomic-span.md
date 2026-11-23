# The bases read as the stop codon where an intron splits it or an RNA edit completes it

**The defect in one sentence.** To decide whether an insertion or deletion overlaps the stop codon, Ensembl VEP 116.3 reads the last three genomic positions of the coding region. On a stop codon that an intron splits, those positions are intron bases beside the codon's last base. On a transcript whose stop codon an `_rna_edit` completes, they are the stop bases the genome holds with the sense bases before them. On the rows where the two engines differ, Ensembl VEP writes `coding_sequence_variant` and vep-rs writes `stop_retained_variant`, both beside `splice_region_variant`.

**Why this is a defect.** The comment above the routine calls it the same test as the cDNA routine beside it, with an insertion's length considered. The cDNA routine reads the last three coding positions of the transcript, which are the stop codon. On a split stop codon the genomic routine misses the stop bases that end the exon before the intron. A deletion of one of them fails the overlap, and the row cannot carry `stop_retained_variant`, whatever the edited sequence holds. On the record of this page, Ensembl VEP's own `stop_lost` test reads cDNA coordinates and finds a stop codon in that edited sequence.

**What vep-rs does instead.** vep-rs reads the overlap of a deletion or a replacement on the genomic bases of the last three coding positions, the codon's own bases, wherever the exons place them. Where the span overlaps one of them, it makes the stop codon test on the edited sequence, and writes `stop_retained_variant` where a stop codon stands at the old stop position. On a transcript whose stop codon an `_rna_edit` completes, it reads the stop bases the genome holds and no sense base. An insertion it reads on the coding sequence: the overlap holds where the insertion lies before or inside the stop codon and every codon before that codon keeps its residue. `stop_codon_overlap_insertion_range` describes that reading.

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `ensembl` at release/116, the modules the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1418-1422, `_overlaps_stop_codon_cil`: the routine takes the transcript's `coding_region_start` and `coding_region_end`. It reads the three genomic positions at the coding region's end as the stop codon: `$t_end - 2` to `$t_end` on the forward strand, `$t_end` to `$t_end + 2` on the reverse.
- The same routine, lines 1403 and 1424-1427: the variant's genomic start and end are tested for overlap with those three positions. Lines 1408-1415 move an insertion's start first, the subject of `stop_codon_overlap_insertion_range`.
- The same file, line 1388, the comment above the routine: `# same as overlaps_stop_codon but "consider insertion length" to see overlap`. Line 1402 gives the reason for genomic coordinates: "because of intron we need to check overlap based on genomic coords when adding inserted seq length".
- The same file, lines 1378-1381, `_overlaps_stop_codon`: the cDNA routine tests the variant's cDNA span against `cdna_coding_end - 2` to `cdna_coding_end`, the stop codon's own three positions.
- `Bio/EnsEMBL/Transcript.pm`, lines 1144-1145, `coding_region_end`: on the forward strand the value is the start of the translation's end exon plus the translation's end minus one. It is counted inside the exon that holds the codon's last base. Lines 1101-1102 do the same for `coding_region_start` on the reverse strand.
- The same file, lines 1038-1057: a transcript's sequence edits move `cdna_coding_end`, and `coding_region_end` reads none of them.
- `VariationEffect.pm`, lines 1313 and 1318, `stop_retained`: for an insertion or deletion whose alternate peptide is undefined, empty or holds an `X`, the term needs the genomic overlap and a stop codon test that finds the codon unaltered. Line 1499: the stop codon test, `_ins_del_stop_altered_cil`, asks the same overlap first. Without it the test returns 0, the answer it also gives for a codon it finds unaltered, and reads no base. Lines 1528-1535 translate the three bases at the old stop position.
- The same file, line 1635, `coding_unknown`: the routine reads `stop_retained`, and where that is false on such a row it writes `coding_sequence_variant`.

**What correcting it would change.** Read on the stop codon's own bases, the overlap holds for a deletion of a stop base before the intron. With that reading alone in Ensembl VEP's routine, `Consequence` changes on 115 rows of the genome-wide ClinVar dataset on GRCh38 and 17 on GRCh37.

On 114 and 16 of them the row would read `stop_retained_variant` only because the stop codon test returns before it tests anything, the defect of `stop_retained_without_stop_codon_test`. The test returns that way where the record runs from the stop codon across the exon's end into the intron, and so has no cDNA position at one end, or where its allele holds a letter that names no base. With that defect corrected as well, they keep `coding_sequence_variant`. On one row of each assembly the record lies inside the exon and has its cDNA positions, so the test is made. It finds a stop codon, and the row reads `stop_retained_variant`, as vep-rs writes it.

On a whole stop codon the two readings are the same three bases.

**The class it produces.** `stop_codon_overlap_genomic_span`: a `Transcript` row of an insertion or deletion on which Ensembl VEP writes `coding_sequence_variant` and vep-rs `stop_retained_variant`, both beside `splice_region_variant`. Every other term is equal. Both rows have both cDNA positions, one for each end of the variant's span and a single value where the span is one base, and a CDS start. The golden corpora record each of their rows of the shape under the class's name. The format comparator sets a row of the class aside in a comparison that names the class; 1.7 says which rows and fields.

**How often the shape occurs.** Tier 1, see [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data carries the shape, and the divergence needs no flag or plugin. A record shows this much of the shape without its transcript: it is an insertion or deletion. The page linked above, the index of the entries, counts each entry's shape over five sets of chromosome 21 records: three cut from population datasets and two from ClinVar. Of those records, 137,713 per million have an allele of another length than the reference: an insertion, a deletion, or a replacement of the removed bases by another number of bases.

The transcript decides whether the divergence occurs. The Ensembl VEP caches of release 116 hold 2,327 transcripts on GRCh38, and 436 on GRCh37, whose stop codon an intron splits and whose coding sequence is complete at its 3' end.

Of the rows the two engines write for the genome-wide ClinVar dataset, the class sets aside 9 rows on GRCh38 and 3 on GRCh37. One row of each assembly is a deletion inside a split stop codon, the record of this page. The others are rows of insertions, which 1.7 names.

The second shape is carried by six RefSeq mitochondrial transcripts on GRCh38: ND1, ND2, COX3, ND3, ND4 and CYTB. On those six no row changes its term through this reading, and 1.4 says why.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. The second base of a split stop codon deleted, on ENST00000673897

**In short.** `3 36993662 2673709 TG T` on GRCh38 deletes one base, 3:36993663, the second base of the stop codon of the forward-strand transcript ENST00000673897. An intron of 7,291 bases splits the codon after that base. Ensembl VEP 116.3 reads 3:37000953-37000955 as the stop codon, two intron bases and the codon's third base, finds no overlap and writes `coding_sequence_variant`. vep-rs reads the codon's own three bases, finds a stop codon at the old stop position of the edited sequence, and writes `stop_retained_variant`. Both rows carry `splice_region_variant` and `NMD_transcript_variant`, and the other 73 rows of the record are equal.

### 1.1 What a reader needs to know first

A transcript's coding sequence is read from its exons joined end to end. A coding base is counted from the first base of that sequence, and codon n is coding bases 3n - 2 to 3n. The stop codon is the last codon: it ends the protein and adds no residue to it. A cDNA position counts the same bases from the transcript's first base, the 5' UTR included. The 3' UTR is the part of the transcript after the stop codon.

An intron can lie between two bases of one codon. The codon's bases then stand in two exons, and on the genome the whole intron separates them. This page calls such a stop codon a split stop codon, and a stop codon inside one exon a whole stop codon.

Three terms matter here. `stop_lost` says that a change leaves no stop codon where the stop codon stood. `stop_retained_variant` says that a change touches the stop codon and a stop codon stands in its place afterwards. `coding_sequence_variant` is the term Ensembl VEP's `coding_unknown` writes on these rows when neither of the two holds.

For an insertion or deletion whose changed codon cannot be translated, Ensembl VEP decides `stop_retained_variant` in two steps. An overlap test asks whether the variant touches the stop codon. A stop codon test then writes the change into the coding sequence with the 3' UTR after it, and translates the three bases at the old stop position, the place of the stop codon's three coding positions. The term is written when the overlap holds and those three bases are a stop codon. This page is about the overlap test.

### 1.2 The record

ClinVar record 2673709, a deletion of one base, against the Ensembl GRCh38 cache of release 116:

```
3	36993662	2673709	TG	T	.	.	.
```

Both engines drop the shared first base, so the record is one base, `G`, deleted at 3:36993663. ENST00000673897.1 is a forward-strand transcript of MLH1 with 18 exons and the biotype `nonsense_mediated_decay`. Three of its parts matter here:

| Part | Genomic span | Bases |
| --- | --- | --- |
| Exon 1 | 3:36993531-36993663 | 133 |
| Intron 1 | 3:36993664-37000954 | 7,291 |
| Exon 2 | 3:37000955-37001053 | 99 |

The coding sequence opens at cDNA 18 and ends at cDNA 134: 134 - 18 + 1 = 117 bases, 38 codons and the stop codon, 3 × 38 + 3 = 117. The stop codon is coding bases 115 to 117, cDNA 132 to 134.

Exon 1 holds cDNA 1 to 133. The codon's first two bases are therefore the exon's last two, 36993531 + 131 = 36993662 and 36993531 + 132 = 36993663. Its third base, cDNA 134, is the first base of exon 2, 37000955. The codon reads `T`, `G`, `A`. The record deletes the `G`, cDNA 133 and coding base 133 - 18 + 1 = 116.

The base after the stop codon, cDNA 135 at 37000956, is an `A`, the first base of the 3' UTR. No flag beyond the cache is needed.

### 1.3 What each engine writes

The row on ENST00000673897 in the default format. Ensembl VEP 116.3:

```
2673709	3:36993663	-	ENSG00000076242	ENST00000673897	Transcript	splice_region_variant,coding_sequence_variant,NMD_transcript_variant	133	116	39	*/X	tGa/ta	-	IMPACT=LOW;STRAND=1
```

vep-rs:

```
2673709	3:36993663	-	ENSG00000076242	ENST00000673897	Transcript	splice_region_variant,stop_retained_variant,NMD_transcript_variant	133	116	39	*/X	tGa/ta	-	IMPACT=LOW;STRAND=1
```

The rows differ in one term of the `Consequence` column: `coding_sequence_variant` in Ensembl VEP's row and `stop_retained_variant` in vep-rs's. Every other column is equal. Both rows place the deletion at cDNA 133, coding base 116 and residue 39, the stop codon. Both write the codon change as `tGa/ta` and the amino acids as `*/X`.

The record writes 74 rows in each engine, and the other 73 are equal in every column. Twelve of them are rows of other transcripts of the gene whose stop codon is split after the same two bases:

| Third stop base | Transcripts | Both engines |
| --- | --- | --- |
| 3:36996624 | 8 | `stop_lost` |
| 3:36995389 | 1 | `stop_lost` |
| 3:36994830 | 3 | `stop_lost` |

Each of the twelve rows reads `stop_lost,splice_region_variant,NMD_transcript_variant` in both engines, with the same codon change and amino acids as the row above.

### 1.4 Why Ensembl VEP writes that

The deletion leaves two bases of the codon, `ta` in the row's `Codons` column, and the row's alternate amino acid is `X`. An alternate peptide that holds an `X` sends `stop_retained` to its second arm (lines 1313 and 1318): the term needs the genomic overlap test and the stop codon test.

`stop_lost` is asked first (line 1297). Its test (lines 1433-1484) reads cDNA coordinates and finds a stop codon at the old stop position, so the row carries no `stop_lost`.

The overlap test takes the variant's genomic start and end, both 36993663 (line 1403). `coding_region_end` is the start of the translation's end exon, exon 2, plus the translation's end there, 1, minus one: 37000955 + 1 - 1 = 37000955 (`Transcript.pm`, lines 1144-1145). Line 1421 reads the three positions from 37000955 - 2 = 37000953 to 37000955 as the stop codon.

Intron 1 ends at 37000954. The positions 37000953 and 37000954 are its last two bases, and 37000955 is the codon's third base. The deleted base lies 37000953 - 36993663 = 7,290 positions before them, so lines 1424-1427 find no overlap and the routine returns 0.

With the overlap at 0, line 1318 makes `stop_retained` false. The stop codon test is not reached: it asks the same overlap first and returns (line 1499). `frameshift` returns 0 as well, since the reference peptide begins with the stop (line 1554). `coding_unknown` reads `stop_retained` (line 1635) and writes `coding_sequence_variant`.

The cDNA routine answers the other way for the same record. It tests cDNA 133 against 134 - 2 = 132 to 134 (lines 1378-1381), and the overlap holds. Ensembl VEP's two routines, one called the same test as the other, return 0 and 1 for one deletion.

The second shape is read from the same lines and has no record on this page. Six RefSeq mitochondrial transcripts on GRCh38 end their last exon one or two bases into the stop codon, and an `_rna_edit`, a sequence edit on the transcript, supplies the rest. `coding_region_end` is the exon's last base, since no edit moves it. On ND1 the exon ends at MT:4262 and the routine reads MT:4260-4262: the last sense base and the two stop bases the genome holds. On ND2, COX3, ND3, ND4 and CYTB it reads two sense bases and the one stop base the genome holds.

For a deletion of MT:4260 on ND1, Ensembl VEP's genomic routine returns 1 and its cDNA routine 0. Both engines write `frameshift_variant` for it. No row of the six transcripts changes its term through the reading. No base follows their stop codon, so such a deletion shortens the edited sequence, and the stop codon test finds the codon altered whichever bases the overlap reads.

All six are also forward-strand transcripts whose stop codon is a stop codon in both the standard and the mitochondrial codon table. Read from the source lines this section cites, the two readings of the overlap give Ensembl VEP a different term on a reverse-strand transcript, or for a stop codon that the mitochondrial table alone holds. The six are neither.

### 1.5 Why vep-rs's row is the correct one

The stop codon of ENST00000673897 is cDNA 132 to 134. Read through the exons, those are 36993662 and 36993663 in exon 1 and 37000955 in exon 2. The record deletes 36993663, so it removes a base of the stop codon. An overlap test that answers no for it is answering about other bases.

The three positions Ensembl VEP reads are 37000953 to 37000955. Two of them are bases of intron 1, 36993664 to 37000954, and have no cDNA position: 37000954 - 36993664 + 1 = 7,291 bases that the transcript does not hold. They are no part of any codon.

The edited sequence shows what the deletion does. The 114 coding bases before the stop codon, 3 × 38 = 114, are untouched, so the 38 residues are the same. At the old stop position, coding bases 115 to 117, the transcript holds the codon's `T`, its `A`, and the base that followed the codon, the `A` at 37000956. The three bases read `TAA`, a stop codon. The protein ends where it ended. A change inside the stop codon that leaves a stop codon in its place is what `stop_retained_variant` says.

Ensembl VEP reads that sequence the same way in its `stop_lost` test, which finds the stop codon at the old stop position and writes no `stop_lost`. Its row then says neither that the stop codon is lost nor that it is kept.

The twelve other rows show the two engines agreeing on a split stop codon. On those transcripts the codon is the same `T`, `G`, `A`, split after the `G`, and the base that follows it is a `T` or a `C`. With the `G` deleted the three bases read `TAT` or `TAC`, no stop codon, and both engines write `stop_lost`. Ensembl VEP reads those split codons whole, through cDNA coordinates. On ENST00000673897 alone the three bases are a stop codon, the term to write is `stop_retained_variant`, and only that term's overlap test reads genomic positions.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it, from the repository's root:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache tests/golden/116/GRCh38/json_cache \
  --species homo_sapiens --assembly GRCh38
grep ENST00000673897 out.txt
```

The committed GRCh38 cache holds the transcript. vep-rs writes `splice_region_variant,stop_retained_variant,NMD_transcript_variant` on the row. Ensembl VEP 116.3 takes the record with the cache options of [Reproducing these examples](reproducing.md) and its GRCh38 cache, and writes `splice_region_variant,coding_sequence_variant,NMD_transcript_variant`.

The entry's reproduction, `scripts/divergences/reproductions/stop_codon_overlap_genomic_span.toml`, runs the record on both engines. It holds the row of ENST00000673897 to both recorded values in six columns: `Location`, `Feature`, `Consequence`, `cDNA_position`, `CDS_position` and `IMPACT`.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) compares one Ensembl VEP output with one vep-rs output, every field of every row. A comparison case is one such comparison: an input, the flags and the classes it applies. Where a comparison case names the class `stop_codon_overlap_genomic_span`, the comparator sets a row of the class aside and counts it once under that name.

The golden corpora record the class's rows. A corpus manifest lists each one under `divergences` with the class's name and the terms vep-rs writes. On a listed row the golden harness checks vep-rs's row against those terms in place of Ensembl VEP's. Five manifests hold 201 such entries, all rows of insertions at four sites. No committed corpus holds a deletion of the class; the record of this page is ClinVar's.

The class's predicate is `stop_codon_overlap_genomic_span_pair`. The class takes a `Transcript` row of an insertion or deletion on four conditions. The two sets of terms are equal but for Ensembl VEP's `coding_sequence_variant` against vep-rs's `stop_retained_variant`. Both sets hold `splice_region_variant`, which every exon base within three of an exon's end carries. Both rows have both cDNA positions and a CDS start. The allele is `-` or bases `A`, `C`, `G` and `T` alone.

A row away from an exon's end, a row without its positions and the two terms the other way round are no rows of the class. The comparator compares those rows like any other and reports a difference in their `Consequence`.

Of Ensembl VEP 116.3's rows against vep-rs's for the genome-wide ClinVar dataset, the class sets aside 9 rows on GRCh38 and 3 on GRCh37. One row of each assembly is the deletion of this page, on ENST00000673897 on GRCh38 and on ENST00000442249 on GRCh37.

Seven rows on GRCh38 and two on GRCh37 are rows of insertions on reverse-strand transcripts whose stop codon an intron splits. vep-rs writes the term on them through its reading of an insertion on the coding sequence, which belongs to `stop_codon_overlap_insertion_range`, and not through the reading of a deletion that this page describes. The class sets them aside because their rows have its shape.

The ninth row on GRCh38 is ClinVar 2839751 on ENST00001038946, an insertion beside a whole stop codon that opens its exon at 5:138930833-138930835. Its mechanism is the one of `stop_codon_overlap_insertion_range`. The class sets the row aside for its `splice_region_variant`.

On a row the class sets aside, the comparator also sets aside the fields that follow the two terms, by the rule `term_pair_following_cell`. `IMPACT` follows where each side's value is the impact of its own most severe term. A field of the Downstream, NMD or LOFTEE plugin follows where it is written on the side whose terms let the plugin run and left out on the other. An `HGVSp` follows where one side alone writes a frameshift's notation, and that side is the one without `stop_retained_variant`. The rule names no other field. Every other field of the row is compared, and a difference in it is reported.

On the row of this page no field follows the terms: `IMPACT` is `LOW` in both engines.

The SNP/indel comparator (`tools/vep-outputs-compare`) compares tuples: a row's Location, Allele, Feature, Feature_type and set of `Consequence` terms. Two tuples, one from each engine, that differ in the set of terms alone are a pair. The comparator's 116.3 registry holds a pair rule under the class's name, with the same predicate, `stop_codon_overlap_genomic_span_pair`. The rule sets aside both tuples of a pair whose sets differ by Ensembl VEP's `coding_sequence_variant` against vep-rs's `stop_retained_variant` alone.

A tuple carries no cDNA or CDS position, so the rule reads the rows of the two outputs that carry the pair's two tuples. The Allele, the record's alleles in the first column and the Location show whether the pair is an insertion's or a deletion's. `cDNA_position` and `CDS_position` show a row's cDNA positions and its CDS start. The pair is set aside when it is an insertion's or a deletion's, both sets hold `splice_region_variant`, and every such row has both cDNA positions, a CDS start and an allele of `-` or bases alone.

The structural-variant comparator (`tools/vep-sv-compare`) applies it to the pairs of the structural-variant sets.
