# The HGVSp of an insertion or deletion shifted from an intron into the coding sequence

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 writes the HGVSc of the shifted change and no HGVSp for an insertion or deletion written in an intron inside the coding region, or before the coding start, that the 3' shift carries wholly into the coding sequence. It fills the HGVSp only for one written in such an intron within 12 bases of an exon, on a transcript that has an intron of 13 bases or fewer. For such an insertion or deletion vep-rs writes the protein notation of the shifted change wherever the codons the transcript holds give it.

**Why this is a defect.** HGVS describes an insertion or deletion at its most 3' position, and both engines shift the allele there before writing a notation. On these rows that position lies in the coding sequence. The row's HGVSc describes the change there, and its `HGVS_OFFSET` gives the number of bases the shift moved it. Ensembl VEP reads the translation start and end, the first and last residue of the change, at the shifted position. The one test that empties the field asks whether the allele is coding, and Ensembl VEP makes it on the allele as the record writes it. For the same change given at the shifted position it writes the protein notation.

The protein notation is the predicted consequence of the change that the HGVSc describes. An empty HGVSp contradicts no HGVS recommendation: the row omits a description that the nomenclature gives. The defect is the inconsistency inside Ensembl VEP's own row, between the shift it applies and the coding test it makes on the allele as written.

On a transcript that has an intron of 13 bases or fewer, which Ensembl VEP's source calls a frameshift intron, the test widens the allele's span by 12 bases before it looks for an exon. An allele written in an intron inside the coding region, within 12 bases of an exon, then passes. Ensembl VEP writes its HGVSp when the translation start and end exist at the shifted position: the method computes the notation wherever its test lets the row through.

**What vep-rs does instead.** vep-rs writes the protein notation of the shifted change. It is the predicted consequence that the nomenclature gives for the description in the row's HGVSc, the one Ensembl VEP omits. Except at the stop codon and at [a start codon that stands again](defect-b.md#3-start_lost-and-start_retained_variant-emitted-together), it is the value both engines write for the change given at the shifted position.

vep-rs writes it where three things hold. The record writes the allele in an intron inside the coding region, or before the coding start: in the 5' UTR or an intron of it. The shifted change lies wholly in the coding sequence, as a deletion of coding bases or an insertion between two coding bases. The codons the transcript holds give the notation. Where the first two hold and Ensembl VEP's wider test fills the HGVSp, vep-rs writes the same notation.

Those codons always give it where the coding sequence ends on its stop codon. A coding sequence ends on its stop codon when it is a whole number of codons, the last of them a stop codon, and the transcript does not carry the flag `cds_end_NF`, which marks a coding sequence whose end is not found. One that ends inside a codon does not end on a stop codon, whatever its last three bases spell: the codons are counted from the start of the sequence.

The flag `cds_start_NF`, which marks a coding sequence whose start is not found, does not enter that test. Such a sequence can open inside a codon: the phase of its first exon pads it into frame, and the codons of the test and the residue numbers are counted from the first padded base. Both engines number the residues so for a change given in an exon, and vep-rs writes that notation for the shifted change.

On a transcript flagged `cds_end_NF`, with or without `cds_start_NF`, and on any other whose coding sequence does not end on its stop codon, vep-rs writes the notation that the whole codons of the sequence give. A frameshift is described from the first residue that differs, so that residue and the one in its place must both be whole codons of known bases. The changed frame is read through the coding sequence and any bases the transcript holds after it, and the new stop is counted on them. Where the frame meets none before the transcript's last base, the count is the unknown one, `fsTer?`, as both engines write it for a change given in an exon.

A change of whole codons must change a residue of the sequence's whole codons. Replaced residues always do. Deleted, duplicated or inserted residues do where the sequence after them does anything but repeat them in the same order.

Where no whole codon of the sequence differs, vep-rs leaves the HGVSp empty as Ensembl VEP does, unless Ensembl VEP's wider test fills it. A frameshifting deletion in the last whole codon of a sequence with no base after it leaves no whole codon for the residue in its place, and a frameshift whose new frame keeps every residue to the sequence's end changes none. Residues deleted, duplicated or inserted where the sequence ends, or where it only repeats them to its end, change none either: HGVS describes such a change at its last place in the protein, which can lie past that end. vep-rs reads no residue from a codon the sequence does not complete.

A change of the stop codon is read apart. On a transcript flagged `cds_end_NF` whose coding sequence is a whole number of codons, the last of them a stop codon, it means a shifted change that keeps every residue before the stop codon. vep-rs gives such a change the notation it writes for the same change given at the shifted position, beside Ensembl VEP's empty field.

On a transcript flagged `cds_end_NF` Ensembl VEP reads no insertion or deletion as overlapping the stop codon (`Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1373 and 1400). vep-rs follows that reading for the change given at the shifted position, so the notation it writes there can differ from the one the transcript would carry without the flag. The class does not take the notation of such a change on a row whose `FLAGS` field names `cds_end_NF`, so that row stays a difference in a comparison of outputs that carry that field. Where they do not, the class reads the row as one without a flag.

A shifted change that keeps every residue before the stop codon is written as [the page of the stop codon](hgvsp-stop-codon-protein-comparison.md) describes. Where a stop still stands at the stop's residue, vep-rs writes the stop as unchanged. Where residues follow the kept protein, it writes an extension named at the stop. For such a change given at the shifted position the two engines' values can differ, in the ways that page lists.

Two more kinds of row keep an empty HGVSp in vep-rs. The first is an insertion that the record writes inside an intron and that the shift leaves between the intron's last base and the exon's first base. On a transcript that has no intron of 13 bases or fewer, Ensembl VEP writes no HGVSp for such an insertion either. Where Ensembl VEP does write one, the row is a difference that the class does not take. The second is an allele written before the coding start whose shifted change is not wholly coding: both engines leave its HGVSp empty.

Where vep-rs writes the notation, nothing else in the row changes: the consequence terms and the empty position columns stay those of the allele as written, as in Ensembl VEP's row.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1657-1822, `hgvs_protein`: line 1711 shifts the allele 3'. Lines 1730-1738 return nothing unless the `coding` predicate holds (line 1731) and a translation start and end exist (lines 1732-1733). The start and end are read at the shifted position; the predicate is not.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 337-398, `_pre_consequence_predicates`: line 340 computes the predicates only when the allele has none, and line 394 keeps them. They are computed for the row's consequence terms, which are filled before the HGVS fields, so they exist before line 1711 shifts the allele.
- The same file, lines 454-571, `_bvfo_preds`: line 460 takes a shift only from an allele that was shifted when the predicates are computed, which `TranscriptVariation.pm`, lines 132-135, does only under `--shift_3prime 1`. The span of line 462 is therefore the written one. Lines 511-516 set `exon` from the exons that span overlaps, and lines 535-537 set `non_coding` when it overlaps none; `coding` is then never set. Nor is it for a span before the coding region, which lines 523-525 give `utr`.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 988-1040, `_create_intron_trees`: lines 1032-1035 set `_has_frameshift_intron` on a transcript where an intron's end lies at most 12 above its start, an intron of 13 bases or fewer.
- The same file, lines 848-870, `_overlapped_exons`: on such a transcript line 861 sets a stretch of 12, and lines 863-866 fetch the exons that overlap the span widened by it. An allele written in an intron within 12 bases of an exon then gets `exon` (lines 511-516) and `coding` (lines 540-557) in `_bvfo_preds`, and line 1731 of `hgvs_protein` passes.

**What correcting it would change.** With the coding test made on the shifted change, these rows would carry the notation Ensembl VEP writes when the record gives the change at the shifted position. The engines' rows differ in the HGVSp alone, so no other column would change.

**The class it produces.** `hgvsp_shifted_indel_written_span_gate` sets aside a protein notation in vep-rs's `HGVSp` against an empty `HGVSp` in Ensembl VEP's. The row is that of an insertion or deletion, annotated inside the transcript, with no CDS or protein position, and with a cDNA position only beside the term `5_prime_UTR_variant`. Both engines' rows hold the same HGVSc and the same non-zero `HGVS_OFFSET`, negative on a reverse-strand transcript. The HGVSc places the shifted change in the coding sequence: a deletion of coding bases alone, a duplication that ends on a coding base, or an insertion between two coding bases.

A row without a cDNA position is written in an intron, and the table gives what shows that the class holds for it.

| Row | What shows it |
| --- | --- |
| a cDNA position beside `5_prime_UTR_variant` | nothing more is asked |
| no cDNA position, and a notation that gives the first base of the exon the change is shifted into as coding base 2 or later | the row alone: the intron lies inside the coding region |
| no cDNA position, and a duplication that begins before the coding start | the row alone: the record is written before the coding start, or in an intron that the duplicated bases span |
| no cDNA position, and any other notation of coding bases | the exon spans of vep-rs's cache, which the comparison takes as `--json-cache`: the allele as the row's Location writes it lies wholly in the intron before the exon that holds the shifted change, whether that intron lies inside the coding region or in the 5' UTR; without the cache the row is not set aside |

A duplication whose first coordinate lies in an intron inside the coding region names the exon's first base. Any other notation of coding bases gives a lower bound for it: the first coding base of a deletion, or the coding base after a duplication or an insertion, less the offset's size, plus 1.

Three more kinds of row are not in the class, and a protein notation on one side alone is a difference on each. The first is an insertion that the shift leaves between an intron's last base and the exon's first base, whatever bases it inserts and however long the shift. Its HGVSc is a duplication that ends on the intron's last base, or an insertion between that base and the exon's first.

The second is a row on which vep-rs's notation describes a change at the transcript's own stop codon in any form but two. The class takes the stop written as unchanged and an extension named at the stop, each of which keeps every residue before the stop. It takes no notation that replaces or deletes the stop without an extension, names it as a flank of an insertion, or ends a duplicated, deleted or replaced range at it.

The third is a row whose `FLAGS` field names `cds_end_NF` in either engine's output and whose notation is one of those two forms at the stop codon. Every other notation that the class takes on a row without a flag, it takes on such a row, and a row flagged `cds_start_NF` alone is read as a row without a flag. The class reads the flags from the row's `FLAGS` field alone, so an output written without that field shows it none. A coding sequence that does not end on its stop codon shows in no field of a row that carries no flag, and the class asks nothing of it.

A notation on Ensembl VEP's side against none on vep-rs's and two notations that differ are differences as well. A row whose `SHIFT_LENGTH` field is positive was annotated at the shifted position, as under `--shift_3prime 1`, and is not in the class either.

The class is built into the format comparator, `tools/vep-formats-compare`, which holds Ensembl VEP's output against vep-rs's and applies the class row by row. The concordance figures that compare vep-rs with Ensembl VEP are scored on sets of consequence terms, so the HGVSp never enters them.

Eight worked examples follow. In them, angle brackets hold a notation described in words.

## 1. The last base of an intron of ENST00000265857, deleted

**In short.** On GRCh37, the record at chromosome 7, position 934969, has a two-base reference allele and a one-base alternate. It deletes the base at 934970, the last base of intron 8 of the forward-strand transcript ENST00000265857. The next two bases repeat that base, so the 3' rule carries the deletion two bases on: it becomes the deletion of coding base 897, the second base of exon 9. Under `--hgvs` each engine's row holds `splice_acceptor_variant`, no cDNA, CDS or protein position, `HGVSc=<the deletion of coding base 897>` and `HGVS_OFFSET=2`. The HGVSp is empty in Ensembl VEP's row and `<a frameshift at residue 300 with a new stop 3 residues on>` in vep-rs's.

A record at 934971 with a two-base reference allele and a one-base alternate gives the same deletion at the shifted position. Both engines annotate it `frameshift_variant,splice_region_variant` at `Protein_position` 299, with the HGVSp that vep-rs writes for the record in the intron. Coding base 897 is the third base of codon 299, and the notation names residue 300, the first residue the frameshift changes. In the written record the deleted base lies in no exon, so lines 511-516 leave `exon` at 0 and lines 535-537 set `non_coding`. The test of line 1731 fails, although lines 1732-1733 find residue 299 at the shifted position.

## 2. Thirty-one bases inserted 16 bases before an exon of ENST00000164305

**In short.** On GRCh37, the record at chromosome 15, position 55642880, has a one-base reference allele and a 32-base alternate. It inserts 31 bases between 55642880 and 55642881, 17 and 16 bases before exon 10 of the forward-strand transcript ENST00000164305. The inserted bases repeat the 31 reference bases after them, so the 3' rule carries the insertion 31 bases on, to stand between coding bases 1138 and 1139. Under `--hgvs` each engine's row holds `splice_polypyrimidine_tract_variant,intron_variant`, no cDNA, CDS or protein position, `HGVSc=<a duplication from 16 bases before coding base 1124 to coding base 1138>` and `HGVS_OFFSET=31`. The HGVSp is empty in Ensembl VEP's row and `<a frameshift at residue 380 with a new stop 44 residues on>` in vep-rs's.

For the insertion given at the shifted position, both engines write `frameshift_variant`, `Protein_position` 380 and the HGVSp that vep-rs writes for the record in the intron. The duplicated bases are the last 16 of the intron and the first 15 of exon 10. The notation's first coordinate therefore lies in the intron and names coding base 1124 as the exon's first base, which is the base the class reads: the exon opens later than coding base 2.

## 3. One deletion on two transcripts, one of them with an intron of one base

**In short.** On GRCh37, the record at chromosome 7, position 21934500, deletes the base at 21934501, the last base of an intron in two forward-strand transcripts of one gene. The 3' rule carries the deletion four bases on, into the exon that follows: it deletes coding base 12958 of ENST00000328843 and coding base 12937 of ENST00000409508. On ENST00000328843, which has an intron one base long, both engines write `HGVSp=<a frameshift at residue 4320 with a new stop 15 residues on>`. On ENST00000409508 the HGVSp is empty in Ensembl VEP's row and `<a frameshift at residue 4313 with a new stop 15 residues on>` in vep-rs's.

The record holds a two-base reference allele and a one-base alternate. On both rows each engine writes `splice_acceptor_variant`, the HGVSc of the deletion of that coding base and `HGVS_OFFSET=4`. The one-base intron of ENST00000328843 makes lines 1032-1035 set `_has_frameshift_intron` on it, and lines 861-866 widen the span by 12 bases: the deleted base, one base before the exon, gets `exon` and `coding`, and the test of line 1731 passes. In a run of Ensembl VEP with that stretch set to 0, the HGVSp on ENST00000328843 is empty as well. ENST00000409508 has no intron of 13 bases or fewer, and the test fails there as in example 1.

For the deletion given at the shifted position, both engines write `frameshift_variant` and those two notations, one on each transcript. Ensembl VEP's wider test is the subject of [another page](hgvsp-coding-test-widened-by-a-short-intron.md); this example shows only that the row on ENST00000409508 is of this class.

## 4. Two bases deleted on ENST00000397591, whose coding sequence does not end on its stop codon

**In short.** On GRCh37, the record at chromosome 15, position 39874391, has a three-base reference allele and a one-base alternate. It deletes the bases at 39874392 and 39874393, the last two bases of intron 2 of the forward-strand transcript ENST00000397591. The transcript carries the flag `cds_end_NF`, and no stop codon closes its coding sequence. The 3' rule moves the deletion two bases on, to coding bases 68 and 69 in exon 3. Both engines write the HGVSc of the deletion of those two coding bases. The HGVSp is empty in Ensembl VEP's row and `<a frameshift at residue 23 with a new stop 9 residues on>` in vep-rs's.

The transcript has no intron of 13 bases or fewer, so Ensembl VEP's test fails as in example 1. Residue 23, the residue in its place and the new stop are whole codons of the 94 the sequence holds, so vep-rs writes the notation. Given at the shifted position, the deletion gets the same HGVSp from both engines. The notation names no stop codon of the transcript, and the class takes the row.

Two records of GRCh38 show the rule's two edges, each a two-base reference allele and a one-base alternate on a transcript flagged `cds_end_NF`. The record at chromosome 5, position 138904347, deletes the last base of intron 4 of ENST00000518381, and the shift carries it one base on, to coding base 187 of 225. The frame meets no stop in the bases left, and vep-rs writes `<a frameshift at residue 63 with an unknown count to the new stop>` beside Ensembl VEP's empty field.

The record at chromosome 17, position 7221950, deletes the last base of intron 5 of ENST00000579886, and the shift carries it three bases on, to coding base 463 of 465, the first base of the last codon, which is no stop codon. The two bases left are no whole codon, and both engines leave the HGVSp empty.

## 5. A deletion that changes the stop codon of ENST00000995548

**In short.** On GRCh38, the record at chromosome 1, position 12030026, has a two-base reference allele and a one-base alternate. It deletes the base at 12030027, the last base of intron 7 of the forward-strand transcript ENST00000995548. The next two bases repeat that base, so the 3' rule carries the deletion two bases on, to coding base 816: the last base of codon 272, one base before the stop codon. Each engine's row holds `splice_acceptor_variant,NMD_transcript_variant`, `HGVSc=<the deletion of coding base 816>` and `HGVS_OFFSET=2`. The HGVSp is empty in Ensembl VEP's row and `<an extension at the stop of residue 273, a glutamate first, with a new stop 21 residues on>` in vep-rs's.

Without coding base 816, codon 272 takes the stop codon's first base as its last base and keeps its residue. The stop's residue becomes a glutamate, and the frame runs on through the 3' UTR to a stop codon 21 residues on. Every residue before the stop is kept, so the change is an extension, and the class takes the row in that form.

On ENST00000235332, a transcript of the same gene, the same deletion is that of coding base 847, and vep-rs writes `<the valine at residue 283 replaced by a stop>` beside Ensembl VEP's empty field. That notation is written at a residue before the transcript's stop codon, and the class takes the row as it takes the row of example 1.

A stop can also stay a stop. On GRCh38, the record at chromosome 1, position 109690269, with a two-base reference allele and a one-base alternate, deletes the base at 109690270, the last base of intron 5 of the forward-strand transcript ENST00000997889. The 3' rule carries the deletion one base on. Intron 5 splits that transcript's stop codon, coding bases 358 to 360: the exon's first base is coding base 359, the stop codon's second base. Without coding base 359 the codon reads a stop with the first base of the 3' UTR. vep-rs writes `<the stop at residue 120, unchanged>` beside Ensembl VEP's empty field, and the class takes the row in that form.

## 6. An insertion of 311 bases carried 310 bases on, into exon 4 of ENST00000356607

**In short.** On GRCh38, the record at chromosome 1, position 10598981, has a one-base reference allele and a 312-base alternate. It inserts 311 bases between 10598981 and 10598982, in intron 3 of the forward-strand transcript ENST00000356607, 257 and 256 bases before exon 4. The first 310 inserted bases repeat the 310 reference bases after them, so the 3' rule carries the insertion 310 bases on, to stand between coding bases 223 and 224. Each engine's row holds `intron_variant`, `HGVSc=<an insertion of 311 bases between coding bases 223 and 224>` and `HGVS_OFFSET=310`. The HGVSp is empty in Ensembl VEP's row and `<a frameshift at residue 75 with a new stop 5 residues on>` in vep-rs's.

The exon opens at coding base 170, and the notation alone cannot show it. Its bound is the coding base after the insertion, 224, less the offset of 310, plus 1, a number below 2. The class reads the exon spans of the transcript in vep-rs's cache instead. Coding base 224 lies 54 exon bases past the exon's first base, which is therefore coding base 170, and the row is set aside. In a comparison that names no cache it stays a difference.

## 7. One deletion on three transcripts of one gene, two of them flagged

**In short.** On GRCh37, a record at chromosome 7, position 780454, with a two-base reference allele and a one-base alternate, deletes the base at 780455. That base is the last base of an intron in three forward-strand transcripts of one gene, and the 3' rule carries the deletion one base on, into the exon that follows. On each of the three rows both engines write `splice_acceptor_variant`, the HGVSc of the deletion of the exon's first coding base and `HGVS_OFFSET=1`. Ensembl VEP's HGVSp is empty on all three.

ENST00000297440 carries neither flag, and vep-rs writes `<a frameshift at residue 261 with a new stop 64 residues on>` on its row, as in example 1. ENST00000440747 carries the flag `cds_start_NF`: its coding sequence opens at the third base of a codon, behind two bases of padding, and ends on its stop codon. The deleted base is coding base 185 there, padded base 187 in codon 63, and vep-rs writes `<a frameshift at residue 63 with a new stop 64 residues on>`. It is the same change of the same reading frame, numbered on the codons of the shorter sequence.

ENST00000437419 carries both flags, and its coding sequence holds the new stop: vep-rs writes `<a frameshift at residue 34 with a new stop 64 residues on>`. For the deletion given at the shifted position, both engines write `frameshift_variant` and, on each of the three transcripts, the notation above. The class takes all three rows.

## 8. Seven bases added in the 5' UTR of ENST00000409508

**In short.** On GRCh37, a made-up record at chromosome 7, position 21582861, adds seven bases in the 5' UTR of the forward-strand transcript ENST00000409508. They repeat the seven bases after them, the last two of the 5' UTR and coding bases 1 to 5, so the 3' rule carries them seven bases on. Lines 523-525 give the written span `utr`, and the test of line 1731 fails as in example 1.

| Record | Engine | Terms, HGVSc, `HGVS_OFFSET` | HGVSp |
| --- | --- | --- | --- |
| the made-up record | Ensembl VEP | `5_prime_UTR_variant` at cDNA position 29-30, `<a duplication from two bases before the coding start to coding base 5>`, 7 | empty |
| the made-up record | vep-rs | the same | `<a frameshift at residue 3 with a new stop 105 residues on>`: the changed codon 2 keeps its residue, and the class takes the row |
| the same seven bases given after coding base 5, at position 21582868 | both | `frameshift_variant`, the same HGVSc, no offset | that frameshift |
| genome-wide ClinVar, GRCh38, chromosome 3, position 37014474: one of four equal bases deleted, the last of them coding base 1 of ENST00000441265 | Ensembl VEP | `5_prime_UTR_variant` at cDNA position 817, `<the deletion of coding base 1>`, 3 | empty |
| the same ClinVar record | vep-rs | the same | the notation it writes for the deletion given at coding base 1, a start codon that stands again; the class takes the row |
