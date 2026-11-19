# The HGVSc of a deletion straddling a transcript's first base under the 3' shift

**The defect in one sentence.** Under `--hgvs`, when a deletion straddling a transcript's first base shifts 3' by k bases, Ensembl VEP 116.2 starts its notation at transcript base k+1, not at the first transcript base the shifted deletion removes, and so describes fewer bases than the row deletes. A twenty-one-base deletion shifted thirty-three bases into ENST00000785034 is written `n.34_43del`, ten bases, where the shifted span is `n.23_43del`, the notation vep-rs writes.

**Why this is a defect.** A deletion's notation names the first and the last deleted base of the transcript. The shifted span of that deletion runs from n.23 to n.43, twenty-one bases, and Ensembl VEP writes `n.23_43del` for the same allele given at the shifted position. `n.34_43del` describes the loss of ten bases, a change the record does not carry.

**What vep-rs does instead.** vep-rs describes the span the shift moved the deletion to, clamped onto the transcript where that span still runs over the first base: `n.23_43del` for the deletion above, `n.1_123del` for one whose shifted span starts before n.1.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2708-2713, `_var2transcript_slice_coords`: the span as written (`unshifted_start` and `unshifted_end` once the variant was shifted) is mapped onto the transcript slice, and lines 2732-2733 clamp a start before the first base onto position 1.
- The same file, lines 1477-1478, `hgvs_transcript`: the shift's offset k is added to both clamped coordinates, so the start moves from 1 to k+1 while the deleted bases before the first base, which the shift carried inside, are left out.

**What correcting it would change.** The notation would name the whole span the shifted deletion covers inside the transcript, as Ensembl VEP's does for the same allele given at the shifted position.

**The class it produces.** `hgvs_shifted_deletion_clamped_span`: two `del` notations of one transcript and one numbering that end at the same cDNA coordinate, Ensembl VEP's starting later than vep-rs's, on a deletion's row both engines give a non-zero `HGVS_OFFSET` and, where printed, a positive `SHIFT_LENGTH`. The class is built into the format comparator and applied once per row; an empty `HGVSc` on either side or another end stays a difference. The concordance figures score consequence sets and never see the field.

A run that prints `SHIFT_LENGTH` places the row's cDNA position at the shifted span, and there the class also holds vep-rs's notation to that position. Its `n.` coordinates must equal the position's, and a `c.` span inside one region (the 5' UTR, the coding sequence or the 3' UTR) must be as long as it, so a span off by one base stays a difference. A run without the field places the row at the written span, and the class makes no such test.

On a row of the class the HGVSIntronOffset plugin's `HGVS_IntronStartOffset` is read the same way: each engine writes the intron offset of its own notation's start coordinate, so the two values differ as the notations do, while `HGVS_IntronEndOffset` agrees.

## 1. A twenty-one-base deletion over the first base of ENST00000785034

**In short.** `21 13037979 . TTGCCTAACGGTTCTGCTCAGC T` deletes 13037980-13038000 on GRCh38, twenty-one bases that repeat the twenty-one after them, and ENST00000785034 starts at 13037991, inside the deleted span. The shift carries the deletion to 13038013-13038033, thirty-three bases on. Under `--hgvs`, with or without `--shift_3prime 1`, both engines write `non_coding_transcript_exon_variant` and `HGVS_OFFSET=33`; Ensembl VEP writes `HGVSc=ENST00000785034.1:n.34_43del`, vep-rs `HGVSc=ENST00000785034.1:n.23_43del`. The same allele given at the shifted position, `21 13038012 . TCTGCTCAGCTGCCTAACGGTT T`, is `n.23_43del` on both engines. On ENST00000785037, which starts three bases earlier, the row is `n.34_46del` against `n.26_46del`.

The written span maps to slice coordinates (-10, 10): eleven deleted bases before the first base and ten inside, and the eleventh base of the transcript follows. Line 2732 clamps the start to 1, lines 1477-1478 add 33 to both ends, and `hgvs_variant_notation` reads the ten bases 34-43 as the deletion.

The reproduction `scripts/divergences/reproductions/hgvs_shifted_deletion_clamped_span.toml` runs this record on both engines on GRCh38 under `--hgvs` and records the consequence, `HGVSc` and `HGVS_OFFSET` of each engine's row on ENST00000785034. For vep-rs it names the JSON cache of the committed corpus `tests/golden/116/GRCh38-release116`.

## 2. A deletion whose shifted span still starts before the first base

**In short.** A gnomAD v4.1 record on GRCh38 deletes 242 bases of a tandem repeat, 21:46193278-46193519, over the first base of ENST00000474319, a reverse-strand transcript whose first base is 46193370, so the written span covers n.1 to n.93. The shift moves it thirty bases towards lower coordinates, to 46193248-46193489, over n.1 to n.123. Under `--hgvs` both engines write `HGVS_OFFSET=-30`; Ensembl VEP writes `n.31_123del`, vep-rs `n.1_123del`, the one hundred and twenty-three bases of the transcript the shifted deletion covers. The reverse strand changes nothing in the mechanism: the clamp lands on position 1 and the offset moves it to 31.
