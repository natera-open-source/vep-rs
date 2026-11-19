# The HGVSc of an insertion shifted into a transcript's first exon under `--shift_3prime 1`

**The defect in one sentence.** Under `--hgvs --shift_3prime 1`, Ensembl VEP 116.2 writes an insertion anchored on the base before a transcript's first base, shifted k bases into exon 1, as a substitution or a one-base delins at n.k+1. A `T` inserted after fourteen `T`s is written `n.15A>T`, where the shifted insertion is the duplication `n.14dup`, the notation vep-rs writes.

**Why this is a defect.** An insertion has no reference base of its own: its notation names the two bases it falls between, or the bases it duplicates. At the position the shift carried it to, the inserted `T` repeats the base before it, n.14, so the notation the sequence admits is `n.14dup`. `n.15A>T` describes a change of the base n.15 into `T`, a substitution the record does not carry, and `n.15delinsTT` a replacement of that base by two.

**What vep-rs does instead.** vep-rs writes the notation of the span the row was annotated at: the duplication of the bases the shift passed over, or an insertion between n.k and n.k+1 where they are not a repeat, with `HGVS_OFFSET` k.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2689-2736, `_var2transcript_slice_coords`: the unshifted insertion, outside the transcript by one base on its 5' side, maps to the slice coordinates (1, 0) on either strand. Line 2729 returns nothing only when both are below 1, and line 2733 clamps the end to 1, so the insertion becomes the one-base span (1, 1).
- The same file, lines 1474-1483, `hgvs_transcript`: the shift k is added to both ends, and `hgvs_variant_notation` (`Utils/Sequence.pm`, lines 493-619) reads the base at n.k+1 as the reference of a substitution or a delins.
- The same file, lines 2231-2242 and 2271-2278, `_clip_alleles`: where the inserted sequence ends with that base, the base is stripped and the remainder is written as an insertion between n.k and n.k+1.

**What correcting it would change.** The row would carry the duplication Ensembl VEP itself writes for the same allele given at the shifted position.

**The class it produces.** `hgvs_shifted_insertion_clamped_end`: a differing `HGVSc` on an insertion's row whose `SHIFT_LENGTH` is positive, whose `cDNA_position` starts at that length and whose `HGVS_OFFSET` is set. Ensembl VEP's notation sits at one position P, a substitution, a delins, or an insertion of the preceding bases between P-1 and P; vep-rs's is a duplication ending at P-1 or an insertion between P-1 and P, of as many bases as the allele. The class is applied once per row, built into the format comparator, as the oracle's population cases under `--shift_3prime 1` apply it; an empty `HGVSc` on either side stays a difference. The concordance figures score consequence sets and never see the field.

## 1. A `T` inserted before the first base of ENST00000617755

**In short.** `21 25928753 . C CT` inserts a `T` after the `C` that precedes the first base of ENST00000617755, a forward-strand transcript opening with fourteen `T`s and then an `A`. Under `--hgvs` alone both engines write `upstream_gene_variant` with `DISTANCE=0`; Ensembl VEP writes no HGVSc there, and vep-rs `HGVSc=ENST00000617755.1:n.14dup;HGVS_OFFSET=14`, a row of [`hgvs_added_bases_written_outside_the_transcript`](hgvs-added-bases-written-outside-the-transcript.md). Under `--hgvs --shift_3prime 1 --shift_length` both write `non_coding_transcript_exon_variant`, `cDNA_position` `14-15` and `SHIFT_LENGTH=14`; Ensembl VEP adds `HGVSc=ENST00000617755.1:n.15A>T;HGVS_OFFSET=14`, vep-rs `HGVSc=ENST00000617755.1:n.14dup;HGVS_OFFSET=14`. The same allele given at the shifted position, `21 25928767 . T TT`, is `n.14dup` on both engines; the two-base insertion `C CTT` at the first position is `n.15delinsTT` against `n.13_14dup`.
