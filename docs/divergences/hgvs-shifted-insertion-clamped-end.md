# The HGVSc of an insertion shifted into a transcript's first exon under `--shift_3prime 1`

**The defect in one sentence.** Under `--hgvs --shift_3prime 1`, Ensembl VEP 116.3 writes an insertion anchored on the base before a transcript's first base, shifted k bases into exon 1, as a substitution or a one-base delins at n.k+1. A `T` inserted after fourteen `T`s is written `n.15A>T`, where the shifted insertion is the duplication `n.14dup`, the notation vep-rs writes.

**Why this is a defect.** An insertion has no reference base of its own: its notation names the two bases it falls between, or the bases it duplicates. At the position the shift carried it to, the inserted `T` repeats the base before it, n.14, so the notation the sequence admits is `n.14dup`. `n.15A>T` describes a change of the base n.15 into `T`, a substitution the record does not carry, and `n.15delinsTT` a replacement of that base by two.

**What vep-rs does instead.** vep-rs writes the notation of the span the row was annotated at: the duplication of the bases the shift passed over, or an insertion between n.k and n.k+1 where they are not a repeat, with `HGVS_OFFSET` k.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2689-2736, `_var2transcript_slice_coords`: the unshifted insertion, outside the transcript by one base on its 5' side, maps to the slice coordinates (1, 0) on either strand. Line 2729 returns nothing only when both are below 1, and line 2733 clamps the end to 1, so the insertion becomes the one-base span (1, 1).
- The same file, lines 1474-1483, `hgvs_transcript`: the shift k is added to both ends, and `hgvs_variant_notation` (`Utils/Sequence.pm`, lines 493-619) reads the base at n.k+1 as the reference of a substitution or a delins.
- The same file, lines 2231-2242 and 2271-2278, `_clip_alleles`: where the inserted sequence ends with that base, the base is stripped and the remainder is written as an insertion between n.k and n.k+1.

**What correcting it would change.** The row would carry the duplication Ensembl VEP itself writes for the same allele given at the shifted position.

**The class it produces.** `hgvs_shifted_insertion_clamped_end`: a differing `HGVSc` on an insertion's row whose `SHIFT_LENGTH` is positive, whose `cDNA_position` starts at that length and whose `HGVS_OFFSET` is set. Ensembl VEP's notation sits at one position P, a substitution, a delins, or an insertion of the preceding bases between P-1 and P; vep-rs's is a duplication ending at P-1 or an insertion between P-1 and P, of as many bases as the allele. The class is applied once per row, built into the format comparator, as the oracle's population cases under `--shift_3prime 1` apply it; an empty `HGVSc` on either side stays a difference. The concordance figures score consequence sets and never see the field.

On a model that carries an `_rna_edit` attribute, Ensembl VEP respells the alleles of the clamped notation. The substitution is printed with `-` as its reference base, and the insertion names the whole allele. Section 2 shows both forms and says where the class reads them.

## 1. A `T` inserted before the first base of ENST00000617755

**In short.** `21 25928753 . C CT` inserts a `T` after the `C` that precedes the first base of ENST00000617755, a forward-strand transcript opening with fourteen `T`s and then an `A`. Under `--hgvs` alone both engines write `upstream_gene_variant` with `DISTANCE=0`; Ensembl VEP writes no HGVSc there, and vep-rs `HGVSc=ENST00000617755.1:n.14dup;HGVS_OFFSET=14`, a row of [`hgvs_added_bases_written_outside_the_transcript`](hgvs-added-bases-written-outside-the-transcript.md). Under `--hgvs --shift_3prime 1 --shift_length` both write `non_coding_transcript_exon_variant`, `cDNA_position` `14-15` and `SHIFT_LENGTH=14`; Ensembl VEP adds `HGVSc=ENST00000617755.1:n.15A>T;HGVS_OFFSET=14`, vep-rs `HGVSc=ENST00000617755.1:n.14dup;HGVS_OFFSET=14`. The same allele given at the shifted position, `21 25928767 . T TT`, is `n.14dup` on both engines; the two-base insertion `C CTT` at the first position is `n.15delinsTT` against `n.13_14dup`.

## 2. A copy of the first base inserted before NR_135516.1, a model that carries an edit

**In short.** `21 26393634 . C CA` on GRCh38 inserts an `A` after the `C` that precedes the first base of NR_135516.1, a forward-strand RefSeq model that opens with `ACC` and carries one `_rna_edit` attribute. Under `--refseq --hgvs --shift_3prime 1 --shift_length` both engines write `non_coding_transcript_exon_variant`, `cDNA_position` `1-2`, `HGVS_OFFSET=1`, `SHIFT_LENGTH=1` and `BAM_EDIT=OK`. Ensembl VEP writes `HGVSc=NR_135516.1:n.2->A`, vep-rs `HGVSc=NR_135516.1:n.1dup`. The inserted `A` repeats n.1 at the place the shift carried it to, and a substitution with `-` for the base it replaces is no HGVS notation.

`21 18759812 . G GCG` inserts two bases before the first base of NR_109925.1, a reverse-strand model with one `_rna_edit` attribute, and the shift carries them three bases in. Both engines write `cDNA_position` `3-4`, `HGVS_OFFSET=-3` and `SHIFT_LENGTH=3`. Ensembl VEP writes `HGVSc=NR_109925.1:n.3_4insCG`, an `ins` of the whole allele, and vep-rs `HGVSc=NR_109925.1:n.2_3dup`: rotated by the shift, the inserted bases repeat n.2 and n.3.

The steps on the two rows, by line of `TranscriptVariationAllele.pm`:

| Step | Lines | `C CA` on NR_135516.1 | `G GCG` on NR_109925.1 |
| --- | --- | --- | --- |
| The insertion as written maps to the slice coordinates (1, 0), and the end is clamped to 1 | 2712, 2729, 2733 | (1, 1) | the same |
| The shift is added, and the notation is typed against the base read there | 1446-1447, 1474-1483 | Shift 1: a substitution at 2, the base n.2 against the inserted `A` | Shift 3: a `delins` at 4, one base against the two inserted bases, rotated by the shift |
| `_clip_alleles` takes the base the two alleles share | 1492, 2229-2237, 2283 | Nothing: the two bases differ | The last base: an `ins` of one base, with start 4 and end 3 |
| On a model with an `_rna_edit` attribute, both alleles are written over from `feature_seq` | 1508-1514 | The reference becomes `-`, the row's reference allele of an insertion | The alternate becomes the row's whole allele, `CG` |
| The positions are numbered and ordered, and the notation printed | 1524-1525, 1535-1538; `Utils/Sequence.pm`, 654-656 and 670-673 | `n.2->A` | `n.3_4insCG` |

On a model without the attribute lines 1508-1514 write nothing over, and the notation keeps the form of section 1. The three-base insertion `C CACC` at the first record's position passes the same lines and is `n.4delinsACC` against `n.1_3dup`, the form of section 1, since a `delins` prints the alternate alone.

Lines 1508-1514 are the step that [`refseq_edited_model_unminimised_allele_hgvsc`](refseq-edited-model-unminimised-allele-hgvsc.md) describes. Here they change only how the alleles are spelt. The position and the kind come from the lines this page cites, and the difference would stand without the respelling.

| What the class reads | On the two rows |
| --- | --- |
| The substitution printed with `-` | Ensembl VEP's `n.2->A` is read as a substitution at P. Every other condition of the class stands: a printed positive `SHIFT_LENGTH`, a `cDNA_position` that starts at it, and on vep-rs's side a `dup` ending at P-1 or an `ins` between P-1 and P. |
| The `ins` of the whole allele | On a row that prints `BAM_EDIT=OK`, Ensembl VEP's `ins` between P-1 and P names the row's `Allele` on the transcript's strand, against vep-rs's `dup` of as many bases ending at P-1. |
| What it cannot tell | For the substitution, that the model carries an edit. A row of a model that carries an edit and prints no `BAM_EDIT` is not read for the `ins`. |
| What stays a difference | Either pair with the sides exchanged, a whole-allele `ins` on a row without `BAM_EDIT=OK`, and both forms in a run that prints no `SHIFT_LENGTH`. |
