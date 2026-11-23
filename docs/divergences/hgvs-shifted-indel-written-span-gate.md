# The HGVSc of an insertion or deletion shifted into or over a transcript's edge

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 writes no HGVSc and no `HGVS_OFFSET` for an insertion or deletion whose 3' shift carries it into a transcript from wholly before its first base, or over its last base, on a row annotated inside the transcript. On such a row vep-rs writes the notation of the shifted span, with the `HGVS_OFFSET` of the shift.

**Why this is a defect.** The row describes the allele at the position the shift carried it to: its consequence, its `cDNA_position` and its `SHIFT_LENGTH` are the shifted span's. At that span the deleted bases are bases of the transcript, and the inserted bases fall between two of them, so a transcript-level notation exists, and Ensembl VEP itself writes it for the same allele given at the shifted position. An empty field says the row has no such description, which is false.

**What vep-rs does instead.** vep-rs writes the notation of the shifted span with the `HGVS_OFFSET` of the shift: the duplication or insertion of the inserted bases, or the deletion of the deleted ones. Where the shifted span still runs over an edge, the deletion is described over the bases inside the transcript (`n.1_16del`, `n.348del`).

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2689-2736, `_var2transcript_slice_coords`: lines 2708-2713 map the span as written (`unshifted_start` and `unshifted_end` once the variant was shifted) onto the transcript slice, and line 2729 returns nothing when both ends lie before the first base or both past the last. A span wholly before the transcript is refused here whatever the shift.
- The same file, lines 1446-1447 and 1472, `hgvs_transcript`: the shift's offset is read from the shift hash, and the method returns nothing when the written end, clamped onto the last base at line 2733, plus the offset passes the transcript's length. A deletion the shift carries over the last base is refused here.
- Lines 1477-1478 add the offset to the slice coordinates of the written span, for the rows that pass.

**What correcting it would change.** The rows would carry the notation and offset Ensembl VEP writes for the same allele given at the shifted position.

**The class it produces.** `hgvs_shifted_indel_written_span_gate`: an empty `HGVSc` and `HGVS_OFFSET` on Ensembl VEP's side against an insertion or deletion notation of the row's transcript and a non-zero offset on vep-rs's, on an insertion's or deletion's row annotated inside the transcript. Where the run prints a `SHIFT_LENGTH`, it is positive; the two fields of the HGVSIntronOffset plugin, which follow the notation, are read the same way. The class is applied once per row, built into the format comparator, as the oracle's population cases under `--hgvs` apply it; a notation on both sides is another class or a difference. The concordance figures score consequence sets and never see the field.

On a RefSeq model edited to the genome, the row Ensembl VEP leaves without a notation also lacks the `REFSEQ_OFFSET` that vep-rs writes. The class reads that field with the notation, as section 4 shows.

## 1. A `T` inserted two bases before the first base of ENST00000602620

**In short.** `21 16420340 . C CT` inserts a `T` after the `C` at 16420340, two bases before the first base of ENST00000602620, a forward-strand transcript opening with fifteen `T`s and then an `A`. Under `--hgvs` alone both engines write `upstream_gene_variant`; Ensembl VEP writes no HGVSc there, and vep-rs `HGVSc=ENST00000602620.2:n.15dup;HGVS_OFFSET=17`, a row of [`hgvs_added_bases_written_outside_the_transcript`](hgvs-added-bases-written-outside-the-transcript.md). Under `--hgvs --shift_3prime 1 --shift_length` both write `non_coding_transcript_exon_variant`, `cDNA_position` `15-16` and `SHIFT_LENGTH=17`; Ensembl VEP writes no HGVSc and no offset, vep-rs `HGVSc=ENST00000602620.2:n.15dup;HGVS_OFFSET=17`. The same allele given at the shifted position, `21 16420357 . T TT`, is `n.15dup` on both engines.

The written insertion maps to slice coordinates (-1, -2), both below 1, so line 2729 returns nothing. An insertion anchored on the base before the first base maps to (1, 0) and passes; its row is the class [`hgvs_shifted_insertion_clamped_end`](hgvs-shifted-insertion-clamped-end.md).

## 2. A deletion carried into ENST00000447861, and over the first base of ENST00000815447

**In short.** `21 13037979 . TTGCCTAACGGTTCTGCTCAGC T` deletes twenty-one bases, 13037980-13038000, that repeat the twenty-one after them; ENST00000447861 starts at 13038009. Under `--hgvs --shift_3prime 1 --shift_length` both engines write `non_coding_transcript_exon_variant`, `cDNA_position` `5-25` and `SHIFT_LENGTH=33`; Ensembl VEP writes no HGVSc, vep-rs `HGVSc=ENST00000447861.6:n.5_25del;HGVS_OFFSET=33`. The same allele given at the shifted position, `21 13038012 . TCTGCTCAGCTGCCTAACGGTT T`, is `n.5_25del` on both engines.

`21 8815408 . GCCGGCGCTGTGGGCCTCTCTGCGCCTGCGCCGGCGCTGTGGGCCTCTCTGCGCCTGCC G` deletes fifty-eight bases before ENST00000815447, whose first base is 8815477, and the shift moves them to 8815435-8815492, over that base: `cDNA_position` `?-16`, `SHIFT_LENGTH=26`, no HGVSc on Ensembl VEP, `n.1_16del` with `HGVS_OFFSET=26` on vep-rs. Written at the shifted position the deletion straddles the first base, and both engines clamp it to `n.1_16del`.

## 3. A deletion carried over the last base of ENST00000624130

**In short.** `21 9810835 . GAC G` deletes an `AC` of the run `ACACACACACACAC` at 9810836-9810849; ENST00000624130, 348 bases long, ends at 9810848. Under `--hgvs` alone both engines annotate the written span, `cDNA_position` `336-337`; under `--hgvs --shift_3prime 1 --shift_length` the shifted one, `348-?` with `SHIFT_LENGTH=12`. In both runs Ensembl VEP writes no HGVSc, where vep-rs writes `HGVSc=ENST00000624130.1:n.348del;HGVS_OFFSET=12`, the one deleted base inside the transcript. The same allele given at the shifted position, `21 9810847 . CAC C`, straddles the last base as written and is `n.348del` on both engines.

The written span maps to slice coordinates (336, 337), inside the transcript; the clamp leaves it, and line 1472 refuses the end 337 plus the offset 12, 349, against the length 348.

## 4. The last two bases of NM_001331018.1 deleted, a model edited to the genome

**In short.** `21 37223419 . ATA A` on GRCh38 deletes 21:37223420-37223421, the last two bases of NM_001331018.1, a reverse-strand RefSeq model of the committed `GRCh38-merged` corpus cache. Under `--merged --hgvs` both engines write `3_prime_UTR_variant`, `cDNA_position` `3457-3458` and `BAM_EDIT=OK`. Ensembl VEP writes no HGVSc, no `HGVS_OFFSET` and no `REFSEQ_OFFSET`; vep-rs writes `HGVSc=NM_001331018.1:c.*2124del`, `HGVS_OFFSET=-1` and `REFSEQ_OFFSET=1`. The offset belongs to the allele's cDNA start, 3457, whatever the notation: under `--shift_hgvs 0` both engines print `REFSEQ_OFFSET=1` for the same record.

On the transcript's strand the base that follows the deletion, 21:37223419, repeats the first deleted base, 21:37223421, so the shift carries the deletion one base on, over the model's last base. The written span ends at slice position 44689, the last of the slice, and line 1472 refuses that end plus the offset 1, as in section 3. vep-rs writes the one deleted base inside the transcript, `c.*2124`.

Ensembl VEP computes `REFSEQ_OFFSET` inside the routine that writes the notation. Line 1517 of `TranscriptVariationAllele.pm` calls `get_misalignment_offset` (lines 1594-1624), and line 1621 stores the value that `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1778, prints. A row that leaves the routine at line 1472 never reaches line 1517, so the column is lost with the notation.

The model carries two `_rna_edit` attributes: `647 646` adds one base, and `3458 3457` adds eighteen after the last base. The first starts before the allele's cDNA start and the second after it, so the summed length change is 1. [RefSeq and merged cache columns](../output/refseq-and-merged.md) states the column by the allele's cDNA start and the model's edits alone.

The class reads an empty `REFSEQ_OFFSET` on Ensembl VEP's side against a non-zero whole number on vep-rs's, on a row that prints `BAM_EDIT=OK` in both engines and meets the class's conditions on `HGVSc` and `HGVS_OFFSET`. It cannot tell from the row that the value is the sum of the model's edits before the allele: the cache holds them. A `REFSEQ_OFFSET` that Ensembl VEP alone writes stays a difference.
