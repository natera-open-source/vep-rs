# The start-codon test of a boundary deletion under `--shift_3prime 1`

A difference in the terms, and so in the HGVSp, of a deletion that starts in an intron and is shifted into or along the exon under `--shift_3prime 1`. Ensembl VEP tests the start codon from the unshifted cDNA coordinates plus the shift, so its row can gain `start_retained_variant`, or `start_lost` on a transcript with no 5' UTR. vep-rs tests the shifted span's own cDNA coordinates, so its terms and HGVSp are the shifted span's. The format comparator sets aside the keys a manifest documents under the class `shift_start_codon_test`, as the oracle's population cases under `--shift_3prime 1` do for the deletions of this shape; nothing is set aside in the concordance figures, since no measured dataset runs the flag.

Under `--shift_3prime 1`, Ensembl VEP tests whether a shifted allele overlaps the start codon from the allele's unshifted cDNA coordinates plus the shift (`Utils/VariationEffect.pm` lines 977-981).

A deletion that starts in an intron has no unshifted cDNA start, which the sum reads as zero. On a transcript not flagged `cds_start_NF` (line 975) the test therefore holds for such a deletion that ends on an exon base past the start codon where the shift is no greater than the start codon's last cDNA position. One that ends in an intron or past the transcript's last base has no unshifted cDNA end either. For it the test holds only where the shift itself is one of the start codon's cDNA positions, as on a transcript with no 5' UTR and a shift of one to three bases.

A deletion that starts before the transcript's first exon and ends on an exon base past the start codon has no unshifted cDNA start either. The test holds for it on a transcript not flagged `cds_start_NF` (line 975), where the shift is no greater than the start codon's last cDNA position and leaves the deletion's start before the first exon and its end inside the transcript. The shifted span then still has no cDNA start (line 1048), and Ensembl VEP's row carries `start_retained_variant` although the deletion removes the start codon.

Where the shift leaves the deletion's start inside the intron, short of the acceptor site, the shifted span still has no cDNA start (line 1048). A deletion that starts in an intron and ends on an exon base past the start codon therefore gains `start_retained_variant` on Ensembl VEP's row in the same way, on three conditions. The transcript is not flagged `cds_start_NF` (line 975), the shifted span overlaps the coding region, so that the predicate `coding` is set, and the shift is no greater than the start codon's last cDNA position.

One that ends in a later intron or past the transcript's last base gains it only where the shift itself is one of the start codon's cDNA positions, as on a transcript with no 5' UTR and a shift of one to three bases.

`start_retained_variant` is then the negation of a second test (lines 947-963), which fails where the shifted span still has no cDNA start (line 1048) or lies downstream of an intact 5' UTR and `ATG` (line 1065). On a transcript with no 5' UTR the second test holds instead, since the deletion shortens the sequence (line 1069). That row is `start_lost` (line 870), with the HGVSp of a lost start, `?` for the new peptide (`TranscriptVariationAllele.pm` lines 2161-2164).

A deletion of the whole transcript is `transcript_ablation` without the flag. Under the flag it stays `transcript_ablation`, with no start term, while the shifted span still covers the whole transcript. The shifted span does so where the deletion starts before the transcript's first exon and the shift does not carry its start past the transcript's first base.

Where the deletion runs from the transcript's first base to its last, the shift carries its start into the transcript and its end past the last base. On a transcript not flagged `cds_start_NF` (line 975), where 1 plus the shift is no greater than the start codon's last cDNA position, the first test then holds on the unshifted coordinates plus the shift. The second fails because the shifted span has no cDNA end (line 1048), and Ensembl VEP's row then carries `start_retained_variant` beside the terms of the shifted span.

vep-rs tests the shifted span's own cDNA coordinates against the start codon, so the terms and the HGVSp are the shifted span's.

The records below sit at the acceptor of exon 6 of CLIC6, where the deletion shifts one base into the exon and becomes a frameshift, and at an exon end of ENSG00000154719, where it straddles the boundary before and after the shift:

```
21	34716319	del_exon_boundary_coding_ok	AGATT	A	.	.	.
21	25585752	del_exon_end_rev_ok	TACCT	T	.	.	.
```

Under `--shift_3prime 1 --hgvs`, the two engines write for the first record:

| Engine | Feature | Consequence | HGVSp |
| --- | --- | --- | --- |
| Ensembl VEP | ENST00000349499 | `frameshift_variant,splice_region_variant,start_retained_variant` | `p.Ile634TrpfsTer14` |
| vep-rs | ENST00000349499 | `frameshift_variant,splice_region_variant` | `p.Ile634TrpfsTer14` |
| Ensembl VEP | ENST00000360731 | `frameshift_variant,start_lost,splice_region_variant` | `p.Ile652_?653` |
| vep-rs | ENST00000360731 | `frameshift_variant,splice_region_variant` | `p.Ile652TrpfsTer14` |

ENST00000360731 is the transcript with no 5' UTR. For the second record Ensembl VEP writes `splice_acceptor_variant,start_retained_variant` on the fifteen transcripts whose coding sequence the deletion reaches, `stop_retained_variant` added on the two whose stop codon it reaches; vep-rs writes `splice_acceptor_variant,coding_sequence_variant` and `splice_acceptor_variant,stop_retained_variant`. The rows of the transcripts the deletion meets in a 3' UTR or an intron agree, as do positions, codons and the other columns, and the rows of both records without the flag.

To reproduce: Ensembl VEP 116.3 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--shift_3prime 1 --hgvs`; vep-rs with the same flags over the same cache and FASTA.
