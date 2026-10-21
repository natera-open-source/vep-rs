# The start-codon test of a boundary deletion under `--shift_3prime 1`

A difference in the terms, and so in the HGVSp, of a deletion that starts in an intron and is shifted into or along the exon under `--shift_3prime 1`; nothing is set aside, since no measured dataset runs the flag.

Under `--shift_3prime 1`, Ensembl VEP tests whether a shifted allele overlaps the start codon from the allele's unshifted cDNA coordinates plus the shift (`Utils/VariationEffect.pm` lines 977-981). A deletion that starts in an intron has no unshifted cDNA start, which the sum reads as zero, so the test holds for any such deletion that ends past the start codon.

`start_retained_variant` is then the negation of a second test (lines 947-963), which fails where the shifted span still has no cDNA start (line 1048) or lies downstream of an intact 5' UTR and `ATG` (line 1065). On a transcript with no 5' UTR the second test holds instead, since the deletion shortens the sequence (line 1069). That row is `start_lost` (line 870), with the HGVSp of a lost start, `?` for the new peptide (`TranscriptVariationAllele.pm` lines 2161-2164).

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

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--shift_3prime 1 --hgvs`; vep-rs with the same flags over the same cache and FASTA.
