# A record whose alleles shift by different lengths under `--shift_3prime 1`

A difference in the terms, positions, codons and HGVSp of a record whose alleles shift by different lengths under `--shift_3prime 1`: a multi-allelic insertion or deletion, or the reference row of `--process_ref_homs` beside a shifted alternate. Nothing is set aside, since no measured dataset runs the flag.

Under `--shift_3prime 1`, Ensembl VEP shifts each insertion or deletion allele along its transcript's strand and reads every consequence predicate at the allele's shifted span (`TranscriptVariationAllele.pm` lines 109-150 and 411-466, `BaseVariationFeatureOverlapAllele.pm` lines 269-282). The alleles of one record on one transcript share one transcript-variation object, and two of its caches are filled once for all of them. `Utils/VariationEffect.pm` lines 806-807 name the design: the reference allele's shift is not precomputed, for variants with several alternate alleles of potentially different shift lengths.

The span predicates say whether the allele lies within the transcript, in an exon, in an intron or in the coding sequence. The first allele evaluated stores them at its shifted span (`BaseVariationFeatureOverlapAllele.pm` line 363, with the offset of line 460), and they choose each later allele's candidate terms, which are then tested at that allele's own span. The cDNA, CDS and protein coordinates are recomputed by every allele in turn before output, each with its own shift (`OutputFactory.pm` lines 306-323, `TranscriptVariationAllele.pm` lines 502-528). The coding predicates and the output read them back cached, whatever offset is passed (`BaseTranscriptVariation.pm` line 142, `Utils/VariationEffect.pm` line 1551), so a row reads coordinates another allele of the record left there.

So a second allele whose shifted span lies in an exon, after a first allele's lay in an intron, is tested against the intron's candidate terms. Its exonic terms, positions, amino acids and codons are not computed, and its HGVSp is absent. Its `Allele` stays as written where [a coding row](../output/vep-and-tab.md#allele-on-a-coding-row-after-a-shift) prints it rotated; the HGVS notations, which shift on their own, agree. With the alleles in the other order the exonic allele stores the predicates, but its coding predicates read the intronic allele's empty coordinates, so its row loses the coding term and the positions; the intronic allele tests exonic terms at its own span.

Under `--process_ref_homs` the reference allele comes first (`BaseVariationFeatureOverlap.pm` line 370) and is never shifted itself (`TranscriptVariationAllele.pm` line 116), so it stores the predicates at the unshifted span. Its row keeps that span's non-coding terms, and each alternate's candidate terms are the unshifted span's, tested at the shifted span. A deletion shifted from two bases upstream of a transcript into its 5' UTR therefore satisfies neither the upstream nor the downstream test and prints `intergenic_variant`. A coding reference row agrees with vep-rs's, since both engines read its positions and codons from the coordinates the last alternate left.

vep-rs reads each allele's predicates and coordinates at its own shifted span, so a row depends neither on the record's other alleles nor on `--process_ref_homs`. The reference row of `--process_ref_homs` is read at the span the record's last alternate shifts to, as [`SHIFT_LENGTH`](../output/vep-and-tab.md#hgvs_offset-and-shift_length) describes. Its positions, codons and coding terms therefore agree with Ensembl VEP's, and its non-coding terms are the shifted span's.

The records below are two insertions after a run of bases, on the forward-strand HUNK and the reverse-strand NRIP1, and a deletion two bases before the start of most USP25 transcripts:

```
21	31974714	multi_two_motifs_fwd	A	AG,AGT	.	.	.
21	15064747	multi_two_motifs_rev	C	CC,CTC	.	.	.
21	15729922	transcript_edge_hp_del	CG	C	.	.	.
```

In the first record the `G` allele shifts four bases along a run of `G` into intron 7, and the `GT` allele one base, staying in exon 7. Under `--shift_3prime 1 --shift_length --hgvs`, the two engines write on ENST00000270112:

| Engine | Allele | Consequence | Positions, Amino_acids and HGVSp |
| --- | --- | --- | --- |
| Ensembl VEP | `G` | `splice_donor_variant` | none |
| vep-rs | `G` | `splice_donor_variant` | none |
| Ensembl VEP | `GT` | `splice_region_variant` | none |
| vep-rs | `TG` | `frameshift_variant,splice_region_variant` | `1826-1827`, `1171-1172`, `391`, `G/VX`, `p.Gly391ValfsTer19` |

Both engines print `c.1171_1172insTG` with `HGVS_OFFSET=1` and `SHIFT_LENGTH=1` on the `GT` rows. Five of the record's ten rows differ, every one the second allele's. Written as `A AGT,AG`, nine of ten differ. Ensembl VEP's `GT` rows then print `splice_region_variant`, the positions `-` and an HGVSp, and its `G` rows `splice_donor_variant` with a Protein_position (`391` on ENST00000270112); vep-rs prints the same rows as before.

In the second record the `C` allele shifts four bases along a run of `C` into the donor site of exon 1, and the `TC` allele one base, staying in the exon's 5' UTR. The `C` rows agree on every transcript. On ENST00000318948 the `TC` row is `splice_region_variant` with no cDNA position in Ensembl VEP and `splice_region_variant,5_prime_UTR_variant` at `254-255` in vep-rs. Eight of the record's thirty rows differ, the `TC` rows of the transcripts whose first exon holds its span.

The third record's deletion shifts three bases, from two bases upstream of ENST00000285679 to the second base of its 5' UTR. Under `--process_ref_homs --shift_3prime 1 --shift_length` the two engines write:

| Engine | Row | Consequence | Extra |
| --- | --- | --- | --- |
| Ensembl VEP | reference | `upstream_gene_variant` | `DISTANCE=2;SHIFT_LENGTH=0` |
| Ensembl VEP | deletion | `intergenic_variant` | `SHIFT_LENGTH=3` |
| vep-rs | reference | `5_prime_UTR_variant`, cDNA position `2` | `SHIFT_LENGTH=0` |
| vep-rs | deletion | `5_prime_UTR_variant`, cDNA position `2` | `SHIFT_LENGTH=3` |

Without `--process_ref_homs` both engines write the deletion's vep-rs row. Over the record's 36 transcripts Ensembl VEP's reference row is `upstream_gene_variant` on 35, and its deletion row `intergenic_variant` on the 20 the shift carries the deletion into. vep-rs writes the deletion's rows as without the flag, and the reference row with the same terms.

Two more reference rows are of the kind. For `21 10521699 GC G` Ensembl VEP writes `splice_donor_region_variant,intron_variant` on ENST00000427445 where vep-rs writes `splice_region_variant,intron_variant`; for `21 15064746 CC C` it writes `splice_region_variant` on ENST00000318948 where vep-rs writes `splice_donor_variant`. The coding reference row of `21 10561150 AT A` agrees, `synonymous_variant` at `649`, `297` and `99` with the codons `ttT/ttT`.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--shift_3prime 1 --shift_length --hgvs`, and `--process_ref_homs` for the third record; vep-rs with the same flags over the same cache and FASTA.
