# A symbolic record with two or more ALTs under `--minimal`

A difference in the span, the marks and the rows written under `--minimal` for a structural record with more than one symbolic ALT; nothing is set aside, since no measured dataset runs the flag.

Under `--minimal`, Ensembl VEP's parser minimises every feature it builds (`Parser.pm` line 852, `minimise_alleles` at lines 964-1019). A structural record's feature carries its ALTs joined by `/` as its allele string (`Parser/VCF.pm` lines 469 and 572), not a reference and its alternates. A record with one symbolic ALT has no `/` in that string and is left alone.

A record with two is read as a reference and one alternate and trimmed as two sequences (`trim_sequences`, `Utils/Sequence.pm` lines 965-1038). `<DEL>/<DUP>` loses the `<D` the two names share at the front and the `>` at the back, becoming `EL/UP` over a span that starts two bases later and ends one earlier, and the feature is marked minimised. vep-rs minimises sequence records alone and annotates a structural record as written.

The records below overlap exon 6 of ENST00000427445 and an intron of two other transcripts of the gene; `multi_delstar` is a record of the `GRCh38-vcf-minimal` corpus:

```
21	10542400	lone_del	N	<DEL>	.	.	SVTYPE=DEL;END=10542450
21	10542400	del_dup	N	<DEL>,<DUP>	.	.	SVTYPE=DEL;END=10542450
21	10542414	multi_delstar	A	AG,<DEL:*>	.	.	.
```

Both engines write `lone_del` alike, nine rows at `21:10542401-10542450` with the Allele `deletion`, and `multi_delstar` alike, nine rows of the `G` allele at `21:10542414-10542415` with `UPLOADED_ALLELE` `A/AG/<DEL:*>` and no `MINIMISED` mark, since the record's alleles `-/G/<DEL:*>` trim nothing against each other. For `del_dup` both write 18 rows, two alleles on nine transcripts, with the Allele `copy_number_variation`; the span, the positions and the marks differ. On ENST00000427445 the two engines write:

| Engine | Location | cDNA_position | Extra |
| --- | --- | --- | --- |
| Ensembl VEP | `21:10542403-10542449` | `426-?` | `MINIMISED=1;OverlapBP=47` |
| vep-rs | `21:10542401-10542450` | `424-?` | `OverlapBP=50` |

`CDS_position` and `Protein_position` move with the cDNA position (`74-?` and `25-?` against `72-?` and `24-?`), `OverlapPC` is `0.10` against `0.11` on ENST00000612957, and the VCF `CSQ` entries carry the same positions and mark.

A record with three symbolic ALTs is split in the input buffer instead (`InputBuffer.pm` lines 638-718), its first feature `EL/UP`. The rejoin copies a merged feature's `transcript_variations` and not the `transcript_structural_variations` a structural feature holds (`OutputFactory.pm` lines 2311-2318), so the merged features' rows are lost. For `21 10542400 N <DEL>,<DUP>,<INV>` with the same `INFO`, Ensembl VEP writes 18 rows, `ALLELE_NUM` 1 and 2, at the record's span `21:10542401-10542450` with `OverlapBP=50`, `MINIMISED=1` and the positions of the trimmed span, `426-?`. vep-rs writes 27 rows, `ALLELE_NUM` 1 to 3, with the record's positions, `424-?`, and no mark.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--minimal --allele_number --uploaded_allele --variant_class`, and `--vcf` for the VCF shape; vep-rs with the same flags over the same cache.
