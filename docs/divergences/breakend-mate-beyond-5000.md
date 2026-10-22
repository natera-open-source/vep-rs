# A breakend mate between 5,000 bases and `--distance` of a transcript

A difference in which records are written, under a `--distance` above 5,000 on a breakend whose mate lies between 5,000 bases and that distance from a transcript; nothing is set aside, since no measured run names the flag.

`--distance` rewrites the two package variables that size a transcript's flanks (`BaseRunner.pm` lines 511-512). Ensembl VEP collects the transcripts of a buffer through those flanks (`AnnotationType/Transcript.pm` lines 96-108), and a bracket breakend sits in the buffer's tree at its mate's coordinate as well as its own (`InputBuffer.pm` lines 387-389). So a transcript within `--distance` of the mate collects the record and is given a `TranscriptStructuralVariation` for it (`AnnotationType/Transcript.pm` lines 145-151).

The alleles of that object are admitted one by one through `_close_to_feature` (`StructuralVariationOverlap.pm` lines 90-102). That test widens the transcript by `MAX_DISTANCE_FROM_TRANSCRIPT` (line 157), a constant of 5,000 bases (`Utils/VariationEffect.pm` line 60) that the flag never changes. A mate more than 5,000 bases from the transcript fails it, so the object holds no allele. The record is not intergenic either, since the intergenic row is built only when no transcript object exists (`StructuralVariationFeature.pm` line 679), and the output maps the alleles of its transcript objects (`OutputFactory.pm` line 560). The record writes no row at all.

vep-rs admits the mate within the larger of the two distances the flag gives, so the bracket allele's row is written with the mate's distance to the transcript. At the default distance the two windows coincide, and every measured set runs at the default. Up to 5,000 bases the two engines write the row and differ in its terms alone ([class 11](defect-e.md#11-breakend-mate-rows-read-from-the-other-chromosomes-coordinate)).

The records below sit 13 kb above ENST00000286808 (reverse strand, 30165565 to 30166805 on GRCh38), their own position intergenic. The first record's mate lies 5,500 bases above the transcript's end, the second's 4,500:

```
21	30180000	far_mate_5500	N	N[21:30172305[	.	.	SVTYPE=BND
21	30180000	near_mate_4500	N	N[21:30171305[	.	.	SVTYPE=BND
```

Under `--distance 6000 --allele_number` over the GRCh38 cache, the two engines write for the bracket allele (`ALLELE_NUM=1`) on ENST00000286808 the rows below: the Ensembl VEP column follows from the lines above, and the reproduction at the end runs both engines:

| Record | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `far_mate_5500` | no row, and no intergenic row either | `upstream_gene_variant`, `DISTANCE=5500` |
| `near_mate_4500` | `intergenic_variant` on the transcript row (class 11) | `upstream_gene_variant`, `DISTANCE=4500` |

Under the default `--distance 5000` both engines write `far_mate_5500` as two `intergenic_variant` rows, one per allele, and `near_mate_4500` as above. A mate exactly 5,000 bases from the transcript's end passes the test at line 157 and is written by both engines under either distance; one base farther fails it, so under `--distance 6000` Ensembl VEP writes nothing where vep-rs writes `DISTANCE=5001`.

The flag's purpose is to widen the window, and a record that disappears from the output when it widens is a defect of the admission test, not a reading vep-rs reproduces.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--distance 6000 --allele_number`; vep-rs with the same flags over the same cache.
