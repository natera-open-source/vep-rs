# The motif score of a multi-base allele on a reverse-strand motif

**In short.** `MOTIF_SCORE_CHANGE` is the change in a motif's binding score when the allele is written over the motif's sequence. For a multi-base allele on a reverse-strand motif, Ensembl VEP writes the allele starting at the motif position of the variant's genomic start. Along a reverse-strand motif that is the allele's last position, so the allele lands one or more positions too far along the motif. At the motif's end it runs past the sequence and no score is written. vep-rs writes the allele over the positions it covers and scores that. A one-base allele occupies one position and the two engines agree on it; nothing is set aside in the concordance figures, since no measured dataset runs `--regulatory`.

## What a reader needs to know first

A `MotifFeature` row of a sequence variant carries `MOTIF_POS`, the variant's start as a position in the motif, `HIGH_INF_POS`, whether that position carries more than 1.5 bits of information, and `MOTIF_SCORE_CHANGE` ([Motif fields](../output/vep-and-tab.md#motif-fields)). The positions of a motif run along its own strand: on a reverse-strand motif position 1 is the motif's highest genomic coordinate. Ensembl VEP keeps a motif's sequence on its own strand in the cache, and compares an allele with it on that strand, reverse-complemented where the strands differ (`VariationFeatureOverlapAllele.pm` lines 246 to 264).

A variant covering genomic coordinates `s` to `e` inside a reverse-strand motif of length `L` starting at `m` covers motif positions `L - (e - m)` to `L - (s - m)`. `motif_start`, the position of the genomic start `s`, is the larger of the two (`MotifFeatureVariationAllele.pm` lines 99 to 103), and it is what `MOTIF_POS` prints on both engines.

## The records

Four records of chromosome 21 of GRCh38 over motifs of the release 116 cache. Three are multi-base alleles over the forward motif ENSM00000574814 (ENSPFM0290, 5090718 to 5090731) and the reverse motif ENSM00000574815 (ENSPFM0484, 5090723 to 5090738). The fourth lies over the NRF1 motifs ENSM00000588223 (reverse) and ENSM00000588224 (forward), both 46699248 to 46699259.

```
21	5090720	fwd_mnv_inside	GG	AA	.	.	.
21	5090723	rev_mnv_start	CA	TT	.	.	.
21	5090730	rev_mnv_inside	CT	AG	.	.	.
21	46699250	nrf1_rev_mnv	CG	TA	.	.	.
```

## What each engine writes

The motif rows of the default format, `Extra` column alone, every other column equal between the engines. The forward-strand rows agree.

| Record | Feature | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `fwd_mnv_inside` | ENSM00000574814 (+) | `MOTIF_POS=3;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.108` | the same |
| `rev_mnv_start` | ENSM00000574814 (+) | `MOTIF_POS=6;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.111` | the same |
| `rev_mnv_start` | ENSM00000574815 (-) | `MOTIF_POS=16;HIGH_INF_POS=N` | `MOTIF_POS=16;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=0.044` |
| `rev_mnv_inside` | ENSM00000574814 (+) | `MOTIF_POS=13;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.016` | the same |
| `rev_mnv_inside` | ENSM00000574815 (-) | `MOTIF_POS=9;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.092` | `MOTIF_POS=9;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.136` |
| `nrf1_rev_mnv` | ENSM00000588223 (-) | `MOTIF_POS=10;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.126` | `MOTIF_POS=10;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.143` |
| `nrf1_rev_mnv` | ENSM00000588224 (+) | `MOTIF_POS=3;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE=-0.149` | the same |

Each row's `Extra` column opens with `IMPACT=MODIFIER;STRAND=<strand>;MOTIF_NAME=<matrix>` and closes with `TRANSCRIPTION_FACTORS=<factors>`, equal on both sides.

## Why Ensembl VEP writes that

`motif_score_delta` (`MotifFeatureVariationAllele.pm` lines 183 to 265) takes the motif's sequence, scores it, writes the allele over it and scores the result. The write is `substr($mf_seq, $mf_start - 1, $var_len) = $allele_seq` (line 245): the allele, on the motif's strand, starting at `motif_start`. On a forward-strand motif `motif_start` is the allele's first position and the write is right. On a reverse-strand motif `motif_start` is the position of the genomic start, the allele's last position along the motif. So a two-base allele that covers positions 8 and 9 is written over 9 and 10: `rev_mnv_inside` on ENSM00000574815 puts `CT` over `GA` instead of over `AG`, scored `-0.092`.

When the allele covers the motif's last positions, the write starts at the last position and runs past the sequence's end, which Perl's `substr` lengthens the string for. The length check that follows (lines 248 to 251) then finds the sequence longer than the motif and returns no score, so `MOTIF_SCORE_CHANGE` is left out (`rev_mnv_start` on ENSM00000574815, positions 15 and 16).

## Why vep-rs's rows are the correct ones

The allele covers the positions `motif_end` to `motif_start` of a reverse-strand motif, both computed by Ensembl VEP itself (lines 99 to 103 and 130 to 134), and its reverse complement reads along the motif from `motif_end`. Writing it there changes the bases the variant changes and no others: `rev_mnv_inside` replaces the `AG` at positions 8 and 9 of `TGACCCCAGACCTCTG` (the genomic `CT` at 5090730 and 5090731) with `CT`, and the score change of that sequence is `-0.136`. Writing it from `motif_start` leaves position 8 as it was, puts the wrong base at position 9 and changes position 10, a base the variant does not touch.

`MOTIF_POS` keeps Ensembl VEP's value on both engines, the genomic start's position, a defensible definition of the field.

## Try it yourself

With a converted GRCh38 release 116 cache holding its regulatory tree ([Regulatory and motif features](../regulatory.md)), write the four records above to `motifs.vcf` and run:

```
vep --offline --json_cache <cache> --regulatory -i motifs.vcf -o motifs.txt
```

The seven motif rows carry the vep-rs column of the table. Ensembl VEP 116.3 on the same records with `--offline --cache --regulatory`, on the release 116 cache from the Ensembl FTP site, writes its column.

## How the comparator treats this difference

The comparators score consequence sets and never read `MOTIF_SCORE_CHANGE`, so nothing is set aside in any figure. The field comparator (`tools/vep-fields-compare`) compares the field exactly; a run that compares it on records with a multi-base allele over a reverse-strand motif names each such row as a field divergence of this page. The format comparator (`tools/vep-formats-compare`) takes the class `motif_score_reverse_strand_mnv` under `--masks`. A corpus manifest documents each such row's value. On an input with no manifest the comparator accepts a differing or added `MOTIF_SCORE_CHANGE` on a reverse-strand motif row whose allele spans two or more bases, and nothing else of the row.
