# A `--minimal` split record under `--coding_only`

A difference in which rows are written, and in the `REF_ALLELE` they print, under `--minimal --coding_only` for a multi-allelic record whose alleles fall on both sides of a CDS boundary; nothing is set aside, since no measured dataset runs the flags.

Under `--minimal`, Ensembl VEP splits a multi-allelic record whose alleles trim into one feature per allele (`InputBuffer.pm` lines 638-718), builds each feature's rows over its own span, and merges the features back into the record at output (`OutputFactory.pm` lines 2272-2354). The `--coding_only` test, `affects_cds` of each transcript overlap (`OutputFactory.pm` lines 490-502, `BaseTranscriptVariation.pm` lines 662-665), then runs again over the rejoined record.

The default and tab factories clear the CDS coordinates cached on each allele's own span after the rejoin (`OutputFactory.pm` lines 297-323), so the test maps the record's restored span. A row then survives when the record reaches the CDS, whether or not its own allele does. The VCF and JSON factories clear them before the rejoin (`OutputFactory/VCF.pm` lines 305-308, `OutputFactory/JSON.pm` lines 217-220), so each allele is judged on its own span. vep-rs judges every row on its allele's own span in every format.

The records below sit at the acceptor site of exon 6 of ENST00000427445: 10542393 is the second intronic base before the exon and 10542394-10542397 reach its coding sequence. Each splits into an SNV at 10542393 (`splice_acceptor_variant`) and a second allele, a deletion of `TACG` over 10542394-10542397 (`splice_acceptor_variant,coding_sequence_variant`) or a second SNV at 10542393:

```
21	10542393	straddle_a	GTACG	CTACG,G	.	.	.
21	10542393	straddle_b	GTACG	CTACG,TTACG	.	.	.
21	10542393	straddle_c	GTACG	G,CTACG	.	.	.
```

In the default and tab formats, Ensembl VEP writes 14 rows for each record, both alleles on the seven coding transcripts of the site, every row with the record's Location `21:10542393-10542397`; the two transcripts the record meets in an intron alone write none. vep-rs writes the deletion's seven rows for `straddle_a` and `straddle_c` and no row for `straddle_b`. On ENST00000427445 the two engines write for `straddle_a`:

| Engine | Allele | Consequence | REF_ALLELE |
| --- | --- | --- | --- |
| Ensembl VEP | `C` | `splice_acceptor_variant` | `GTACG` |
| Ensembl VEP | `-` | `splice_acceptor_variant,coding_sequence_variant` | `GTACG` |
| vep-rs | `-` | `splice_acceptor_variant,coding_sequence_variant` | `G` |

and for `straddle_b`:

| Engine | Allele | Consequence | REF_ALLELE |
| --- | --- | --- | --- |
| Ensembl VEP | `C` | `splice_acceptor_variant` | `GTACG` |
| Ensembl VEP | `T` | `splice_acceptor_variant` | `GTACG` |

`REF_ALLELE` differs for a second reason. Ensembl VEP reads a record's reference once, on the first row it builds, and keeps the value (`OutputFactory.pm` line 1357, `VariationFeature.pm` lines 1533-1539). The rejoin builds the first allele's feature's rows before it restores the record's allele string (`OutputFactory.pm` lines 2289-2291), so every row of a split record prints that feature's trimmed reference. When the feature kept no row before the rejoin, the value is read after the restore and is the record's own.

The SNV's feature of `straddle_a` and `straddle_b` keeps no row on its own span, so Ensembl VEP prints `GTACG` on every row of the record; vep-rs prints the first allele's trimmed reference, `G`. With the deletion first, `straddle_c`, that feature keeps its rows and both engines print `TACG`.

In the VCF and JSON formats the two engines write the same rows, the deletion's seven for `straddle_a` and `straddle_c` and none for `straddle_b`, and `REF_ALLELE` (the JSON `ref_allele`) differs as above: `GTACG` against `G` on `straddle_a`, `TACG` on both sides for `straddle_c`. Without `--coding_only` the engines agree in every format, 18 rows per record on the nine transcripts.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--minimal --coding_only --allele_number --show_ref_allele --uploaded_allele --variant_class`, and `--tab`, `--vcf` or `--json` for the other formats; vep-rs with the same flags over the same cache.
