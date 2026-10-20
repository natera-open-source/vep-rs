# A repeated alternate allele the looked-up reference equals

A difference in which rows are written under `--lookup_ref --minimal` for a VCF record that repeats an alternate allele; nothing is set aside, since no measured dataset runs the flags.

Ensembl VEP reads a record's alternate alleles as written. It replaces the record's reference with the FASTA's bases under `--lookup_ref` (`Parser.pm` lines 635-646) and splits the record afterwards (`InputBuffer.pm` lines 638-718). A record whose repeated alternate allele the looked-up reference equals therefore splits into one feature per repeat, each `-/-`. Every feature builds its reference allele (`VariationFeatureOverlap.pm` lines 469-474), so under `--process_ref_homs` each repeat writes a `-` reference row numbered 0. Alleles after the first that trim to one allele string share an output key in Ensembl VEP, the later replacing the earlier ([two alleles that trim alike](minimal-same-trimmed-allele.md)); every record below repeats its first alternate allele once. At a site whose reference base is `A`, the record

```
21	10542415	rep	T	A,A	.	.	.
```

becomes `A/A/A`, and Ensembl VEP writes two `-` rows per transcript. In a VCF record whose written alleles do not trim, vep-rs reads a repeated alternate allele as one allele before the reference is looked up. The record therefore has one feature `-/-` and writes one `-` row per transcript. For `T A,G,A` at the same site Ensembl VEP writes `-`, `A`, `G`, `-` per transcript and vep-rs `A`, `G`, `-`. Every other field of the rows both engines write agrees, `REF_ALLELE` `-` and `UPLOADED_ALLELE` `A/A/A` included.

The same reading drops a repeat that trims against the looked-up reference without equalling it. `T GA,GA` at the site is `A/GA/GA`, two features `-/G` in Ensembl VEP, which writes `G` numbered 1 and 2 per transcript and, under `--process_ref_homs`, a `-` reference row for each; vep-rs writes `G` numbered 1 and one `-` row.

A record none of whose alleles equals the looked-up reference or trims against it (`T G,G` at the site, `A/G/G`) is not split, and its repeated allele is one allele under its first number in both engines. A VCF record whose written alleles trim reads alike in both engines: the repeat is annotated when the looked-up alleles trim too (`TC AC,AC`, `AA/AC/AC`), and the repeat once when they do not (`TC TG,TG`, `AA/TG/TG`). An Ensembl-format record reads alike in both engines, since vep-rs keeps its repeats through the parse and reads them after the reference is looked up, as Ensembl VEP does.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache, a FASTA reading `A` at 21:10542415 and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--fasta <file> --lookup_ref --minimal --process_ref_homs --allele_number --show_ref_allele --uploaded_allele`; vep-rs with the same flags over the same cache and FASTA.
