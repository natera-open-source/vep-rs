# An individual's copy of a `--minimal` record under `--safe`

A difference in the alleles annotated for an individual under a flag pair Ensembl VEP accepts only with `--safe`; nothing is set aside, since no measured dataset runs the flags.

Ensembl VEP lists `--individual` against `--minimal` as incompatible (`Config.pm` line 626) and skips every incompatibility check under `--safe` (line 961, `unless($config->{safe})`), a flag whose stated purpose is to die when a plugin fails to compile (line 272). vep-rs skips the same check under `--safe`, so `--safe --individual all --minimal` runs in both engines.

Ensembl VEP post-processes the record before it copies it per individual (`Parser/VCF.pm` lines 362-364). `minimise_alleles` (`Parser.pm` lines 852 and 964-1019) trims a record with one alternate allele and marks it; a record with two or more passes untouched, its written string kept as `original_allele_string` (lines 978-984). Each copy is then a feature of its own in the input buffer, where `split_variants` (`InputBuffer.pm` lines 638-718) splits a feature only when its allele string names two alternate alleles (line 648) and one of them trims.

A copy carrying one alternate allele of a multi-allelic record, or two genotype bits naming one allele, is therefore annotated on the record's anchor-chopped alleles. Nothing minimises it, it carries no `MINIMISED`, and its `UPLOADED_ALLELE` is the record's string. vep-rs reads such a copy the same way.

The records below, with the samples `S1`, `S2` and `S3`, reach nine transcripts of the GRCh38 cache:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
21	10542414	bi_trim	AAAG	AG	.	.	.	GT	0/1	1/1	0/0
21	10542419	multi_trim	AGCAA	CGCAA,A	.	.	.	GT	1/2	0/1	2/2
21	10542424	multi_dup	CCG	CGG,CGG	.	.	.	GT	1/2	0/2	1/1
```

Both engines write the `0/1` carrier of `multi_trim` as `CGCAA` over `AGCAA`, a `substitution` with `ALLELE_NUM` 1 and `UPLOADED_ALLELE` `AGCAA/CGCAA/A`, and its `2/2` carrier as `A` over `AGCAA`, an `indel`. Every carrier of `multi_dup` is `CGG` over `CCG`, a `substitution`, the `1/2` genotype read as `HOM` since both bits name one allele. Those five copies agree column for column on every transcript in the default, tab, VCF and JSON outputs, codons and peptides included. The list form `--individual S1` reads `S1` as `all` reads it.

## A copy of two alleles that trim

The `1/2` carrier of `multi_trim` names two alternate alleles that trim, so both engines split its copy and rejoin the rows under the record's span, `21:10542419-10542423`. Ensembl VEP lists the copy's alleles in Perl's hash order ([an individual carrying two alternate alleles](individual-allele-order.md)). Under `--minimal` the first allele of that order also decides the copy's `REF_ALLELE`, its trimmed reference, and the `VARIANT_CLASS` of every row (`InputBuffer.pm` lines 694-697, `OutputFactory.pm` lines 1357-1360). `UPLOADED_ALLELE` and the JSON `allele_string` print the copy's string in that order.

Under a fixed hash seed Ensembl VEP lists the deletion first; vep-rs lists the record's ALT order. On ENST00000427445 the two engines write:

| Engine | Allele | ALLELE_NUM | Consequence |
| --- | --- | --- | --- |
| Ensembl VEP | `-` | 1 | `frameshift_variant` |
| Ensembl VEP | `C` | 2 | `synonymous_variant` |
| vep-rs | `C` | 1 | `synonymous_variant` |
| vep-rs | `-` | 2 | `frameshift_variant` |

| Engine | REF_ALLELE | VARIANT_CLASS | UPLOADED_ALLELE |
| --- | --- | --- | --- |
| Ensembl VEP | `GCAA` | `deletion` | `AGCAA/A/CGCAA` |
| vep-rs | `A` | `SNV` | `AGCAA/CGCAA/A` |

Each allele's consequence terms, codons and peptides agree, and both engines print `MINIMISED` 1.

## The copies of an indel minimised past the anchor chop

`bi_trim` is a record with one alternate allele, so `minimise_alleles` trims it from the chopped `AAG/G` to `AA/-` at 21:10542415-10542416 before the copies are made. Ensembl VEP then builds each copy from the chopped alleles against the minimised reference `AA` ([an individual's copy of an indel the anchor chop leaves trimmable](individual-trimmable-indel.md)). The `0/1` carrier's copy is `AA/AAG/G`, three alleles, which the split trims to `-/G` and `AA/G`; the `1/1` carrier's is `AA/G`; the `0/0` individual's is `AA/AAG`, read as `HOM`. vep-rs compares chopped alleles with the chopped reference, so both carriers annotate the minimised `AA/-` and the `0/0` individual is reference-homozygous. On ENST00000427445:

| Engine | Genotype | Rows (Allele, ALLELE_NUM, ZYG) | REF_ALLELE, VARIANT_CLASS |
| --- | --- | --- | --- |
| Ensembl VEP | `0/1` | `G` 1 `HET`; `G` 2 `HET` | `-`, `insertion` |
| Ensembl VEP | `1/1` | `G` 1 `HOM` | `AA`, `indel` |
| Ensembl VEP | `0/0` | `AAG` 1 `HOM` | `AA`, `indel` |
| vep-rs | `0/1` | `-` 1 `HET` | `AA`, `deletion` |
| vep-rs | `1/1` | `-` 1 `HOM` | `AA`, `deletion` |
| vep-rs | `0/0` | none | none |

Every row of both engines reads `frameshift_variant` on this transcript and prints `MINIMISED` 1. `UPLOADED_ALLELE` is `AAAG/AG` on every row but the `0/1` carrier's in Ensembl VEP, which print the split copy's own string, `AA/AAG/G`. In the VCF output the `0/0` individual's line carries a `CSQ` in Ensembl VEP and none in vep-rs; in the JSON output its object has `allele_string` `AA/AAG` and `variant_class` `indel` in Ensembl VEP, and `AA`, `sequence_alteration` and `most_severe_consequence` `?` in vep-rs.

Over the nine transcripts Ensembl VEP writes 99 default rows to vep-rs's 81: the second `G` row of the `0/1` carrier and the `0/0` individual's `AAG` row, nine each. Of the remaining 81 rows, 45 agree in every column and 36, the `1/2` carrier of `multi_trim` and the two carriers of `bi_trim`, differ as the tables show.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--safe --individual all --minimal --allele_number --show_ref_allele --uploaded_allele --variant_class`, the input written with its `FORMAT` and sample columns; vep-rs with the same flags over the same cache. The order of the `1/2` copy's alleles in Ensembl VEP changes with Perl's hash seed (`PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0` fixes it).
