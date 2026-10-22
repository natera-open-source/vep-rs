# HGVSg on a record whose ALT list repeats an allele or names the reference

A difference in the `HGVSg` text of a multi-allelic record whose ALT column repeats an allele or names the reference base, not in any tuple; nothing is set aside, since the comparators score consequence sets and never see the field.

Under `--hgvsg`, Ensembl VEP computes a notation for each alternate allele of the record's allele string, every written ALT kept (`Parser/VCF.pm` line 342), and stores it under the written allele (`VariationFeature.pm` line 2039). A row prints the notation stored under one entry of that string: the second entry when the string has two, else the entry at the row's allele number (`OutputFactory.pm` lines 1386-1393). The allele number is counted over a list from which a repeated allele (`VariationFeatureOverlap.pm` line 463) and an allele equal to the reference (line 481) are dropped first. So the number and the string fall out of step at every allele written after such an entry, and the row prints another allele's notation.

The allele that equals the reference gets a notation although no row is written for it. The multi-allelic trim (`VariationFeature.pm` line 1970) empties `A/A` to a span of no bases after 34775850; the pair is neither a deletion nor an insertion (lines 1981-1982), so the untrimmed `A` stays the sequence checked (line 1946). `hgvs_variant_notation` finds the base before the site equal to it and writes a duplication (`Utils/Sequence.pm` lines 570-589). A duplication of the reference base describes no change. The lookup indexes the string by the row's allele number, so its intent is the row's own allele; the value it reaches here is an accident of two lists counted differently, not a reading of the record.

vep-rs computes each row's notation from the row's own allele at its trimmed span, so a repeated or reference-equal entry changes no other row's value. Both engines count `ALLELE_NUM` alike, and neither writes a row for the allele that equals the reference.

The records sit on GRCh37 at the first base of a start codon of ENSG00000159128's transcripts, where the reference reads `ATG`:

```
21	34775850	refeq	A	A,C	.	.	.
21	34775850	dupalt	A	C,C,G	.	.	.
```

Under `--hgvsg --allele_number`, the two engines write on each of the five transcripts of that gene:

| Record and Allele | ALLELE_NUM | Ensembl VEP HGVSg | vep-rs HGVSg |
| --- | --- | --- | --- |
| `refeq`, `C` | 1 | `21:g.34775850dup` | `21:g.34775850A>C` |
| `dupalt`, `C` | 1 | `21:g.34775850A>C` | `21:g.34775850A>C` |
| `dupalt`, `G` | 2 | `21:g.34775850A>C` | `21:g.34775850A>G` |

Every other column of the fifteen rows agrees. A record whose ALT list carries no repeat and no reference-equal entry prints the same `HGVSg` in both engines, whatever its allele count.

To reproduce: Ensembl VEP 116.2 with the GRCh37 cache and FASTA and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--hgvsg --allele_number`; vep-rs with the same flags over the same cache and FASTA.
