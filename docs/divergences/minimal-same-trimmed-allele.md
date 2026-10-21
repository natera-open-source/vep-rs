# Two `--minimal` alleles that trim to the same allele string

A difference in which rows are written under `--minimal` for one multi-allelic shape; nothing is set aside, since no measured dataset runs the flag.

Under `--minimal`, Ensembl VEP splits a multi-allelic record whose alleles trim into one feature per allele (`InputBuffer.pm` lines 638-718) and, at output, merges every allele after the first back into the record's feature under the key `<trimmed allele string>_<transcript key>` (`OutputFactory.pm` line 2316). Two alternate alleles that trim to the same allele string at different positions share that key, and the later allele's rows replace the earlier allele's. The record below has three alternate alleles, which trim to `A/C` at 10542414, `A/T` at 10542414 and `A/T` at 10542415:

```
21	10542414	collide	AA	CA,TA,AT	.	.	.
```

Ensembl VEP writes the rows of the first and third alleles alone, the second allele's `stop_gained` among the rows lost. vep-rs sorts the rows by the same key and merges nothing, so every allele keeps its rows: on the nine transcripts of the cache slice `tests/golden/116/GRCh38-identifiers/json_cache`, 18 rows against 27. On ENST00000622113 the two engines write:

| Engine | Allele | Consequence | ALLELE_NUM |
| --- | --- | --- | --- |
| Ensembl VEP | `T` | `missense_variant` | 3 |
| Ensembl VEP | `C` | `missense_variant` | 1 |
| vep-rs | `T` | `stop_gained` | 2 |
| vep-rs | `T` | `missense_variant` | 3 |
| vep-rs | `C` | `missense_variant` | 1 |

Every row of both engines prints the record's Location, `21:10542414-10542415`, and `UPLOADED_ALLELE=AA/CA/TA/AT`. A repeated alternate allele after the first in a record that splits is one key too: `GAA GGT,GGA,GGA` trims to `AA/GT`, `A/G` and `A/G`, and Ensembl VEP writes the rows of alleles 1 and 3 where vep-rs writes all three. Under `--pick` each engine ranks the rows it keeps, so on `AAA AGA,AGA,AGA` Ensembl VEP picks the row with `ALLELE_NUM=3` and vep-rs the row with `ALLELE_NUM=2`. A record none of whose alleles trims is not split, and both engines read a repeated allele as one allele under its first number.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--minimal --allele_number --uploaded_allele`; vep-rs with the same flags over the same cache.
