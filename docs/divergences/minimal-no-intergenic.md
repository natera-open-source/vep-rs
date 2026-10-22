# A `--minimal` record with one intergenic allele under `--no_intergenic`

A difference in which rows are written under `--minimal --no_intergenic` for one multi-allelic shape; nothing is set aside, since no measured dataset runs the flags.

Under `--minimal`, Ensembl VEP splits a multi-allelic record whose alleles trim into one feature per allele (`InputBuffer.pm` lines 638-718) and annotates each feature on its own. An allele that overlaps no transcript gets an intergenic overlap of its own (`VariationFeature.pm` lines 822-840). At output the split features are merged back into the record's first feature (`OutputFactory.pm` lines 2272-2354).

The merge copies a later feature's intergenic overlap onto the first feature when the first has none (lines 2335-2338, a branch commented `this probably won't happen`), and leaves the first feature's own in place when it has one. `--no_intergenic` then drops every overlap of a feature that carries an intergenic overlap (line 482), so the rejoined record writes nothing, the genic allele's transcript row included.

vep-rs judges each allele alone: `--no_intergenic` removes the intergenic row of the allele that overlaps no transcript and keeps the other allele's transcript row, its CSQ entry and its JSON `transcript_consequences`. Without `--minimal` the two engines agree, since a record's alleles then share one span and one transcript set. Under `--minimal` alone they agree too: the transcript row Ensembl VEP drops under the flag is one it writes without it.

The record below has two alternate alleles that trim to different spans. `GA` trims to `T/A` at 8885976, the last base of the downstream window of ENST00000624291, and `GTC` to an insertion of `C` after that base, one base outside the window:

```
21	8885975	str_ENST00000624291_8885976	GT	GA,GTC	.	.	.
```

Under `--minimal --allele_number` both engines write two rows: the `A` allele's `downstream_gene_variant` row on ENST00000624291 with `DISTANCE=5000` and `ALLELE_NUM=1`, and the `C` allele's `intergenic_variant` row with `ALLELE_NUM=2`. Under `--minimal --no_intergenic` they write:

| Format | Ensembl VEP | vep-rs |
| --- | --- | --- |
| default, tab | no row | the `A` row: ENST00000624291, `downstream_gene_variant`, `DISTANCE=5000` |
| VCF | the input line with INFO `.` and no `CSQ` entry | the input line with `CSQ=A\|downstream_gene_variant\|MODIFIER\|\|ENSG00000278931\|Transcript\|ENST00000624291\|...` |
| JSON | the record with `"most_severe_consequence":"?"`, `"minimised":1` and no consequence list | the record with `"most_severe_consequence":"downstream_gene_variant"` and the `A` allele's `transcript_consequences` entry |

The same holds with the alleles in the other order, where the first feature is the intergenic one and keeps its own overlap through the merge. For `21 38728001 rstr_ENST00000785411_38728002 GT AT,GTC` Ensembl VEP writes nothing under the flags and vep-rs the `C` allele's `upstream_gene_variant` row.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--minimal --no_intergenic`, and `--tab`, `--vcf` or `--json` for the other formats. vep-rs takes the same flags over the same cache, or over the cache slice `tests/golden/116/GRCh38-pick-summary/json_cache`, which holds both transcripts.
