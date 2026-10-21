# SpliceAI: GRCh38 SNV only

SpliceAI adds deep-learning splice-site predictions to a transcript consequence.
It is an annotation-store plugin: it reads bgzipped, tabix-indexed VCFs and
matches by position, allele and gene symbol. This page covers its parameters,
the one data file that is supported from open data, and how a reverse-strand
indel is matched.

The supported SpliceAI data is the openly published GRCh38 SNV score file from
the Ensembl MANE release. GRCh37 scores and the indel files for both assemblies
are distributed solely through Illumina BaseSpace, which requires an account
signup, so they are not supported. They drop in unchanged for a user who holds
credentials: same parameters, same tabix contract.

## Parameters and data {#spliceai}

| Property      | Value                                                                                                                                                                                                          |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh38, SNV scores                                                                                                                                                                                             |
| Matches       | position, allele and gene symbol                                                                                                                                                                               |
| Output fields | `SpliceAI_pred_DS_AG`, `SpliceAI_pred_DS_AL`, `SpliceAI_pred_DS_DG`, `SpliceAI_pred_DS_DL`, `SpliceAI_pred_DP_AG`, `SpliceAI_pred_DP_AL`, `SpliceAI_pred_DP_DG`, `SpliceAI_pred_DP_DL`, `SpliceAI_pred_SYMBOL` |
| Data source   | [Ensembl MANE release](https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/) (GRCh38 SNV)                                                                                             |

| Parameter        | Meaning                                                                      |
| ---------------- | ---------------------------------------------------------------------------- |
| `snv=<path>`     | The SNV score VCF.                                                           |
| `indel=<path>`   | The indel score VCF; a BaseSpace-gated file, for a user who has one.         |
| `cutoff=<float>` | Report a record only when its largest delta score reaches this value.        |

```bash
--plugin SpliceAI,\
snv=/path/to/spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz
--plugin SpliceAI,snv=/path/to/spliceai_snv.vcf.gz,\
indel=/path/to/spliceai_indel.vcf.gz,cutoff=0.5
```

**Preparing the data.** The GRCh38 SNV score VCF, bgzipped with its `.tbi`,
comes from the Ensembl MANE release at
[ftp.ensembl.org](https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/).
Remote tabix works against it, so a chr21 subset can be sliced without the full
download.

## Reverse-strand indels

An Ensembl-format or region record given on strand -1 is matched to the score
file on its forward-strand alleles, as every tabix-backed plugin matches it. A
reverse-strand indel therefore scores the row its forward-strand twin scores.
Ensembl VEP 116.2's SpliceAI plugin complements such an indel twice: it reads
the forward-strand alleles `to_VCF_record` writes (`VariationFeature.pm` line
2376), then reverse-complements them again (`SpliceAI.pm` lines 259-271 of the
`VEP_plugins` release/116 branch). The alleles it compares are the
reverse-strand ones, and no row of the file matches.

For `21 26021981 26021980 -/CTT -1` against a file holding the forward-strand
insertion `G GAAG` at 26021980, Ensembl VEP 116.2 writes no SpliceAI field on
any of the record's 57 rows, and vep-rs scores that row on every transcript row.
A reverse-strand SNV is complemented once in both engines and matches alike. The
[index of differences](../intended-divergences.md#differences-outside-the-classes)
lists it.
