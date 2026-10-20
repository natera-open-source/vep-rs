# SpliceAI: GRCh38 SNV only

SpliceAI adds deep-learning splice-site predictions to a transcript consequence.
It is an annotation-store plugin: it reads bgzipped, tabix-indexed VCFs and
matches by position, allele and gene symbol. This page covers its parameters and
the one data file that is supported from open data.

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
