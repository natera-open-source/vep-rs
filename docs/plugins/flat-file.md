# Flat-file plugins

A flat-file plugin loads its whole data file at init and looks records up in
memory, by gene symbol or by position. This page documents LoFtool, pLI and
GWAS. LoFtool and pLI are keyed by gene symbol, so one data file serves both
assemblies. GWAS is GRCh38-only because the NHGRI-EBI catalog publishes
`CHR_POS` in GRCh38 coordinates only.

## LoFtool {#loftool}

LoFtool adds a gene's intolerance to functional variation to a transcript
consequence.

| Property      | Value                                                                                                                             |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                                    |
| Matches       | gene symbol lookup                                                                                                                |
| Output fields | `LoFtool_percentile`                                                                                                              |
| Data source   | `LoFtool_scores.txt` from the [Ensembl VEP_plugins](https://github.com/Ensembl/VEP_plugins) repository; vep-rs does not bundle it |

| Parameter | Meaning                                              |
| --------- | ---------------------------------------------------- |
| `<path>`  | The scores file, as the only parameter. Required.    |

```bash
--plugin LoFtool,/path/to/LoFtool_scores.txt
```

**Preparing the data.** Fetch `LoFtool_scores.txt` from the VEP_plugins
repository and pass its path.

## pLI {#pli}

pLI adds a gene's probability of loss-of-function intolerance to a transcript
consequence.

| Property      | Value                                                                                                                                                                                      |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Builds        | GRCh37, GRCh38                                                                                                                                                                             |
| Matches       | gene symbol lookup                                                                                                                                                                         |
| Output fields | `pLI_gene_value`                                                                                                                                                                           |
| Data source   | `pLI_values.txt` from the [Ensembl VEP_plugins](https://github.com/Ensembl/VEP_plugins) repository (ExAC r0.3 gene pLI); any uncompressed TSV with `gene` and `pLI` header columns is accepted |

| Parameter | Meaning                                              |
| --------- | ---------------------------------------------------- |
| `<path>`  | The gene constraint file, as the only parameter.     |

```bash
--plugin pLI,/path/to/pLI_values.txt
```

## GWAS {#gwas}

GWAS adds NHGRI-EBI GWAS Catalog associations to a variant.

| Property      | Value                                                                                                                                                                             |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh38                                                                                                                                                                            |
| Matches       | position and risk allele                                                                                                                                                          |
| Output fields | `GWAS_accessions`, `GWAS_associated_gene`, `GWAS_beta_coef`, `GWAS_odds_ratio`, `GWAS_p_value`, `GWAS_pmid`, `GWAS_risk_allele`, `GWAS_study` (the Perl GWAS plugin's field set) |
| Data source   | [GWAS Catalog downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads)                                                                                                          |

| Parameter     | Meaning                                                      |
| ------------- | ------------------------------------------------------------ |
| `file=<path>` | The catalog associations file; a bare path is accepted too.  |

```bash
--plugin GWAS,file=/path/to/gwas_catalog_associations.tsv
```

**Preparing the data.** Download the associations TSV from
[GWAS Catalog downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads). The
plugin reads it plain or gzipped and locates its columns by header name, so no
reformatting is needed.
