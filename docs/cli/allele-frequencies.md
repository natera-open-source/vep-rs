# Allele frequencies and co-located variants

These flags add population allele frequencies, PubMed identifiers and gene-phenotype associations to the output. All of them are switched on by `--everything`. `--check_existing`, which reports the co-located known variants themselves, is described under [Output control](output-control.md#opt_check_existing).

## Summary

| Flag | Default | What it does |
| --- | --- | --- |
| [`--af`](#opt_af) | off | Global allele frequency from 1000 Genomes Phase 3 |
| [`--af_1kg`](#opt_af_1kg) | off | Continental frequencies: AFR, AMR, EAS, EUR, SAS |
| [`--af_gnomade`](#opt_af_gnomade) | off | gnomAD exome frequencies |
| [`--af_gnomadg`](#opt_af_gnomadg) | off | gnomAD genome frequencies |
| [`--max_af`](#opt_max_af) | off | Highest frequency across all populations |
| [`--pubmed`](#opt_pubmed) | off | PubMed IDs for co-located variants |
| [`--gene_phenotype`](#opt_gene_phenotype) | off | Gene-phenotype associations |

## Options

### `--af` {#opt_af}

**`--af`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the global allele frequency from 1000 Genomes Phase 3.

### `--af_1kg` {#opt_af_1kg}

**`--af_1kg`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the continental frequencies: AFR, AMR, EAS, EUR and SAS.

### `--af_gnomade` {#opt_af_gnomade}

**`--af_gnomade`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the gnomAD exome frequencies.

### `--af_gnomadg` {#opt_af_gnomadg}

**`--af_gnomadg`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the gnomAD genome frequencies.

### `--max_af` {#opt_max_af}

**`--max_af`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the highest frequency across all populations.

### `--pubmed` {#opt_pubmed}

**`--pubmed`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the PubMed IDs for co-located variants.

### `--gene_phenotype` {#opt_gene_phenotype}

**`--gene_phenotype`**, boolean, default off. Ensembl VEP: same flag and default.

Adds the gene-phenotype associations.
