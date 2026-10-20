# Allele frequencies and co-located variants

These flags report the population allele frequencies, PubMed identifiers and gene-phenotype associations of co-located known variants, and filter the input by frequency. `--check_existing`, which reports the co-located variants themselves, is listed under [Output control](output-control.md); the frequency flags `--everything` switches on are listed under [option sets](../cli-reference.md#option-sets-configpm351-593).

## Flags

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--af` | bool | false | Add global allele frequency from 1000 Genomes Phase 3 |
| `--af_1kg` | bool | false | Add continental allele frequencies from 1000 Genomes |
| `--af_gnomade` | bool | false | Add gnomAD exome allele frequencies |
| `--af_gnomad` | bool | false | Add gnomAD exome allele frequencies under the older spelling of `--af_gnomade`; the two spellings cannot be combined, and this one alone fails on a cache without gnomAD exome data |
| `--af_gnomadg` | bool | false | Add gnomAD genome allele frequencies |
| `--max_af` | bool | false | Add highest allele frequency across populations |
| `--pubmed` | bool | false | Add PubMed IDs for co-located variants |
| `--old_maf` | bool | false | Rewrite the AFR, AMR, ASN and EUR frequency values of each co-located variant as a bare minor-allele frequency: the allele is stripped, the first value kept and a value above 0.5 inverted; the allele-keyed frequency columns then find no allele in them and print empty |
| `--check_frequency` | bool | false | Filter the input by the frequency of co-located known variants in the population `--freq_pop` names, against the threshold `--freq_freq` on the side `--freq_gt_lt` gives, dropping the alleles that pass or keeping the records whose alleles all pass as `--freq_filter` says; switches on `--check_existing` and adds the FREQS field (the frequencies the filter read, NA where an allele has none) |
| `--freq_pop` | String | `1KG_ALL` | The population `--check_frequency` reads: 1KG\_ALL (the 1000 Genomes AF column of the cache), a frequency column of the cache in any case (AF, AFR, gnomADe\_NFE; the gnomAD genome columns are refused, as Ensembl VEP refuses them), or 1KG\_ followed by such a column |
| `--freq_freq` | f64 | 0.01 | The frequency `--check_frequency` compares against |
| `--freq_gt_lt` | `gt`, `lt` | `gt` | Whether an allele passes `--check_frequency` with a frequency at or above (gt) or at or below (lt) `--freq_freq` |
| `--freq_filter` | `exclude`, `include` | `exclude` | Whether `--check_frequency` drops the alleles that pass (exclude) or keeps a record only when every allele passes (include) |
| `--filter_common` | bool | false | Drop the alleles with a 1000 Genomes frequency of 0.01 or more: `--check_frequency` with `--freq_pop` 1KG\_ALL `--freq_freq` 0.01 `--freq_gt_lt` gt `--freq_filter` exclude |
| `--gene_phenotype` | bool | false | Add the GENE\_PHENO column: 1 on a transcript whose gene the cache associates with a phenotype, disease or trait, empty otherwise |

