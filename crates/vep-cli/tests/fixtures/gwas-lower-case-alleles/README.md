# GWAS lower-case risk allele fixture

A curated-catalog file of the header line and two rows of the GRCh38 GWAS
corpus's catalog slice
(`tests/golden/116/GRCh38-plugin-gwas/plugin_data/gwas_catalog_slice.tsv`, its
second and twentieth lines): the row of `rs2834259` at `21:33760786` with its
risk allele rewritten in lower case (`rs2834259-t`), and the row of
`rs2801270` at `21:13297803` as the slice writes it (`rs2801270-G`).
`variants.vcf` holds the two records the rows describe, `lower_case_row` and
`upper_case_row`. `vepdir/` holds the processed catalogue Ensembl VEP's plugin
writes under `--dir` on its first run (`GWAS.pm` lines 106 to 118), so a run
whose `--dir` is that directory, mounted read-only, finds it there.

`GWAS.pm` lines 182 to 194 hand the row's risk allele to
`get_matched_variant_alleles`, whose keys are compared as text
(`Utils/Sequence.pm` lines 1207 and 1232), while the engine has upper-cased
the input's alleles (`Parser.pm` line 554).

Ensembl VEP 116.3 wrote, on this fixture under `file=<file>` with a writable
`--dir`: the GWAS fields, `GWAS_risk_allele=G` and `GWAS_study=GCST008758`
among them, on the row of `upper_case_row`, and no GWAS item on any row of
`lower_case_row`. vep-rs writes the lower-case row's fields on the rows of
`lower_case_row` as well, `GWAS_risk_allele=t` as the file writes it, the
documented class `tabix_allele_case_insensitive_match`.
