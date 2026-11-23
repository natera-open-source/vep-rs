# OpenTargets lower-case alleles fixture

A file in the layout of the Open Targets export the plugin reads (a `#` header
line naming 21 columns, the first four the chromosome, the position and the
reference and alternate alleles), bgzip-compressed and indexed with
`tabix -s 1 -b 2 -e 2`. Its two records are fabricated, lead GWAS associations
at positions of the GRCh38 OpenTargets corpus: one at `21:13980370` whose
allele cells are written in lower case, naming the gene `ENSG00000000008`, and
one at `21:14281918` written in capitals, naming `ENSG00000000009`.
`variants.vcf` holds the two records they describe, `lower_case_row` and
`upper_case_row`.

`OpenTargets.pm` lines 376 to 388 hand the record's cells to
`get_matched_variant_alleles`, whose keys are compared as text
(`Utils/Sequence.pm` lines 1207 and 1232), while the engine has upper-cased
the input's alleles (`Parser.pm` line 554).

Ensembl VEP 116.3 wrote, on this fixture under `file=<file>`: the five default
fields, `OpenTargets_gwasGeneId=ENSG00000000009` among them, on every row of
`upper_case_row`, and no OpenTargets item on any row of `lower_case_row`.
vep-rs writes the lower-case record's fields on the rows of `lower_case_row`
as well, the documented class `tabix_allele_case_insensitive_match`.
