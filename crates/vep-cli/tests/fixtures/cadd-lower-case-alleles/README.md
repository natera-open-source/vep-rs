# CADD lower-case alleles fixture

A CADD-shaped SNV file (a `##` line naming GRCh38, the `#Chrom` header line and
the columns `Chrom`, `Pos`, `Ref`, `Alt`, `RawScore`, `PHRED`), bgzip-compressed
and indexed with `tabix -s 1 -b 2 -e 2`. Its two rows are fabricated, at
positions of the GRCh38 CADD corpus: a row at `21:19268555` whose `Ref` and
`Alt` cells are written in lower case, scored `7.654`, and a row at
`21:20204701` written in capitals, scored `6.543`. `variants.vcf` holds the two
records the rows describe, `lower_case_row` and `upper_case_row`.

`CADD.pm` lines 239 to 251 hand the row's cells to
`get_matched_variant_alleles`, whose keys are compared as text
(`Utils/Sequence.pm` lines 1207 and 1232), while the engine has upper-cased
the record's alleles (`Parser.pm` line 554).

Ensembl VEP 116.3 wrote, on this fixture under `snv=<file>`:
`CADD_PHRED=6.543` and `CADD_RAW=0.401234` on the row of `upper_case_row`, and
no CADD item on the row of `lower_case_row`. vep-rs writes the lower-case row's
scores on `lower_case_row` as well, the documented class
`tabix_allele_case_insensitive_match`.
