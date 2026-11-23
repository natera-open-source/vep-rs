# LOEUF trailing space fixture

A table of 94 columns whose header names the three columns the plugin reads
where its v4 layout expects them (`gene_id` second, `transcript` third,
`lof.oe_ci.upper` ninety-fourth and last, `LOEUF.pm` lines 89 to 92), with the
coordinate columns `chrom`, `start` and `end` fourth to sixth, bgzip-compressed
and indexed with `tabix -s 4 -b 5 -e 6`. Its two rows are fabricated, for two
transcripts of one gene of the GRCh38 LOEUF corpus: `ENST00000346798` with the
value `5.0000e-01` followed by one space, and `ENST00000348990` with
`7.0000e-01`. `variants.vcf` holds one record inside both rows' span,
`trailing_space`.

`BaseVepTabixPlugin.pm` lines 605 and 606 hand `parse_data` the line as the
library returns it, `LOEUF.pm` line 215 splits it on tabs and line 181 returns
the cell, so the space after the last cell is printed; line 192 reads such a
cell as a number and line 204 returns the gene's minimum as its cell is
written.

Ensembl VEP 116.3 wrote, on this fixture in the default format: under
`match_by=transcript`, `LOEUF=5.0000e-01` with the space after it on the row of
`ENST00000346798`, `LOEUF=7.0000e-01` on the row of `ENST00000348990` and no
LOEUF item on the record's other rows; under `match_by=gene`, the value with
its space on every row of the record. vep-rs writes the same.
