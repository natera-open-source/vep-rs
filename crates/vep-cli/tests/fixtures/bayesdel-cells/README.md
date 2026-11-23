# BayesDel cells fixture

A BayesDel-shaped file of five columns (`#Chr`, `Pos`, `Ref`, `Alt`,
`BayesDel_addAF`), bgzip-compressed and indexed with `tabix -s 1 -b 2 -e 2`,
the preparation `BayesDel.pm` describes (lines 54 to 62). Its four rows are fabricated, at
positions of the GRCh37 BayesDel fixture corpus: a row at `21:26078748` for
`T>C`, an intergenic record, whose fifth cell is empty and is followed by a
sixth cell; a row at `21:36252953` for `C>T` of four cells, with no score
cell; a row at `21:40762624` for `A>T` scored `-0.55399`; and a row at
`21:43505455` for `T>G` of five cells, the fifth empty (the line ends in a
tab). `variants.vcf` holds the four records.

`BayesDel.pm` line 151 splits a row into five variables, so Ensembl VEP keeps a
trailing empty cell as the empty string and leaves a missing cell undefined,
and line 146 returns the matched row's result either way; the plugin's result
prints on every row of the record, the intergenic row included (line 117).

Ensembl VEP 116.3 wrote, on this fixture: for the two rows whose fifth cell is
empty, `BayesDel=` in the default format, an empty cell in the tab format, an
empty CSQ field in the VCF format and `"bayesdel": ""` in JSON; for the
four-cell row, `BayesDel=` in the default format, `-` in the tab format, an
empty CSQ field and no JSON key (`VEP_output.pm` lines 117 to 119; `Tab.pm`
line 128; `JSON.pm` line 327); for the scored row, the score in every format.
