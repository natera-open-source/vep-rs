# ClinPred cells fixture

A ClinPred-shaped file of five columns (`#chr`, `Start`, `Ref`, `Alt`,
`ClinPred_Score`), bgzip-compressed and indexed with `tabix -s 1 -b 2 -e 2`,
the preparation `ClinPred.pm` describes (lines 36 to 44). Its four rows are fabricated, at
positions of the GRCh38 ClinPred fixture corpus: a row at `21:30282172` for
`T>C` whose fifth cell is empty and is followed by a sixth cell; a row at
`21:30282928` for `G>A` of five cells, the fifth empty (the line ends in a tab);
a row at `21:30396087` for `A>G` of four cells, with no score cell; and a row at
`21:30396091` for `C>T` scored `0.411650`. `variants.vcf` holds the four records.

`ClinPred.pm` line 127 splits a row into five variables, so Ensembl VEP keeps a
trailing empty cell as the empty string and leaves a missing cell undefined,
and lines 116 to 122 return the matched row's result either way.

Ensembl VEP 116.3 wrote, on this fixture: for the two rows whose fifth cell is
empty, `ClinPred=` in the default format, an empty cell in the tab format, an
empty CSQ field in the VCF format and `"clinpred": ""` in JSON; for the
four-cell row, `ClinPred=` in the default format, `-` in the tab format, an
empty CSQ field and no JSON key (`VEP_output.pm` lines 117 to 119; `Tab.pm`
line 128; `JSON.pm` line 327); for the scored row, the score in every format.
