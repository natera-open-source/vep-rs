# AVI cells fixture

An AVI-shaped file of six columns (chromosome, position, reference allele,
alternate allele, raw score, PHRED score) with no header, bgzip-compressed and
indexed with `tabix -s 1 -b 2 -e 2`. Its four rows are fabricated, at positions
of the GRCh38 AVI corpus: a row at `chr21:5030192` for `A>G` whose sixth cell is
empty and is followed by a seventh cell; a row at `chr21:5030633` for `C>T` of
five cells, with no PHRED cell; a row at `chr21:5030664` for `A>T` scored
`1.446`; and a row at `chr21:5131441` for `T>C` of six cells, the sixth empty
(the line ends in a tab). `variants.vcf` holds the four records, each
intergenic in the corpus cache.

`AVI.pm` line 161 splits a row into six variables with a limit of seven, so
Perl keeps a trailing empty cell as the empty string, and line 163 returns
nothing for a row lacking its sixth cell, which is then never matched.

Ensembl VEP 116.2 wrote, on this fixture under `raw=1,links=1`: for the two
rows whose sixth cell is empty, `AVI_PHRED=` in the default format, an empty
cell in the tab format, an empty CSQ field in the VCF format and
`"avi_phred": ""` in JSON, with the raw score and the link beside it; for the
five-cell row, no AVI item in the default format, `-` in the tab format, an
empty CSQ field and no JSON key (`VEP_output.pm` lines 117 to 119; `Tab.pm`
line 128; `JSON.pm` line 327); for the scored row, the score in every format.
