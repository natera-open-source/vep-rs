# DosageSensitivity cells fixture

A DosageSensitivity-shaped table of three columns (`#gene`, `pHaplo`,
`pTriplo`), gzip-compressed, the file `DosageSensitivity.pm` reads whole. Its
four lines are fabricated, for genes of the GRCh38 DosageSensitivity corpus: a
line for `APP` of two cells, with no `pTriplo` cell; a line for `BACH1` whose
`pHaplo` cell is empty; a line for `SOD1` whose `pTriplo` cell is empty (the
line ends in a tab); and a line for `KRTAP20-4` of one cell, the gene alone.
`variants.vcf` holds one corpus record on each gene: `app_exon`, `bach1_exon`,
`sod1_exon2` and `krtap20_4_snv`; the `SOD1` and `KRTAP20-4` records also
overlap transcripts of neighbouring genes the table lacks, whose rows carry no
field in either engine.

`DosageSensitivity.pm` line 138 assigns the split line to three variables, so
Ensembl VEP keeps an empty cell, trailing or not, as the empty string and leaves a
missing cell undefined, and lines 141 to 142 store both keys for the gene.

Ensembl VEP 116.3 wrote, on this fixture: for a missing cell, `pTriplo=` in the
default format, `-` in the tab format, an empty CSQ field in the VCF format and
no JSON key; for an empty cell, `pHaplo=` in the default format, an empty tab
cell, an empty CSQ field and `"phaplo": ""` in JSON (`VEP_output.pm` lines 117
to 119; `Tab.pm` line 128; `JSON.pm` line 327); for a present score, the score
in every format.
