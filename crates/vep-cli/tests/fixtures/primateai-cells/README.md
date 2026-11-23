# PrimateAI cells fixture

A made-up PrimateAI-shaped score file, bgzip-compressed and indexed with
`tabix -s 1 -b 2 -e 2` as `PrimateAI.pm` describes (lines 57 to 69): eleven
comment lines, then the column header behind a `#` on line twelve. The header
is the one the tree's other made-up PrimateAI tables carry, which names the
three columns the module reads (`pos`, `alt`, `primateDL_score`, line 124),
with one more column, `note`, after `primateDL_score`, so a row's score cell is
not the last cell of its line.

Its three rows are fabricated, at positions of the GRCh38 PrimateAI fixture
corpus: a row at `chr21:19268555` and a row at `chr21:20204701` whose score
cell is empty and whose `note` cell is filled, and a row at `chr21:30282169`
scored `0.6125`. `variants.vcf` holds the three records, each a substitution
of one base to its row's alternate allele: `empty_score_intergenic`, which
overlaps no transcript in the corpus cache, `empty_score` and `scored`.

`PrimateAI.pm` line 188 splits a row with no limit and line 196 reads the score
at the column the header names, so an empty score cell with a filled cell after
it is a defined empty string: it passes the check of line 198 and lines 204 to
207 return it. The plugin's result prints on every row of the record, the
intergenic row included.

Ensembl VEP 116.3 wrote, on this fixture: for the two records whose score cell
is empty, `PrimateAI=` in the default format, an empty cell in the tab format,
an empty CSQ field in the VCF format and `"primateai": ""` in JSON
(`VEP_output.pm` lines 117 to 119; `Tab.pm` line 128; `JSON.pm` line 327); for
the scored record, the score in every format.
