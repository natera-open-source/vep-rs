# dbNSFP cells fixture

A dbNSFP-shaped file of thirteen columns (`#chr`, `pos(1-based)`, `ref`, `alt`,
`aaref`, `aaalt`, `rs_dbSNP`, `hg19_chr`, `hg19_pos(1-based)`,
`Ensembl_transcriptid`, `SIFT_score`, `SIFT_pred`, `GERP++_RS`),
bgzip-compressed and indexed with `tabix -s 1 -b 2 -e 2`, the GRCh38
preparation of the plugin's documentation (`dbNSFP.pm` lines 69 to 70). Its two
rows are fabricated: a row at `1:1045000`, inside the
region `dbNSFP.pm` reads to find the `;`-joined columns (lines 444 to 483), whose `SIFT_pred`
holds one value so that column is not transcript-specific; and a row at
`21:33296341` for the `A>G` change `Q/R` naming two transcripts, whose
`SIFT_score` list is `.;0.297` and whose `SIFT_pred` cell is empty.
`variants.vcf` is that one record.

Under `transcript_match=1` the first transcript's `SIFT_score` entry is `.`,
which Ensembl VEP's plugin leaves out, and the empty `SIFT_pred` cell is a
value the row has: `SIFT_pred=` in the default format, an empty cell in the tab
format and `""` in JSON (`dbNSFP.pm` lines 388, 407 and 408; `Tab.pm` line
128; `JSON.pm` line 327).
