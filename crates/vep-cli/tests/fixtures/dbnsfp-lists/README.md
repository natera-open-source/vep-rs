# dbNSFP list fixture

A dbNSFP-shaped file of twelve columns (`#chr`, `pos(1-based)`, `ref`, `alt`,
`aaref`, `aaalt`, `rs_dbSNP`, `hg19_chr`, `hg19_pos(1-based)`,
`Ensembl_transcriptid`, `SIFT_score`, `SIFT_pred`), bgzip-compressed and indexed
with `tabix -s 1 -b 2 -e 2`, the GRCh38 preparation of the plugin's
documentation (`dbNSFP.pm` lines 69 to 70). Its two rows are
fabricated: a row at `1:1045000`, inside the region `dbNSFP.pm` reads to find
the `;`-joined columns (lines 444 to 483), and a row at `21:33296341` for the `A>G` change `Q/R`
naming six transcripts whose `SIFT_score` list holds two entries and whose
`SIFT_pred` list holds six. `variants.vcf` is that one record.

Under `transcript_match=1` the six-entry list is read per transcript and the
two-entry list is malformed: Ensembl VEP's plugin leaves the whole value and
reports every transcript row on standard error (`dbNSFP.pm` lines 378-385).
