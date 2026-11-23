# Shared plugin key fixture

Four made-up tables for two pairs of plugins that write one output key, each
bgzip-compressed and indexed with `tabix -s 1 -b 2 -e 2`, and `variants.vcf`,
one record at `21:10521702` of the GRCh38 TranscriptAnnotator corpus, a
substitution of one base that overlaps nine transcripts in the corpus cache.
Every table holds one fabricated row, for that record:

- `BayesDel_shared_key_synthetic.txt.gz`: five columns (`#Chr`, `Pos`, `Ref`,
  `Alt`, `BayesDel_addAF`); the row's fifth cell is present and empty (the line
  ends in a tab).
- `transcript_annotator_BayesDel.tsv.gz`: a TranscriptAnnotator table with one
  value column, headed `BayesDel`; the row names `ENST00000618007` and ends
  after its five leading cells.
- `CADD_shared_key_GRCh38_synthetic.tsv.gz`: a CADD-shaped table (`#Chrom`,
  `Pos`, `Ref`, `Alt`, `RawScore`, `PHRED`); the row's PHRED cell is `23.1`.
- `transcript_annotator_PHRED.tsv.gz`: a TranscriptAnnotator table with one
  value column, headed `PHRED`; the row names `ENST00000618007` and its cell is
  `7.5`.

Ensembl VEP runs its plugins on each output line in load order and assigns
every key of a plugin's result over the line's hash, an undefined value
included (`run_plugins`, `OutputFactory.pm` lines 2239 to 2241), so under a key
two plugins write, the plugin loaded later stands on the line.

Ensembl VEP 116.3 wrote, on this fixture:

- Under `--plugin BayesDel,file=<scores>` then
  `--plugin TranscriptAnnotator,file=<BayesDel table>,prefix=0`: on the row of
  `ENST00000618007`, `BayesDel=` in the default format, `-` in the tab format,
  an empty CSQ field and no JSON key, TranscriptAnnotator's undefined value
  (`TranscriptAnnotator.pm` lines 245 to 253) over BayesDel's empty string
  (`BayesDel.pm` line 151); on the other eight rows, `BayesDel=`, an empty tab
  cell, an empty CSQ field and `"bayesdel": ""`.
- Under `--plugin TranscriptAnnotator,file=<PHRED table>,prefix=CADD` then
  `--plugin CADD,<scores>`: `23.1` under `CADD_PHRED` on all nine rows, in
  every format.
- Under the same two flags with CADD first: `7.5` on the row of
  `ENST00000618007` and `23.1` on the other eight, in every format.

Its tab header and its CSQ format list the key once per plugin, both columns
holding the row's one value; vep-rs lists the key once
(`duplicate_output_column_collapsed` in `docs/intended-divergences.md`).
