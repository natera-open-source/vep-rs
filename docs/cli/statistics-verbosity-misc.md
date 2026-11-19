# Statistics, verbosity and miscellaneous

These flags cover the statistics report written beside the output, the logging level, and the flags that change how a record is handled.

## Statistics

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--no_stats` | bool | false | Write no statistics report and keep no counts for one, whatever `--stats_text`, `--stats_html` and `--stats_file` say; the annotation output is unchanged |
| `--stats_file`, `--sf` | String | none | Name of the statistics report; default \<output\_file\>\_summary with the report's extension (.txt under `--stats_text`, .html otherwise), so -o out.txt writes out.txt\_summary.html and -o STDOUT writes STDOUT\_summary.html in the working directory. A name ending .txt or .html takes the report's extension (x.html names x.txt for the text report); any other extension and a bare name gain it (x.foo names x.foo.html). An existing report is refused before the run unless `--force_overwrite`, and is checked before the output file |
| `--stats_text` | bool | false | Write the statistics report as plain text, \<output\_file\>\_summary.txt. It opens with the run's settings and times, the cache's data versions, and the counts of records read, annotated and filtered, of known variants, and of the genes, transcripts and regulatory features overlapped. The charts follow, each a key and count per line under its bracketed title: the variant classes, the consequences (most severe per record, all, coding), the SIFT and PolyPhen terms, the records per chromosome, their distribution along each chromosome when a FASTA gives the lengths, and the coding variants' position along their proteins in deciles. Together with `--stats_html` both reports are written |
| `--stats_html` | bool | false | Write the statistics report as an HTML page, \<output\_file\>\_summary.html, with the sections and charts of the text report; the default when neither `--stats_text` nor `--stats_html` is given |

## Verbosity

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--verbose`, `-v` | bool | false | Verbose output (debug-level logging) |
| `--quiet`, `-q` | bool | false | Suppress the status and progress lines and the standard-error copy of each warning; the warnings file is still written |
| `--no_progress` | bool | false | Suppress the progress line, `Processed N variants...` on standard error every 10,000 variants; the status lines, the summary line and the output are unchanged |

## Miscellaneous

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--dont_skip` | bool | false | Keep and annotate a record that fails a validation check (an unknown chromosome, an invalid allele string, a reference allele `--check_ref` rejects, a chromosome `--chr` leaves out) instead of skipping it; the skip warning is still written, except for a chromosome `--chr` leaves out, which is dropped without a warning and kept without one. Two records are skipped under the flag too, each with its warning, where Ensembl VEP keeps the feature: a record whose start is past its end plus one (Ensembl VEP's run then ends on the first transcript the record reaches), and a region, SPDI or HGVS transcript record whose span exceeds `--max_sv_size`, refused as its line is parsed. A line that cannot be parsed is skipped either way |
| `--allow_non_variant` | bool | false | Keep a VCF record whose first ALT is `.` (a reference-only call): the VCF output carries its line with any earlier CSQ removed unless `--keep_csq`, the JSON output writes its object with no consequence, and the default and tab outputs write no row for it; under `--individual` a copy whose genotype names a later ALT of the record is annotated as that allele |
| `--mirna` | bool | false | Add the miRNA column: the secondary-structure elements of a miRNA transcript the variant overlaps, as the SO terms miRNA\_stem and miRNA\_loop (a native cache carries the structure only when built with `--core-dumps`) |
| `--lrg` | bool | false | Refused with an error: `--lrg` is not supported: LRG mapping needs the Ensembl core database (Ensembl VEP: Cannot map to LRGs in offline mode) |
| `--tmpdir` | String | none | Accepted without effect: no module of Ensembl VEP reads it, and vep-rs writes no temporary file |

