# Statistics, verbosity and miscellaneous

These flags cover the statistics file Ensembl VEP writes and vep-rs does not, the logging level, and the flags that change how a record is handled.

## Statistics

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--no_stats` | bool | false | Accepted for VEP compatibility; not implemented (no statistics file is written, with or without it) |
| `--stats_file`, `--sf` | String | none | Accepted for VEP compatibility; not implemented (no statistics file is written) |
| `--stats_text` | bool | false | Accepted for VEP compatibility; not implemented (no statistics file is written) |

## Verbosity

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--verbose`, `-v` | bool | false | Verbose output (debug-level logging) |
| `--quiet`, `-q` | bool | false | Suppress the status and progress lines and the standard-error copy of each warning; the warnings file is still written |

## Miscellaneous

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--dont_skip` | bool | false | Keep and annotate a record that fails a validation check (an unknown chromosome, an invalid allele string, a reference allele `--check_ref` rejects, a chromosome `--chr` leaves out) instead of skipping it; the skip warning is still written, except for a chromosome `--chr` leaves out, which is dropped without a warning and kept without one. Two records are skipped under the flag too, each with its warning, where Ensembl VEP keeps the feature: a record whose start is past its end plus one (Ensembl VEP's run then ends on the first transcript the record reaches), and a region, SPDI or HGVS transcript record whose span exceeds `--max_sv_size`, refused as its line is parsed. A line that cannot be parsed is skipped either way |
| `--allow_non_variant` | bool | false | Keep a VCF record whose first ALT is `.` (a reference-only call): the VCF output carries its line with any earlier CSQ removed unless `--keep_csq`, the JSON output writes its object with no consequence, and the default and tab outputs write no row for it; under `--individual` a copy whose genotype names a later ALT of the record is annotated as that allele |
| `--mirna` | bool | false | Add the miRNA column: the secondary-structure elements of a miRNA transcript the variant overlaps, as the SO terms miRNA\_stem and miRNA\_loop (a native cache carries the structure only when built with `--core-dumps`) |
| `--lrg` | bool | false | Refused with an error: `--lrg` is not supported: LRG mapping needs the Ensembl core database (Ensembl VEP: Cannot map to LRGs in offline mode) |
| `--tmpdir` | String | none | Accepted without effect: no module of Ensembl VEP reads it, and vep-rs writes its Parquet intermediate beside the output file |

