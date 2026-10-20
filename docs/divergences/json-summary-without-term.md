# A JSON summary row with no consequence

A difference in whether a run completes, under two flags no measured dataset is run with; nothing is set aside.

Under `--most_severe` or `--summary`, Ensembl VEP's `summary_only` returns a record's hash whatever rows the record has, writing `WARNING: Unable to assign consequence type` when the record has none (`OutputFactory.pm` lines 454-458, the warning at line 455). A record with no consequence, a reference-homozygous individual's copy under `--individual` or a non-variant record under `--allow_non_variant`, therefore writes one row with no term in the default and tab formats. Ensembl VEP's JSON factory then dereferences the hash's absent consequence list and dies (`OutputFactory/JSON.pm` line 374), so a `--json` run over such a record ends with `Can't use an undefined value as an ARRAY reference`.

vep-rs writes the default and tab row and the warning as Ensembl VEP does. In JSON it writes the record's object with no consequence entry and `most_severe_consequence` `?`, the object it writes for the record without the two flags.
