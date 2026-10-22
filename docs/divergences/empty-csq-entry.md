# An empty `CSQ` entry in an input VCF's INFO

A difference in the text of a VCF output's INFO column, not in any tuple: the comparators score consequence sets and never read another INFO field, and nothing is set aside.

A VCF output rewrites each input record with the annotation appended to INFO, and an earlier annotation under the same key (`CSQ`, or the `--vcf_info_field` name) is removed first unless `--keep_csq` is given. Ensembl VEP removes it with one substitution over the INFO text: the key and its `=`, then one or more non-blank characters, up to the next `;` or the end of the field (`OutputFactory/VCF.pm` lines 327-331).

An entry with no value, `CSQ=`, gives those characters nothing of their own to match, so they take the `;` and the field after it: `CSQ=;DP=4` loses `DP=4`, and `DP=4;CSQ=;AF=0.5` loses `AF=0.5`. What is left of the field, with the new entry appended, is the output's INFO; a record with no consequence and nothing left prints `.`. An empty entry that ends the field, `DP=4;CSQ=`, matches nothing and stays, so the output carries two `CSQ` entries.

vep-rs splits INFO on `;`, removes the first entry under the key whatever its value, and keeps every other field. For `CSQ=;DP=4`, `DP=4;CSQ=;AF=0.5` and `DP=4;CSQ=` it writes `DP=4;CSQ=...`, `DP=4;AF=0.5;CSQ=...` and `DP=4;CSQ=...`. An entry with a value is removed alike by both engines: `CSQ=G|old;DP=4` becomes `DP=4;CSQ=...` in each.

Under `--keep_csq` neither engine removes anything: both write the INFO as read, `CSQ=;DP=4`, and append the new entry after a `;` when the record has consequences.

The records are `21 5030602 empty_csq A G . . CSQ=;DP=4`, intergenic on GRCh38, and the non-variant line `21 30643020 nonvar_empty TTG . . . CSQ=;DP=4`. Both engines run with `--vcf --allow_non_variant` against the release-116 GRCh38 cache. The INFO column each engine writes, the new entry's empty fields elided:

```
Ensembl VEP 116.2   empty_csq     CSQ=G|intergenic_variant|MODIFIER|...
vep-rs              empty_csq     DP=4;CSQ=G|intergenic_variant|MODIFIER|...
Ensembl VEP 116.2   nonvar_empty  .
vep-rs              nonvar_empty  DP=4
```
