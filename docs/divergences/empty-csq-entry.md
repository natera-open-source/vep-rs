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

## The same substitution under `--vcf_info_field`

Under `--vcf_info_field <name>` the entry removed is the one under that name, through the same substitution (`OutputFactory/VCF.pm` lines 312 and 327-331). So an empty entry under any name takes the field after it in Ensembl VEP. Under `--vcf_info_field ANN` the INFO `DP=13;ANN=;CSQ=A|old_csq,C|old_csq;AF=0.5` loses its `CSQ` entry, another key that the removal otherwise leaves in place. An empty `ANN=` that ends the field matches nothing and stays, and the new entry is appended after it, so the output carries two `ANN` keys. vep-rs removes the empty entry alone in both shapes, as it does under the default name.

The record is the first of the `GRCh38-writers` corpus, `21 10541165 ENST00000618007_multi G A,C`, whose `A` allele is a missense variant of TPTE. Both engines run with `--vcf --vcf_info_field ANN` against the release-116 GRCh38 cache. The INFO column each engine writes for each input INFO, the new entry's fields elided:

```
input        DP=13;ANN=;CSQ=A|old_csq,C|old_csq;AF=0.5
Ensembl VEP  DP=13;AF=0.5;ANN=A|missense_variant...
vep-rs       DP=13;CSQ=A|old_csq,C|old_csq;AF=0.5;ANN=A|missense_variant...
input        DP=13;CSQ=A|old_csq;ANN=
Ensembl VEP  DP=13;CSQ=A|old_csq;ANN=;ANN=A|missense_variant...
vep-rs       DP=13;CSQ=A|old_csq;ANN=A|missense_variant...
```

A `;` that ends the column right after the entry under the name goes with the entry in both engines, since the one or more non-blank characters take it: `DP=87;ANN=T|old;` becomes `DP=87;ANN=...`. A `;` ending the column anywhere else stays, so `ANN=x;DP=87;` becomes `DP=87;;ANN=...` in both.

## The name read as a regular expression

Ensembl VEP interpolates the field name unescaped into the header test (`OutputFactory/VCF.pm` line 221) and the substitution over INFO (lines 328-329), so a name holding a metacharacter matches other keys. Under `--vcf_info_field A.N` the input lines `##INFO=<ID=ANN,...>` and `##contig=<ID=ANN,...>` are dropped and the entry `ANN=G|old_ann` is removed, as if the name were `ANN`. vep-rs matches the name literally: those lines and that entry stay, and the new entry is appended under `A.N`.

A dot is the one metacharacter a VCF INFO key admits. A value that is not a key (a letter or underscore, then letters, digits, underscores or dots, or the reserved `1000G`), the empty string and `0` included, is an argument syntax error here (exit 2). Ensembl VEP accepts it, declares the given name in the `##INFO` line and, for a name Perl reads as false, keys the entry `CSQ` (lines 219 and 312).

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--vcf --allow_non_variant` for the first two records and `--vcf --vcf_info_field ANN` for the third. The dotted name takes `--vcf --vcf_info_field A.N`; vep-rs runs with the same flags over the same cache.
