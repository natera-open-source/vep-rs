# An INFO column of exactly `0` in an input VCF

A difference in a VCF output's INFO column for an input record whose INFO column is the one character `0`, not in any tuple; nothing is set aside. Ensembl VEP joins its annotation onto the `0` with no separator, and writes `.` in the column's place when the record has no consequence. vep-rs keeps the `0` in both cases and writes `;` between it and the annotation.

A VCF output rewrites each input record with the annotation appended to its INFO column. Ensembl VEP appends in two branches (`OutputFactory/VCF.pm` lines 347-353). For a record with consequences, line 348 adds a `;` to the column if the column's text is true, and line 349 adds the entry. For a record with none, line 352 replaces the column with `.` unless its text is true.

Perl reads the string `0` as false, as it reads the empty string. So both branches treat a column that holds `0` as a column that holds nothing: the first adds no `;`, and the second writes `.` over it.

vep-rs asks whether the column holds any text. Both engines read a `.` as an empty column (line 323). vep-rs keeps any other text, `0` included, and writes `;` and the new entry after it when the record has one.

The records are `21 5030602 with_row A G . . 0`, intergenic on GRCh38, and `21 5030610 no_row T * . . 0`. The one ALT of the second is `*`, so neither engine writes a consequence for it. Both engines run with `--vcf` against the release-116 GRCh38 cache. The INFO column each engine writes, the new entry's empty fields elided:

```
Ensembl VEP 116.3   with_row  0CSQ=G|intergenic_variant|MODIFIER|...
vep-rs              with_row  0;CSQ=G|intergenic_variant|MODIFIER|...
Ensembl VEP 116.3   no_row    .
vep-rs              no_row    0
```

The input column names one key, `0`, with no value. Ensembl VEP's first line names one key too, `0CSQ`: the input's key is gone, and the record has no entry under `CSQ`, the name the output's `##INFO` line declares. A reader that looks up `CSQ` finds no annotation on the record. Its second line has lost the record's own INFO.

vep-rs's lines hold the input's key as written, and the first holds the annotation under its declared name. Every other column of the two lines agrees byte for byte in the two engines.

`0` is the only text of one or more characters that Perl reads as false, so the difference stops at this one value. The columns `00`, `0.0` and `0;DP=4` are kept by both engines, each with `;` before the new entry.

To reproduce: Ensembl VEP 116.3 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--vcf`, over a VCF of the two records above; vep-rs with the same flags over the same cache.
