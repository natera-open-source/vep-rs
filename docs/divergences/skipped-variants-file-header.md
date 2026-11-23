# The header and record column of `--skipped_variants_file`

A difference in the text of the skipped-variants file and of the skip warnings, not in any row of the annotation. Four things differ: the header's name for the input, the width of the record column on a line with a multi-byte character, the cut of a quoted line whose byte 46 falls inside one, and the blanks around a VCF line read under `--delimiter`. Ensembl VEP names a file handle under `--input_data`, counts the column and the cut in bytes and collapses the raw line. vep-rs names `STDIN` there, counts 40 characters, cuts at a character boundary and re-joins VCF fields with tabs first. Nothing is set aside, since the comparators never read the file or the warnings.

Ensembl VEP writes the file when the run ends (`Runner.pm` line 176): a header naming the input, then one line per skipped record (`BaseRunner.pm` lines 371-376). The header prints the `input_file` parameter as typed: `STDIN` by default (`Config.pm` line 301), `stdin` for `-i stdin`, and under `--input_data` the in-memory file handle the string is read through, which prints as `*Bio::EnsEMBL::VEP::Runner::IN` (`Runner.pm` lines 818-820).

The record column is `%-40.40s` of the skipped line (`BaseRunner.pm` line 375), after its whitespace runs are collapsed to one space and a line over 50 bytes is cut to its first 46 bytes plus `...` (`BaseVEP.pm` lines 849-850). The skip warning quotes that cut line in full. The input is opened with no decoding layer (`ensembl-io`, `TextParser.pm` line 70), so `length`, `substr` and `%-40.40s` count bytes. A line carrying a two-byte character such as `é` is cut one character earlier, or padded one space less, than its character count gives.

When byte 46 falls inside a character, the cut splits it and the warning quotes the character's lone lead byte, `\xC3` for `é`, before the `...`: a byte sequence that is not UTF-8.

vep-rs writes the same file in the same layout, with the same collapse, and four things differ. The quoted line is cut by bytes too, at the last character boundary at or before byte 46, so a line whose byte 46 falls inside a character loses that whole character and the stray lead byte is never written. The header names the input as `-i` spelled it, `STDIN` by default, and under `--input_data` it names `STDIN`, or the `--input_file` value when one stands beside the string.

The record column is the line's first 40 characters, padded to 40 characters, so a line with a multi-byte character pads to a different width. On a line of ASCII characters under `-i <path>`, `-i STDIN` in any letter case or no `-i`, the two files are byte-identical.

Under a `--delimiter` other than a tab, the skip warning and the record column quote a VCF line's fields re-joined with tabs before the collapse. The delimiter then prints as a space, and a blank the split drops at either end of the line is lost. An Ensembl-format line is quoted as read, as Ensembl VEP quotes it. Ensembl VEP collapses the raw line's whitespace runs and keeps every other character (`BaseVEP.pm` line 849), so a VCF line opening with a space under `--delimiter ' '` quotes as `( ZZ 30282169 ...)` there and `(ZZ 30282169 ...)` here.

The records below are a sequence record whose alleles are no bases, and a record on a chromosome the cache lacks whose identifier carries an `é` and runs past the column:

```
21	30165901	invalid_alleles	X	Y	.	.	.
Z	1000	idé_accent_then_a_long_identifier_tail	A	G	.	.	K=v
```

Under `--format vcf --skipped_variants_file skipped.txt`, the header each engine writes for the input given three ways:

| Input | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `-i in.vcf` | `from in.vcf` | `from in.vcf` |
| `-i stdin`, the file piped | `from stdin` | `from stdin` |
| `--input_data` with the first record | `from *Bio::EnsEMBL::VEP::Runner::IN` | `from STDIN` |

The first record's line is identical in both files: `Line 1`, the record padded to 40, then `Invalid allele string X/Y or possible parsing error`. The second record's line reads `Line 2`, the record column below, then `Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features`:

| Engine | Record column | Width |
| --- | --- | --- |
| Ensembl VEP | `Z 1000 idé_accent_then_a_long_identifie` | 40 bytes, 39 characters |
| vep-rs | `Z 1000 idé_accent_then_a_long_identifier` | 40 characters, 41 bytes |

To reproduce: Ensembl VEP 116.3 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--skipped_variants_file skipped.txt`; vep-rs with the same flags over the same cache. Run each engine once with the file, once with the file piped under `-i stdin`, and once with the first record as `--input_data`.
