# vep filter

`vep filter` filters an annotated output line by line, as Ensembl VEP's `filter_vep` does, over a VCF or tab output of `vep` or of Ensembl VEP. It takes the same options, read as Getopt::Long reads them: after one dash or two, in any case, and abbreviated to any prefix that names one option. A repeated option keeps its last value, and every `--filter` given is kept. The short forms are `-i`, `-o`, `-f`, `-s`, `-c`, `-l` and `-y`.

An abbreviation is resolved against the options of the release that reads it, so a later release that adds an option can give a prefix to that option or make it name two options, which is refused (exit 2). A script that must keep working writes each option's full name.

```bash
# The lines with a moderate-impact missense entry
vep filter -i annotated.vcf \
  --filter "Consequence is missense_variant and IMPACT is MODERATE"

# Only the passing entries, kept in the CSQ field
vep filter -i annotated.vcf --only_matched -f "SYMBOL is BRCA2" -o brca2.vcf

# The field names the input's header block offers
vep filter -i annotated.tab --list
```

A filter is `<field> <operator> <value>`. Several filters are joined by `and` and `or`, evaluated left to right with no precedence, grouped by parentheses and inverted by `not`; repeated `--filter` flags are joined by `and`. The operators are `is` (`eq`, `=`), `ne` (`!=`), `gt` (`>`), `lt` (`<`), `gte` (`>=`), `lte` (`<=`), `exists` (`ex`, `defined`), `nex`, `match` (`matches`, `regex`, `re`), `nre`, `in` and `is_child`. `match` is a case-insensitive, unanchored regular expression; `in` takes a comma-separated list or a file with one value per line.

A field alone tests that it is defined, `ne` without a value is `nex`, and a value starting with `#` names another field. A value shaped `text(number)` or `text:number` compares through its number against a numeric value, and against an `in` list when the list's first written entry is a number, where `filter_vep` picks a hash-order entry ([the class](../divergences/transcript-filter-in-list-numeric-by-hash-order.md)). Write an `in` list's entries in one shape, `HGNC:2559,HGNC:1132` or `2559,1132`, and both tools keep the same lines. A value holding several entries joined by `,` or `&` after a digit passes when any entry passes.

A field name resolves to a key of the line exactly, then in lower case, upper case, with a leading `_`, then as the unique case-insensitive prefix of one key; a prefix of several keys ends the run with an error.

A VCF line's fields are its columns, its INFO keys and, for each entry of the consequence field, that entry's fields; an INFO flag without `=` is defined and empty. The line passes when any entry passes. Every INFO item under the field's name is read, so an output written with [`--keep_csq`](input-output.md#output) is filtered on its earlier annotation and its new one, and `--only_matched` keeps each item's passing entries. A tab line's fields are its columns, with an `Extra` column exploded into its pairs; a word without `=` among them is undefined. An empty value is undefined, and so is `-` under every field but `Allele`; an INFO value on a line without a consequence entry stands as written.

The header block is written through, and a line passes unchanged unless `--only_matched` or `--soft_filter` rewrites it. The input is read as bytes and every line, header or data, is written as the bytes it was read from, so a line that is not UTF-8 is filtered as `filter_vep` filters it, each byte outside UTF-8 one character in a comparison.

A command-line word is read as bytes too, so a value, a field name or a file name outside UTF-8 compares byte for byte.

`match` runs its pattern over the bytes of the value as Perl runs it over a byte string: `.` is one byte, and `\w`, `\s`, `\d` and the case fold are ASCII. A `\p{}` class, a `\x{}`, `\o{}` or octal escape above 255 or a `\N{U+}` escape switches the pattern to Unicode rules over the bytes read as Latin-1, and a `u`, `a` or `aa` flag switches its own group, as in Perl. Under those rules the case fold is full: `ss` matches `ß` and `\x{DF}` matches `ss`.

`--list` sorts the names by byte. An error naming several keys of a line lists them as the line holds them, the columns, the INFO keys, then the consequence fields; `filter_vep` lists them in its hash's order, which changes from run to run.

`--ontology` rewrites `Consequence is <term>` to `Consequence is_child <term>`, which passes an entry whose term is the parent or lies below it in the Sequence Ontology. The answer comes from the built-in table of consequence terms rather than an Ensembl database, and the connection flags are refused by name. Where `vep filter` reads a line differently from `filter_vep` is stated on [Abbreviated fields, `is_child` and the regex dialect](../divergences/transcript-filter-fields.md#vep-filter).

Every error is one line `ERROR: [filter] <message>` on standard error and exit 1; a command-line syntax error is exit 2.

A line that standard error refuses, on a full device or a pipe without a reader, is dropped and changes no exit status. A run that ends with the warning of lines without the consequence field exits 0 with its output whole, and a refusal exits 1. `filter_vep` completes on a full device too, a refusal still ending non-zero; on a pipe without a reader it is killed by `SIGPIPE` (exit 141) at its warning or its error, with none of its output written.

## Input and output

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--input_file`, `-i` | String | none | The annotated output to filter, a VCF or tab file, read from standard input when absent; a `.gz` file, or any gzip stream, is decompressed |
| `--format` | String | none | The input format, `vcf` or `tab`, detected from the first data line when absent |
| `--gz` | bool | false | Read the input as gzip-compressed whatever its name |
| `--output_file`, `-o` | String | `stdout` | Where the filtered lines go; `stdout` (in any case) is standard output |
| `--force_overwrite` | bool | false | Overwrite an existing output file, which is otherwise an error |
| `--vcf_info_field` | String | `CSQ` | The INFO key the consequence annotations of a VCF input sit under |

## Filters

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--filter`, `-f` | String (repeatable) | none | A filter; repeated flags are joined by `and` |
| `--list`, `-l` | bool | false | Print the field names the input's header block names, sorted, and stop |
| `--count`, `-c` | bool | false | Print the number of passing lines alone |
| `--only_matched` | bool | false | In a VCF input, keep only the passing entries of the consequence field |
| `--soft_filter` | bool | false | In a VCF input, write every line and mark its FILTER column `filter_vep_pass` or `filter_vep_fail` |
| `--no_missing_consequence_warning` | bool | false | Do not warn about VCF lines without the consequence field |

## Run control

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--start`, `-s` | i64 | 1 | Skip the first N-1 passing lines |
| `--limit` | i64 | none | Write at most N passing lines |
| `--test` | i64 | none | Read at most N data lines |

## Ontology

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--ontology`, `-y` | bool | false | Match `Consequence is <term>` against the term and every Sequence Ontology term below it, from the built-in table |
| `--host` | String | none | Refused: Ensembl VEP's `filter_vep` connects to an Ensembl database for `--ontology`; the built-in table answers it here |
| `--user` | String | none | Refused: the database user of `--host` |
| `--pass`, `--password` | String | none | Refused: the database password of `--host` |
| `--port` | String | none | Refused: the database port of `--host` |
| `--version` | String | none | Refused: the Ensembl database release `--host` connects to |
| `--registry` | String | none | Refused: an Ensembl registry file naming the database |

