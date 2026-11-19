# A filter unit whose first word is `0`

**The defect in one sentence.** Ensembl VEP 116.2 tests whether a filter unit has its field yet with Perl's truth test, and the string `0` is false to Perl. So a unit whose first word is `0` loses that word: the operator word becomes the field, and the next word, when there is one, is read as the operator, refused unless it names one. vep-rs reads the first word as the unit's field whatever it is, so `0` names a field undefined everywhere and a unit with a value keeps nothing.

**Why this is a defect.** The parser fills a unit's three slots in order, deciding which slot a word fills by testing whether the earlier slots are set. The field slot is tested as `!$current->{field}`, true for an empty slot and true again for one holding `0`, so the word after `0` overwrites it. In `0 is 1` the field becomes `is`, and `1` is read as the operator and refused as `No such operator "1"`, an error naming the value. In `0 ne` the field becomes `ne`, and the unit tests whether a field named `ne` is defined. In `0 biotype is protein_coding` every word lands one slot late, and the unit runs as `biotype is protein_coding`.

The value slot is tested with `defined`, so only the field slot loses a `0`.

**What vep-rs does instead.** vep-rs reads the first word of a unit as its field whatever the word is, so `0` names a field called `0`, which no transcript key and no output column is named or begins with, undefined everywhere. `0 is 1` and `0 in <list>` keep nothing, since a unit with a value fails on an undefined field, `0 ne` and `0 nex` keep everything, and `0 biotype is protein_coding` is refused, since `biotype` names no operator.

`0 is`, with no value, compares the undefined field with no value and keeps everything, as `filter_eq` compares two undefined values. The same rule governs `vep filter`, whose `filter_vep` counterpart shares the grammar.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/FilterSet.pm`):

- Line 223, `parse_filters`: the words of a filter string, split at blanks and parentheses.
- Lines 238-255: an empty word is skipped, `and` and `or` close the unit, `not` inverts it.
- Lines 260-262: the field slot is filled while `!$current->{field}` holds; `0` is false to Perl, so the next word fills the slot again.
- Lines 265-281: the operator slot is filled while `!$current->{operator}` holds; a word that names no predicate and no synonym dies with `No such operator`.
- Lines 284-286: the value slot is filled while `!defined($current->{value})` holds, so a value spelled `0` stays.
- Line 374, `finish_filter_node`: a unit without an operator tests that its field is defined.
- Lines 481-483, `evaluate`, and 538-572, `get_input`: a field no key matches exactly, in another case, with a leading `_` or as a prefix is undefined, and a unit with a value fails on it.

**What correcting it would change.** With the field slot tested as `!defined($current->{field})`, as the value slot is, `0 is 1` would name the field `0`, undefined on every transcript, and keep nothing; `0 ne` would keep everything. A unit whose first word is anything but `0` reads the same either way.

**The class it produces.** `transcript_filter_zero_field_read_as_unset`: the outcome of a `--transcript_filter` or `vep filter` run whose filter holds a unit whose first word is `0`. Ensembl VEP reads the unit's words one slot late, refusing `0 is 1`, running `0 biotype is protein_coding` as `biotype is protein_coding` and testing `0 ne` as a field named `ne`; vep-rs reads `0` as the field. The class is documented only: no committed flag set and no committed tool run passes such a unit.

**How often the shape occurs.** Tier 3: the shape is a value of `--transcript_filter`, or of `vep filter`'s `--filter`, that no ordinary command line carries. No ordinary field is named `0` either. A `--custom` source given `short_name=0` is named after its file on both engines, since `AnnotationSource/File.pm:154-157` falls back to the file's last path segment for a short name false to Perl, and a VCF INFO key may not begin with a digit. So the shape arises only in a filter string as typed. Of the 392 committed flag sets, 9 pass `--transcript_filter` and none passes a unit whose first word is `0`; none of the 28 committed `vep filter` runs does.

## 1. A unit `0 is 1` and a unit `0 ne` over 41 transcripts

**In short.** The record `tf_17608184`, `21 17608184 A>G`, lies within reach of 41 transcripts. Under `--transcript_filter '0 is 1'`, Ensembl VEP 116.2 refuses the run with `ERROR: No such operator "1"` and writes nothing; vep-rs runs and writes the record's one `intergenic_variant` row, every transcript dropped. Under `--transcript_filter '0 ne'`, Ensembl VEP writes that intergenic row and vep-rs the 41 transcript rows.

### 1.1 What a reader needs to know first

`--transcript_filter` keeps the transcripts a filter passes before any consequence is computed; a record whose every transcript is dropped is written as `intergenic_variant`. A filter string is units of `field operator value` joined by `and` and `or`. A unit without an operator tests that its field is defined, and `ne` without a value tests that it is not. A field name resolves to a transcript key exactly, in another case, with a leading `_`, or as the unique prefix of one key; a name no key matches is undefined, and a unit with a value fails on an undefined field.

### 1.2 The record

The input is the `tests/golden/116/GRCh38-transcript-filter` corpus's own record:

```
21	17608184	tf_17608184	A	G	.	.	.
```

Both engines run over the GRCh38 release 116 cache, with the filter quoted below added.

### 1.3 What each engine writes

Ensembl VEP 116.2 with `--transcript_filter '0 is 1'` writes no output and ends with exit status 2 and this line in its exception block:

```
MSG: ERROR: No such operator "1"
```

vep-rs with the same filter writes one row:

```
tf_17608184	21:17608184	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

With `0 in ENST00000284971,ENST00000400454`, Ensembl VEP ends with `No such operator "ENST00000284971,ENST00000400454"`, and with `0 in <file>` it names the file's path as the operator; vep-rs writes the intergenic row for either.

Ensembl VEP with `--transcript_filter '0 ne'` writes the intergenic row above. vep-rs writes the 41 transcript rows of the unfiltered run, the first two:

```
tf_17608184	21:17608184	G	ENSG00000154640	ENST00000339775	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
tf_17608184	21:17608184	G	ENSG00000154640	ENST00000348354	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
```

The corpus holds 32 records, 952 transcript rows unfiltered, and the counts below hold in each of the four output formats. Under `0 is 1` and the two `0 in` shapes, Ensembl VEP refuses and vep-rs writes 35 intergenic rows, one per allele. Under `0 ne`, Ensembl VEP writes the 35 intergenic rows and vep-rs the 952. Under `0 exists` and under `0` alone, both write the 35 intergenic rows. Under `0 biotype is protein_coding`, Ensembl VEP writes the 365 rows of `biotype is protein_coding`, byte for byte, and vep-rs refuses the run with `ERROR: [filter] No such operator "biotype"`.

Under `0 is`, Ensembl VEP reads `is` as a field, the prefix of `is_canonical`, and writes 108 rows: the 109 rows of `is_canonical` less the first transcript the prefix resolves on, which is [the abbreviated-fields class](transcript-filter-fields.md). vep-rs writes the 952.

Through the tools over the unfiltered output, `filter_vep -filter '0 is 1'` and `-filter '0 in <list>'` refuse with the same messages, and `-filter '0 is'` keeps no line, since no column is named `is`. `vep filter` keeps no line under `0 is 1` or `0 in <list>` and every line under `0 is`.

### 1.4 Why Ensembl VEP writes that

The tokenizer yields `0`, `is` and `1`. The field slot takes `0`. For `is`, the test `!$current->{field}` is `!"0"`, true, so `is` overwrites the field. For `1`, the field is set and the operator slot is empty, and `1` names no predicate and no synonym: the run dies with `No such operator "1"`. In `0 biotype is protein_coding` the word after `0` is `biotype`, which becomes the field, then `is` names an operator and `protein_coding` is the value, so the unit runs.

For `0 ne`, the unit ends with the field `ne` and no operator, which `finish_filter_node` makes `ex`. `get_input` finds no transcript key named `ne` in any case, none named `_ne` and none beginning with `ne`, so the field is undefined and the test fails on every transcript. For `0 is`, the field `is` is a prefix of the one key `is_canonical`, so the unit keeps the canonical transcripts.

### 1.5 Why vep-rs's row is the correct one

The grammar names the first word of a unit as its field, and nothing in it exempts a word spelled `0`. Ensembl VEP's own value slot keeps a `0`, so `biotype is 0` compares with `0` on both engines, and only the field slot loses the word. The refusal names `1`, a word written in the value's place, as an operator; the two-word unit tests a field named after an operator word the user never meant as a field. vep-rs reads every unit in the documented order and treats `0` as a field name like any other, undefined on every transcript, so a unit with a value keeps nothing and a negated existence test keeps everything.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --transcript_filter '0 is 1'
```

vep-rs writes the record's intergenic row. Ensembl VEP 116.2 with the same filter over the GRCh38 release 116 cache ends with `ERROR: No such operator "1"`. With `--transcript_filter '0 ne'`, vep-rs writes the 41 transcript rows and Ensembl VEP the intergenic row.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/transcript_filter_zero_field_read_as_unset.toml` holds the record, the filter `0 is 1`, Ensembl VEP's refusal and exit status, and vep-rs's clean run. The class covers the outcome of a `--transcript_filter` or `vep filter` run whose filter holds a unit whose first word is `0`. No committed flag set or tool run reaches the shape, so no manifest names the class and the class is documented only.
