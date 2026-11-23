# A `--custom` key that spells a field Ensembl VEP holds as a plain value

**The defect in one sentence.** Ensembl VEP 116.3 files a `--custom` source's values under the source's key by pushing them onto a list. Where that key spells a flag field it has already set as one plain value, such as `VARIANT_CLASS` under `--variant_class` or `IMPACT`, the push fails and the run ends on a Perl error. vep-rs writes every row and, under the key, the field's value followed by the source's values. The push fails the same way under a key that holds a statistic of a source, a plain number; there vep-rs writes the statistic, then the values of the source records that name the row's allele.

Four terms recur. A source's key is its short name, under which the names of its records print, or its short name, `_` and a field it asks for. A flag field is a field of the output beside the thirteen fixed columns, most of them added by a flag. A plain value is one value that is not a list. A push adds a value to the end of the list under a key. Fourteen flag fields are plain values when a push is made; [section 2](#2-the-fourteen-fields) lists them.

A source that asks for `summary_stats` also adds a statistic, a figure over the records it matches at an input record, such as their `count`. Ensembl VEP holds it as a plain value too ([section 3](#3-a-key-that-spells-a-sources-statistic)).

**Why this is a defect.** A source's short name is any string without commas, and Ensembl VEP accepts each of these names: the run starts, writes its header and annotates. It dies at the first push onto a plain value, on a message of the Perl interpreter that names neither the option nor the input record. The same command completes in JSON and with the field's flag off.

Ensembl VEP's own rows show what a push means where the field is a list. On the input record of section 1, under the source of that section named `Consequence`, its row of the allele `C` holds the row's term and then the source record's name: the `Consequence` column reads `intergenic_variant,exact_2`.

**What vep-rs does instead.** The key holds both values. vep-rs writes the field's value, then each value Ensembl VEP pushes after it, in the order of its pushes. The default and tab formats join them with a comma and a `CSQ` entry with `&`. The key heads one tab column and one `CSQ` field ([A key listed twice by the tab and VCF writers](duplicate-output-column.md)). A key that one source declares as a field and as a statistic has one description line in the tab and VCF headers, which carries both texts ([section 3.6](#36-the-description-of-a-key-one-source-declares-twice)).

Ensembl VEP makes two pushes. The first takes the source records that name no allele and belongs to the input record. The second is made on each row and takes the source records that name the row's allele. A field set before both pushes is followed by the values of both. A field set between them replaces the values of the first push and is followed by the values of the second alone.

A source record that names an allele prints on the rows of that allele ([Custom annotation](../custom-annotation.md#1-an-exact-custom-record-repeated-once-per-row-beside-a-whole-variant-record-of-the-same-name)). In JSON both engines write the field's own value under its key and the source's records under `custom_annotations`. Where Ensembl VEP completes a run under a key that spells one of the fourteen fields, vep-rs writes Ensembl VEP's cells.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/OutputFactory.pm` unless a path says otherwise):

- Line 65: `use strict`.
- Lines 1297-1312, `_add_custom_annotations_to_hash`: line 1301 pushes a record's name onto the list under the short name, and line 1303 a field's value onto the list under the short name, `_` and the field. A key that holds a plain string is read as the name of a list, which `strict` forbids.
- Lines 883-1025, `VariationFeature_to_output_hash`: the input record's fields are set at lines 894-999, the first push is made at lines 1004-1014, and `CHECK_REF` is set at line 1020.
- Line 391: each row starts from a copy of the input record's fields.
- Lines 1333-1446, `VariationFeatureOverlapAllele_to_output_hash`: the row's fields are assigned at lines 1341-1410 and the second push is made at lines 1426-1435.
- `OutputFactory/JSON.pm`, lines 249-258: the JSON writer takes the custom annotations off the input record before any field is set, so it makes no push.

**What correcting it would change.** A push that kept a plain value as the first item of the list would complete the run and write the field's value before the source's, as vep-rs does.

**The class it produces.** `custom_key_named_as_a_plain_value_field`: every row of a run in the default, tab or VCF format from the first push onto a plain value. The reference, Ensembl VEP's side of a comparison, ends there, with exit status 255 under `--fork` and 25 without it. The candidate, vep-rs's side, writes every row. No comparator codes the class: the reference holds no such row to compare, and no committed run meets the shape. In a run the reference completes, the class holds the header lines and the second column of a key that one source declares twice ([section 3.6](#36-the-description-of-a-key-one-source-declares-twice)).

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a `--custom` source whose key spells one of the fourteen fields, a statistic Ensembl VEP assigns a source or one of the three columns of section 4. Of the 511 committed flag sets, 19 pass `--custom`, and none names such a source. No record of the cuts and no golden corpus carries the shape, which is a property of the command.

## 1. A source named `VARIANT` under `--variant_class`

**In short.** The input record is a substitution with two alternate alleles and no transcript near it. A VCF source named `VARIANT` that asks for its `CLASS` field declares the key `VARIANT_CLASS`, the name of the field `--variant_class` adds. The source holds a record on the first alternate allele. Ensembl VEP 116.3 dies with `Can't use string ("SNV") as an ARRAY ref`, its header written and no row. vep-rs writes both rows, `SNV,exact` on the first allele's and `SNV` on the other's.

### 1.1 What a reader needs to know first

`--custom` annotates each input record from a file of the user's own. Its `short_name` names the source and makes its keys: `VARIANT` for the names of its records, which for a VCF source are their identifiers, and `VARIANT_CLASS` for their `CLASS` field. A VCF source of `type=exact` matches an input record's alleles one by one unless either record is structural, and a record so matched names one allele. Any other match is by position; its record names none. [Custom annotation](../custom-annotation.md) describes the option.

`--variant_class` adds the field `VARIANT_CLASS`, the Sequence Ontology class of the input record, `SNV` for a substitution of one base.

### 1.2 The record

The input record is line 9 of the input of the `GRCh38-custom` golden corpus, read alone against that corpus's cache (`tests/golden/116/GRCh38-custom/`):

```
21	5380280	rs1444031311	A	C,G	.	.	.
```

The source is the corpus's `custom/sites.vcf.gz`, whose record at the position is named `exact_2`, carries `CLASS=exact` and names the allele `C`. Both engines run offline with:

```
--custom file=<corpus>/custom/sites.vcf.gz,short_name=VARIANT,format=vcf,type=exact,fields=AF%CLASS
--variant_class --tab
```

### 1.3 What each engine writes

Ensembl VEP 116.3 writes its header and no row, and this line to standard error:

```
Can't use string ("SNV") as an ARRAY ref while "strict refs" in use at /opt/vep/src/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm line 1303, <__ANONIO__> line 3.
```

Run with `--fork 1`, as the reproduction runs it, the line `Died in forked process` and a process number follow, and the exit status is 255. Run without `--fork` it writes the line alone and exits with status 25.

vep-rs writes the two rows. Their last eight columns are `IMPACT`, `DISTANCE`, `STRAND`, `FLAGS`, `VARIANT_CLASS`, `SOURCE`, `VARIANT` and `VARIANT_AF`:

```
MODIFIER	-	-	-	SNV,exact	-	exact_2	0.004
MODIFIER	-	-	-	SNV	-	-	-
```

In the default format the first row's `Extra` column reads `IMPACT=MODIFIER;VARIANT_CLASS=SNV,exact;VARIANT=exact_2;VARIANT_AF=0.004;VARIANT_FILTER=PASS`, the last item being the source record's `FILTER`, which that format prints for a VCF source. In the VCF the first allele's `CSQ` entry holds `SNV&exact`.

The flag and the format decide the outcome. Without `--variant_class` both engines write `exact` on the first row and `-` on the second. With `--json` both write `"variant_class":"SNV"` and the source's record under `custom_annotations`. Under `type=overlap`, where the source record names no allele, Ensembl VEP dies the same way and vep-rs writes `SNV,exact` on both rows.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP builds one set of fields for the input record and copies it for each row. Under `--variant_class` that set holds `VARIANT_CLASS` as the string `SNV`. The row of the allele `C` then takes the source records that name that allele, and the push of their `CLASS` value goes to the key `VARIANT_CLASS`.

A push expects a list under the key. It makes one where the key is unset, which is how the same key takes `exact` with the flag off. Here the key holds a string, Perl reads the string as the name of a list, and `strict` ends the process.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/OutputFactory.pm`):

- Line 967: `VARIANT_CLASS` is set on the input record under `--variant_class`.
- Line 391: the row's fields start as a copy of the input record's.
- Lines 1426-1435: the source records that name the row's allele are handed to the push.
- Line 1303: the push that fails.

### 1.5 Why vep-rs's row is the correct one

The input record's alleles are single bases, so its class is `SNV`: the value both engines write for it under any other source name. The source's record at position 5380280 names the allele `C`, so its `CLASS`, `exact`, belongs to the row of `C` and to no other: the value both engines write there with the flag off.

The command asks for both under one key, and Ensembl VEP sets the class before it pushes the source's value. The row of `C` therefore holds two values under `VARIANT_CLASS`, the class and then the source's value, and the row of `G` holds the class alone. An exit that writes no row answers neither.

### 1.6 Try it yourself

Save the record of 1.2 as `record.vcf` and run vep-rs:

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --variant_class --tab \
  --custom file=<corpus>/custom/sites.vcf.gz,short_name=VARIANT,format=vcf,type=exact,fields=AF%CLASS
```

vep-rs writes the rows of 1.3. Ensembl VEP 116.3 takes the same flags with the cache options of [Reproducing these examples](reproducing.md), and ends with the line of 1.3.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/custom_key_named_as_a_plain_value_field.toml` runs the record of 1.2 on both engines and records the `VARIANT_CLASS` cell of each row: no row and exit status 255 for Ensembl VEP, `SNV,exact` and `SNV` for vep-rs.

No comparator sets the class aside: no committed flag set names a source whose key spells one of the fourteen fields.

## 2. The fourteen fields

Each field is a plain value when the push the third column names is made. Ensembl VEP 116.3 dies at that push when a source's key spells the field, and vep-rs writes the field's value followed by the values of that push.

| Field | Set at line | Push that fails | Flag |
| --- | --- | --- | --- |
| `VARIANT_CLASS` | 967 | the first or the second | `--variant_class` |
| `IND` | 971 | the first or the second | `--individual` |
| `ZYG` | 976 | the first or the second | `--individual` |
| `MINIMISED` | 992 | the first or the second | `--minimal` |
| `AMBIGUITY` | 999 | the first or the second | `--ambiguity` |
| `CHECK_REF` | 1020 | the second | `--check_ref --dont_skip` |
| `IMPACT` | 1344 | the second | none |
| `ALLELE_NUM` | 1354 | the second | `--allele_number` |
| `REF_ALLELE` | 1357 | the second | `--show_ref_allele` |
| `UPLOADED_ALLELE` | 1360 | the second | `--uploaded_allele` |
| `PICK` | 1363 | the second | `--flag_pick` and its variants |
| `HGVSg` | 1393 | the second | `--hgvsg` |
| `SPDI` | 1401 | the second | `--spdi` |
| `GA4GH_VRS` | 1410 | the second | `--ga4gh_vrs` |

The first five are set before both pushes, so a source record of either kind reaches a failing push. vep-rs writes the field's value, then the names or values of the records that name no allele, then those of the records that name the row's allele. Under `--individual_zyg` the field `ZYG` is a list, which a push appends to, and `IND` is not set: that list is not one of the fourteen, and vep-rs writes it followed by the source's values.

The other nine are set after the first push and replace its values. Under a source whose records name no allele Ensembl VEP therefore completes, and both engines write the field's value alone on a row that sets the field. On the input record of section 1 a source named `IMPACT` of `type=overlap` leaves `MODIFIER` on both rows. One named `PICK` under `--flag_pick`, which marks the row of `C`, leaves `1` on that row and `exact_2` on the row of `G`, where the field is not set.

Under the same two sources of `type=exact` Ensembl VEP dies at the second push. vep-rs writes `MODIFIER,exact_2` and `1,exact_2` on the row of `C`, and `MODIFIER` and `-` on the row of `G`.

A `--summary` or `--most_severe` row is the input record's fields alone: it takes the first push and no second, and of the fourteen it carries the first five and `CHECK_REF`. Under a source whose records each name an allele both engines write the field's value alone there. Where one names none, Ensembl VEP dies at the first push when the field is one of the first five, and vep-rs writes the field's value and the source's after it.

Three of the thirteen fixed columns are plain values at a push as well: `Uploaded_variation` and `Location`, set at lines 894-932, and `Allele`, set at line 1347. Ensembl VEP's run ends the same way under a key that spells one of them, the first two as the first five fields and `Allele` as the other nine. They are not among the fourteen: [section 4](#4-the-three-columns-that-name-a-row) gives each engine's cells under them.

## 3. A key that spells a source's statistic

**In short.** `summary_stats` asks a source for statistics over the records it matches at an input record: `count`, the number of those records, and `min`, `max`, `mean` and `sum` of their scores ([Custom annotation](../custom-annotation.md#records-and-statistics-per-input-record)). Each statistic has a key of its own: the source's short name, `_` and the statistic's name. Ensembl VEP 116.3 assigns the statistic on the input record as one plain number. A source record pushed under the same key afterwards ends the run as in section 1. vep-rs writes the statistic, then the values of the source records that name the row's allele.

Two kinds of record are pushed under a statistic's key. A second source whose short name is the key lists its record names there. A record field whose name is a statistic's name (`count`, `min`, `max`, `mean` or `sum`) is listed under its source's short name, `_` and that name, which is the key of that source's statistic of the same name.

A key holds a statistic only where Ensembl VEP assigns one. `count` is assigned for every source that asks for it and matched a record. `min`, `max`, `mean` and `sum` are assigned where a matched record carries a score, the fifth column of a BED record or the `QUAL` of a VCF record ([Custom annotation](../custom-annotation.md#values)). The records of the source file of section 1, `custom/sites.vcf.gz`, carry none: their `QUAL` is `.`. So where a source `S` over that file asks for `summary_stats=min`, no `min` is assigned. A second source named `S_min` in the same command lists its record's name alone in both engines, and Ensembl VEP completes.

Take a BED source `R` whose one record at the input record's position has the score 4, and ask it for all five statistics. Ensembl VEP assigns each: `4` for `min`, `max`, `mean` and `sum`, and `1` for `count`. A second source of `type=exact` over the source file of section 1, named `R_min`, `R_max`, `R_mean`, `R_count` or `R_sum`, ends Ensembl VEP's run at line 1301. Under the key that second source is named as, vep-rs writes the statistic and then `exact_2` on the row of `C`, `4,exact_2` or for `R_count` `1,exact_2`, and the statistic alone on the row of `G`.

### 3.1 The record and the command

The input record and the source file, the corpus's `custom/sites.vcf.gz`, are those of section 1. The file is given twice, as a source named `S` that asks for the count of its records and as a source named `S_count`. Neither asks for `fields`, and a VCF source's `FILTER` item prints only with `fields`, so the rows carry no `FILTER` item:

```
--custom file=<corpus>/custom/sites.vcf.gz,short_name=S,format=vcf,\
type=exact,summary_stats=count
--custom file=<corpus>/custom/sites.vcf.gz,short_name=S_count,format=vcf,\
type=exact
--tab
```

### 3.2 What each engine writes

Ensembl VEP 116.3 writes its header and no row. It ends on `Can't use string ("1") as an ARRAY ref while "strict refs" in use`, at line 1301 of `OutputFactory.pm`, with the exit statuses of section 1. vep-rs writes the two rows. Their last two columns are `S` and `S_count`:

```
exact_2	1,exact_2
-	1
```

In the default format the first row's `Extra` column reads `IMPACT=MODIFIER;S=exact_2;S_count=1,exact_2` and the second's `IMPACT=MODIFIER;S_count=1`. In the VCF the first allele's `CSQ` entry holds `1&exact_2` and the second's `1`. Ensembl VEP completes in JSON, where both engines write the statistic, the object `{"count":1}`, as the first item of the list of `S` under `custom_annotations`. It also completes on a `--most_severe` row, where both write `1` under `S_count`.

One source is enough where a field of its records is named as a statistic. Give `S`, as `type=exact`, a file whose record on the allele `C` carries the INFO field `count=7`, and ask for `fields=count,summary_stats=count`. Ensembl VEP 116.3 dies at line 1303, the push of a field's value, and vep-rs writes `1,7` under `S_count` on the row of `C` and `1` on the row of `G`.

### 3.3 Why Ensembl VEP writes that

At the first push Ensembl VEP calls `_add_custom_annotations_to_hash` once for each short name that matched a record, with that short name's records that name no allele and with its statistics. The sub pushes the records, if any, and then assigns each statistic under the short name, `_` and the statistic's name: the input record's fields hold `S_count` as the number `1`. At the second push it calls the same sub with the records that name the row's allele, and with no statistics. On the row of `C` the sub pushes `exact_2` under `S_count`, onto the number.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/OutputFactory.pm`):

- Line 1012: the first push hands the statistics of the short name to the sub.
- Lines 1307-1309: the sub assigns each statistic as a plain value.
- Lines 1426-1435: the second push, which hands no statistics.
- Lines 1301 and 1303: the pushes that fail.

### 3.4 Why vep-rs's row is the correct one

Source `S` matches one record of the file at the position, so its count is `1`. A statistic belongs to the whole input record, so both rows carry it. The record `exact_2` names the allele `C`, so its name belongs under `S_count` on the row of `C` and on no other: both engines write it there when `S` asks for no statistic.

The command asks for both under one key, and Ensembl VEP assigns the statistic before it pushes the record's name. The row of `C` therefore holds the statistic and then the name, and the row of `G` the statistic alone.

### 3.5 A record that names no allele

A record that names no allele is pushed at the first push. Under a key that holds a statistic, the statistic has replaced it wherever Ensembl VEP completes. A record of the statistic's own source is pushed before the statistics are assigned, inside one call of the sub. Give the source with the `count` field as `type=overlap`: Ensembl VEP's run completes, both engines write `1` on both rows, and the field's value `7` prints on neither.

A record of another source is pushed before or after the assignment. Ensembl VEP walks the short names as the keys of a hash (line 1004), in Perl's hash order, which changes between runs unless the interpreter's `PERL_HASH_SEED` and `PERL_PERTURB_KEYS` are fixed. With the two sources of 3.1 given as `type=overlap` and with `PERL_PERTURB_KEYS=0`, Ensembl VEP 116.3 completes with `1` on both rows under `PERL_HASH_SEED=11` and dies at line 1301 under `PERL_HASH_SEED=0`.

Of those two orders vep-rs takes the one in which Ensembl VEP completes. In the other order the names would follow the statistic. It writes `1` on both rows, the cell of every run Ensembl VEP completes, whichever of the two sources the command names first. The record names of the source whose short name is the key then print in no column of either engine's row; the JSON record holds them.

### 3.6 The description of a key one source declares twice

The source of 3.5 declares `S_count` twice, as the record field `count` and as the statistic `count`. A header describes a field's key as `<Description>. <field> field from <file>` and a statistic's as `<stat> data from <file>` ([Custom annotation](../custom-annotation.md#output)). Ensembl VEP tests a source's keys against the keys of the earlier sources alone, so it lists this key twice, each time with its own text. Its tab header reads a key's text from a hash, in which the statistic's text replaces the field's, and its VCF header prints an `##INFO` line for each listing.

Name the file `statfield.vcf.gz` and let its header describe the field as `A made-up field named as a statistic`. With `type=overlap`, where Ensembl VEP completes, the headers read as follows; `F` stands for `[PATH]/statfield.vcf.gz`:

| Header | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Tab | `## S_count : count data from F`, twice | `## S_count : count data from F; A made-up field named as a statistic. count field from F`, once |
| VCF | two `##INFO=<ID=S_count,...>` lines, with the field's text and then with the statistic's | one line, with the statistic's text, `; ` and the field's |
| Default | a `## S_count :` line with the field's text, then one with the statistic's | the same two lines |
| Columns | two `S_count` tab columns and two `CSQ` fields, each with the row's one value | one column and one field |

vep-rs lists the key once, as one column and one `CSQ` field ([A key listed twice by the tab and VCF writers](duplicate-output-column.md)), and its one line in the tab and VCF headers carries both texts, the statistic's first.

What the cell can hold decides the text. A header is written before any row, and the cell under such a key holds one of three things: the statistic alone, the field's values alone, or both. It holds the field's values alone where Ensembl VEP assigns no statistic. Where it assigns one, the cell holds the statistic, followed by the field's value on a row of an allele that an `exact` record names, the `1,7` of 3.2. No one text is true of all three, and the line follows the order of the values in the cell.

A statistic that needs a score shows the field alone. Give the file's record the field `min=7` in place of `count`, with `.` in its `QUAL` column, and ask for `fields=min,summary_stats=min`. No `min` is assigned, and both engines print the field's `7` in the cell. They do so under `type=exact` as well, on the row of `C`, where Ensembl VEP then completes. Ensembl VEP's tab header still reads `## S_min : min data from ...` on both of its lines. With a `QUAL` of 30 and `type=overlap` both engines print the statistic, `30`.

With `count` under `type=overlap` the statistic always replaces the field's value, which the JSON record alone then holds. The field's half of the line names that value all the same: the line is one text for every type and every file.

From version 4.3 on, the VCF specification requires the `ID` of an `##INFO` line to be unique among them; the text of version 4.2 does not state the rule. Both engines write the `##fileformat` line of the input record's file, `VCFv4.2` here. Given the same record in a file that declares `VCFv4.3`, Ensembl VEP's output declares that version over the same two lines. A reader that keeps the first of two `##INFO` lines of one `ID`, as `bcftools` 1.21 does without a warning, reads the field's text alone from Ensembl VEP's header, over a cell that holds the statistic.

The default format's header lists a source's declarations and no columns, so vep-rs keeps Ensembl VEP's two lines there. Under a `--fields` list that names the key, Ensembl VEP's tab header has one line for it, with the statistic's text, and vep-rs's has its joined line.

A field or a statistic named twice (`fields=count%count`) declares its key twice with one text: Ensembl VEP repeats the line, and vep-rs prints it once in the tab and VCF headers. Two sources that share a short name are another case. There the later source's file joins, after a comma, the text of the declaration that source repeats ([the repeated short name's header lines](pc-description-text.md#2-the-header-line-that-receives-a-repeated-custom-short-names-file)).

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `OutputFactory.pm`, lines 2444-2512, `get_custom_headers`: line 2451 reads the keys listed so far, once a source, before the source's own are pushed; lines 2482 and 2499 test a field's and a statistic's key against that reading; lines 2493 and 2503-2506 push the key with its text.
- `OutputFactory/Tab.pm`, line 149: the hash of the texts, in which the later text of a key stands; line 153: one description line per listed column.
- `OutputFactory/VCF.pm`, line 276: one `##INFO` line for each listed key. `OutputFactory/VEP_output.pm`, line 167: one key line for each.
- `AnnotationSource/File.pm`, line 346: `count` is set for every source that asks for it; lines 424-434: `min`, `max` and the sum behind `mean` and `sum` are set only where a matched record carries a score.

## 4. The three columns that name a row

Each table row is a command on the record and source of section 1, the source named as the column, under `--tab` unless the row says otherwise.

| Key and source | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Location`, `type=exact` | ends at line 1301 | `21:5380280,exact_2` on the row of `C`, `21:5380280` on the row of `G` |
| `Location`, `type=overlap` | ends at line 1301, under `--most_severe` or `--summary` too | `21:5380280,exact_2` on every row, the `--most_severe` or `--summary` row included |
| `Location`, `type=exact`, `--most_severe` or `--summary` | `21:5380280`: line 387 returns that row before the rows of the alleles are built, so it takes no second push | `21:5380280` |
| `Location`, either type, under `--shift_3prime 1 --shift_genomic 1`, on the row of any sequence variant against a transcript: lines 1711-1715 assign the row's location again, after the call that makes the second push | ends at line 1301 | the location that assignment gives, alone, since it replaces the list |
| `Location`, either type, under the same two flags, on the intergenic row of an inserted or deleted alternate allele: lines 2001-2007 assign the row's location again, before the call that makes the second push | ends at line 1301 | the location that assignment gives, then the names of the records that name the row's allele, and none of a record that names no allele |
| `Uploaded_variation`, either type | ends or completes as in the first three table rows | as in those three, with `rs1444031311` in place of `21:5380280`; no row assigns the identifier again |
| `Uploaded_variation` as the key of a source named `Uploaded` whose records carry a field `variation` that it asks for, either type | ends at line 1303 where a matched record carries the field; completes where none does, and under `type=exact` on a `--most_severe` or `--summary` row | `rs1444031311`, then the field's values, on the output rows that hold `exact_2` in the first three table rows |
| `Allele`, `type=exact` | ends at line 1301 | `C,exact_2` on the row of `C`, `G` on the row of `G` |
| `Allele`, `type=exact`, `--most_severe` or `--summary` | `-`, since that row sets no allele and takes no second push | `-` |
| `Allele`, `type=exact,num_records=0`, a limit under which the source lists no record but one named `...`, three dots, that names no allele | the row's allele; the name `...` on a `--most_severe` or `--summary` row | the same |
| `Allele`, `type=overlap` | the row's allele; `exact_2` on a `--most_severe` or `--summary` row | the same |
| any of the three, `--json` | the input record's own `id`, coordinates and alleles, and the source's record under `custom_annotations` | the same |

`Uploaded_variation` and `Location` are set at lines 895 and 931, before both pushes, and behave as the first five fields of section 2, `Location` on a row that does not assign it again. `Allele` is set at line 1347, between the pushes, and behaves as the other nine. Every other cell is as under another source name.
