# A `--custom` key named `ENSP` under `--protein_version`

**The defect in one sentence.** Under `--protein_version` without `--protein`, Ensembl VEP 116.3 appends a translation's version to the value under the key `ENSP`, which is the list a `--custom` source of that name filed there. Its cell then reads `ARRAY(0xaaab21c79598).2`, the list's memory address and the version, in place of the source's values. vep-rs prints the source's values on every row.

Three terms recur. A source's key is its short name, under which the names of its records print. The field `ENSP` is the protein identifier of a row's transcript, which `--protein` adds. `--protein_version` adds the version of the transcript's translation to that identifier, `ENSP00000300278.2` for `ENSP00000300278`.

**Why this is a defect.** The cell is neither the source's value nor a protein identifier. It is the text Perl gives a reference to a list: an address in the memory of the process, different in every run. The value the source gave the row is lost, and nothing warns. A source's short name is any string without commas, and Ensembl VEP accepts `ENSP`.

Ensembl VEP's own rows show the value the cell should hold. With neither flag it prints the source's record name, `multi_7`, on every row of the record of section 1. Under `--protein_version` it prints the same name on the three rows whose transcript has no translation.

**What vep-rs does instead.** The key holds the source's values, as it does in a run with neither flag. With `--protein` off the row has no protein identifier, so `--protein_version` has nothing to add a version to and changes no cell. With both flags on, a row whose transcript has a protein identifier holds that identifier with its version in place of the source's values, in both engines. The key heads one tab column and one `CSQ` field ([A key listed twice by the tab and VCF writers](duplicate-output-column.md)).

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/OutputFactory.pm` unless a path says otherwise):

- `AnnotationSource/File/VCF.pm`: lines 283-297 give a record of a VCF source its identifier as its name, and lines 227-246 give a record of `type=exact` the allele it matched.
- Lines 1297-1312, `_add_custom_annotations_to_hash`: line 1301 pushes a source record's name onto the list under the source's short name.
- Lines 1004-1014: the source records that name no allele are pushed onto the input record's fields. Line 391: each row starts from a copy of those fields, and the copies share each list.
- Lines 1426-1435: on each row, the source records that name the row's allele are pushed onto the row's fields.
- Lines 1696-1697: a transcript's row makes that push first and then runs lines 1565-1574.
- Lines 1565-1568: where `--protein` is on and the transcript has a protein identifier, `ENSP` is assigned that identifier.
- Lines 1569-1574: under `--protein_version`, where the key holds a true value that does not end in a dot and digits and the transcript's translation has a version, `.=` appends a dot and that version to the key's value.
- Line 387: `--summary` and `--most_severe` return before any row of a transcript is built.
- `OutputFactory/JSON.pm`, lines 249-258: the JSON writer moves the custom annotations from the input record to the output's `custom_annotations` before any row is built. Both pushes read them from the input record, so no row holds the list.

**What correcting it would change.** An append made only where line 1565 assigned an identifier would leave the source's list as it stands. Ensembl VEP would then print the source's values on every row, as vep-rs does.

**The class it produces.** `custom_key_named_ensp_under_protein_version`: the cell under `ENSP` on each row that holds a record of the source and whose transcript's translation has a version, in the default, tab and VCF formats. The run names a `--custom` source `ENSP` and passes `--protein_version` without `--protein`. The reference, Ensembl VEP's side of a comparison, holds the address and the version there. The candidate, vep-rs's side, holds the source's values. No comparator holds a rule for the class, so a comparison reports each such cell as a difference ([1.7](#17-how-the-comparator-treats-this-class)).

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a `--custom` source whose short name is `ENSP`, in a run that passes `--protein_version` and not `--protein`. No committed flag set names a `--custom` source whose short name is `ENSP`. No record of the cuts and no golden corpus carries the shape, which is a property of the command.

## 1. A source named `ENSP` under `--protein_version`

**In short.** The input record is a substitution inside a gene, with rows on sixteen transcripts. Thirteen of them have a translation. A VCF source named `ENSP` holds, at the input record's position, one record that names the input record's allele. Its name is `multi_7`. Under `--protein_version` alone, Ensembl VEP 116.3 prints `multi_7` on the three rows without a translation and an address with a version, such as `ARRAY(0xaaab21c79598).2`, on the other thirteen. vep-rs prints `multi_7` on all sixteen.

### 1.1 What a reader needs to know first

`--custom` annotates each input record from a file of the user's own. Its `short_name` names the source and is the key under which the names of its records print, which for a VCF source are their identifiers. A VCF source of `type=exact` matches an input record's alleles one by one unless either record is structural, and a record so matched names one allele. Any other match is by position; its record names none. [Custom annotation](../custom-annotation.md) describes the option.

A key that spells a field of the output shares that field's column. [A `--custom` key that spells a field Ensembl VEP holds as a plain value](custom-key-named-as-a-plain-value-field.md) describes the fields Ensembl VEP sets as one value before a source's values arrive. Ensembl VEP assigns `ENSP` after them, and only under `--protein`.

### 1.2 The record

The input record is line 42 of `variants.vcf`, the input of the `GRCh38-custom` golden corpus, read alone against that corpus's cache. `<corpus>` below stands for the corpus's directory, `tests/golden/116/GRCh38-custom`:

```
21	33549625	rs2145814066	C	T	.	.	.
```

The source is the corpus's `custom/sites.vcf.gz`, whose record at the position is named `multi_7` and names the input record's allele. Both engines run offline with:

```
--custom file=<corpus>/custom/sites.vcf.gz,short_name=ENSP,format=vcf,type=exact,fields=AF%CLASS
--protein_version --tab
```

### 1.3 What each engine writes

Both engines write sixteen rows, and every cell but the one under `ENSP` is equal. The last eight columns are `IMPACT`, `DISTANCE`, `STRAND`, `FLAGS`, `SOURCE`, `ENSP`, `ENSP_AF` and `ENSP_CLASS`. On the row of `ENST00000300278` Ensembl VEP 116.3 writes the first line below and vep-rs the second:

```
HIGH	-	1	-	-	ARRAY(0xaaab21c79598).2	0.021	multi_allelic
HIGH	-	1	-	-	multi_7	0.021	multi_allelic
```

The cell under `ENSP` on each row, from one run of Ensembl VEP:

| Transcript | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `ENST00000290239` | `multi_7` | `multi_7` |
| `ENST00000300278` | `ARRAY(0xaaab21c79598).2` | `multi_7` |
| `ENST00000356577` | `ARRAY(0xaaab21cbf5f0).4` | `multi_7` |
| `ENST00000381679` | `ARRAY(0xaaab21cd0160).4` | `multi_7` |
| `ENST00000381692` | `ARRAY(0xaaab21cd3e50).2` | `multi_7` |
| `ENST00000436227` | `ARRAY(0xaaab21cd6b30).1` | `multi_7` |
| `ENST00000455528` | `ARRAY(0xaaab21cddbd8).1` | `multi_7` |
| `ENST00000475072` | `multi_7` | `multi_7` |
| `ENST00000492229` | `multi_7` | `multi_7` |
| `ENST00000695559` | `ARRAY(0xaaab21ce7ac0).1` | `multi_7` |
| `ENST00000704334` | `ARRAY(0xaaab21ceb520).1` | `multi_7` |
| `ENST00001036731` | `ARRAY(0xaaab21cef678).1` | `multi_7` |
| `ENST00001036732` | `ARRAY(0xaaab21cf35e8).1` | `multi_7` |
| `ENST00001129604` | `ARRAY(0xaaab21cf5630).1` | `multi_7` |
| `ENST00001129605` | `ARRAY(0xaaab21cf5a20).1` | `multi_7` |
| `ENST00001131645` | `ARRAY(0xaaab21cfdb98).1` | `multi_7` |

The digits after `0x` change with every run of Ensembl VEP, and the thirteen cells of one run hold thirteen addresses. Two runs of this command shared no address. The three transcripts whose rows read `multi_7` have no translation. The default format holds the same cells in its `Extra` column, and the VCF format in each `CSQ` entry.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP builds each row as a table of keys and values. The source's record names the row's allele, so line 1301 pushes `multi_7` onto a list under the key `ENSP` on every row (lines 1426-1435). The lines that write the protein identifier run after that push (lines 1696-1697).

Line 1565 assigns the identifier only where `--protein` is on. Without that flag the key still holds the source's list when line 1569 is reached. Of the key, its condition asks only that it hold a true value that does not end in a dot and digits, and a reference to a list passes both tests. So on a row whose transcript's translation has a version, `.=` appends a dot and the version.

Appending to a reference turns it into text first. Perl writes a reference to a list as `ARRAY(`, the list's address and `)`, and the key then holds that text with the version after it. The list itself is not changed. On a row whose transcript has no translation line 1569 does nothing, and the writer prints the list's one value.

### 1.5 Why vep-rs's row is the correct one

The source's file holds one record at position 33549625 of chromosome 21. Its name is `multi_7` and one of its alternate alleles is `T`, the input record's allele. The annotation that source gives a row of that allele, under its short name, is therefore `multi_7`, on each of the sixteen rows.

No flag of the run asks for a protein identifier. `--protein_version` adds a version to an identifier, and with `--protein` off the row has none. Nothing else of the run writes under the key, so the cell is the source's value. An address in the memory of one process is no annotation of the record.

### 1.6 Try it yourself

Save the record of 1.2 as `record.vcf` and run vep-rs:

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --protein_version --tab \
  --custom file=<corpus>/custom/sites.vcf.gz,short_name=ENSP,format=vcf,type=exact,fields=AF%CLASS
```

vep-rs writes `multi_7` under `ENSP` on the sixteen rows. Ensembl VEP 116.3 takes the same flags with the cache options of [Reproducing these examples](reproducing.md), and writes the cells of 1.3 with addresses of its own.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/custom_key_named_ensp_under_protein_version.toml` runs the record of 1.2 on both engines and records the `ENSP` cell of each row beside its transcript. It reads each address as one fixed word, since the digits change between runs.

No comparator sets the class aside: a comparison of such a run reports each of these cells as a difference. The reference's cell holds an address where the source's values stood, so a rule that set the cell aside could check the candidate's cell against nothing in the two engines' rows. No committed flag set meets the shape.

## 2. The same source under other flags, types and formats

The engines were run on the record of 1.2 under each row of this table, which gives the cell under `ENSP`. A row that names no type was run under `type=exact` and under `type=overlap`:

| Run | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `--protein_version`, `type=exact`, default, tab or VCF | an address on each of 13 rows | `multi_7` on 16 |
| `--protein_version`, `type=overlap`, default, tab or VCF | one address on 13 rows | `multi_7` on 16 |
| `--protein_version --pick`, `type=exact`, tab | an address on the row of `ENST00000356577` | `multi_7` |
| `--protein_version`, JSON | no `protein_id` | the same |
| neither flag, tab | `multi_7` on 16 rows | the same |
| `--protein`, tab | the identifier on 13 rows | the same |
| `--protein --protein_version`, tab | the identifier and version on 13 rows | the same |
| `--protein_version` with `--summary` or `--most_severe`, `type=exact`, tab | `-` on the one row | the same |
| `--protein_version` with `--summary` or `--most_severe`, `type=overlap`, tab | `multi_7` on the one row | the same |

A source of `type=overlap` files its record on the input record, and the rows share one list. Each row still turns its own copy of the reference into text, so the thirteen cells of such a run hold one address.

Under `--protein`, with or without `--protein_version`, line 1565 writes the identifier over the list on the thirteen rows, and the other three keep `multi_7`. In both of those runs Ensembl VEP's column line names `ENSP` twice, the two columns holding the same cell, and vep-rs's names it once. In JSON neither engine writes a `protein_id` without `--protein`, and the source's record stands under the input record's `custom_annotations`.

`--summary` and `--most_severe` write one row for the input record and build no row of a transcript, so line 1569 is not reached. That row holds the values filed on the input record: the name of a record that names no allele, and nothing of a record that names one, which is filed on the rows of its allele. Both engines write the row so, and it is no part of this difference.
