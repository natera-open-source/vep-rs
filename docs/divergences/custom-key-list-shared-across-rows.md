# A `--custom` key's list that every row of a record shares

**The defect in one sentence.** A `--custom` key that spells a column or a flag field can hold a list in the hash of the input record. Ensembl VEP 116.3 builds every row of the record over that one list, so what one row adds to the key prints on other rows too. vep-rs lists on each row the source's values and the row's own items, each once.

**Why this is a defect.** A cell states what holds for its row: the frequency of the row's allele, the known variants at the record's site, the sets the row's transcript belongs to. A source whose key spells such a cell's name shares the cell with it, and both engines then list the source's values beside the row's own. Ensembl VEP's cell does not follow the row. It prints one allele's frequency on the rows of another ([1.3](#13-what-each-engine-writes)). It lists a known variant's identifier once for each row of the record ([2.2](#22-what-each-engine-writes)). It prints `MANE_Select` on the rows of transcripts outside that set ([3.2](#32-what-each-engine-writes)).

The same output contradicts the cell. Each population column beside a shared frequency key holds the frequency of the row's own allele alone, and the column `MANE_SELECT`, which holds the set's RefSeq match, is filled on the row of the set's transcript alone. Where the record's hash holds no list under the key, Ensembl VEP writes vep-rs's cells: sections 1 to 3 give such controls for every key but `NEAREST`.

**What vep-rs does instead.** vep-rs builds each row's cell from the row: the values the source's key carries on it and what the row's own field holds, each item once. Under a frequency flag the field holds each frequency of the row's own allele that none of the source's values spells. Under `--check_existing` it holds the co-located identifiers, and under `--var_synonyms` the synonyms. Under `--mane_select` it holds `MANE_Select` on the row of a transcript of that set. Under `--nearest` the nearest names stand first and the records of the row's allele after them.

The values are joined as the output format joins a list, by a comma in the default and tab outputs and by `&` in a VCF; `VAR_SYNONYMS` joins its items by two hyphens.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3; `OutputFactory.pm` is `Bio/EnsEMBL/VEP/OutputFactory.pm`):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/VCF.pm`, lines 227-263: a record is built with the `allele` it matched only where the match returned the list of matched alleles (line 228), which the `exact` match of lines 405-421 returns. The `overlap` match of lines 389-404 returns a truth value and a percentage, so an `overlap` source's record names no allele.
- `OutputFactory.pm`, lines 1003-1014, `VariationFeature_to_output_hash`: the records without an `allele` go into the hash of the input record, which is built once for the record (line 347). They are handed to `_add_custom_annotations_to_hash` (lines 1297-1312), whose line 1303 pushes each field's value onto a list under the short name, an underscore and the field.
- `OutputFactory.pm`, lines 389-404: each row's hash starts as `my %copy = %$hash;` (line 391). The copy is one level deep, so the key of every row names the one list of the record's hash.
- `OutputFactory.pm`, lines 1441-1443, call `add_colocated_frequency_data` for each row. Line 1200 finds the known variant's allele that matches the row's and lines 1234-1242 read its frequency. Line 1247 adds it to the list under the column's key, `merge_arrays($hash->{$out_key} ||= [], [$f])`, the one statement every frequency column is written through.
- `Bio/EnsEMBL/VEP/Utils.pm`, lines 303-310, `merge_arrays`: each value the list lacks is pushed onto the list itself.
- `OutputFactory.pm`, lines 266-272 and 288-304, turn the rows into lines once every row has its hash; the VCF writer does the same for each record (`Bio/EnsEMBL/VEP/OutputFactory/VCF.pm`, lines 342-345). So each row prints the list as the record's last row left it.
- `Bio/EnsEMBL/VEP/OutputFactory/JSON.pm`, lines 250-251, 278-286 and 513-515: the JSON factory moves the source's records to `custom_annotations` and reads each allele's frequencies into a hash of its own, so no list is shared there.
- `OutputFactory.pm`, line 1301, in the same `_add_custom_annotations_to_hash`: each record's name is pushed onto a list under the short name alone.
- `OutputFactory.pm`, lines 1437-1438, call `add_colocated_variant_info` for each row. Its line 1086 pushes each co-located variant's name onto the list under `Existing_variation`, and line 1090 its synonyms onto the list under `VAR_SYNONYMS`. Line 1141 then joins the `VAR_SYNONYMS` list into a string, which replaces the key of the row's own hash.
- `OutputFactory.pm`, lines 1596-1602: on the row of a transcript that carries the attribute, line 1601 pushes `MANE_Select` onto the list under `MANE`.
- `OutputFactory.pm`, line 1017, sets `NEAREST` in the record's hash to the list the nearest search left on the input record (`Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, line 172). Lines 1426-1435 push the records that name the row's allele, through line 1301.

**What correcting it would change.** A copy at line 391 that gave each row its own lists would have Ensembl VEP write vep-rs's rows. The same copy is behind [`custom_exact_source_allele_record_accumulation`](../custom-annotation.md#1-an-exact-custom-record-repeated-once-per-row-beside-a-whole-variant-record-of-the-same-name), where another statement adds to a shared list: a source's own record, pushed once for each row.

**The class it produces.** `custom_key_list_shared_across_rows` is the cell of a `--custom` key on a default or tab row and, for the frequency key of section 1, a VCF row. The key spells a column or a flag field, and the hash of the input record holds a list under it: the values of a source's record of the whole input record, or the nearest names. Ensembl VEP's cell holds the source's values and what the record's rows added to that list. A key printed from the list holds every row's items; `VAR_SYNONYMS`, which each row joins as it is built, holds those of the rows built so far. vep-rs's cell holds the source's values and the row's own items.

Sections 1 to 3 give the five keys that were run in both engines. The class is documented only: no committed flag set names such a source beside the flag that fills the field, so no manifest names the class and no comparator applies it.

**How often the shape occurs.** Tier 3 of [the index's order](../intended-divergences.md#differences-outside-the-classes): ordinary data does not carry the shape. It needs a source whose key spells a column or a flag field, and none of the committed flag sets passes such a source beside the flag that fills the field.

## 1. An `overlap` source named `gnomADg` that is asked for its `AF` field

**In short.** Three records of two alternate alleles each are run under `--af_gnomadg` with a VCF file as an `overlap` source named `gnomADg` that asks for `AF`. At the first site the cache holds a gnomAD genomes frequency for the allele `C` alone. Ensembl VEP 116.3 writes `gnomADg_AF=0.004,7.994e-06` on the row of `C` and on the row of `G`, where vep-rs writes `0.004` on the row of `G`. At the second site each allele holds a frequency. Ensembl VEP writes both on all four rows, and vep-rs writes on each row its own allele's.

### 1.1 What a reader needs to know first

A `--custom` source prints each field it is asked for under a key made of its short name, an underscore and the field's name. It prints the names of its records under the short name alone, and a VCF source's `FILTER` value under the short name and `_FILTER`. A source of `type=overlap` gives its record to every row of an input record the record overlaps. A source of `type=exact` gives it to the rows of the allele the record names.

`--af_gnomadg` prints eleven columns from the co-located known variants of the cache: `gnomADg_AF`, the frequency of the row's allele over the gnomAD genomes, and one column for each of ten populations. A source named `gnomADg` that asks for `AF` therefore shares the key `gnomADg_AF` with the flag. Both engines then list the source's values first and the frequency after them.

vep-rs lists the source's values and then each frequency of the row's own allele that none of them spells. A row of an allele without a frequency holds the source's values alone. The difference needs a record whose rows have different frequencies to add, one of them possibly none. Without the flag the engines agree: the key holds the source's values alone. They agree as well where the source's record names an allele, as an `exact` source's does, since such a record is listed on the rows of its allele alone. The JSON output is the same in both engines: the frequencies stand under `colocated_variants`, by allele, and the source's record under `custom_annotations`.

### 1.2 The records

The first and the third record are lines 9 and 10 of `tests/golden/116/GRCh38-custom/variants.vcf`. The second is the first with its alternate alleles in the other order:

```
21	5380280	rs1444031311	A	C,G	.	.	.
21	5380280	alleles_reversed	A	G,C	.	.	.
21	8987421	rs766964238	G	C,T	.	.	.
```

The source is `tests/golden/116/GRCh38-custom/custom/sites.vcf.gz`, with its `.tbi` index beside it. It holds one record at each site, and the record's name and `CLASS` value are the corpus's own labels:

```
21	5380280	exact_2	A	C	.	PASS	AF=0.004;CLASS=exact
21	8987421	overlap_2	G	A	.	PASS	AF=0.005;CLASS=overlap
```

Read as `overlap`, each record is matched to an input record by its span alone. vep-rs reads the known variants from the cache of `tests/golden/116/GRCh38-everything` and Ensembl VEP from the release's GRCh38 cache, which holds the same line for each of the two. The gnomAD genomes column of `rs1444031311` reads `C:7.994e-06`: the cache holds no frequency for `G`. That of `rs766964238` reads `C:3.593e-05,T:1.794e-05`. Both engines run:

```
--af_gnomadg --custom file=sites.vcf.gz,short_name=gnomADg,format=vcf,\
type=overlap,fields=AF%CLASS
```

### 1.3 What each engine writes

The cell `gnomADg_AF` of each row:

| Record | Allele | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `rs1444031311` | `C` | `0.004,7.994e-06` | `0.004,7.994e-06` |
| `rs1444031311` | `G` | `0.004,7.994e-06` | `0.004` |
| `alleles_reversed` | `G` | `0.004,7.994e-06` | `0.004` |
| `alleles_reversed` | `C` | `0.004,7.994e-06` | `0.004,7.994e-06` |
| `rs766964238` | `C` | `0.005,3.593e-05,1.794e-05` | `0.005,3.593e-05` |
| `rs766964238` | `T` | `0.005,3.593e-05,1.794e-05` | `0.005,1.794e-05` |

The first two records have one row for each allele. Each allele of `rs766964238` has two rows, one for each of two transcripts, and both hold the cell of the table. The default output holds the cell as an item of `Extra`, the tab output as a column and the VCF output as a `CSQ` field whose values `&` joins.

On the row of `G` of `rs1444031311`, the `Extra` column of Ensembl VEP's default output reads `IMPACT=MODIFIER;gnomADg_AF=0.004,7.994e-06;gnomADg=exact_2;gnomADg_CLASS=exact;gnomADg_FILTER=PASS`. vep-rs's reads `IMPACT=MODIFIER;gnomADg_AF=0.004;gnomADg=exact_2;gnomADg_CLASS=exact;gnomADg_FILTER=PASS`. The `CSQ` field of that row is `0.004&7.994e-06` in Ensembl VEP's VCF and `0.004` in vep-rs's.

Apart from the key's second listing, described below, every other value of the eight rows agrees between the engines. On the rows of `rs766964238` that includes the ten population columns: `gnomADg_AFR_AF` is `2.801e-05` on the rows of `C` and `0` on the rows of `T`, and `gnomADg_SAS_AF` is `0` and `0.001082`. The JSON records of the three input records are equal between the engines, apart from the `vep_rs_version` key that vep-rs writes on every record.

Ensembl VEP lists the key twice in the tab and VCF outputs: two columns, and two `CSQ` fields, that hold the same cell. That is [`duplicate_output_column_collapsed`](duplicate-output-column.md), and vep-rs lists the key once.

Two controls on the first record show that the difference needs both the flag and a record of the whole input record. Without `--af_gnomadg`, both engines write `0.004` on both rows. With the same file read as `type=exact`, both write `0.004,7.994e-06` on the row of `C` and no value on the row of `G`.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP builds one Perl hash for the input record and one for each row. The record's hash takes what belongs to every row, the `overlap` source's record included: `gnomADg_AF` becomes a list holding `0.004`. Each row's hash is a copy of the record's. Perl holds a list in a hash as a reference, and the copy copies the reference, so `gnomADg_AF` is one list under every row.

The frequency routine then runs for each row and adds the frequency of the row's allele to the list its key holds, unless the list holds that text already. On `rs1444031311` the row of `C` adds `7.994e-06` and the row of `G` adds nothing. On `rs766964238` the first row of `C` adds `3.593e-05` and the first row of `T` adds `1.794e-05`.

The lines are written after every row of a record is built, so each row prints the whole list. The second record shows it: its row of `G` is built before the row of `C` adds the frequency, and prints it too.

A column that no key of a source spells has no list in the record's hash. The routine then starts one in the row's own hash, which is why the population columns hold each row's own allele.

### 1.5 Why vep-rs's rows are the correct ones

Each value follows from the command and the two files. Read as `overlap`, the source's record is matched to the input record's span and to no allele of it, so every row owes its `AF`: `0.004` at the first site and `0.005` at the second. The cache line of `rs1444031311` gives `C` the frequency `7.994e-06` and gives `G` none, so the row of `G` owes no frequency. The cache line of `rs766964238` gives `C` and `T` one frequency each, so a row owes its own allele's. Those are vep-rs's cells.

Ensembl VEP writes the same values wherever its rows share no list. Read as `exact`, the source leaves the row of `G` without a frequency. The population columns never hold another allele's value, and the JSON output keys each frequency to its allele. The order inside the cell is Ensembl VEP's own for a key that a source and the frequency routine share: the source's values first.

### 1.6 Try it yourself

With vep-rs, where `records.vcf` holds the records of 1.2 under a VCF header, `<cache>` is the `json_cache` directory of `tests/golden/116/GRCh38-everything` and `sites.vcf.gz` is the source of 1.2:

```
vep -i records.vcf -o out.txt --offline --json_cache <cache> \
  --assembly GRCh38 --af_gnomadg --custom file=sites.vcf.gz,\
short_name=gnomADg,format=vcf,type=overlap,fields=AF%CLASS
```

vep-rs writes the cells of 1.3. Ensembl VEP 116.3 takes the same `--af_gnomadg` and `--custom` options with `--offline --cache` and the release's GRCh38 cache, and writes every allele's frequency on every row. Naming the source so that no key of it is a frequency column makes the engines agree.

## 2. An `overlap` source named `Existing_variation` under `--check_existing`

**In short.** The three records of section 1 are run under `--check_existing` with the same file as an `overlap` source named `Existing_variation`. The source's record names share that column with the identifiers of the co-located known variants. Ensembl VEP 116.3 lists an identifier once for each row of its record: twice on the two rows of the first record, four times on the four rows of the third. vep-rs lists it once.

### 2.1 The command

`--check_existing` fills the column `Existing_variation` with the identifiers of the known variants the cache holds at a record's site. A source's record names print under its short name, so a source named `Existing_variation` lists them in that column, before the identifiers. The cache holds one known variant at each of the two sites, `rs1444031311` and `rs766964238`. Both engines run the records and the source of [1.2](#12-the-records):

```
--check_existing --custom file=sites.vcf.gz,short_name=Existing_variation,\
format=vcf,type=overlap,fields=AF%CLASS
```

### 2.2 What each engine writes

The cell `Existing_variation`, which every row of a record holds:

| Record | Rows | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `rs1444031311` | 2 | `exact_2,rs1444031311,rs1444031311` | `exact_2,rs1444031311` |
| `alleles_reversed` | 2 | `exact_2,rs1444031311,rs1444031311` | `exact_2,rs1444031311` |
| `rs766964238` | 4 | `overlap_2,rs766964238,rs766964238,rs766964238,rs766964238` | `overlap_2,rs766964238` |

The default output holds the cell in its `Existing_variation` column, and the tab output in the column of that name. Ensembl VEP's tab output lists that column twice, as in [1.3](#13-what-each-engine-writes). Every other cell of the eight rows agrees between the engines.

Three controls on the first record read the same in both engines. Without `--check_existing` the cell is `exact_2` on both rows. Under `--summary`, with `--check_existing` kept, the one row's cell is `exact_2`: a summary row lists no identifier in either engine. With the same file read as `type=exact` the cell is `exact_2,rs1444031311` on the row of `C`, the allele the source's record names, and `rs1444031311` on the row of `G`.

### 2.3 Why Ensembl VEP writes that, and why vep-rs's cell is the correct one

The `overlap` source's record goes into the record's hash, where `Existing_variation` becomes a list that holds its name. Every row's hash holds that one list ([1.4](#14-why-ensembl-vep-writes-that)). The routine that reads the co-located variants runs for each row and pushes each variant's identifier onto the list its key holds. A record of four rows pushes it four times, and each row prints the whole list.

The inputs hold one source record and one known variant at each site, so a row owes one name and one identifier. Ensembl VEP's count of the identifier is the number of rows the record has, which neither the cache nor the source states. Read as `exact`, the source leaves no list in the record's hash, each row starts its own, and Ensembl VEP lists the identifier once.

## 3. A source named `MANE`, `VAR_SYNONYMS` or `NEAREST`

**In short.** Line 42 of `tests/golden/116/GRCh38-custom/variants.vcf` is one substitution with sixteen rows, one for each of sixteen transcripts. Three runs give the corpus's source file a short name that spells a flag field that is a list. Ensembl VEP 116.3 writes `MANE_Select` on all sixteen rows, repeats the variant's synonyms once more on each later row, and lists the allele's record sixteen times on every row. vep-rs writes each item once, on the row it belongs to.

### 3.1 The record and the commands

The input record, and the source's one record at its site:

```
21	33549625	rs2145814066	C	T	.	.	.
21	33549625	multi_7	C	T,A	.	PASS	AF=0.021,0.521;CLASS=multi_allelic
```

The first alternate allele of the source's record is the input record's. The cache holds one known variant at the site, `rs2145814066`, with one string of synonyms, `ClinVar::RCV002273474,VCV001699617`. It gives the MANE Select attribute to one of the sixteen transcripts, `ENST00000356577`, and the nearest gene's symbol is `SON`. The three runs:

```
--mane_select --custom file=sites.vcf.gz,short_name=MANE,format=vcf,\
type=overlap,fields=AF%CLASS
--check_existing --var_synonyms --custom file=sites.vcf.gz,\
short_name=VAR_SYNONYMS,format=vcf,type=overlap,fields=AF%CLASS
--nearest symbol --custom file=sites.vcf.gz,short_name=NEAREST,format=vcf,\
type=exact,fields=AF%CLASS
```

### 3.2 What each engine writes

| Key | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `MANE` | `multi_7,MANE_Select` on all sixteen rows | `multi_7,MANE_Select` on the row of `ENST00000356577`, `multi_7` on the other fifteen |
| `VAR_SYNONYMS` | `multi_7`, then the synonyms once on the first row, twice on the second and sixteen times on the sixteenth | `multi_7` and the synonyms once on every row |
| `NEAREST` | `SON`, then `multi_7` sixteen times, on all sixteen rows | `SON,multi_7` on every row |

The cells are the tab output's, and the default output holds each as an item of `Extra`. vep-rs's `VAR_SYNONYMS` cell reads `multi_7--ClinVar::RCV002273474,VCV001699617`. Ensembl VEP's reads the same on the first row, and on the second row `multi_7--ClinVar::RCV002273474,VCV001699617--ClinVar::RCV002273474,VCV001699617`.

In each run every other cell of the sixteen rows agrees between the engines. `MANE_SELECT` is one of them: it names `NM_138927.4` on the row of `ENST00000356577` alone. Ensembl VEP's tab output also lists each key's column twice, as in [1.3](#13-what-each-engine-writes). With the same file read as `type=exact`, both engines write vep-rs's cells of the table under `MANE` and under `VAR_SYNONYMS`.

### 3.3 Why Ensembl VEP writes that, and why vep-rs's cells are the correct ones

Under `MANE` and `VAR_SYNONYMS` the `overlap` source's record puts a list in the record's hash, as in sections 1 and 2. The row of `ENST00000356577` pushes `MANE_Select` onto it, and every row prints it. Each row pushes the synonyms onto the list under `VAR_SYNONYMS` and then joins the list into a string of its own. A row therefore holds what the rows built before it pushed, and its own.

Under `NEAREST` the list is the input record's own: the names the nearest search left on it, which the record's hash takes as they are. An `exact` source's record names an allele, so it is pushed for each row of that allele, sixteen times here, onto the one list. The nearest names stand first in both engines, since the record's hash holds them before a row pushes anything.

The other keys of this page need a source's record of the whole input record to open the list; `NEAREST` does not. In [`custom_exact_source_allele_record_accumulation`](../custom-annotation.md#1-an-exact-custom-record-repeated-once-per-row-beside-a-whole-variant-record-of-the-same-name) a whole-variant record of the same short name opens the list that an `exact` record is pushed onto. Here the nearest search opens it.

`MANE_Select` is an attribute of one transcript, and the cache gives it to `ENST00000356577` alone. The cache holds one string of synonyms for the one known variant, and the source one record for the row's allele. A row owes each once, and `MANE_Select` only on its own transcript's row: those are vep-rs's cells. Ensembl VEP's counts follow the number of the record's rows and the order they are built in.

## 4. How the comparator treats this class

The reproduction `scripts/divergences/reproductions/custom_key_list_shared_across_rows.toml` holds the records, the commands and both engines' values: the run of section 1 and, as further runs, those of sections 2 and 3. The class covers the cell of a `--custom` key that spells a column or a flag field, on a row of a record whose hash holds a list under that key. A flag set is one committed run of a test corpus and a manifest lists that corpus's documented rows; [How to read the counts](../intended-divergences.md#how-to-read-the-counts) defines both and the comparators. No committed flag set names such a source beside the flag that fills the field, so no manifest names the class and the class is documented only.

The runs of sections 2 and 3 were made in the default and tab formats. The class therefore names a VCF row for the frequency key alone, and the JSON output of those runs was not compared.
