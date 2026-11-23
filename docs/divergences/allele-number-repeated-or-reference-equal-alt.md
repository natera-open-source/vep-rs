# `ALLELE_NUM` on a record whose ALT list repeats an allele or names the reference

**The defect in one sentence.** Under `--allele_number`, Ensembl VEP 116.3 counts a record's alternate alleles after it has dropped a repeated allele and an allele equal to the reference, so every allele written after one prints a number below its place in the record. vep-rs prints the place.

**Why this is a defect.** `ALLELE_NUM` exists to lead a reader from an output row back to an allele of the input record. The header line both engines print says so: `Allele number from input; 0 is reference, 1 is first alternate etc`. For the record `C G,G,T`, Ensembl VEP numbers the row of `T` 2. The second alternate of that record is the repeated `G`, and the row of `G` carries another consequence. A reader who follows the number reads the wrong allele.

**What vep-rs does instead.** vep-rs numbers every alternate that has a row by its place in the record's ALT column: 1 for the alternate listed first, 2 for the one listed second, 3 for the third. A repeated allele and an allele equal to the reference have no row in either engine, and each still holds its place. So does a `*`, which both engines number and give no row. An Ensembl-format line is numbered the same way, by the places of its allele string.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/`, and `ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `Parser/VCF.pm`, line 342: the feature's allele string joins every written alternate, a repeated one included.
- `VariationFeatureOverlap.pm`, lines 459-467: the alleles of the string are made unique. An allele seen before is dropped at line 463, and the reference is the first one seen.
- `VariationFeatureOverlap.pm`, lines 477-491: the alleles left are counted. The reference's own entry is passed over at line 481 and every other entry takes the next number at line 483.
- `OutputFactory.pm`, line 1354: the row prints that number as `ALLELE_NUM`.

**What correcting it would change.** With the count taken over the written allele string, or the number set before the unique filter, Ensembl VEP would print the place and the two engines would agree. Its own lookup of a row's `HGVSg`, which indexes the written string with the number, would then reach the row's own allele ([HGVSg on a record whose ALT list repeats an allele or names the reference](hgvsg-allele-lookup.md)).

**The class it produces.** `allele_number_repeated_or_reference_equal_alt`: the `ALLELE_NUM` of a row whose record lists a repeated or reference-equal alternate before the row's allele. vep-rs's number is the allele's place in the record's allele list. Ensembl VEP's is that place less the count of such entries before it. The golden tests hold a corpus's rows to the values its manifest documents. The format comparator, which compares two outputs row by row, accepts the pair under the class's name where both rows print the list.

**How often the shape occurs.** Tier 3 of [the index's order](../intended-divergences.md#differences-outside-the-classes): the shape is a record whose ALT column repeats an allele, or names the reference, before another allele. No record of the cuts, the index's five sets of chromosome 21 records, carries it. Of the committed flag sets of the golden corpora, 21 pass `--allele_number` or `--minimal`, which switches it on, and 2 corpora hold a record of the shape: `GRCh38-vcf-edge-cases`, whose runs pass `--allele_number`, and `GRCh38-vcf-minimal`, whose runs pass `--minimal` and split it.

## 1. An alternate written after a repeated one

**In short.** The record `dup_alt` of the `GRCh38-vcf-edge-cases` corpus writes three alternates on one base, and its second repeats its first. Under `--allele_number`, Ensembl VEP 116.3 numbers the row of the third alternate 2 and vep-rs numbers it 3. Two made-up records show the other forms: an alternate equal to the reference before another one, and a repeat after a `*`.

### 1.1 What a reader needs to know first

A VCF record lists its alternate alleles in the ALT column, separated by commas. The order is part of the record: a genotype names an allele by its place, 0 for the reference, 1 for the first alternate, 2 for the second.

`--allele_number` adds `ALLELE_NUM` to every row: an `Extra` item in the default format, a column in the tab format, a `CSQ` field in the VCF format and the key `allele_num` in the JSON format ([Default and tab output](../output/vep-and-tab.md)). `--uploaded_allele` adds the record's allele list, reference first, as `UPLOADED_ALLELE`.

`--minimal` trims the bases an alternate shares with the reference. Where at least one alternate of a record with several trims, the record is split into one feature for each alternate. The flag also switches `--allele_number` on.

### 1.2 The records

Three records on one base, read against the cache of the `GRCh38-vcf-edge-cases` golden corpus (`tests/golden/116/GRCh38-vcf-edge-cases/json_cache`). The first is line 15 of the corpus's `variants.vcf`. The other two are made up for this page:

```
21	30289183	dup_alt	C	G,G,T	.	.	.
21	30289183	after_reference	C	C,T	.	.	.
21	30289183	after_star	C	*,G,G,T	.	.	.
```

`21:30289183` is the middle base of the stop codon of ENST00000416044, a reverse-strand transcript of KRTAP25-1. The alternate `G` turns the codon into a serine codon, and the alternate `T` into another stop codon. Both engines take two flags:

```
--allele_number --uploaded_allele
```

### 1.3 What each engine writes

Each engine writes five rows on ENST00000416044, and the rows agree in every column but `ALLELE_NUM`:

| Record and Allele | Consequence | Ensembl VEP | vep-rs |
| --- | --- | --- | --- |
| `dup_alt`, `G` | `stop_lost` | 1 | 1 |
| `dup_alt`, `T` | `stop_retained_variant` | 2 | 3 |
| `after_reference`, `T` | `stop_retained_variant` | 1 | 2 |
| `after_star`, `G` | `stop_lost` | 2 | 2 |
| `after_star`, `T` | `stop_retained_variant` | 3 | 4 |

The rows of one record carry one `UPLOADED_ALLELE`: `C/G/G/T`, `C/C/T` and `C/*/G/G/T`. Neither engine writes a row for the repeated `G`, for the alternate `C` that equals the reference, or for the `*`.

Under `--minimal` the engines agree on a record the flag splits: Ensembl VEP then builds one feature for each written alternate and numbers each by its place. The record `multi_dup_plus_ins` of the `GRCh38-vcf-minimal` corpus repeats its first alternate, and its third alternate prints 3 in both engines. The flag does not split `dup_alt`, none of whose alternates trims, so its rows keep the numbers of the table.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP keeps every written alternate in the record's allele string, `C/G/G/T`. When it builds the alleles it annotates, it first makes the entries of that string unique. The reference is the first entry it sees, so a later entry equal to the reference is dropped, and so is the second copy of a repeated alternate.

It then counts what is left. For `dup_alt` the list left is `C`, `G`, `T`: `G` takes 1 and `T` takes 2. For `after_reference` the list left is `C`, `T`, and `T` takes 1. A `*` is not dropped: it takes its number and writes no row, so `after_star` numbers `G` 2 and `T` 3.

Source (`ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/`):

- `VariationFeatureOverlap.pm`, lines 459-467: the unique filter.
- `VariationFeatureOverlap.pm`, lines 477-491: the count. Line 484 leaves a `*` its number and gives it no row.

### 1.5 Why vep-rs's row is the correct one

The ALT column of `dup_alt` reads `G`, `G`, `T`. The row in question is the row of `T`, and `T` stands once in that column, as its third entry. Its number is 3. The same row prints `UPLOADED_ALLELE=C/G/G/T`, whose entry 2 is `G`.

Ensembl VEP states the same meaning in four places. Its header line calls the number the allele number from input. Its accessor documents it as the order of the allele in the feature's allele string (`BaseVariationFeatureOverlapAllele.pm`, lines 217-218). Its `HGVSg` lookup indexes the written string with it (`OutputFactory.pm`, lines 1386-1390). And under `--minimal` it numbers a split record's alternates by their places (`InputBuffer.pm`, lines 656 and 684).

In `after_reference` the alternate `T` is the second written alternate, so its number is 2. In `after_star` the alternates `G` and `T` stand at places 2 and 4.

### 1.6 Try it yourself

```
vep -i records.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 \
  --allele_number --uploaded_allele
```

vep-rs writes the numbers of 1.3. Ensembl VEP 116.3 takes the same two flags with the cache options of [Reproducing these examples](reproducing.md) and writes its own.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/allele_number_repeated_or_reference_equal_alt.toml` runs the three records of 1.2 on both engines and records the five rows of 1.3 for each. A further run of it gives `dup_alt` to both engines under `--minimal` and records 2 for Ensembl VEP and 3 for vep-rs on its `T` row.

The `GRCh38-vcf-edge-cases` corpus names the class in its `masks`. Six of its records write a `T` row on ENST00000416044 at this position. Two of them repeat an alternate: `dup_alt` and `trailing_tabs`, the same record under another name. The manifest documents the `ALLELE_NUM` of those two rows in one entry, with both engines' values, for each of the corpus's five runs. The entry states `rows_alike` 4, the count of rows on its key on which both engines print one number. So the golden tests hold the two rows that differ to vep-rs's number, compare the other four like any row, and hold the count of those four.

The format comparator compares two outputs row by row. Given a corpus's manifest, it holds the rows the manifest documents as the golden tests do. On a pair no manifest documents, with or without a manifest, it accepts an `ALLELE_NUM` pair under the class only where three things hold:

- both rows print the same `UPLOADED_ALLELE`;
- vep-rs's number is the place of the row's allele in that list;
- Ensembl VEP's number is that place less the repeated and reference-equal entries before it.

On an insertion or deletion record a row spells its allele without the first base the record's alleles share, so the comparator finds the allele's place by the list's entries with that base taken off. An allele that fits one entry as written and another with that base taken off, at two places, leaves the pair a difference. So does a pair on the entry that equals the reference: neither engine writes a row for that entry.

A pair on rows that print no list, as in a run without `--uploaded_allele`, stays a difference. So does any other difference of the field.
