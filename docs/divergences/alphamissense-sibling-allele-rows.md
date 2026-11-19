# The AlphaMissense prediction of a sibling allele on a multi-allelic record

**The defect in one sentence.** Ensembl VEP's AlphaMissense plugin matches the prediction table against every alternate allele of a record and returns the first row whose amino-acid change equals the row's own, so an allele of a multi-allelic record takes a sibling allele's prediction when both give one amino acid. vep-rs matches each allele's rows against that allele alone, so every row carries its own allele's prediction or none.

**Why this is a defect.** Each table row predicts one substitution: a chromosome, a position, a reference base, an alternate base, the protein variant and the transcript. The plugin tests the row's own amino-acid letters and, under `transcript_match=1`, its own transcript, so its contract is the prediction for this allele on this transcript, and on a one-allele record it reads exactly that row. Only the record's allele list widens the match to a sibling's row, and the table's row order then decides which allele's prediction every allele gets.

**What vep-rs does instead.** vep-rs matches each allele's rows against that allele alone, so every row carries its own allele's prediction or none, whatever the record's other alleles and whatever the table's row order. The plugin's parameters and fields are described under [AlphaMissense](../plugins/annotation-store.md#alphamissense).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-variation` at release/116):

- `AlphaMissense.pm`, line 245, `run`: the variant side of the match is `$tva->base_variation_feature->alt_alleles`, every alternate allele of the record (`VariationFeature.pm`, line 1555, returns every allele of the allele string but the reference), not the allele being annotated.
- `AlphaMissense.pm`, lines 250-253: each table row at the position is taken in file order and kept when the first and last letters of its `protein_variant` equal the row's reference and alternate amino acids (`_aminoacid_changes_match`, lines 221-236); lines 256-261 add the transcript test under `transcript_match`.
- `AlphaMissense.pm`, lines 263-275: the kept row's reference, alternate and position are matched by `get_matched_variant_alleles` (`Utils/Sequence.pm`, lines 1144-1247) against the record's reference and the whole allele list, and lines 280-283 return the first row that matches any allele of the list.

**What correcting it would change.** Each allele's rows would carry the table's prediction for that allele, as they do when the record is split into one line per allele, and swapping two table rows would change no allele's value.

**The class it produces.** `alphamissense_sibling_allele_rows`: an `am_` field (the JSON `alphamissense` object) that differs between the engines on a row of a multi-allelic record whose alleles give one amino-acid change. The format comparator sets the field aside where a comparison case names the class, under the test of section 1.7. Each AlphaMissense golden corpus holds one two-allele record, whose alleles give different amino acids, and its table carries no two rows of one site and one amino-acid change with different predictions, so no committed row has the shape. Where no case names the class, a disagreement of this shape is reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): 6,347 of the cuts' 15,646,481 records (405.6 per million) carry two or more sequence alleles; the codon decides whether two give one amino acid, and the table whether their rows differ. 64 corpora carry the shape and 6 flag sets run the plugin.

The published tables predict the protein variant, so a transcript's rows for one amino-acid change carry one value, and the rows that differ belong to two genes' transcripts overlapping at one site. The GRCh38 table's rows form 5,842,071 groups of one site and one amino-acid change with two or more alternate alleles, and one group carries two values; the GRCh37 table's form 5,688,744, and one carries two values. Neither site is on chromosome 21.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Two alleles of one codon that both give leucine, against a table whose rows for them differ

**In short.** The record `21 13610388 A C,T` on the forward-strand transcript `ENST00000299443`, against the GRCh38 AlphaMissense golden corpus's table with one row changed: the `A T` row's prediction set to `0.9000`, `likely_pathogenic`, where the `A C` row keeps `0.1079`, `likely_benign`. Both alleles give the amino-acid change `M/L`. Both engines write the `A C` row on the `C` rows. On the `T` rows Ensembl VEP 116.2 writes the `A C` row's values, `am_class=likely_benign;am_pathogenicity=0.1079`, in either ALT order; vep-rs writes the `A T` row's, `am_class=likely_pathogenic;am_pathogenicity=0.9000`. Run as two one-allele records, Ensembl VEP writes vep-rs's rows.

### 1.1 What a reader needs to know first

The published AlphaMissense tables have one row per site, alternate allele and transcript. Each row names the protein variant it predicts, `M54L`, with its score and class. A multi-allelic VCF record lists several alternate alleles on one line, and both engines write one set of rows per allele. The plugin matches a table row to an output row by the first and last letters of the protein variant, the reference and alternate amino acids. Two alleles of one codon that give the same amino acid therefore both match the rows of the other allele.

### 1.2 The record

A two-allele record at the site of the `GRCh38-plugin-alphamissense` golden corpus's table (`tests/golden/116/GRCh38-plugin-alphamissense/plugin_data/`), in both ALT orders:

```
21	13610388	sibling_ct	A	C,T	.	.	.
21	13610388	sibling_tc	A	T,C	.	.	.
```

The table's rows at the site, after the `A T` row's score and class are changed from the corpus's `0.1079`, `likely_benign` (the value the published table gives both rows) and the table is re-indexed:

```
chr21	13610388	A	C	hg38	Q86YR6	ENST00000299443.5	M54L	0.1079	likely_benign
chr21	13610388	A	T	hg38	Q86YR6	ENST00000299443.5	M54L	0.9000	likely_pathogenic
chr21	13610388	A	G	hg38	Q86YR6	ENST00000299443.5	M54V	0.0772	likely_benign
```

Both engines run `--plugin AlphaMissense,file=<table>`.

### 1.3 What each engine writes

The rows on `ENST00000299443` in the default format. Ensembl VEP 116.2 on `C,T`:

```
sibling_ct	21:13610388	C	ENSG00000166351	ENST00000299443	Transcript	missense_variant	612	160	54	M/L	Atg/Ctg	-	IMPACT=MODERATE;STRAND=1;am_class=likely_benign;am_pathogenicity=0.1079
sibling_ct	21:13610388	T	ENSG00000166351	ENST00000299443	Transcript	missense_variant	612	160	54	M/L	Atg/Ttg	-	IMPACT=MODERATE;STRAND=1;am_class=likely_benign;am_pathogenicity=0.1079
```

Ensembl VEP 116.2 on `T,C`:

```
sibling_tc	21:13610388	T	ENSG00000166351	ENST00000299443	Transcript	missense_variant	612	160	54	M/L	Atg/Ttg	-	IMPACT=MODERATE;STRAND=1;am_class=likely_benign;am_pathogenicity=0.1079
sibling_tc	21:13610388	C	ENSG00000166351	ENST00000299443	Transcript	missense_variant	612	160	54	M/L	Atg/Ctg	-	IMPACT=MODERATE;STRAND=1;am_class=likely_benign;am_pathogenicity=0.1079
```

vep-rs on `C,T`:

```
sibling_ct	21:13610388	C	ENSG00000166351	ENST00000299443	Transcript	missense_variant	612	160	54	M/L	Atg/Ctg	-	IMPACT=MODERATE;STRAND=1;am_class=likely_benign;am_pathogenicity=0.1079
sibling_ct	21:13610388	T	ENSG00000166351	ENST00000299443	Transcript	missense_variant	612	160	54	M/L	Atg/Ttg	-	IMPACT=MODERATE;STRAND=1;am_class=likely_pathogenic;am_pathogenicity=0.9000
```

vep-rs on `T,C`:

```
sibling_tc	21:13610388	T	ENSG00000166351	ENST00000299443	Transcript	missense_variant	612	160	54	M/L	Atg/Ttg	-	IMPACT=MODERATE;STRAND=1;am_class=likely_pathogenic;am_pathogenicity=0.9000
sibling_tc	21:13610388	C	ENSG00000166351	ENST00000299443	Transcript	missense_variant	612	160	54	M/L	Atg/Ctg	-	IMPACT=MODERATE;STRAND=1;am_class=likely_benign;am_pathogenicity=0.1079
```

Each record writes four rows in each engine, two per allele; the rows on `ENST00000620442` follow the same pattern. The `am_` fields by ALT order:

| ALT order | Allele | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- | --- |
| `C,T` | `C` | `likely_benign`, `0.1079` | `likely_benign`, `0.1079` |
| `C,T` | `T` | `likely_benign`, `0.1079` | `likely_pathogenic`, `0.9000` |
| `T,C` | `T` | `likely_benign`, `0.1079` | `likely_pathogenic`, `0.9000` |
| `T,C` | `C` | `likely_benign`, `0.1079` | `likely_benign`, `0.1079` |

Given alone, `A T` gets `likely_pathogenic`, `0.9000` and `A C` gets `likely_benign`, `0.1079` on every row in both engines. With the two table rows swapped, the `A T` row written before the `A C` row, Ensembl VEP writes `likely_pathogenic`, `0.9000` on the `C` rows as well as the `T` rows; vep-rs's rows do not change. Against the corpus's unchanged table, whose two rows carry one value, the engines agree in both orders.

### 1.4 Why Ensembl VEP writes that

Line 245 of `AlphaMissense.pm` takes the record's allele list, `C` and `T`, for every row of the record. On the `T` allele's row the loop at line 250 reads the table rows at the site in file order. The `A C` row comes first, its letters `M` and `L` equal the row's amino acids (line 253), and `get_matched_variant_alleles` finds its alternate `C` in the list (lines 263-275), so lines 281-282 return the `A C` row's columns. The `A T` row is never read.

In the order `T,C` the list is `T` and `C`, which holds the same two alleles, so the same row is returned. With the table rows swapped the `A T` row comes first, its alternate `T` is in the list, and both alleles' rows take it.

### 1.5 Why vep-rs's row is the correct one

The table's `A T` row is the prediction for the `T` substitution: its score and class belong to that allele, and the `A C` row's to `C`. A `T` row written with the `C` row's values states a prediction the table does not make for `T`. Ensembl VEP itself writes the `A T` row's values when `T` is the only allele on the line, and its own `C` rows keep their row in both ALT orders. vep-rs's rows inside the multi-allelic record are those one-allele rows, in either ALT order and whichever table row comes first.

### 1.6 Try it yourself

Write the table of 1.2 from the corpus's, then annotate the record with each engine:

```
zcat <corpus>/plugin_data/alphamissense.tsv.gz \
  | awk -F'\t' -v OFS='\t' '$2 == "13610388" && $4 == "T" { $9 = "0.9000"; $10 = "likely_pathogenic" } { print }' \
  | bgzip -c > am_edited.tsv.gz && tabix -s 1 -b 2 -e 2 am_edited.tsv.gz
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --dir_plugins <corpus>/plugin_data \
  --plugin AlphaMissense,file=am_edited.tsv.gz
```

Ensembl VEP's `T` rows carry `am_class=likely_benign;am_pathogenicity=0.1079` and vep-rs's `am_class=likely_pathogenic;am_pathogenicity=0.9000`, whichever ALT comes first; the `C` rows agree. Splitting the record into `A C` and `A T` makes Ensembl VEP write vep-rs's rows. Against the unchanged corpus table the engines agree, since its two rows carry one value.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies `alphamissense_sibling_allele_rows` when a comparison case names it. An `am_` field, or the JSON `alphamissense` object, is set aside, counted once per row, where Ensembl VEP's value is the one its row of another allele of the record carries on the same feature beside the same `Amino_acids`. vep-rs leaves the field empty there or fills it differently. Every other difference stays one: a field vep-rs alone fills, an Ensembl VEP value that none of its rows of another allele carries on that feature beside those amino acids, and any value on a record of one allele.

The golden harness compares the two AlphaMissense corpora's `am_` fields verbatim, which their one two-allele record each, giving different amino acids per allele, never puts to the test, so the class applies to none of their rows.
