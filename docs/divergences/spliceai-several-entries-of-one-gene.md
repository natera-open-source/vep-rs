# A gene named by several entries of one SpliceAI score record

A worked example of `spliceai_multi_entry_cell`: a SpliceAI score record whose entries name one gene more than once, and each engine's rows for it.

This page continues [Two genes in one record](../plugins/spliceai.md#two-genes-in-one-record) on the SpliceAI page, which states the defect, where it is in Ensembl VEP's source and the class it produces. Tuple, pair, set aside and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 3. A gene named by several entries of one record

**In short.** `21 17620500 G A` on GRCh38 lies inside BTG3-AS1 and CXADR. The GRCh38 SpliceAI golden corpus's synthetic SNV file scores it with one score row of three entries: one of BTG3-AS1, one of CXADR and a second of BTG3-AS1. vep-rs writes both BTG3-AS1 entries on the three BTG3-AS1 rows, in the score row's order, and CXADR's entry on the CXADR row. Ensembl VEP 116.3 writes the whole cell, CXADR's entry inside it, on the BTG3-AS1 rows and no SpliceAI field on the CXADR row.

### 3.1 What a reader needs to know first

A `SpliceAI` INFO cell is a comma-separated list of entries, each naming its allele, its gene and eight values. The layout keeps no gene to one entry, so a score row can name one gene twice, with equal scores or with others. vep-rs's rule for a row of such a gene has four parts:

- The row carries every entry of the score row that names its gene, in the score row's order, and no entry of another gene.
- `SpliceAI_pred`, the one field of the default and tab formats, holds those entries as the score row spells them. The first stands without its allele and each later one with it, joined by commas.
- Each of the nine `SpliceAI_pred_` fields of the VCF and JSON formats and of `split_output` holds its values over those entries, joined by commas. A CSQ field writes the comma as `&`, and a JSON member of several values is one text.
- `SpliceAI_cutoff` is one verdict a row: `PASS` where any of those entries reaches the cutoff by its own largest delta score, `FAIL` where none does.

A gene that one entry names keeps the row it had. Where several score rows stand at an allele, the last one that names the gene is the row's score row.

### 3.2 The record

The record `snv_gene_around_other_gene` of the GRCh38 SpliceAI golden corpus, against the corpus's own synthetic files (`tests/golden/116/GRCh38-plugin-spliceai-fixture/plugin_data/`), and its score row in the SNV file:

```
21	17620500	snv_gene_around_other_gene	G	A	.	.	.
```

```
21	17620500	.	G	A	.	.	SpliceAI=A|BTG3-AS1|0.09|0.20|0.26|0.04|37|6|-32|31,A|CXADR|0.48|0.25|0.72|0.29|-2|33|26|-29,A|BTG3-AS1|0.19|0.06|0.14|0.13|-36|-1|-26|28
```

Each SpliceAI golden corpus carries six made-up records of this kind, one hundred bases apart, in the span its two genes share: 17620100 to 17620600 on GRCh38 (BTG3-AS1, CXADR) and 15482100 to 15482600 on GRCh37 (AP001347.6, LIPI). The table names the first gene of each pair `G` and the second `H`.

| Record | Entries of its score row | Flags under `cutoff=0.5` |
| --- | --- | --- |
| `snv_gene_two_entries_equal` | `G`, `G`, equal scores | `G` rows: `PASS` in both engines |
| `snv_gene_two_entries_differ` | `G`, `G`, other scores | `G` rows: vep-rs `PASS`, Ensembl VEP `FAIL` |
| `snv_gene_four_entries` | `G` four times | `G` rows: vep-rs `PASS`, Ensembl VEP `FAIL` |
| `snv_gene_twice_then_other_gene` | `G`, `G`, `H` | `G` rows: `PASS` in both; `H` row: vep-rs `FAIL`, Ensembl VEP none |
| `snv_gene_around_other_gene` | `G`, `H`, `G` | `G` rows: `FAIL` in both; `H` row: vep-rs `PASS`, Ensembl VEP none |
| `snv_other_gene_then_gene_twice` | `H`, `G`, `G` | `G` rows: vep-rs `PASS`, Ensembl VEP none; `H` row: `FAIL` in both |

Both engines run `--plugin SpliceAI,snv=<snv file>,indel=<indel file>` with the corpus's reference FASTA.

### 3.3 What each engine writes

The record has four transcript rows: `ENST00000626994`, `ENST00000701920` and `ENST00000807308` of BTG3-AS1, and `ENST00001128578` of CXADR. Ensembl VEP 116.3 writes the score row's whole cell on a BTG3-AS1 row and no SpliceAI field on the CXADR row:

```
snv_gene_around_other_gene	21:17620500	A	ENSG00000280594	ENST00000626994	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SpliceAI_pred=BTG3-AS1|0.09|0.20|0.26|0.04|37|6|-32|31,A|CXADR|0.48|0.25|0.72|0.29|-2|33|26|-29,A|BTG3-AS1|0.19|0.06|0.14|0.13|-36|-1|-26|28
snv_gene_around_other_gene	21:17620500	A	ENSG00000154639	ENST00001128578	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
```

vep-rs writes the BTG3-AS1 row both BTG3-AS1 entries and the CXADR row CXADR's:

```
snv_gene_around_other_gene	21:17620500	A	ENSG00000280594	ENST00000626994	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SpliceAI_pred=BTG3-AS1|0.09|0.20|0.26|0.04|37|6|-32|31,A|BTG3-AS1|0.19|0.06|0.14|0.13|-36|-1|-26|28
snv_gene_around_other_gene	21:17620500	A	ENSG00000154639	ENST00001128578	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SpliceAI_pred=CXADR|0.48|0.25|0.72|0.29|-2|33|26|-29
```

In the VCF format Ensembl VEP 116.3 writes the BTG3-AS1 rows the first entry's nine values, the ninth as `31&A`. vep-rs writes each field its two values:

```
SpliceAI_pred_SYMBOL=BTG3-AS1&BTG3-AS1  SpliceAI_pred_DS_AG=0.09&0.19  SpliceAI_pred_DP_DL=31&28
```

In JSON Ensembl VEP's `DP_DL` is the text `31,A` beside eight members of one value. vep-rs's members are `"SYMBOL":"BTG3-AS1,BTG3-AS1"`, `"DS_AG":"0.09,0.19"` and so on to `"DP_DL":"31,28"`. Under `cutoff=0.5` both engines flag the BTG3-AS1 rows `FAIL`. vep-rs flags the CXADR row `PASS`, where Ensembl VEP writes no flag.

The record `snv_gene_two_entries_differ` (`21 17620200 C T`) has a score row that names BTG3-AS1 alone, in two entries of other scores. Both engines write `SpliceAI_pred` as the score row's value behind its first allele, byte for byte. Their nine fields differ as above, and under `cutoff=0.5` Ensembl VEP flags the rows `FAIL` from the first entry's largest score, 0.49. vep-rs flags them `PASS`: the second entry's largest score is 0.84.

### 3.4 Why Ensembl VEP writes that

`SpliceAI.pm`, at the `VEP_plugins` release/116 branch, splits a score row's cell once into its first allele, its first gene and everything after it (lines 328-332). It files that remainder under the first gene (line 305). The BTG3-AS1 rows therefore print the whole remainder, CXADR's entry included. A CXADR transcript finds nothing under its symbol (lines 313-316).

For the nine fields the plugin splits the remainder on `|` and keeps the first nine parts (lines 276-287). The ninth runs on to the next entry's allele, and every later entry is dropped. The cutoff reads the largest of the remainder's first four scores (lines 336-338), the first entry's.

### 3.5 Why vep-rs's row is the correct one

Each entry of the score row names its own allele and gene. Two of the three name BTG3-AS1 and one names CXADR. A BTG3-AS1 row therefore carries the two BTG3-AS1 entries, in the order the score row gives them, and none of CXADR's scores. The CXADR row carries CXADR's entry.

A row's flag follows from the same entries: the score row gives BTG3-AS1 no score that reaches 0.5 and gives CXADR 0.72. Where a score row names one gene alone, the one field vep-rs writes is the value Ensembl VEP writes.

### 3.6 Try it yourself

```
vep -i <corpus>/variants.vcf -o out.txt --offline \
  --json_cache <corpus>/json_cache --assembly GRCh38 \
  --fasta <corpus>/reference.fa.gz \
  --plugin SpliceAI,snv=<snv file>,indel=<indel file>,cutoff=0.5
```

On the rows of `snv_gene_around_other_gene` vep-rs writes the two values above, flagged `FAIL` and `PASS`. Ensembl VEP 116.3, run with the same plugin words, writes the whole cell on the BTG3-AS1 rows alone.

### 3.7 How the comparator treats this class

The comparison tools read every score row at a differing row's allele, in file order. They set the row aside under `spliceai_multi_entry_cell` where each of these holds:

- An entry there names the row's gene, and a score row that holds such an entry holds several entries.
- vep-rs writes every entry of the last score row that names the gene which names it, in that score row's order: as the one field, or as the nine fields with each field's values joined by commas.
- Ensembl VEP writes the whole value of the last score row the gene's entry leads, or no SpliceAI field where it leads none.
- The two cutoff flags are equal, or one cutoff gives both.

The tools do not know a run's cutoff. They accept Ensembl VEP's `FAIL` beside vep-rs's `PASS` where the largest delta score of the first entry of Ensembl VEP's score row is below the largest among the gene's entries. A cutoff between the two gives both flags. A row on which vep-rs writes one entry of a gene that several name stays a difference.

Where a score row names the row's gene alone, both engines write one `SpliceAI_pred` and the default and tab formats show no difference. The corpus builder's `classify` then asks the class about the row's nine fields, which the VCF and JSON formats write. Each corpus's manifest documents the rows of the six records under the class.
