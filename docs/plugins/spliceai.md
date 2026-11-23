# SpliceAI

SpliceAI adds deep-learning splice-site predictions to a transcript consequence.
It is an annotation-store plugin: it reads the bgzipped, tabix-indexed SpliceAI
score VCFs and matches a record by position, allele and gene symbol. This page
covers its parameters, the data files and the shape its fields take in each
output format. It also covers the places where it differs from Ensembl VEP's
plugin.

The plugin's documentation lists the Illumina SNV and indel score files for
both assemblies on BaseSpace (`SpliceAI.pm` lines 101-112). They are for
academic use only under Illumina's terms, so the project does not redistribute
them. Ensembl's GRCh38 MANE SNV file (terms not stated by the plugin) is passed
under `snv=`.

## Parameters and data {#spliceai}

| Property      | Value                                                                                            |
| ------------- | ------------------------------------------------------------------------------------------------ |
| Builds        | GRCh37 and GRCh38, the file of the run's assembly                                                |
| Matches       | position, allele and gene symbol                                                                 |
| Output fields | `SpliceAI_pred` in the default and tab formats; `SpliceAI_pred_SYMBOL`, `SpliceAI_pred_DS_AG`, `SpliceAI_pred_DS_AL`, `SpliceAI_pred_DS_DG`, `SpliceAI_pred_DS_DL`, `SpliceAI_pred_DP_AG`, `SpliceAI_pred_DP_AL`, `SpliceAI_pred_DP_DG`, `SpliceAI_pred_DP_DL` in the VCF format and under `split_output`; `SpliceAI_cutoff` under `cutoff` |
| Data source   | Illumina BaseSpace (both assemblies, SNV and indel files); the [Ensembl MANE release](https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/) (GRCh38 SNV) |

| Parameter          | Meaning                                                                                       |
| ------------------ | --------------------------------------------------------------------------------------------- |
| `snv=<path>`       | The SNV score VCF.                                                                            |
| `indel=<path>`     | The indel score VCF. At least one of the two files is given; both are read for every record. |
| `cutoff=<number>`  | A value in 0 to 1 as Ensembl VEP reads it; `cutoff=` sets no cutoff. Adds `SpliceAI_cutoff`: `PASS` when the largest delta score of the row's entries reaches the cutoff, `FAIL` otherwise. No record is dropped. `cutoff=0` adds the flag too, `PASS` on every scored record, as `cutoff=0.0` does in both engines. Ensembl VEP 116.3 accepts the literal `0` in its range check and then reads it as no cutoff, since the string `0` is false in Perl (`SpliceAI.pm`, lines 165-170, 222, 295 and 335 of the `VEP_plugins` release/116 branch): it writes neither the `SpliceAI_cutoff` header line nor the flag. No committed run passes the value (tier 3 of [Intended divergences from Ensembl VEP](../intended-divergences.md); 8 flag sets run the plugin). |
| `split_output=1`   | Print the nine fields in every output format.                                                 |

Every parameter is `key=value`; a file given without its key is refused, as
Ensembl VEP's plugin refuses it.

```bash
--plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.hg38.vcf.gz,\
indel=/path/to/spliceai_scores.raw.indel.hg38.vcf.gz
--plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.hg38.vcf.gz,\
indel=/path/to/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5
--plugin SpliceAI,\
snv=/path/to/spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz
```

**Preparing the data.** vep-rs reads each file bgzipped beside its `.tbi`
index. The plugin's documentation has the Illumina files indexed with
`tabix -p vcf` (`SpliceAI.pm` lines 114-122).

## Output shape per format

The plugin prints what Ensembl VEP's `SpliceAI.pm` prints for the run's output
format, except where this page records a difference. The plugin reads a score
row's value as
`ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`
(`SpliceAI.pm` lines 276-287 and 328-332). The parts are the gene
symbol, four delta scores (acceptor gain, acceptor loss, donor gain, donor
loss) and the four positions of the largest change relative to the variant.

- The default and tab formats print one field, `SpliceAI_pred`, holding the
  value after its allele: `SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`.
- The VCF format, whose CSQ field separates its own parts with `|`, prints the
  nine `SpliceAI_pred_*` fields, as does every format under `split_output=1`.
- The JSON format nests the nine parts under `spliceai` in each consequence
  object, their names without the `SpliceAI_pred_` prefix, the scores and
  positions as numbers. The flag sits in the same object as `SpliceAI_cutoff`.
- Parquet, a vep-rs format with no Ensembl VEP counterpart, takes the tab
  format's shape: one column per field.

Each header description opens with the tool version the first file's `##INFO`
line names, which the plugin reads at `SpliceAI.pm` lines 195-203
(`SpliceAI (v1.3)` for the committed files); a file whose header names
no version gets the descriptions without it. There Ensembl VEP 116.3 prints the
parentheses empty: `SpliceAI () predicted effect on splicing` for the merged
field. In the VCF format and under `split_output=1` it prints
`SpliceAI () gene symbol` and the eight other split descriptions the same way.
It warns
`Use of uninitialized value $spliceai_version` on standard error for each such
description. The committed score files name a version (tier 3 of
[Intended divergences from Ensembl VEP](../intended-divergences.md); 8 flag sets
run the plugin).

## Gene matching

A score file can hold scores for several genes that overlap a variant, and
Ensembl VEP can annotate several genes at one site (`SpliceAI.pm` lines 83-84).
A transcript
row takes the score row whose `SYMBOL` equals its gene symbol exactly. A row of
a gene the file does not name, and a row whose transcript has no symbol, take
nothing. When several score rows of one gene match a variant, the last one in
the file wins. The Ensembl VEP plugin's hash of score rows by gene keeps the
same one.

A row carries every entry of its score row that names its gene
([worked example](../divergences/spliceai-several-entries-of-one-gene.md)).
Where score rows at an allele stand again or differ, the row's score row is
the last one that names its gene.

## Two genes in one record

**In short.** Ensembl VEP 116.3's SpliceAI plugin reads a score row's INFO
cell as one entry. A cell that joins two genes' entries with a comma gives the
first gene the whole cell as its value and the second gene nothing. vep-rs
reads each entry under its own gene.

The committed score files declare their `SpliceAI` INFO field `Number=.`, and VCF 4.2 writes
the values of a field declared so as a comma-separated list, so one score row
can hold an entry for each gene at its site. The plugin's own description (`SpliceAI.pm`
lines 83 and 40 of `VEP_plugins` release/116) says: "SpliceAI can contain scores for
multiple genes that overlap a variant" and "The plugin works with any file
outputed by the SpliceAI tool".

Both SpliceAI golden corpora, synthetic fixtures, carry the shape in both score
files. The indel file of each writes the two entries of its two-gene deletion
as one score row. The SNV file writes the two entries of its two-gene
substitution as one score row for each alternate base.

### What each engine writes for a two-gene deletion

In `GRCh37-plugin-spliceai-fixture` the record
`21 15482000 del_two_genes GCA G` lies inside AP001347.6 and LIPI. The
synthetic score row `21 15482000 . GCA G . .` has the INFO cell
`SpliceAI=G|AP001347.6|0.03|0.12|0.15|0.11|-37|-17|-20|-13,G|LIPI|0.06|0.10|0.16|0.09|-22|33|12|35`.
The two engines write different rows.

From that synthetic cell, Ensembl VEP 116.3 writes
`SpliceAI_pred=AP001347.6|0.03|0.12|0.15|0.11|-37|-17|-20|-13,G|LIPI|0.06|0.10|0.16|0.09|-22|33|12|35`
on `ENST00000428809` and `ENST00000432621`, the AP001347.6 rows, and no
SpliceAI field on `ENST00000344577` and `ENST00000536861`, the LIPI rows. In
the VCF format the AP001347.6 rows carry `SpliceAI_pred_DP_DL` as `-13&G`, the
CSQ field's spelling of `-13,G`, and the LIPI rows carry nine empty fields. In
JSON `DP_DL` is the string `-13,G` and the LIPI rows have no `spliceai` object.
Under `cutoff=0.5` the AP001347.6 rows are flagged `FAIL` from the first gene's
scores and the LIPI rows carry no flag.

vep-rs writes the synthetic cell's
`SpliceAI_pred=AP001347.6|0.03|0.12|0.15|0.11|-37|-17|-20|-13` on
the AP001347.6 rows and `SpliceAI_pred=LIPI|0.06|0.10|0.16|0.09|-22|33|12|35`
on the LIPI rows in every format, `DP_DL` as `-13` and `35`. Under
`cutoff=0.5` it flags all four rows `FAIL`, the largest delta scores being 0.15
and 0.16.

In `GRCh38-plugin-spliceai-fixture` the record `21 17620000 del_two_genes CTT C` lies
inside BTG3-AS1 and CXADR, and its synthetic score row joins BTG3-AS1's entry and
CXADR's. Ensembl VEP 116.3 writes the whole cell on the three BTG3-AS1 rows,
`ENST00000626994`, `ENST00000701920` and `ENST00000807308`, and no SpliceAI
field on `ENST00001128578`, the CXADR row. vep-rs writes each gene's own entry,
and under `cutoff=0.5` flags the BTG3-AS1 rows `PASS` and the CXADR row `FAIL`.

### What each engine writes for a two-gene substitution

A two-gene substitution takes the same course. In `GRCh37-plugin-spliceai-fixture` the
record `21 15481500 snv_two_genes T C` lies inside the same two genes. The
synthetic score row `21 15481500 . T C . .` has the INFO cell
`SpliceAI=C|AP001347.6|0.15|0.08|0.06|0.17|-6|-41|-17|-9,C|LIPI|0.02|0.12|0.10|0.05|-25|-23|37|17`.
Ensembl VEP 116.3 writes
`SpliceAI_pred=AP001347.6|0.15|0.08|0.06|0.17|-6|-41|-17|-9,C|LIPI|0.02|0.12|0.10|0.05|-25|-23|37|17`
on the AP001347.6 rows and no SpliceAI field on the LIPI rows. vep-rs writes
`SpliceAI_pred=AP001347.6|0.15|0.08|0.06|0.17|-6|-41|-17|-9` on the AP001347.6
rows and `SpliceAI_pred=LIPI|0.02|0.12|0.10|0.05|-25|-23|37|17` on the LIPI
rows.

In `GRCh38-plugin-spliceai-fixture` the record `21 17619500 snv_two_genes A G` has a
synthetic score row that joins BTG3-AS1's entry and CXADR's. Ensembl VEP 116.3 writes
the whole cell on the three BTG3-AS1 rows and no SpliceAI field on the CXADR
row. vep-rs writes each gene's own entry, and under `cutoff=0.5` flags the
BTG3-AS1 rows `FAIL` and the CXADR row `PASS`.

### Why Ensembl VEP gives the first gene the whole cell

`SpliceAI.pm`, lines 328-332 of the `VEP_plugins` release/116 branch, removes
`SpliceAI=` from the whole INFO cell. It splits the rest on `|` into at most
three parts: the allele, the gene and everything after the gene. The value it
prints is the gene joined to that remainder. The synthetic cell's second entry
sits inside the value's ninth part, `-13,G`. The parts after it are dropped when the
value is split into nine fields (lines 276-287).

The score row is stored under the first entry's gene alone (line 305). A LIPI
transcript therefore finds no score row (lines 313-316). A cutoff reads the
first gene's four scores (lines 336-338).

An INFO key written after the `SpliceAI` value lands inside the ninth part the
same way. Ensembl VEP 116.3 prints the key with the ninth value, its
`;` escaped (`SpliceAI.pm` lines 328-331 and 287); vep-rs, reading the
`SpliceAI` field alone, prints the value.

### Why each gene's own entry is the correct value

Each entry of the cell names its own allele and gene, and the `Number=.`
declaration admits the list. The plugin's own description, quoted above, says
a file can hold scores for several genes. The score row's own values give
AP001347.6 the first entry and LIPI the second, and those are the rows vep-rs
writes. A site scored by one score row per gene gives identical rows in both
engines, as each corpus's record `snv_two_genes_record_per_gene` shows.

### How the comparison tools treat a cell of several entries

The shape is a property of the score file, not of an input record, so the
comparison tools read the files themselves. The format comparator, the plugin
scorer and the corpus builder's `classify` take the run's two score files after
`--spliceai-table`. They read them through their tabix indexes at a differing
row's variant.

They set the row aside under `spliceai_multi_entry_cell` only where each of
these holds:

- An entry there names the row's gene.
- A score row that holds such an entry holds several entries.
- vep-rs writes the gene's entries of the last score row that names it.
- Ensembl VEP writes the value of the last score row the gene's entry leads,
  that score row's cell behind its first allele. One cutoff gives both
  flags. Where the gene's entry leads no score row, it writes no SpliceAI
  field.

The row's gene is its `SYMBOL` where the row writes one, else its transcript's
symbol in the candidate cache. The comparator and the scorer take the cache
after `--json-cache`, and `classify` reads the corpus's own. A row whose gene
neither its `SYMBOL` nor the cache names is not set aside. Without the files
nothing is set aside.

### How often a cell holds several entries

86 of the 177 rows of the four committed score files hold several entries.
The four files are synthetic, and a row of a synthetic table is not ordinary
data, so the entry stands at tier 3 of
[Intended divergences from Ensembl VEP](../intended-divergences.md); 8 flag
sets run the plugin and both corpora carry such rows. Each corpus's
manifest documents the four rows of `del_two_genes`, of `snv_two_genes`, of
`snv_two_genes_line_twice`, of `snv_two_genes_line_thrice`, of
`snv_two_genes_lines_differ` and of the worked example's six records under it,
and the plugin gate's eight SpliceAI runs hand the scorer the score
files. The
[index of differences](../intended-divergences.md#differences-outside-the-classes)
lists it.

### 2. A record line that stands more than once at an allele

**In short.** `21 17619600 A G` on GRCh38 lies inside BTG3-AS1 and CXADR. The
GRCh38 SpliceAI golden corpus's synthetic SNV file holds the score row of the two
genes' entries twice for that allele. vep-rs writes each gene's rows its own
entry once. Ensembl VEP 116.3 writes the whole cell once on the three BTG3-AS1
rows and no SpliceAI field on the CXADR row, as it does where the score row
stands once.

#### 2.1 What a reader needs to know first

The plugin reads every score row that the tabix query returns for a variant,
in file order. Each SpliceAI golden corpus's SNV file holds some score rows
more than once at an allele. Each of those score rows joins several genes'
entries with commas in one INFO cell, the shape
[Two genes in one record](#two-genes-in-one-record) describes.

#### 2.2 The record

The record `snv_two_genes_line_twice` of the GRCh38 SpliceAI golden corpus,
against the corpus's own synthetic files
(`tests/golden/116/GRCh38-plugin-spliceai-fixture/plugin_data/`). The SNV file
holds the record's score row twice, one copy after the other:

```
21	17619600	snv_two_genes_line_twice	A	G	.	.	.
```

```
21	17619600	.	A	G	.	.	SpliceAI=G|BTG3-AS1|0.01|0.04|0.16|0.19|-38|17|-47|-16,G|CXADR|0.15|0.01|0.02|0.14|47|-44|26|-12
21	17619600	.	A	G	.	.	SpliceAI=G|BTG3-AS1|0.01|0.04|0.16|0.19|-38|17|-47|-16,G|CXADR|0.15|0.01|0.02|0.14|47|-44|26|-12
```

Both engines run `--plugin SpliceAI,snv=<snv file>,indel=<indel file>` with the
corpus's reference FASTA.

#### 2.3 What each engine writes

The record has four transcript rows: `ENST00000626994`, `ENST00000701920` and
`ENST00000807308` of BTG3-AS1, and `ENST00001128578` of CXADR. Ensembl VEP
116.3 writes the score row's whole cell on a BTG3-AS1 row and no SpliceAI field
on the CXADR row:

```
snv_two_genes_line_twice	21:17619600	G	ENSG00000280594	ENST00000626994	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SpliceAI_pred=BTG3-AS1|0.01|0.04|0.16|0.19|-38|17|-47|-16,G|CXADR|0.15|0.01|0.02|0.14|47|-44|26|-12
snv_two_genes_line_twice	21:17619600	G	ENSG00000154639	ENST00001128578	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
```

vep-rs writes each row its gene's entry, once:

```
snv_two_genes_line_twice	21:17619600	G	ENSG00000280594	ENST00000626994	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SpliceAI_pred=BTG3-AS1|0.01|0.04|0.16|0.19|-38|17|-47|-16
snv_two_genes_line_twice	21:17619600	G	ENSG00000154639	ENST00001128578	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SpliceAI_pred=CXADR|0.15|0.01|0.02|0.14|47|-44|26|-12
```

The record `snv_two_genes_line_thrice` (`21 17619800 A G`) has a synthetic score row of three
entries that stands three times: BTG3-AS1's, CXADR's and one naming `DECOY`, a
symbol no transcript carries. Ensembl VEP 116.3 writes the three-entry cell
once on the BTG3-AS1 rows and nothing on the CXADR row. vep-rs writes
`SpliceAI_pred=BTG3-AS1|0.63|0.75|0.78|0.67|-30|-31|7|32` on the BTG3-AS1 rows
and `SpliceAI_pred=CXADR|0.11|0.05|0.08|0.16|-2|16|-14|24` on the CXADR row.
The `DECOY` entry reaches no row in vep-rs, and in Ensembl VEP it stands inside
the cell the BTG3-AS1 rows carry. The GRCh37 corpus's records of the same names
behave the same way on AP001347.6 and LIPI.

The record `snv_two_genes_lines_differ` (`21 17619900 G A`) is scored by two
score rows of the same two genes whose scores differ. Ensembl VEP 116.3 writes
the second score row's whole cell on the BTG3-AS1 rows and nothing on the CXADR
row. vep-rs writes each gene's rows its entry of the second score row.

#### 2.4 Why Ensembl VEP writes that

`SpliceAI.pm`, at the `VEP_plugins` release/116 branch, reads each score row's
cell as one entry (lines 328-332). It files the entry under the score row's
first gene, with everything after the first allele as its value. Line 305 keeps
one value for each gene in a hash. The second copy of the score row therefore
takes the first copy's place with the same value, and the row prints that value
once. A CXADR transcript finds no value under its symbol (lines 313-316),
because CXADR leads no score row.

#### 2.5 Why vep-rs's row is the correct one

Each copy of the score row states the same two entries, one naming BTG3-AS1 and
one naming CXADR, each with its own allele and scores. Reading the score row
once or twice gives each gene the same entry. The BTG3-AS1 rows therefore take
the first entry and the CXADR row takes the second, as vep-rs writes them. A
copy of a score row adds no entry and no gene.

#### 2.6 Try it yourself

```
vep -i <corpus>/variants.vcf -o out.txt --offline \
  --json_cache <corpus>/json_cache --assembly GRCh38 \
  --fasta <corpus>/reference.fa.gz \
  --plugin SpliceAI,snv=<snv file>,indel=<indel file>
```

On the rows of `snv_two_genes_line_twice` vep-rs writes the two values above.
Ensembl VEP 116.3, run with the same plugin words, writes the whole cell on the
BTG3-AS1 rows alone.

#### 2.7 How the comparator treats this class

The comparison tools read every score row at a differing row's allele, in file
order. They set the row aside under `spliceai_multi_entry_cell` where an entry
names the row's gene. A score row that holds such an entry must hold several
entries. vep-rs must write the gene's entries of the last such score row.
Ensembl VEP's value must be that of the last score row the gene's entry leads,
or no field where it leads none.

A row whose gene no entry names stays a difference, and so does a row on which
vep-rs writes an earlier score row's entry. When the score files are given, the
scorer's summary counts the rows the class was asked about. It also counts the
rows the class did not set aside, by the ground of each refusal.

## Indels and the reference FASTA

Ensembl VEP's plugin compares an indel through a VCF record it rebuilds from the
variant, which needs the base before the change. It reads that base from the
reference (`BaseVariationFeature.pm` `_get_prev_base`, `N` without one). A run
without `--fasta` therefore builds a record no score row equals and scores no
indel, and a run with `--fasta` scores them. vep-rs compares the input record's
own anchor base and scores indels with or without a reference FASTA. The golden corpora run both engines with
a reference, so the comparison covers the indel records.

## Reverse-strand indels

An Ensembl-format or region record given on strand -1 is matched to the score
file on its forward-strand alleles, as every tabix-backed plugin matches it. A
reverse-strand indel therefore scores the row its forward-strand twin scores.
Ensembl VEP 116.3's SpliceAI plugin complements such an indel twice. It reads
the forward-strand alleles `to_VCF_record` writes (`VariationFeature.pm` line
2376), then reverse-complements them again (`SpliceAI.pm` lines 259-271 of the
`VEP_plugins` release/116 branch). The alleles it compares are the
reverse-strand ones, and no row of the file matches.

Take `21 26021981 26021980 -/CTT -1` against a made-up score file holding the forward-strand
insertion `G GAAG` at 26021980. Ensembl VEP 116.3 writes no SpliceAI field on
any of the record's 57 rows, and vep-rs scores that row on every transcript row.

A reverse-strand SNV is complemented once by vep-rs, by Ensembl VEP's SpliceAI
(`SpliceAI.pm` lines 268-271) and by every Ensembl VEP plugin that matches
through `get_matched_variant_alleles`. Ensembl VEP's ClinPred compares the
record's as-written base with the table row's alternate allele (`ClinPred.pm`
lines 109 and 117) and scores another substitution's row. The
[index of differences](../intended-divergences.md#differences-outside-the-classes)
lists it.

## A substitution beside a deletion in one record

Both engines read a VCF record whose alternate alleles all share the reference's
first base and include a deletion or insertion, such as `GC G,GA`, with that
base removed. The example becomes `C/-/A` one position on (`Bio/EnsEMBL/VEP/Parser/VCF.pm`, lines
300-323). Ensembl VEP 116.3's SpliceAI plugin, seeing the `-`, rebuilds the
whole record with its anchor base before matching (`SpliceAI.pm`, lines 259-266,
VEP_plugins at release/116; `Bio/EnsEMBL/Variation/VariationFeature.pm`, lines
2400-2420). The substitution allele is then compared as `GC GA` at the anchor's
position, which a score row that writes it as `C A` one base later does not
equal (`SpliceAI.pm` line 272). A reference FASTA changes nothing: the rebuilt alleles stay two bases
long.

vep-rs compares each allele at its own position against the score row trimmed
of its shared leading base. The substitution therefore scores, and the deletion
scores when the indel file holds its row.

For `21 26965499 GC G,GA` against the GRCh38 SpliceAI golden corpus's synthetic files,
Ensembl VEP 116.3 writes no SpliceAI field on any of the record's twelve rows.
vep-rs writes `SpliceAI_pred=ADAMTS5|0.05|0.71|0.72|0.66|-47|-14|3|17`, the
file's `26965500 C A` row, on the six `A` rows. The deletion rows carry nothing
in both engines, since the indel file holds no row for it. A record whose
alleles share no first base, such as `C A,G` or `A G,AC`, is kept unchopped and
matches alike in both engines.

The shape is rare. 1 of the 15,646,481 records of the chromosome 21
cuts is a multi-allelic record whose alleles share the anchor and, chopped, mix
an empty allele with a substitution. That is tier 2 of
[Intended divergences from Ensembl VEP](../intended-divergences.md), and 8 flag
sets run the plugin. The
[index of differences](../intended-divergences.md#differences-outside-the-classes)
lists it.

The format comparator (`tools/vep-formats-compare`) applies the
class `spliceai_chopped_multiallelic` when a comparison case names it. It sets
aside a `SpliceAI_` field, or the JSON `spliceai` object, where vep-rs alone
fills it on a substitution allele's row. The record's rows must carry an allele
written `-`. Any other difference stays one, a field on the row of
the `-` allele itself among them: the class covers the substitution's rows
alone.

## A score row spelled with a shared base

**The defect in one sentence.** Ensembl VEP's SpliceAI plugin compares a score
row's REF and ALT as written, so a row whose alleles share a leading or trailing
base, apart from the one base anchoring an insertion or deletion, scores nothing
for a record spelling the change without it.

**Why this is a defect.** A base the REF and ALT share at the same end is a
reference base the change leaves as it is. `TG` to `AG` at 18353100 and `T` to
`A` at 18353100 produce one sequence. The reference carries `G` at
18353101 and the row keeps it.

Ensembl VEP 116.3 scores the row for a record spelled as the row is, so the
row's values are the values of that change. The same change spelled with its
one changed base finds no row. The plugin
resolves a second spelling for an insertion or deletion, whose record it
rebuilds with its anchor base before comparing, and for no substitution.

**What vep-rs does instead.** vep-rs trims the bases a row's REF and ALT share
at either end before comparing them with the record's alleles. It moves the
row's position past the trimmed leading bases. It scores the row when the
trimmed alleles equal the record's at its position. The row's values are written as
they stand, as both engines write them for a record spelled as the row is.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):
`SpliceAI.pm`, line 272, takes a row when its start, REF and ALT each equal the
record's with `==` and `eq`. The REF and ALT are the row's fourth and fifth
columns as `parse_data` splits the row (lines 326-333). Lines 259-266 rebuild
only an allele string carrying `-`, through `to_VCF_record`, so a substitution
is compared in the record's own spelling.

**What correcting it would change.** A record would take a row that spells its
change beside a shared base, as it takes the row spelled without one.

### The class for a score row with a shared base and how often it occurs

**The class it produces.** `spliceai_trimmed_table_row_match`: a `SpliceAI_pred`
field, its nine split fields or the `SpliceAI_cutoff` flag that vep-rs writes
where Ensembl VEP writes none. vep-rs writes it from a score row whose REF and
ALT share a leading or trailing base beyond an indel's anchor.

No comparison tool sets a row of the shape aside. Every committed record spells
its change as its rows do, so no committed run reaches it. A disagreement of
the shape is reported as a difference and charged to vep-rs. The concordance
figures score consequence sets and never read a plugin field.

**How often the shape occurs.** Tier 3, see
[Intended divergences from Ensembl VEP](../intended-divergences.md): 2 of the 177
rows of the four committed score files, which are synthetic, share a base at an
end, beyond the one leading base that anchors an insertion or deletion. Each is
the substitution
allele of a multi-allelic record, written as that record spells it: `TG AG`
beside `TG T` at 18353100 on GRCh38, `CA GA` beside `CA C` at 18919100 on
GRCh37. 8 flag sets run the plugin and the 2 corpora hold such a row. A file
whose rows spell every substitution with its changed bases alone holds none.

### 1. A substitution spelled beside its unchanged neighbour

**In short.** `21 18353100 T A` on GRCh38 lies in an intron of TMPRSS15, at the
site of the GRCh38 SpliceAI golden corpus's record `TG T,AG`. The corpus's synthetic
indel file scores that record's `AG` allele as the row `TG AG`, the same change with
the reference `G` at 18353101 appended. On ENST00000284885, vep-rs writes the
row's `TMPRSS15|0.82|0.64|0.86|0.67|-25|-36|-29|-20`. Ensembl VEP 116.3 writes no
SpliceAI field, and writes the row's values once the record is spelled `TG AG`.

#### 1.1 What a reader needs to know first

A score file is a VCF whose rows name a position, a REF, an ALT and the gene
the scores belong to (`SpliceAI.pm` lines 326-332). The committed indel files,
which are synthetic, hold one row per allele that is not a single base,
spelled as the corpus's record spells it, beside the record's own rows.

#### 1.2 The record

A substitution at the site of the GRCh38 SpliceAI golden corpus's record
`multi_deletion_unanchored` (`21 18353100 TG T,AG`), against the corpus's own
synthetic files (`tests/golden/116/GRCh38-plugin-spliceai-fixture/plugin_data/`). The reference is
`TG` at 18353100-18353101, and the indel file holds the record's two rows:

```
21	18353100	shared_trailing_base	T	A	.	.	.
```

```
21	18353100	.	TG	AG	.	.	SpliceAI=AG|TMPRSS15|0.82|0.64|0.86|0.67|-25|-36|-29|-20
21	18353100	.	TG	T	.	.	SpliceAI=T|TMPRSS15|0.98|0.48|0.81|0.42|-11|42|39|13
```

Both engines run `--plugin SpliceAI,snv=<snv file>,indel=<indel file>` with the
corpus's reference FASTA.

#### 1.3 What each engine writes

The record has one transcript row, ENST00000284885 of TMPRSS15. Ensembl VEP
116.3 writes it without a SpliceAI field:

```
shared_trailing_base	21:18353100	A	ENSG00000154646	ENST00000284885	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
```

vep-rs writes the `TG AG` row's values on it:

```
shared_trailing_base	21:18353100	A	ENSG00000154646	ENST00000284885	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;SpliceAI_pred=TMPRSS15|0.82|0.64|0.86|0.67|-25|-36|-29|-20
```

The record `21 18353100 row_spelling TG AG` is the same change spelled as the
row is. It gets the synthetic row's `SpliceAI_pred=TMPRSS15|0.82|0.64|0.86|0.67|-25|-36|-29|-20` on that
transcript row in both engines, the rows identical. On GRCh37, `21 18919100 C G`
against the GRCh37 corpus's `CA GA` row behaves the same way on the six CXADR
rows. They carry `SpliceAI_pred=CXADR|0.20|0.08|0.14|0.08|29|46|27|-19` in
vep-rs and no field in Ensembl VEP 116.3.

#### 1.4 Why Ensembl VEP writes that

The tabix query over 18353100 returns both rows. For each, line 272 compares
the record's start, REF and ALT with the row's. `18353100` equals the row's
start, but `T` is not `TG` and `A` is not `AG`. So neither row is taken and the
row set is empty. The `-`-allele branch of lines 259-266, which would rebuild
the record's spelling, is not entered for a substitution.

#### 1.5 Why vep-rs's row is the correct one

The row's REF `TG` and ALT `AG` share their `G`, the reference base at
18353101. The row therefore states the change of 18353100 from `T` to `A` and nothing
else. The record states that change. Ensembl VEP 116.3 scores the row when the
record repeats the row's spelling, so the values are agreed to belong to the
change. vep-rs writes them for the change in either spelling.

The four positions are the row's own, relative to its position, as both engines
write them for every row they take. A trailing base leaves that position at the
record's.

#### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --fasta <corpus>/reference.fa.gz \
  --plugin SpliceAI,snv=<snv file>,indel=<indel file>
```

vep-rs writes the synthetic row's
`SpliceAI_pred=TMPRSS15|0.82|0.64|0.86|0.67|-25|-36|-29|-20` on
the record's one row; Ensembl VEP 116.3 writes no SpliceAI field. Spelling the
record `TG AG` makes Ensembl VEP write vep-rs's row.

#### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator
(`tools/vep-formats-compare`) has no class of this shape, and
the golden harness compares the two SpliceAI corpora's fields as written, which
their records, spelled as their rows are, never put to the test.

## Where next

- [Plugins](../plugins.md): every plugin's support status and how one is
  enabled.
- [Reproducing the plugin concordance](reproducing-concordance.md): both
  engines' commands and the scorer for a supported plugin.
- [Intended divergences from Ensembl VEP](../intended-divergences.md): the
  index of every documented difference between the engines.
