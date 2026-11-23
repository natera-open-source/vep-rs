# A dbscSNV file whose header lacks the `Ensembl_gene` column

**The defect in one sentence.** Ensembl VEP's dbscSNV plugin dies as it loads when its file's header lacks the `Ensembl_gene` column, a column no line of the plugin reads, so the run goes on without the plugin and no row carries a score. vep-rs reads the file and writes `ada_score` and `rf_score`.

**Why this is a defect.** A check that a column is present protects a read of that column, and the plugin makes none: it matches a file row on its position and alternate allele and writes that row's `ada_score` and `rf_score` cells. The test that compared the row's `Ensembl_gene` cell with the transcript's gene is commented out, and the gene identifier fetched for it is read nowhere else. The check outlived the test: a file cut down to the six columns the plugin reads holds everything a score needs. It is turned away for a column whose content changes nothing, since a file whose `Ensembl_gene` cells all hold a made-up identifier gets the whole file's scores.

**What vep-rs does instead.** vep-rs requires the columns it reads: `alt`, and the assembly's position column, `pos` under GRCh37 and `hg38_pos` under GRCh38. It reads `ada_score` and `rf_score` by name where the header carries them. A header without `Ensembl_gene` is read, and each record takes the scores of the first file row at its position with its alternate allele, as it does from the whole file. The plugin's parameters and fields are described under [dbscSNV](../plugins/annotation-store.md#dbscsnv).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-vep` at release/116.3):

- `dbscSNV.pm`, lines 102-108: the header is the last line the plugin's `tabix -fh` call on the file prints, split on whitespace.
- Lines 110-113: `foreach my $h(qw(alt Ensembl_gene))` dies with `ERROR: Could not find required column $h in $file` when the header lacks either name. The comment above the loop names a third column, `Ensembl_transcriptid`.
- Line 158, `run`: the gene of the row's transcript is fetched into `$g_id`.
- Lines 163-173: a file row matches when its position cell equals the variant's start and its `alt` cell equals the allele. The two lines that would test `$tmp_data->{Ensembl_gene}` against `$g_id`, 169 and 170, are commented out, and no other line reads `$g_id`.
- Lines 120-123 and 179-183: of the matched row's cells the plugin returns `ada_score` and `rf_score`, each unless it is `.`.
- Line 196, `parse_data`: every cell is stored under its header name. No line reads the `Ensembl_gene` key.
- `Runner.pm`, lines 946-954: a plugin whose constructor dies is reported as `Failed to instantiate plugin` in the warnings file and skipped; under `--safe` the message ends the run.

**What correcting it would change.** The plugin would require the columns it reads and no other, and a file without `Ensembl_gene` would load and score as the whole file does.

**The class it produces.** `dbscsnv_file_without_gene_column`: the `ada_score` and `rf_score` fields, and their two header lines, that vep-rs writes and Ensembl VEP does not when the dbscSNV file's header lacks `Ensembl_gene`. No comparison tool sets a row of the shape aside: no committed flag set passes such a file, so a disagreement of this shape would be reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the fields.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a dbscSNV file whose header lacks the column, which the plugin's own check names required. The two committed made-up files each name eighteen columns, `Ensembl_gene` among them, and each of the 4 flag sets that run the plugin passes one of them. The rows the shape reaches are the transcript rows of every single-base substitution the file scores: on the 7 records of the GRCh37 corpus, 57 of 67 rows.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A splice donor substitution against a file of six columns

**In short.** `21 10942709 splice_donor_a A C` on GRCh37 is a single-base substitution at a splice donor of ENST00000298232. Against the corpus's made-up file both engines write `ada_score=0.644219` and `rf_score=0.058437` on that row. Against the same file cut to `chr pos ref alt ada_score rf_score`, Ensembl VEP 116.3 drops the plugin and writes neither field; vep-rs writes both.

### 1.1 What a reader needs to know first

dbscSNV writes two splice-site predictions of a single-base substitution on each of its transcript rows, `ada_score` and `rf_score`. Ensembl VEP's plugin reads them from a bgzip-compressed, tabix-indexed table whose header line names the columns, and finds each cell by its column's name. Its documentation prepares the table under the header line of the download (`dbscSNV.pm`, lines 51-59). A plugin that dies as it loads is skipped with a warning, and the run goes on.

### 1.2 The record

A single-base substitution against the Ensembl 116 GRCh37 cache, the first record of the `GRCh37-plugin-dbscsnv-fixture` corpus:

```
21	10942709	splice_donor_a	A	C	.	.	.
```

The whole file is the corpus's made-up `dbscSNV1.1_synthetic_GRCh37.txt.gz`, eighteen columns with fabricated scores. The cut file keeps its columns 1 to 4, 17 and 18 and is indexed as the plugin documents for GRCh37 (`dbscSNV.pm`, lines 58-59):

```
bgzip -dc dbscSNV1.1_synthetic_GRCh37.txt.gz | cut -f1-4,17,18 | bgzip -c > dbscSNV1.1_synthetic_six_columns_GRCh37.txt.gz
tabix -s 1 -b 2 -e 2 -c c dbscSNV1.1_synthetic_six_columns_GRCh37.txt.gz
```

Its header line and its three rows at the record's position:

```
chr	pos	ref	alt	ada_score	rf_score
21	10942709	A	C	0.644219	0.058437
21	10942709	A	G	0.987025	0.700332
21	10942709	A	T	0.784800	0.037445
```

Both engines run `--plugin dbscSNV,<file>` under `--assembly GRCh37`.

### 1.3 What each engine writes

The row on ENST00000298232 in the default format. Against the whole file, both engines:

```
splice_donor_a	21:10942709	C	ENSG00000166157	ENST00000298232	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1;ada_score=0.644219;rf_score=0.058437
```

Against the cut file, Ensembl VEP 116.3:

```
splice_donor_a	21:10942709	C	ENSG00000166157	ENST00000298232	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1
```

vep-rs:

```
splice_donor_a	21:10942709	C	ENSG00000166157	ENST00000298232	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=-1;ada_score=0.644219;rf_score=0.058437
```

The record writes 5 rows in each engine, and vep-rs writes both fields on all 5. Over the corpus's 7 records vep-rs writes a score on 57 of 67 rows, each equal to the row Ensembl VEP writes from the whole file; the 10 rows without one are the rows of the corpus's deletion. Ensembl VEP's output lacks the `## ada_score` and `## rf_score` header lines, and its warnings file holds one line:

```
WARNING: Failed to instantiate plugin dbscSNV: ERROR: Could not find required column Ensembl_gene in <file>
```

vep-rs writes no warning. Under `--safe` Ensembl VEP ends the run on that message with exit status 2 and writes no output; vep-rs writes the same rows with the flag as without it.

### 1.4 Why Ensembl VEP writes that

The header read of lines 102-108 gives six names. Line 111 walks `alt` and `Ensembl_gene`: the first is among the six, the second is not, and line 112 dies. `Runner.pm` catches the failure (lines 946-954), writes the warning and annotates without the plugin, so no row is handed to its `run`.

### 1.5 Why vep-rs's row is the correct one

The file row at 10942709 with the record's alternate allele carries `0.644219` and `0.058437`, and the position and the alternate allele are all the plugin's match reads (lines 165-168). The gene test is commented out, so no row of any file is kept or dropped by its `Ensembl_gene` cell. Ensembl VEP shows it on a file of seven columns, the six and an `Ensembl_gene` column holding the made-up identifier `MADE_UP_GENE` on every row: it writes the 57 scored rows it writes from the whole file. vep-rs's row is the row the plugin's own match gives.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh37 \
  --plugin dbscSNV,dbscSNV1.1_synthetic_six_columns_GRCh37.txt.gz
```

vep-rs writes the two header lines and `ada_score=0.644219;rf_score=0.058437` on each of the record's 5 rows. Ensembl VEP 116.3, given the same record and file, writes the 5 rows without either field and the warning above. Given the corpus's whole file, both engines write vep-rs's rows.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`tools/vep-formats-compare`) has no class of this shape, and the golden harness compares the `ada_score` and `rf_score` cells of the two dbscSNV corpora under their committed flag sets, each of which passes the corpus's eighteen-column file. A disagreement of the shape would be reported as a difference and charged to vep-rs.
