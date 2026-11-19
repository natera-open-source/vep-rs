# A GWAS summary-statistics file with an empty cell before `hm_pos`

**The defect in one sentence.** Under `type=sstate`, Ensembl VEP 116.2's GWAS plugin filters and sorts a harmonised summary-statistics file with `awk` and `sort` counting blank-separated fields, then indexes the sorted copy with `tabix` counting tab-separated columns. A row with an empty cell before `hm_pos` is filtered and sorted on the wrong columns, so when positions differ in their number of digits the copy is sorted as text, the index build fails and the plugin is dropped.

**Why this is a defect.** The file is tab-separated, and the plugin reads its column numbers from the tab-separated header. An empty cell is a column like any other; only the two text tools, left at their default field splitting, read it as no field and shift every later column by one. The sort key the plugin asks for, the position as a number, becomes the position as text, under which `14482906` precedes `9053097`, and the index the plugin then builds on the tab columns meets positions out of order. The file's rows are well formed; the processing is not.

**What vep-rs does instead.** vep-rs reads a summary-statistics file by its tab-separated columns throughout: the `hm_chrom` and `hm_pos` cells by the header's column numbers, an `NA` in either skipping the row, and the rows of a chromosome ordered by position as a number. An empty cell before `hm_pos` is read as the empty cell it is, so the file loads and both records are matched.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-vep` at release/116.2):

- `GWAS.pm`, lines 434-451, `parse_sstate_header`: the header is split on tabs and `hm_chrom` and `hm_pos` are found by name; lines 127-131 turn their indices into the one-based column numbers `$c` and `$p`.
- `GWAS.pm`, line 460: the rows are filtered with `awk '$c != "NA" && $p != "NA"'`, which, with no field separator set, splits a line on runs of blanks, tabs included, and skips leading blanks. A line opening with an empty cell starts with a tab, so `awk`'s field `$p` is the column after `hm_pos`.
- `GWAS.pm`, line 464: the filtered rows are sorted with `sort -k$c,$c -k$p,$pn`, which splits fields at the transition from non-blank to blank and keeps the leading blanks with the following field. The same shift makes `-k$c,$c` the position compared as text and `-k$p,$pn` the next column compared as a number.
- `GWAS.pm`, line 471: `tabix -s $c -b $p -e $p -f` indexes the sorted copy on its tab-separated columns, where column `$p` is the position. Positions sorted as text are out of numeric order whenever they differ in digit count, `hts_idx_push` reports `Unsorted positions on sequence #1`, and the plugin dies at line 472, `Failed to create index for <file>.sorted.gz: 256`.
- `Bio/EnsEMBL/VEP/Runner.pm`, lines 946-951: the constructor's death is caught, the plugin is dropped with `Failed to instantiate plugin GWAS: <error>` in the warnings file, and the run continues without it.

**What correcting it would change.** With the two tools splitting on tabs (`awk -F'\t'`, `sort -t$'\t'`), the filter and the sort would read the columns the header names and the copy would be sorted by position as a number. The index would then build and the plugin would annotate the records as it does when the cell is filled.

**The class it produces.** `gwas_sstate_index_failure`: the eight `GWAS_` fields, and the `## GWAS_` header lines, that vep-rs alone writes under `type=sstate` on a harmonised file whose rows carry an empty cell before `hm_pos` and whose positions differ in digit count. Ensembl VEP's warnings file carries `Failed to instantiate plugin GWAS`. The class is documented only. The catalog's harmonised files write `NA`, never an empty cell, for a missing value, and the committed file has none.

## 1. Two rows with an empty `hm_variant_id` and positions of seven and eight digits

**In short.** A harmonised file under the committed file's name, `17463246-GCST000028-EFO_0001360.h.tsv`, holds two rows whose first cell, `hm_variant_id`, is empty and whose positions are `9053097` and `14482906`. Run with `--plugin GWAS,type=sstate,file=<file>` on the two matching records, Ensembl VEP 116.2 fails to build the index of its sorted copy, drops the plugin and writes no `GWAS_` field and no `## GWAS_` header line. vep-rs writes the eight header lines and the eight fields on every row of both records. With the first cell filled, both engines write vep-rs's rows.

### 1.1 What a reader needs to know first

Under `type=sstate` the plugin reads one study's harmonised summary-statistics file, named `<PubMed id>-<study accession>-<trait accession>.h.tsv`, and writes the file-name parts as `GWAS_pmid`, `GWAS_study` and `GWAS_accessions`. Ensembl VEP's plugin first writes a processed copy of the file: the rows whose `hm_chrom` and `hm_pos` are not `NA`, sorted by chromosome and position, compressed and indexed with tabix. A record matches a row at its position whose `hm_effect_allele` equals one of the record's alternate alleles.

### 1.2 The record

The file, ten tab-separated columns, every data row opening with an empty cell:

```
hm_variant_id	hm_rsid	hm_chrom	hm_pos	hm_other_allele	hm_effect_allele	hm_beta	hm_odds_ratio	p_value	ci_upper
	rs552964654	21	9053097	G	C	NA	NA	0.5	NA
	rs2822690	21	14482906	C	T	NA	NA	0.641	NA
```

The input, the two records the rows name:

```
21	9053097	rs552964654	G	C	.	.	.
21	14482906	rs2822690	C	T	.	.	.
```

Both engines run `--plugin GWAS,type=sstate,file=<file>` with an empty `--dir` for the processed copy.

### 1.3 What each engine writes

Ensembl VEP 116.2 writes the records' rows with no `GWAS_` field and no `## GWAS_` header line:

```
rs552964654	21:9053097	C	ENSG00000280436	ENST00000623201	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=25;STRAND=1
rs2822690	21:14482906	T	ENSG00000155307	ENST00000285670	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2322;STRAND=-1
```

Its standard error and warnings file read:

```
[E::hts_idx_push] Unsorted positions on sequence #1: 14482906 followed by 9053097
tbx_index_build failed: 17463246-GCST000028-EFO_0001360.h.tsv.sorted.gz
WARNING: Failed to instantiate plugin GWAS: Failed to create index for 17463246-GCST000028-EFO_0001360.h.tsv.sorted.gz: 256
```

vep-rs writes the eight header lines and, on the same rows:

```
rs552964654	21:9053097	C	ENSG00000280436	ENST00000623201	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=25;STRAND=1;GWAS_accessions=EFO_0001360;GWAS_associated_gene=;GWAS_beta_coef=NA;GWAS_odds_ratio=NA;GWAS_p_value=0.5;GWAS_pmid=17463246;GWAS_risk_allele=C;GWAS_study=GCST000028
rs2822690	21:14482906	T	ENSG00000155307	ENST00000285670	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2322;STRAND=-1;GWAS_accessions=EFO_0001360;GWAS_associated_gene=;GWAS_beta_coef=NA;GWAS_odds_ratio=NA;GWAS_p_value=0.641;GWAS_pmid=17463246;GWAS_risk_allele=T;GWAS_study=GCST000028
```

The tab output carries the eight columns in vep-rs alone, and the VCF output the eight `CSQ` fields and `##GWAS_` header lines in vep-rs alone. In the JSON output each of vep-rs's consequence entries carries the eight `gwas_` keys, `"gwas_p_value": 0.5` and `"gwas_study": "GCST000028"` among them, that Ensembl VEP's entries lack. With `21_9053097_G_C` and `21_14482906_C_T` in the first cell, both engines write vep-rs's rows and header lines, and neither warns.

### 1.4 Why Ensembl VEP writes that

The header names `hm_chrom` as column 3 and `hm_pos` as column 4. Each data row opens with a tab, so `awk` and `sort`, splitting on blanks, number `hm_rsid` as field 1, `hm_chrom` as 2 and `hm_pos` as 3. The filter tests fields 3 and 4, the position and `hm_other_allele`, neither `NA`, and keeps both rows. The sort orders by field 3 as text, the position, under which `14482906` precedes `9053097`; its numeric second key is the other allele. The sorted copy is:

```
	rs2822690	21	14482906	C	T	NA	NA	0.641	NA
	rs552964654	21	9053097	G	C	NA	NA	0.5	NA
```

`tabix` reads tab columns 3 and 4, chromosome `21` and the positions in that order, and refuses the index. The plugin dies in its constructor, and the runner drops it.

### 1.5 Why vep-rs's row is the correct one

The file states, in its tab-separated columns, two associations at `21:9053097` and `21:14482906` with effect alleles `C` and `T`, and the two records carry those alleles at those positions. An empty cell in an unrelated column changes none of that, and Ensembl VEP's own plugin writes the same fields when the cell is filled. The rows vep-rs writes are the ones the plugin's rule produces from the file as written.

### 1.6 Try it yourself

```
mkdir -p vepdir
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --dir vepdir \
  --plugin GWAS,type=sstate,file=17463246-GCST000028-EFO_0001360.h.tsv
```

vep-rs writes the eight `GWAS_` fields on every row of both records. Ensembl VEP 116.2 writes the index-build error and no field. Filling the first cell of both rows, or writing `NA` there as the catalog does, makes the engines agree.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/gwas_sstate_index_failure.toml` holds the file, the two records and both engines' rows. The class covers the eight `GWAS_` fields and the `## GWAS_` header lines that vep-rs alone writes under `type=sstate` on a file whose rows carry an empty cell before `hm_pos` and whose positions differ in digit count. No committed dataset reaches the shape, so the class is documented only.
