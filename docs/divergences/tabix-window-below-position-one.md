# No annotation from seven tabix-backed plugins at the start of a chromosome

**The defect in one sentence.** Seven of Ensembl VEP's tabix-backed plugins read their table rows from one or two bases before a record's start. A record at position 1, or at position 2 under a two-base window, asks for a window that starts at or below 0, and the read dies. The record carries none of the plugin's fields and a `Plugin '<Name>' went wrong` warning. vep-rs clamps the query at position 1, so a table row at position 1 or 2 annotates its record.

**Why this is a defect.** The window reaches one or two bases before the record's start. A record at the first positions of a chromosome has fewer bases before it than the window asks for, and the correct window is the bases that exist. The base class every such plugin reads its file through passes the arithmetic result on: a start of 0 is false in Perl and refused as no start at all, and a start of -1 is an unparsable region for the tabix reader. Either way the exception is caught by the plugin runner and reported, and a table row at position 1 or 2 can annotate nothing.

**What vep-rs does instead.** vep-rs clamps every tabix query at position 1 and matches the record against the rows it holds, so a table row at position 1 or 2 annotates its record as a row anywhere else does.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-variation` at release/116, `ensembl-vep` at release/116.2, Bio::DB::HTS 2.11):

- The two-base window, `get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})` in each plugin's `run`: `CADD.pm`, line 227; `AlphaMissense.pm`, line 248; `OpenTargets.pm`, line 350; `SpliceVault.pm`, line 158; `TranscriptAnnotator.pm`, line 218. The one-base window, `$vf->{start} - 1`, which reaches 0 for a record at position 1 alone: `dbNSFP.pm`, line 314; `dbscSNV.pm`, line 163. `CADD.pm`, line 220, takes a structural variant's start one base earlier before the same subtraction, so a structural record at position 1, 2 or 3 under a CADD structural-variant table reaches the shape.
- `Bio/EnsEMBL/Variation/Utils/BaseVepTabixPlugin.pm`, line 354, `get_data`: `die("ERROR: No start specified\n") unless $s`, which a start of 0 fails; lines 596-602, `_get_data_hts`: the region is formatted `chr:start-end`, so a start of -1 gives `21:-1-1`.
- `Bio/DB/HTS/Tabix.pm`, lines 102-104, `query`: the region must match `chr`, `chr:start` or `chr:start-end` with unsigned digits, and `21:-1-1` is refused with `You must specify a region in the format chr, chr:start or chr:start-end`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2243-2245: the exception from the plugin's `run` becomes the warning `Plugin '<Name>' went wrong: <reason>` and the row is written without the plugin's fields.
- `GWAS.pm`, line 176, reads its processed catalog through the same base class with the same two-base window; its class is [No GWAS annotation at position 1 or 2 of a chromosome](gwas-window-below-position-one.md).

**What correcting it would change.** A record whose window reaches below position 1 would be matched against the rows from position 1 to its end, as vep-rs matches it, and the warning would go.

**The class it produces.** `tabix_window_below_position_one`: a field of the seven plugins (the JSON keys of the same names) that vep-rs alone fills on a record whose start is 1 or 2; dbNSFP and dbscSNV reach it at 1 alone, a structural record under CADD at 1 to 3. In `--warnings` mode it also covers the reference's `Plugin '<Name>' went wrong` lines on such records. The class is documented only: no manifest names it and no comparator applies it. The committed records whose start is 1 or 2, a mitochondrial deletion at position 1 of the fly corpus and the deletion and insertion written at position 0 of the VCF edge-case corpus, run under none of the seven plugins.

## 1. Three records at the start of chromosome 21 under CADD

**In short.** With the `GRCh38-plugin-cadd` corpus's synthetic SNV table carrying its three rows at 19268555 moved to positions 1 to 3 with the alternate `G`, the records `21 1 pos1 C G`, `21 2 pos2 C G` and `21 3 pos3 C G` are run under CADD. Ensembl VEP writes `CADD_PHRED=3.312;CADD_RAW=1.239896` on the record at 3 and no CADD field on the other two. vep-rs writes each record's own row: `CADD_PHRED=39.801;CADD_RAW=9.984677` at 1, `CADD_PHRED=22.363;CADD_RAW=7.465190` at 2 and `CADD_PHRED=3.312;CADD_RAW=1.239896` at 3.

### 1.1 What a reader needs to know first

A CADD table row is read under the column names `#Chrom`, `Pos`, `Ref` and `Alt` (`CADD.pm` lines 265-268). Ensembl VEP's plugin reads the rows from two bases before the record's start to its end (`CADD.pm`, line 227) and takes the first row whose alleles match. The three records here are intergenic, at the first three bases of chromosome 21, so each writes one row and the only cell that can differ is the CADD field.

### 1.2 The record

Three substitutions at positions 1, 2 and 3, against the Ensembl 116 GRCh38 cache, with a CADD table carrying one row at each of the three positions.

```
21	1	pos1	C	G	.	.	.
21	2	pos2	C	G	.	.	.
21	3	pos3	C	G	.	.	.
```

The table is the committed `tests/golden/116/GRCh38-plugin-cadd/plugin_data/whole_genome_SNVs_synthetic.tsv.gz`, the corpus's synthetic fixture, with its three rows at 19268555 moved to positions 1, 2 and 3 and their alternate set to `G`, compressed and indexed again. Both engines run `--plugin CADD,snv=<that table>`.

### 1.3 What each engine writes

The rows in the default format. Ensembl VEP 116.2:

```
pos1	21:1	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
pos2	21:2	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
pos3	21:3	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;CADD_PHRED=3.312;CADD_RAW=1.239896
```

Its warnings file carries two lines behind its `WARNING:` prefix, the first for the record at 1 and the second for the record at 2:

```
Plugin 'CADD' went wrong: You must specify a region in the format chr, chr:start or chr:start-end at /usr/local/lib/aarch64-linux-gnu/perl/5.34.0/Bio/DB/HTS/Tabix.pm line 104, <__ANONIO__> line 5.
Plugin 'CADD' went wrong: ERROR: No start specified
```

vep-rs:

```
pos1	21:1	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;CADD_PHRED=39.801;CADD_RAW=9.984677
pos2	21:2	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;CADD_PHRED=22.363;CADD_RAW=7.465190
pos3	21:3	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;CADD_PHRED=3.312;CADD_RAW=1.239896
```

The tab, VCF and JSON formats carry the same values in their own columns and keys, and vep-rs writes no warning about the plugin. Both engines exit 0.

### 1.4 Why Ensembl VEP writes that

For the record at 1 the plugin asks for the window -1 to 1. The start is not 0, so the base class formats the region `21:-1-1` and hands it to the tabix reader, whose region pattern allows only unsigned digits after the colon; the reader dies with its format message. For the record at 2 the window is 0 to 2, and a start of 0 is false in Perl, so the base class dies with `No start specified` before any read. The plugin runner catches each exception, writes it to the warnings file and writes the row without the plugin's fields. The record at 3 asks for 1 to 3, a valid region, and is scored.

### 1.5 Why vep-rs's row is the correct one

The row at position 1 names the record at 1, and the row at position 2 names the record at 2. The window is a search range, not a property of the record, and a range that reaches before the first base of a chromosome means the bases that exist. vep-rs clamps the range at 1 and reads the rows there. Ensembl VEP scores the same row for the same alleles anywhere the window stays at or above 1, as the record at 3 shows.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin CADD,snv=cadd_pos_synthetic.tsv.gz
```

vep-rs writes a `CADD_PHRED` and a `CADD_RAW` value on all three rows. Ensembl VEP 116.2 writes them on the row at 3 alone and two `Plugin 'CADD' went wrong` lines to its warnings file. Moving the three table rows to positions 3, 4 and 5 and the records with them makes the engines agree.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/tabix_window_below_position_one.toml` holds the three records, the table's construction and both engines' rows. The class covers a field of the seven plugins that vep-rs alone fills on a record whose window reaches below position 1. No committed dataset reaches the shape, so no manifest names the class and the class is documented only.
