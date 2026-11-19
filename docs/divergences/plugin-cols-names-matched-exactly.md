# A `cols` name read as a pattern: AlphaMissense, OpenTargets and TranscriptAnnotator

**The defect in one sentence.** Ensembl VEP's AlphaMissense, OpenTargets and TranscriptAnnotator plugins test a `cols` name against the file's column names as an anchored regular expression, so a name holding a metacharacter heads a column of undefined values; vep-rs matches names exactly and refuses one the header lacks.

**Why this is a defect.** `cols` names the columns of the data file a plugin writes, `cols=gwasGeneId:gwasDiseases` or `cols=pathogenicity`. Each module checks the list with `grep(/^$col$/, ...)` over the column names: the name becomes a regular expression anchored at both ends and is tested against every column name. A name spelt exactly passes. A name holding a metacharacter, `gwasGeneI.` for `gwasGeneId` or `SIFT_.core` for `SIFT_score`, passes too, because the pattern matches a real column.

The module then keeps the name as typed, so the header line, the value lookup and the output key all use the pattern's text, which names no column. Every value is undefined: the plugin prints a header line with an empty description and an empty cell, `NA` or `-` on every row, and Perl warns `Use of uninitialized value` once per cell. A name whose text is not a valid pattern, `gwasGeneId(`, dies inside the test with a regex error instead of the column message.

**What vep-rs does instead.** vep-rs splits `cols` on `:` and compares each name with the column names as text, so a name passes only when a column is spelt exactly so. A name the header lacks is refused as the plugin loads, with the valid names listed: `OpenTargets: columns not found in the header of <file>: gwasGeneI.`, or `AlphaMissense: the following columns were not found in file header: am_pathogenicit.; valid columns are: ...`. TranscriptAnnotator refuses with `no valid columns selected. Valid column names: ...` when no name passes. The refusal is a `Failed to instantiate plugin` warning; the plugin is dropped and the run goes on without its header line or field, as both engines do when a constructor dies.

The plugins are described under [AlphaMissense](../plugins/annotation-store.md#alphamissense), [OpenTargets](../plugins/annotation-store.md#opentargets) and [TranscriptAnnotator](../plugins/annotation-store.md#transcriptannotator).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `OpenTargets.pm`, lines 219-224: `_parse_colnames` splits `cols` on `:`, stores the names as typed and pushes a name to `@invalid_cols` unless `grep(/^$col$/, @{ $self->{colnames} })` finds a column. Lines 227-228 die when the list is non-empty. Line 308 filters the header by the stored names and line 318 reads each record's hash by them, so a pattern's text finds no description and no value; line 339 writes `NA` for the undefined value.
- `AlphaMissense.pm`, lines 149-163: the same test on the `am_`-prefixed names, the refusal listing the valid columns. Line 216 filters the header and line 279 the record by the stored names.
- `TranscriptAnnotator.pm`, lines 134-150: `grep /^$col$/` over the column names keeps a name as typed and warns about the others; the plugin dies only when no name passes. Line 199 filters the header by the stored names.
- `MaveDB.pm`, line 123, and `PromoterAI.pm`, line 117, validate `cols` with the same test; `EVE.pm`, line 136, tests `class_number` and `GO.pm`, line 131, tests `match` the same way. None of the four is ported.
- The idiom runs the other way in `GWAS.pm`, line 446, where each summary-statistics header name is the pattern and the required names are tested against it, and in `TranscriptAnnotator.pm`, line 115, where fixed names are tested against the header as unanchored, case-insensitive patterns. There the text comes from the file, not the command line; vep-rs reads those header names exactly too, and a header spelt with a metacharacter is a malformed table outside this class.
- The remaining ports that take a column list compare text: `dbNSFP.pm`, line 239, with `eq`; `dbscSNV.pm`, lines 120-123 and 137-140, and `PrimateAI.pm`, lines 136-138, read fixed names; `Geno2MP.pm`, lines 118-119, keeps the names in a hash without any test.

**What correcting it would change.** Each module would compare a `cols` name with the column names as text, `grep { $_ eq $col }` as `dbNSFP.pm` does, and die on a name the header lacks, so `gwasGeneI.` would be refused where it is given, as vep-rs refuses it.

**The class it produces.** `plugin_cols_names_matched_exactly`: a header line, a field or a JSON key Ensembl VEP writes under a `cols` name holding a regular-expression metacharacter, where vep-rs drops the plugin. No comparison tool sets a row of the shape aside: no committed flag set passes such a name, so a disagreement of this shape would be reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the fields.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a metacharacter in a `cols` name, which no ordinary run passes, and no committed flag set carries one. The 15 committed flag sets that run the three plugins pass `cols=all`, column names spelt exactly, or no `cols`.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. OpenTargets under `cols=gwasGeneI.`

**In short.** `21 14824378 ot_gwas_and_eqtl_leads C T` on GRCh38 leads three GWAS credible sets and one eQTL association in the committed export. Under `cols=gwasGeneId` both engines write `OpenTargets_gwasGeneId=ENSG00000180530,ENSG00000180530,ENSG00000180530,NA` on all 17 of its rows. Under `cols=gwasGeneI.`, Ensembl VEP 116.2 writes a header line and a field named `OpenTargets_gwasGeneI.` with four `NA` values on every row; vep-rs drops the plugin.

### 1.1 What a reader needs to know first

OpenTargets reads the Open Targets Platform's tabix-indexed export and writes one field per selected column, `OpenTargets_<column>`, with the matched records' values `,`-joined in file order and `NA` where a record has no value. `cols` selects the columns by name; `all` selects every column.

### 1.2 The record

A substitution against the Ensembl 116 GRCh38 cache, with the corpus's export slice.

```
21	14824378	ot_gwas_and_eqtl_leads	C	T	.	.	.
```

Both engines run `--plugin OpenTargets,file=<export>,cols=gwasGeneI.`; the ordinary run is `cols=gwasGeneId`.

### 1.3 What each engine writes

The row on ENST00000432230 in the default format. Under `cols=gwasGeneId`, both engines print `## OpenTargets_gwasGeneId : column from open_targets_vep.slice.tsv.gz` and write:

```
ot_gwas_and_eqtl_leads	21:14824378	T	ENSG00000281903	ENST00000432230	Transcript	non_coding_transcript_exon_variant	394	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;OpenTargets_gwasGeneId=ENSG00000180530,ENSG00000180530,ENSG00000180530,NA
```

Under `cols=gwasGeneI.`, Ensembl VEP 116.2 prints the header line `## OpenTargets_gwasGeneI. : ` with an empty description and writes, on every one of the 17 rows:

```
ot_gwas_and_eqtl_leads	21:14824378	T	ENSG00000281903	ENST00000432230	Transcript	non_coding_transcript_exon_variant	394	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;OpenTargets_gwasGeneI.=NA,NA,NA,NA
```

Its standard error carries `Use of uninitialized value in concatenation (.) or string at .../OutputFactory/VEP_output.pm line 167.` The tab output has an `OpenTargets_gwasGeneI.` column holding `NA,NA,NA,NA`, the VCF a `CSQ` field holding `NA&NA&NA&NA`, and the JSON an `opentargets` list of four objects, each `{"OpenTargets_gwasGeneI.": null}`. vep-rs prints no `## OpenTargets_` header line, writes the row without the field, and warns:

```
WARNING: Failed to instantiate plugin OpenTargets: OpenTargets: columns not found in the header of <file>: gwasGeneI.
```

Under `cols=gwas.*`, a pattern matching every GWAS column, Ensembl VEP writes one column, `OpenTargets_gwas.*`, with every value `NA`; the pattern selects no column, it only passes the test. Under `cols=gwasGeneId(`, a text that is not a valid pattern, Ensembl VEP's test itself dies, and its warnings file reads:

```
WARNING: Failed to instantiate plugin OpenTargets: Unmatched ( in regex; marked by <-- HERE in m/^gwasGeneId( <-- HERE $/ at .../OpenTargets.pm line 224.
```

vep-rs refuses that name as a column the header lacks; both engines drop the plugin and their rows agree.

### 1.4 Why Ensembl VEP writes that

Line 224 interpolates `gwasGeneI.` into `/^gwasGeneI.$/`, whose `.` matches the `d` of `gwasGeneId`, so `grep` finds one column and nothing is pushed to `@invalid_cols`. Line 220 has stored the name as typed, and line 286 prefixes it to `OpenTargets_gwasGeneI.`. Line 308 reads `$header{OpenTargets_gwasGeneI.}`, a key the header hash lacks, so the description is undefined and the header line ends after the colon. Line 318 reads each matched record's hash by the same key and line 339 writes `NA` for the undefined value, once per matched record.

### 1.5 Why vep-rs's row is the correct one

`gwasGeneI.` is the name of no column: the export has `gwasGeneId`, and a run that passes the test with a pattern reads no value from it. An output that carries the pattern's text as a field name, with every value `NA`, looks like a run whose records lacked the column, and it hides the misspelling behind a header line with no description. vep-rs refuses the name where it is given, lists it, and drops the plugin, the outcome Ensembl VEP itself produces for a name that matches no column at all.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --dir_plugins <corpus> --plugin OpenTargets,file=<export>,cols=gwasGeneI.
```

vep-rs warns `Failed to instantiate plugin OpenTargets: OpenTargets: columns not found in the header of <export>: gwasGeneI.` and writes no `## OpenTargets_` line. Ensembl VEP 116.2 writes `## OpenTargets_gwasGeneI. : ` and `OpenTargets_gwasGeneI.=NA,NA,NA,NA` on every row. Under `cols=gwasGeneId` both write the three gene identifiers and `NA` for the eQTL record.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`scripts/concordance/compare_vep_formats.py`) has no class of this shape, and the golden harness runs the `GRCh38-plugin-opentargets` corpus under the default columns, `cols=all` and `cols=qtlGeneId:qtlBiosampleName:studyId:studyType`, every name spelt exactly. A disagreement of the shape would be reported as a difference and charged to vep-rs.

## 2. AlphaMissense under `cols=am_pathogenicit.`

**In short.** `21 13610379 two_alts_13610379 G A,T` on GRCh38 is a missense substitution of ENST00000299443 by either allele, with predictions `0.0926` and `0.1494` in the committed table. Under `cols=am_pathogenicity` both engines write `am_pathogenicity=0.0926` on the `A` row. Under `cols=am_pathogenicit.`, Ensembl VEP 116.2 writes `am_pathogenicit.=` with no value on every missense row; vep-rs drops the plugin.

### 2.1 What a reader needs to know first

AlphaMissense reads a tabix-indexed prediction table and writes, on each missense transcript row, the selected columns of the matching prediction, each prefixed `am_`. `cols` selects the columns by their output names or the file's own; the default is `am_pathogenicity:am_class`.

### 2.2 The record

A two-allele substitution against the Ensembl 116 GRCh38 cache, with the corpus's prediction table.

```
21	13610379	two_alts_13610379	G	A,T	.	.	.
```

Both engines run `--plugin AlphaMissense,file=<table>,cols=am_pathogenicit.`; the ordinary run is `cols=am_pathogenicity`.

### 2.3 What each engine writes

The `A` row on ENST00000299443 in the default format. Under `cols=am_pathogenicity`, both engines write:

```
two_alts_13610379	21:13610379	A	ENSG00000166351	ENST00000299443	Transcript	missense_variant	603	151	51	D/N	Gac/Aac	-	IMPACT=MODERATE;STRAND=1;am_pathogenicity=0.0926
```

Under `cols=am_pathogenicit.`, Ensembl VEP 116.2 prints the header line `## am_pathogenicit. : ` with an empty description and writes, on each of the four missense rows:

```
two_alts_13610379	21:13610379	A	ENSG00000166351	ENST00000299443	Transcript	missense_variant	603	151	51	D/N	Gac/Aac	-	IMPACT=MODERATE;STRAND=1;am_pathogenicit.=
```

Its warnings file holds `Use of uninitialized value $data in substitution (s///) at .../OutputFactory/VEP_output.pm line 118` and the matching `concatenation` line for line 119, one pair per row. The tab output has an `am_pathogenicit.` column holding `-`, the VCF an empty `CSQ` field, and the JSON an `alphamissense` object `{"am_pathogenicit.": null}`. vep-rs prints no `## am_` header line, writes the row without the field, and warns:

```
WARNING: Failed to instantiate plugin AlphaMissense: AlphaMissense: the following columns were not found in file header: am_pathogenicit.; valid columns are: am_genome, am_uniprot_id, am_transcript_id, am_protein_variant, am_pathogenicity, am_class
```

### 2.4 Why Ensembl VEP writes that

Line 152 prefixes the name to `am_pathogenicit.` and line 158 interpolates it into `/^am_pathogenicit.$/`, whose `.` matches the `y` of `am_pathogenicity`, so the test passes and line 161 does not die. Line 153 has stored the name as typed. Line 216 reads `$header{am_pathogenicit.}`, a key the header hash lacks, so the description is undefined. Line 279 reads the matched prediction's hash by the same key, and the undefined value reaches the writer, which prints the key with nothing after `=` and warns on each row.

### 2.5 Why vep-rs's row is the correct one

`am_pathogenicit.` names no column of the table. The prediction `0.0926` is in the row the plugin matched, under `am_pathogenicity`, and the run asked for a name one character away from it. Writing `am_pathogenicit.=` on every row reports a column that does not exist, with no value, under a header line that describes nothing. vep-rs refuses the name where it is given, lists the six columns the table has, and drops the plugin, as Ensembl VEP's own refusal does for a name that matches no column.

### 2.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --dir_plugins <corpus>/plugin_data \
  --plugin AlphaMissense,file=<table>,cols=am_pathogenicit.
```

vep-rs warns `Failed to instantiate plugin AlphaMissense: AlphaMissense: the following columns were not found in file header: am_pathogenicit.; valid columns are: ...` and writes no `## am_` line. Ensembl VEP 116.2 writes `## am_pathogenicit. : ` and `am_pathogenicit.=` on every missense row. Under `cols=am_pathogenicity` both write `am_pathogenicity=0.0926` on the `A` row and `0.1494` on the `T` row.

### 2.7 How the comparator treats this class

As in [1.7](#17-how-the-comparator-treats-this-class): nothing is set aside, and the golden harness runs the AlphaMissense corpora of both assemblies under the default columns, `cols=all` and `transcript_match=1`.

## 3. TranscriptAnnotator under `cols=SIFT_.core`

**In short.** `21 10521702 rs463608 C T` on GRCh38 has a table row for ENST00000618007 with `SIFT_score` `0.61`. Under `cols=SIFT_score` both engines write `transcript_annotator_SIFT_score=0.61` on that row. Under `cols=SIFT_.core`, Ensembl VEP 116.2 writes `transcript_annotator_SIFT_.core=` with no value there; vep-rs drops the plugin.

### 3.1 What a reader needs to know first

TranscriptAnnotator reads a tabix-indexed table of per-transcript annotations and writes the selected columns on the row of the transcript the table names, each prefixed with the file's name. `cols` selects the columns; every column after the fifth is written by default. A name the header lacks is a warning in Ensembl VEP and the plugin dies only when no name passes.

### 3.2 The record

A substitution against the Ensembl 116 GRCh38 cache, with the corpus's table, whose header is `#Chrom Pos Ref Alt Transcript SIFT_score SIFT_pred Comment`.

```
21	10521702	rs463608	C	T	.	.	.
```

Both engines run `--plugin TranscriptAnnotator,file=<table>,cols=SIFT_.core`; the ordinary run is `cols=SIFT_score`.

### 3.3 What each engine writes

The row on ENST00000618007 in the default format. Under `cols=SIFT_score`, both engines print `## transcript_annotator_SIFT_score : column from <table>` and write:

```
rs463608	21:10521702	T	ENSG00000274391	ENST00000618007	Transcript	splice_region_variant,intron_variant	-	-	-	-	-	-	IMPACT=LOW;STRAND=1;transcript_annotator_SIFT_score=0.61
```

Under `cols=SIFT_.core`, Ensembl VEP 116.2 prints the header line `## transcript_annotator_SIFT_.core : ` with an empty description and writes:

```
rs463608	21:10521702	T	ENSG00000274391	ENST00000618007	Transcript	splice_region_variant,intron_variant	-	-	-	-	-	-	IMPACT=LOW;STRAND=1;transcript_annotator_SIFT_.core=
```

The other eight rows of the record, whose transcripts the table does not name, carry no field in either engine. The tab output has a `transcript_annotator_SIFT_.core` column holding `-`, the VCF an empty `CSQ` field, and the JSON no key for the column on the row. vep-rs prints no `## transcript_annotator_` header line, writes the row without the field, and warns:

```
WARNING: Failed to instantiate plugin TranscriptAnnotator: TranscriptAnnotator: no valid columns selected. Valid column names: SIFT_score, SIFT_pred, Comment
```

### 3.4 Why Ensembl VEP writes that

Line 136 interpolates `SIFT_.core` into `/^SIFT_.core$/`, whose `.` matches the `s` of `SIFT_score`, so the test passes and line 137 pushes the name as typed; no name is invalid, so neither the warning of line 149 nor the death of line 146 fires. Line 180 prefixes it to `transcript_annotator_SIFT_.core`. Line 199 reads the header hash by that key, so the description is undefined, and the row's values are read by the same key, so the field is written with nothing after `=`.

### 3.5 Why vep-rs's row is the correct one

`SIFT_.core` names no column of the table; its score, `0.61`, sits under `SIFT_score`. Ensembl VEP writes a field named by the pattern's text and never a value, and its own check would have warned had the pattern matched nothing. vep-rs matches the name to `SIFT_score`, `SIFT_pred` and `Comment` as text, finds no column, and, with no valid name left, refuses the plugin naming the three columns the table has.

### 3.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin TranscriptAnnotator,file=<table>,cols=SIFT_.core
```

vep-rs warns `Failed to instantiate plugin TranscriptAnnotator: TranscriptAnnotator: no valid columns selected. Valid column names: SIFT_score, SIFT_pred, Comment` and writes no `## transcript_annotator_` line. Ensembl VEP 116.2 writes `## transcript_annotator_SIFT_.core : ` and `transcript_annotator_SIFT_.core=` on the ENST00000618007 row. Under `cols=SIFT_score` both write `transcript_annotator_SIFT_score=0.61` there.

### 3.7 How the comparator treats this class

As in [1.7](#17-how-the-comparator-treats-this-class): nothing is set aside, and the golden harness runs the TranscriptAnnotator corpora under `cols=SIFT_score:Comment`, `prefix=0`, `trim=0` and the default columns.
