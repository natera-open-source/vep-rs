# A repeated `QTL` study type and an OpenTargets export without `isLead`

**The defect in one sentence.** Ensembl VEP 116.2's OpenTargets plugin stores the `QTL` token's selection in the variable it then tests as a list, so a `study_types` list naming `QTL` twice, with only `gwas` or `QTL` before the second, dies with a Perl type error and the plugin is dropped. The same plugin filters on the `isLead` and `studyType` cells without checking that the export has those columns, so an export lacking the column under the filter loads, prints the plugin's header lines and annotates nothing. vep-rs reads the repeated `QTL` as one selection and annotates, and refuses such an export as the plugin loads, writing no OpenTargets header line.

**Why this is a defect.** A study type named twice names one set of studies, as `gwas:gwas` does in the same plugin, and the error the plugin dies with is not a refusal written for the user but the interpreter's. A filter on a column the export does not carry can keep no record, so the header lines announce fields no row will carry, and the only trace is a stream of interpreter warnings about an undefined value. In both shapes the engine's answer, the plugin silently absent or silently empty, is not an answer to the request.

**What vep-rs does instead.** vep-rs reads a repeated `QTL` as the one selection it names, every study type ending in `qtl`, and annotates with the QTL columns. It refuses, as the plugin loads, an export whose header lacks `isLead` under `lead_variants_only=1` or lacks `studyType` under a `study_types` filter, with a warning naming the column and the parameter that lifts the filter, and writes no OpenTargets header line.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-vep` at release/116.2):

- `OpenTargets.pm`, lines 147-200, `_get_studytypes`: the tokens of `study_types` are read in order. Line 170 stores the `QTL` token's selection as a compiled pattern in `$qtl_studytypes`, the variable lines 161 and 183 test as a list. A second `QTL` reaches line 161, `scalar(@$qtl_studytypes)`, with the pattern in the variable, and the dereference dies: `Not an ARRAY reference`. A named type after `QTL` reaches line 183, which tests the variable's kind first and dies with the plugin's own conflict message, the refusal both engines share.
- `Bio/EnsEMBL/VEP/Runner.pm`, lines 946-951: a plugin whose constructor dies is dropped with the warning `Failed to instantiate plugin OpenTargets: <error>` and the run goes on without it.
- `OpenTargets.pm`, lines 105-124, `_get_colnames`: the columns are the export header's from the fifth on; nothing checks that `isLead` or `studyType` is among them. Line 354 filters the records at a site by `$_->{result}->{'OpenTargets_isLead'} eq 'true'`, and lines 362-366 by the `OpenTargets_studyType` cell. A column the header lacks reads as an undefined value, which never equals `true` and matches no pattern, so every record is filtered out, with a `Use of uninitialized value` warning per record.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2407-2421, `get_plugin_headers`: the header lines come from the plugin's `get_header_info`, which names the selected columns whether or not a record can reach them.

**What correcting it would change.** A repeated `QTL` would select the QTL study types once and the plugin would load and annotate. An export without the filtered column would be refused as it loads, or the filter lifted, so that no run prints OpenTargets header lines over rows that carry no field.

**The class it produces.** `opentargets_duplicate_study_type_and_missing_islead`: an `OpenTargets_` field, or the JSON `opentargets` object, that vep-rs alone fills under a `study_types` list naming `QTL` twice with only `gwas` or `QTL` before the second. It also covers the `## OpenTargets_` header lines (the `##OpenTargets_` lines and `CSQ` fields of the VCF output, the columns of the tab output) that Ensembl VEP alone prints on an export whose header lacks the filtered column. That is `isLead` under `lead_variants_only=1` or `studyType` under a `study_types` filter, the data rows agreeing. The class is documented only. No committed flag set lists a study type twice, and the one committed export carries both columns.

## 1. `study_types=QTL:QTL` on the committed export

**In short.** Against the OpenTargets dataset's export (`tests/golden/116/GRCh38-plugin-opentargets/`), `--plugin OpenTargets,file=<export>,study_types=QTL:QTL` drops the plugin in Ensembl VEP 116.2 with `Failed to instantiate plugin OpenTargets: Not an ARRAY reference at OpenTargets.pm line 161`, so no row carries an `OpenTargets_` field and no header line names one. vep-rs loads the plugin as under `study_types=QTL`: two header lines, and `OpenTargets_qtlBiosampleName` and `OpenTargets_qtlGeneId` on every row of a record with a lead QTL association.

### 1.1 What a reader needs to know first

`study_types` is a `:`-separated list. `gwas` names the GWAS studies, `QTL` every study type ending in `qtl`, a name such as `eqtl` or `sqtl` one QTL type, and `all` both families. A QTL type listed beside `QTL` is refused by both engines as a conflict, since `QTL` already selects it. The tokens are read in order: `all` ends the reading, and a QTL type or an unknown word is refused where it stands, so a second `QTL` is reached only when every token before it is `gwas` or `QTL`. Under the default `lead_variants_only=1` only associations whose `isLead` is `true` are written. Without `cols`, the columns written are each selected family's defaults: `qtlGeneId` and `qtlBiosampleName` for the QTL family.

### 1.2 The record

The record `ot_sqtl_lead_gwas_nonlead_empty_cells` of the dataset, whose site carries eleven associations in the export: nine GWAS associations and one eQTL association that are not lead variants, and one lead sQTL association:

```
21	15551809	ot_sqtl_lead_gwas_nonlead_empty_cells	G	A	.	.	.
```

```
21	15551809	G	A	21_15551809_G_A	brainseq_leafcutter_brain_21_15526723_15538583_clu_54996_-	sqtl	c6bacb413e74738b23496ed3afaba9d0	4.244	-7	0.504602	true	true	0.05042927658870074	SuSie	true				ENSG00000229425	dorsolateral prefrontal cortex
```

Both engines run `--plugin OpenTargets,file=<export>,study_types=QTL:QTL` against the dataset's cache.

### 1.3 What each engine writes

Ensembl VEP 116.2 writes the record's twenty transcript rows without an `OpenTargets_` field and no `## OpenTargets_` header line; its warnings file reads:

```
WARNING: Failed to instantiate plugin OpenTargets: Not an ARRAY reference at /plugins/OpenTargets.pm line 161.
```

vep-rs writes two header lines and the QTL columns on all twenty rows; the first two rows:

```
## OpenTargets_qtlBiosampleName : column from open_targets_vep.slice.tsv.gz
## OpenTargets_qtlGeneId : column from open_targets_vep.slice.tsv.gz
ot_sqtl_lead_gwas_nonlead_empty_cells	21:15551809	A	ENSG00000229425	ENST00000634644	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;OpenTargets_qtlBiosampleName=dorsolateral prefrontal cortex;OpenTargets_qtlGeneId=ENSG00000229425
ot_sqtl_lead_gwas_nonlead_empty_cells	21:15551809	A	ENSG00000229425	ENST00000634708	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;OpenTargets_qtlBiosampleName=dorsolateral prefrontal cortex;OpenTargets_qtlGeneId=ENSG00000229425
```

Ensembl VEP's rows are the same two lines ending at `STRAND=-1`. The tab output carries the two columns in vep-rs alone, and the VCF output the two `CSQ` fields and `##OpenTargets_` header lines in vep-rs alone. In the JSON output each of vep-rs's transcript entries carries an `opentargets` list, `[{"OpenTargets_qtlBiosampleName": "dorsolateral prefrontal cortex", "OpenTargets_qtlGeneId": "ENSG00000229425"}]`, that Ensembl VEP's entries lack. Under `study_types=QTL` both engines write vep-rs's rows, header lines and object.

### 1.4 Why Ensembl VEP writes that

The first `QTL` token stores the QTL pattern in `$qtl_studytypes` (line 170). The second reaches `if(scalar(@$qtl_studytypes) > 0)` (line 161), the test for a named QTL type listed before it, which dereferences the variable as a list. The variable holds a compiled pattern at that point, the dereference dies with the interpreter's `Not an ARRAY reference`, and the constructor's death reaches the runner's `eval`, which drops the plugin with the warning and continues. The run has no plugin to write header lines or fields.

### 1.5 Why vep-rs's row is the correct one

`QTL` names the set of every study type ending in `qtl`; naming the set twice does not change it, as the same plugin reads `gwas:gwas` as `gwas`. The record's lead association is an sQTL, a member of that set, so the QTL columns belong on its rows, and both engines write them under `study_types=QTL`. The interpreter's type error is not a decision about the request but a fault in the test meant to catch a conflicting type, which a repeated `QTL` is not.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --dir_plugins <corpus> \
  --plugin OpenTargets,file=<corpus>/open_targets_vep.slice.tsv.gz,study_types=QTL:QTL
```

vep-rs writes the two QTL columns on the record's rows. Ensembl VEP 116.2 writes none and the warning above. Replacing `QTL:QTL` by `QTL` makes the engines agree.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/opentargets_duplicate_study_type_and_missing_islead.toml` holds this record and run with both engines' rows. The class covers an `OpenTargets_` field, or the JSON `opentargets` object, that vep-rs alone fills under a `study_types` list naming `QTL` twice with only `gwas` or `QTL` before the second. No committed dataset reaches the shape, so the class is documented only.

## 2. An export without `isLead` under the default filter

**In short.** The committed export with its `isLead` column removed, run under the default `lead_variants_only=1`, loads in Ensembl VEP 116.2. It prints the five `## OpenTargets_` header lines, writes no field on any row and warns `Use of uninitialized value in string eq at OpenTargets.pm line 354` once per association read. vep-rs drops the plugin as it loads, with a warning naming the missing column, and prints no header line. The data rows are the same in both engines. An export without `studyType` under a `study_types` filter behaves the same way, from line 364.

### 2.1 What a reader needs to know first

The export's columns from the fifth on are the plugin's record columns, named by the header. `lead_variants_only=1`, the default, keeps the associations whose `isLead` cell is `true`; `lead_variants_only=0` keeps them all. The header lines of the default output, one `## OpenTargets_<column> : column from <file>` per selected column, announce the fields the rows may carry; the VCF output declares them in its `CSQ` header and the tab output as columns.

### 2.2 The record

The export of [section 1](#1-study_typesqtlqtl-on-the-committed-export) with its sixteenth column, `isLead`, cut from the header and every row, compressed with bgzip and indexed on its first two columns, so that its header reads:

```
#chromosome	position	referenceAllele	alternateAllele	variantId	studyId	studyType	studyLocusId	pValueMantissa	pValueExponent	beta	is95CredibleSet	is99CredibleSet	posteriorProbability	finemappingMethod	gwasGeneId	gwasLocusToGeneScore	gwasDiseases	qtlGeneId	qtlBiosampleName
```

The input is the record of section 1 and the dataset's record `ot_gwas_insertion_lead`:

```
21	15551809	ot_sqtl_lead_gwas_nonlead_empty_cells	G	A	.	.	.
21	14979687	ot_gwas_insertion_lead	G	GT	.	.	.
```

Both engines run `--plugin OpenTargets,file=<export without isLead>` with no other parameter.

### 2.3 What each engine writes

Ensembl VEP 116.2's header block carries the five default columns' lines:

```
## OpenTargets_gwasDiseases : column from no_islead.tsv.gz
## OpenTargets_gwasGeneId : column from no_islead.tsv.gz
## OpenTargets_gwasLocusToGeneScore : Locus-to-gene (L2G) scores to predict causal genes at GWAS loci; column from no_islead.tsv.gz
## OpenTargets_qtlBiosampleName : column from no_islead.tsv.gz
## OpenTargets_qtlGeneId : column from no_islead.tsv.gz
```

and none of its thirty-five rows carries an `OpenTargets_` field; its standard error repeats, once per association read:

```
Use of uninitialized value in string eq at /plugins/OpenTargets.pm line 354, <__ANONIO__> line 4.
```

vep-rs writes no `## OpenTargets_` line, the same thirty-five rows, and in its warnings file:

```
WARNING: Failed to instantiate plugin OpenTargets: OpenTargets: lead_variants_only=1 needs an isLead column, which the header of no_islead.tsv.gz lacks; pass lead_variants_only=0
```

In the tab output Ensembl VEP has the five columns, every cell `-`, and vep-rs has no such column. In the VCF output Ensembl VEP has five `##OpenTargets_` lines and five `CSQ` fields, every one empty, and vep-rs none. The JSON outputs are the same, no record carrying an `opentargets` object. Under `lead_variants_only=0` both engines load the export, print the five header lines and write every association of each site, lead or not, on its record's rows.

### 2.4 Why Ensembl VEP writes that

The plugin takes its columns from the header (lines 105-124) and its header lines from the selected columns (lines 295-311), with no test that `isLead` is among the columns. At each site the filter `grep { $_->{result}->{'OpenTargets_isLead'} eq 'true' } @data` (line 354) reads a cell the records do not have: an undefined value, which the interpreter warns about and compares as the empty string, never `true`. Every association is filtered out before the allele match, so no row receives a field, while the header lines, decided before any record is read, stand. A `study_types` filter on an export without `studyType` fails the same way at lines 362-366, where an undefined cell matches no pattern and equals no name.

### 2.5 Why vep-rs's row is the correct one

The data rows agree; the difference is the header, and a header line announces a field the output can carry. Under a filter on a column the export does not have, no row can carry the field, so the five lines describe nothing in the file, and the one sign that the filter cannot be applied is an interpreter warning per association on standard error. Refusing the plugin as it loads names the column and the parameter that lifts the filter, and leaves the header true to the rows. The filter itself is an unanswerable request, not a column selection.

### 2.6 Try it yourself

```
zcat <corpus>/open_targets_vep.slice.tsv.gz | cut -f 1-15,17-21 | bgzip -c > no_islead.tsv.gz
tabix -s 1 -b 2 -e 2 no_islead.tsv.gz
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --dir_plugins <corpus> --plugin OpenTargets,file=no_islead.tsv.gz
```

Ensembl VEP 116.2 prints the five `## OpenTargets_` header lines and no field. vep-rs prints the warning and no header line. Adding `lead_variants_only=0` to the plugin's parameters makes the engines agree on every line.

### 2.7 How the comparator treats this class

The reproduction of [section 1.7](#17-how-the-comparator-treats-this-class) reproduces the first shape; this one is read from the engines' header blocks. The class covers the `## OpenTargets_` header lines, the VCF `##OpenTargets_` lines and empty `CSQ` fields, and the tab columns that Ensembl VEP alone prints on an export lacking `isLead` under `lead_variants_only=1` or lacking `studyType` under a `study_types` filter. The warnings comparison counts the two engines' warnings under it. The committed export carries both columns, so the class is documented only.
