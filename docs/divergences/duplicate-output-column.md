# A key listed twice by the tab and VCF writers

**The defect in one sentence.** Ensembl VEP's tab and VCF writers append every plugin key and every `--custom` short name to the column list after the flag columns, without testing whether the list already carries the name. So a source named `SOURCE`, a dbNSFP column named `TSL` under `--tsl`, or a plugin named twice heads two columns of the tab output and two fields of the VCF's `CSQ` list, each carrying the one value the row holds. vep-rs lists each such name once, as one column and one field.

**Why this is a defect.** A row holds one value per key: the writers read each listed column from one hash of the row, so the second column repeats the first byte for byte, and the tab header's description line for the key repeats with it. A reader that indexes the table by column name, `filter_vep` and `vep filter` among them, sees one column; a reader that counts columns sees one more than the key list promises.

The flag columns are already made unique among themselves, so two flags that both add `SOURCE` (`--merged` and `--custom`) list it once. The plugin keys were left out of that test, against the flag columns, against each other and against the custom keys, and the custom keys against the flag columns. Under `--fields` the writers print the list as given and append nothing, so no name repeats there.

**What vep-rs does instead.** vep-rs lists each column once. A plugin key or a custom short name that equals a flag column takes the flag's place under its flag, and a plugin key handed twice is one column: one tab column, one description line and one `CSQ` field, carrying the value the row holds. The default and JSON outputs of both engines agree, since they write the row's keys rather than a column list.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/OutputFactory/Tab.pm`, lines 187-202, `fields`: the 13 default columns (line 193), the flag columns (194), then every plugin header key and every custom header key (195), appended with no test against the names already listed. `OutputFactory.pm`, lines 2407-2420, `get_plugin_headers`, hands over every loaded plugin's keys, a plugin loaded twice twice.
- `Bio/EnsEMBL/VEP/OutputFactory/Tab.pm`, lines 125-129, `output_hash_to_line`: one cell per listed column, read from the row's hash by name, so a column listed twice prints its one value twice. Lines 144-155, `description_headers`: one `## <key> : <text>` line per listed column (153), the text taken from the flag's description table first, so the repeated line repeats the flag's text.
- `Bio/EnsEMBL/VEP/OutputFactory/VCF.pm`, lines 434-455, `fields`: the flag columns are filtered against the VCF's own fixed columns alone (445-447), and the plugin and custom keys are appended untested (449). Lines 266-269 print the list as the `CSQ` header's `Format`, and lines 386-387 write one `CSQ` field per listed name.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2516-2537, `flag_fields`: the flag columns are made unique among themselves (2530-2534); the plugin and custom keys never pass through this test.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1004-1014 and 1417-1427, with line 1301: a custom source's record names are pushed into the row's hash under the short name. Line 1583 then sets `SOURCE` on a transcript that names a source cache, a merged cache's, so on such a transcript the flag's value stands and on every other row the names do. Lines 2225-2231: a plugin's keys are written over the hash after the core fields, so the plugin's value stands.

**What correcting it would change.** The class vanishes: the tab output and the `CSQ` list carry the key once in both engines, with the value they already share.

**The class it produces.** `duplicate_output_column_collapsed`: a tab column line, a `## <key> : <text>` description line or a `CSQ` `Format` list in which Ensembl VEP lists a plugin key or custom short name a second time, and the row cells under the second listing; vep-rs lists each once. The class is documented only: no committed flag set names a `--custom` source `SOURCE`, a plugin column under its own flag or a plugin twice, so the golden harness and the format comparator meet no row of it. The concordance figures score consequence sets and never read a column list.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A `--custom` source named `SOURCE`

**In short.** `--custom` adds the `SOURCE` column to every run, the cache a merged cache's transcript came from. A source whose short name is `SOURCE` adds the same name a second time. Ensembl VEP 116.2's tab output has two `SOURCE` columns and two `## SOURCE : Source of transcript` lines, and its VCF lists `SOURCE|SOURCE` at the end of the `CSQ` format, both cells holding the source's record names. vep-rs has one column, one line and one field, with the same names. The default and JSON outputs are identical between the engines.

### 1.1 What a reader needs to know first

Every `--custom` source adds its short name as a key of the row, listing the names of the records that matched ([Custom annotation](../custom-annotation.md#output)). Giving a run any `--custom` source also switches on the `SOURCE` flag column, which names `Ensembl` or `RefSeq` on a merged cache's transcript and is empty elsewhere. The two names meet when the short name is `SOURCE`: one key of the row, holding the source's record names wherever the transcript names no source cache, and two entries in the writers' column lists.

The tab output prints one cell per listed column and one `## <column> : <description>` line per listed column. The VCF output prints one `CSQ` field per listed name and announces the list as the `Format` of its `##INFO=<ID=CSQ,...>` line. The default output prints the row's keys as `key=value` items of its `Extra` column, and the JSON prints the record's `custom_annotations` object, so neither reads a column list.

### 1.2 The record

One deletion, `21 5028997 . TCA T`, run against the `GRCh38-custom` golden corpus (`tests/golden/116/GRCh38-custom/`), whose `custom/regions.bed.gz` holds two records that overlap it, named `exact_1` and `overlap_1`:

```
21	5028997	.	TCA	T	.	.	.
```

Both engines run with the source named `SOURCE`, then, as the control, with the same file named `regions`:

```
--custom file=<corpus>/custom/regions.bed.gz,short_name=SOURCE,format=bed,type=overlap
```

### 1.3 What each engine writes

The tab format's column line and row. Ensembl VEP:

```
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACT	DISTANCE	STRAND	FLAGS	SOURCE	SOURCE
21_5028997_TCA/T	21:5028998-5028999	-	-	-	-	intergenic_variant	-	-	-	-	-	-	MODIFIER	-	-	-	exact_1,overlap_1	exact_1,overlap_1
```

vep-rs:

```
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACT	DISTANCE	STRAND	FLAGS	SOURCE
21_5028997_TCA/T	21:5028998-5028999	-	-	-	-	intergenic_variant	-	-	-	-	-	-	MODIFIER	-	-	-	exact_1,overlap_1
```

Ensembl VEP's tab header carries `## SOURCE : Source of transcript` twice, vep-rs's once. In the VCF, Ensembl VEP's `CSQ` format ends `...|SYMBOL_SOURCE|HGNC_ID|SOURCE|SOURCE` and the record's entry ends `|exact_1&overlap_1|exact_1&overlap_1`; vep-rs's ends `|SOURCE` and `|exact_1&overlap_1`. Removing the second `SOURCE` column from Ensembl VEP's tab output, and the second field from its VCF, gives vep-rs's output byte for byte, the version and time lines of the headers aside.

The default output agrees between the engines: `IMPACT=MODIFIER;SOURCE=exact_1,overlap_1` in `Extra`, under a header that names the flag's column, `## SOURCE : Source of transcript`, and the source's file, `## SOURCE : [PATH]/regions.bed.gz`, one line each. The JSON records are equal. Under the control name `regions`, every format agrees: 19 tab columns and 25 `CSQ` fields on both sides.

### 1.4 Why Ensembl VEP writes that

The tab writer's `fields` list is the 13 default columns, then the flag columns, then every plugin and custom header key, in that order (`Tab.pm` lines 193-195). The flag columns were made unique among themselves (`OutputFactory.pm` lines 2530-2534), so `SOURCE` stands once among them; the custom keys are appended afterwards with no test, so `SOURCE` is listed again.

The row is written one cell per listed column from the row's hash (`Tab.pm` line 128), which holds one `SOURCE` value: the source's record names, pushed under the key (`OutputFactory.pm` lines 1004-1014 and 1301). Line 1583 replaces them with the transcript's source cache only on a merged cache's transcript. The description lines follow the same list (line 153), and the flag's description table supplies `Source of transcript` for both.

The VCF writer builds its list the same way: the flag columns filtered against the VCF's fixed columns alone, then the plugin and custom keys appended (`VCF.pm` lines 445-449). The list is printed as the `CSQ` `Format` (lines 266-269) and read field by field for every entry (lines 386-387).

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/OutputFactory/Tab.pm`, lines 193-195: the column list; line 128: one cell per listed column; line 153: one description line per listed column.
- `Bio/EnsEMBL/VEP/OutputFactory/VCF.pm`, lines 445-449: the `CSQ` field list; lines 266-269 and 386-387: its header and the entries.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2530-2534: the flag columns made unique; lines 1004-1014 and 1301: the custom names pushed into the row's hash; line 1583: `SOURCE` set from a merged cache's transcript afterwards.
- `Bio/EnsEMBL/VEP/Constants.pm`, line 101: the `custom` flag's column, `SOURCE`.

### 1.5 Why vep-rs's row is the correct one

The row holds one `SOURCE` value, and both engines agree on it: the source's record names, or a merged cache's transcript source where the transcript has one, printed in the default output's single `SOURCE=` item and in the JSON once.

A column list is the promise that each name heads one column. A name listed twice over one value adds a cell that no key distinguishes, so a reader indexing the table by name reads one column and a reader counting columns finds a cell the header cannot explain. Ensembl VEP itself keeps the flag columns unique for this reason. vep-rs applies the same rule to the plugin and custom keys and writes the column once, with the value both engines hold.

### 1.6 Try it yourself

Run the record of 1.2 through both engines with the corpus's `custom/regions.bed.gz`. vep-rs's command is below; Ensembl VEP takes the same `--custom` string with the cache options of [Reproducing these examples](reproducing.md):

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --fasta <corpus>/reference.fa.gz --species homo_sapiens --assembly GRCh38 \
  --custom file=<corpus>/custom/regions.bed.gz,short_name=SOURCE,format=bed,type=overlap \
  --tab
```

| Output | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Tab columns | 19, `SOURCE` at 18 and 19 | 18, `SOURCE` at 18 |
| Tab description lines | `## SOURCE : Source of transcript` twice | once |
| `CSQ` format | `...\|HGNC_ID\|SOURCE\|SOURCE` | `...\|HGNC_ID\|SOURCE` |
| Default `Extra` and JSON | `SOURCE=exact_1,overlap_1`; equal records | the same |

Renaming the source `regions` makes the two engines agree in every format, which is the control.

### 1.7 How the comparator treats this class

The class is documented only. No committed flag set names a source `SOURCE`, so the golden harness (`crates/vep-cli/tests/common/mod.rs`) meets no reference output with a repeated column, and the format comparator (`scripts/concordance/compare_vep_formats.py`) has no rule for the class. A comparison case that names such a source compares its tab column lines and `CSQ` formats as text, and the second listing is reported as a difference under this class's name from the case's documentation.

## 2. A dbNSFP column named as an active flag's column

**In short.** The dbNSFP plugin writes the columns a run asks of its file under their own names, and the `GRCh38-plugin-dbnsfp` corpus's synthetic file carries a `TSL` column and an `APPRIS` column, the names of the `--tsl` and `--appris` flag columns. Under `--tsl` with `TSL` asked of the plugin, Ensembl VEP 116.2's tab output has two `TSL` columns and its `CSQ` list two `TSL` fields, each carrying the plugin's value on a row the plugin annotated and the transcript's own support level elsewhere. vep-rs has one `TSL` column with the same values. Without `--tsl`, or with `--tsl` and a plugin column list that leaves `TSL` out, every format agrees.

### 2.1 What a reader needs to know first

The dbNSFP plugin takes the columns a run names after its file, or every column under `ALL`, and writes each under the file's own column name unless `output_prefix=1` prefixes them. `--tsl` adds the `TSL` flag column, the transcript's support level. The two meet when the plugin's list names `TSL`: one key of the row, the plugin's value written over the transcript's, and two entries in the writers' column lists.

### 2.2 The record and the runs

The record `21 29067663 missense_rev_two_changes T C` of the `GRCh38-plugin-dbnsfp` golden corpus (`tests/golden/116/GRCh38-plugin-dbnsfp/`), whose `plugin_data/dbNSFP4.9a_grch38_synthetic.gz` annotates its missense rows. The run under test:

```
--tsl --plugin dbNSFP,<corpus>/plugin_data/dbNSFP4.9a_grch38_synthetic.gz,TSL,SIFT_score
```

The two controls drop `--tsl`, or keep it and ask the plugin for `SIFT_score` alone.

### 2.3 What each engine writes

The tab column line, then one of the record's 38 rows, a non-coding transcript the plugin does not annotate. Ensembl VEP:

```
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACT	DISTANCE	STRAND	FLAGS	TSL	SIFT_score	TSL
missense_rev_two_changes	21:29067663	C	ENSG00000156261	ENST00000470450	Transcript	non_coding_transcript_exon_variant	348	-	-	-	-	-	MODIFIER	-	-1	-	1	-	1
```

vep-rs:

```
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACT	DISTANCE	STRAND	FLAGS	TSL	SIFT_score
missense_rev_two_changes	21:29067663	C	ENSG00000156261	ENST00000470450	Transcript	non_coding_transcript_exon_variant	348	-	-	-	-	-	MODIFIER	-	-1	-	1	-
```

On a row the plugin annotates, both of Ensembl VEP's `TSL` cells hold the plugin's list (`1,NA` on the `ENST00000540844` row, beside `SIFT_score` `0.235,0.723`), and vep-rs's one cell holds the same list. Ensembl VEP's header reads `## TSL : Transcript support level`, `## SIFT_score : SIFT_score from dbNSFP file`, `## TSL : Transcript support level`; vep-rs's the first two. The `CSQ` format ends `...|HGNC_ID|TSL|SIFT_score|TSL` in Ensembl VEP and `...|HGNC_ID|TSL|SIFT_score` in vep-rs, and each entry's cells follow. Removing the second `TSL` from Ensembl VEP's tab and VCF outputs gives vep-rs's, the version and time lines aside; the default and JSON outputs are equal as written.

### 2.4 Why Ensembl VEP writes that

The `--tsl` flag lists `TSL` among the flag columns (`Constants.pm` line 90), and the plugin's header lists `TSL` among its keys, so the writers' column lists hold the name twice (`Tab.pm` line 195, `VCF.pm` line 449). The plugin's value is written into the row's hash over the transcript's own (`OutputFactory.pm` lines 2225-2231), and the one hash value is printed in both cells. Without `--tsl` the name is listed once, as a plugin key; with `--tsl` and no `TSL` plugin column it is listed once, as a flag column.

### 2.5 Why vep-rs's row is the correct one

The row holds one `TSL` value, the plugin's where the plugin wrote one and the transcript's where it did not, and the engines agree on it in every cell. vep-rs lists the column once, under the flag, carrying that value, which is the column a reader asking for `TSL` expects and the one the default and JSON outputs of both engines already print.

### 2.6 Try it yourself

Run the record of 2.2 through both engines with the corpus's dbNSFP file and `--tsl`:

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --tsl --dir_plugins <corpus>/plugin_data \
  --plugin dbNSFP,<corpus>/plugin_data/dbNSFP4.9a_grch38_synthetic.gz,TSL,SIFT_score --tab
```

| Output | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Tab columns | 20, `TSL` at 18 and 20 | 19, `TSL` at 18 |
| Tab description lines | `## TSL : Transcript support level` twice | once |
| `CSQ` format | `...\|HGNC_ID\|TSL\|SIFT_score\|TSL` | `...\|HGNC_ID\|TSL\|SIFT_score` |
| Default and JSON | equal records | the same |

The same file's `APPRIS` column under `--appris` takes the same path. The committed `dbnsfp-all` flag set asks the plugin for every column, `TSL` and `APPRIS` among them, with neither flag on, so its reference lists each name once and the engines agree.

### 2.7 How the comparator treats this class

As in 1.7: documented only. No committed flag set asks a plugin for a column named as a flag it switches on, so the golden harness meets no repeated `TSL` or `APPRIS` column, and a comparison case that does names the class from its documentation.

## 3. A plugin named twice

**In short.** A plugin given twice on the command line is loaded twice, and each instance hands the writers its header keys. With REVEL named twice, Ensembl VEP 116.2's tab output has two `REVEL` columns and two `## REVEL : Rare Exome Variant Ensemble Learner` lines, and its `CSQ` list ends `|REVEL|REVEL`, both cells holding the one score; vep-rs has one column, one line and one field. Two plugins sharing a key take the same path.

### 3.1 What a reader needs to know first

Every loaded plugin declares its header keys, and the writers list them after the flag columns in load order ([Plugins](../plugins.md)). A plugin named twice is two loaded instances with the same keys; a key both instances write lands in the one row hash once, the later instance's value standing.

### 3.2 The record

The record `21 26965413 4637665 T C` of the `GRCh38-plugin-revel` golden corpus (`tests/golden/116/GRCh38-plugin-revel/`), whose `plugin_data/revel.tsv.gz` scores its missense row:

```
--plugin REVEL,file=<corpus>/plugin_data/revel.tsv.gz --plugin REVEL,file=<corpus>/plugin_data/revel.tsv.gz
```

### 3.3 What each engine writes

The tab column line and the scored row. Ensembl VEP:

```
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACT	DISTANCE	STRAND	FLAGS	REVEL	REVEL
4637665	21:26965413	C	ENSG00000154736	ENST00000284987	Transcript	missense_variant	1676	979	327	S/G	Agc/Ggc	-	MODERATE	-	-1	-	0.176	0.176
```

vep-rs:

```
#Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACT	DISTANCE	STRAND	FLAGS	REVEL
4637665	21:26965413	C	ENSG00000154736	ENST00000284987	Transcript	missense_variant	1676	979	327	S/G	Agc/Ggc	-	MODERATE	-	-1	-	0.176
```

Ensembl VEP's header carries `## REVEL : Rare Exome Variant Ensemble Learner` twice, vep-rs's once; the `CSQ` format ends `|HGNC_ID|REVEL|REVEL` against `|HGNC_ID|REVEL`, and each of the six entries carries the cell twice. Removing the second `REVEL` column and field from Ensembl VEP's outputs gives vep-rs's, the version and time lines aside.

### 3.4 Why Ensembl VEP writes that

`get_plugin_headers` (`OutputFactory.pm` lines 2407-2420) walks every loaded plugin and pushes each one's keys, so two instances of REVEL push `REVEL` twice, and the writers append the pushed list without a test (`Tab.pm` line 195, `VCF.pm` line 449). Each instance writes its score into the row's hash under the same key (`OutputFactory.pm` lines 2225-2231), so both cells print the one value.

### 3.5 Why vep-rs's row is the correct one

As in 1.5: one key, one value, one column. The second instance adds nothing the row can hold, and vep-rs lists the key once whether the duplicate comes from a flag column or from an earlier plugin.

### 3.6 Try it yourself

Run the record of 3.2 through both engines with the corpus's REVEL table named twice:

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 --dir_plugins <corpus>/plugin_data \
  --plugin REVEL,file=<corpus>/plugin_data/revel.tsv.gz \
  --plugin REVEL,file=<corpus>/plugin_data/revel.tsv.gz --tab
```

| Output | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Tab columns | 19, `REVEL` at 18 and 19 | 18, `REVEL` at 18 |
| Tab description lines | `## REVEL : ...` twice | once |
| `CSQ` format | `...\|HGNC_ID\|REVEL\|REVEL` | `...\|HGNC_ID\|REVEL` |

### 3.7 How the comparator treats this class

As in 1.7: documented only. No committed flag set names a plugin twice or two plugins sharing a key.
