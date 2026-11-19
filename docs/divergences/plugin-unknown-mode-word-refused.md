# A mode word outside a plugin's set: LOEUF, pLI and Conservation

**The defect in one sentence.** Ensembl VEP's LOEUF, pLI and Conservation plugins test a mode word against their known words with no branch for any other, so a word outside the set loads the plugin and silently selects nothing, or the other method; vep-rs refuses the plugin.

**Why this is a defect.** A mode word chooses what the plugin computes: LOEUF's `match_by` is `gene` or `transcript`, pLI's second parameter is `gene` or `transcript`, and Conservation's closing `MAX` chooses the maximum over the average. Each module compares the word with string equality and does nothing when no comparison holds. LOEUF and pLI then load, write their header line and annotate no row, so a mistyped word, `match_by=genes` or `transcripts`, yields a run that looks complete and carries no value. Conservation reads its method from its last parameter alone, so a word after `MAX` silently selects the average under a header that says so.

**What vep-rs does instead.** vep-rs refuses each plugin when it loads: `LOEUF: match_by must be gene or transcript, not '<word>'`, `pLI: mode '<word>' is not gene or transcript`, and `Conservation: MAX must be the last parameter`. The refusal is a `Failed to instantiate plugin` warning, in the warnings file beside the output and on standard error unless `--quiet`. The plugin is dropped and the run goes on without its header line or field, as both engines go on when a plugin's constructor dies. The three plugins are described under [LOEUF](../plugins/annotation-store.md#loeuf), [pLI](../plugins/flat-file.md#pli) and [Conservation](../plugins/conservation.md).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-vep` at release/116.2):

- `LOEUF.pm`, lines 146-149 store any `match_by` value; line 152 dies only when the key is absent. In `run`, line 177 tests `eq 'transcript'`, line 184 `eq 'gene'`, and lines 207-209 return `{}` for every other word, so every row gets nothing.
- `pLI.pm`, line 109 takes the second parameter as the mode and line 111 stores it. Lines 138-148 read the scores and set the `pLI_gene_value` header from the file's `gene` column whatever the mode. In `run`, line 218 tests `eq "gene"` and line 224 `eq "transcript"`; with neither true the sub falls off its end and returns the false value of the last test, which is not a hash. `OutputFactory.pm` line 2235 then warns `Plugin 'pLI' did not return a hashref, output ignored!` and writes the row without a value.
- `Conservation.pm`, lines 82-89: `for(@{$self->params})` sets the method to `MAX` when the parameter `eq 'MAX'` and to `AVERAGE` for any other parameter, the file path included, so the last parameter decides. Line 140 names that method in the header and line 241 reads it.

**What correcting it would change.** Each module would die on a word outside its set, as LOEUF already dies when `match_by` is absent, and Ensembl VEP would drop the plugin with a `Failed to instantiate plugin` warning, as vep-rs does.

**The class it produces.** `plugin_unknown_mode_word_refused`: a LOEUF, pLI or Conservation field Ensembl VEP writes, or a header line it prints, under a mode word outside the plugin's set, where vep-rs drops the plugin. No comparison tool sets a row of the shape aside: no committed flag set passes such a word, so a disagreement of this shape would be reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the fields.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a word outside the set in a plugin's parameters, which no ordinary run passes, and no committed flag set carries one. The committed flag sets pass LOEUF `match_by=gene` and `match_by=transcript`, pLI's `transcript` and its default, and Conservation's `MAX` as the last parameter or no method word.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. LOEUF under `match_by=exon`

**In short.** `21 26112000 app_exon G A` on GRCh38 is a synonymous substitution in APP, whose gene has a LOEUF record of `2.5883e-01` in the committed table. Under `match_by=gene` both engines write `LOEUF=2.5883e-01` on all 55 of its rows. Under `match_by=exon`, Ensembl VEP 116.2 writes the `## LOEUF` header line and no value on any row; vep-rs drops the plugin.

### 1.1 What a reader needs to know first

LOEUF reads a tabix-indexed constraint table and writes one field, `LOEUF`, on each transcript row: under `match_by=gene` the smallest numeric value of the gene's records, under `match_by=transcript` the value of the record naming the transcript. The key is required, and its value is one of those two words.

### 1.2 The record

A substitution against the Ensembl 116 GRCh38 cache, with the corpus's constraint table.

```
21	26112000	app_exon	G	A	.	.	.
```

Both engines run `--plugin LOEUF,file=<table>,match_by=exon`; the ordinary run is `match_by=gene`.

### 1.3 What each engine writes

The row on ENST00000346798 in the default format. Under `match_by=gene`, both engines:

```
app_exon	21:26112000	A	ENSG00000142192	ENST00000346798	Transcript	synonymous_variant	354	204	68	G	ggC/ggT	-	IMPACT=LOW;STRAND=-1;LOEUF=2.5883e-01
```

Under `match_by=exon`, Ensembl VEP 116.2 prints the header line `## LOEUF : Loss-of-function observed/expected upper bound fraction` and writes:

```
app_exon	21:26112000	A	ENSG00000142192	ENST00000346798	Transcript	synonymous_variant	354	204	68	G	ggC/ggT	-	IMPACT=LOW;STRAND=-1
```

vep-rs prints no `## LOEUF` header line, writes the same row, and warns:

```
WARNING: Failed to instantiate plugin LOEUF: LOEUF: match_by must be gene or transcript, not 'exon'
```

Over the corpus's 9 records, Ensembl VEP writes 244 rows with a value on 138 of them under `match_by=gene` and on none under `match_by=exon`, with no warning; vep-rs writes the same 244 rows under both words, the field under `gene` alone.

### 1.4 Why Ensembl VEP writes that

Lines 146-149 store `exon` without a test and the plugin loads, so `get_header_info` prints the header line. On each row, `run` reaches line 177, `exon eq 'transcript'` is false, line 184, `exon eq 'gene'` is false, and lines 207-209 return the empty hash, which adds no field.

### 1.5 Why vep-rs's row is the correct one

`exon` is not a mode the plugin has: neither branch of `run` can fire, and the table is never read. A run that loads the plugin, prints its header and writes no value hides the mistyped word behind an output that looks complete and empty. vep-rs refuses the word where it is given, names it, and drops the plugin, the outcome Ensembl VEP itself produces for a `match_by` that is absent.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin LOEUF,file=<table>,match_by=exon
```

vep-rs warns `Failed to instantiate plugin LOEUF: LOEUF: match_by must be gene or transcript, not 'exon'` and writes no `## LOEUF` line. Ensembl VEP 116.2 writes the header line and no `LOEUF` value on any row. Under `match_by=gene` both write `LOEUF=2.5883e-01` on the record's rows.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`scripts/concordance/compare_vep_formats.py`) has no class of this shape, and the golden harness runs the `GRCh38-plugin-loeuf` and `GRCh37-plugin-loeuf` corpora under `match_by=gene` and `match_by=transcript` alone. A disagreement of the shape would be reported as a difference and charged to vep-rs.

## 2. pLI under the mode word `exon`

**In short.** `21 17561425 4453310 G T` on GRCh38 is a missense substitution in a gene whose symbol has the pLI value `0.34` in the committed values file. Under `pLI,<file>,gene` both engines write `pLI_gene_value=0.34` on all 33 of its rows. Under `pLI,<file>,exon`, Ensembl VEP 116.2 writes the `## pLI_gene_value` header line, no value on any row and one warning; vep-rs drops the plugin.

### 2.1 What a reader needs to know first

pLI takes a values file as its first parameter and a mode as its second, `gene` by default: under `gene` each transcript row gets the score of the gene's symbol, under `transcript` the score of its stable id. The header line names the mode's field.

### 2.2 The record

A substitution against the Ensembl 116 GRCh38 cache, with the corpus's gene values file, whose header is `gene pLI`.

```
21	17561425	4453310	G	T	.	.	.
```

Both engines run `--plugin pLI,<file>,exon`; the ordinary run is `--plugin pLI,<file>,gene`.

### 2.3 What each engine writes

The row on ENST00000284878 in the default format. Under `gene`, both engines:

```
4453310	21:17561425	T	ENSG00000154639	ENST00000284878	Transcript	missense_variant	869	782	261	R/L	cGt/cTt	-	IMPACT=MODERATE;STRAND=1;pLI_gene_value=0.34
```

Under `exon`, Ensembl VEP 116.2 prints the header line `## pLI_gene_value : pLI value by gene` and writes:

```
4453310	21:17561425	T	ENSG00000154639	ENST00000284878	Transcript	missense_variant	869	782	261	R/L	cGt/cTt	-	IMPACT=MODERATE;STRAND=1
```

Its warnings file holds one line, `WARNING: 11 : WARNING: Plugin 'pLI' did not return a hashref, output ignored!`. vep-rs prints no `## pLI_gene_value` header line, writes the same row, and warns:

```
WARNING: Failed to instantiate plugin pLI: pLI: mode 'exon' is not gene or transcript
```

Over the corpus's 29 records, Ensembl VEP writes 652 rows with a value on 551 of them under `gene` and on none under `exon`; vep-rs writes the same 652 rows under both words, the field under `gene` alone.

### 2.4 Why Ensembl VEP writes that

Line 109 takes `exon` as the mode and line 111 stores it. Lines 138-148 find the file's `gene` column, set the `pLI_gene_value` header and read every score, so the constructor returns and the plugin loads. On each row, line 218 tests `exon eq "gene"` and line 224 `exon eq "transcript"`; both are false, `run` ends without a `return`, and its value is the false result of the last test. `OutputFactory.pm` line 2235 finds a defined value that is not a hash, warns once, and writes the row without a field.

### 2.5 Why vep-rs's row is the correct one

`exon` is not a mode the plugin has: neither branch of `run` can fire, so no score can reach a row. Ensembl VEP's own output says as much in its warning, once, in the warnings file, while the output carries the header of a mode the user did not ask for. vep-rs refuses the word where it is given, names it, and drops the plugin.

### 2.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin pLI,<file>,exon
```

vep-rs warns `Failed to instantiate plugin pLI: pLI: mode 'exon' is not gene or transcript` and writes no `## pLI_gene_value` line. Ensembl VEP 116.2 writes the header line, no `pLI_gene_value` on any row and the `did not return a hashref` warning. Under `gene` both write `pLI_gene_value=0.34` on the record's rows.

### 2.7 How the comparator treats this class

As in [1.7](#17-how-the-comparator-treats-this-class): nothing is set aside, and the golden harness runs the pLI corpora under `transcript` and the default mode alone.

## 3. Conservation under `MAX` followed by another word

**In short.** `21 18383705 4764123 GTTC G` on GRCh38 deletes three bases whose GERP scores in the committed bigWig average `2.406` and peak at `3.240`. Under `Conservation,<file>,MAX` both engines write `Conservation=3.240` on its 3 rows. Under `Conservation,<file>,MAX,foo`, Ensembl VEP 116.2 writes `Conservation=2.406`, the average, under a header that says average; vep-rs drops the plugin.

### 3.1 What a reader needs to know first

Conservation reads a bigWig of per-base scores and writes one field, `Conservation`, on every row of a record: the value of a single covered base, or the average of several, or under `MAX` their maximum. `MAX` is a positional word after the file. The header line names the method.

### 3.2 The record

A three-base deletion against the Ensembl 116 GRCh38 cache, with the corpus's bigWig.

```
21	18383705	4764123	GTTC	G	.	.	.
```

Both engines run `--plugin Conservation,<file>,MAX,foo`; the ordinary run is `--plugin Conservation,<file>,MAX`.

### 3.3 What each engine writes

The row on ENST00000284885 in the default format. Under `MAX`, both engines print `## Conservation : The maximum conservation score for this site` and write:

```
4764123	21:18383706-18383708	-	ENSG00000154646	ENST00000284885	Transcript	inframe_deletion	578-580	415-417	139	E/-	GAA/-	-	IMPACT=MODERATE;STRAND=-1;Conservation=3.240
```

Under `MAX,foo`, Ensembl VEP 116.2 prints `## Conservation : The average conservation score for this site` and writes:

```
4764123	21:18383706-18383708	-	ENSG00000154646	ENST00000284885	Transcript	inframe_deletion	578-580	415-417	139	E/-	GAA/-	-	IMPACT=MODERATE;STRAND=-1;Conservation=2.406
```

vep-rs prints no `## Conservation` header line, writes the row without the field, and warns:

```
WARNING: Failed to instantiate plugin Conservation: Conservation: MAX must be the last parameter
```

Over the corpus's 87 records, Ensembl VEP's rows under `MAX,foo` are row for row its rows under no method word: 3,043 of 3,171 rows carry a value, and 1,626 of those values differ from the run under `MAX`. vep-rs writes the 3,171 rows with no field.

### 3.4 Why Ensembl VEP writes that

Lines 82-89 walk the parameters in order: the file path is not `MAX`, so the method is set to `AVERAGE`; `MAX` sets it to `MAX`; `foo` is not `MAX`, so the method is set back to `AVERAGE`. Line 140 prints the average header and line 241 averages the span's scores. Nothing reads `foo` again.

### 3.5 Why vep-rs's row is the correct one

The run asked for the maximum. Under Ensembl VEP the word after `MAX` is never read as a parameter, yet it decides the method, so the output carries the other method's values and says so only in a header line. A `MAX` that is not the last parameter is a parameter list the plugin cannot read as given; vep-rs refuses it where it is given and drops the plugin. A word the plugin does not know in a place that does not follow `MAX` is read by neither engine: `Conservation,<file>,max` selects the average in both, `Conservation,<file>,foo,MAX` the maximum in both, and vep-rs reports the word it did not read.

### 3.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin Conservation,<file>,MAX,foo
```

vep-rs warns `Failed to instantiate plugin Conservation: Conservation: MAX must be the last parameter` and writes no `## Conservation` line. Ensembl VEP 116.2 writes the average header line and `Conservation=2.406` on the record's rows. Under `MAX` alone both write `Conservation=3.240`.

### 3.7 How the comparator treats this class

As in [1.7](#17-how-the-comparator-treats-this-class): nothing is set aside, and the golden harness runs the `GRCh38-plugin-conservation` corpus under `MAX` as the last parameter, under `file=`, and under no method word.
