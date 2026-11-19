# An `in` list mixing a number and other text

**The defect in one sentence.** Ensembl VEP 116.2 decides whether an `in` list compares a `text:number` or `text(number)` value through its number by one entry it picks from the list's hash. So a list that mixes a number and other text, `gene_hgnc_id in HGNC:2559,1132`, keeps one gene's rows under one Perl hash seed and another gene's under another. vep-rs decides the reading by the list's first written entry and keeps the same rows in every run.

**Why this is a defect.** The `in` operator's list is kept as a hash of its entries. Before the comparison, the value read from the transcript is reduced to its number when the filter's value is a number, and for a list that filter value is the first key the hash yields. A hash yields its keys in an order that depends on the process's hash seed, which Perl draws at random unless `PERL_HASH_SEED` fixes it.

When the picked entry is `1132`, the transcript's `HGNC:1132` is read as `1132`, which is in the list, and `HGNC:2559` as `2559`, which is not. When the picked entry is `HGNC:2559`, both values are compared whole: `HGNC:2559` is in the list and `HGNC:1132` is not. The same command over the same cache therefore keeps the 34 BTG3 rows of the record below in one run and the 2 CXADR rows in another. Output that changes with a hash seed is not a reading of the command line.

**What vep-rs does instead.** vep-rs decides the reading by the list's first written entry, the same in every run: `gene_hgnc_id in 1132,HGNC:2559` reads every value through its number and keeps the BTG3 rows, and `gene_hgnc_id in HGNC:2559,1132` compares every value whole and keeps the CXADR rows. Each of those is one of the two readings Ensembl VEP gives, chosen by the written order rather than by the seed. A list whose entries are all numbers, or none of them, reads the same on both engines in every run, and so does a list read from a file. The same rule governs `vep filter`, whose `filter_vep` counterpart shares the grammar.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/FilterSet.pm`):

- Lines 391-392, `finish_filter_node`: a value with text on both sides of a comma is split on commas into a hash, one key per entry.
- Lines 589-594, `get_input`: a defined value is matched by `([\w\.\-]+)?\:?\(?([\-\d\.e]*)\)?`, and the numeric reading runs when the number group is not empty; `CELL_TYPE` and the HGVS fields are exempt.
- Lines 596-599: when the filter's value is a hash, `(keys %{$value})[0]` replaces it by the first key the hash yields, under the comment `Assumes that all the values are of similar format`.
- Lines 600-602: when that key is true to Perl and matches `^[\-\d\.e]+$`, the value becomes its number, or its text when the text is itself a number.
- Lines 603-604: otherwise the value is returned whole.
- Line 718, `filter_in`: the comparison is `exists` on the hash.

**What correcting it would change.** With the picked entry taken from the list as written, or the numeric test made per entry, the rows kept would be the same in every run of the same command. A list whose entries are all numbers or all other text is unaffected either way.

**The class it produces.** `transcript_filter_in_list_numeric_by_hash_order`: the rows a `--transcript_filter` keeps, and the lines or consequence entries a `vep filter` run keeps, under an `in` list some of whose entries are numbers and some other text, on a value shaped `text:number` or `text(number)`. Ensembl VEP keeps the rows of one reading or the other by hash seed; vep-rs keeps the rows of the reading the first written entry selects. The class is documented only: no committed flag set and no committed tool run passes such a list.

**How often the shape occurs.** Tier 3: the shape is a value of `--transcript_filter`, or of `vep filter`'s `--filter`, that no ordinary command line carries. Of the 392 committed flag sets, 9 pass `--transcript_filter` and none passes an `in` list mixing a number and other text; of the 28 committed `vep filter` runs, one passes an `in` list, and its three entries are all consequence terms.

## 1. Two genes named one by HGNC identifier and one by its number

**In short.** The record `tf_17608184`, `21 17608184 A>G`, lies within reach of 41 transcripts: 34 of BTG3 (`HGNC:1132`), 5 of BTG3-AS1 (`HGNC:53145`) and 2 of CXADR (`HGNC:2559`). Under `--transcript_filter 'gene_hgnc_id in HGNC:2559,1132'`, vep-rs writes the 2 CXADR rows. Ensembl VEP 116.2 writes the 34 BTG3 rows under `PERL_HASH_SEED=0` and the 2 CXADR rows under `PERL_HASH_SEED=11`, and the reversed list `1132,HGNC:2559` gives it the same two answers under the same two seeds, where vep-rs writes the 34 BTG3 rows.

### 1.1 What a reader needs to know first

`--transcript_filter` keeps the transcripts a filter passes before any consequence is computed; a record whose every transcript is dropped is written as `intergenic_variant`. The transcript's `_gene_hgnc_id` key holds the gene's HGNC identifier as `HGNC:<number>`, and `gene_hgnc_id` names that key through the leading-`_` rule. The filter grammar reads a value shaped `text:number` through its number when the filter compares it with a number, so `gene_hgnc_id is 1132` keeps the BTG3 transcripts on both engines. The `in` operator takes a comma-separated list or a file with one entry per line.

### 1.2 The record

The input is the `tests/golden/116/GRCh38-transcript-filter` corpus's own record:

```
21	17608184	tf_17608184	A	G	.	.	.
```

Both engines run `--symbol` over the GRCh38 release 116 cache, with the filter quoted below added.

### 1.3 What each engine writes

Ensembl VEP 116.2 under `PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0`, with `--transcript_filter 'gene_hgnc_id in HGNC:2559,1132'`, writes 34 rows, every one on a BTG3 transcript; the first two:

```
tf_17608184	21:17608184	G	ENSG00000154640	ENST00000339775	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;SYMBOL=BTG3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:1132
tf_17608184	21:17608184	G	ENSG00000154640	ENST00000348354	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;SYMBOL=BTG3;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:1132
```

Ensembl VEP 116.2 under `PERL_HASH_SEED=11`, with the same filter, writes 2 rows:

```
tf_17608184	21:17608184	G	ENSG00000154639	ENST00001092433	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=426;STRAND=1;SYMBOL=CXADR;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:2559
tf_17608184	21:17608184	G	ENSG00000154639	ENST00001128578	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SYMBOL=CXADR;SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:2559
```

Of the seeds 0 to 20, sixteen give the 34 BTG3 rows and five (11, 13, 15, 17 and 19) the 2 CXADR rows. Of 24 runs without a seed, over the four output formats and both written orders, 13 gave 34 rows and 11 gave 2. The reversed list, `gene_hgnc_id in 1132,HGNC:2559`, gives the same answer as the written one under every one of those seeds. The default, tab, VCF and JSON outputs agree with each other in every run: 34 or 2 transcript rows, consequence entries or `transcript_consequences` members.

vep-rs with `gene_hgnc_id in HGNC:2559,1132` writes the 2 CXADR rows above, byte for byte, in every format; with `gene_hgnc_id in 1132,HGNC:2559` it writes the 34 BTG3 rows of the first block, byte for byte.

The same list through the filter tools, over the run's unfiltered 41-row output with `HGNC_ID in HGNC:2559,1132`: `filter_vep` keeps 34 rows under `PERL_HASH_SEED=0` and 2 under `PERL_HASH_SEED=11`, with either written order; `vep filter` keeps 2 rows for that order and 34 for `HGNC_ID in 1132,HGNC:2559`. Over the VCF output under `--only_matched`, the consequence entries kept are the same counts.

### 1.4 Why Ensembl VEP writes that

The list becomes a two-key hash, `HGNC:2559` and `1132`. For each transcript, `get_input` matches the value `HGNC:1132` or `HGNC:2559` and finds a number after the colon, so it picks one key of the hash to stand for the filter's value. When the picked key is `1132`, a number, the value is reduced to its number: `1132` is in the hash and `2559` is not, so the BTG3 transcripts pass. When the picked key is `HGNC:2559`, not a number, the value is compared whole: `HGNC:2559` is in the hash and `HGNC:1132` is not, so the CXADR transcripts pass.

Which key the hash yields first is fixed for a process by its hash seed, so one run keeps 34 rows and another 2, and the written order of the list plays no part. `filter_vep` runs the same `get_input` over the output's `HGNC_ID` field.

### 1.5 Why vep-rs's row is the correct one

The command line names two entries in an order, and the cache holds the same 41 transcripts in every run. A reading of that command over that cache gives one set of rows, the same each time; the two sets Ensembl VEP writes are the two readings of its own rule, and a seed drawn at random chooses between them. vep-rs's rule picks the entry the user wrote first, so the command `gene_hgnc_id in HGNC:2559,1132` means the same thing in every run, and a user who wants the numeric reading writes the number first. Where a list's entries are all of one kind, both engines agree in every run.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --symbol \
  --transcript_filter 'gene_hgnc_id in HGNC:2559,1132'
```

vep-rs writes the 2 CXADR rows; with `1132,HGNC:2559` the 34 BTG3 rows. Ensembl VEP 116.2 with the same filter over the GRCh38 release 116 cache writes 34 rows under `PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0` and 2 under `PERL_HASH_SEED=11 PERL_PERTURB_KEYS=0`, with either written order. Writing the list as `2559,1132` or as `HGNC:2559,HGNC:1132` makes the engines agree under every seed, 36 rows each, since every entry then reads the same way.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/transcript_filter_in_list_numeric_by_hash_order.toml` holds the record, the filter and both engines' row counts, Ensembl VEP's under the seed of the committed references. The class covers the rows a `--transcript_filter` keeps, and the lines or consequence entries a `vep filter` run keeps, under an `in` list mixing a number and other text on a `text:number` or `text(number)` value. No committed flag set or tool run reaches the shape, so no manifest names the class and the class is documented only.
