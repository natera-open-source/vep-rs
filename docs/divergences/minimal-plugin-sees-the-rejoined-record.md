# The record a plugin reads for a `--minimal` split allele

**The defect in one sentence.** Under `--minimal`, Ensembl VEP's output pass hands every plugin the rejoined record, its span and allele string restored, so a plugin that reads the record's span or alleles values every split allele by the record and not by the allele's own trimmed form. vep-rs hands the plugin each split allele's own trimmed span and alleles.

**Why this is a defect.** `--minimal` exists to annotate each allele of a multi-allelic record on its own trimmed span, and Ensembl VEP computes every split allele's consequences that way. It then rejoins the record for output and runs the plugins over the rejoined feature, whose span and alleles are the record's. A plugin that asks "is this a one-base substitution" or "which alternate alleles does this record carry" is answered for the record, not for the allele whose row it is writing. AVI scores no allele of a split record, its span being longer than one base, and OpenTargets writes an association one allele matched on every allele's rows. Written alone, the same allele gets the opposite answer.

**What vep-rs does instead.** vep-rs hands each split allele its own trimmed span and alleles in every output format, so a split record's rows carry what each allele would carry as a record of its own.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `VEP_plugins` at release/116):

- `Bio/EnsEMBL/VEP/InputBuffer.pm`, lines 638-718, `split_variants`: a multi-allelic record one of whose alleles trims becomes one feature per allele on its trimmed span, the first marked `first`, the others `merge_with`.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 2281-2363, `rejoin_variants_in_InputBuffer`: lines 2298 and 2312 build each split feature's rows once, a pass whose plugin results are discarded. Lines 2300-2302 restore the first feature's allele string, start and end to the record's. Line 2323 points every later allele's transcript, regulatory and motif overlaps at that feature, and lines 2332-2347 its intergenic overlap alleles.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 2235, `run_plugins`: the output pass hands each plugin the overlap allele, whose variation feature is now the rejoined one.
- `AVI.pm`, lines 129-134, `run`: the plugin reads `start`, `end` and `ref_allele_string` from that feature and returns nothing unless the span is one base.
- `OpenTargets.pm`, lines 348-381, `run`: the plugin reads the feature's `alt_alleles` (line 349), its window (line 350) and its `ref_allele_string` and `start` (lines 378-380), and matches an export row against them.
- The same read sits in the other ported plugins. In `AlphaMissense.pm` it is lines 245-268, `BayesDel.pm` lines 127-144, `CADD.pm` lines 203-213, `ClinPred.pm` lines 108-120, `dbNSFP.pm` lines 294-336, `dbscSNV.pm` lines 146-166, `DosageSensitivity.pm` lines 156-163 and `LOEUF.pm` lines 171-187. In `MaxEntScan.pm` it is lines 259-265, `PrimateAI.pm` lines 158-178, `REVEL.pm` lines 164-181, `SpliceAI.pm` lines 231-268, `SpliceVault.pm` lines 155-170, `TranscriptAnnotator.pm` lines 212-229 and `UTRAnnotator.pm` lines 147-153.
- `gnomADc.pm`, lines 159-163, `gnomADMt.pm`, lines 263-267, `Conservation.pm`, lines 206-210, `TSSDistance.pm`, lines 86-89, `NearestExonJB.pm`, lines 225-226, and `GWAS.pm`, lines 176-186, read the feature's span the same way, and `GWAS.pm`, line 185, its `alt_alleles`.

**What correcting it would change.** The output pass would hand each plugin the split allele's own feature, so every plugin's field on a split allele's rows would equal the field the allele carries as a record of its own, in Ensembl VEP as in vep-rs.

**The class it produces.** `minimal_plugin_sees_the_rejoined_record`: a plugin field (the VCF `CSQ` field of that name, the JSON key or object) that differs between the engines on a row of a `MINIMISED` multi-allelic record under a plugin whose `run` reads the record's span or alleles. AVI and OpenTargets are measured below, and the class covers every plugin the source list names; gnomADc, gnomADMt, Conservation, NearestExonJB, TSSDistance and GWAS read the same feature, and vep-rs hands them the allele's own span and alleles too. One committed flag set runs a covered plugin under `--minimal`, gnomADc over the coverage dataset (`tests/golden/116/GRCh38-plugin-gnomadc/`), and the golden harness sets the coverage fields of its two split records aside under this class.

Two sibling classes describe the same rejoin through other reads: [`minimal_plugin_value_follows_first_allele`](../plugins/ancestral-allele.md#1-the-ancestral-allele-printed-for-a-later-allele-of-a---minimal-split-record), a value cached on the first allele's feature, and [`nmd_minimal_rejoined_end`](nmd-minimal-rejoined-end.md), the rejoined end a rule tests.

**How often the shape occurs.** Tier 2: 1,422 of the cuts' records (91 per million) are multi-allelic records one of whose alleles trims beyond the anchor, so that `--minimal` splits and rejoins them; the plugin's table or rule decides whether the record's reading and the allele's differ. 22 corpora carry the shape and one committed flag set runs a covered plugin under `--minimal`, the coverage dataset's gnomADc run, which reaches it on two records.

## 1. AVI on a split record of two one-base substitutions

**In short.** The record `21 5030191 avi_aa_ag_at AA AG,AT`, written over the AVI dataset's synthetic table (`tests/golden/116/GRCh38-plugin-avi-fixture/`), has two alleles that each trim to a one-base substitution at 21:5030192 under `--minimal`, `A>G` and `A>T`. The table scores both. vep-rs writes each allele's row with its score and link; Ensembl VEP 116.3 writes both rows without an AVI field, since the rejoined record spans two bases. Without `--minimal` neither engine scores the record, its alleles being two bases long, and a bi-allelic `AA AG` under `--minimal` is scored by both.

### 1.1 What a reader needs to know first

AVI scores a one-base substitution whose reference and alternate are each one of A, C, G and T, from a table row at the same chromosome, position, reference and alternate (`AVI.pm:135-136,139-151`); every other allele gets no field. A two-base record whose alleles share their first base, such as `AA AG,AT`, is not chopped when it is read, since no allele is an insertion or deletion, so without `--minimal` its alleles are `AG` and `AT` over two bases. Under `--minimal` each allele trims to one base at 21:5030192 and the record is split, annotated and rejoined, every row carrying `MINIMISED=1` and its `ALLELE_NUM`. The dataset's synthetic table holds `A>G` (raw 1.8865, PHRED 3.993) and `A>T` (raw 0.7105, PHRED 2.071) at 21:5030192.

### 1.2 The record

```
21	5030191	avi_aa_ag_at	AA	AG,AT	.	.	.
```

Both engines run `--minimal --plugin AVI,file=<table>,raw=1,links=1` against the dataset's cache, the dataset's own plugin options.

### 1.3 What each engine writes

Ensembl VEP 116.3, default format:

```
avi_aa_ag_at	21:5030191-5030192	G	-	-	-	intergenic_variant	-	-	-	-	-	-	ALLELE_NUM=1;IMPACT=MODIFIER;MINIMISED=1
avi_aa_ag_at	21:5030191-5030192	T	-	-	-	intergenic_variant	-	-	-	-	-	-	ALLELE_NUM=2;IMPACT=MODIFIER;MINIMISED=1
```

vep-rs:

```
avi_aa_ag_at	21:5030191-5030192	G	-	-	-	intergenic_variant	-	-	-	-	-	-	ALLELE_NUM=1;IMPACT=MODIFIER;MINIMISED=1;AVI_LINK=https://deepmind.google.com/science/alphagenome/atlas?q=chr21:5030192:A>G&m=variant;AVI_PHRED=3.993;AVI_RAW=1.8865
avi_aa_ag_at	21:5030191-5030192	T	-	-	-	intergenic_variant	-	-	-	-	-	-	ALLELE_NUM=2;IMPACT=MODIFIER;MINIMISED=1;AVI_LINK=https://deepmind.google.com/science/alphagenome/atlas?q=chr21:5030192:A>T&m=variant;AVI_PHRED=2.071;AVI_RAW=0.7105
```

In the tab output Ensembl VEP's `AVI_LINK`, `AVI_PHRED` and `AVI_RAW` cells are `-` on both rows and vep-rs's hold the synthetic table's values. In the VCF output Ensembl VEP's two `CSQ` entries end in three empty fields and vep-rs's in the link, the PHRED and the raw score. In the JSON output Ensembl VEP's two `intergenic_consequences` entries carry no `avi_` key and vep-rs's carry `avi_link`, `avi_phred` and `avi_raw`, the scores as numbers. Without `--minimal` both engines write the two rows with no AVI field, and on the bi-allelic record `21 5030191 AA AG` under `--minimal` both write `AVI_PHRED=3.993` and `AVI_RAW=1.8865` with the link.

### 1.4 Why Ensembl VEP writes that

The input buffer splits the record into a feature per allele, `A/G` and `A/T` at 21:5030192. The rejoin builds each feature's rows, then writes the record's string `AA/AG/AT`, start 5030191 and end 5030192 back onto the first feature and points the second allele's overlap at it. The output pass runs the plugin on each row with that feature: `start` is 5030191, `end` 5030192, and the test at line 134, `$start == $end`, fails for both alleles, so `run` returns an empty hash and no field is written. The bi-allelic record is minimised in place, with no split and no rejoin, so its feature keeps the one-base span and the plugin scores it.

### 1.5 Why vep-rs's row is the correct one

Under `--minimal` the `G` allele of this record is the substitution of `A` by `G` at 21:5030192, the row's `Allele` column says so, and the table's row for that substitution is its score. Both engines score the same allele when it is the record's only alternate. The two-base span the plugin reads after the rejoin is a property of the way the two alleles were written beside each other, not of either allele, and the plugin's own rule, applied to the allele the row names, gives the score. vep-rs applies the rule to the allele.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --dir_plugins <corpus> --minimal \
  --plugin AVI,file=<corpus>/avi_scores.synthetic.tsv.gz,raw=1,links=1
```

From the synthetic table, vep-rs writes `AVI_PHRED=3.993` on the `G` row and `AVI_PHRED=2.071` on the `T` row. Ensembl VEP 116.3 writes neither. Dropping `--minimal`, or writing the record with one alternate allele, makes the engines agree.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/minimal_plugin_sees_the_rejoined_record.toml` holds this record and run with both engines' `AVI_PHRED` values. The class covers a plugin field that differs on a row of a `MINIMISED` multi-allelic record under a plugin whose `run` reads the record's span or alleles. The format comparator has no predicate for the class and applies it from a corpus manifest alone. The golden harness sets aside the gnomADc coverage fields of the two split records in the coverage dataset's `--minimal` run, the one committed flag set that reaches it.

## 2. OpenTargets on a split record with one matching allele

**In short.** The record `21 14640467 ot_gaat_g_ga GAAT G,GA`, written over the OpenTargets dataset's export (`tests/golden/116/GRCh38-plugin-opentargets/`), carries the export's lead GWAS deletion `GAAT>G` as its first allele and a shorter deletion as its second. Under `--minimal` Ensembl VEP 116.3 writes the association's five fields on all six rows of the record; vep-rs writes them on the first allele's three rows alone. Written as a record of its own, the second allele gets no field from either engine, with or without the flag.

### 2.1 What a reader needs to know first

OpenTargets matches an export row by chromosome, position, reference and alternate against every alternate allele of the record, so a row matched by one allele is written on every allele's rows in both engines when the record is not split. The export holds one association at this site, `21 14640467 GAAT G`: study GCST90618078, gene ENSG00000155307, locus-to-gene score 0.5601661801338196, trait OBA_2041707, a lead variant. The record's three alleles share their first base, so both engines read it as `AAT/-/A` at 21:14640468-14640470. Under `--minimal` the `A` allele trims to the deletion of `AT` at 21:14640469-14640470, so the record is split and rejoined.

Under the default `lead_variants_only=1` and column set, five fields are written: `OpenTargets_gwasDiseases`, `OpenTargets_gwasGeneId`, `OpenTargets_gwasLocusToGeneScore`, `OpenTargets_qtlBiosampleName` and `OpenTargets_qtlGeneId`, the two QTL cells `NA`.

### 2.2 The record

```
21	14640467	ot_gaat_g_ga	GAAT	G,GA	.	.	.
```

Both engines run `--minimal --plugin OpenTargets,file=<export>` against the dataset's cache.

### 2.3 What each engine writes

Ensembl VEP 116.3 writes six rows, three transcripts for each allele, every one carrying the five fields; the second allele's row on ENST00000647101:

```
ot_gaat_g_ga	21:14640468-14640470	-	ENSG00000155307	ENST00000647101	Transcript	intron_variant	-	-	-	-	-	-	ALLELE_NUM=2;IMPACT=MODIFIER;STRAND=-1;MINIMISED=1;OpenTargets_gwasDiseases=OBA_2041707;OpenTargets_gwasGeneId=ENSG00000155307;OpenTargets_gwasLocusToGeneScore=0.5601661801338196;OpenTargets_qtlBiosampleName=NA;OpenTargets_qtlGeneId=NA
```

vep-rs writes the same six rows, the five fields on the three `ALLELE_NUM=1` rows alone; its second allele's row on the same transcript:

```
ot_gaat_g_ga	21:14640468-14640470	-	ENSG00000155307	ENST00000647101	Transcript	intron_variant	-	-	-	-	-	-	ALLELE_NUM=2;IMPACT=MODIFIER;STRAND=-1;MINIMISED=1
```

Both alleles trim to a deletion, so `Allele` is `-` on every row and `ALLELE_NUM` tells them apart. In the tab output Ensembl VEP fills the five columns on all six rows and vep-rs on the first allele's three, the second allele's cells `-`. In the VCF output Ensembl VEP's six `CSQ` entries carry the five fields and vep-rs's `ALLELE_NUM=2` entries leave them empty. In the JSON output every one of Ensembl VEP's six `transcript_consequences` entries carries an `opentargets` list and vep-rs's second-allele entries carry none.

Without `--minimal` both engines write the five fields on all six rows. On the bi-allelic records `21 14640467 GAAT G` and `21 14640467 GAAT GA`, with or without `--minimal`, the engines agree: the first carries the fields on every row, the second on none.

### 2.4 Why Ensembl VEP writes that

The input buffer splits the record into the features `AAT/-` at 21:14640468 and `AT/-` at 21:14640469. The rejoin builds each feature's rows, restores the record's string `AAT/-/A`, start 14640468 and end 14640470 on the first feature and points the second allele's transcript overlaps at it. The output pass then runs the plugin on every row with that feature: `alt_alleles` (line 349) is `-` and `A`, the reference `AAT` and the start 14640468, and the export row `AAT/-` at 14640468 matches the first of them. The plugin returns the association for the row whatever allele the row names, so the second allele's rows carry it too.

### 2.5 Why vep-rs's row is the correct one

The second allele, `GAAT>GA`, is the deletion of `AT` at 21:14640469. The export's association belongs to `GAAT>G`, the deletion of `AAT` at 21:14640468, another variant. Ensembl VEP itself writes no field for `GAAT>GA` as a record of its own, with or without `--minimal`, and `--minimal` is the flag that asks for each allele to be judged as its own variant. The association reaches the second allele's rows only through the rejoined record's allele list, which is the first allele's. vep-rs matches the export against the allele the row names, so its rows say of each allele what the export says of it.

### 2.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --dir_plugins <corpus> --minimal \
  --plugin OpenTargets,file=<corpus>/open_targets_vep.slice.tsv.gz
```

vep-rs writes `OpenTargets_gwasGeneId=ENSG00000155307` on the three `ALLELE_NUM=1` rows and no OpenTargets field on the `ALLELE_NUM=2` rows. Ensembl VEP 116.3 writes it on all six. Dropping `--minimal` makes the engines agree on all six rows; writing `GAAT GA` alone makes both engines write no OpenTargets field on any row.

### 2.7 How the comparator treats this class

The reproduction of [section 1.7](#17-how-the-comparator-treats-this-class) holds the AVI shape; this one is the same class read through the plugin's allele list. The class covers an OpenTargets field, or the JSON `opentargets` list, that Ensembl VEP alone writes on a later allele's rows of a `MINIMISED` record whose first allele the export names. No committed flag set runs OpenTargets under `--minimal`.
