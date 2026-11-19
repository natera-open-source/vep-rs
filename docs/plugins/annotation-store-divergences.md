# Where the annotation-store plugins differ from Ensembl VEP

Where vep-rs's AlphaMissense, CADD, SpliceVault, TranscriptAnnotator, gnomADMt,
LOEUF and OpenTargets output differs from Ensembl VEP 116.2's: eight classes,
each with a worked example.

In each the plugin annotates a row, or describes a field, in a way the record or
the data file does not support. vep-rs writes what the record supports. Six are
classes the comparison tools set aside. Section 8, OpenTargets' cell of `0`, is
a documented class no comparison tool codes, since no committed row falls under
it.

[How to read the counts](../intended-divergences.md#how-to-read-the-counts) on
[Intended divergences from Ensembl VEP](../intended-divergences.md) defines the
format comparator, the golden harness and the corpus manifest the sections below
name. The concordance figures score consequence sets and never read a plugin
field, so nothing is set aside in them. The engines, caches and commands the
examples use are on [Reproducing these examples](../divergences/reproducing.md);
every row is in the default output format.

## 1. AlphaMissense: a two-digit transcript version under `transcript_match=1`

**In short.** Under `transcript_match=1`, Ensembl VEP's AlphaMissense plugin
strips one digit of a prediction's transcript version. A prediction on a
transcript whose version has two digits therefore never matches, and the row
gets no `am_` field. vep-rs strips the whole version and matches it. On the
GRCh38 chromosome 21 ClinVar records, 3,168 rows of 3,168 records differ, every
one a prediction on a two-digit version, and no value differs where both engines
match.

### 1.1 What a reader needs to know first

`transcript_match=1` writes a prediction only when its transcript is the row's.
The prediction file names transcripts with their version,
`ENST00000285679.10`; the row's `Feature` is the stable id without one. The
versions of the release 116 transcripts run from 1 to over 20.

### 1.2 The record

The ClinVar GRCh38 record `21 15730419 4602719 A G` against the published
GRCh38 prediction file, whose row for the site is:

```
chr21	15730419	A	G	hg38	Q9UHP3	ENST00000285679.10	Q9R	0.196	likely_benign
```

Both engines run `--plugin AlphaMissense,file=<file>,transcript_match=1`.

### 1.3 What each engine writes

Ensembl VEP's row on the transcript carries no `am_` field:

```
4602719	21:15730419	G	ENSG00000155313	ENST00000285679	Transcript	missense_variant	495	26	9	Q/R	cAg/cGg	-	IMPACT=MODERATE;STRAND=1
```

vep-rs's row carries the prediction:

```
4602719	21:15730419	G	ENSG00000155313	ENST00000285679	Transcript	missense_variant	495	26	9	Q/R	cAg/cGg	-	IMPACT=MODERATE;STRAND=1;am_class=likely_benign;am_pathogenicity=0.196
```

### 1.4 Why Ensembl VEP writes that

`AlphaMissense.pm`, line 259, removes the version with `s/\.[0-9]//g`, one dot
and one digit. `ENST00000285679.10` therefore becomes `ENST000002856790`, which
equals no transcript's stable id, and the comparison at line 260 fails. A
one-digit version is removed whole and matches.

### 1.5 Why vep-rs's row is the correct one

The prediction names the row's transcript: `ENST00000285679` at version 10. The
`transcript_match` test asks whether the prediction's transcript is the row's,
and the version is not part of the row's identifier. The plugin's own match on
a one-digit version shows the intent; the two-digit case fails on the pattern
alone.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin AlphaMissense,file=AlphaMissense_hg38.tsv.gz,transcript_match=1
```

Ensembl VEP writes no `am_` field on `ENST00000285679`; vep-rs writes
`am_class=likely_benign;am_pathogenicity=0.196`. Without `transcript_match=1`
both engines write the prediction. The chromosome 21 counts, with the inputs'
digests and both commands, are in
[the comparison's record](../concordance-provenance/alphamissense-transcript-match-clinvar-grch38-chr21.json).

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies
`alphamissense_two_digit_transcript_version` when a comparison case names it in
`masks`. An `am_` field, or the JSON `alphamissense` object, that the candidate
fills where the reference has none is set aside, counted once per row. A
differing value where both engines fill the field, or a field the reference
alone fills, stays a difference. The golden corpora carry one-digit versions
alone, so the class applies to none of their rows.

## 2. CADD: structural alleles outside insertion, deletion and duplication

**In short.** Ensembl VEP's CADD plugin matches a structural variant to the
CADD-SV file by a class word it derives from the record. A class outside
`insertion`, `deletion` and `duplication` gives no word. Every row of such a
record over a span the file covers then throws inside the plugin. Ensembl VEP
writes the exception once per run as a `Plugin 'CADD' went wrong` warning and
the rows without a score.

A deletion written on the reverse strand in Ensembl's own input format has its
class word reverse-complemented and matches no row. Its rows get no score.
vep-rs writes no warning, and scores the reverse-strand deletion as it scores
the forward-strand one.

### 2.1 What a reader needs to know first

The plugin matches a CADD-SV row by class word (`DEL`, `INS`, `DUP`; `CADD.pm`
lines 75-79) and the base before the start (line 220).
A VCF record's class comes from its
symbolic allele; a record in Ensembl's default input format (chromosome, start,
end, allele, strand, identifier) carries a class word in its allele column and a
strand of its own. Ensembl VEP writes each distinct warning text once per run.

### 2.2 The record

A hand-written `<INV>` record over the span the `GRCh38-plugin-cadd-sv` golden
corpus's synthetic structural-variant file scores (`tests/golden/116/GRCh38-plugin-cadd-sv/plugin_data/`), and the same span as a deletion:

```
21	5029000	inv	T	<INV>	.	.	SVTYPE=INV;END=5047000
21	5029000	del	T	<DEL>	.	.	SVTYPE=DEL;END=5047000
```

and the same deletion written twice in Ensembl's default format, tab-separated,
once on each strand, read with `--format ensembl`:

```
21	5029001	5047000	DEL	+	del_fwd
21	5029001	5047000	DEL	-	del_rev
```

Both engines run `--plugin CADD,sv=<corpus>/plugin_data/1000G_phase3_SVs_synthetic.tsv.gz`.

### 2.3 What each engine writes

Both engines write the inversion's row without a `CADD_` field and the VCF
deletion's row with the synthetic file's score:

```
inv	21:5029001-5047000	inversion	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
del	21:5029001-5047000	deletion	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;CADD_PHRED=12.475;CADD_RAW=1.769367361985296
```

On the Ensembl-format records, Ensembl VEP scores the forward-strand deletion
and not the reverse-strand one:

```
del_fwd	21:5029001-5047000	deletion	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;CADD_PHRED=12.475;CADD_RAW=1.769367361985296
del_rev	21:5029001-5047000	deletion	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER
```

vep-rs scores both:

```
del_fwd	21:5029001-5047000	deletion	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;CADD_PHRED=12.475;CADD_RAW=1.769367361985296
del_rev	21:5029001-5047000	deletion	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;CADD_PHRED=12.475;CADD_RAW=1.769367361985296
```

Ensembl VEP's warnings file for the VCF records carries, for the inversion, two
Perl warnings about the undefined allele and the plugin's failure with its
exception block:

```
WARNING: 13 : Use of uninitialized value $orig_a_alt in hash element at /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/Sequence.pm line 1193, <__ANONIO__> line 4.
Use of uninitialized value $alt in numeric gt (>) at /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/Sequence.pm line 1252, <__ANONIO__> line 4.
WARNING: Plugin 'CADD' went wrong: 
-------------------- EXCEPTION --------------------
MSG: Missing reference or alternate sequence
STACK Bio::EnsEMBL::Variation::Utils::Sequence::trim_sequences /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/Sequence.pm:968
STACK Bio::EnsEMBL::Variation::Utils::Sequence::get_matched_variant_alleles /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/Sequence.pm:1203
STACK CADD::run /plugins/CADD.pm:251
```

The block runs on through the output factory's stack to a closing dashed line.
An inversion over a span with eight transcript rows throws on each of the eight
and writes the block once. Under `--fork` the block appears twice, the forked
process's own copy and the parent's relay of it. vep-rs writes no warning for
any of the records.

### 2.4 Why Ensembl VEP writes that

`CADD.pm`, lines 222-223, looks the record's Sequence Ontology class up in a
three-entry table (`insertion`, `deletion`, `duplication`). An inversion, a
tandem duplication or a copy-number class therefore gets an undefined allele
string. `trim_sequences` in `ensembl-variation`'s `Utils/Sequence.pm` (line 968)
throws on it for every row whose span the file covers. `OutputFactory.pm`, line
2244, catches each exception as the warning, which `BaseVEP.pm`, lines 809-813,
write once per distinct text. A record over a span with no file row never
reaches the match and throws nothing.

For an Ensembl-format record the plugin passes the record's strand to the allele
match (`CADD.pm`, line 244). `Utils/Sequence.pm`, lines 1182 and 1196,
reverse-complement the alleles of a record whose strand differs from the file
row's. `DEL` reversed is `LED`, and the complement table maps `D` to `H`, so the
key compared with the row's `DEL` is `LEH`, which matches nothing. `INS` becomes
`SNI` and `DUP` becomes `PUH`.

### 2.5 Why vep-rs's row is the correct one

A class the file does not score has no score, and the engine should write
nothing for it without an error. The record is well formed and every other
column of its rows is right. A class word is a label, not a sequence, so it has
no reverse complement. The record's strand changes nothing about the span
deleted, `21:5029001-5047000`. The plugin reads a row's class and start, no
strand (`CADD.pm` lines 265-278).

### 2.6 Try it yourself

```
vep -i records.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --plugin CADD,sv=<corpus>/plugin_data/1000G_phase3_SVs_synthetic.tsv.gz
vep -i records.txt -o ens.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --format ensembl \
  --plugin CADD,sv=<corpus>/plugin_data/1000G_phase3_SVs_synthetic.tsv.gz
```

Compare `out.txt_warnings.txt` from each engine: Ensembl VEP's carries the
block above for the inversion, vep-rs's nothing. In `ens.txt`, Ensembl VEP's
`del_rev` row ends `IMPACT=MODIFIER` and vep-rs's
`IMPACT=MODIFIER;CADD_PHRED=12.475;CADD_RAW=1.769367361985296`, the synthetic
file's score; the `del_fwd` rows agree.

### 2.7 How the comparator treats this class

The format comparator applies `cadd_structural_allele_match` when a comparison
case names it. A `CADD_` field the candidate fills on a structural allele's row
where the reference has none is set aside, counted once per row. In
`--warnings` mode the reference's `Plugin 'CADD' went wrong` lines are left out
of the comparison with their exception blocks and the two undefined-allele
warnings before them, each block counted once. A sequence allele's row, and a
value both engines fill differently, stay differences. The golden corpora place
their `<INV>` and `<DUP:TANDEM>` records where no CADD-SV row lies and carry no
Ensembl-format structural record, so the class applies to none of their rows.

## 3. SpliceVault: the table row of a sibling allele

**In short.** Ensembl VEP's SpliceVault plugin matches a table row against every
alternate allele of the record and takes the first row that matches any of them.
On a multi-allelic record every allele's rows therefore carry the first matching
row. vep-rs matches each row's own allele: an allele with a row of its own gets
that row, one without gets nothing.

### 3.1 What a reader needs to know first

A SpliceVault row names a site, an alternate allele and a transcript
(`SpliceVault.pm` line 201).
A multi-allelic VCF record lists several alternate alleles on one line, and
both engines write one set of rows per allele.

### 3.2 The record

A two-allele record on the site of the `GRCh38-plugin-splicevault-fixture`
golden corpus's record `donor_fwd_two_transcripts`, against the corpus's
synthetic table
(`tests/golden/116/GRCh38-plugin-splicevault-fixture/data/`), which holds a
row for each allele on `ENST00000299443`:

```
21	13618457	sibling	T	C,G	.	.	.
```

```
chr21	13618457	T	C	ENST00000299443	Acceptor_loss	chr21:13618451	0.51	Frameshift:0/2	Top1:CA;+449;0.493%;inFrame|Top2:CD;+44;30%;inFrame	44461	191
chr21	13618457	T	G	ENST00000299443	Acceptor_loss	chr21:13618443	0.24	Frameshift:2/4	Top1:ES;8-9;2 492%;Frameshift|Top2:CD;-94;22%;inFrame|Top3:CD;-416;17%;Frameshift|Top4:CD;-366;37%;inFrame	157839	122
```

Both engines run `--plugin SpliceVault,file=<table>`.

### 3.3 What each engine writes

The `G` allele's row on `ENST00000299443`. Ensembl VEP writes the `C` row's
values:

```
sibling	21:13618457	G	ENSG00000166351	ENST00000299443	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SpliceVault_SpliceAI_delta=0.51;SpliceVault_out_of_frame_events=Frameshift:0/2;SpliceVault_site_max_depth=191;SpliceVault_site_pos=chr21:13618451;SpliceVault_site_sample_count=44461;SpliceVault_site_type=Acceptor_loss;SpliceVault_top_events=Top1:CA:+449:0.493%:inFrame,Top2:CD:+44:30%:inFrame
```

vep-rs writes the `G` row's:

```
sibling	21:13618457	G	ENSG00000166351	ENST00000299443	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SpliceVault_SpliceAI_delta=0.24;SpliceVault_out_of_frame_events=Frameshift:2/4;SpliceVault_site_max_depth=122;SpliceVault_site_pos=chr21:13618443;SpliceVault_site_sample_count=157839;SpliceVault_site_type=Acceptor_loss;SpliceVault_top_events=Top1:ES:8-9:2_492%:Frameshift,Top2:CD:-94:22%:inFrame,Top3:CD:-416:17%:Frameshift,Top4:CD:-366:37%:inFrame
```

The `C` allele's rows agree between the engines.

### 3.4 Why Ensembl VEP writes that

`SpliceVault.pm`, line 156, takes the record's alternate alleles as a list and
tests each table row against the list. It returns the first row that matches any
of them. The row is then written for the allele being annotated, whichever
allele it named. On a one-allele record the list holds the row's allele alone.

### 3.5 Why vep-rs's row is the correct one

The synthetic table states the site's splicing events per alternate allele.
The `G` row's
events, depth and sample count describe the `G` allele, and the `C` row's
describe `C`. A row written for `G` with the `C` allele's numbers states
something the table does not say about `G`. Ensembl VEP writes the `G` row
itself when `G` is the only allele on the line.

### 3.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --dir_plugins <corpus>/data \
  --plugin SpliceVault,file=<corpus>/data/SpliceVault_data_GRCh38_synthetic.tsv.gz
```

On the synthetic table, Ensembl VEP's `G` rows carry
`SpliceVault_SpliceAI_delta=0.51` and
`SpliceVault_site_pos=chr21:13618451`; vep-rs's carry `0.24` and
`chr21:13618443`. Splitting the record into `T C` and `T G` makes Ensembl VEP
write vep-rs's rows.

### 3.7 How the comparator treats this class

The format comparator applies `splicevault_sibling_allele_rows` when a
comparison case names it. A `SpliceVault_` field, or the JSON `splicevault`
object, is set aside, counted once per row, where the reference's value is the
one its row of another allele of the record carries on the same feature. A field
the candidate alone fills, and any value on a record of one allele, stay
differences. The corpora carry one allele per record, so the class applies to
none of their rows.

## 4. SpliceVault: the JSON object of a record's later transcript rows

**In short.** In the JSON output format, Ensembl VEP's SpliceVault plugin writes
the `splicevault` object on the first row that reads a cached table row. It
writes none on the later rows reading the same cached row, with a
`Plugin 'SpliceVault' went wrong: Not an ARRAY reference` warning.

Under `transcript_match=0` those are a record's later transcript rows. Under the
default match they are the rows of a multi-allelic record's later alleles, and
of a later record at the same site with the same allele. Its default, tab and
VCF outputs carry the seven fields on every row. vep-rs writes the object on
every row, with the values its own text formats carry.

### 4.1 What a reader needs to know first

`transcript_match=0` drops the transcript test of the table lookup. Every
transcript row of the record takes the first table row at the site whose allele
matches, whatever transcript that row names. `transcript_match=1`, the default,
gives a row the table rows naming its own transcript. In the JSON format the
seven `SpliceVault_` fields nest under `splicevault`, and
`SpliceVault_top_events` becomes an object keyed by the event's rank.

Under the default match, too, several rows can read one table row. A
multi-allelic record's alleles share the first allele's row
([section 3](#3-splicevault-the-table-row-of-a-sibling-allele)), and a later
record at the site reads the parsed rows the plugin kept.

### 4.2 The record

The record `donor_fwd_two_transcripts` of the `GRCh38-plugin-splicevault-fixture`
golden corpus, against the corpus's synthetic table
(`tests/golden/116/GRCh38-plugin-splicevault-fixture/data/`), whose first `T`
to `C` row at the site is:

```
21	13618457	donor_fwd_two_transcripts	T	C	.	.	.
```

```
chr21	13618457	T	C	ENST00000299443	Acceptor_loss	chr21:13618451	0.51	Frameshift:0/2	Top1:CA;+449;0.493%;inFrame|Top2:CD;+44;30%;inFrame	44461	191
```

The record has three transcript rows: `ENST00000299443` and `ENST00000620442`
of `ENSG00000166351`, and the upstream snRNA `ENST00000384745`. Both engines
run `--plugin SpliceVault,file=<table>,transcript_match=0 --json`.

### 4.3 What each engine writes

Ensembl VEP's JSON record carries `splicevault` on the first transcript entry
alone; its second and third entries are:

```
{"transcript_id":"ENST00000384745","gene_id":"ENSG00000207476","impact":"MODIFIER","strand":1,"distance":3120,"variant_allele":"C","consequence_terms":["upstream_gene_variant"]}
{"transcript_id":"ENST00000620442","gene_id":"ENSG00000166351","impact":"HIGH","strand":1,"variant_allele":"C","consequence_terms":["splice_donor_variant"]}
```

and its warnings file reads:

```
WARNING: Plugin 'SpliceVault' went wrong: Not an ARRAY reference at /plugins/SpliceVault.pm line 183, <__ANONIO__> line 3.
```

vep-rs writes the object on all three entries; its third is:

```
{"consequence_terms":["splice_donor_variant"],"gene_id":"ENSG00000166351","impact":"HIGH","splicevault":{"SpliceVault_SpliceAI_delta":0.51,"SpliceVault_out_of_frame_events":"Frameshift:0/2","SpliceVault_site_max_depth":191,"SpliceVault_site_pos":"chr21:13618451","SpliceVault_site_sample_count":44461,"SpliceVault_site_type":"Acceptor_loss","SpliceVault_top_events":{"Top1":{"frame":"inFrame","percent_of_supporting_samples":"0.493%","site_type":"CA","transcript_impact":449},"Top2":{"frame":"inFrame","percent_of_supporting_samples":"30%","site_type":"CD","transcript_impact":44}}},"strand":1,"transcript_id":"ENST00000620442","variant_allele":"C"}
```

The object equals Ensembl VEP's first-entry object member for member. In the
default format both engines write the same three rows from the synthetic table,
each ending
`SpliceVault_SpliceAI_delta=0.51;...;SpliceVault_top_events=Top1:CA:+449:0.493%:inFrame,Top2:CD:+44:30%:inFrame`,
the rows the committed `expected_splicevault-transcript-match-0/default.txt.gz`
carries.

### 4.4 Why Ensembl VEP writes that

The plugin's table rows are parsed once per site and cached
(`BaseVepTabixPlugin.pm`, lines 361-369, return the same parsed rows for the
same region key). In the JSON format, `SpliceVault.pm`, lines 180-191, split the
cached row's `SpliceVault_top_events`, a list of
`rank:type:impact:percent:frame` strings (line 216), into an object keyed by
rank. They write that object back onto the cached row (line 190).

The next transcript row of the record reads the same cached row. Line 183,
`for my $event (@$top_events)`, meets the object where it expects a list and
dies. `OutputFactory.pm`, line 2244, catches the error as the warning and writes
the row without the object. The text formats return the cached row unchanged
(lines 192-194), so every row carries the fields there.

Two more readers meet the rewritten row under the default match. Lines 156-178
match a table row against every alternate allele, so on a record `T C,G` the `G`
rows take the `T` to `C` row the `C` rows rewrote, and die. The base class keeps
the last region's parsed rows (`BaseVepTabixPlugin.pm`, lines 558-570), so a
second `T C` record at the site dies on every transcript. A `T G` record there
reads its own rows.

### 4.5 Why vep-rs's row is the correct one

Under `transcript_match=0` the plugin's own rule hands the first table row
matching the site and allele to every transcript row. Ensembl VEP's default,
tab and VCF outputs write the seven fields on all three rows, and its JSON
writes the object on the first. The object the later rows lose is the same
object. vep-rs builds it from each row's own fields, so its JSON says what its
text formats say.

### 4.6 Try it yourself

```
vep -i record.vcf -o out.json --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --dir_plugins <corpus>/data \
  --plugin SpliceVault,file=<corpus>/data/SpliceVault_data_GRCh38_synthetic.tsv.gz,transcript_match=0 \
  --json
```

Ensembl VEP's record carries `splicevault` on `ENST00000299443` alone and its
`out.json_warnings.txt` the line above; vep-rs's record carries it on all three
transcript entries. Without `--json` both engines write the seven fields on
every row. Under `transcript_match=1` both write the object on the two
transcripts the table names for this one-allele record.

With the ALT written `C,G` under `transcript_match=1`, Ensembl VEP writes the
object on the `C` entries alone, with the same warning. A second `T C` record at
the site loses it on every entry. vep-rs writes each entry's own object.

### 4.7 How the comparator treats this class

The format comparator applies `splicevault_json_later_rows` when a comparison
case names it. A JSON `splicevault` object the candidate fills on a transcript
entry where the reference has none is set aside, counted once per entry. An
object both engines fill with different values, or one the reference alone
fills, stays a difference, and the text formats compare as usual. The corpus's
`transcript_match=0` flag set compares the default, tab and VCF formats and
leaves the JSON format out. The class therefore applies to none of the golden
corpora's rows.

## 5. TranscriptAnnotator: the table row of a sibling allele

**In short.** Ensembl VEP's TranscriptAnnotator plugin matches a table row
against every alternate allele of the record and returns the first row that
matches any of them. An allele without a row of its own therefore takes a
sibling allele's row. vep-rs matches each row's own allele.

### 5.1 What a reader needs to know first

The table names a chromosome, position, reference and alternate allele and a
transcript per row, then the columns to append. The plugin writes those columns
on the rows whose transcript the table names. The fields carry the file's name
as their prefix unless `prefix=` names another.

### 5.2 The record

A two-allele record on the site of the `GRCh38-plugin-transcript-annotator`
golden corpus's record `rs463608`, against the corpus's table, which holds a
row for the `T` allele alone:

```
21	10521702	rs463608	C	T,A	.	.	.
```

```
21	10521702	C	T	ENST00000618007	0.61	Tolerated	Intronic
```

Both engines run `--plugin TranscriptAnnotator,file=<table>`.

### 5.3 What each engine writes

The `A` allele's row on `ENST00000618007`. Ensembl VEP writes the `T` row's
columns:

```
rs463608	21:10521702	A	ENSG00000274391	ENST00000618007	Transcript	splice_region_variant,intron_variant	-	-	-	-	-	-	IMPACT=LOW;STRAND=1;transcript_annotator_Comment=Intronic;transcript_annotator_SIFT_pred=Tolerated;transcript_annotator_SIFT_score=0.61
```

vep-rs writes none:

```
rs463608	21:10521702	A	ENSG00000274391	ENST00000618007	Transcript	splice_region_variant,intron_variant	-	-	-	-	-	-	IMPACT=LOW;STRAND=1
```

The `T` allele's row agrees between the engines.

### 5.4 Why Ensembl VEP writes that

`TranscriptAnnotator.pm`, line 215, takes the record's alternate alleles as a
list. Line 227 hands the list to the allele match of each table row. Line 237
returns the first row that matched any of them on the right transcript.

### 5.5 Why vep-rs's row is the correct one

The table's row is about the `T` allele: its SIFT score and comment were
computed for that substitution. Written on the `A` allele's row they describe
a change the row does not report. Ensembl VEP writes nothing for `A` when `A`
is the only allele on the line.

### 5.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 \
  --plugin TranscriptAnnotator,file=<corpus>/data/transcript_annotator.tsv.gz
```

Ensembl VEP's `A` row carries the three `transcript_annotator_` fields;
vep-rs's carries none.

### 5.7 How the comparator treats this class

The golden harness (`crates/vep-cli/tests/common/mod.rs`) and the format
comparator apply `transcript_annotator_sibling_allele_rows` from a corpus
manifest. A `field_divergences` entry of the class names the row and the field
with the value vep-rs prints, and the comparison checks vep-rs's value against
it. No predicate decides the class from the outputs, since the field names carry
the user's prefix. The corpora carry one record per allele at the shared
positions, so the class applies to none of their rows.

## 6. gnomADMt: the `##INFO` descriptions of the file's fields

**In short.** Under `fields=all`, Ensembl VEP's header describes five gnomADMt
fields with their `##INFO` descriptions cut at the first `=` and stripped of
every `<` and `>`, because its VCF header parser mangles the lines. vep-rs
writes the file's text.

### 6.1 What a reader needs to know first

Each gnomADMt field's header line carries the description the file's `##INFO`
line gives it. Five of the file's descriptions hold a comparison sign:
`AC_het`, `AC_hom`, `AF_het` and `AF_hom` end `heteroplasmy level >= 0.95` or
`>= 0.10 and < 0.95`, and `common_low_heteroplasmy` holds `> 0 and < 0.50` and
`<0.01`.

### 6.2 The record

The first record of the `GRCh38-plugin-gnomadmt` golden corpus
(`tests/golden/116/GRCh38-plugin-gnomadmt/variants.vcf`), against the corpus's
slice of the gnomAD v3.1 chrM sites file, under
`--plugin gnomADMt,file=<corpus>/data/gnomad.v3.1.chrM.slice.vcf.bgz,fields=all`:

```
MT	3	.	T	C	.	.	.
```

### 6.3 What each engine writes

Ensembl VEP's header line for `gnomAD_AC_hom` and for
`gnomAD_common_low_heteroplasmy`:

```
## gnomAD_AC_hom : Allele count restricted to variants with a heteroplasmy level 
## gnomAD_common_low_heteroplasmy : Present if variant is found at an overall frequency of .001 across all samples with a heteroplasmy level  0 and  0.50 (includes variants 0.01 heteroplasmy which are subsequently filtered)
```

vep-rs's:

```
## gnomAD_AC_hom : Allele count restricted to variants with a heteroplasmy level >= 0.95
## gnomAD_common_low_heteroplasmy : Present if variant is found at an overall frequency of .001 across all samples with a heteroplasmy level > 0 and < 0.50 (includes variants <0.01 heteroplasmy which are subsequently filtered)
```

The record's first row is the same in both engines; with the five fields
selected (`fields=AC_het:AC_hom:AF_het:AF_hom:common_low_heteroplasmy`) it
reads:

```
MT_3_T/C	MT:3	C	ENSG00000198888	ENST00000361390	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=3304;STRAND=1;gnomAD_AC_het=1;gnomAD_AC_hom=19;gnomAD_AF_het=1.77198e-05;gnomAD_AF_hom=0.000336676;gnomAD_common_low_heteroplasmy=
```

Under `fields=all` every data row of the corpus is identical between the
engines.

### 6.4 Why Ensembl VEP writes that

The plugin reads the file's header through `ensembl-io`'s `BaseVCF4.pm`. Its
line 135 removes every `<` and `>` of a `##INFO` line, and its line 146 splits
each `key=value` part at `=`. A description is therefore kept up to its first
`=` and loses the characters the two lines strip. `heteroplasmy level >= 0.95`
is cut at the `=` of `>=`; `> 0 and < 0.50` and `<0.01` lose their signs.

### 6.5 Why vep-rs's row is the correct one

The header line's text is the file's description, which the gnomAD file states
in full; a sentence cut inside its comparison names no threshold. vep-rs's line
is the `##INFO` description as the file writes it.

### 6.6 Try it yourself

```
vep -i variants.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 \
  --plugin gnomADMt,file=<corpus>/data/gnomad.v3.1.chrM.slice.vcf.bgz,fields=all
```

Compare the `## gnomAD_` header lines of the two outputs: five differ as above,
every other line and every data row agree. Under the default field list the
five fields are not written and the headers agree.

### 6.7 How the comparator treats this class

No comparison of this repository reads a header line. The format comparator and
the golden harness compare the entries of the data rows, so nothing is set aside
and every row compares as usual. The `fields-all` flag set of the corpus
compares the JSON format, whose records carry no header, and keeps the
reference's text outputs beside it. `gnomadmt_header_descriptions` names what a
comparison of header lines would set aside.

## 7. LOEUF: a gene's minimum of `0`

**In short.** Under `match_by=gene`, Ensembl VEP's LOEUF plugin reports the
smallest score over the gene's transcripts. It drops a minimum whose text is
`0`, since the string `0` is false in Perl; vep-rs prints `LOEUF=0`. The
published gnomAD tables write a zero as `0.0000e+00`, which is true in Perl.
Only a table of the user's own that writes a bare `0` therefore reaches the
shape.

### 7.1 What a reader needs to know first

Under `match_by=gene` the plugin takes every table row whose gene is the
transcript's and keeps the ones whose LOEUF cell is a number. It reports the
smallest, as the table wrote it. Under `match_by=transcript` it reports the cell
of the row naming the transcript, whatever its value. Perl treats two strings as
false, the empty string and `0`; `0.0`, `0.00` and `0.0000e+00` are true.

### 7.2 The record

The record `app_exon` of the `GRCh38-plugin-loeuf` golden corpus
(`tests/golden/116/GRCh38-plugin-loeuf/variants.vcf`), on APP
(`ENSG00000142192`), against a copy of the corpus's gnomAD v4 table in which one
cell is changed: the `lof.oe_ci.upper` column (column 94) of the row for
`ENST00000359726`, the gene's minimum, from `2.5883e-01` to `0`. The copy is
re-indexed with `tabix -s 114 -b 115 -e 116`, the plugin's documented columns.

```
21	26112000	app_exon	G	A	.	.	.
```

Both engines run `--plugin LOEUF,file=<table>,match_by=gene`.

### 7.3 What each engine writes

The row on `ENST00000359726`. Ensembl VEP:

```
app_exon	21:26112000	A	ENSG00000142192	ENST00000359726	Transcript	synonymous_variant	291	99	33	G	ggC/ggT	-	IMPACT=LOW;STRAND=-1
```

vep-rs:

```
app_exon	21:26112000	A	ENSG00000142192	ENST00000359726	Transcript	synonymous_variant	291	99	33	G	ggC/ggT	-	IMPACT=LOW;STRAND=-1;LOEUF=0
```

Every one of the record's 55 APP transcript rows differs the same way: the
rows agree in every column, and vep-rs's end `LOEUF=0`. Against the unchanged
table both engines write `LOEUF=2.5883e-01` on all 55 rows.

### 7.4 Why Ensembl VEP writes that

`LOEUF.pm`, lines 191-203, keep the smallest numeric LOEUF over the gene's rows
as the string the table wrote (`looks_like_number`, line 192, admits `0`).
Line 204 returns it through `return $min_loeuf_score ? {LOEUF =>
$min_loeuf_score} : {}`, a truth test on that string. `0` is false, so the
plugin returns an empty result and the row gets no `LOEUF`; `0.0000e+00` is
true and prints. Under `match_by=transcript`, line 181 returns the matched
row's result without testing the value, so the same cell prints `LOEUF=0`
there.

### 7.5 Why vep-rs's row is the correct one

The table states the transcript's upper bound as `0`, the cell passes the
plugin's own numeric test, and it is the smallest value among the gene's rows.
The plugin's contract under `match_by=gene` is the gene's most severe score, and
a bound of `0` is the most severe a table can state. Ensembl VEP itself prints
the cell under `match_by=transcript`; only the gene branch's truth test loses
it.

### 7.6 Try it yourself

Write the table copy of 7.2 and annotate the record with each `match_by`:

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --plugin LOEUF,file=loeuf_zero.tsv.gz,match_by=gene
```

| Table cell, `match_by` | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `0`, `gene` | no `LOEUF` | `LOEUF=0` |
| `0`, `transcript` | `LOEUF=0` on `ENST00000359726` | the same |
| `0.0000e+00`, `gene` | `LOEUF=0.0000e+00` | the same |

### 7.7 How the comparator treats this class

The format comparator applies `loeuf_zero_minimum` when a comparison case names
it. A `LOEUF` the candidate fills with `0` where the reference has none is set
aside, counted once per row. Any other value, or a `LOEUF` the reference alone
fills, stays a difference. The corpus's table writes no bare `0`, so the class
applies to none of the golden corpora's rows.

## 8. OpenTargets: a cell of `0`

**The defect in one sentence.** In the default, tab and VCF formats, Ensembl
VEP's OpenTargets plugin writes `NA`, its marker for a cell the record leaves
empty, for a selected cell whose text is `0`.

**Why this is a defect.** The plugin joins a field's values across the matched
records with `$res->{$_} || "NA"`, a truth test on the cell's text. `0` is one
of the two strings Perl treats as false. The same module's JSON branch copies
the cell as written, so one export cell reads `0` in the JSON object and `NA` in
the text fields of the same run. `NA` stands for an association without a value
in that column; a locus-to-gene score or a p-value exponent of `0` is a value.

**What vep-rs does instead.** vep-rs writes `NA` for an empty or absent cell
alone and copies every other cell as written, `0` included, in every output
format.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `OpenTargets.pm`, lines 336-340, `_join_results`: the text-format branch
  pushes `$res->{$_} || "NA"` for each selected column; lines 330-334: the JSON
  branch pushes the record's cells unchanged.

**What correcting it would change.** A text-format field would read `0`
wherever the export writes `0`, as the JSON object of the same row does.

**The class it produces.** `opentargets_zero_cell`: an `OpenTargets_` field of
the default, tab or VCF format whose `,`-joined items agree but for an item the
reference writes as `NA` and the candidate as `0`; the JSON `opentargets`
object agrees and is not read. No committed row: the committed export slice
(`tests/golden/116/GRCh38-plugin-opentargets/`, 69 associations) read by the
corpus's four runs and by the four comparison cases carries no cell equal to
`0`. Its empty cells, 165 in the five default columns, print `NA` on both
engines.

**How often the shape occurs.** Tier 3 under the rule of [Intended divergences
from Ensembl VEP](../intended-divergences.md): ordinary data does not carry the
shape. None of the 69 rows of the committed export has a bare `0` in an output
column. 4 of the 392 committed flag sets run the plugin.

### 8.1 A locus-to-gene score of `0` on an insertion

**In short.** Against a copy of the committed export in which the
`gwasLocusToGeneScore` cell of one record, `0.84708172082901`, is changed to
`0`, the record `ot_gwas_insertion_lead` (`21 14979687 G GT`) has 15 transcript
rows. Ensembl VEP 116.2 writes `OpenTargets_gwasLocusToGeneScore=NA` on every
one; vep-rs writes `OpenTargets_gwasLocusToGeneScore=0`, and both write `0` in
the row's JSON object. Every other field of the rows agrees, the two empty QTL
cells printing `NA` on both engines; against the unchanged export both write
`0.84708172082901`.

## Where next

- [Annotation-store
  plugins](annotation-store.md#where-these-plugins-differ-from-ensembl-vep):
  three more classes, on REVEL, TranscriptAnnotator and gnomADMt, and the list
  of the five classes with pages of their own.
- [The AlphaMissense prediction of a sibling
  allele](../divergences/alphamissense-sibling-allele-rows.md): one of those
  five pages.
- [Intended divergences from Ensembl VEP](../intended-divergences.md): the index
  of every class.
