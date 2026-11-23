# The NMD escape of a `--minimal` split allele

**The defect in one sentence.** Under `--minimal`, Ensembl VEP's NMD plugin judges every allele of a split multi-allelic record by the rejoined record's end, and in the default and tab formats by the record's CDS coordinates as well, so an allele's escape depends on the alleles written beside it. vep-rs judges each split allele over its own trimmed span in every output format, so the answer equals the one the allele gets written alone.

**Why this is a defect.** The plugin's rules are positional facts of one allele: where its end lies against the transcript's last two exons, and how far into the coding sequence it reaches. `--minimal` exists to annotate each allele of a multi-allelic record on its own trimmed span, and Ensembl VEP computes the allele's terms, codons and amino acids there. The plugin then runs after the rejoin and reads the record's span, which the other alleles set. The same allele written alone on a record gets the other answer.

**What vep-rs does instead.** vep-rs judges each split allele over its own trimmed span in every output format: its own genomic end against the last exon and the penultimate exon's window, and its own CDS coordinates: both ends defined, and the end against coding base 101. The alleles beside it change nothing, and the answer equals the one the allele gets written alone. The plugin's rules are described under [NMD](../plugins.md).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `NMD.pm`, line 137: `my $vf_end = $vf->seq_region_end;` is the position tested against the last exon (line 148) and the penultimate exon's window (lines 159-161 and 165-167). Line 88 reads `$tv->cds_end`, which line 98 tests against 101, and lines 103-105 return nothing unless both `$tv->cds_start` and `$tv->cds_end` are defined.
- `OutputFactory.pm`, lines 2293-2325: the rejoin annotates each split feature on its own span (lines 2298 and 2312), then writes the record's span back onto the first feature (lines 2300-2302) and points every merged feature's TranscriptVariations at it (line 2323). The output pass that runs the plugin follows, so `$vf->seq_region_end` is the record's end in every format.
- `OutputFactory.pm`, lines 294-302: the default and tab formats rejoin, then clear every cached CDS, cDNA and protein position (`clear_shifting_variables`, `TranscriptVariationAllele.pm` lines 509-517), so the output pass maps `$tv->cds_start` and `$tv->cds_end` from the record's span (`BaseTranscriptVariation.pm` lines 505-513). `OutputFactory/VCF.pm` lines 305-308 and `OutputFactory/JSON.pm` lines 217-220 clear first and rejoin after, so those formats keep the CDS coordinates cached on the allele's own span.

**What correcting it would change.** The plugin would read each allele's own end and CDS coordinates, and a split allele's `NMD` value would equal the value it gets on a record of its own, in every format.

**The class it produces.** `nmd_minimal_rejoined_end`: an `NMD` value that differs between the engines on a row of a multi-allelic record whose alleles trim to different spans under `--minimal`. No comparison tool sets a row of the shape aside: no committed flag set runs the plugin with `--minimal`, so a disagreement of this shape would be reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape the record alone shows is a multi-allelic record whose sequence alleles trim to different spans, 1,894 of the cuts' 15,646,481 records (121.0 per million), in 27 corpora. The transcript decides whether an allele's own span and the record's differ in what the plugin reads, the genomic end or the CDS coordinates, and whether the two fall on different sides of a rule. Of the committed flag sets, none runs the plugin with `--minimal`.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion two bases short of the penultimate exon's window

**In short.** `21 36071003 acgta_agta_acgtc ACGTA AGTA,ACGTC` on GRCh38 lies in the second of CBR1's three exons on ENST00000290349, whose penultimate exon ends at 36071057. Under `--minimal`, `AGTA` trims to the deletion of `C` at 36071004, two bases before the plugin's window, the 52 positions within 51 bases of that end, starts at 36071006, and `ACGTC` to a substitution at 36071007, inside it. Ensembl VEP 116.3 writes `NMD_escaping_variant` on the deletion's frameshift rows of 11 transcripts; vep-rs on one, whose last exon holds 36071004.

### 1.1 What a reader needs to know first

The NMD plugin marks a stop-gained, frameshift or splice-site row `NMD_escaping_variant` when the variant's end lies in the transcript's last exon or within 51 bases of the penultimate exon's 3' end, or when its CDS end is at or before coding base 101. A row whose CDS position is open at either end is never marked. Under `--minimal` Ensembl VEP trims each allele of a multi-allelic record, annotates the trimmed alleles one at a time and rejoins their rows under the record's span. The `Location` column prints the record's span.

### 1.2 The record

A multi-allelic record against the Ensembl 116 GRCh38 cache. Both engines drop the shared first base, so the record spans 36071004-36071007, coding bases 344-347 of ENST00000290349. The first alternate allele is the deletion of coding base 344 at 36071004, a frameshift; the second changes coding base 347 at 36071007, a missense change.

```
21	36071003	acgta_agta_acgtc	ACGTA	AGTA,ACGTC	.	.	.
```

Both engines run `--minimal --plugin NMD`.

### 1.3 What each engine writes

The deletion's rows on ENST00000290349 and ENST00000439427 in the default format. Ensembl VEP 116.3:

```
acgta_agta_acgtc	21:36071004-36071007	-	ENSG00000159228	ENST00000290349	Transcript	frameshift_variant	436-439	344-347	115-116	F/X	tTc/tc	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;MINIMISED=1;NMD=NMD_escaping_variant
acgta_agta_acgtc	21:36071004-36071007	-	ENSG00000159228	ENST00000439427	Transcript	frameshift_variant	613-616	344-347	115-116	F/X	tTc/tc	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;MINIMISED=1;NMD=NMD_escaping_variant
```

vep-rs:

```
acgta_agta_acgtc	21:36071004-36071007	-	ENSG00000159228	ENST00000290349	Transcript	frameshift_variant	436-439	344-347	115-116	F/X	tTc/tc	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;MINIMISED=1
acgta_agta_acgtc	21:36071004-36071007	-	ENSG00000159228	ENST00000439427	Transcript	frameshift_variant	613-616	344-347	115-116	F/X	tTc/tc	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;MINIMISED=1
```

The record writes 38 rows in each engine, 19 per allele, equal in every column but `NMD`. The substitution's rows are missense and carry no value in either engine. The deletion's `NMD` cell by transcript, in every output format:

| Transcript | Penultimate exon's window | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| ENST00000290349 and nine others | 36071006-36071057 | `NMD_escaping_variant` | none |
| ENST00000439427 | its last exon spans 36070950-36071583 | `NMD_escaping_variant` | `NMD_escaping_variant` |
| ENST00001125702, ENST00001125703 | 36071010-36071061 | none | none |

The nine transcripts with the first window are ENST00000399191, ENST00000530908, ENST00000993894, ENST00000993895, ENST00001121351, ENST00001125699, ENST00001125700, ENST00001125701 and ENST00001131789. Ensembl VEP marks 11 rows and vep-rs one.

### 1.4 Why Ensembl VEP writes that

The rejoin writes the record's span back onto the first allele's feature and points the second allele's TranscriptVariations at it, so when the plugin runs, `$vf->seq_region_end` is 36071007 on every row of both alleles. 36071007 lies inside the window 36071006-36071057 of ten transcripts and inside ENST00000439427's last exon, so the deletion's eleven frameshift rows there escape; the two transcripts whose window starts at 36071010 do not. The substitution's rows carry no gate term, so the plugin returns nothing on them whatever the end.

### 1.5 Why vep-rs's row is the correct one

The deletion's end is a fact of its own span: it removes base 36071004 and nothing after it. Against the ten transcripts whose window starts at 36071006 it lies two bases short, and only ENST00000439427's last exon contains it. Ensembl VEP agrees when the deletion is given alone: the bi-allelic record `21 36071003 ACGTA AGTA` under `--minimal` is trimmed in place to the same deletion and writes `NMD_escaping_variant` on the ENST00000439427 row alone in both engines. vep-rs's rows inside the multi-allelic record are those rows.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --minimal --plugin NMD
```

vep-rs writes `NMD=NMD_escaping_variant` on the deletion's ENST00000439427 row alone. Ensembl VEP 116.3 writes it on eleven of the deletion's rows, the ENST00000290349 row included. Writing the deletion on a record of its own gives one row in both engines.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`tools/vep-formats-compare`) has no class of this shape, and the golden harness compares the `NMD` cell of every row of the four corpora that run the plugin, whose single-allele records never put the shape to the test. A disagreement of the shape would be reported as a difference and charged to vep-rs.

## 2. A deletion of coding base 101 beside a deletion of bases 101-102

**In short.** `21 36070215 gtg_gg_g GTG GG,G` on GRCh38 lies in CBR1's first exon. Under `--minimal`, `GG` trims to the deletion of `T` at 36070216, coding base 101 of ENST00000290349, and `G` to the deletion of `TG` at 36070216-36070217, coding bases 101-102. The plugin's first rule marks a CDS end at or before 101. Ensembl VEP 116.3 marks the one-base deletion's 14 frameshift rows in the VCF and JSON formats and none of them in the default and tab formats; vep-rs marks the 14 in every format, and neither engine marks the two-base deletion.

### 2.1 What a reader needs to know first

Ensembl VEP computes a split allele's CDS position on its trimmed span and caches it. Its default and tab formats clear that cache after the rejoin and recompute the position over the record's span, which is what their `CDS_position` column prints; its VCF and JSON formats clear the cache before the rejoin and print the allele's own position. vep-rs prints the same columns as each format of Ensembl VEP. The plugin reads whichever CDS end its engine holds when it runs.

### 2.2 The record

```
21	36070215	gtg_gg_g	GTG	GG,G	.	.	.
```

Both engines drop the shared first base, so the record spans 36070216-36070217, coding bases 101-102. Both engines run `--minimal --plugin NMD`, and `--tab`, `--vcf` or `--json` for the other formats.

### 2.3 What each engine writes

The two alleles' rows on ENST00000290349 in the default format. Ensembl VEP 116.3:

```
gtg_gg_g	21:36070216-36070217	-	ENSG00000159228	ENST00000290349	Transcript	frameshift_variant	193-194	101-102	34	V/X	gTg/gg	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;MINIMISED=1
gtg_gg_g	21:36070216-36070217	-	ENSG00000159228	ENST00000290349	Transcript	frameshift_variant	193-194	101-102	34	V/X	gTG/g	-	ALLELE_NUM=2;IMPACT=HIGH;STRAND=1;MINIMISED=1
```

vep-rs:

```
gtg_gg_g	21:36070216-36070217	-	ENSG00000159228	ENST00000290349	Transcript	frameshift_variant	193-194	101-102	34	V/X	gTg/gg	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;MINIMISED=1;NMD=NMD_escaping_variant
gtg_gg_g	21:36070216-36070217	-	ENSG00000159228	ENST00000290349	Transcript	frameshift_variant	193-194	101-102	34	V/X	gTG/g	-	ALLELE_NUM=2;IMPACT=HIGH;STRAND=1;MINIMISED=1
```

The same two alleles' `CSQ` entries for the transcript in the VCF format, where the two engines agree:

```
-|frameshift_variant|HIGH|CBR1|ENSG00000159228|Transcript|ENST00000290349|protein_coding|1/3||||193|101|34|V/X|gTg/gg||1||1||1|HGNC|HGNC:1548|NMD_escaping_variant
-|frameshift_variant|HIGH|CBR1|ENSG00000159228|Transcript|ENST00000290349|protein_coding|1/3||||193-194|101-102|34|V/X|gTG/g||2||1||1|HGNC|HGNC:1548|
```

The record writes 38 rows in each engine, 19 per allele. The one-base deletion (`ALLELE_NUM=1`) has 14 frameshift rows and the two-base deletion (`ALLELE_NUM=2`) 14; the rows carrying `NMD_escaping_variant`, by format, with the `CDS_position` both engines print for the allele:

| Format and allele | `CDS_position` | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| default and tab, one base | `101-102` | none | 14 rows |
| VCF and JSON, one base | `101` | 14 rows | 14 rows |
| every format, two bases | `101-102` | none | none |

### 2.4 Why Ensembl VEP writes that

The one-base deletion's own CDS end is 101, computed and cached when the rejoin annotates the split feature. In the default and tab formats the output pass clears that cache after the rejoin and maps the record's span, 36070216-36070217, to coding bases 101-102. `$tv->cds_end` is then 102 and the first rule does not fire; the record's end lies in the first exon, so neither does any other. In the VCF and JSON formats the cache is cleared before the rejoin and never again, so the plugin reads the cached 101 and the rows escape. The two-base deletion's own span is the record's, so its CDS end is 102 in every format.

### 2.5 Why vep-rs's row is the correct one

The one-base deletion removes coding base 101 and reaches no further; a CDS end of 102 belongs to the allele beside it. The plugin's own rule marks a CDS end at or before 101, and Ensembl VEP applies it to this allele's own end in its VCF and JSON formats and to the bi-allelic record `21 36070215 GTG GG` in every format. vep-rs's rows inside the multi-allelic record are those rows: the one-base deletion escapes and the two-base deletion does not, whatever the output format.

### 2.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --minimal --plugin NMD
```

vep-rs writes `NMD=NMD_escaping_variant` on the one-base deletion's 14 frameshift rows and nothing on the two-base deletion's, in every format. Ensembl VEP 116.3 writes the same 14 under `--vcf` and `--json` and none in the default and tab formats. Writing `GTG GG` on a record of its own gives the 14 rows in both engines.

### 2.7 How the comparator treats this class

As in [1.7](#17-how-the-comparator-treats-this-class): nothing is set aside, and a disagreement of either shape would be reported as a difference and charged to vep-rs.
