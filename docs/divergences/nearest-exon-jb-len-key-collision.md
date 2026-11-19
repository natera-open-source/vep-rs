# Every overlapped exon reported as `len` under NearestExonJB's `max_range=len`

**The defect in one sentence.** Ensembl VEP's NearestExonJB plugin files each exon's distances and its length in one hash, the length under the key `len`, and tests the range value as a key of that hash, so `max_range=len` reports every overlapped exon with the word `len` as its distance. vep-rs reads the value as 0, warns that it is not a number, and writes no `NearestExonJB` value on any row.

**Why this is a defect.** `max_range` is a number of bases: the plugin reports the exons whose nearer junction lies within it. Ensembl VEP's numeric comparison reads the word `len` as 0 and warns that it is not numeric, so no exon can lie within the range. The lookup that chooses the exons to report is a hash-key test on the same value, and `len` is the key the plugin itself stores every exon's length under. The item it then prints is `<exon>|len|<length>|<length>`: the word in the distance place, the exon's length in the junction place and again in its own.

**What vep-rs does instead.** vep-rs reads `max_range=len` as Ensembl VEP's numeric comparison reads it, as 0, and warns that the value is not a number, in the warnings file beside the output and on standard error unless `--quiet`. A range of 0 admits no distance below it, and an exon at exactly the range is reported only when the range is written as a plain whole number. Under `max_range=len` no row carries a `NearestExonJB` value; the header line is written and the run goes on. The plugin's parameters and fields are described under [NearestExonJB](../plugins/transcript-model.md#nearestexonjb).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `NearestExonJB.pm`, line 186: `my $min = $CONFIG{max_range}`, the parameter's text, `len`.
- Lines 225-239: for each overlapped exon the plugin stores its distance under `$dists{$exon_id}{$startD}` or `{$endD}` (lines 228, 232 and 236) and its length under `$dists{$exon_id}{len}` (lines 229, 233 and 237). It then lowers `$min` with `$min = $startD if $min > $startD` (lines 230, 234 and 238). The comparison numifies `len` to 0, warns `Argument "len" isn't numeric in numeric gt (>)`, and never lowers it, since no distance is below 0.
- Lines 267-271: `if (exists $dists{$exon}{$min})` tests the text `len` as a key, which every exon's hash holds, and pushes `$exon.$char_sep.$min.$char_sep.$dists{$exon}{$min}.$char_sep.$dists{$exon}{len}`, the exon, the word, its length and its length again.
- Lines 244-264, the `intronic=1` branch with two flanking exons, read the `len` key by name and are not affected.

**What correcting it would change.** The plugin would keep the lengths apart from the distances, or test the range as a number, and `max_range=len` would report no exon, as `max_range=abc` already does in both engines.

**The class it produces.** `nearest_exon_jb_len_key_collision`: a `NearestExonJB` value Ensembl VEP writes under `max_range=len` on a row with an overlapped exon, where vep-rs writes none. No comparison tool sets a row of the shape aside: no committed flag set passes the value, so a disagreement of this shape would be reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is the parameter value `max_range=len` itself, a word no ordinary run passes, and no committed flag set carries it. The rows it reaches are every row with an overlapped exon: on the 47 records of the committed corpus that runs the plugin, 961 of 1,606 rows.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion 23 bases into an exon under `max_range=len`

**In short.** `21 36896882 1069237 TG T` on GRCh38 deletes the base at 36896883, 23 bases from the start of ENSE00001282176, an exon of 272 bases on ENST00000336648. Under the default range both engines write `ENSE00001282176|23|start|272` on that row. Under `--plugin NearestExonJB,max_range=len`, Ensembl VEP 116.2 writes `ENSE00001282176|len|272|272` and vep-rs writes nothing.

### 1.1 What a reader needs to know first

NearestExonJB reports, for a row whose variant overlaps an exon, the exon's stable id, the distance to its nearer junction, which junction that is and the exon's length, as `ExonID|distance|start/end|length`. `max_range` bounds the distance, 10000 bases by default, and only the exons at the smallest distance are reported. A value that is not a number is read as 0 by Ensembl VEP's own comparison, with a warning per row in the warnings file, and vep-rs reads it the same way with one warning.

### 1.2 The record

A one-base deletion against the Ensembl 116 GRCh38 cache. The reference is `TG` at 36896882-36896883; the deleted `G` at 36896883 is 23 bases from the start of ENSE00001282176, in genomic order, on the reverse-strand transcript ENST00000336648.

```
21	36896882	1069237	TG	T	.	.	.
```

Both engines run `--plugin NearestExonJB,max_range=len`; the default range is `--plugin NearestExonJB`.

### 1.3 What each engine writes

The row on ENST00000336648 in the default format. Under the default range, both engines:

```
1069237	21:36896883	-	ENSG00000159267	ENST00000336648	Transcript	frameshift_variant	1899	1428	476	P/X	ccC/cc	-	IMPACT=HIGH;STRAND=-1;NearestExonJB=ENSE00001282176|23|start|272
```

Under `max_range=len`, Ensembl VEP 116.2:

```
1069237	21:36896883	-	ENSG00000159267	ENST00000336648	Transcript	frameshift_variant	1899	1428	476	P/X	ccC/cc	-	IMPACT=HIGH;STRAND=-1;NearestExonJB=ENSE00001282176|len|272|272
```

vep-rs:

```
1069237	21:36896883	-	ENSG00000159267	ENST00000336648	Transcript	frameshift_variant	1899	1428	476	P/X	ccC/cc	-	IMPACT=HIGH;STRAND=-1
```

The record writes 17 rows in each engine. Ensembl VEP writes an item on all 17, each with `len` in the distance place and the exon's length twice; vep-rs writes one on none. Over the corpus's 47 records, Ensembl VEP writes a value on the 961 rows that carry one under the default range, 954 of them one item and 7 a deletion's three overlapped exons, where the default range names the two nearer ones. Its warnings file holds one `Argument "len" isn't numeric in numeric gt (>)` line per row with an overlapped exon, the first comparison of each row's own range value, 961 in all; vep-rs writes one line in its warnings file:

```
WARNING: NearestExonJB: max_range='len' is not a number; it is read as 0, as Ensembl VEP reads it
```

### 1.4 Why Ensembl VEP writes that

Line 186 sets `$min` to the text `len`. For the overlapped exon, line 228 stores the start distance 23 under `$dists{ENSE00001282176}{23}` and line 229 stores the length under `$dists{ENSE00001282176}{len}`; line 230 compares `len > 23` as numbers, reads `len` as 0, warns, and leaves `$min` as `len`. Line 268 then asks whether the exon's hash has a key `len`, which it does for every exon, and line 269 prints the exon, `$min`, the value under that key and the length: `ENSE00001282176|len|272|272`. A deletion over three exons prints all three, since each exon's hash holds the key.

### 1.5 Why vep-rs's row is the correct one

The range is a number of bases, and the only number Ensembl VEP reads from `len` is 0, the one its comparison uses. An exon 23 bases away lies outside a range of 0, so the plugin's own rule reports nothing, as both engines report nothing under `max_range=abc`, a word that names no key. The item Ensembl VEP writes is not a junction distance: its distance place holds a word and its junction place holds the exon's length. vep-rs's row, with no value, is the row the plugin's rule gives.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin NearestExonJB,max_range=len
```

vep-rs writes the `## NearestExonJB` header line, one warning and no `NearestExonJB` value on any row. Ensembl VEP 116.2 writes `NearestExonJB=ENSE00001282176|len|272|272` on the row of ENST00000336648 and an item of the same shape on every row with an overlapped exon. Running both engines with `max_range=abc` shows them agree: no value on any row.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`scripts/concordance/compare_vep_formats.py`) has no class of this shape, and the golden harness compares the `GRCh38-plugin-computed` corpus's `NearestExonJB` cells under its committed flag sets, none of which passes `max_range=len`. A disagreement of the shape would be reported as a difference and charged to vep-rs.
