# The NearestExonJB exon set of a shifted indel under `--shift_3prime 1`

**The defect in one sentence.** Under `--shift_3prime 1`, Ensembl VEP 116.2's NearestExonJB plugin takes the exons an insertion or deletion overlaps from the shifted span while it measures every distance from the record's own span. An indel the shift carries off an exon's last base therefore gets no value, and one carried onto an exon's first base gets that exon at a distance measured from outside it. vep-rs reads the exon set at the record's span, the span its distances measure from, so the first indel keeps its value and the second gets none.

**Why this is a defect.** The plugin's rule is one span's: the exons the variant overlaps, each with the distance from the variant to its nearer junction. Everything the plugin reads for itself is the record's span: its distances, the key of its result cache, the exon set it computes when no consequence test ran before it, and the exon set its `gff3=` option reads from a file. The exon set the consequence tests leave behind is the only shifted reading, and it belongs to another variant than the distances. Given the same exons through `gff3=`, Ensembl VEP 116.2 writes the value on every row its cache run leaves empty.

**What vep-rs does instead.** vep-rs reads the exon set at the record's span, the span its distances measure from, so a row's `NearestExonJB` describes one variant, the record's. The row's terms and positions describe the shifted span, as Ensembl VEP's do, and `SHIFT_LENGTH` relates the two. The plugin's parameters and fields are described under [NearestExonJB](../plugins/transcript-model.md#nearestexonjb).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-variation` at release/116, `ensembl-vep` at release/116.2):

- `NearestExonJB.pm`, line 184: the exon set is `$trv->_overlapped_exons`, called without a span. Lines 225-226 measure `$vf->start` and `$vf->end`, the record's span; line 179 keys the result cache by the same span; lines 279-288, `_gff3_overlapped_exons`, read the `gff3=` exons over `$vf->{start}` and `$vf->{end}`.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 848-870, `_overlapped_exons`: computed at the first call and memoised on the TranscriptVariation (line 851); a call without a span reads the feature's `start` and `end` (lines 853-855), the record's span.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 454-463, `_bvfo_preds`: the span the consequence tests read adds the shift to the feature's `start` and `end` (lines 460-462), and line 514 calls `_overlapped_exons($min_vf, $max_vf)` with it, filling the memo at the shifted span. Lines 269-282: each predicate's call moves the feature's `start` and `end` by the shift and back, so the feature the plugin reads keeps the record's span.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 395-398: a row's consequence terms are computed before its plugins run, so the memo the plugin reads at `NearestExonJB.pm` line 184 is the one the tests filled.

**What correcting it would change.** The plugin would read its exon set at the span its distances measure from, and a shifted indel would get the value the same record gets without the flag, which its `gff3=` option writes with it.

**The class it produces.** `nearest_exon_jb_shifted_span`: a `NearestExonJB` cell that differs between the engines on a transcript row whose `SHIFT_LENGTH` is positive, which only an insertion or deletion is given: a value in one engine against none in the other, or exons the other engine's cell lacks. The format comparator sets the cell aside where a comparison case names the class. No committed flag set runs the plugin with `--shift_3prime 1`, so no committed row has the shape. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape the record alone shows is an insertion or deletion, 2,154,547 of the cuts' 15,646,481 records (137,702 per million), in 104 corpora. The transcript decides whether the shift carries it across an exon's edge. No committed flag set runs the plugin with `--shift_3prime 1`.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. The last base of a DYRK1A exon deleted, the shift carrying the deletion into the intron

**In short.** `21 37480825 shift_probe AG A` on GRCh38 deletes the `G` at 37480826, the last base of ENSE00003584551, the 189-base exon at 37480638-37480826 of ENST00000338785, a forward-strand transcript of DYRK1A. The bases at 37480823-37480831 are `ACAGGTAAT`, so the first intron base, 37480827, is also `G`, and under `--shift_3prime 1` the deletion moves one base, onto the intron's first base and its donor site. Without the flag both engines write 27 identical rows, 25 with a value. Under `--plugin NearestExonJB --shift_3prime 1 --shift_length` Ensembl VEP 116.2 writes a value on 2 of the 27 rows and vep-rs on 25: the 23 `splice_donor_variant` rows carry a value in vep-rs, `ENSE00003584551|0|end|189` on 22 and `ENSE00003675684|0|end|189` on ENST00000643355, and nothing in Ensembl VEP.

### 1.1 What a reader needs to know first

NearestExonJB reports, for every exon a variant overlaps, the exon's stable id, the distance from the variant to its nearer junction, that junction and the exon's length. A row with no overlapped exon gets no value unless `intronic=1` is set and its terms include `intron_variant`. Under `--shift_3prime 1` an insertion or deletion is moved as far along the transcript's 3' direction as the repeated sequence allows before its terms, positions and HGVS are read; `Location` keeps the record's span and `SHIFT_LENGTH` says how far the row's span moved. A deletion of an exon's last base whose next base matches moves off the exon, and its row becomes a splice-donor variant of the intron.

### 1.2 The record

A deletion against the Ensembl 116 GRCh38 cache. The reference is `AG` at 37480825-37480826; 37480826 is the last base of ENSE00003584551 and 37480827, the next `G`, the first base of the intron.

```
21	37480825	shift_probe	AG	A	.	.	.
```

Both engines run `--plugin NearestExonJB --shift_3prime 1 --shift_length` with the reference FASTA, which the shift reads.

### 1.3 What each engine writes

The rows on ENST00000338785 and ENST00000999350 in the default format. Ensembl VEP 116.2:

```
shift_probe	21:37480826	-	ENSG00000157540	ENST00000338785	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1
shift_probe	21:37480826	-	ENSG00000157540	ENST00000999350	Transcript	frameshift_variant	2139	490	164	V/X	Gta/ta	-	IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1;NearestExonJB=ENSE00004532916|188|start|1703
```

vep-rs:

```
shift_probe	21:37480826	-	ENSG00000157540	ENST00000338785	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1;NearestExonJB=ENSE00003584551|0|end|189
shift_probe	21:37480826	-	ENSG00000157540	ENST00000999350	Transcript	frameshift_variant	2139	490	164	V/X	Gta/ta	-	IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1;NearestExonJB=ENSE00004532916|188|start|1703
```

Without the two shift flags, both engines write the ENST00000338785 row as:

```
shift_probe	21:37480826	-	ENSG00000157540	ENST00000338785	Transcript	frameshift_variant,splice_region_variant	2266	516	172	Q/X	caG/ca	-	IMPACT=HIGH;STRAND=1;NearestExonJB=ENSE00003584551|0|end|189
```

The record writes 27 rows in each engine under the shift flags, equal in every column but `NearestExonJB`:

| Rows | Consequence | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- | --- |
| 22 | `splice_donor_variant` | none | `ENSE00003584551\|0\|end\|189` |
| 1, ENST00000643355 | `splice_donor_variant,non_coding_transcript_variant` | none | `ENSE00003675684\|0\|end\|189` |
| 1, ENST00000999350 | `frameshift_variant` | `ENSE00004532916\|188\|start\|1703` | the same |
| 1, ENST00000462274 | `non_coding_transcript_exon_variant` | `ENSE00003828235\|188\|start\|1012` | the same |
| 2 | `intron_variant,non_coding_transcript_variant`; `downstream_gene_variant` | none | none |

On ENST00000999350 the shift moved `cDNA_position` from 2138 to 2139, `CDS_position` from 489 to 490 and `Protein_position` from 163 to 164, while the distance stays 188: from 37480826, the record's position, to 37480638, the exon's start.

### 1.4 Why Ensembl VEP writes that

The row's consequence terms are computed first, over the shifted span 37480827, and that computation fills the TranscriptVariation's memoised exon set with the exons overlapping 37480827. For ENST00000338785 none does: 37480827 is the first base of the intron, so the plugin's list is empty and it writes nothing. ENSE00004532916 of ENST00000999350 and ENSE00003828235 of ENST00000462274 start at 37480638, like ENSE00003584551, and end past 37480827, so their rows keep an exon. The plugin then measures that exon's distance from the record's span, 188 from 37480826 to the exon's start, while the same row's positions describe 37480827. The distance is the one the row gets without the flag.

### 1.5 Why vep-rs's row is the correct one

An item names the distance from the variant to the junction of an exon the variant overlaps, so its exon set and its distance describe one span. The plugin measures from the record's span, so the exon set is the record's, and Ensembl VEP's own plugin agrees whenever its exons come from anywhere but the memo. Given the same exons through `gff3=` (a GFF3 of the transcripts' `exon` lines, each with `ID=exon:<exon id>;Parent=transcript:<transcript id>`, compressed and indexed for tabix), Ensembl VEP 116.2 under the same flags writes `ENSE00003584551|0|end|189` on the 22 rows and `ENSE00003675684|0|end|189` on ENST00000643355: all 27 rows equal vep-rs's. Without the flag both engines write those values from the cache.

A deletion the shift moves inside the exon shows which span the distance has. `21 37480756 inside_exon TA T` deletes the `A` at 37480757, and the shift carries it five bases along `AAAAAA`, to 37480762, still inside ENSE00003584551. Both engines write `ENSE00003584551|69|end|189` on ENST00000338785's row under `SHIFT_LENGTH=5`: 69 is the record's distance, from 37480757 to 37480826, where the shifted span's would be 64.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --fasta <reference> --plugin NearestExonJB --shift_3prime 1 --shift_length
```

vep-rs writes `NearestExonJB=ENSE00003584551|0|end|189` on the `splice_donor_variant` row of ENST00000338785; Ensembl VEP 116.2 writes no `NearestExonJB` on it. Without the two shift flags both engines write the value, and Ensembl VEP with `--plugin NearestExonJB,gff3=<exons>` writes it under them.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies `nearest_exon_jb_shifted_span` when a comparison case names it. A `NearestExonJB` cell on a transcript row whose `SHIFT_LENGTH` is positive is set aside, counted once per row, where one engine fills it and the other leaves it empty, or where the two values name different exons. Any other difference stays one: two values that name the same exons, and any cell on a row whose `SHIFT_LENGTH` is zero or missing. The comparator reads the shift from the row's own `SHIFT_LENGTH`, which a run prints under `--shift_length`, so the class applies to runs with that flag; under `--shift_3prime 1` alone the same difference stays one.

No committed flag set runs the plugin with `--shift_3prime 1`; the golden harness compares the `GRCh38-plugin-computed` corpus's `NearestExonJB` cells as multisets on runs without the flag, so the class applies to none of its rows.

## 2. The acceptor's last base deleted, the shift carrying the deletion onto the exon

**In short.** `21 37480636 onto_exon AG A` deletes the `G` at 37480637, the last intron base before ENSE00003584551 and the second base of its acceptor site `AG`. The exon's first base, 37480638, is also `G`, so under `--shift_3prime 1` the deletion moves onto it. Without the flag both engines write 27 rows and no value, since the record's span overlaps no exon. Under `--plugin NearestExonJB --shift_3prime 1 --shift_length` Ensembl VEP 116.2 writes a value on 26 of the 27 rows, `ENSE00003584551|1|start|189` on ENST00000338785's `frameshift_variant,splice_region_variant` row, and vep-rs on none.

### 2.1 What a reader needs to know first

The 3' shift of a forward-strand transcript moves a deletion to higher coordinates. A deletion of the last intron base before an exon, when the exon's first base repeats it, therefore moves onto the exon, and its row becomes a frameshift of the exon's first base with `SHIFT_LENGTH=1`. The record's span stays in the intron, one base before the exon.

### 2.2 The record

The reference is `AG` at 37480636-37480637, the acceptor site; 37480638 is the exon's first base.

```
21	37480636	onto_exon	AG	A	.	.	.
```

Both engines run `--plugin NearestExonJB --shift_3prime 1 --shift_length` with the reference FASTA.

### 2.3 What each engine writes

The row on ENST00000338785 in the default format. Ensembl VEP 116.2:

```
onto_exon	21:37480637	-	ENSG00000157540	ENST00000338785	Transcript	frameshift_variant,splice_region_variant	2078	328	110	V/X	Gtt/tt	-	IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1;NearestExonJB=ENSE00003584551|1|start|189
```

vep-rs:

```
onto_exon	21:37480637	-	ENSG00000157540	ENST00000338785	Transcript	frameshift_variant,splice_region_variant	2078	328	110	V/X	Gtt/tt	-	IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1
```

Ensembl VEP writes an item on 26 rows: `ENSE00003584551|1|start|189` on 22, and one each of `ENSE00004532916|1|start|1703`, `ENSE00003828235|1|start|1012`, `ENSE00003675684|1|start|189` and `ENSE00001803088|1|start|32`, every exon one that starts at 37480638. The 27th row, the `intron_variant` row of ENST00000839928 on the reverse strand, has `SHIFT_LENGTH=0` and no value in either engine. vep-rs writes no item on any of the 27 rows. Every other column of the 27 rows agrees, except on the row of ENST00000398956, where the two engines' start-codon tests give different consequence terms, a difference [the start-codon class](shift-start-codon-test.md) documents.

### 2.4 Why Ensembl VEP writes that

The consequence tests fill the memoised exon set at the shifted span, 37480638, with the exons that start there. The plugin then measures each exon's distance from the record's span: 1, from 37480637 to the exon's start, which lies one base outside the exon the item names. The item places the variant inside an exon and at distance 1 from the exon's start, from the wrong side.

### 2.5 Why vep-rs's row is the correct one

At the record's span the variant overlaps no exon, and the plugin's rule writes nothing for such a row without `intronic=1`, as both engines do for this record without the flag. Given the same exons through `gff3=`, Ensembl VEP 116.2 under the shift flags writes no item on any of the 27 rows, as vep-rs does. A value measured from the record's span and naming the exons of another span describes no variant of the record.

### 2.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --fasta <reference> --plugin NearestExonJB --shift_3prime 1 --shift_length
```

vep-rs writes no `NearestExonJB` on any row; Ensembl VEP 116.2 writes `NearestExonJB=ENSE00003584551|1|start|189` on the row of ENST00000338785 and an item on 25 more. Without the two shift flags neither engine writes a value.

### 2.7 How the comparator treats this class

As in [1.7](#17-how-the-comparator-treats-this-class): the comparator sets the cell of either shape aside where a comparison case names the class.
