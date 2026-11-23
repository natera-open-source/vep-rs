# The NearestExonJB value of a later allele or record at the same span

**The defect in one sentence.** Ensembl VEP's NearestExonJB plugin keeps one finished item string per transcript and feature span for the whole run, under a key that names no allele, so a later row at that span takes the first row's value or emptiness, whatever its own exons and terms. vep-rs values every row from its own span and terms, so the ALT order and the records before a row change nothing.

**Why this is a defect.** The plugin's own rules are allele-specific: the exon set it reads for a `--minimal` allele is the one its trimmed span overlaps, and under `intronic=1` it reports the flanking exons only on a row whose terms include `intron_variant`. The cache then serves every later row with the same transcript and span, whichever allele or record it belongs to. For the same allele of the same record, Ensembl VEP 116.3 writes a value in one ALT order and nothing in the other, and writes for a record the value of the record before it.

**What vep-rs does instead.** vep-rs values every row from its own span and terms: a `--minimal` allele's exon set and distances from its trimmed span ([a documented class](minimal-plugin-sees-the-rejoined-record.md)), and the `intronic=1` flanking exons on the rows whose terms carry `intron_variant`. The ALT order and the records before a row change nothing. The plugin's parameters and fields are described under [NearestExonJB](../plugins/transcript-model.md#nearestexonjb).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `NearestExonJB.pm`, line 179: the key is `sprintf("%s:%s-%i-%i", $trv->transcript_stable_id, $vf->{chr} || $vf->seq_region_name, $vf->{start}, $vf->{end})`, the transcript and the feature's span and no allele. Line 181 computes only on a miss, line 274 stores the joined items or `undef`, line 276 serves the stored value. `_cache` is written at line 274 and read at lines 181 and 276 alone, so it spans the run.
- `NearestExonJB.pm`, line 184: the exon set is `$trv->_overlapped_exons`, memoised on the TranscriptVariation at `BaseTranscriptVariation.pm` lines 851-855 from the feature's span at its first call. Lines 190-194 gate the `intronic` branch on the allele's own terms, `$consequences =~ /intron/`.
- Without `--minimal`, `OutputFactory.pm` lines 389-398 run the plugin once per allele of a record in the order the alleles were added, the alternate alleles in ALT order (`BaseVariationFeatureOverlap.pm` lines 317-318), so the first alternate allele's string is handed to the second.
- Under `--minimal`, `InputBuffer.pm` lines 638-718 split a record into one feature per trimmed allele (line 665 trims, lines 674-676 give it the trimmed string and span). `OutputFactory.pm` lines 2298 and 2312 run each split feature on its own span, lines 2300-2302 restore the first feature's span to the record's, and lines 2323-2325 move the merged feature's TranscriptVariations onto the first feature under the key `<trimmed allele string>_<transcript key>` (`BaseVariationFeature.pm` lines 310-321). The output pass runs them in `sort keys` order (`VariationFeature.pm` line 666), so a merged allele whose key sorts before the first feature's `ENST` key runs first and fills the record-span key.

**What correcting it would change.** The plugin would compute every row from its own exon set and terms, and a row's value would not depend on the ALT order or on the records before it.

**The class it produces.** `nearest_exon_jb_shared_span_cache`: a `NearestExonJB` value that differs between the engines on a row of a multi-allelic record under `--minimal` or `intronic=1`, or on a record at a span an earlier record of the run shares. No comparison tool sets a row of the shape aside: the committed corpus that runs the plugin carries no such record, so a disagreement of this shape is reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape the record alone shows is a multi-allelic record whose sequence alleles trim to different spans, 1,894 of the cuts' 15,646,481 records (121.0 per million), in 27 corpora. The transcript decides whether the alleles' exon sets or `intron_variant` terms differ there. Of the committed flag sets, one runs the plugin with `--minimal` or `intronic=1`. Two records at one span in one run, the second shape, are counted by no predicate.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A `--minimal` record whose alleles trim to the two sides of an exon boundary

**In short.** `21 36897131 cc_ct_gc CC CT,GC` on GRCh38 sits on the last exon base and the first intron base of ENSE00001282176, an exon of ENST00000336648 on the reverse strand. Under `--minimal`, `CT` trims to a `T` at 36897132, in the intron, and `GC` to a `G` at 36897131, in the exon. On ENST00000336648, Ensembl VEP 116.3 writes `ENSE00001282176|1|end|272` on both alleles' rows in the ALT order `CT,GC` and on neither in the order `GC,CT`. vep-rs writes `ENSE00001282176|0|end|272` on the `G` row and nothing on the `T` row, in both orders.

### 1.1 What a reader needs to know first

NearestExonJB reports, for every exon a variant overlaps, the exon's stable id, the distance to its nearer junction, that junction and its length. A row with no overlapped exon gets no value unless `intronic=1` is set. Under `--minimal` Ensembl VEP trims each allele of a multi-allelic record, annotates the trimmed alleles one at a time and rejoins their rows under the record's span. The `Location` column prints the record's span, each allele's exon set is the one its trimmed span overlaps, and Ensembl VEP measures the item's distance from the span the row prints, where vep-rs measures it from the allele's trimmed span ([a documented class](minimal-plugin-sees-the-rejoined-record.md)).

### 1.2 The record

A multi-allelic substitution against the Ensembl 116 GRCh38 cache. The reference is `CC` at 36897131-36897132; 36897131 is the last base of ENSE00001282176 in genomic order and 36897132 the first base of the intron after it.

```
21	36897131	cc_ct_gc	CC	CT,GC	.	.	.
```

Both engines run `--minimal --plugin NearestExonJB`. The two ALT orders are the records `cc_ct_gc` (`CT,GC`) and `cc_gc_ct` (`GC,CT`).

### 1.3 What each engine writes

The rows on ENST00000336648 in the default format. Ensembl VEP 116.3 on `CT,GC`:

```
cc_ct_gc	21:36897131-36897132	G	ENSG00000159267	ENST00000336648	Transcript	missense_variant,splice_region_variant	?-1651	?-1180	?-394	E/Q	Gaa/Caa	-	ALLELE_NUM=2;IMPACT=MODERATE;STRAND=-1;MINIMISED=1;NearestExonJB=ENSE00001282176|1|end|272
cc_ct_gc	21:36897131-36897132	T	ENSG00000159267	ENST00000336648	Transcript	splice_acceptor_variant	-	-	-	-	-	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1;MINIMISED=1;NearestExonJB=ENSE00001282176|1|end|272
```

Ensembl VEP 116.3 on `GC,CT`:

```
cc_gc_ct	21:36897131-36897132	T	ENSG00000159267	ENST00000336648	Transcript	splice_acceptor_variant	-	-	-	-	-	-	ALLELE_NUM=2;IMPACT=HIGH;STRAND=-1;MINIMISED=1
cc_gc_ct	21:36897131-36897132	G	ENSG00000159267	ENST00000336648	Transcript	missense_variant,splice_region_variant	?-1651	?-1180	?-394	E/Q	Gaa/Caa	-	ALLELE_NUM=1;IMPACT=MODERATE;STRAND=-1;MINIMISED=1
```

vep-rs on `CT,GC`:

```
cc_ct_gc	21:36897131-36897132	G	ENSG00000159267	ENST00000336648	Transcript	missense_variant,splice_region_variant	?-1651	?-1180	?-394	E/Q	Gaa/Caa	-	ALLELE_NUM=2;IMPACT=MODERATE;STRAND=-1;MINIMISED=1;NearestExonJB=ENSE00001282176|0|end|272
cc_ct_gc	21:36897131-36897132	T	ENSG00000159267	ENST00000336648	Transcript	splice_acceptor_variant	-	-	-	-	-	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1;MINIMISED=1
```

vep-rs on `GC,CT`:

```
cc_gc_ct	21:36897131-36897132	T	ENSG00000159267	ENST00000336648	Transcript	splice_acceptor_variant	-	-	-	-	-	-	ALLELE_NUM=2;IMPACT=HIGH;STRAND=-1;MINIMISED=1
cc_gc_ct	21:36897131-36897132	G	ENSG00000159267	ENST00000336648	Transcript	missense_variant,splice_region_variant	?-1651	?-1180	?-394	E/Q	Gaa/Caa	-	ALLELE_NUM=1;IMPACT=MODERATE;STRAND=-1;MINIMISED=1;NearestExonJB=ENSE00001282176|0|end|272
```

The record writes 34 rows in each engine, 17 per allele, with the same consequence terms on every row. The `NearestExonJB` cell by ALT order:

| ALT order | Allele | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `CT,GC` | `G` | `ENSE00001282176\|1\|end\|272` | `ENSE00001282176\|0\|end\|272` |
| `CT,GC` | `T` | `ENSE00001282176\|1\|end\|272` | none |
| `GC,CT` | `G` | none | `ENSE00001282176\|0\|end\|272` |
| `GC,CT` | `T` | none | none |

The other 16 rows of each allele follow the same pattern with their own transcripts' exon ids. Ensembl VEP writes a value on the 17 `T` rows in the first order and on none of the 34 rows in the second; vep-rs writes the 17 `G` rows alike in both, each at distance 0 from its exon's end.

### 1.4 Why Ensembl VEP writes that

The split leaves two features: the first ALT's, which keeps the record's span after the rejoin, and the merged allele's, whose TranscriptVariations move onto the first feature under keys prefixed with the trimmed allele string. In the order `CT,GC` the merged allele is `G`, keyed `C/G_ENST00000336648_…`, which sorts before the first feature's `ENST00000336648_…`, so the `G` row runs first. Its memoised exon set is the exon at 36897131 and its distances are measured from the record's span, 36897131-36897132, so it writes `ENSE00001282176|1|end|272` and stores that string under `ENST00000336648:21-36897131-36897132`. The `T` row runs next with the same transcript and span, line 181 finds the key, and line 276 serves the `G` row's string.

In the order `GC,CT` the merged allele is `T`, keyed `C/T_…`, and runs first. Its memoised exon set, read over 36897132 alone, is empty, so line 274 stores `undef` under the record-span key. The `G` row then finds the key and gets nothing, although its own span overlaps the exon.

### 1.5 Why vep-rs's row is the correct one

Each allele's exon set is a fact of its own trimmed span: `G` changes the exon's last base and `T` the intron's first. Ensembl VEP agrees when each allele is given alone. The bi-allelic `CC GC` under `--minimal` writes a value on every `G` row in both engines, `ENSE00001282176|0|end|272` on ENST00000336648's, zero because a bi-allelic record is trimmed in place and measured from the trimmed base; `CC CT` writes nothing on any `T` row in both.

Without `--minimal` the untrimmed alleles both overlap the exon and the engines write 34 identical rows, each with a distance of 1 to its exon's end. vep-rs's `G` rows inside the multi-allelic record are the rows `CC GC` writes alone: the exon set and the distance of the allele's own span, in either ALT order.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --minimal --plugin NearestExonJB
```

vep-rs writes `NearestExonJB=ENSE00001282176|0|end|272` on the `G` row of ENST00000336648 and no value on the `T` row, whichever ALT comes first. Ensembl VEP 116.3 writes it on both rows when `CT` comes first and on neither when `GC` does. Swapping the two ALT alleles in the record shows the difference.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`tools/vep-formats-compare`) has no class of this shape, and the golden harness compares the `GRCh38-plugin-computed` corpus's `NearestExonJB` cells as multisets, which the corpus's 47 single-allele records never put to the test. A disagreement of the shape is reported as a difference and charged to vep-rs. A `NearestExonJB` value on a row of a `MINIMISED` split record belongs to [`minimal_plugin_sees_the_rejoined_record`](minimal-plugin-sees-the-rejoined-record.md), the class the two mechanisms share on that shape.

## 2. An `intronic=1` record whose alleles differ in `intron_variant`, and two records at one span

**In short.** `21 37480827 gtaa_gtag_gcaa GTAA GTAG,GCAA` on GRCh38 lies in the first four bases of a DYRK1A intron, after exon ENSE00003584551 of ENST00000338785, which ends at 37480826. `GTAG` changes the fourth intron base and is `splice_donor_region_variant,intron_variant`; `GCAA` changes the second and is `splice_donor_variant` alone. Under `--plugin NearestExonJB,intronic=1`, without `--minimal`, the untrimmed alleles share the span 37480827-37480830. Ensembl VEP 116.3 writes the two flanking exons on both alleles' rows in the order `GTAG,GCAA` and on neither in the order `GCAA,GTAG`; vep-rs writes them on the `GTAG` rows alone in both orders. Two bi-allelic records at the span behave the same way in Ensembl VEP, the second taking the first's value; vep-rs values each record as when it is alone.

### 2.1 What a reader needs to know first

Under `intronic=1` a row with no overlapped exon reports the two exons flanking its intron, each with its own distance, when its terms include `intron_variant`. A `splice_donor_variant` row carries no `intron_variant`, so the plugin's own gate reports nothing on it. Without `--minimal` every allele of a record is annotated on the record's untrimmed span, so the alleles share the key and differ only in their terms.

### 2.2 The record

```
21	37480827	gtaa_gtag_gcaa	GTAA	GTAG,GCAA	.	.	.
```

The reference is `GTAA` at 37480827-37480830. Both engines run `--plugin NearestExonJB,intronic=1`. The two ALT orders are the records `gtaa_gtag_gcaa` and `gtaa_gcaa_gtag`; the shape with two records is one file holding `21 37480827 first_gcaa GTAA GCAA` and then `21 37480827 second_gtag GTAA GTAG`, and the same two records in the other order.

### 2.3 What each engine writes

The rows on ENST00000338785 in the default format. Ensembl VEP 116.3 on `GTAG,GCAA`:

```
gtaa_gtag_gcaa	21:37480827-37480830	GTAG	ENSG00000157540	ENST00000338785	Transcript	splice_donor_region_variant,intron_variant	-	-	-	-	-	-	IMPACT=LOW;STRAND=1;NearestExonJB=ENSE00001033299|5640|start|148,ENSE00003584551|4|end|189
gtaa_gtag_gcaa	21:37480827-37480830	GCAA	ENSG00000157540	ENST00000338785	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1;NearestExonJB=ENSE00001033299|5640|start|148,ENSE00003584551|4|end|189
```

Ensembl VEP 116.3 on `GCAA,GTAG`:

```
gtaa_gcaa_gtag	21:37480827-37480830	GCAA	ENSG00000157540	ENST00000338785	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
gtaa_gcaa_gtag	21:37480827-37480830	GTAG	ENSG00000157540	ENST00000338785	Transcript	splice_donor_region_variant,intron_variant	-	-	-	-	-	-	IMPACT=LOW;STRAND=1
```

vep-rs on `GTAG,GCAA`:

```
gtaa_gtag_gcaa	21:37480827-37480830	GTAG	ENSG00000157540	ENST00000338785	Transcript	splice_donor_region_variant,intron_variant	-	-	-	-	-	-	IMPACT=LOW;STRAND=1;NearestExonJB=ENSE00003584551|4|end|189,ENSE00001033299|5640|start|148
gtaa_gtag_gcaa	21:37480827-37480830	GCAA	ENSG00000157540	ENST00000338785	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
```

vep-rs on `GCAA,GTAG`:

```
gtaa_gcaa_gtag	21:37480827-37480830	GCAA	ENSG00000157540	ENST00000338785	Transcript	splice_donor_variant	-	-	-	-	-	-	IMPACT=HIGH;STRAND=1
gtaa_gcaa_gtag	21:37480827-37480830	GTAG	ENSG00000157540	ENST00000338785	Transcript	splice_donor_region_variant,intron_variant	-	-	-	-	-	-	IMPACT=LOW;STRAND=1;NearestExonJB=ENSE00003584551|4|end|189,ENSE00001033299|5640|start|148
```

The record writes 54 rows in each engine, 27 per allele. Given alone, `GTAG` gets a value on 26 of its 27 rows (the flanking exons on its 24 `intron_variant` rows, an overlapped exon on 2 exonic rows) and `GCAA` on 3 of its 27 (one `intron_variant` row, the same 2 exonic rows). The multi-allelic record's valued rows per allele:

| ALT order | Allele | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `GTAG,GCAA` | `GTAG` | 26 of 27 | 26 of 27 |
| `GTAG,GCAA` | `GCAA` | 26 of 27 | 3 of 27 |
| `GCAA,GTAG` | `GCAA` | 3 of 27 | 3 of 27 |
| `GCAA,GTAG` | `GTAG` | 3 of 27 | 26 of 27 |

The two items of a row are the same in both engines; their order is Ensembl VEP's hash order, which the comparison tools read as a multiset. With the two records in one file, Ensembl VEP values 3 of the second record's 27 rows when `first_gcaa` comes first and 26 when `first_gtag` does. vep-rs values each record as when it is alone, 26 of 27 for `GTAG` and 3 of 27 for `GCAA`, in either order.

### 2.4 Why Ensembl VEP writes that

Without `--minimal` the plugin runs once per allele in ALT order on one feature, the record's span. In the order `GTAG,GCAA` the `GTAG` row runs first: no exon overlaps the span, its terms include `intron_variant`, and lines 190-217 collect the two flanking exons, stored under `ENST00000338785:21-37480827-37480830`. The `GCAA` row has the same transcript and span, so line 276 serves it the stored string although line 194 would have forbidden the items on its own terms. In the order `GCAA,GTAG` the `GCAA` row runs first, line 194 finds no `intron` in `splice_donor_variant`, line 274 stores `undef`, and the `GTAG` row finds the key and gets nothing.

Two bi-allelic records at the span behave the same way across records: the cache is never cleared, so the second record's row on each transcript takes the first record's stored value.

### 2.5 Why vep-rs's row is the correct one

The flanking exons are a fact of the row's own terms: `GTAG` is an `intron_variant` and `GCAA` a `splice_donor_variant` on ENST00000338785, and the plugin's rule reports the exons on the first and nothing on the second. Given alone, `GTAA GTAG` writes a value on 26 rows in both engines (the two flanking exons on 24) and `GTAA GCAA` on 3. vep-rs's rows for each allele inside the multi-allelic record equal Ensembl VEP's rows for that allele alone, as item sets, in both orders. A record's rows in vep-rs do not change with the records before it.

### 2.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --assembly GRCh38 \
  --plugin NearestExonJB,intronic=1
```

vep-rs writes `NearestExonJB=ENSE00003584551|4|end|189,ENSE00001033299|5640|start|148` on the `GTAG` row of ENST00000338785 and no value on the `GCAA` row, whichever ALT comes first. Ensembl VEP 116.3 writes the two items on both rows when `GTAG` comes first and on neither when `GCAA` does. Splitting the record into two bi-allelic records in one file shows the second taking the first's value.

### 2.7 How the comparator treats this class

As in [1.7](#17-how-the-comparator-treats-this-class): nothing is set aside, and a disagreement of either shape is reported as a difference and charged to vep-rs.
