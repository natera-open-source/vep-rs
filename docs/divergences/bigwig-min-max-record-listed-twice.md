# A bigWig record listed twice under `min` and `max`

**The defect in one sentence.** For a bigWig source of `type=overlap` with `summary_stats` of `min` and `max` alone and no other matching option, Ensembl VEP 116.2 runs its per-record walk twice over every variant of at most 65,536 bases. Each such variant lists every record it matches twice. vep-rs lists each record once, and the minimum and the maximum are the same in both engines.

**Why this is a defect.** A source's list holds the records a variant matches, and an interval the file holds once is one record. The same variant lists each record once when the source asks a further statistic, as in `min%max%count`, or none at all. The reader's own description of its span rule says that a small variant falls through to the base class with its behaviour unchanged. The minimum and the maximum are the same either way; the list alone is doubled, and under a record limit each pass lists up to the limit and closes with `...`.

**What vep-rs does instead.** vep-rs annotates every variant once on that path: a variant of at most 65,536 bases by the per-record walk, a longer one by the span rule.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/AnnotationSource/`):

- `File/BigWig.pm`, lines 181-182, `annotate_InputBuffer`: a source the span rule does not take goes to the base class once. Lines 239-249, `_can_use_zoom_summary`, decide it: `min` and `max` alone, `type` overlap, and no `overlap_cutoff`, `distance`, `reciprocal`, `same_type` or `PC` field.
- Lines 184-192: a variant of at most `zoom_min_span` bases, 65,536 (lines 86-88), is set aside for the base class, and a longer one is annotated by the span rule.
- Lines 195-204: the buffer is pointed at the variants set aside, and `SUPER::annotate_InputBuffer` is called on it inside an `eval`, lines 197-201. Line 202 calls it on the same buffer again.
- `File.pm`, lines 436-440, `annotate_VariationFeature`: a pass pushes each matching record onto the variant's list while its own count is below `num_records`, so the second pass appends the first pass's records again.
- `File/BigWig.pm`, lines 157-167: the method's description, by which small variants "fall through to the base class, so their behaviour is unchanged".

**What correcting it would change.** Without the second call each record would be listed once, as under every other bigWig source. `<short_name>_min` and `<short_name>_max` would not move.

**The class it produces.** `bigwig_min_max_record_listed_twice`: the `<short_name>` list of a bigWig source the span rule takes, for a variant of at most 65,536 bases. The reference's list is the candidate's list written twice. The class is documented only: no comparator codes it, since no committed flag set reaches it.

**How often the shape occurs.** Tier 3: the shape is a bigWig source of `type=overlap` whose `summary_stats` name `min`, `max` or both and nothing else, with no other matching option. Of the 392 committed flag sets, 4 name a bigWig source and none names one of that shape: the one committed `overlap` source with statistics asks `min%max%mean%count%sum`.

## 1. Two intervals listed four times

**In short.** The bigWig `scores.bw` of the `GRCh38-custom` golden corpus holds one-base intervals at `21:5028997` and `21:5028999`, with the values `0.001` and `0.002`. Under `--custom file=scores.bw,short_name=cons,format=bigwig,type=overlap,summary_stats=min%max`, the record `mnv`, `21 5028997 TCA>AGT`, overlaps both. Ensembl VEP 116.2 prints `cons=0.00100000004749745,0.0020000000949949,0.00100000004749745,0.0020000000949949`; vep-rs prints `cons=0.00100000004749745,0.0020000000949949`. Both print the same `cons_min` and `cons_max`.

### 1.1 What a reader needs to know first

A `--custom` source lists, under its short name, the records a variant matches, and a bigWig record is one interval of the file, listed by its value. `summary_stats` adds one item per statistic asked, `<short_name>_min` and `<short_name>_max` here, computed over the matched values ([bigWig sources](../output/custom-bigwig.md#bigwig-sources)).

Ensembl VEP reads a bigWig source by its per-record walk, which visits each interval of the variant's window and lists the ones that match. A source of `type=overlap` that asks `min` and `max` alone, with no other matching option, also has a span rule. A variant of more than 65,536 bases takes the rule, which reads the extrema over the whole span; a shorter variant is left to the walk.

### 1.2 The record

The bigWig is `tests/golden/116/GRCh38-custom/custom/scores.bw`. It names its chromosome `21`, and its first three intervals, written here as bedGraph lines, are:

```
21	5028996	5028997	0.001
21	5028998	5028999	0.002
21	5029000	5029001	0.042
```

The input is one VCF record, a three-base substitution over positions `5028997` to `5028999`, which the first two intervals overlap:

```
21	5028997	mnv	TCA	AGT	.	.	.
```

Both engines read the corpus's cache and take one flag:

```
--custom file=scores.bw,short_name=cons,format=bigwig,type=overlap,summary_stats=min%max
```

### 1.3 What each engine writes

The row in the default format. Ensembl VEP 116.2:

```
mnv	21:5028997-5028999	AGT	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;cons=0.00100000004749745,0.0020000000949949,0.00100000004749745,0.0020000000949949;cons_max=0.0020000000949949;cons_min=0.00100000004749745
```

vep-rs:

```
mnv	21:5028997-5028999	AGT	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;cons=0.00100000004749745,0.0020000000949949;cons_max=0.0020000000949949;cons_min=0.00100000004749745
```

The rows agree in every other column and in both statistics. With `summary_stats=min%max%count`, or with no `summary_stats`, both engines list `cons=0.00100000004749745,0.0020000000949949`. With `summary_stats=max` Ensembl VEP lists the four values again. With `num_records=1` added to the `min%max` source, Ensembl VEP lists `cons=0.00100000004749745,...,0.00100000004749745,...` and vep-rs `cons=0.00100000004749745,...`.

### 1.4 Why Ensembl VEP writes that

The bigWig reader overrides the base class's annotation loop for a source its span rule takes. It sends each variant of more than 65,536 bases to the rule and collects the others. It then points the buffer at the collected variants and calls the base class's loop on them inside an `eval`. The next statement calls the same loop on the same buffer once more, and only then is the full buffer restored.

Each call walks the variant's window and pushes every matching record onto the variant's list, so the list holds each record once per call. Each call also counts its own records from zero, so a record limit is applied to each pass by itself. The minimum and the maximum are taken over the values of both passes, which repeat one another, so neither moves.

Source (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/AnnotationSource/`):

- `File/BigWig.pm`, lines 184-192: the variants are divided at 65,536 bases.
- `File/BigWig.pm`, lines 197-202: the two calls of the base class's loop.
- `File.pm`, lines 332-341 and 436-440: a pass counts from zero and pushes each matching record.

### 1.5 Why vep-rs's row is the correct one

The file holds one interval at position `5028997` and one at `5028999`, and the record overlaps each once, so its list has two records. Ensembl VEP lists those two for the same record whenever the source is outside the span rule, with a third statistic or with none, and its reader describes the short variant's behaviour under the rule as unchanged. Asking for the minimum and the maximum should add two items to the row; it should not change which records the row lists. vep-rs lists the same records whatever statistics are asked.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 \
  --custom file=<corpus>/custom/scores.bw,short_name=cons,format=bigwig,type=overlap,summary_stats=min%max
```

vep-rs writes `cons=0.00100000004749745,0.0020000000949949`. Ensembl VEP 116.2 takes the same `--custom` flag with the cache options of [Reproducing these examples](reproducing.md) and writes the two values twice.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/bigwig_min_max_record_listed_twice.toml` runs the record of 1.2 on both engines and records the `cons` item of its row. The class covers the list of a bigWig source the span rule takes, for a variant of at most 65,536 bases. No committed flag set names such a source, so no comparator codes the class.
