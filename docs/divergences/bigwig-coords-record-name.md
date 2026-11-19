# The name of a bigWig record under `coords`

**The defect in one sentence.** Under `coords`, Ensembl VEP 116.2's per-record walk names a bigWig record by the file's chromosome name less a leading `chr` and by the interval's 0-based start with its end. A one-base interval at position `5029000` of `chr21` prints as `21:5028999-5029000`. vep-rs prints the same interval as `chr21:5029000-5029000`, the file's chromosome name with the interval's 1-based start and end.

**Why this is a defect.** The name has the form of a region, `<chromosome>:<start>-<end>`, and a region Ensembl VEP prints is 1-based and inclusive: the `Location` column is, and so is a BED record's name under `coords`. A bigWig stores an interval 0-based and half-open, and the walk prints the stored start beside the end. The start printed is one below the interval's first base, so a one-base interval reads as two bases and two adjacent intervals read as sharing one. The reader's matching code adds one to that start before it tests the interval, and its span rule prints the start plus one under the file's chromosome name. One source so names its records in two forms.

**What vep-rs does instead.** vep-rs names every record by the file's chromosome name with the interval's 1-based start and end, the form Ensembl VEP's span rule prints, whatever the variant's length. The values, the matching and the statistics agree; the name alone differs.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/AnnotationSource/`; `ensembl-io` at release/116, `Bio/EnsEMBL/IO/`):

- `File/BigWig.pm`, lines 416-428, `_get_record_name`: under `coords` the name is `sprintf('%s:%i-%i', ...)` over the parser's `get_seqname`, `get_start` and `get_end`.
- `BigFileParser.pm`, lines 105-108: the file's chromosomes are keyed by their names with `s/^chr//` applied, each key mapped to the file's own name.
- `Parser/BigWig.pm`, lines 60 and 68, `seek`: an interval is stored under the key its chromosome was looked up by, with its start less one. Lines 107-110, `get_start`: the stored start plus one, which is the interval's 0-based start.
- `Parser/Wig.pm`, lines 169-173, `get_seqname`: the stored name with `s/^chr//` applied, which leaves the key as it is.
- `File.pm`, lines 529-532: for a bigWig the matching code adds one to that start, under the comment `BigWig format is 0-based half-open`.
- `File/BigWig.pm`, lines 391-393, `_collect_span_records`: the span rule names a record `sprintf('%s:%i-%i', $seq_id, $i->start + 1, $i->end)`, the file's own chromosome name and the 1-based start.

**What correcting it would change.** Printing the start plus one in `_get_record_name`, under the name the span rule prints, would give every record one name whatever the variant's length. No value, match or statistic would move.

**The class it produces.** `bigwig_coords_record_name`: the `<short_name>` list of a bigWig source read with `coords` on, for a variant the per-record walk reads. The reference names each record with a start one below the candidate's, and without the leading `chr` the candidate prints when the file's chromosome names carry one. The class is documented only: no comparator codes it, since no committed flag set reaches it.

**How often the shape occurs.** Tier 3: the shape is a bigWig source whose `coords` option is on, and every record such a source lists on the walk is named alike. Of the 392 committed flag sets, 4 name a bigWig source and none passes `coords` to one.

## 1. A one-base interval named as a two-base region

**In short.** The bigWig `named.bw` holds three one-base intervals on `chr21`, at positions `5029000` to `5029002`. Under `--custom file=named.bw,short_name=named,format=bigwig,type=overlap,coords=1`, the record `snv`, `21 5029000 T>G`, overlaps the first. Ensembl VEP 116.2 prints `named=21:5028999-5029000`; vep-rs prints `named=chr21:5029000-5029000`. The record `pair`, a two-base substitution at `21:5029001-5029002`, lists `21:5029000-5029001,21:5029001-5029002` in Ensembl VEP and `chr21:5029001-5029001,chr21:5029002-5029002` in vep-rs.

### 1.1 What a reader needs to know first

A `--custom` source lists, under its short name, the records a variant matches. A bigWig record is one interval of the file. It is listed by its value, or by its region when the source's `coords` option is on ([bigWig sources](../output/custom-bigwig.md#bigwig-sources)). A bigWig stores an interval as a chromosome name, a 0-based start and an end, so the one-base interval that covers position `5029000` is stored as `5028999` to `5029000`.

Ensembl VEP reads a bigWig source in one of two ways. Its per-record walk visits each interval of the variant's window and lists the ones that match. Its span rule serves a variant of more than 65,536 bases under `type=overlap` with `summary_stats` of `min` and `max` alone and no other matching option, and lists the intervals of the span's first 32,768 bases. Every other variant, and every variant of another source, takes the walk.

### 1.2 The record

The bigWig is `scripts/divergences/reproductions/inputs/bigwig_coords_record_name/named.bw`. It names its chromosome `chr21` and holds three intervals, written here as bedGraph lines:

```
chr21	5028999	5029000	0.25
chr21	5029000	5029001	0.5
chr21	5029001	5029002	0.75
```

The input is two VCF records, a substitution at the first interval and a two-base substitution over the other two:

```
21	5029000	snv	T	G	.	.	.
21	5029001	pair	AC	GT	.	.	.
```

Both engines read the cache of the `GRCh38-custom` golden corpus (`tests/golden/116/GRCh38-custom/`) and take one flag:

```
--custom file=named.bw,short_name=named,format=bigwig,type=overlap,coords=1
```

### 1.3 What each engine writes

The two rows in the default format. Ensembl VEP 116.2:

```
snv	21:5029000	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;named=21:5028999-5029000
pair	21:5029001-5029002	GT	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;named=21:5029000-5029001,21:5029001-5029002
```

vep-rs:

```
snv	21:5029000	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;named=chr21:5029000-5029000
pair	21:5029001-5029002	GT	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;named=chr21:5029001-5029001,chr21:5029002-5029002
```

The rows agree in every other column. Without `coords` both engines list the values, `named=0.25` for `snv` and `named=0.5,0.75` for `pair`.

### 1.4 Why Ensembl VEP writes that

The walk builds a record's name from three getters of the bigWig parser. `get_seqname` returns the name the interval is stored under. That is the key the parser looks a chromosome up by, the file's name with a leading `chr` removed, so `chr21` prints as `21`. `get_start` returns the interval's start as the file stores it, 0-based: the parser's `seek` subtracts one from the start its library reports, and `get_start` adds it back. `get_end` returns the end, which for a half-open interval is its last base. The interval that covers position `5029000` is therefore named `21:5028999-5029000`.

The same reader treats that start as 0-based where it matches. Before it tests an interval against the variant it adds one to `get_start`, under the comment `BigWig format is 0-based half-open`. The name is built from the start as stored.

Source (`ensembl-vep` at release/116.2, `ensembl-io` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/BigWig.pm`, lines 416-428: `_get_record_name`.
- `Bio/EnsEMBL/IO/BigFileParser.pm`, lines 105-108: the chromosome keys.
- `Bio/EnsEMBL/IO/Parser/Wig.pm`, lines 169-173: `get_seqname`.
- `Bio/EnsEMBL/IO/Parser/BigWig.pm`, lines 60, 68 and 107-110: the stored key, the stored start and `get_start`.
- `Bio/EnsEMBL/VEP/AnnotationSource/File.pm`, lines 529-532: the start raised by one for matching.

### 1.5 Why vep-rs's row is the correct one

The interval `chr21 5028999 5029000` covers one base, position `5029000`, the base the record `snv` changes. Ensembl VEP writes a region 1-based and inclusive, as in the `Location` column of the same row, so this interval's region is `5029000-5029000`. The name `21:5028999-5029000` spans two bases, and the interval does not cover the first. The next interval's name, `21:5029000-5029001`, spans position `5029000` as well, so two intervals that share no base read as overlapping.

Ensembl VEP's own span rule prints the form vep-rs prints. With `summary_stats=min%max` added to the source, the deletion `21 5028990 N <DEL>` with `END=5100000` takes the span rule, and both engines list `chr21:5029000-5029000,chr21:5029001-5029001,chr21:5029002-5029002,...` for it. In that run Ensembl VEP still names the first interval `21:5028999-5029000` for the record `snv`, twice over, which is the documented [`bigwig_min_max_record_listed_twice`](bigwig-min-max-record-listed-twice.md) difference. vep-rs gives an interval one name on both paths.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 \
  --custom file=named.bw,short_name=named,format=bigwig,type=overlap,coords=1
```

vep-rs writes `named=chr21:5029000-5029000` for `snv`. Ensembl VEP 116.2 takes the same `--custom` flag with the cache options of [Reproducing these examples](reproducing.md) and writes `named=21:5028999-5029000`.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/bigwig_coords_record_name.toml` runs the two records of 1.2 on both engines and records the `named` item of each row. The class covers the name of a bigWig record under `coords` for a variant the per-record walk reads. No committed flag set passes `coords` to a bigWig source, so no comparator codes the class.
