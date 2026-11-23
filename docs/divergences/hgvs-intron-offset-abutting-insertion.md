# The HGVSIntronOffset offsets of an insertion directly beside a transcript

**The defect in one sentence.** For an insertion directly beside a transcript's first or last base, Ensembl VEP 116.3's HGVSIntronOffset plugin writes `HGVS_IntronStartOffset` and `HGVS_IntronEndOffset` on a row that lies outside the transcript and carries no HGVSc. The row reads `upstream_gene_variant` or `downstream_gene_variant` with `DISTANCE=0`, and the two values are the offsets of a notation Ensembl VEP never prints. vep-rs writes the two fields beside an HGVSc and nowhere else: a row without an HGVSc carries neither.

**Why this is a defect.** The two fields are the intron parts of an HGVSc's start and end coordinates. Such an insertion falls between two neighbouring positions of which exactly one is a base of the transcript, its first or its last, so the record writes the inserted bases outside the transcript. The row's own term and `DISTANCE=0` say so, and Ensembl VEP gives the row no HGVSc. A row without a notation has no notation coordinates, and so no intron offsets of them. The values belong to another change: on most rows a substitution or a replacement of one base of the transcript, a base the record leaves as it is.

**What vep-rs does instead.** vep-rs writes the two fields beside an HGVSc and nowhere else, as [the plugin's page](../plugins/transcript-model.md#hgvsintronoffset) states. It gives such a row an HGVSc only where the HGVS 3' rule places the inserted bases on the transcript, which [`hgvs_added_bases_written_outside_the_transcript`](hgvs-added-bases-written-outside-the-transcript.md) describes. Every other such row carries no HGVSc and neither field.

**Where it is in Ensembl VEP's source** (`VEP_plugins`, `ensembl-vep` and `ensembl-variation` at release/116):

- `HGVSIntronOffset.pm`, line 62: the plugin asks the allele for its transcript notation on every Transcript row. Ensembl VEP's own output asks only on a row inside the feature (`Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1730 and 1760-1763), which is why the row has no HGVSc.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2689-2736, `_var2transcript_slice_coords`: an insertion's start is one above its end. Directly before the first base it maps to the slice coordinates (1, 0); directly after the last base of a slice of L bases, to (L+1, L). Line 2729 returns nothing only when both lie before the first base or both past the last, and lines 2732-2733 clamp the pair onto the edge base: (1, 1) or (L, L).
- The same file, `hgvs_transcript`: lines 1446 and 1477-1478 add the length k of the insertion's 3' shift to both coordinates, and lines 1474-1483 describe the inserted bases against the one base there. That base is the edge base where k is 0, and base 1+k where the shift carries an insertion before the first base into the transcript. Lines 1543 and 1545 store the intron offsets of the notation's coordinates, which the plugin reads at its lines 64-65.
- The same file, lines 2188-2248, `_clip_alleles`: the bases the notation's two alleles share at their start and at their end are stripped, and the notation's start or end moves by as many positions. An inserted allele that ends in the edge base and does not begin with it loses that base at lines 2231-2238, which move the notation's end one position down. At a forward-strand transcript's first base that is from slice coordinate 1 to 0, the position before the transcript.
- The same file, lines 2802-2826, `_get_cDNA_position`: the position after a reverse-strand transcript's last base and the position before a forward-strand transcript's first base both lie before the first exon in genomic order. The loop reaches its intron branch at index 0, where line 2802 reads the exon at index -1, the transcript's other end, as the exon before. The edge exon is the nearer of the two, so lines 2822-2826 build the coordinate off it: "last base +1" on the reverse strand and "first base -1" on the forward strand.

**What correcting it would change.** The plugin would write its fields only where the row carries an HGVSc, and a row outside the transcript without one would carry neither field.

**The class it produces.** `hgvs_intron_offset_abutting_insertion`: the two fields on Ensembl VEP's side alone, on an insertion's row to which both engines give the one term `upstream_gene_variant` or `downstream_gene_variant`, `DISTANCE=0`, no HGVSc and no `HGVS_OFFSET`. The pair is 0 and 0, or 1 and 0 on the downstream row of a reverse-strand transcript, or 0 and -1 on the upstream row of a forward-strand transcript; a row that prints no `STRAND` is read as either strand. The format comparator holds the rule. Any other pair, a value on vep-rs's side, a row with a term inside the transcript and the fields on vep-rs's side alone stay differences.

The same clamp gives [`hgvs_shifted_insertion_clamped_end`](hgvs-shifted-insertion-clamped-end.md) its notation: under `--shift_3prime 1` a shifted insertion's row moves into the transcript and prints it. Here the row stays outside, and only the plugin asks for a notation.

**How often the shape occurs.** Tier 2, the plugin being needed. On the genome-wide ClinVar GRCh37 file 344 rows are an insertion directly beside a transcript. Ensembl VEP writes the fields on 142 of them, from 135 records and 107 transcripts: 0 and 0 on 125, 1 and 0 on 17. On the GRCh38 file it writes them on 429 of 734 such rows, from 260 records and 257 transcripts: 0 and 0 on 403, 1 and 0 on 26. Neither file holds the pair 0 and -1, which section 4 shows on a gnomAD record.

On the other 202 and 305 rows the plugin's call returns no notation. On 35 of the 142 and 167 of the 429, all upstream rows of forward-strand transcripts, the 3' shift carries the insertion into the transcript and the notation sits that many bases in. No committed corpus carries the shape.

vep-rs writes neither field on 78 of the 142, from 75 records and 56 transcripts, and on 160 of the 429, from 101 records and 111 transcripts: the rows of this class. Ensembl VEP's pair there is 0 and 0 on 72 and 158, 1 and 0 on 6 and 2. On 56 and 105 the insertion lies before a first base and has no 3' shift; on 22 and 55 it lies past a last base and is no whole copy of the transcript's last bases.

On the other 64 and 269, those 35 and 167 among them, vep-rs writes an HGVSc with the two fields of that notation, as it does on 63 and 128 of the rows where the plugin's call returns nothing. [`hgvs_added_bases_written_outside_the_transcript`](hgvs-added-bases-written-outside-the-transcript.md) describes those rows.

## 1. A `T` inserted just before the first base of ENST00000480892

**In short.** `2 27670443 1365957 G GT`, a ClinVar record on GRCh37, inserts a `T` between 27670443 and 27670444. ENST00000480892, a reverse-strand transcript, runs from 27668808 to 27670443: its first base is 27670443, and the inserted base lies just outside, on its 5' side. Both engines write `upstream_gene_variant`, `DISTANCE=0` and no HGVSc. Ensembl VEP 116.3 adds `HGVS_IntronStartOffset=0` and `HGVS_IntronEndOffset=0`. vep-rs writes neither field.

### 1.1 What a reader needs to know first

The HGVSIntronOffset plugin writes two fields beside a transcript notation, the intron parts of its start and end coordinates: `45` for `c.123+45`, `-10` for `c.124-10`, `0` for an exonic coordinate. Its own description says it "returns HGVS intron start and end offsets" and is "To be used with --hgvs option".

An insertion sits between two neighbouring bases. Ensembl VEP gives it a start one above its end, here 27670444 and 27670443, and a distance of 0 to a transcript that one of the two bases belongs to.

ENST00000480892 (ENSG00000138002) has two exons on the reverse strand of chromosome 2: 1,636 bases of genome, 27668808 to 27670443, and 477 spliced bases. Its first exon is 27670382 to 27670443, so n.1 is the base at 27670443.

### 1.2 The record

```
2	27670443	1365957	G	GT	.	.	.
```

Both engines run `--hgvs --plugin HGVSIntronOffset` on GRCh37 with a reference FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.3, default format:

```
1365957	2:27670443-27670444	T	ENSG00000138002	ENST00000480892	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=-1;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

vep-rs:

```
1365957	2:27670443-27670444	T	ENSG00000138002	ENST00000480892	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=-1
```

Every other column is equal. Ensembl VEP writes the same pair on the row when the plugin runs without `--hgvs`. On ENST00000260570, a longer transcript of the gene, the same record is the duplication `c.4597dup`, and both engines write that HGVSc with the two fields, 0 and 0.

### 1.4 Why Ensembl VEP writes that

The transcript is on the reverse strand, so line 2708 maps the insertion's end and start, 27670443 and 27670444, to 27670443 less each, plus 1: the slice coordinates 1 and 0. Line 2729 asks whether both are below 1. One is not, so the method goes on, and line 2733 clamps the 0 to 1.

The record has no 3' shift, so lines 1477-1478 add nothing. The call then describes the span (1, 1), the transcript's first base, against the inserted base: a substitution of n.1. Both of its coordinates are exonic, so lines 1543 and 1545 store 0 and 0, and the plugin, having received a notation at line 62, writes the pair. The output's own code never asks, the row not being inside the feature, and prints no HGVSc.

### 1.5 Why vep-rs's row is the correct one

The record inserts a base between 27670443 and 27670444. Of the two, only 27670443 is a base of the transcript, and it is the first: the inserted base lies before n.1, outside the transcript, and differs from n.1, so the HGVS 3' rule leaves it there. Both engines say so with `upstream_gene_variant` and `DISTANCE=0`, and neither prints an HGVSc. The offsets of an HGVSc's coordinates do not exist where there is no HGVSc. The pair Ensembl VEP writes belongs to a substitution of n.1, and the record leaves that base as it is.

### 1.6 Try it yourself

Save the record of 1.2 as `record.vcf` and annotate it against a GRCh37 cache with the reference FASTA. vep-rs's command is below; Ensembl VEP takes the same record with its own cache and the FASTA under `--offline --cache`:

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh37 --hgvs --plugin HGVSIntronOffset
```

vep-rs's row for ENST00000480892 ends `DISTANCE=0;STRAND=-1`. Ensembl VEP 116.3 adds `HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0` to it, with `--hgvs` and without.

### 1.7 How the comparator treats this class

The format comparator's rule (`tools/vep-formats-compare`, `intron_offset_abutting_insertion`) accepts both fields of the row, `hgvs_intronstartoffset` and `hgvs_intronendoffset` in JSON: Ensembl VEP's side alone carries them, the pair is 0 and 0, the row is an insertion's, and both engines give it `upstream_gene_variant`, `DISTANCE=0` and no HGVSc. The same pair on a row inside a transcript and a value on vep-rs's side stay differences. So does any pair but the three of the class, and so does 1 and 0 or 0 and -1 on a row other than the one the class names for it.

## 2. Four bases inserted just past the last base of ENST00000465204

**In short.** `2 48033890 801710 C CTTCT`, a ClinVar record on GRCh37, inserts four bases between 48033890 and 48033891. ENST00000465204, on the reverse strand, runs from 48033891 to 48038800 over five exons and 3,158 spliced bases: 48033891 is its last base, and the inserted bases lie just past it. Both engines write `downstream_gene_variant`, `DISTANCE=0` and no HGVSc. Ensembl VEP writes `HGVS_IntronStartOffset=1` and `HGVS_IntronEndOffset=0`; vep-rs writes neither field.

The insertion maps to the slice coordinates (L+1, L) and is clamped onto the last base. The call's notation is an insertion between n.3158 and the position after it, 48033890 on the genome, which lies before the transcript's first exon in genomic order. The loop of `_get_cDNA_position` reaches its intron branch at index 0, line 2802 reads the exon at index -1, and lines 2822-2826 build "3158+1". The start offset 1 is the intron part of that coordinate, at a position where no intron is.

On a forward-strand transcript the position after the last base falls through the loop and the call returns nothing, so every row carrying 1 and 0 is the downstream row of a reverse-strand transcript. ClinVar holds six insertions at this position, of one to five bases. On this transcript's row Ensembl VEP writes 1 and 0 for four of them, 0 and 0 for one and no field for one: what the call returns follows the inserted bases, where the insertion's place is the same.

The other five insertions repeat the transcript's last bases. vep-rs writes each as the duplication of those bases with its two fields, 0 and 0, a row of [`hgvs_added_bases_written_outside_the_transcript`](hgvs-added-bases-written-outside-the-transcript.md). This record's four bases repeat only the transcript's last two, a partial copy, so the row keeps no notation and stays this class's.

## 3. A `G` whose 3' shift runs seven bases into ENST00000478223

**In short.** `1 2179434 3895388 C CG`, a ClinVar record on GRCh37, inserts a `G` just before 2179435, the first base of ENST00000478223 on the forward strand. Both engines write `upstream_gene_variant` and `DISTANCE=0`, and Ensembl VEP no HGVSc. Ensembl VEP adds the two fields, 0 and 0; vep-rs writes `HGVSc=ENST00000478223.2:n.7dup`, `HGVS_OFFSET=7` and the two fields, 0 and 0. The row is not this class's: [`hgvs_added_bases_written_outside_the_transcript`](hgvs-added-bases-written-outside-the-transcript.md) holds it.

The insertion's 3' shift is seven bases long. The pair (1, 0) is clamped to (1, 1), lines 1477-1478 add 7 to both, and the call's notation is a substitution of n.8, an exonic base of the 76-base first exon. It is the notation [`hgvs_shifted_insertion_clamped_end`](hgvs-shifted-insertion-clamped-end.md) describes, which a row prints once `--shift_3prime 1` has moved it into the transcript. Under `--hgvs` alone the row stays upstream, Ensembl VEP prints no HGVSc, and the plugin's two fields are all that shows of it. vep-rs describes the inserted base where the HGVS 3' rule places it, as the duplication of n.7, and its two fields are that notation's.

## 4. An inserted allele that ends in the first base of ENST00001091162

**In short.** `21 10521561 . A AAG`, a gnomAD v4.1 record on GRCh38, inserts `AG` between 10521561 and 10521562. ENST00001091162, a forward-strand transcript, runs from 10521562 to 10605709: its first base is 10521562, and the inserted bases lie just outside, on its 5' side. Both engines write `upstream_gene_variant`, `DISTANCE=0` and no HGVSc. Ensembl VEP 116.3 adds `HGVS_IntronStartOffset=0` and `HGVS_IntronEndOffset=-1`. vep-rs writes neither field.

### 4.1 The record

```
21	10521561	.	A	AAG	.	.	.
```

ENST00001091162 (ENSG00000274391) has its first exon at 10521562 to 10521694 and its start codon 461 bases in, so its first base is c.-461. That base is a `G`, the last base of the inserted allele. Both engines run `--hgvs` and the HGVSIntronOffset plugin on GRCh38 with a reference FASTA, in a run that loads further plugins.

### 4.2 What each engine writes

Ensembl VEP 116.3, default format:

```
21_10521561_A/AAG	21:10521561-10521562	AG	ENSG00000274391	ENST00001091162	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1;AA=-;HGVS_IntronEndOffset=-1;HGVS_IntronStartOffset=0;LoFtool=0.0451;NearestExonJB=ENSE00004715058|0|start|133;TSSDistance=1;gnomAD_100x_cov=0.0000;gnomAD_10x_cov=0.0014;gnomAD_15x_cov=0.0004;gnomAD_1x_cov=0.0901;gnomAD_20x_cov=0.0001;gnomAD_25x_cov=0.0000;gnomAD_30x_cov=0.0000;gnomAD_50x_cov=0.0000;gnomAD_5x_cov=0.0123;gnomAD_mean_cov=0.2242;gnomAD_median_approx_cov=0.0000;gnomAD_total_DP_cov=163894.0000;pHaplo=0.139986193997075;pLI_gene_value=0.00;pTriplo=0.0764328113393634
```

vep-rs:

```
21_10521561_A/AAG	21:10521561-10521562	AG	ENSG00000274391	ENST00001091162	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1;AA=-;LoFtool=0.0451;NearestExonJB=ENSE00004715058|0|start|133;TSSDistance=1;gnomAD_100x_cov=0.0000;gnomAD_10x_cov=0.0014;gnomAD_15x_cov=0.0004;gnomAD_1x_cov=0.0901;gnomAD_20x_cov=0.0001;gnomAD_25x_cov=0.0000;gnomAD_30x_cov=0.0000;gnomAD_50x_cov=0.0000;gnomAD_5x_cov=0.0123;gnomAD_mean_cov=0.2242;gnomAD_median_approx_cov=0.0000;gnomAD_total_DP_cov=163894.0000;pHaplo=0.139986193997075;pLI_gene_value=0.00;pTriplo=0.0764328113393634
```

The two offsets are the only difference: the fields of the run's other plugins stand on both rows with the same values. On ENST00001091161, a transcript of the same gene that starts two bases further on, the same record's row reads `DISTANCE=2` and neither engine writes a field.

### 4.3 Why Ensembl VEP writes that

Line 2712 maps the insertion's start and end, 10521562 and 10521561, to the slice coordinates 1 and 0, and line 2733 clamps the 0 to 1. The record has no 3' shift. The call describes the span (1, 1), the transcript's first base `G`, against the inserted `AG`. `_clip_alleles` strips the `G` the two share at their end and moves the notation's end from 1 to 0 (lines 2231-2238), and lines 2271-2284 leave an insertion of `A` with start 1 and end 0.

`_get_cDNA_position` reads the start as c.-461, an exonic base. The end, slice coordinate 0, is 10521561 on the genome, which lies before the first exon. The loop reaches its intron branch at index 0, line 2802 reads the exon at index -1, and lines 2822-2826 build "1-1", the first base with an intron offset of -1, which the start-codon shift of lines 2862-2867 prints as "-461-1". Lines 1543 and 1545 store 0 and -1, and the plugin writes the pair.

### 4.4 Why vep-rs's row is the correct one

The record inserts two bases between 10521561 and 10521562. Only 10521562 is a base of the transcript, and it is the first: the inserted bases lie before c.-461, outside the transcript, and the record leaves that base as it is. The first of the two differs from c.-461, so the HGVS 3' rule leaves them there, and neither engine writes a notation of this transcript for them: both rows read `upstream_gene_variant`, `DISTANCE=0` and no HGVSc. The end offset -1 is counted at 10521561, where the transcript has no intron: its first exon starts at its first base. A row without a notation has no notation coordinates, and so no intron offsets.

### 4.5 How the comparator treats the pair

The rule accepts 0 and -1 on the upstream row of a forward-strand transcript under the class's other conditions: Ensembl VEP's side alone carries the two fields, the row is an insertion's, and both engines give it `upstream_gene_variant`, `DISTANCE=0`, no HGVSc and no `HGVS_OFFSET`. The same pair on a downstream row, on a reverse-strand transcript's row and on a row inside a transcript stays a difference.
