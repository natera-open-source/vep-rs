# A structural insertion point just outside a RefSeq model whose edits renumber its end base

**The defect in one sentence.** On a RefSeq model whose `_rna_edit` attributes renumber its outermost genomic base, Ensembl VEP 116.3 writes `feature_elongation` for a symbolic insertion point lying between that base and the base beyond the model, beside the `upstream_gene_variant` or `downstream_gene_variant` that vep-rs writes alone, with `IMPACT=MODIFIER`.

**Why this is a defect.** The point's outer neighbour is not a base of the model: no exon holds it, and the row's own `upstream_gene_variant` or `downstream_gene_variant` says the variant begins beyond the model's end. The row prints no cDNA position, where a `feature_elongation` row of an insertion point inside a model prints one. On a model without edits the same point reads `upstream_gene_variant` or `downstream_gene_variant` alone, and so it does on an edited model unless its first base is numbered above 1 or its last base below the sum of its exon lengths. `feature_elongation` has HIGH impact, so the row's `IMPACT` changes with the term.

**What vep-rs does instead.** vep-rs reads an insertion point by its two neighbouring bases. A point with a neighbour outside the model's span lies outside the model, and its row carries `upstream_gene_variant` or `downstream_gene_variant` alone, with `IMPACT=MODIFIER`, on an edited model as on any other.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation`, `ensembl` and `ensembl-io` at release/116, the modules the Ensembl VEP 116.3 image holds):

- `Bio/EnsEMBL/IO/Parser/BaseVCF4.pm`, lines 262-264, `get_start`: a record with a symbolic allele starts one base after its `POS`. Lines 302-309, `get_end`: its end is counted from `SVLEN` where the record has one, and is its `END` otherwise. A symbolic record with no `SVLEN` whose `END` is its `POS` therefore has its start one above its end.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 337-348, `feature_elongation`: the term needs `within_cdna`, `complete_within_feature` and an insertion or a gain. `Bio/EnsEMBL/Variation/Utils/Constants.pm`, lines 810-813: its filter does not ask for a variant within the feature.
- `VariationEffect.pm`, lines 163-166, `complete_within_feature`: the variant's start is compared with the feature's start and its end with the feature's end. An insertion point has its start one above its end, so a point straddling either end of the feature passes.
- `VariationEffect.pm`, lines 670-694, `within_cdna`: a mapped coordinate counts when its end is above 0 and its start is at most the transcript's length (line 680). `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, line 488, `cdna_coords`: the coordinates are the variant's start and end, mapped through the model's mapper. `Bio/EnsEMBL/Mapper.pm`, lines 485-580, `map_insert`: the point's one mapped neighbour, at cDNA position c, becomes the coordinate (c, c - 1) at the 5' end and (c + 1, c) at the 3' end.
- `Bio/EnsEMBL/TranscriptMapper.pm`, lines 207-211, `_load_mapper`: an edit that adds bases moves the cDNA numbers after it up, and one that removes bases moves them down. `Bio/EnsEMBL/Transcript.pm`, lines 1666-1682, `length`: the sum of the exon lengths, without the edits on a cached model, whose `edits_enabled` the cache build drops (`Bio/EnsEMBL/VEP/Pipeline/DumpVEP/Dumper/Core.pm`, line 164).
- The length that line 680 compares with is also behind [`refseq_edited_model_tail_term`](refseq-edited-model-tail-term.md), where a base numbered above it fails the test and the row of a sequence allele loses its term.

**What correcting it would change.** `feature_elongation` would ask for both neighbours of an insertion point inside the feature, as the `within_feature` test does. A point just outside an edited model's end would carry `upstream_gene_variant` or `downstream_gene_variant` alone, with `IMPACT=MODIFIER`, as it does on a model without edits.

**The class it produces.** `refseq_edited_model_edge_insertion_elongation`: a `Transcript` row of a symbolic insertion point whose Consequence is `feature_elongation` with `upstream_gene_variant` or `downstream_gene_variant` in Ensembl VEP and the upstream or downstream term alone in vep-rs, with the `IMPACT` that follows the term. No comparison tool sets a row of the shape aside: a disagreement of this shape is reported as a difference and charged to vep-rs. The concordance datasets are annotated on the Ensembl caches, which hold no `_rna_edit`, and no corpus record is such a point at such a model's end, so the entry is documented on this page and the comparator carries no rule for it.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a symbolic insertion written as a point, with `END` at its `POS` and no `SVLEN`, between the end base of a RefSeq model and the base beyond it. The model is one whose edits number its first base above 1 or its last base below the sum of its exon lengths. The two Ensembl caches the concordance datasets run on hold no `_rna_edit` attribute, so no model of theirs is renumbered: 0 of 669,092 on GRCh38 and 0 of 238,823 on GRCh37. No concordance dataset is run with `--refseq` or `--merged`, so the datasets hold 0 rows of the shape.

The 32 files of the two structural-variant sets hold 10,212 symbolic insertions and no point. Each of their 3,527 symbolic records with `END` at `POS` carries an `SVLEN` that makes it a span.

The committed caches hold 160 models with an `_rna_edit` attribute, each counted once in every cache that holds it, in ten corpora: eight RefSeq and merged corpora, `AaegL5-aedes` and `AgamP4-anopheles`. On 14 of them the last base is numbered below the sum of the exon lengths: 6 in `AaegL5-aedes` and 8 on chromosome 21, `NM_130444.3` and `NM_030582.4` in `GRCh37-refseq`, `GRCh38-refseq` and `GRCh38-merged`, and `NM_198689.3` in `GRCh37-refseq` and `GRCh37-merged`. On none is the first base numbered above 1.

The inputs of those ten corpora hold 1,508 records, 22 of them symbolic and none a point. Of the 8,113 lines of the committed Ensembl VEP outputs that carry `feature_elongation`, none carries it beside an upstream or downstream term.

A set of 39 hand-written points on the RefSeq GRCh38 cache of release 116 gives 455 rows in each engine under `--refseq`. On 282 of them the point lies just outside the row's own model. On 9 of the 282 the model's first base is numbered above 1 or its last base below the sum of its exon lengths, and there Ensembl VEP writes `feature_elongation` on all 9 and vep-rs on none. On the other 273 neither engine writes it, and the engines' other 446 rows are equal in every column. Eight of the 9 lie at a model's 5' end and one at a 3' end.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. An insertion point just beyond the first base of RYBP's model `NM_012234.7`

**In short.** `3:72446623 <INS>` with `END=72446623`, on the RefSeq GRCh38 cache under `--refseq`, is the point between the model's first genomic base and the base beyond it. Ensembl VEP writes `feature_elongation,upstream_gene_variant` with `IMPACT=HIGH`; vep-rs writes `upstream_gene_variant` with `IMPACT=MODIFIER`.

### 1.1 What a reader needs to know first

A RefSeq model can hold bases the genome lacks. The cache records them as `_rna_edit` attributes and numbers the model's cDNA with them. `NM_012234.7`, on the reverse strand, has 213 such bases before its first genomic base, so that base, `3:72446623`, is cDNA 214 and not cDNA 1.

A symbolic insertion whose `END` equals its `POS`, with no `SVLEN`, is a point between two neighbouring bases. Both engines print this one as `3:72446623-72446624`.

A cDNA position counts a model's bases from its first base, in the direction of transcription, so the first genomic base of a reverse-strand model is its highest genomic position. The model's mapper pairs each unbroken stretch of cDNA positions with the genomic positions it stands on. A genomic base that no pair holds has no cDNA position.

`feature_elongation` is the term for a structural variant that lengthens a feature. `upstream_gene_variant` and `downstream_gene_variant` are the terms for a variant outside a model, before its first base or after its last.

### 1.2 The record

```
3	72446623	.	T	<INS>	.	.	SVTYPE=INS;END=72446623
```

The model spans `3:72374597-72446623`. The point's lower neighbour, `3:72446623`, is the model's first base; its higher neighbour, `3:72446624`, is outside the model.

The cache holds, for the model, the `_rna_edit` attribute `1 0` followed by 213 bases: the model holds 213 bases before the position that would be its first without them. Its first mapper pair sets cDNA 214-341 on `3:72446496-72446623`, and its four exons sum to 4,445 bases.

### 1.3 What each engine writes

The row of `NM_012234.7`, by column:

| Field | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| Location | `3:72446623-72446624` | `3:72446623-72446624` |
| Consequence | `feature_elongation,upstream_gene_variant` | `upstream_gene_variant` |
| cDNA_position | `-` | `-` |
| IMPACT | `HIGH` | `MODIFIER` |
| DISTANCE | `0` | `0` |

The two rows differ in `Consequence` and `IMPACT` and in nothing else.

Eight models of the RefSeq GRCh38 primary assembly hold bases before their first genomic base: `NM_012234.7`, `NM_053005.5`, `NM_182705.3`, `NM_001394353.1`, `NM_001394354.1`, `NM_001394355.1`, `NR_131204.1` and `NR_165399.1`. The point just outside the 5' end reads the same two ways on each.

At a 3' end the same happens where the edits remove bases, so that the model's last base is numbered below the sum of its exon lengths. The last base of `NR_185482.1`, `9:62474853`, is numbered 7871 against 7,878 exon bases: for the point just beyond it Ensembl VEP writes `feature_elongation,downstream_gene_variant` and vep-rs `downstream_gene_variant`.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP starts a symbolic record one base after its `POS` and ends it at its `END`, so this record has the start 72446624 and the end 72446623.

The point's start, 72446624, is at least the model's start, and its end, 72446623, is at most the model's end, so `complete_within_feature` holds. The mapper holds `3:72446623` at cDNA 214 and nothing at `3:72446624`, so the point's coordinate is (214, 213). Its end is above 0 and its start is below the model's 4,445 exon bases, so `within_cdna` holds. On a model without the added bases the coordinate would be (1, 0), and the test would refuse it.

In numbers, 72446624 is at least 72374597 and 72446623 is at most 72446623 for the first test, and 213 is above 0 and 214 is at most 4,445 for the second. Both tests hold and the record is an insertion, so `feature_elongation` is written.

`upstream_gene_variant` is written beside it because the point's start, 72446624, is above the model's end, 72446623 (`VariationEffect.pm`, lines 413-414, `_after_end`). On a reverse-strand model that side is upstream (lines 420-422).

### 1.5 Why vep-rs's row is the correct one

`3:72446624` is not a base of the model, and the 213 added bases are not the genome's bases beyond the model. Compared base by base, they agree at 66 of 213 positions, and the added base next to the model differs from the base at `3:72446624`. The point is next to the model, as `upstream_gene_variant` on Ensembl VEP's own row says.

The coordinate (214, 213) places the point between the model's first genomic base and cDNA 213, the last of the added bases. The record's own coordinates place it between `3:72446623` and `3:72446624`, a base that no exon and no mapper pair holds. vep-rs asks whether both neighbours lie inside the model's span, `3:72374597-72446623`. 72446624 is above 72446623, so the point is outside, and the row carries `upstream_gene_variant` alone.

### 1.6 Try it yourself

The commands of [Reproducing these examples](reproducing.md#reproducing-these-examples), with `--refseq` and the RefSeq GRCh38 cache, print this row; `grep NM_012234.7` selects it.

vep-rs runs on the conversion of that cache ([Cache conversion](../cache-conversion.md)). With `--use_given_ref` in place of a reference FASTA, each engine prints its row of 1.3. Section 2 gives a record that runs on a cache the repository holds.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`tools/vep-formats-compare`) has no class of this shape, and no committed corpus holds a row of it. A disagreement of the shape is reported as a difference in `Consequence` and `IMPACT` and charged to vep-rs.

Rows beside the shape are no part of this entry. A point inside a model gets `feature_elongation` beside a cDNA position from both engines. A point just outside a model whose end base is not renumbered gets the upstream or downstream term alone from both.

## 2. An insertion point just beyond the last base of two of COL18A1's models

**In short.** `21:45513720 <INS>` with `END=45513720`, on the cache of the committed `GRCh38-refseq` corpus under `--refseq`, is the point between the last base of three COL18A1 models and the base beyond it. Two of them, `NM_130444.3` and `NM_030582.4`, carry an edit that removes nine bases, and on their rows Ensembl VEP writes `feature_elongation,downstream_gene_variant` with `IMPACT=HIGH`. vep-rs writes `downstream_gene_variant` with `IMPACT=MODIFIER` there, as both engines do on the row of the third model, `NM_001379500.1`, which carries no edit.

### 2.1 What a reader needs to know first

An `_rna_edit` attribute can also say that a model lacks bases the genome holds. It then names the cDNA positions the missing bases would have and gives no bases: `4092 4100` says that nine bases are missing, 4100 - 4092 + 1 = 9. Every base after them is numbered nine lower than it would be without the edit.

Ensembl VEP's `within_cdna` test compares a coordinate with the model's length. The length a cached model returns is the sum of its exon lengths without the edits. On a model that lacks nine bases, that sum is nine above the number of the model's last base.

### 2.2 The record

The record is hand-written on the cache and reference of the `GRCh38-refseq` golden corpus (`tests/golden/116/GRCh38-refseq/`):

```
21	45513720	col18a1_end_point	A	<INS>	.	.	SVTYPE=INS;END=45513720
```

Three forward-strand models of COL18A1 in that cache end at `21:45513720`, the last base of the last exon they share, `21:45512188-45513720`. `NM_130444.3` and `NM_030582.4` span `21:45455522-45513720`, and `NM_001379500.1` spans `21:45405165-45513720`.

| Model | `_rna_edit` | Exon bases | Number of the last base |
| --- | --- | --- | --- |
| `NM_130444.3` | `4092 4100` | 6,596 | 6,587 |
| `NM_030582.4` | `3387 3395` | 5,891 | 5,882 |
| `NM_001379500.1` | none | 5,408 | 5,408 |

On the two edited models the missing bases are `21:45504526-45504534`, nine genomic bases inside an exon: 6,596 - 9 = 6,587 and 5,891 - 9 = 5,882. The last of them lies 45513720 - 45504534 = 9,186 bases before the models' last base.

The point's lower neighbour, `21:45513720`, is the last base of all three models; its higher neighbour, `21:45513721`, is outside them. Both engines run `--refseq` with the corpus's reference FASTA.

### 2.3 What each engine writes

The three models' rows in the default format. Ensembl VEP 116.3:

```
col18a1_end_point	21:45513720-45513721	insertion	80781	NM_001379500.1	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1
col18a1_end_point	21:45513720-45513721	insertion	80781	NM_030582.4	Transcript	feature_elongation,downstream_gene_variant	-	-	-	-	-	-	IMPACT=HIGH;DISTANCE=0;STRAND=1;BAM_EDIT=OK
col18a1_end_point	21:45513720-45513721	insertion	80781	NM_130444.3	Transcript	feature_elongation,downstream_gene_variant	-	-	-	-	-	-	IMPACT=HIGH;DISTANCE=0;STRAND=1;BAM_EDIT=OK
```

vep-rs:

```
col18a1_end_point	21:45513720-45513721	insertion	80781	NM_001379500.1	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1
col18a1_end_point	21:45513720-45513721	insertion	80781	NM_030582.4	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1;BAM_EDIT=OK
col18a1_end_point	21:45513720-45513721	insertion	80781	NM_130444.3	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1;BAM_EDIT=OK
```

The rows of the two edited models differ in `Consequence` and `IMPACT` and in nothing else. `BAM_EDIT=OK` marks a model that the cache's alignment edited ([RefSeq and merged cache columns](../output/refseq-and-merged.md)).

The record gives 28 rows in each engine, and the other 26 are equal: the row of `NM_001379500.1` and 25 rows of SLC19A1's models on the other strand, which the point lies inside. On 22 of those both engines write `feature_elongation,3_prime_UTR_variant` beside a cDNA position, such as `2847-2848` on `NM_001205206.4`. On the other 3 both write `intron_variant`.

### 2.4 Why Ensembl VEP writes that

The record has the start `POS` + 1 = 45513721 and the end `END` = 45513720. `complete_within_feature` holds on all three models: 45513721 is at least each model's start, 45455522 or 45405165, and 45513720 is at most its end, 45513720.

The mapper holds `21:45513720` at the number of the model's last base and nothing at `21:45513721`. On `NM_130444.3` the point's coordinate is therefore (6588, 6587). Its end is above 0 and its start, 6588, is at most the 6,596 exon bases, so `within_cdna` holds and the term is written. On `NM_030582.4` the coordinate is (5883, 5882), and 5883 is at most 5,891.

On `NM_001379500.1` the coordinate is (5409, 5408). Its start is above the model's 5,408 exon bases, `within_cdna` refuses it, and Ensembl VEP writes `downstream_gene_variant` alone.

### 2.5 Why vep-rs's row is the correct one

`21:45513721` is not a base of any of the three models: their last exon ends at `21:45513720`. The point lies after the models' last base, as `downstream_gene_variant` on all three of Ensembl VEP's own rows says.

The position Ensembl VEP's test accepts on `NM_130444.3`, 6588, is not a position of the model. The cache holds the model's sequence with 6,587 bases, the 6,596 exon bases less the nine it lacks, and its last base is number 6,587. The bound the test applies, 6,596, counts nine bases the model does not have.

Ensembl VEP's term turns on the edit and not on the record. The three models share their last exon and their last base, and the record touches none of the nine missing bases, which lie 9,186 bases away. On the model without the edit Ensembl VEP writes what vep-rs writes on all three.

### 2.6 Try it yourself

From the repository's root, the command as vep-rs takes it:

```
printf '##fileformat=VCFv4.2\n' > point.vcf
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >> point.vcf
printf '21\t45513720\tcol18a1_end_point\tA\t<INS>\t.\t.\tSVTYPE=INS;END=45513720\n' >> point.vcf
vep -i point.vcf -o out.txt --offline --cache \
  --dir_cache tests/golden/116/GRCh38-refseq/cache_root --cache_version 116 \
  --species homo_sapiens --assembly GRCh38 --refseq \
  --fasta tests/golden/116/GRCh38-refseq/reference.fa.gz --force_overwrite --no_stats --quiet
grep -E 'NM_(130444|030582|001379500)' out.txt
```

vep-rs prints the three rows of 2.3 with `downstream_gene_variant` alone. Ensembl VEP 116.3 takes the same record and `--refseq` with its RefSeq GRCh38 cache of release 116 and the cache options of [Reproducing these examples](reproducing.md). It prints `feature_elongation,downstream_gene_variant` on the two edited models' rows.

The entry's reproduction, `scripts/divergences/reproductions/refseq_edited_model_edge_insertion_elongation.toml`, runs this record on both engines and holds the three rows to their recorded values.

### 2.7 How the comparator treats this class

Nothing is set aside here either. The two edited models' rows are reported as differences in `Consequence` and `IMPACT` and charged to vep-rs. The row of `NM_001379500.1` and the 25 rows of SLC19A1's models are equal in both engines and are no part of this entry.
