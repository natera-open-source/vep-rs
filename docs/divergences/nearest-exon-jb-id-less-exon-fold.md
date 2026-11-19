# The NearestExonJB item of two exons without stable ids

**The defect in one sentence.** Ensembl VEP's NearestExonJB plugin files each exon's distance and length under the exon's stable id, so on a transcript whose exons carry no id every exon of a row shares one empty key and one item is printed from the parts of different exons. vep-rs writes one item per reported exon, the exon's own distance, junction and length, with the empty id where the exon has none.

**Why this is a defect.** An item, `ExonID|distance|start/end|length`, describes one exon: its nearer junction, how far away it is and how long it is. Under a shared key the length printed is the last exon's filed, and the distance is the nearest exon's or, under `intronic=1`, the one Perl's hash order reads last. A deletion across an intron then prints one exon's junction with the other's length, an item no exon has, and two exons equally near print one item instead of two. Under `intronic=1` the two flanking exons print one item: the first exon's distance with the second's length, or the second's item alone. The transcripts of the same gene whose exons carry ids print each exon's own item.

**What vep-rs does instead.** vep-rs writes one item per reported exon, the exon's own distance, junction and length, with the empty id where the exon has none, so a transcript without exon ids prints the items its id-carrying twin prints, less the ids. The plugin's parameters and fields are described under [NearestExonJB](../plugins/transcript-model.md#nearestexonjb).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-variation` at release/116):

- `NearestExonJB.pm`, line 220: `$exon_id` is the exon's `stable_id`. Lines 228-237 file the exon's distance and side under `$dists{$exon_id}{$distance}` and its length under `$dists{$exon_id}{len}`. An exon without an id gives the key the empty string, so every id-less exon of a row writes into one entry, and a later write of `len`, or of an equal distance, replaces the one before.
- `NearestExonJB.pm`, lines 267-271: the default branch prints one item per key, the key's `$min` entry and its `len`, so the item carries the nearest exon's junction and the length of the exon filed last.
- `NearestExonJB.pm`, lines 244-263: the `intronic=1` branch walks the key's inner hash, `for my $internal_key (keys %{$inner_hash})`, and the last distance key read sets the printed distance and side, so the item follows Perl's hash order, which changes between runs and from row to row.
- `BaseTranscriptVariation.pm`, line 865: the overlapped exons come from the transcript's interval tree in the tree's own order, so which of two exons is filed last differs between transcripts of the same exons. The flanking exons of the `intronic=1` branch come from `_sorted_exons` (`NearestExonJB.pm`, lines 197-215), the downstream exon last.

**What correcting it would change.** The plugin would key the exons by the exon itself, or print an item for each exon it measured, so a transcript whose exons carry no id would print the items its id-carrying twin prints.

**The class it produces.** `nearest_exon_jb_id_less_exon_fold`: a `NearestExonJB` value whose items carry no exon id and that differs between the engines: Ensembl VEP writes one item where vep-rs writes two or more, or where the variant begins in one exon and ends in another. The format comparator sets the value aside where a comparison case names the class, under the test of section 1.7. The committed corpus that runs the plugin annotates from an Ensembl cache, whose exons all carry ids, so no committed row has the shape. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a transcript with two or more exons lacking a stable id, 675 of the 77,601 transcripts of the committed caches, in 2 corpora. Those are the merged and RefSeq caches, whose aligned RefSeq-source models (biotype `aligned_transcript`, written under `--all_refseq`) carry no exon ids. A record shows it when it overlaps two such exons, or lies in such a transcript's intron under `intronic=1`. Of the committed flag sets, none runs the plugin on a merged or RefSeq cache; the two that run it annotate from an Ensembl cache.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion across an intron of a transcript whose exons carry no id

**In short.** `21 10578400 idless_two_exons ATAATG…GTTT A` on GRCh38 deletes 10578401-10578530: the last 18 bases of the exon at 10578337-10578418 (82 bases), the 98-base intron and the first 14 bases of the exon at 10578517-10578605 (89 bases) of `compmerge.2.Heart.chr21`, a model whose exons carry no stable id. Both engines run `--merged --all_refseq --use_given_ref --plugin NearestExonJB`. Ensembl VEP 116.2 writes `|64|start|89`, the first exon's distance and junction with the second exon's length; vep-rs writes `|64|start|82`. On the Ensembl transcripts of the same gene, whose exons carry ids, both engines write `ENSE00001659869|64|start|82`.

### 1.1 What a reader needs to know first

NearestExonJB reports, for each exon a variant overlaps, the exon's stable id, the distance to its nearer junction, that junction and the exon's length; only the exons at the smallest distance are reported, all of them when two are equally near. The merged and RefSeq caches carry, under `--all_refseq`, RefSeq-source models aligned to the genome whose exons have no stable id. Both engines write such an exon's item with the id empty, `|75|end|89` for a row overlapping one exon alone.

### 1.2 The record

A deletion against the Ensembl 116 GRCh38 merged cache. The reference allele is the anchor base at 10578400 and the 130 deleted bases; the second record extends the deletion by 11 bases, to 10578541, so that both exons end up 64 bases from their nearer junction.

```
21	10578400	idless_two_exons	ATAATGTCCCCACTCTACAGTAAGTTTTATGCTAAATTGACTATCAGAAGAGGGCTAAATATATTGGGCTTGATTTTTTGAAAATATTTTAACTAAAAATCTTAAACTTGAAATCAGTCAGATGGTGGTTT	A	.	.	.
21	10578400	idless_two_exons_tie	ATAATGTCCCCACTCTACAGTAAGTTTTATGCTAAATTGACTATCAGAAGAGGGCTAAATATATTGGGCTTGATTTTTTGAAAATATTTTAACTAAAAATCTTAAACTTGAAATCAGTCAGATGGTGGTTTTCACCAAGGAA	A	.	.	.
```

Both engines run `--merged --all_refseq --use_given_ref --plugin NearestExonJB`.

### 1.3 What each engine writes

The rows on `compmerge.2.Heart.chr21` and on `ENST00000427445`, whose exons `ENSE00001659869` and `ENSE00002240871` span the same two exons with ids, in the default format. Ensembl VEP 116.2 on the first deletion:

```
idless_two_exons	21:10578401-10578530	-	compmerge.2.Heart.chr21	compmerge.2.Heart.chr21	Transcript	splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,non_coding_transcript_exon_variant,intron_variant	1266-1297	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SOURCE=RefSeq;NearestExonJB=|64|start|89
idless_two_exons	21:10578401-10578530	-	ENSG00000274391	ENST00000427445	Transcript	splice_acceptor_variant,splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant,intron_variant	1159-1190	807-838	269-280	HNVPTLHQMVVF/HX	caTAATGTCCCCACTCTACATCAGATGGTGGTTTTC/catc	-	IMPACT=HIGH;STRAND=1;SOURCE=Ensembl;NearestExonJB=ENSE00001659869|64|start|82
```

vep-rs on the first deletion:

```
idless_two_exons	21:10578401-10578530	-	compmerge.2.Heart.chr21	compmerge.2.Heart.chr21	Transcript	splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,non_coding_transcript_exon_variant,intron_variant	1266-1297	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SOURCE=RefSeq;NearestExonJB=|64|start|82
idless_two_exons	21:10578401-10578530	-	ENSG00000274391	ENST00000427445	Transcript	splice_acceptor_variant,splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant,intron_variant	1159-1190	807-838	269-280	HNVPTLHQMVVF/HX	caTAATGTCCCCACTCTACATCAGATGGTGGTTTTC/catc	-	IMPACT=HIGH;STRAND=1;SOURCE=Ensembl;NearestExonJB=ENSE00001659869|64|start|82
```

Ensembl VEP 116.2 on the second deletion, the tie:

```
idless_two_exons_tie	21:10578401-10578541	-	compmerge.2.Heart.chr21	compmerge.2.Heart.chr21	Transcript	splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,non_coding_transcript_exon_variant,intron_variant	1266-1308	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SOURCE=RefSeq;NearestExonJB=|64|end|89
idless_two_exons_tie	21:10578401-10578541	-	ENSG00000274391	ENST00000427445	Transcript	splice_acceptor_variant,splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant,intron_variant	1159-1201	807-849	269-283	HNVPTLHQMVVFTKE/X	caTAATGTCCCCACTCTACATCAGATGGTGGTTTTCACCAAGGAA/ca	-	IMPACT=HIGH;STRAND=1;SOURCE=Ensembl;NearestExonJB=ENSE00001659869|64|start|82,ENSE00002240871|64|end|89
```

vep-rs on the tie:

```
idless_two_exons_tie	21:10578401-10578541	-	compmerge.2.Heart.chr21	compmerge.2.Heart.chr21	Transcript	splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,non_coding_transcript_exon_variant,intron_variant	1266-1308	-	-	-	-	-	IMPACT=HIGH;STRAND=1;SOURCE=RefSeq;NearestExonJB=|64|start|82,|64|end|89
idless_two_exons_tie	21:10578401-10578541	-	ENSG00000274391	ENST00000427445	Transcript	splice_acceptor_variant,splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant,intron_variant	1159-1201	807-849	269-283	HNVPTLHQMVVFTKE/X	caTAATGTCCCCACTCTACATCAGATGGTGGTTTTCACCAAGGAA/ca	-	IMPACT=HIGH;STRAND=1;SOURCE=Ensembl;NearestExonJB=ENSE00001659869|64|start|82,ENSE00002240871|64|end|89
```

Each deletion writes 94 rows on models whose exons carry no id, the same 94 in both engines; 90 of them overlap both exons and 4 the second exon alone. The `NearestExonJB` cell on those rows:

| Record | Rows | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- | --- |
| first deletion, both exons | 90 | `\|64\|start\|82` on 48, `\|64\|start\|89` on 42 | `\|64\|start\|82` on 90 |
| first deletion, second exon alone | 4 | `\|75\|end\|89` | `\|75\|end\|89` |
| tie, both exons | 90 | `\|64\|start\|82` on 48, `\|64\|end\|89` on 42 | `\|64\|start\|82,\|64\|end\|89` on 90 |
| tie, second exon alone | 4 | `\|64\|end\|89` | `\|64\|end\|89` |

The rows on models with exon ids agree between the engines on both deletions, the tie's two items in either order. These values are the same in every run of Ensembl VEP: an id-less row has one key, so the default branch reads only that key's entry at the minimum (line 268), and the order of the keys plays no part.

### 1.4 Why Ensembl VEP writes that

Both exons overlap the deletion, so line 219 files both. The first exon is 64 bases from its start and 112 from its end (10578530 against 10578418), the second 116 from its start and 75 from its end. Lines 228-229 write `$dists{''}{64} = 'start'` and `$dists{''}{len} = 82` for the first, lines 232-233 `$dists{''}{75} = 'end'` and `$dists{''}{len} = 89` for the second, each `len` replacing the one before. Line 268 finds the key's entry at the minimum, 64, and line 269 prints it with the key's `len`: `|64|start|89` when the second exon was filed last, `|64|start|82` when the first was.

Which exon is filed last is the interval tree's order (`BaseTranscriptVariation.pm`, line 865), fixed for a transcript but not the same across transcripts. On 42 of the 90 models the second exon comes last and the length is wrong; on 48 the first does and the item is the first exon's. On the tie both exons are 64 away, so the second write of `$dists{''}{64}` replaces the first's side as well: the row prints `|64|end|89` or `|64|start|82`, one exon's item where the model with ids prints two.

### 1.5 Why vep-rs's row is the correct one

The exon at 10578337-10578418 is 82 bases long and its start is 64 bases before the deletion's first base; the exon at 10578517-10578605 is 89 bases long and its end 75 bases past the deletion's last base. `|64|start|89` names no exon: no exon of 89 bases has a junction 64 bases from the deletion. Ensembl VEP agrees wherever the exons carry ids: `ENSE00001659869`, the same 10578337-10578418, gets `ENSE00001659869|64|start|82` on `ENST00000427445` in both engines, and on the tie both engines write that exon and `ENSE00002240871|64|end|89`. vep-rs's rows on the id-less model are those rows with the ids blank.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <merged cache> --assembly GRCh38 \
  --merged --all_refseq --use_given_ref --plugin NearestExonJB
```

vep-rs writes `NearestExonJB=|64|start|82` on the `compmerge.2.Heart.chr21` row of the first deletion and `|64|start|82,|64|end|89` on the tie's. Ensembl VEP 116.2 writes `|64|start|89` and `|64|end|89`. The `ENST00000427445` rows of the same record carry the per-exon items with their ids in both engines.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies `nearest_exon_jb_id_less_exon_fold` when a comparison case names it. A `NearestExonJB` value is set aside, counted once per row, where both engines fill it with items of an empty exon id and Ensembl VEP's one item has the distance and the side of one of vep-rs's items. The row is one vep-rs gives two or more items, or one whose consequence carries both splice-site terms and whose cDNA position names two positions, a variant that begins in one exon and ends in another.

Any other difference stays one: a value with an exon id on either side, an Ensembl VEP item whose distance and side none of vep-rs's items has, and a deletion over one whole exon. That deletion's row carries both terms, but its cDNA position names no two positions, and it overlaps that exon alone. The golden harness compares the `GRCh38-plugin-computed` corpus's `NearestExonJB` cells as multisets, which that corpus's Ensembl cache, every exon with an id, never puts to the test.

## 2. An intronic record of such a transcript under `intronic=1`

**In short.** `21 10521800 idless_intronic T A` on GRCh38 lies in the first intron of `compmerge.2.Heart.chr21`, 106 bases past exon 1 (10521560-10521694, 135 bases) and 2,780 bases before exon 2 (10524580-10524688, 109 bases). Under `--plugin NearestExonJB,intronic=1`, Ensembl VEP 116.2 writes one item, `|106|end|109`, exon 1's distance and junction with exon 2's length; vep-rs writes `|106|end|135,|2780|start|109`. Which exon's distance Ensembl VEP prints follows Perl's hash order: under a fixed hash seed every such row printed exon 1's, and in three runs under the default order 55, 54 and 54 of the 114 rows printed exon 1's and the rest `|2780|start|109`.

### 2.1 What a reader needs to know first

Under `intronic=1` a row with no overlapped exon whose terms include `intron_variant` reports the two exons flanking its intron, each with its own distance, junction and length, two items. The flanking exons are read from the model in transcript order, the downstream exon second. Perl walks a hash's keys in an order that changes between runs and, from one hash to the next, within a run.

### 2.2 The record

```
21	10521800	idless_intronic	T	A	.	.	.
```

Both engines run `--merged --all_refseq --use_given_ref --plugin NearestExonJB,intronic=1`.

### 2.3 What each engine writes

The rows on `compmerge.2.Heart.chr21` and on `ENST00000427445` in the default format. Ensembl VEP 116.2:

```
idless_intronic	21:10521800	A	compmerge.2.Heart.chr21	compmerge.2.Heart.chr21	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq;NearestExonJB=|106|end|109
idless_intronic	21:10521800	A	ENSG00000274391	ENST00000427445	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl;NearestExonJB=ENSE00002206977|2780|start|109,ENSE00003729301|106|end|142
```

vep-rs:

```
idless_intronic	21:10521800	A	compmerge.2.Heart.chr21	compmerge.2.Heart.chr21	Transcript	intron_variant,non_coding_transcript_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq;NearestExonJB=|106|end|135,|2780|start|109
idless_intronic	21:10521800	A	ENSG00000274391	ENST00000427445	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl;NearestExonJB=ENSE00003729301|106|end|142,ENSE00002206977|2780|start|109
```

Ensembl VEP writes 114 rows on models whose exons carry no id, every one `|106|end|109` under a fixed hash seed. In three runs under Perl's default hash order the same rows split 55 and 59, 54 and 60, 54 and 60 between `|106|end|109` and `|2780|start|109`, a different set of rows each run. vep-rs ran on a cache slice that holds 112 of those models and lacks the other two. On each of the 112 it writes `|106|end|<exon 1's length>,|2780|start|109`, exon 1's length by the model (66 to 168 bases), and the same two items whichever run. The 13 rows on models with exon ids carry both exons with their own lengths in both engines, in either order.

### 2.4 Why Ensembl VEP writes that

No exon overlaps the record and the row's terms include `intron_variant`, so lines 190-215 take exon 1 and exon 2 from the sorted exon list. Line 219 files exon 1 first, `$dists{''}{106} = 'end'` and `$dists{''}{len} = 135`, then exon 2, `$dists{''}{2780} = 'start'` and `$dists{''}{len} = 109`, so the key holds both distances and exon 2's length. Two exons under `intronic=1` take the branch at line 244, which prints one item per key: lines 251-259 walk the key's inner hash, `106`, `2780` and `len`, in hash order, and the last distance key read sets the item's distance and side. The item is `|106|end|109` or `|2780|start|109` by that order, which Perl draws afresh for each hash and each run.

### 2.5 Why vep-rs's row is the correct one

The two flanking exons are facts of the model: exon 1 ends 106 bases before the record and is 135 bases long, exon 2 starts 2,780 bases after it and is 109 bases long. `|106|end|109` pairs exon 1's junction with exon 2's length, and `|2780|start|109` is exon 2's item alone, so neither run of Ensembl VEP reports exon 1 as it is. On `ENST00000427445`, whose exons carry ids, both engines write both exons with their own lengths, `ENSE00003729301|106|end|142` and `ENSE00002206977|2780|start|109`. vep-rs's row on the id-less model is that row with the ids blank and exon 1's own length.

### 2.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <merged cache> --assembly GRCh38 \
  --merged --all_refseq --use_given_ref --plugin NearestExonJB,intronic=1
```

vep-rs writes `NearestExonJB=|106|end|135,|2780|start|109` on the `compmerge.2.Heart.chr21` row in every run. Ensembl VEP 116.2 writes `|106|end|109` or `|2780|start|109`, and running it twice shows rows changing between the two.

### 2.7 How the comparator treats this class

As in [1.7](#17-how-the-comparator-treats-this-class): the comparator sets the value aside where a comparison case names the class.
