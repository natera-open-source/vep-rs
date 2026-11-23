# The HGVSc of inserted bases written outside a transcript

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 writes no HGVSc and no `HGVS_OFFSET` for inserted bases that a record writes outside a transcript, directly past its last base or before its first, where the HGVS nomenclature places the same bases on the transcript. Its HGVSIntronOffset plugin writes on such a row the two offsets of a notation read off one base of the transcript, or none. On such a row vep-rs writes the duplication of the transcript's last bases, or the duplication or insertion at the most 3' placement of a run that crosses the first base, and the plugin's two fields beside it.

**Why this is a defect.** Inserted bases that repeat the bases beside them have several placements that give one changed sequence, and a record writes one of them. The HGVS nomenclature describes the change, whichever placement the record writes: its 3' rule assigns the change to the most 3' position, and a copy inserted directly after its original is a duplication. Ensembl VEP 116.3 applies both on a row inside a transcript. `2 48033890 801709 C CT` copies the base at 48033891, and on ENST00000405808, where that base lies in an intron, both engines write `c.169+1443dup`.

The same base is the last base of ENST00000465204, n.3158, and there Ensembl VEP writes no HGVSc. The copy is written directly past that base, between it and the next, so the row reads `downstream_gene_variant`, and Ensembl VEP asks for no notation outside a transcript. The duplicated base is a base of the transcript all the same. The plugin asks on every row and receives a notation read off the last base alone. That is none for this record. For the records that insert two and three bases at the same point it is a duplication of one base, with the offsets 0 and 0, and an insertion at a position past the last base, with 1 and 0.

**What vep-rs does instead.** On a row outside a transcript whose allele adds bases, vep-rs reads the run of placements at which the bases give the same sequence, along the reference bases on either side of the written point. Where the most 3' placement stands directly past the last base and the bases repeat the transcript's last bases, it writes the duplication of those bases, with no `HGVS_OFFSET`, nothing having moved. Where the run starts at or before the first base and ends between two bases of the transcript, it writes the notation at that point.

That notation is a duplication where the copied bases all lie in the transcript, and an insertion between the two bases where the added bases copy nothing. `HGVS_OFFSET` is the count of bases from the written point to the described one, negative on a reverse strand. The row keeps its term, its `DISTANCE` and its empty position columns, and has no HGVSp. The HGVSIntronOffset plugin's two fields are the intron parts of the notation's coordinates, as on every row with an HGVSc.

Every other row outside a transcript carries no HGVSc: a single placement outside, bases that repeat the last bases in part, a run that goes on past the last base, and a copy whose original starts before the first base. Under `--shift_3prime 1` the row is annotated at the shifted span. A record the shift carries inside is then a row inside the transcript, and a row the shift leaves outside carries no HGVSc on either engine. Under `--shift_hgvs 0`, without a reference FASTA and on a RefSeq model that carries an `_rna_edit`, a row outside carries none either.

**Where it is in Ensembl VEP's source** (`ensembl-vep`, `ensembl-variation` and `VEP_plugins` at release/116):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1730 and 1759-1764: HGVSc and `HGVS_OFFSET` are written only on a row within the feature. `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 460-481, sets that flag from the overlap of the transcript with the record's span, as written or, under `--shift_3prime 1`, as shifted. An insertion written directly beside a transcript, or before it, does not overlap it.
- `HGVSIntronOffset.pm`, line 62: the plugin asks the allele for its transcript notation on every Transcript row, and lines 64-72 write the two offsets of whatever notation the call returns.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2726-2733, `_var2transcript_slice_coords`: line 2729 refuses a span wholly before or wholly past the transcript, so the call returns nothing for an insertion written further before the first base. An insertion directly beside the transcript has one end inside, and lines 2732-2733 clamp the other onto the first or last base.
- The same file, `hgvs_transcript`, lines 1446 and 1474-1483: the length of the 3' shift is added to both clamped ends, and the whole inserted allele is read against the one base there. `Bio/EnsEMBL/Variation/Utils/Sequence.pm`, line 526, takes the reference allele from that span. Lines 529-531 return nothing for an inserted allele equal to that base, lines 549-555 name a substitution for one other base, lines 602-611 a `dup` for the base twice and a repeat for more copies, and line 616 a `delins` for any other allele.
- `TranscriptVariationAllele.pm`, lines 2204-2242 and 2271-2284, `_clip_alleles`, which line 1492 calls for every kind but a `dup`: the bases the two alleles share are stripped, and bases left against an empty reference become an `ins`. Past a reverse-strand transcript's last base that `ins` starts at the position after the base, which lines 2802 and 2820-2827 of `_get_cDNA_position` number "last base +1". Lines 1543 and 1545 store the intron parts the plugin reads.

**What correcting it would change.** Ensembl VEP would write, on a row outside a transcript, the HGVSc and `HGVS_OFFSET` it writes for the same change on a row inside one. That is the duplication of the transcript's last bases, or the duplication or insertion at the most 3' placement of a run that crosses the first base. The plugin's two fields would be the intron parts of that notation's coordinates, and a row with no notation would carry neither.

**The class it produces.** `hgvs_added_bases_written_outside_the_transcript`: a differing HGVSc, `HGVS_OFFSET` or plugin field on a row to which both engines give the one term `upstream_gene_variant` or `downstream_gene_variant`. The row's allele is an insertion's, or holds more bases than its Location spans. Ensembl VEP's side prints no HGVSc and no `HGVS_OFFSET`, and carries both plugin fields or neither. vep-rs's side prints a `dup` or an `ins` of the row's transcript, of as many bases as the allele holds or adds to its span, and the plugin's two fields, the intron parts of that notation's coordinates in their printed order.

On a downstream row the notation is a `dup` beside no `HGVS_OFFSET`. On an upstream row `HGVS_OFFSET` is a non-zero integer, negative on a reverse-strand transcript. On an insertion's row it is larger in size than the row's `DISTANCE`. On the row of an allele that holds more bases than its Location spans, whose `DISTANCE` counts from that span, it is no smaller in size than the `DISTANCE`. In a run without the plugin neither side carries its fields, and the HGVSc and `HGVS_OFFSET` are the class.

A row with a term inside the transcript stays a difference, and so does an HGVSc or `HGVS_OFFSET` on Ensembl VEP's side, or one plugin field there without the other. On vep-rs's side a row without a notation, another kind of notation, fields that are not its notation's, an `HGVS_OFFSET` on a downstream row and none on an upstream row stay differences. So does a vep-rs row without both fields beside an Ensembl VEP row that carries them.

The format comparator reads the bases a row adds from the record's alleles where its entry carries them: a VCF or JSON entry, and a default or tab row whose first column spells them. A pure insertion's notation counts the whole allele, and that of an allele that keeps the bases it replaces counts the bases added. A deletion written as a sequence, a substitution and a replacement that keeps neither end of its bases stay differences.

A row under an identifier carries no alleles and is read from its Location: two neighbouring positions read as an insertion, an allele longer than its Location spans as added bases. Two bases replaced by a sequence share an insertion's Location, so a notation that counts the allele's length is accepted there, as is one that counts a longer insertion's surplus over two positions. An allele that adds to the two bases it keeps reads as an insertion: beside a notation in an intron, an `HGVS_OFFSET` equal to the `DISTANCE` stays a difference.

The class is built into the format comparator and the golden harness's classifier and applied once per row. The concordance figures score consequence sets and never see the fields. The committed corpus `GRCh37-hgvs` holds four rows of the class, upstream rows at two Locations on two transcripts each, and its manifest records vep-rs's HGVSc and `HGVS_OFFSET` under it.

A row outside a transcript on which vep-rs writes no notation is not this class. Where Ensembl VEP's plugin writes its two fields on such a row, directly beside a transcript, the row is the one [`hgvs_intron_offset_abutting_insertion`](hgvs-intron-offset-abutting-insertion.md) describes.

**How often the shape occurs.** Tier 2, `--hgvs` being needed. On the genome-wide ClinVar GRCh37 file vep-rs writes the notation on 338 rows outside a transcript, from 246 records and 198 transcripts. Of these, 92 are downstream rows of reverse-strand transcripts, each a duplication of the last bases. The other 246 are upstream rows of forward-strand transcripts, 159 duplications and 87 insertions: 35 of a record written directly before the first base and 211 of one written further before it. Ensembl VEP writes no HGVSc and no `HGVS_OFFSET` on any of the 338. Its plugin writes no field on 274 of them, 0 and 0 on 53 and 1 and 0 on 11.

On the GRCh38 file vep-rs writes the notation on 850 such rows, from 378 records and 468 transcripts. Of these, 230 are downstream rows of reverse-strand transcripts, all duplications. The other 620 are upstream rows of forward-strand transcripts, 351 duplications and 269 insertions: 167 of a record written directly before the first base and 453 of one written further before it. Ensembl VEP writes no HGVSc and no `HGVS_OFFSET` on any of the 850. Its plugin writes no field on 581 of them, 0 and 0 on 245 and 1 and 0 on 24.

Thirteen files of three population datasets were run with the plugin: eight chromosomes of 1000 Genomes Phase 3 on GRCh37, chromosome 6 of 1000 Genomes high-coverage, and a part of each of four gnomAD v4.1 chromosomes on GRCh38. There vep-rs writes the notation on 3,834 rows outside a transcript. 1,851 are downstream rows of reverse-strand transcripts, all duplications. 1,983 are upstream rows of forward-strand transcripts, 1,236 duplications and 747 insertions, 17 of them of an allele that keeps a base it shares with its reference.

Ensembl VEP's plugin writes no field on 2,162 of the 3,834 rows, 0 and 0 on 811 and 1 and 0 on 853. On 8 upstream rows, whose notation lies in an intron, it writes another pair.

A record written at the leftmost of its placements stands at the most 3' one on a reverse-strand transcript and at the most 5' one on a forward-strand transcript. That is why the downstream rows of these files are reverse-strand rows and the upstream rows forward-strand rows.

## 1. A base written past the last base of ENST00000465204 that repeats it

**In short.** `2 48033890 801709 C CT`, a ClinVar record on GRCh37, inserts a `T` between 48033890 and 48033891. ENST00000465204, on the reverse strand, runs from 48033891 to 48038800: 48033891 is its last base, n.3158, an `A` on its strand, and the inserted base is a copy of it written just past it. Both engines write `downstream_gene_variant` and `DISTANCE=0`. Ensembl VEP 116.3 writes no HGVSc and neither plugin field. vep-rs writes `HGVSc=ENST00000465204.1:n.3158dup` with `HGVS_IntronStartOffset=0` and `HGVS_IntronEndOffset=0`.

### 1.1 What a reader needs to know first

An insertion sits between two neighbouring bases. Ensembl VEP gives it a start one above its end, here 48033891 and 48033890, and annotates it outside a transcript when that span does not overlap the transcript's. Such a row carries `upstream_gene_variant` or `downstream_gene_variant`, a `DISTANCE` and no position columns.

Under `--hgvs` a row inside a transcript carries the HGVSc of its change, and `HGVS_OFFSET` where the 3' shift moved the change before it was described. The HGVSIntronOffset plugin adds the intron parts of the notation's start and end coordinates, 0 for an exonic coordinate.

ENST00000465204 (ENSG00000138081) has five exons on the reverse strand of chromosome 2 and 3,158 spliced bases. Its last exon is 48033891 to 48035386, so n.3158 is the base at 48033891. The reference holds a `C` at 48033890 and eighteen `T` from 48033891 to 48033908. On the transcript's strand its last eighteen bases, n.3141 to n.3158, are `A`.

### 1.2 The record

```
2	48033890	801709	C	CT	.	.	.
```

Both engines run `--hgvs --plugin HGVSIntronOffset` on GRCh37 with a reference FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.3, default format:

```
801709	2:48033890-48033891	T	ENSG00000138081	ENST00000465204	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=-1
```

vep-rs:

```
801709	2:48033890-48033891	T	ENSG00000138081	ENST00000465204	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=-1;HGVSc=ENST00000465204.1:n.3158dup;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

Every other column is equal. On ENST00000405808, a transcript of the same gene with 48033891 in an intron, both engines write `HGVSc=ENST00000405808.1:c.169+1443dup` and the two fields, 1443 and 1443.

Two more ClinVar records write copies of the same base at the same point, `2 48033890 89510 C CTT` and `2 48033890 89511 C CTTT`. Ensembl VEP writes no HGVSc on any of the three rows of ENST00000465204:

| Record | Inserted | Ensembl VEP's two fields | vep-rs |
| --- | --- | --- | --- |
| 801709 | `T` | none | `n.3158dup`, 0 and 0 |
| 89510 | `TT` | 0 and 0 | `n.3157_3158dup`, 0 and 0 |
| 89511 | `TTT` | 1 and 0 | `n.3156_3158dup`, 0 and 0 |

### 1.4 Why Ensembl VEP writes that

The row lies outside the transcript, so `OutputFactory.pm` line 1730 skips the HGVSc and `HGVS_OFFSET`. The plugin calls `hgvs_transcript` itself. The transcript spans 4,910 bases of the genome, and line 2708 maps the insertion to the slice coordinates 4911 and 4910. One of them lies inside, so line 2729 lets the pair pass, and line 2732 clamps 4911 to 4910, the last base. The record stands at its most 3' placement on this strand, so the shift adds nothing.

The call then reads the inserted allele against that one base. For 801709 the allele is the base itself, `A` on the transcript's strand, and `Utils/Sequence.pm` lines 529-531 return nothing: the call has no notation, and the plugin writes no field. For 89510 the allele is the base twice, lines 602-606 name a `dup` of the one base, n.3158, and the plugin writes its offsets, 0 and 0.

For 89511 the allele is the base three times, a repeat, which `_clip_alleles` clips. It strips the base the two alleles share and leaves an `ins` of two bases that starts at slice position 4911. That position lies before the first exon in genomic order, and `_get_cDNA_position` numbers it `3158+1`. The plugin writes 1 and 0: an intron offset at a position where the transcript has no intron.

### 1.5 Why vep-rs's row follows

The record adds a `T` after 48033890, and the reference holds a `T` at 48033891: the added base is a copy of that base, placed directly beside it. On the transcript's strand the copy is an `A` directly 3' of n.3158, itself an `A`. No placement lies further 3': the base beyond the copy is the `C` at 48033890.

The HGVS Nomenclature (`hgvs-nomenclature.org/stable/`) names the kind. Its Duplication page, `recommendations/DNA/duplication/`, defines a duplication as "a sequence change where, compared to a reference sequence, a copy of one or more nucleotides is inserted directly 3' of the original copy of that sequence". Its notes add: "when a variant can be described as a duplication, it must be described as a duplication and not as, e.g., an insertion".

The original copy is n.3158, a base of the transcript, so the description names no position outside it: `n.3158dup`. The duplicated base is exonic, which gives the plugin's fields 0 and 0. The record does not move, so the row has no `HGVS_OFFSET`. It is the duplication both engines write as `c.169+1443dup` on ENST00000405808, of the same base. The row keeps `downstream_gene_variant` and `DISTANCE=0`: the term says where the record writes the base, and the notation says what changes.

### 1.6 Try it yourself

Save the record of 1.2 as `record.vcf` and annotate it against a GRCh37 cache with the reference FASTA. vep-rs's command is below; Ensembl VEP takes the same record with its own cache and the FASTA under `--offline --cache`:

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh37 --hgvs --plugin HGVSIntronOffset
```

vep-rs's row for ENST00000465204 carries `HGVSc=ENST00000465204.1:n.3158dup` and the two fields. Ensembl VEP 116.3's ends `DISTANCE=0;STRAND=-1`. The reproduction `scripts/divergences/reproductions/hgvs_added_bases_written_outside_the_transcript.toml` runs the six records of this page on both engines.

### 1.7 How the comparator treats this class

The format comparator's rule (`tools/vep-formats-compare`, `added_bases_written_outside_the_transcript`) accepts the HGVSc and the two plugin fields of this row, `hgvsc`, `hgvs_intronstartoffset` and `hgvs_intronendoffset` in JSON. Both engines give the row the one term `downstream_gene_variant`. Ensembl VEP's side prints no HGVSc and no `HGVS_OFFSET`. vep-rs's prints a `dup` of the row's transcript that names as many bases as the allele holds, no `HGVS_OFFSET`, and the two fields of that notation's coordinates.

On the row of 89510 the two fields are equal on both sides and the HGVSc alone is the class. On the row of 89511 the HGVSc and the start offset are.

## 2. A base written before the first base of ENST00000478223 that a run carries inside

**In short.** `1 2179434 3895388 C CG`, a ClinVar record on GRCh37, inserts a `G` between 2179434 and 2179435, the first base of ENST00000478223 on the forward strand. The transcript opens with seven `G`, so the inserted base gives the same sequence after each of them, and its most 3' placement is after n.7. Both engines write `upstream_gene_variant` and `DISTANCE=0`. Ensembl VEP 116.3 writes no HGVSc, no `HGVS_OFFSET` and the plugin's fields 0 and 0. vep-rs writes `HGVSc=ENST00000478223.2:n.7dup`, `HGVS_OFFSET=7` and the same two fields.

### 2.1 The record

```
1	2179434	3895388	C	CG	.	.	.
```

ENST00000478223 (ENSG00000157933) has three exons on the forward strand, the first from 2179435 to 2179510, so n.1 is the base at 2179435. The reference holds a `C` at 2179434, seven `G` from 2179435 to 2179441 and a `C` at 2179442.

### 2.2 What each engine writes

Ensembl VEP 116.3, default format:

```
3895388	1:2179434-2179435	G	ENSG00000157933	ENST00000478223	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

vep-rs:

```
3895388	1:2179434-2179435	G	ENSG00000157933	ENST00000478223	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1;HGVSc=ENST00000478223.2:n.7dup;HGVS_OFFSET=7;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

On ENST00000378536, a transcript of the same gene with these bases in an intron, both engines write `HGVSc=ENST00000378536.4:c.969+18267dup` and `HGVS_OFFSET=7`.

### 2.3 Why Ensembl VEP writes that

The row lies outside the transcript, and line 1721 skips the HGVSc and `HGVS_OFFSET`. In the plugin's call, line 2712 maps the insertion to the slice coordinates 1 and 0, and line 2733 clamps the 0 to 1. The record's 3' shift is seven bases long, and lines 1477-1478 add 7 to both ends. The call describes the base at slice position 8, the `C` at 2179442, against the inserted `G`: a substitution of n.8, a base the record leaves as it is. Its coordinate is exonic, and the plugin writes 0 and 0.

### 2.4 Why vep-rs's row follows

The inserted `G` gives one changed sequence at eight placements, from before 2179435 to after 2179441: eight `G` where the reference holds seven. The HGVS Nomenclature's general recommendations, `recommendations/general/`, state the 3' rule: "for all descriptions, the most 3' position possible of the reference sequence is arbitrarily assigned to have been changed". They add that the rule "applies to ALL descriptions (genome, gene, transcript, and protein) of a given variant".

The most 3' placement is after n.7. There the added base is a copy of n.7 inserted directly 3' of it, a duplication by the sentences of 1.5: `n.7dup`. `HGVS_OFFSET=7` counts the bases from the written point to that placement. The duplicated base is exonic, so the plugin's fields are 0 and 0. Ensembl VEP carries the same record the same seven bases on ENST00000378536 and writes the duplication of the same base, 2179441, there.

Under `--shift_3prime 1` both engines annotate the record at the shifted span, inside the first exon, and vep-rs writes the same `n.7dup` on that row. What Ensembl VEP writes there is the subject of [`hgvs_shifted_insertion_clamped_end`](hgvs-shifted-insertion-clamped-end.md).

## 3. Two bases written before the first base of ENST00000408108 that copy nothing

**In short.** `8 41517958 802400 A AGG`, a ClinVar record on GRCh37, inserts `GG` between 41517958 and 41517959, the first base of ENST00000408108, a forward-strand miRNA transcript that opens `GT`. The two bases give the same sequence one base further on, between n.1 and n.2, and there they copy nothing. Both engines write `upstream_gene_variant` and `DISTANCE=0`. Ensembl VEP 116.3 writes no HGVSc, no `HGVS_OFFSET` and the plugin's fields 0 and 0. vep-rs writes `HGVSc=ENST00000408108.1:n.1_2insGG`, `HGVS_OFFSET=1` and the same two fields.

### 3.1 The record

```
8	41517958	802400	A	AGG	.	.	.
```

ENST00000408108 (ENSG00000221035) is one exon of 68 bases, 41517959 to 41518026. The reference holds an `A` at 41517958, a `G` at 41517959 and a `T` at 41517960.

### 3.2 What each engine writes

Ensembl VEP 116.3, default format:

```
802400	8:41517958-41517959	GG	ENSG00000221035	ENST00000408108	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

vep-rs:

```
802400	8:41517958-41517959	GG	ENSG00000221035	ENST00000408108	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=0;STRAND=1;HGVSc=ENST00000408108.1:n.1_2insGG;HGVS_OFFSET=1;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

### 3.3 Why each engine writes that

In the plugin's call the insertion is clamped onto the first base and the 3' shift, one base long, is added: the call describes the `T` at n.2 against the inserted `GG`, a `delins` of a base the record leaves as it is. Its coordinate is exonic, and the plugin writes 0 and 0.

Written before n.1 the record gives `A`, `GG`, `G`, `T`. Written between n.1 and n.2 it gives `A`, `G`, `GG`, `T`, the same bases, and no placement lies further 3'. The 3' rule assigns the change to that point. The two bases 5' of it are the `A` at 41517958 and the `G` at n.1, so the inserted `GG` is no copy of them.

The nomenclature's Insertion page, `recommendations/DNA/insertion/`, defines an insertion as a change in which bases are inserted "and where the insertion is not a copy of a sequence immediately 5'". It asks for "two flanking nucleotides, e.g., `123_124`". Both flanks, n.1 and n.2, are bases of the transcript: `n.1_2insGG`, with `HGVS_OFFSET=1` for the one base the rule moved the change.

Where the first exon is shorter than the run, the placement lies in the first intron, and the plugin's fields show it. A 1000 Genomes high-coverage record on GRCh38 inserts 33 bases between 6:170160615 and 6:170160616, directly before the first base of ENST00000801843, whose first exon has 24 bases. A run carries them 32 bases. vep-rs writes an `ins` of the 33 bases between `n.24+8` and `n.24+9`, `HGVS_OFFSET=32` and the fields 8 and 9. Ensembl VEP 116.3 writes no HGVSc and the fields 9 and 9, the intron part of the one base its call describes.

## 4. A base written one position before ENST00000494977 that a run carries inside

**In short.** `12 32735311 637079 A AG`, a ClinVar record on GRCh37, inserts a `G` between 32735311 and 32735312. ENST00000494977, on the forward strand, starts at 32735313, so the record's span lies one base before the transcript: both engines write `upstream_gene_variant` and `DISTANCE=1`. The reference holds four `G` from 32735312 to 32735315, the last three of them the transcript's first three bases, and the most 3' placement of the inserted base is after c.3. Ensembl VEP 116.3 writes no HGVSc, no `HGVS_OFFSET` and neither plugin field. vep-rs writes `HGVSc=ENST00000494977.1:c.3dup`, `HGVS_OFFSET=4` and the fields 0 and 0.

### 4.1 The record

```
12	32735311	637079	A	AG	.	.	.
```

ENST00000494977 (ENSG00000139132) carries `FLAGS=cds_start_NF`: its coding sequence opens at its first base, so the base at 32735313 is c.1. The reference holds an `A` at 32735311, four `G` from 32735312 to 32735315 and a `C` at 32735316.

### 4.2 What each engine writes

Ensembl VEP 116.3, default format:

```
637079	12:32735311-32735312	G	ENSG00000139132	ENST00000494977	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=1;STRAND=1;FLAGS=cds_start_NF
```

vep-rs:

```
637079	12:32735311-32735312	G	ENSG00000139132	ENST00000494977	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=1;STRAND=1;FLAGS=cds_start_NF;HGVSc=ENST00000494977.1:c.3dup;HGVS_OFFSET=4;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0
```

On ENST00000427716, a transcript of the same gene that holds the four `G` in an exon, both engines write `HGVSc=ENST00000427716.2:c.514dup` and `HGVS_OFFSET=4`.

### 4.3 Why each engine writes that

Line 2712 maps the insertion to the slice coordinates 0 and -1. Both lie before the first base, line 2729 returns nothing, and the plugin's call has no notation: Ensembl VEP writes no field.

The inserted `G` gives one changed sequence at five placements, from before 32735312 to after 32735315. The most 3' is after c.3, where the added base is a copy of c.3 inserted directly 3' of it: `c.3dup`, by the 3' rule and the Duplication page's sentences. `HGVS_OFFSET=4` counts the four bases from the written point to that placement, one of them before the transcript, which is why the offset is larger than the row's `DISTANCE`. It is the duplication of the same base, 32735315, that both engines write as `c.514dup` on ENST00000427716.

### 4.4 How the comparator treats an upstream row

The rule accepts the HGVSc, the `HGVS_OFFSET` and the plugin's fields of an upstream row under the class's conditions. Both engines give the row the one term `upstream_gene_variant`. Ensembl VEP's side prints no HGVSc and no `HGVS_OFFSET`. vep-rs's prints a `dup` or an `ins` of the row's transcript, a non-zero `HGVS_OFFSET` larger in size than the row's `DISTANCE`, this row being an insertion's, and the two fields of that notation's coordinates. On an insertion's row an `HGVS_OFFSET` no larger in size than the `DISTANCE` stays a difference, and on any row so does one of the wrong sign for the row's `STRAND`.
