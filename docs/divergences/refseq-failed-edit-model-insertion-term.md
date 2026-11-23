# An insertion at the twice-read position of a failed-edit RefSeq model

**The defect in one sentence.** On a RefSeq model whose alignment correction failed (`BAM_EDIT=FAILED`) and whose coding sequence reads one position twice, Ensembl VEP reads an insertion written at that position as a replacement of one cached base, so one inserted base is `synonymous_variant`. vep-rs reads the insertion between the two genomic bases around that position, so one inserted base is `frameshift_variant`.

**Why this is a defect.** The RefSeq record of such a model equals the genome's exons joined, base for base. Its coding sequence is a join of two spans that share one position, under the qualifier `ribosomal_slippage`: translation reads that position twice. The cache writes the double reading as one extra position of the model, by an edit. An edit is an `_rna_edit` attribute: a first and a last cDNA position and the bases the model holds there.

The cached coding sequence is one base longer than the record's coding span, and the cached mapper, the map between cDNA and genomic positions, gives the extra position no genomic base. This page calls it the cache's own position. It stands for the shared position read once more and is no base of the transcript.

The two genomic bases around the cache's own position are therefore neighbours in the transcript, and an insertion between them adds bases to the coding sequence. With the join's rule kept, one added base shifts the frame of everything after it. Ensembl VEP's row states a replacement instead: `synonymous_variant`, a `CDS_position` of one base and one `Amino_acids` value, for a base the transcript does not hold.

For an inserted base equal to the base at the position read twice, Ensembl VEP's row contradicts itself: its term is `synonymous_variant` and its own HGVSp is a frameshift. The same base written one position on is `frameshift_variant` in both engines. For another base Ensembl VEP's HGVSp is an unchanged residue, as its term says.

**What vep-rs does instead.** vep-rs reads the status the cache holds for the model and takes the genome as the reference of a model whose status is `failed`. An insertion between two exonic bases that the cached mapper numbers more than one position apart is read at the point after the 5' one of the two. The inserted bases are added there and do not replace the cached position between them. For such an insertion the term, `IMPACT`, the position columns, `Codons` and `Amino_acids` are all read at that point. The HGVS notations start from the same point and move 3' along the genome where its bases repeat the inserted ones.

One inserted base is `frameshift_variant`, with the two positions, the lengthened codon and the two `Amino_acids` halves of an insertion. Two inserted bases are `frameshift_variant` as well, and three are `inframe_insertion`. In the HGVSc, an inserted base equal to the bases that follow the point is a duplication of the last of them. Any other inserted bases are an `ins` between the base 5' of the point and the next position of the coding sequence. Sections 1.3 and 1.6 give the rows.

The rule reads the status alone. A model the alignment corrected (`BAM_EDIT=OK`) is left as it is: its cached sequence, with its edits, is the RefSeq sequence, so the bases it holds and the genome lacks are bases of the transcript, and its rows are read on them. A record that keeps both genomic bases around the cache's own position and adds bases between them gets the term, the codons, the residues, the protein position and the HGVS notations of the same bases written as an insertion. Its `cDNA_position` and `CDS_position` stay the record's span; section 2 describes it.

**Where it is in Ensembl VEP's source** (`ensembl` at commit `0d852313c420a1d5128b6774635125ebbac03350`, `ensembl-variation` at commit `db7394a02523a6e519cb0d89f84b2889bac12737`, `ensembl-vep` at commit `2cb0bbe216bb31c75de8f8000e2da7ff4fb7b451`):

- `ensembl`, `Bio/EnsEMBL/TranscriptMapper.pm`, lines 174-210, `_load_mapper`: for an edit that lengthens the cDNA, the pair before the edit ends at the cDNA position before it. The next pair starts one cDNA position further on, at the next genomic base (lines 179, 199, 202 and 209-210). One cDNA position is left with no genomic base.
- `ensembl-variation`, `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, line 1079, `_mapper`: a run reads the mapper the cache holds, which was built with the model's own edit applied.
- `ensembl`, `Bio/EnsEMBL/Mapper.pm`, lines 268-269, and `map_insert`, lines 489-542: an insertion is mapped as the two bases around it. Where the two lie in two pairs, the first result's start moves on and the second's end moves back (lines 524-528 and 538-542), and both land on the position between the pairs.
- `BaseTranscriptVariation.pm`, lines 263-264, `cds_start`: the row's coding span runs from the first mapped result's start to the last one's end. Here the two are equal, a span of one position, where an insertion between two neighbours of one pair has an empty span.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 851-877, `codon`, with lines 2409-2410 and 2421, `_get_alternate_cds`: the alternate coding sequence is the cached one with the allele spliced in place of the span, here one base in place of one.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1556-1564, `frameshift`: the test compares the allele's length with the span's, one against one, and refuses. Line 1082, `synonymous_variant`, accepts where the two peptides are equal.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 1696, 1742-1746, 1760 and 1769: the term and the row's cells are asked first, with no 3' shift made, and each HGVS notation then makes its own shift. For a base equal to the bases after the point, the notations are built three bases on, where the span is an ordinary empty one.
- `ensembl-vep`, `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, lines 579, 621 and 624-625, `apply_edits`: a model is `failed` when its sequence, with the edit it came with and the alignment's edits applied, is not the RefSeq sequence the alignment holds. The alignment's edits are then removed, the edit the model came with stays, and no cDNA is cached.
- `TranscriptVariationAllele.pm`, lines 1517, 1524-1525 and 1605-1617, `get_misalignment_offset`: the position inside the HGVSc is counted along the exons and moved on by an offset for the model's edit. For an edit that writes one position as two bases, line 1617 adds 2, where the transcript gains no base.

**What correcting it would change.** Every row of an insertion written at such a position would carry the term of its inserted length: `frameshift_variant` for one base and `inframe_insertion` for three, where Ensembl VEP writes `synonymous_variant` and `frameshift_variant`. Its `IMPACT`, `cDNA_position`, `CDS_position`, `Codons` and `Amino_acids` would be an insertion's. On a row of a base that differs from the bases after the point, the HGVSp would be a frameshift.

The position inside the HGVSc is a separate difference with a cause of its own, the offset of the last source line above. Section 1.4 describes it, and [The HGVSc position past the own edit of a failed-edit RefSeq model](refseq-failed-edit-model-hgvsc-offset.md) is its entry. No row of a model the alignment corrected would change.

**The class it produces.** None. `refseq_failed_edit_model_insertion_term` names this entry in the index and in its reproduction, and no comparison tool holds a rule, a class or a test under that name. Nothing is set aside for these rows, so a difference on them is charged to vep-rs in the concordance figures. The reproduction also records the entry's input class, which is another attribute than a comparison class: `unobserved`, the index's word for an input that the format and the documentation allow and that no measured data carries.

Nothing the comparisons measure reaches the entry. The concordance datasets run on caches that hold no `_rna_edit` attribute, and the caches of the committed corpora hold no model with `BAM_EDIT=FAILED`. The RefSeq GRCh38 cache of release 116 holds three such models, `NM_015068.3`, `NM_001184961.1` and `NM_001172437.2` of gene `PEG10`, which share one such position. That cache holds 1,218 pairs of a model and a place where the model holds positions the genome lacks, at 613 places, and these three are the only pairs whose model is `FAILED`. None of the 4,398,127 records of genome-wide ClinVar on GRCh38 is an insertion written at their position.

**What is not shown.** Each of these bounds what the page shows:

- The three models lie on the forward strand. Nothing is shown for a reverse-strand model.
- Every row of this page was written under `--refseq --hgvs` on the RefSeq GRCh38 cache of release 116: `--refseq` reads that cache and `--hgvs` writes the HGVS notations. In section 1.6, 6 rows were also written with `--shift_hgvs 0` added, under which the notations make no 3' shift. The RefSeq cache of GRCh37 was not read.
- Apart from the three records of section 1.3, the rows come from constructed records that the reproduction does not hold.

## 1. The twice-read position of `NM_015068.3`

**In short.** `NM_015068.3` (gene `PEG10`, forward strand, `BAM_EDIT=FAILED`) reads `7:94664513` on GRCh38 twice in its coding sequence. The cache writes the double reading as cDNA 1436, a position with no genomic base. For one base inserted between `7:94664512` and `7:94664513`, Ensembl VEP 116.3 writes `synonymous_variant` at `CDS_position` 957, and vep-rs writes `frameshift_variant` at `CDS_position` 956-957. The same base inserted one position on is `frameshift_variant` in both engines.

### 1.1 What a reader needs to know first

The models of a RefSeq cache were corrected against an alignment of the RefSeq sequences. A model whose corrected sequence equals the RefSeq sequence prints `BAM_EDIT=OK` on its rows. A model whose sequence, with the edit it came with and the alignment's edits applied, does not equal it has the status `failed` and prints `BAM_EDIT=FAILED`. The alignment's edits are then removed, the edit the model came with stays, and no cDNA is cached. Ensembl VEP's documentation says of that case: "If the corrected model does not match the original RefSeq sequence in the BAM, the corrected model is discarded."

This page's model came with one edit, the `_rna_edit` attribute `1436 1436` with two bases: it writes the one cDNA position 1436 as two bases.

On a model whose edit failed, the two engines read as follows. Ensembl VEP has no cached cDNA, and the cDNA it reads for the model at run time is the genome's exons joined. Its mapper and its coding sequence are the cached ones, and both hold the model's own edit. An insertion at the edit's position is therefore read on the cached position. vep-rs takes the genome as the model's reference there and reads the insertion between the two genomic bases.

A cDNA position counts from the model's first base, and a CDS position from the first base of its coding sequence. On `NM_015068.3` the coding sequence starts at cDNA 480, so cDNA 1436 is CDS 957, in the codon of residue 319.

### 1.2 The position in RefSeq's record, on the genome and in the cache

RefSeq's record of `NM_015068.3` gives the coding sequence as `join(480..1436,1436..2605)` with the qualifier `ribosomal_slippage`. The two spans share record position 1436, so the join reads it twice. The record's coding span, positions 480 to 2605, holds 2,126 positions, and the join read span by span is 2,127 bases long. The record equals the genome's exons joined over all 6,618 exonic bases: its positions 1435 and 1436 are `7:94664512` and `7:94664513`, neighbours in the model's second exon.

The place is the point between `7:94664512` and `7:94664513`. The run is the three equal bases the genome holds from `7:94664513` to `7:94664515`; `7:94664516` holds another base.

| Where | The two bases around the place | The run of equal bases |
| --- | --- | --- |
| The genome | `7:94664512` and `7:94664513`, neighbours | 3 bases, `7:94664513` to `7:94664515` |
| RefSeq's record | positions 1435 and 1436, neighbours | 3 bases, from position 1436 |
| The cache | cDNA 1435 and cDNA 1437, with cDNA 1436 between | 4 bases, CDS 957 to 960 |

Two pairs of the cached mapper meet at the place: cDNA 257-1435 on `7:94663334-94664512`, and cDNA 1437-6619 on `7:94664513-94669695`. No pair holds cDNA 1436, the cache's own position.

The cached coding sequence is the join read span by span, 2,127 bases, base for base. It holds the run's base four times, at CDS 957 to 960. The position read twice is `7:94664513`, and the cache spells it as cDNA 1436 and cDNA 1437. In that numbering the position without a genomic base, cDNA 1436, stands first of the two, and cDNA 1437 maps to `7:94664513`.

The models `NM_001184961.1` and `NM_001172437.2` share the place. Their own positions are CDS 1059 and CDS 1185.

### 1.3 The three records and both engines' rows

Three one-base insertions, written for this example as VCF lines, are run under `--refseq --hgvs`. This page prints none of their bases; the reproduction of section 1.7 holds the lines and the rows with their letters.

The first record inserts one copy of the run's base at the place, between `7:94664512` and `7:94664513`. Its row on `NM_015068.3`:

| Cell | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `synonymous_variant` | `frameshift_variant` |
| `cDNA_position` | 1436 | 1435-1436 |
| `CDS_position` | 957 | 956-957 |
| `Protein_position` | 319 | 319 |
| `Amino_acids` | one value | two halves |
| `Codons` | two equal halves of 3 bases | halves of 3 and 4 bases |
| `HGVSc` | a duplication at `c.961` | a duplication at `c.959` |
| `HGVSp` | a frameshift first changing residue 321 | the same notation |
| `HGVS_OFFSET` | 3 | 3 |

The second record inserts the same base one position on, between `7:94664513` and `7:94664514`. Both engines write `frameshift_variant`, `cDNA_position` 1437-1438, `CDS_position` 958-959 and the frameshift HGVSp of the first record. Each writes the duplication it writes on the first record, `c.961` in Ensembl VEP and `c.959` in vep-rs, with `HGVS_OFFSET=2`.

The third record inserts, at the place, one base that is not the run's. Its row on `NM_015068.3`:

| Cell | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `synonymous_variant` | `frameshift_variant` |
| `cDNA_position` | 1436 | 1435-1436 |
| `CDS_position` | 957 | 956-957 |
| `Protein_position` | 319 | 319 |
| `Amino_acids` | one value | two halves |
| `Codons` | halves of 3 bases that differ in one base | halves of 3 and 4 bases |
| `HGVSc` | an `ins` between `c.958` and `c.959` | an `ins` between `c.956` and `c.957` |
| `HGVSp` | residue 319 unchanged | a frameshift first changing residue 321 |
| `HGVS_OFFSET` | none | none |

vep-rs prints `IMPACT=HIGH` on its `frameshift_variant` rows, and both engines print `REFSEQ_OFFSET=2` on all three rows.

Six constructed records insert one base at the place, and on the three models they give 18 rows: 6 rows of the run's base and 12 rows of another base. This page does not describe how the six records are written. Ensembl VEP writes `synonymous_variant` on all 18 and vep-rs `frameshift_variant` on all 18.

### 1.4 The position inside the coding notation

The rows of 1.3 differ in the position inside the HGVSc as well, by a second cause. Ensembl VEP counts that position along the exons and moves it on by an offset for the model's edit, the value it prints as `REFSEQ_OFFSET`. The offset is the number of bases the edit writes, 2. vep-rs counts the position along the exons alone, which are the positions of RefSeq's record.

The duplicated base is the run's last, `7:94664515`. Counted along the exons it is cDNA 1438, `c.959` from the coding start, the position vep-rs writes and the record's own position of that base. With the offset it is cDNA 1440, `c.961`. No count of the transcript's bases gives the run's last base that number, so Ensembl VEP's notation names a base that is not the one duplicated.

For the third record, the two genomic bases around the place are cDNA 1435 and 1436 along the exons, `c.956` and `c.957`, and vep-rs writes the `ins` between them. With the offset they are cDNA 1437 and 1438, `c.958` and `c.959`, so Ensembl VEP's `ins` stands two positions past the insertion point.

On the three models Ensembl VEP's duplication is at `c.961`, `c.1063` and `c.1189`, and vep-rs's at `c.959`, `c.1061` and `c.1187`, on 12 of 12 duplication rows. The 12 are the 6 rows of the run's base at the place and the 6 rows of two records that write the run's base off the place, 3 rows each. One is written one position on, as the second record of section 1.3, and the other further on, at a point this page does not give. On those 6 both engines write `frameshift_variant` and a frameshift HGVSp. The offset accounts for the position on all 12, and for the two positions of Ensembl VEP's `ins` on the 12 rows of another base.

### 1.5 The other reading, and why vep-rs does not take it

One reading makes `synonymous_variant` the right term: the transcript with the insertion read straight through, with no position read twice. RefSeq's join reads the shared position twice, so the coding sequence it writes holds the run's base four times. A transcript with one more copy of the run's base holds four copies. Read straight through, those four give the reference coding sequence itself, and a record of another base changes one base of one codon.

Under that reading the run's base is no change wherever in the run it is written. The second record of 1.3 would then be no change as well, where both engines write `frameshift_variant`.

vep-rs does not take that reading. RefSeq's record states the slip for the transcript, and no source states a reading of the transcript without it. With the join's rule kept, the shared position read twice, the coding sequence with the insertion is 2,128 bases against 2,127, and the frame shifts after the added base. No record says how a transcript one base longer is translated. Keeping the join's rule for it is a rule of vep-rs and not a statement of any source.

### 1.6 Other lengths, records written before the place, and a kept base

Constructed records insert two and three bases at the place, 3 rows each on the three models:

| Inserted at the place | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| two bases | `frameshift_variant`, the cache's own position, `Codons` halves of 3 and 4 bases | `frameshift_variant`, the point's two positions, halves of 3 and 5 bases |
| three bases | `frameshift_variant` | `inframe_insertion`, `IMPACT=MODERATE`, the point's two positions, halves of 3 and 6 bases |

For three bases Ensembl VEP's `frameshift` test compares three bases with a span of one position, and the remainder of two is a frameshift (line 1564). vep-rs writes an `ins` between the base 5' of the point and the next position, and the HGVSp of an in-frame change at residue 319.

Two constructed records are written before the place: one base one position before it, and two bases two positions before it. On the genome, the 3' shift that each HGVS notation makes moves both onto the place, with `HGVS_OFFSET` 1 and 2. On their 6 rows the term and the position cells are those of the written point in both engines. Both engines also write the same HGVSc, the duplication of the base or the two bases before the place, `c.956` and `c.955_956` on `NM_015068.3`.

The HGVSp of those 6 rows differs. vep-rs writes a frameshift on all 6. On `NM_015068.3` it first changes residue 321 for the first record and residue 320 for the second. Ensembl VEP's protein notation makes its own shift onto the place and reads the allele in place of the cache's own position there (`TranscriptVariationAllele.pm`, lines 1711 and 1732-1733). Its notation is an unchanged residue on the first record's 3 rows, beside the row's own `frameshift_variant`, and a frameshift with another distance to the new stop on the second's. With `--shift_hgvs 0` added to `--refseq --hgvs`, Ensembl VEP's notations make no 3' shift, and it writes vep-rs's notation on all 6.

Two further constructed records, each of several alternate alleles, hold a one-base insertion at the place: the run's base in one and another base in the other. In each, the allele keeps a base it shares with the reference. Beside another alternate allele the parse leaves such a base in the allele, where it trims it from a record of one alternate allele.

On their 6 rows on the three models, vep-rs writes the term and the HGVSp of the same base written as a plain insertion. Their HGVSc follows vep-rs's rule for an allele that keeps a shared base, which other entries document. This page states nothing of that cell, of which base is kept or of Ensembl VEP's rows on those records.

### 1.7 How the comparator treats these rows

Nothing is set aside under this entry's name: no comparison tool holds a rule, a class or a test under it. A comparison of the two engines' outputs counts every differing cell of these rows against vep-rs, but for one. The position inside the `HGVSc` of the second record's rows, whose `cDNA_position` the engines write alike, is the offset's, and the format comparator sets it aside where the comparison names that entry's class, `refseq_failed_edit_model_hgvsc_offset`.

The reproduction `scripts/divergences/reproductions/refseq_failed_edit_model_insertion_term.toml` runs the three records of 1.3 on both engines and records nine cells of each row on `NM_015068.3`. No committed corpus's cache holds a model with `BAM_EDIT=FAILED`, so the reproduction converts the one region file that holds the three models from an Ensembl VEP RefSeq GRCh38 cache it is given. It holds no record of sections 1.6 and 2.

## 2. A record that keeps both bases around the cache's own position

**In short.** A record whose reference is the two genomic bases around the place, and whose alternate allele keeps both with other bases between them, states an insertion at the place. Ensembl VEP reads it across the cache's own position, on a span of three positions, so one added base is `synonymous_variant`. vep-rs gives the row the term, the codons, the residues, the protein position and the HGVS notations of the same bases written as an insertion. `cDNA_position` and `CDS_position` stay the record's span in both engines.

### 2.1 What Ensembl VEP reads

An allele keeps both bases through the parse only beside a second alternate allele, so the shape is an allele of a record of several alternates. Its two reference bases, `7:94664512` and `7:94664513`, lie in the two pairs of the cached mapper, at cDNA 1435 and cDNA 1437 on `NM_015068.3`.

The row's span runs from the first mapped base's start to the last one's end (`BaseTranscriptVariation.pm`, lines 263-264). That is three positions, cDNA 1435-1437 and CDS 956-958, with the cache's own position in the middle. `_get_alternate_cds` then splices the whole alternate allele in place of the span (`TranscriptVariationAllele.pm`, lines 2409-2410 and 2421). The two kept bases stand in their own places, and one added base stands where the cache's own position stood, so the row reads as no change of length. Two added bases read as one base more, and three as two more.

### 2.2 One row, and the eighteen

A constructed record has the two reference bases and, as its first alternate, those two with one other base between them. Its row on `NM_015068.3`:

| Cell | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `synonymous_variant` | `frameshift_variant` |
| `cDNA_position` | 1435-1437 | 1435-1437 |
| `CDS_position` | 956-958 | 956-958 |
| `Protein_position` | 319-320 | 319 |
| `Codons` | six bases against six | three bases against four |
| `HGVSc` | an `ins` of three bases between `c.956` and `c.957` | an `ins` of one base between `c.956` and `c.957` |
| `HGVSp` | residues unchanged | a frameshift |

Ensembl VEP's HGVSc is an `ins` of the whole alternate allele, the two kept bases with the added one. The offset of section 1.4 does not reach this row. Ensembl VEP counts an edit into the offset only where the row's cDNA start is at or after the edit's first position (`TranscriptVariationAllele.pm`, lines 1608-1610), and this row starts at cDNA 1435, before the edit at cDNA 1436. Its two positions are the exons' count of the two bases around the place, `c.956` and `c.957`.

Apart from the two position cells, vep-rs's cells in this table are those it writes for the same base as a plain insertion at the place.

Constructed records on the three models give 18 rows of the shape, from six sets of added bases: four of one base, one of two and one of three. Ensembl VEP writes `synonymous_variant` on the 12 rows of one added base. On the other 6 it writes `stop_gained,frameshift_variant` on 3 and `frameshift_variant` on 3. vep-rs writes `frameshift_variant` on 15 rows and `inframe_insertion` on 3. Each of its rows equals the row of the same bases written as an insertion in `Consequence`, `IMPACT`, `Codons`, `Amino_acids`, `Protein_position`, `HGVSc`, `HGVSp` and `HGVS_OFFSET`.

### 2.3 The two position cells

`cDNA_position` and `CDS_position` of these 18 rows are not settled. Both engines write the record's span on the cache's numbering, three positions with the cache's own position inside, and the two cells are equal between the engines on all 18 rows. The same bases written as an insertion carry two positions in vep-rs's row and one position in Ensembl VEP's, as section 1.3 shows. No page of the HGVS Nomenclature decides between the span of three and the span of two for such a record, and the cells stay as both engines write them.

### 2.4 What reaches the shape

No row of a concordance dataset, a committed corpus or a comparison case reaches the shape. The 18 rows are constructed, and the reproduction holds no record of it.
