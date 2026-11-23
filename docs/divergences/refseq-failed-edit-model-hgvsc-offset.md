# The HGVSc position past the own edit of a failed-edit RefSeq model

**The defect in one sentence.** Under `--hgvs`, on a RefSeq model whose alignment correction failed (`BAM_EDIT=FAILED`) and whose own edit changes the length of the cDNA, Ensembl VEP 116.3 numbers an HGVSc position past the edit with the edit counted. The model's own edit is the `_rna_edit` attribute it came with, which the alignment did not write. vep-rs counts every position along the genome's exons, which are the model's transcript.

**Why this is a defect.** A failed model keeps no edit of the alignment, so its transcript is the genome's exons joined. Each of the ten such models is translated with a ribosomal frameshift. RefSeq's record gives its coding sequence as a join that reads one position twice, on the three models of `PEG10`, or passes over one, on the seven of `OAZ1`, `OAZ2` and `OAZ3`. The model's edit marks that place for translation. It adds no base to the transcript and takes none out.

An HGVSc position names a base of the transcript by its count from the start codon or from the stop codon, so the edit moves no position. Ensembl VEP writes one of two numbers, by the kind of notation, and each includes the edit. A deletion whose `CDS_position` is 1247 is written `c.1248del`, while a one-base substitution is written at its `CDS_position`. Neither is the base's position in RefSeq's record.

Genome-wide ClinVar gives 24 rows past such an edit. For 7 of them ClinVar's own title of the record names the row's model, and each of the 7 titles is vep-rs's notation. Sections 1.5 and 2.2 quote two.

**What vep-rs does instead.** On a model whose status is `failed`, vep-rs counts each HGVSc position along the genome's exons: a coding position from the first base of the start codon, and a 3' UTR position from the last base of the stop codon. On the `PEG10` models a coding position past the edit is one below the row's `CDS_position`, and on the `OAZ` models one above. `cDNA_position` and `CDS_position` keep the cached numbering, which both engines print. The two intron offsets of the HGVSIntronOffset plugin follow the notation, and no other cell does. vep-rs prints the `REFSEQ_OFFSET` Ensembl VEP prints, since that column states the sum Ensembl VEP makes.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at commit `db7394a02523a6e519cb0d89f84b2889bac12737`, `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`):

- Lines 1516-1517: for an `NM_` or `XM_` model that carries an `_rna_edit` attribute, the offset is the sum `get_misalignment_offset` returns, the value a row prints as `REFSEQ_OFFSET`.
- Lines 1594-1624, `get_misalignment_offset`: line 1605 splits an edit's value into a first position, a last position and, where the edit writes bases, their text. Line 1606 passes over an edit whose span is as long as its text, and lines 1608-1610 over an edit that starts after the row's cDNA start.
- Line 1617 adds the length of the text of an edit of three fields, whatever span the text replaces. Line 1614 adds the negative length of the span of an edit of two fields.
- Lines 1519-1522: a one-base substitution that has a CDS position is numbered at the row's CDS start, a number of the cached map, which counts the edit.
- Lines 1524-1525: each end of any other notation is moved by the offset along the transcript and counted along the exons from there (`_get_cDNA_position`, lines 2766-2768 and 2787-2794).
- Lines 2837-2852 and 2862-2868: the count is made relative to the model's coding end or coding start, which the cache holds with the edit counted.

**What correcting it would change.** With every position counted along the exons, and the coding start and end found on the genome the same way, each notation would stand at the transcript's own position, as vep-rs writes it. The rows that would change lie past the own edit of the ten models. The intron rows there, and the 3' UTR rows of the `OAZ` models, already stand at that position. No row of a model the alignment corrected (`BAM_EDIT=OK`) is in question: its cached sequence, with its edits, is its transcript.

**The class it produces.** `refseq_failed_edit_model_hgvsc_offset`: an `HGVSc` on a row that carries `BAM_EDIT=FAILED` and the same `cDNA_position` in both outputs, past such an edit. Both notations name the same transcript in coding numbering, `c.`, with no intron offset. The candidate, vep-rs's side of a comparison, counts each coordinate along the genome's exons. The reference, Ensembl VEP's side, writes a one-base substitution that has a CDS position at that CDS position, and any other notation at the candidate's count plus the offset, between the cached coding positions. [Section 1.7](#17-how-the-comparator-treats-this-class) says what the format comparator reads for it.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a row past the own edit of such a model. The cuts are chromosome 21 records and hold no such model, and no golden corpus's cache holds one. The RefSeq caches of release 116 hold ten such models, the same ten on GRCh38 and GRCh37. Genome-wide ClinVar holds nine records past their edits on each assembly, all one-base substitutions, on 24 rows: 18 coding rows and 6 rows in the 3' UTR of the `PEG10` models.

**What is not shown.** Each of these bounds what the page shows:

- Every row of the tables was written under `--refseq --hgvs` on the RefSeq GRCh38 cache of release 116. On the RefSeq GRCh37 cache the nine ClinVar records were run, and on all 24 rows each engine writes the notation it writes on GRCh38. The made-up records were not run there, and the merged caches were not read.
- The made-up records lie on `NM_015068.3` of `PEG10`, on `NM_004152.3` of `OAZ1` and on the reverse-strand `NM_002537.3` of `OAZ2`. The three `PEG10` models lie on the forward strand, so no row shows a reverse-strand model whose edit adds a position.
- No failed model carries more than one edit, an edit that adds bases between two positions, or an edit that writes a span as a shorter text of one base or more. For such an edit the class follows the source lines above, and no row shows it.
- A one-base deletion in the 3' UTR of `NM_015068.3` is not the class. At both places tried, Ensembl VEP's notation stands two positions above vep-rs's, and Ensembl VEP's row alone carries `HGVS_OFFSET=1`, the length of a 3' shift. One of the two positions is this entry's; this page does not trace the other. On `NM_004152.3` both engines write such a deletion alike.
- Every record of the tables lies more than two bases inside its exon. Ensembl VEP moves a position by the offset before it finds the exon, so within the offset's distance of an exon's edge the moved position leaves the exon (lines 2783 and 2834). On `NM_004152.3` a deletion of the first base of an exon past the edit is `c.293-1del` in Ensembl VEP and `c.293del` in vep-rs. A substitution in the last two bases of `NM_015068.3` has no `HGVSc` in Ensembl VEP. The class takes no such row.
- Line 1519 tests the record, not the allele. A one-base substitution allele of a made-up record that also holds an insertion takes the offset's path in Ensembl VEP: `c.1179` against vep-rs's `c.1177` on `NM_015068.3`. The class leaves that row a difference.
- A row that holds the position an edit adds is not the class. [An insertion at the twice-read position of a failed-edit RefSeq model](refseq-failed-edit-model-insertion-term.md) describes those of the `PEG10` models. An edit that takes a position out has no cached position of its own, and the row of the first base past it is a row of the class. A substitution of the base `NM_004152.3` passes over was run once: both engines write the same `HGVSc`, and the row's `Consequence` differs, which this page does not trace.
- The statement for an edit of any span and text is read from lines 1614 and 1617 and the cached coordinates. The rows show two edits: one position written as two bases, and one position taken out.

## 1. Positions past the edit of `NM_015068.3`

**In short.** `NM_015068.3` (gene `PEG10`, forward strand, `BAM_EDIT=FAILED`) carries one edit, which writes cDNA 1436 as two bases. RefSeq's record holds that position once and reads it twice. For ClinVar's record `7 94667457 T A`, in the 3' UTR, Ensembl VEP 116.3 writes `c.*1776` and vep-rs `c.*1775`. For ClinVar's coding record `7 94664733 A C`, Ensembl VEP writes `c.1178` and vep-rs `c.1177`. The base each notation names differs too, by another entry's cause.

### 1.1 What a reader needs to know first

A RefSeq model whose alignment correction failed prints `BAM_EDIT=FAILED` and keeps the edit it came with. [An insertion at the twice-read position of a failed-edit RefSeq model](refseq-failed-edit-model-insertion-term.md) describes these models: the coding sequence of `PEG10` reads one position twice, and the cache writes the double reading as one extra position.

The edit is the `_rna_edit` attribute `1436 1436` with two bases. The cached coordinates hold it: cDNA 1436 has no genomic base, and every cached position past it is numbered one above its count along the exons. The cached coding sequence runs from cDNA 480 to cDNA 2606, and its last base is `7:94665682`.

RefSeq's record of `NM_015068.3` gives the coding sequence as `join(480..1436,1436..2605)` with the qualifier `ribosomal_slippage`. The record holds the twice-read position once, so its coding sequence ends at position 2605, the count along the exons of `7:94665682`.

On a row past the edit Ensembl VEP also reads the reference base one position further along the transcript. The row's `USED_REF`, the reference allele it was annotated against, shows it, and the first section of [A coding transcript's cDNA under `--use_transcript_ref`](use-transcript-ref-cache-cdna.md#1-a-refseq-model-whose-alignment-correction-failed) documents it. That moves the base a notation names and not its position.

### 1.2 The records

The two records are lines 1,698,357 and 1,698,359 of the ClinVar GRCh38 file, variations 4134666 and 4340897, run under `--refseq --hgvs` with the GRCh38 reference FASTA:

```
7	94664733	4134666	A	C	.	.	.
7	94667457	4340897	T	A	.	.	.
```

### 1.3 What each engine writes

The row of the second record on `NM_015068.3`, `3_prime_UTR_variant` in both engines:

| Cell | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `cDNA_position` | 4381 | 4381 |
| `BAM_EDIT` | `FAILED` | `FAILED` |
| `REFSEQ_OFFSET` | 2 | 2 |
| `USED_REF` | `G` | `T` |
| `HGVSc` | `NM_015068.3:c.*1776G>A` | `NM_015068.3:c.*1775T>A` |

The rows on `NM_001184961.1` and `NM_001172437.2` carry the same two `HGVSc` numbers, `*1776` and `*1775`. ClinVar's next record, `7 94667481 C G`, is written `c.*1800` by Ensembl VEP and `c.*1799` by vep-rs on all three models.

The row of the first record on `NM_015068.3`, `missense_variant` in both engines:

| Cell | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `cDNA_position` | 1657 | 1657 |
| `CDS_position` | 1178 | 1178 |
| `USED_REF` | `G` | `A` |
| `HGVSc` | `NM_015068.3:c.1178G>C` | `NM_015068.3:c.1177A>C` |

On `NM_001184961.1` the same record is `c.1280` in Ensembl VEP and `c.1279` in vep-rs, and on `NM_001172437.2` `c.1406` and `c.1405`. The file's line between the two records, `7 94665299 G T`, is a coding record too: `c.1744` in Ensembl VEP and `c.1743` in vep-rs on `NM_015068.3`.

Four made-up records were run on both engines, each on `NM_015068.3`. Each is written where the reference bases beside it differ from it and from each other, so that no 3' shift moves it:

| Record | `CDS_position` | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| a coding deletion before the edit | 744 | `c.744del` | `c.744del` |
| a coding deletion past the edit | 1247 | `c.1248del` | `c.1246del` |
| a coding insertion past the edit | 1345-1346 | `c.1346_1347insA` | `c.1344_1345insA` |
| a 3' UTR insertion | none | `c.*1920_*1921insC` | `c.*1919_*1920insC` |

`cDNA_position` and `CDS_position` are equal between the engines on every row. Before the edit the notations are equal. Past it, Ensembl VEP's coding substitution stands at the row's `CDS_position`, one position above vep-rs's. Its coding deletion and insertion stand one above the row's `CDS_position`, two above vep-rs's. In these tables its 3' UTR notations stand one above vep-rs's.

Three made-up records of two bases past the edit were run as well. A coding substitution reads `c.1448_1449delinsAG` against vep-rs's `c.1446_1447delinsAG`, a coding deletion `c.1547_1548del` against `c.1545_1546del`, and a 3' UTR substitution `c.*2037_*2038delinsCA` against `c.*2036_*2037delinsCA`. On the two coding ones `Consequence` differs as well, which this page does not trace. The reproduction does not hold the three.

### 1.4 Why Ensembl VEP writes that

The 3' UTR record's count along the exons is 4380; the row's `cDNA_position`, on the cached numbering, is 4381. `get_misalignment_offset` gives its sum for the model's one edit. The edit has three fields and its span, one position, is not as long as its text, two bases, so line 1617 adds the text's length, 2. Ensembl VEP moves the record's position on by those 2 and counts along the exons from there: 4382.

The model's coding end is the cached cDNA 2606, which counts the edit's extra position. Lines 2848-2852 subtract it: 4382 less 2606 is 1776, written `*1776`.

The coding deletion's count along the exons is 1725. With the offset it is 1727, and counted from the coding start, cDNA 480, it is `c.1248`. The coding substitution takes the other path: lines 1519-1522 write its CDS start, 1178, the cached cDNA 1657 counted from cDNA 480. Each number counts the edit, the first as two positions and the second as one.

### 1.5 Why vep-rs's row is the correct one

The last coding base is `7:94665682`, and the three reference bases that end there are a stop codon. The 3' UTR record lies in the same exon, 1,775 bases past it. The number after `*` is that count, 1775, and ClinVar titles the record `c.*1775T>A` on `NM_015068.3`.

The first base of the start codon is `7:94663557`, in the same exon as the coding record. The record's base is the 1,177th base from it, `c.1177`. It is position 1656 of RefSeq's record, whose coding sequence starts at position 480, which gives the same number. The deletion's base is the 1,246th.

The residue is not in question. Both engines write `Protein_position` 393 on the coding record's row: the protein is read through the join, which reads the shared position twice, and `CDS_position` 1178 counts that reading.

### 1.6 Try it yourself

The reproduction `scripts/divergences/reproductions/refseq_failed_edit_model_hgvsc_offset.toml` runs both engines on thirteen records. Six lie on `NM_015068.3`: the two records of 1.2 and the four made-up records of 1.3. Seven lie on `NM_004152.3`: the records of section 2.2. To run vep-rs on the two records of 1.2, save their lines as `records.vcf` under a VCF header and run:

```
vep -i records.vcf -o out.txt --offline \
  --json_cache <json_cache> --species homo_sapiens --assembly GRCh38 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --refseq --hgvs
```

`<json_cache>` is the conversion of the RefSeq GRCh38 cache of release 116. The reproduction copies the two region files of that cache that hold the two models and converts them with `tools/vep-storable-to-json`. On the rows of the two models it compares `Consequence`, `cDNA_position`, `CDS_position`, `BAM_EDIT`, `REFSEQ_OFFSET` and `HGVSc`. Its record holds each engine's values, and sections 1.3 and 2.2 give their `HGVSc`.

### 1.7 How the comparator treats this class

The format comparator `tools/vep-formats-compare` holds the class `refseq_failed_edit_model_hgvsc_offset`. It reads the candidate's cache, given as `--json-cache`, for the edits of the row's transcript that change the cDNA's length, and for the model's cached coding start and end. Without the cache it sets nothing aside.

The class reads vep-rs's coordinates as counts along the genome's exons, between the coding start and end counted the same way. For a one-base substitution on a row with one CDS position, Ensembl VEP's coordinate must be that CDS position, and vep-rs's the coding position of the same base counted along the genome's exons. For any other notation, each of Ensembl VEP's coordinates must be vep-rs's count plus the offset summed over the edits before the row, read between the cached coding start and end.

Each notation must also stand at its row. A one-base substitution stands at the row's cDNA position. Any other notation stands within the row's `HGVS_OFFSET` and its allele's length of the row's positions.

What follows the coordinates is the same text in both. For a substitution whose base differs as well, the comparison must also name `use_transcript_ref_cache_cdna`, which reads the base.

Any other distance, a row before the edit, a row that holds a position an edit adds or runs across a base an edit takes out, and a notation with an intron offset stay differences. On the nine ClinVar records the two classes together set aside the 24 `HGVSc` cells, on each assembly, in the default format and in JSON. No golden corpus's cache holds such a model, so no corpus manifest names the class.

## 2. The models whose edit takes a position out

**In short.** Seven models of `OAZ1`, `OAZ2` and `OAZ3` carry an edit that takes one cDNA position out. RefSeq's record holds that position and its coding join passes over it. Past the edit Ensembl VEP's coding positions stand one below vep-rs's, and the two engines' 3' UTR positions are equal.

### 2.1 The model

The seven models are `NM_004152.3` and `NM_001301020.1` of `OAZ1`, `NM_002537.3` and `NM_001301302.1` of `OAZ2`, and `NM_016178.2`, `NM_001134939.1` and `NM_001301371.1` of `OAZ3`.

`NM_004152.3` (forward strand, `BAM_EDIT=FAILED`) carries the `_rna_edit` attribute `318 318` with no text, which takes cDNA 318 out. The cached map gives the genomic base `19:2271443` no cDNA position, so every cached position past it is numbered one below its count along the exons. The cached coding sequence runs from cDNA 114 to cDNA 800.

RefSeq's record gives the coding sequence as `join(114..317,319..801)` with the qualifier `ribosomal_slippage`. Position 318 is a base of the record, and the count along the exons is the record's.

### 2.2 What each engine writes

ClinVar's record `19 2273083 A T` on GRCh38, variation 2648969, has `CDS_position` 662 on `NM_004152.3` in both engines. Ensembl VEP writes `NM_004152.3:c.662G>T` and vep-rs `NM_004152.3:c.663A>T`. ClinVar titles the record `c.663A>T` on that model. As in section 1, the base Ensembl VEP names is the one it read one position off, here one position back along the transcript.

Six made-up records on `NM_004152.3`, written as those of section 1.3, were run on both engines:

| Record | `CDS_position` | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| a coding deletion before the edit | 165 | `c.165del` | `c.165del` |
| a coding deletion past the edit | 312 | `c.312del` | `c.313del` |
| a coding insertion past the edit | 360-361 | `c.360_361insG` | `c.361_362insG` |
| a coding substitution past the edit | 410 | `c.410G>A` | `c.411G>A` |
| an intron substitution past the edit | none | `c.465+12G>A` | `c.465+12G>A` |
| a 3' UTR insertion | none | `c.*92_*93insA` | `c.*92_*93insA` |

Three more lie on `NM_002537.3`, on the reverse strand, whose edit takes cDNA 330 out. A coding deletion at `CDS_position` 327 reads `c.327del` against vep-rs's `c.328del`. A coding insertion at `CDS_position` 275-276 reads `c.275_276insG` against `c.276_277insG`, and a 3' UTR insertion `c.*203_*204insG` in both. The reproduction holds the seven records on `NM_004152.3` and not these three.

### 2.3 Why the rows differ

The edit has two fields, so line 1614 adds the negative length of its span: the offset is -1, and Ensembl VEP writes a coding position past the edit at the count less one. Its one-base substitution stands at the cached CDS start, which is one below the count too. In the 3' UTR the offset and the cached coding end, one below the count of the last coding base, cancel. In an intron the row has no cDNA start, the offset is 0, and both engines count from the exon's end along the exons.

vep-rs's coding position is the base's position in RefSeq's record less the 113 positions before the start codon. The base passed over is one of the record's bases, so a base past it is one further from the start codon than its place in the translation.
