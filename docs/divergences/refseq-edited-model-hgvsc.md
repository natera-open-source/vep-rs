# The HGVSc of a RefSeq model whose alignment inserted or deleted bases

**The defect in one sentence.** On a RefSeq transcript model whose alignment to the genome inserted or deleted bases, Ensembl VEP numbers the HGVSc from the unedited exons. An intronic distance counts from the unedited exon end, and an exonic position on an `NM_` or `XM_` model is slid along the genome by the edit length before it is read. An `NR_` or `XR_` model is never corrected, nor is an Ensembl model whose `_rna_edit` inserts or deletes bases (section 2). A deletion that opens before bases the model lacks and closes after them is written through them, its end too high by their count. vep-rs numbers every HGVSc position against the edited model, the transcript the identifier names.

**Why this is a defect.** HGVS coordinates on a transcript refer to the transcript's own sequence, which for a RefSeq model is the edited one. The cache records each difference from the genome as an `_rna_edit` attribute and stores the model's coordinate mapper with the edits folded in, so the cDNA numbering of every exon past an insertion or deletion is shifted by the edit's length. The cDNA, CDS and protein position columns of the same row use that mapper and agree between the engines; only the HGVSc is read through a mapper rebuilt without the edits.

**What vep-rs does instead.** vep-rs numbers every HGVSc position against the edited model, the transcript the identifier names, on RefSeq, merged and Ensembl-source caches alike.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116, `ensembl` at release/116):

- `Bio/EnsEMBL/VEP/Pipeline/DumpVEP/Dumper/Core.pm`, lines 164 and 169, inside `clean_transcript` (line 150): the dumper deletes the transcript's `edits_enabled` flag and its `transcript_mapper`.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1524-1525 and 2755, `_get_cDNA_position`: each end of the allele is turned into a cDNA coordinate by walking the exons and reading each exon's cDNA start and end through `Exon::cdna_start` (`Bio/EnsEMBL/Exon.pm`, lines 383-410). That method rebuilds a mapper at run time, and with the flag gone `TranscriptMapper.pm`'s `_load_mapper` folds no edit in.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1517 and 1594-1624, `get_misalignment_offset`: for an `NM_` or `XM_` model alone, the length changes of the indel edits at or before the allele's cDNA start are summed (the `REFSEQ_OFFSET` column) and added to the genomic position of both ends. An allele whose start has no cDNA position, an intronic one, counts no edit (lines 1608-1610).
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2787-2794, `_get_cDNA_position`: an end inside an exon is the exon's cDNA start plus the end's genomic distance from the exon's start, which counts bases the model lacks. The offset is one sum: lines 1608-1610 count an edit only where the allele's cDNA start lies at or past it, and lines 1524-1525 add the sum to both ends.

**What correcting it would change.** The class vanishes: every HGVSc on an edited model carries the edited numbering, as the cDNA and CDS position columns of the same row already do.

**The class it produces.** `refseq_edited_model_hgvsc_offset`, applied by the format comparator to the `HGVSc` of a row carrying `BAM_EDIT=OK` whose two notations name the same reference and differ in their coordinate tokens alone. Where both state how many positions a pair of coordinates covers, the counts must agree. A model of an Ensembl-source cache whose `_rna_edit` inserts or deletes bases prints no `BAM_EDIT`, so the comparator's own rule leaves its row a difference and a corpus manifest alone sets it aside. For the golden harness each corpus documents its rows of the class in its manifest; [section 1.7](#17-how-the-comparator-treats-this-class) says how they are named. The concordance figures score consequence sets and never see the field, so nothing is set aside in them.

A deletion over a model's first base that runs across an edit is a row of the entry as well: Ensembl VEP numbers its end on the unedited exon. The two notations then cover different numbers of positions, and the format comparator's rule leaves the pair a difference. [Section 3](#3-a-deletion-over-the-first-base-that-runs-across-an-edit) shows the row.

A deletion that opens inside the coding sequence, before bases the model lacks, and closes after them is a row of the entry too. Its row states both counts of positions, and the comparator accepts the pair as the one exception to its count test ([section 4](#4-a-deletion-across-the-nine-bases-col18a1s-models-lack)).

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion over a donor site downstream of a nine-base edit

**In short.** A 31-base deletion over the donor site of a COL18A1 exon lies on two RefSeq models whose alignment deleted nine bases upstream. Ensembl VEP writes `c.3546_3553+23del`, eight exonic bases none of which is deleted and an intronic distance nine short; vep-rs writes `c.3544+2_3544+32del`, the span the genome has, the same shape the unedited sibling model `NM_001379500.1` gets from both engines.

### 1.1 What a reader needs to know first

A converted RefSeq or merged cache carries the attributes Ensembl VEP writes about each model's alignment to the genome ([RefSeq and merged cache columns](../output/refseq-and-merged.md)). A model the alignment corrected prints `BAM_EDIT=OK` on its rows, and its `_rna_edit` attributes name each insertion or deletion as `<start> <end> <bases>`. `REFSEQ_OFFSET` prints the summed length change of the edits upstream of the allele. Under `--hgvs` both engines shift an insertion or deletion 3' along the transcript before writing its notation, and `HGVS_OFFSET` prints the distance moved.

### 1.2 The record

The record is `21 45505272 rs2037126010 AGGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG A`, record 34 of `variants.vcf` in the `GRCh38-refseq` golden corpus (`tests/golden/116/GRCh38-refseq/`), a gnomAD v4.1 deletion of `21:45505273-45505303`:

```
21	45505272	rs2037126010	AGGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG	A	.	.	.
```

Two models carry a nine-base deletion edit upstream of it: `NM_030582.4` (`_rna_edit` `3387 3395`) and `NM_130444.3` (`_rna_edit` `4092 4100`). Both engines run the corpus's `hgvs` flag set, `--refseq --hgvs`, and shift the deletion seven bases 3' along the transcript, to `21:45505280-45505310` (`HGVS_OFFSET=7`).

### 1.3 What each engine writes

The two models' rows, default format. Ensembl VEP (the corpus's `expected_hgvs/default.txt.gz`):

```
rs2037126010	21:45505273-45505303	-	80781	NM_030582.4	Transcript	splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant	3548-?	3539-?	1180-?	-	-	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=-9;GIVEN_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;USED_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;BAM_EDIT=OK;HGVSc=NM_030582.4:c.3546_3553+23del;HGVS_OFFSET=7
rs2037126010	21:45505273-45505303	-	80781	NM_130444.3	Transcript	splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant	4253-?	4244-?	1415-?	-	-	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=-9;GIVEN_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;USED_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;BAM_EDIT=OK;HGVSc=NM_130444.3:c.4251_4258+23del;HGVS_OFFSET=7
```

vep-rs:

```
rs2037126010	21:45505273-45505303	-	80781	NM_030582.4	Transcript	splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant	3548-?	3539-?	1180-?	-	-	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=-9;GIVEN_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;USED_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;BAM_EDIT=OK;HGVSc=NM_030582.4:c.3544+2_3544+32del;HGVS_OFFSET=7
rs2037126010	21:45505273-45505303	-	80781	NM_130444.3	Transcript	splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant	4253-?	4244-?	1415-?	-	-	-	IMPACT=HIGH;STRAND=1;REFSEQ_OFFSET=-9;GIVEN_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;USED_REF=GGCAGAGTAAGTCAGTGGGGAGTGGGCCCCG;BAM_EDIT=OK;HGVSc=NM_130444.3:c.4249+2_4249+32del;HGVS_OFFSET=7
```

Every column but the HGVSc agrees, `REFSEQ_OFFSET` included: vep-rs keeps Ensembl VEP's value of the column and applies no offset to the notation.

### 1.4 Why Ensembl VEP writes that

`hgvs_transcript` turns each end of the shifted allele into a cDNA coordinate with `_get_cDNA_position`, which walks the transcript's exons and reads each exon's cDNA start and end through `Exon::cdna_start`. That method maps through the transcript's own mapper, built at run time, and the cache dumper deleted the transcript's edit flag and stored mapper, so the rebuilt mapper folds no `_rna_edit` in and numbers every exon by exon lengths alone. The edited numbering survives only in the cached mapper the cDNA, CDS and protein position columns read, which is why those columns agree while the HGVSc does not.

For an `NM_` or `XM_` model, `get_misalignment_offset` then sums the length changes of the indel edits at or before the allele's cDNA start, -9 here, and adds the sum to the genomic position of both ends, a move along the genome by the edit length. An allele whose start has no cDNA position, an intronic one, counts no edit, so an intronic distance is anchored on the unedited exon end. Here the -9 slides the start from `21:45505280` into the exon at `21:45505271`, numbers the exon end 3553 where the edited model's last exon base is 3544, and counts the intronic distance 23 for 32.

Source (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116, `ensembl` at release/116):

- `Bio/EnsEMBL/VEP/Pipeline/DumpVEP/Dumper/Core.pm`, lines 164 and 169: the flag and the mapper deleted.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1524-1525, 2755 and 2881-2905: the cDNA positions read through the rebuilt mapper.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1517 and 1594-1624: the offset, and lines 1608-1610, the intronic start that counts no edit.
- `Bio/EnsEMBL/Exon.pm`, lines 383-410; `Bio/EnsEMBL/Transcript.pm`, lines 2268-2271: the run-time mapper.

### 1.5 Why vep-rs's row is the correct one

The model's exon is `21:45505134-45505278`, edited cDNA 3409-3553 and CDS 3400-3544; the nine bases the edit `3387 3395` removed lie upstream, so the unedited numbering of the exon runs nine higher. The shifted deletion removes `21:45505280-45505310`, every base of it past the exon end `21:45505278` by 2 to 32 bases: `c.3544+2_3544+32del`. Ensembl VEP's `c.3546_3553+23del` names eight exonic bases, c.3546 to c.3553, none of which is deleted, and an intronic run nine bases short. The unedited sibling model `NM_001379500.1` of the same gene, which shares the exon's genomic span, is written by both engines as `c.3013+2_3013+32del`, the shape vep-rs gives the edited models. Both notations describe 31 bases; only vep-rs's places them where the genome has them.

### 1.6 Try it yourself

```
vep -i tests/golden/116/GRCh38-refseq/variants.vcf -o out.txt --offline --cache \
  --dir_cache tests/golden/116/GRCh38-refseq/cache_root --cache_version 116 \
  --species homo_sapiens --assembly GRCh38 --refseq --hgvs \
  --fasta tests/golden/116/GRCh38-refseq/reference.fa.gz --force_overwrite --no_stats --quiet
grep rs2037126010 out.txt | grep -o 'HGVSc=[^;]*'
```

| Model | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `NM_030582.4` | `c.3546_3553+23del` | `c.3544+2_3544+32del` |
| `NM_130444.3` | `c.4251_4258+23del` | `c.4249+2_4249+32del` |
| `NM_001379500.1` (unedited) | `c.3013+2_3013+32del` | `c.3013+2_3013+32del` |

On the `GRCh37-refseq` corpus's cache under `--refseq --hgvs`, `rs1232890876` at `21:38639535` on `NM_001331018.1`, `NM_001331021.1` and `NM_001331022.1` reads `c.56+4A>G` from Ensembl VEP and `c.57+4A>G` from vep-rs, the models carrying the insertion edit `647 646 G`.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) applies `refseq_edited_model_hgvsc_offset` when a comparison case names it in `masks`. A differing `HGVSc` (the JSON `hgvsc`) on a row carrying `BAM_EDIT=OK` is accepted when both notations name the same reference, the text before the colon, and are equal once every cDNA coordinate token is replaced by one marker. So `c.3546_3553+23del` and `c.3544+2_3544+32del` compare equal, and a changed allele text or reference version does not. Each accepted row is counted once; the RefSeq comparison cases under `--hgvs` expect 2 rows on GRCh38 and 3 on GRCh37.

The rule asks one thing more: a pair of coordinates covers as many positions in one notation as in the other, where both state the count. A pair states its count as two plain positions or as two offsets from one exon boundary. A pair that mixes the two, such as `c.3546_3553+23`, states none and is not tested. The row of [section 3](#3-a-deletion-over-the-first-base-that-runs-across-an-edit) fails that test and stays a difference. No field of it holds the length of the edits inside the span, which the corpus tool reads from the cache. The row of [section 4](#4-a-deletion-across-the-nine-bases-col18a1s-models-lack) fails it as well, but `GIVEN_REF` and the cDNA position state both counts there. [Section 4.7](#47-how-the-comparator-treats-this-row) says when such a pair is accepted.

The golden harness reads the documented rows from each corpus's manifest. The corpus tool names the class from the corpus cache, not from the `BAM_EDIT` column. The row's transcript is filed there with an `_rna_edit` that inserts or deletes bases. The two notations name the same reference and differ in their coordinates alone, and a number of one lies a running sum of the edits' length changes from a number of the other.

A documented row of a RefSeq or merged cache carries `BAM_EDIT=OK`, as the two rows above do in the GRCh38 RefSeq corpus's `hgvs` flag set. A documented row of an Ensembl-source cache carries no `BAM_EDIT`, as the row of [section 2](#2-a-substitution-in-an-intron-of-an-ensembl-model-whose-edit-deletes-a-base) does. On every documented row only the `HGVSc` changes.

## 2. A substitution in an intron of an Ensembl model whose edit deletes a base

**In short.** A substitution lies in the intron of the yellow fever mosquito model `AAEL025507-RA`, 763 bases before its last exon. An `_rna_edit` deletes one base of the model's first exon. Ensembl VEP writes `c.310-763G>T`, anchored on the exon start the unedited exons give; vep-rs writes `c.309-763G>T`, the exon's first base on the model's own sequence. The cache is an Ensembl one, so neither row carries `BAM_EDIT`.

### 2.1 What a reader needs to know first

An `_rna_edit` is not a RefSeq attribute alone. An Ensembl-source cache files one on a model whose transcript sequence differs from the genome, as `<start> <end> <bases>` in the cDNA numbering of the spliced exons. Such a cache holds no alignment status and no `NM_` or `XM_` model, so its rows print neither `BAM_EDIT` nor `REFSEQ_OFFSET`, and nothing in a row says that the model carries an edit.

### 2.2 The record

The record is line 86 of `variants.vcf` in the `AaegL5-aedes` golden corpus (`tests/golden/116/AaegL5-aedes/`), a designed substitution on an unplaced scaffold of the Aedes aegypti cache of Ensembl Metazoa:

```
NIGP01001870	3784	intron_12	C	A	.	.	.
```

`AAEL025507-RA` is a reverse-strand model of two exons, `NIGP01001870:6067-6439` and `NIGP01001870:3018-3021`, with the coding start at cDNA base 65. Its one `_rna_edit`, `302 302` with no bases, deletes cDNA base 302, genomic base 6138 of the first exon. Both engines run the corpus's `hgvs-sift` flag set, `--hgvs --sift b`; the `everything` flag set gives the same two notations.

### 2.3 What each engine writes

The model's row, default format. Ensembl VEP (the corpus's `expected_hgvs-sift/default.txt.gz`):

```
intron_12	NIGP01001870:3784	A	AAEL025507	AAEL025507-RA	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=AAEL025507-RA:c.310-763G>T
```

vep-rs:

```
intron_12	NIGP01001870:3784	A	AAEL025507	AAEL025507-RA	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=AAEL025507-RA:c.309-763G>T
```

Every column but the `HGVSc` agrees.

### 2.4 Why Ensembl VEP writes that

The position is intronic, so `_get_cDNA_position` names it by the nearer exon: the last one, whose first base `NIGP01001870:3021` is 763 bases away, against 2,283 to the first exon. It reads that exon's cDNA start through the mapper rebuilt at run time, which folds no edit in. The first exon's 373 genomic bases put the last exon at cDNA 374, and from the coding start at cDNA 65 that is `c.310`. The model is neither `NM_` nor `XM_`, so no offset is computed, and an intronic allele counts no edit in an offset on any model. The source lines are those of [section 1.4](#14-why-ensembl-vep-writes-that).

### 2.5 Why vep-rs's row is the correct one

The model's own sequence holds 376 bases: the first exon's 373 less the deleted one, then the last exon's four. Its coding sequence, bases 65 to 376, is 312 bases, 103 codons and a stop codon, and translates to the 103-residue protein the cache stores. The unedited exons hold 313 bases from the same start, not a whole number of codons, and reach no stop codon. The edited model is therefore the transcript the identifier names. The last exon opens at base 373 of that sequence, coding base 309, and the record lies 763 bases before it: `c.309-763G>T`. Ensembl VEP's `c.310` is the exon's second base, the first base of the stop codon.

### 2.6 Try it yourself

```
vep -i tests/golden/116/AaegL5-aedes/variants.vcf -o out.txt --offline \
  --json_cache tests/golden/116/AaegL5-aedes/json_cache \
  --species aedes_aegypti_lvpagwg --assembly AaegL5 --hgvs --sift b \
  --fasta tests/golden/116/AaegL5-aedes/reference.fa.gz \
  --force_overwrite --no_stats --quiet
grep intron_12 out.txt | grep -o 'HGVSc=[^;]*'
```

| Model | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `AAEL025507-RA` | `c.310-763G>T` | `c.309-763G>T` |

### 2.7 How the comparator treats this row

The row carries no `BAM_EDIT`, so the format comparator's own rule does not accept it: run without the corpus manifest, the comparator reports the `HGVSc` as a difference. Given the manifest, the comparator sets the row aside as the manifest documents it when the comparison case names `refseq_edited_model_hgvsc_offset` in `masks`. It counts the row once in each of the corpus's two flag sets under `--hgvs`, `everything` and `hgvs-sift`. The golden harness reads the same manifest entries.

## 3. A deletion over the first base that runs across an edit

**In short.** On GRCh37, a record at chromosome 21, position 46020492, with a 182-base reference allele and a one-base alternate, deletes 21:46020493-46020673: four bases before the first base of the forward-strand model `NM_198689.3` and 177 genomic bases of its one exon. The deletion runs across the `_rna_edit` `160 174`, fifteen genomic bases that the model lacks. Under `--refseq --hgvs --shift_hgvs 0` both engines write `coding_sequence_variant,5_prime_UTR_variant`, `cDNA_position` `?-162`, `CDS_position` `?-137` and `BAM_EDIT=OK`. Ensembl VEP writes `c.-25_152del` and vep-rs `c.-25_137del`, whose end is the row's own `CDS_position` end.

The steps on the row, by line of `TranscriptVariationAllele.pm`:

| Step | Lines | On this record |
| --- | --- | --- |
| The written span is mapped onto the transcript slice, and its start is clamped to 1 | 2712, 2732-2733 | (-3, 177) becomes (1, 177) |
| The notation is typed over the clamped span | 1474-1483 | A `del` over 1 to 177 |
| `get_misalignment_offset` sums the edits at or before the allele's cDNA start. A span that starts before the transcript has none, which is read as 0, so no edit is counted | 1516-1517, 1594-1624, 1608-1610 | The offset is 0, and the row prints no `REFSEQ_OFFSET` |
| Each end is numbered along the unedited exon | 1524-1525, 2755 | 1 is `c.-25`, and 177 less the 25 bases of the 5' UTR is `c.152` |

The model pairs cDNA 1 to 159 with 21:46020497-46020655 and cDNA 160 onwards with 21:46020671 onwards. The fifteen genomic bases between them have no place on the model. Inside the transcript the record removes 21:46020497-46020673, the model's bases 1 to 162. The coding sequence starts at cDNA 26, so the removed bases are c.-25 to c.137: `c.-25_137del`. Ensembl VEP's `c.152` is the model's base 177, fifteen bases past the last base the record removes.

A model whose edit adds a base shows the opposite length. On GRCh38, a record at chromosome 21, position 37267458, with a 655-base reference allele, deletes 21:37267459-37268112 over the first base of the reverse-strand model `NM_001331018.1`. The model holds one base at cDNA 647 that the genome lacks. Ensembl VEP's end, `c.-166`, is the model's base 650, one base short of the last base the record removes.

| Record and run | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `NM_198689.3`, `--refseq --hgvs --shift_hgvs 0` | `c.-25_152del` | `c.-25_137del` |
| The same record under `--refseq --hgvs`, `HGVS_OFFSET=2` | `c.-23_154del` | `c.-25_139del` |
| `NM_001331018.1`, `--merged --hgvs --shift_hgvs 0`, `cDNA_position` `?-651` | `c.-815_-166del` | `c.-815_-165del` |
| The same record under `--merged --hgvs`, `HGVS_OFFSET=-1` | `c.-814_-165del` | `c.-815_-164del` |

In the two runs with a shift the start differs as well. That is the clamped start of [`hgvs_shifted_deletion_clamped_span`](hgvs-shifted-deletion-clamped-span.md).

| What a comparison reads | On these rows |
| --- | --- |
| The format comparator's rule | It holds each pair of coordinates to one count of positions in both notations. Its one exception, the deletion of [section 4](#4-a-deletion-across-the-nine-bases-col18a1s-models-lack), opens at a coding position on a row whose cDNA position has both sides. `c.-25_152del` covers 177 positions and `c.-25_137del` 162, and this row's cDNA position is `?-162`, so the pair stays a difference. |
| What the row cannot show | The length of the edits inside the span, 15 bases on `NM_198689.3` and 1 on `NM_001331018.1`. The offset is 0, so the row prints no `REFSEQ_OFFSET`. |
| The corpus tool | Its test reads the model's edits from the corpus cache. It asks that a number of one notation lie a running sum of the edits' length changes from a number of the other. These ends lie 15 and 1 apart, each such a sum on its model, so the test holds: the tool names the class for such a row. No golden corpus carries either record, so no manifest documents these rows. |

## 4. A deletion across the nine bases COL18A1's models lack

**In short.** The cache's alignment of `NM_030582.4` and `NM_130444.3` leaves out the nine genomic bases `21:45504526-45504534`. The gnomAD v4.1 deletion `rs1179352929` removes seventeen genomic bases that hold all nine, so it removes eight bases of each model. Ensembl VEP writes `c.3372_3388del`, seventeen positions; vep-rs writes `c.3372_3379del`, the eight the models hold. Both open at the same position, and Ensembl VEP's closes nine positions later.

### 4.1 What a reader needs to know first

A deletion edit names bases the genome has and the model lacks. They lie inside one exon of the model's genomic structure, and no cDNA position belongs to them. A genomic deletion that opens before them and closes after them removes fewer bases from the model than from the genome. The row's `cDNA_position` and `CDS_position` span the model's bases, and `GIVEN_REF` spells the genomic ones.

### 4.2 The record

The record is `rs1179352929` of gnomAD v4.1, a GRCh38 deletion of `21:45504519-45504535`:

```
21	45504518	rs1179352929	CCAGGCCCCCCAGGCCCA	C	.	.	.
```

The models' edits are `3387 3395` (`NM_030582.4`) and `4092 4100` (`NM_130444.3`), each nine positions in the numbering of the unedited exons, inside the exon `21:45504416-45504556`. Both engines run `--refseq --hgvs` and shift the deletion one base 3', to `21:45504520-45504536` (`HGVS_OFFSET=1`).

### 4.3 What each engine writes

The two models' rows, default format, as the commands of [section 4.6](#46-try-it-yourself) print them. Ensembl VEP:

```
rs1179352929	21:45504519-45504535	-	80781	NM_030582.4	Transcript	frameshift_variant	3380-3387	3371-3378	1124-1126	PGP/X	cCAGGCCCC/c	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=AGGCCCCCCAGGCCCAC;BAM_EDIT=OK;HGVSc=NM_030582.4:c.3372_3388del;HGVSp=NP_085059.2:p.Gly1125TrpfsTer136;HGVS_OFFSET=1
rs1179352929	21:45504519-45504535	-	80781	NM_130444.3	Transcript	frameshift_variant	4085-4092	4076-4083	1359-1361	PGP/X	cCAGGCCCC/c	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=AGGCCCCCCAGGCCCAC;BAM_EDIT=OK;HGVSc=NM_130444.3:c.4077_4093del;HGVSp=NP_569711.2:p.Gly1360TrpfsTer136;HGVS_OFFSET=1
```

vep-rs:

```
rs1179352929	21:45504519-45504535	-	80781	NM_030582.4	Transcript	frameshift_variant	3380-3387	3371-3378	1124-1126	PGP/X	cCAGGCCCA/c	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=AGGCCCCCCAGGCCCAC;BAM_EDIT=OK;HGVSc=NM_030582.4:c.3372_3379del;HGVSp=NP_085059.2:p.Gly1125TrpfsTer136;HGVS_OFFSET=1
rs1179352929	21:45504519-45504535	-	80781	NM_130444.3	Transcript	frameshift_variant	4085-4092	4076-4083	1359-1361	PGP/X	cCAGGCCCA/c	-	IMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=AGGCCCCCCAGGCCCAC;BAM_EDIT=OK;HGVSc=NM_130444.3:c.4077_4084del;HGVSp=NP_569711.2:p.Gly1360TrpfsTer136;HGVS_OFFSET=1
```

The `Codons` differ as well; they are the entry [`refseq_edited_model_deleted_span_codons`](refseq-edited-model-deleted-span-codons.md). Every other column agrees. Both engines place the deletion at eight cDNA positions, 3380-3387 and 4085-4092, and neither prints `REFSEQ_OFFSET`.

### 4.4 Why Ensembl VEP writes that

`_get_cDNA_position` numbers an end inside an exon as the exon's cDNA start plus the end's genomic distance from the exon's start. The exon opens before the edit, so the start of the deletion, before the nine bases, is numbered as the models number it. The end lies after them, and its genomic distance counts all nine.

`get_misalignment_offset` exists to take an edit's length off. It counts an edit only where the allele's cDNA start lies at or past the edit's start, and its one sum is added to both ends. The start it tests is the row's own, before the shift: cDNA 3380 and 4085. Both lie before the edits at 3387 and 4092, so the offset is 0 and the end stays nine too high.

Source (`ensembl-variation` at release/116): `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2787-2794, the genomic distance; lines 1608-1610 and 1524-1525, the one offset.

### 4.5 Why vep-rs's row is the correct one

The sequence the record leaves is the model's exons spliced from the genome without the seventeen bases. It is eight bases shorter than `NM_030582.4`'s own sequence, which never held the nine. A deletion of eight model bases gives exactly that sequence when it opens at c.3371 or at c.3372 and at no other position, and HGVS names the most 3' one: `c.3372_3379del`.

Ensembl VEP's `c.3372_3388del` names seventeen bases of the model. The model without them is nine bases shorter than the sequence the record leaves. Ensembl VEP's own row disagrees with it: its `cDNA_position` and `CDS_position` span eight positions.

The unedited sibling model `NM_001379500.1` holds all seventeen bases. Both engines write `c.2832_2848del` for it, seventeen positions over a `CDS_position` of seventeen.

### 4.6 Try it yourself

```
printf '21\t45504518\trs1179352929\t%s\tC\t.\t.\t.\n' \
  CCAGGCCCCCCAGGCCCA > record.vcf
vep -i record.vcf -o out.txt --offline --cache \
  --dir_cache tests/golden/116/GRCh38-refseq/cache_root \
  --cache_version 116 --species homo_sapiens --assembly GRCh38 \
  --refseq --hgvs --force_overwrite --no_stats --quiet \
  --fasta tests/golden/116/GRCh38-refseq/reference.fa.gz
grep -o 'HGVSc=NM_[^;]*' out.txt
```

| Record and model | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `rs1179352929`, `NM_030582.4` | `c.3372_3388del` | `c.3372_3379del` |
| `rs1179352929`, `NM_130444.3` | `c.4077_4093del` | `c.4077_4084del` |
| ClinVar `1338712`, `NM_030582.4` | `c.3376_3393del` | `c.3376_3384del` |
| ClinVar `1338712`, `NM_130444.3` | `c.4081_4098del` | `c.4081_4089del` |
| `rs1179352929`, `NM_001379500.1` (unedited) | `c.2832_2848del` | `c.2832_2848del` |

ClinVar variation `1338712` deletes the eighteen genomic bases `21:45504521-45504538`: the nine the models lack and nine of the models' own. Both engines shift it three bases 3'.

### 4.7 How the comparator treats this row

The two notations cover different numbers of positions, so the row fails the count test of [section 1.7](#17-how-the-comparator-treats-this-class). The format comparator still accepts the pair under `refseq_edited_model_hgvsc_offset` where the rest of that rule holds (`BAM_EDIT=OK`, one reference) and all of these do:

- the row's `STRAND` is 1;
- each notation is one `del` between two coding positions, with no `-`, `*` or intron offset, and the two open at the same position;
- Ensembl VEP's notation covers as many positions as `GIVEN_REF` has bases;
- vep-rs's notation covers as many positions as the row's cDNA position spans, both of its sides stated, and fewer than Ensembl VEP's.

A duplication, an insertion or a `delins` across such bases and a pair that opens at different positions stay differences. So do a reverse-strand row, a position outside the coding sequence and a deletion over bases the model holds and the genome lacks.

The test reads the row, and the row does not say where the lacked bases lie. No row of vep-rs needs it to. Where a shifted span opens or closes on a lacked base, vep-rs names that end from the nearest modelled base with an offset ([HGVSc on an edited model](../output/refseq-and-merged.md#hgvsc-on-an-edited-model)). Such a notation states no count of positions, so the rule of section 1.7 does not test the pair.

No golden corpus carries the record, so the row counts of section 1.7 do not include it.
