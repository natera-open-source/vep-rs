# `--regulatory_gff` motif rows: their order and `MOTIF_POS` without a FASTA

Two behaviours of Ensembl VEP 116.2's `--regulatory_gff` motif rows are not reproduced: it orders a record's rows by a read counter and omits `MOTIF_POS` without a FASTA, where vep-rs orders them by `ID` and writes `MOTIF_POS` either way. Each is a documented class, which no comparison tool applies: no committed dataset carries a row of it. [Intended divergences from Ensembl VEP](../intended-divergences.md) lists it; its [How to read the counts](../intended-divergences.md#how-to-read-the-counts) defines the terms. The rows below are in the default output format, against the Ensembl 116 GRCh38 cache of chromosome 21 and the GRCh38 regulatory GFF dataset's three files. The motif rows themselves and their other fields are the same on both engines ([Regulatory features from GFF files](../regulatory.md#regulatory-features-from-gff-files)).

Each class carries its tier under that page's rule: tier 1 is a shape the chromosome 21 population cuts carry with no flag, tier 2 one that needs a flag or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts are over the cuts' 15,646,481 records, the 111 golden corpora and their 392 flag sets.

## 1. Two motif rows of one record in the order of a read counter

**The defect in one sentence.** Ensembl VEP orders a record's motif rows by the text of a counter it gives each record of a `--regulatory_gff` motif file as it reads it, so the order of two motif rows at one position depends on the records before them in the input.

**Why this is a defect.** A motif read from a GFF file has no database identifier, so Ensembl VEP gives it a counter advanced by one per record read and sorts the rows by that counter as text, `10` before `9`. The counter follows the run's read sequence, each 1 Mb region of the file in the order the input's records first reach it. So the same two motifs at one position come as `ENSM00000574661, ENSM00000574660` for one input, and the reverse for an input with one more record earlier on the chromosome. An order that changes with unrelated records is not a property of the features.

**What vep-rs does instead.** vep-rs orders a GFF motif's rows by its stable id, the `ID` attribute: the key the regulatory rows already use, and the same for every input.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/RegFeat.pm`, line 452, `_build_motif_feature`: `dbID => ++$self->{_rf_dbID}`, one counter per source, advanced per record read; lines 245-266: each region of the file read in the order given. `Bio/EnsEMBL/VEP/AnnotationSource.pm`, lines 119-136 and 235: the regions gathered from the buffer's records in input order.
- `Bio/EnsEMBL/Variation/MotifFeatureVariation.pm`, lines 70-76: the row keyed by the motif's `dbID`; `VariationFeature.pm`, line 937: stored under that key, and line 750: the rows returned in `sort keys` order, a text sort. `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 477-516: no later sort of the motif rows.

**What correcting it would change.** The motif rows of a record would come in one order whatever the input, as vep-rs's do; Ensembl VEP's own order of two motifs at one position would change for some inputs.

**The class it produces.** `regulatory_gff_motif_row_order`: two or more `MotifFeature` rows of one record and allele under `--regulatory_gff` (JSON `motif_feature_consequences` entries) whose fields agree between the engines and whose order differs. The golden harness and the format comparator check that the entries both engines share come in one order, outside a flag set marked `unordered`, so they report such rows as an order difference. The class is documented only. Every committed comparison orders such rows alike: the 6 sequence-variant records of the GRCh38 regulatory GFF dataset that overlap two or more motifs, one of them under two alleles, read their motifs in identifier order on both engines.

**How often the shape occurs.** Tier 2: the shape is two motif records of the file sharing a base, which a variant inside both reaches under `--regulatory_gff`; the run's read order decides whether the counters' text sorts unlike the identifiers. 4,425 of the committed motif file's 13,577 records share a base with the record before them, and 6 flag sets run `--regulatory_gff`.

### 1.1 `21:5014060 G>A` inside the motif pair `ENSM00000574660` and `ENSM00000574661`

**In short.** The record `21 5014060 motif_9_and_10 G A` with the GRCh38 regulatory GFF dataset's three files (`tests/golden/116/GRCh38-regulatory-gff/gff/`) and FASTA, under `--regulatory_gff file=<regulatory>,motifs=<motifs>,emars=<EMARs>`. Both engines write the same five rows: the EMAR and enhancer rows, a `MotifFeature` row each for `ENSM00000574660` (`STRAND=1`, `MOTIF_POS=3`) and `ENSM00000574661` (`STRAND=-1`, `MOTIF_POS=13`), both `TF_binding_site_variant` on matrix `ENSPFM0363`, and the intergenic row. vep-rs writes `ENSM00000574660` first. Ensembl VEP writes `ENSM00000574661` first: the two records are the ninth and tenth motif records read for their region, and `10` sorts before `9` as text.

With `21 12500000 lead_bin12 G A` placed before the record in the input, Ensembl VEP reads the region holding position 12500000 first, the pair take counters 16 and 17, and it writes `ENSM00000574660` first, as vep-rs does in both runs. The JSON and VCF outputs swap the same two entries.

## 2. `MOTIF_POS` omitted from a motif row without a FASTA

**The defect in one sentence.** Without a FASTA, Ensembl VEP writes a `--regulatory_gff` motif row with no `MOTIF_POS`, and four `Use of uninitialized value` warnings, because the stand-in slice it builds for the chromosome has no strand.

**Why this is a defect.** `MOTIF_POS` is the variant's start as a position in the matrix, computed from the record's and the motif's coordinates alone; no sequence is read, and the regulatory source's own comment says its consequences are purely coordinate-based. Without a FASTA or a database the source builds a stand-in slice of 2,000,000,000 bases for the chromosome and leaves its strand undefined, so the motif's genomic start is taken from the slice's end, the position falls outside the matrix and the field is dropped. The same record with a FASTA prints the values vep-rs prints.

**What vep-rs does instead.** vep-rs computes `MOTIF_POS` from the two coordinates and the motif's strand and writes it with or without a FASTA.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `ensembl-variation` and `ensembl` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource/File/RegFeat.pm`, lines 528-549, `_get_or_fake_slice`: the stand-in slice built when `get_slice` returns nothing, with `start`, `end` and `seq_region_length` and no `strand`; line 125, `$FAKE_SLICE_LENGTH`; lines 505-508, the comment that no sequence is required.
- `Bio/EnsEMBL/Feature.pm`, lines 1090-1095, `seq_region_start`: `$slice->strand() == 1` is false for an undefined strand (the warning), so the start is `$slice->end() - $self->end() + 1`.
- `Bio/EnsEMBL/Variation/MotifFeatureVariationAllele.pm`, lines 97-108, `motif_start`: the position from the two starts, flipped on the reverse strand, undefined outside 1 to the matrix length; `Bio/EnsEMBL/VEP/OutputFactory.pm`, line 1928: `MOTIF_POS` written only when defined.

**What correcting it would change.** The stand-in slice would carry a forward strand, and a `--regulatory_gff` motif row would print the same `MOTIF_POS` with or without a FASTA, as vep-rs's does.

**The class it produces.** `regulatory_gff_motif_pos_without_fasta`: a `MOTIF_POS` value (the JSON `motif_pos`) that vep-rs alone writes on a `--regulatory_gff` motif row in a run without a FASTA, and in `--warnings` mode the reference's `Use of uninitialized value in numeric eq (==)` lines naming `Feature.pm` line 1090. The class is documented only. Every committed `--regulatory_gff` comparison runs with its dataset's FASTA.

**How often the shape occurs.** Tier 3: the shape is a `--regulatory_gff` run naming a motif file with no FASTA in reach, which 0 of the 392 flag sets is; 6 flag sets run `--regulatory_gff`.

### 2.1 The record of 1.1 without `--fasta`

**In short.** The record of 1.1 run as there, without `--fasta`. vep-rs writes the five rows of 1.1 unchanged, `MOTIF_POS=3` and `MOTIF_POS=13` included. Ensembl VEP writes the same rows without `MOTIF_POS` on either motif row, and its warnings file carries, after the cache notice, four lines `Use of uninitialized value in numeric eq (==) at .../Feature.pm line 1090`. Every other field agrees.
