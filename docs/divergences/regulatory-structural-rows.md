# A structural variant's regulatory and motif rows: their order and the `--pick` row

One behaviour of Ensembl VEP 116.3's regulatory and motif rows for a structural variant is not reproduced: it writes them in the order the input reaches each 1 Mb region, which also sets the `--pick` row; vep-rs writes one order for every input. It is a documented class, which no comparison tool applies: no committed dataset carries a row of it. [Intended divergences from Ensembl VEP](../intended-divergences.md) lists it; its [How to read the counts](../intended-divergences.md#how-to-read-the-counts) defines the terms. The rows below are in the default output format, against the Ensembl 116 GRCh38 cache of chromosome 21 and the GRCh38 regulatory GFF dataset's three files. The rows themselves and their other fields are the same on both engines ([Regulatory and motif features](../regulatory.md)).

The class carries its tier under that page's rule: tier 1 is a shape the chromosome 21 cuts carry with no flag, tier 2 one that needs a flag or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts are over the cuts' 15,646,481 records, the 152 golden corpora and their 511 flag sets.

## 1. A structural variant's rows in the order its input reaches each 1 Mb region

**The defect in one sentence.** Ensembl VEP writes a structural variant's regulatory and motif rows region by region, in the order the records of its buffer first reach each 1 Mb region of the chromosome. So a record placed before the variant, a breakend whose mate lies ahead of it, or an earlier buffer in a run without `--fork` puts a later region's rows first, and under `--pick` a different row is kept. A breakend whose mate lies in an earlier region puts its own region's rows first by itself.

**Why this is a defect.** The order of an unchanged variant's rows, and the row `--pick` keeps for it, change with the variant's neighbours in the input, with the input's sort, and with `--fork` and `--buffer_size`: a property of the run, not of the variant and the annotation. The deletion of 1.1 alone gets `ENSR21_94D2C` first from Ensembl VEP; behind a breakend whose mate lies in the next region it gets `ENSR21_94DGM` first, with the same three files and flags. An order that changes with unrelated records is not a property of the features.

**What vep-rs does instead.** A `--regulatory_gff` run's rows come source by source (`file`, `emars`, `motifs`) and in file order within a source. A cache's rows come in the cache's order, region by region ascending ([VEP default and tab](../output/vep-and-tab.md#row-order)). Both orders are the same for every input, buffer size and thread count, and on a full `--pick` tie the first row in that order is kept.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/AnnotationSource.pm`, lines 119-134: the features of every region already in the source's memory cache first, then the uncached regions' features; line 136 and lines 392-411: `clean_cache` keeps the current buffer's regions alone. Line 235: the buffer's records in input order; lines 190-200: each record's 1 Mb regions ascending, a region once; lines 250-257: a breakend's mate regions after its own.
- `Bio/EnsEMBL/VEP/AnnotationSource/File/RegFeat.pm`, lines 245-269: the regions read from the file in the order given; lines 291-329: file order within a region; lines 628-640: a record read in two regions kept at its first.
- `Bio/EnsEMBL/VEP/AnnotationType/RegFeat.pm`, line 135: one overlap pushed per feature met; `Bio/EnsEMBL/Variation/StructuralVariationFeature.pm`, lines 495-509: the lists read as pushed, regulatory features before motifs.
- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 724, 789 and 808: under `--pick` every non-transcript row scores alike, the sort per category is stable, and the first of a full tie is kept.
- `Bio/EnsEMBL/VEP/Runner.pm`, lines 387-392 and 470-530: under `--fork` each chunk of records is annotated in a child process, so the memory cache dies with the chunk and the walk is per chunk. `Bio/EnsEMBL/VEP/InputBuffer.pm`, lines 216-221 and 240-250, with `Constants.pm` lines 348-349: 100 records more than 100 bases out of order are accepted.

**What correcting it would change.** The rows of a structural variant would come in one order whatever the input, as vep-rs's do; Ensembl VEP's own order, and its `--pick` row, for a deletion behind a breakend, an unsorted record or an earlier buffer would change.

**The class it produces.** `regulatory_structural_row_order`: the `RegulatoryFeature` and `MotifFeature` rows of one structural record and allele under `--regulatory` or `--regulatory_gff` (JSON `regulatory_feature_consequences` and `motif_feature_consequences` entries) whose fields agree between the engines and whose order within a feature type differs. Under `--pick`, `--pick_allele`, `--flag_pick` and `--flag_pick_allele` the row kept or flagged for a full tie among such rows: the engines keep different rows of one equal-ranked set. The class is documented only. The sibling class [`regulatory_gff_motif_row_order`](regulatory-gff-motif-rows.md#1-two-motif-rows-of-one-record-in-the-order-of-a-read-counter) is the same region walk read through a sequence variant's motif counter.

The golden harness and the format comparator check that the entries both engines share come in one order, outside a flag set marked `unordered`, and compare a kept or flagged row as any other row. They report a differing order as an order difference, a differing flagged row as a field difference, and a differing kept row as two entries, each on one side alone.

Every committed comparison orders such rows alike: the 43 structural records of the GRCh38 regulatory GFF dataset are position-sorted, none spans a 1 Mb region boundary and none is a breakend, so both engines reach their regions in ascending order, which is file order. The one record of the shape in the datasets run under `--regulatory`, a breakend whose mate lies two regions ahead, meets a regulatory feature at its own position alone, so one region's row is written on both engines.

**How often the shape occurs.** Tier 2: the shape is a symbolic structural record spanning a 1 Mb region boundary, or a breakend whose own-chromosome mate lies in another region. A mate in an earlier region puts the breakend's own region first by itself; otherwise the input's order, a neighbour's mate, `--fork` and the buffer fill decide whether a later region is reached first. 13 of the cuts' 15,646,481 records carry the shape, every one a symbolic record of the 1000 Genomes cut spanning a boundary; 8 corpora carry it, and 60 flag sets run `--regulatory` or `--regulatory_gff`.

### 1.1 A 43 kb deletion over the 6 Mb boundary behind a breakend whose mate lies past it

**In short.** The records `21 4500000 bnd_mate_in_region6 N N[21:6500000[ . . SVTYPE=BND` and `21 5972900 sv_del_over_6mb N <DEL> . . SVTYPE=DEL;END=6016000`, position-sorted, with the GRCh38 regulatory GFF dataset's three files (`tests/golden/116/GRCh38-regulatory-gff/gff/`) and FASTA, under `--regulatory_gff file=<regulatory>,motifs=<motifs>,emars=<EMARs>`. The deletion covers the enhancer `ENSR21_94D2C` (5972926-5973418) and its EMAR `ENSR00001793255` in the region 5000001-6000000, and the enhancer `ENSR21_94DGM` (6015540-6015943), its EMAR `ENSR00001793257` and the motifs `ENSM00000575084` and `ENSM00000575085` in the region 6000001-7000000. Both engines write the same seven rows for it: four `RegulatoryFeature` rows with `regulatory_region_ablation,regulatory_region_variant`, two `MotifFeature` rows with `TFBS_ablation,TF_binding_site_variant`, and the intergenic row, every field equal. The regulatory rows differ in order alone:

| Input | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| the deletion alone | `ENSR21_94D2C`, `ENSR21_94DGM`, `ENSR00001793255`, `ENSR00001793257` | the same |
| the deletion, then `21 6500000 snv_in_region6 A C` | the same | the same |
| the breakend, then the deletion | `ENSR21_94DGM`, `ENSR21_94D2C`, `ENSR00001793257`, `ENSR00001793255` | `ENSR21_94D2C`, `ENSR21_94DGM`, `ENSR00001793255`, `ENSR00001793257` |
| `21 6500000 snv_in_region6 A C`, then the deletion | `ENSR21_94DGM`, `ENSR21_94D2C`, `ENSR00001793257`, `ENSR00001793255` | `ENSR21_94D2C`, `ENSR21_94DGM`, `ENSR00001793255`, `ENSR00001793257` |

The breakend's mate at 6500000 puts the later region on the buffer's list before the deletion adds the earlier one, so Ensembl VEP reads the later region's features first in each file. The sequence record at 6500000 does the same from an input one record out of order, which Ensembl VEP accepts. The motif rows lie in one region and keep their order in every run. With `--buffer_size 1` and no `--fork`, the breakend's buffer leaves the later region in the memory cache, which the deletion's buffer lists first: the same flipped order. With `--buffer_size 1 --fork 1` the forked child starts with no cache, and Ensembl VEP writes the file order. vep-rs writes the file order in every run.

### 1.2 The same deletion under `--pick` with the regulatory file alone

**In short.** The two records of 1.1 under `--regulatory_gff file=<regulatory>` and `--pick`. The deletion's two `RegulatoryFeature` rows tie on every pick category, so the first listed is kept. Ensembl VEP keeps `ENSR21_94DGM` (`OverlapBP=404`); vep-rs keeps `ENSR21_94D2C` (`OverlapBP=493`). For the deletion alone both engines keep `ENSR21_94D2C`. With all three files and `--pick` both engines keep the motif `ENSM00000575084`, whose `TFBS_ablation` outranks the regulatory terms, and the runs agree.

### 1.3 A cache's rows: a deletion over the 15 Mb boundary behind a record past it

**In short.** The records `21 15500000 snv_in_region15 A C` and `21 14981353 sv_cn0_like_1kg N <DEL> . . SVTYPE=DEL;END=15016529`, in that order, under `--regulatory` on the Ensembl 116 GRCh38 cache of chromosome 21, with no GFF file. The deletion covers two regulatory features and five motifs in the region 14000001-15000000 and five regulatory features and eight motifs in the region 15000001-16000000. Ensembl VEP writes the later region's rows first, and vep-rs writes the earlier region's first, in ascending order.

Ensembl VEP writes the later region's rows first: `ENSR21_97XQC`, `ENSR21_97XRC`, `ENSR21_B4ZDW`, `ENSR21_97XT7`, `ENSR21_B4ZGH`, then `ENSR21_8FW9H`, `ENSR21_97XP4`; the motifs `ENSM00000577283` to `ENSM00000577290`, then `ENSM00000577278` to `ENSM00000577282`. With the deletion first in the input it writes the earlier region's rows first, `ENSR21_8FW9H`, `ENSR21_97XP4`, `ENSR21_97XQC` and so on, and the motifs from `ENSM00000577278` upwards. vep-rs writes that ascending order for both inputs, every field equal to Ensembl VEP's sorted-input rows.

### 1.4 A breakend alone whose mate lies in an earlier region

**In short.** The one record `21 6500000 bnd_mate_earlier N N[21:5973000[ . . SVTYPE=BND` with the three files and FASTA of 1.1. Its own position lies in the enhancer `ENSR21_CMF6` (6498466-6500510) and its EMAR `ENSR00001793359`, its mate in `ENSR21_94D2C` and `ENSR00001793255` of the earlier region. Ensembl VEP adds the mate's region to the buffer's list after the record's own, so it writes `ENSR21_CMF6`, `ENSR21_94D2C`, `ENSR00001793359`, `ENSR00001793255`; vep-rs writes `ENSR21_94D2C`, `ENSR21_CMF6`, `ENSR00001793255`, `ENSR00001793359`, the file order.

The mate-side rows carry `feature_truncation` in Ensembl VEP and `regulatory_region_variant` in vep-rs, which is [Defect E](defect-e.md). The same breakend written from the other end, `21 5973000 bnd_mate_later N N[21:6500000[`, reaches its regions in ascending order, and both engines write `ENSR21_94D2C`, `ENSR21_CMF6`, `ENSR00001793255`, `ENSR00001793359`.
