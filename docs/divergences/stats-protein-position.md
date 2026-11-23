# The protein-position chart of the statistics report under `--hgvs`

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP's statistics report counts a coding variant's position in the protein from the HGVS-shifted codon while its `Protein_position` column prints the unshifted one. The `[Position in protein]` chart of a `--hgvs` run therefore differs from the chart of the same run without `--hgvs`, although the rows' positions are the same. vep-rs counts the position its `Protein_position` column prints, so its chart is the same in both runs.

**Why this is a defect.** The chart summarises the `Protein_position` column: each coding row's position as a percentile of its protein's length. Both the column and the chart should count the same position. Ensembl VEP's HGVSp routine recomputes the protein position with the 3' shift it uses for the notation and leaves the shifted value cached, and the statistics logger reads that cached value, so the chart moves with a flag that changes no row's position. The run's own column and its own chart disagree.

**What vep-rs does instead.** vep-rs counts the position its `Protein_position` column prints, under `--hgvs` and without it, so its chart is the same in both runs and equals Ensembl VEP's chart for the run without `--hgvs`.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116, `ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1724-1733, `hgvs_protein`: when an allele carries a shift, the TranscriptVariation's `translation_start` and `translation_end` are deleted and recomputed with the shift's offset, and the recomputed values stay cached.
- `Bio/EnsEMBL/VEP/Stats.pm`, line 359, `log_BaseTranscriptVariationAllele`: the chart counts `translation_start` as the object holds it.

**What correcting it would change.** The class vanishes: the chart of a `--hgvs` run is the chart of the run without it, the one the `Protein_position` column supports.

**The class it produces.** `stats_protein_position_chart`: the `[Position in protein]` section of the statistics report of a `--hgvs` run. The concordance figures and the output-format comparisons never read the statistics report, so nothing is set aside in them; the golden statistics corpus runs without `--hgvs`, where the two charts agree.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. One input, two charts

**In short.** On the input of the `GRCh37-hgvs` golden corpus, Ensembl VEP's `[Position in protein]` chart reads `171 75 79 95 90 129 65 78 113 158` across its ten deciles without `--hgvs` and `166 76 79 95 87 126 66 80 112 145` with it. Every other section of the report and every row of the output keep their values. vep-rs's chart is the first vector in both runs.

### 1.1 What a reader needs to know first

`--stats_text` writes a plain-text statistics report beside the output (`--stats_file` names it; the HTML report is the default form). Its `[Position in protein]` section counts the coding rows by the decile of the protein their position falls in: a row at residue 50 of a 500-residue protein is counted under `00-10%`. The count is over rows, one per transcript consequence, so one record contributes as many counts as it has coding rows.

Under `--hgvs`, an insertion or deletion is shifted to the most 3' position along the transcript at which the altered sequence is the same before its notation is written. The shift changes the HGVSc and HGVSp strings; the `Protein_position` column reports the position of the allele as written in the input.

### 1.2 The record

The input is `variants.vcf` of the `GRCh37-hgvs` golden corpus (`tests/golden/116/GRCh37-hgvs/`), 293 ClinVar GRCh37 chromosome 21 records chosen for the HGVS flags, annotated against the corpus's cache and FASTA slice. Of the run's 2,818 rows, 319 are coding rows whose allele the HGVS shift moves (`HGVS_OFFSET` set and a protein position printed). Both engines run twice: with `--stats_text` alone and with `--stats_text --hgvs`.

### 1.3 What each engine writes

The `[Position in protein]` section of each report. Ensembl VEP, without `--hgvs`:

```
[Position in protein]
00-10%	171
10-20%	75
20-30%	79
30-40%	95
40-50%	90
50-60%	129
60-70%	65
70-80%	78
80-90%	113
90-100%	158
```

Ensembl VEP, with `--hgvs`:

```
[Position in protein]
00-10%	166
10-20%	76
20-30%	79
30-40%	95
40-50%	87
50-60%	126
60-70%	66
70-80%	80
80-90%	112
90-100%	145
```

vep-rs writes the first section in both runs. Every other section of the two Ensembl VEP reports is identical between them, and the `Protein_position` column of the `--hgvs` output equals the column of the run without the flag on every row.

### 1.4 Why Ensembl VEP writes that

The statistics logger counts a row's protein position from the TranscriptVariation's `translation_start`. When the row's HGVSp is written, `hgvs_protein` deletes the object's `translation_start` and `translation_end` and recomputes them with the allele's shift offset, so that the notation names the shifted codon, and the recomputed values stay on the object. The logger runs after the notation and reads the shifted position; the `Protein_position` column was filled from the unshifted one. Without `--hgvs` no notation is written, nothing is recomputed, and the logger reads the position the column printed.

The two charts also count different numbers of rows, 1,053 and 1,032, the sums of the two vectors above. A shifted start that maps to no coding position, an insertion or deletion carried past the stop codon into the 3' UTR or across an exon boundary into the intron, is stored as undefined (`BaseTranscriptVariation.pm`, line 375). The logger counts only a row whose `translation_start` is defined (`Stats.pm`, line 359), while a start the shift moves from an intron into the coding sequence is counted for the first time. The 21 rows are the balance of the two.

Source (`ensembl-variation` at release/116, `ensembl-vep` at release/116.3):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1724-1733: the cached positions deleted and recomputed with the shift.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 365-380: a shifted coordinate the mapper returns as a gap stored as an undefined position.
- `Bio/EnsEMBL/VEP/Stats.pm`, lines 359-360: the chart's read of `translation_start`, a defined value alone, and its decile.

### 1.5 Why vep-rs's row is the correct one

The chart is a summary of the output. Counting each coding row's printed `Protein_position` against its protein's length gives `171 75 79 95 90 129 65 78 113 158`, the vector vep-rs prints and the vector Ensembl VEP prints without the flag. Ensembl VEP's `Protein_position` column is the same in both of its runs on every row. The second vector counts the shifted rows at positions no column of the output states and leaves out the rows whose shifted start maps to no coding position. It therefore disagrees with the first on how many rows lie in each decile of the same proteins and on how many rows there are, which no notation flag can change.

### 1.6 Try it yourself

Run the corpus input through both engines twice, with and without `--hgvs`, writing the text report. vep-rs's command is below; Ensembl VEP takes the same flags with the cache and FASTA options of [Reproducing these examples](reproducing.md):

```
vep -i variants.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --fasta <corpus>/reference.fa.gz --species homo_sapiens --assembly GRCh37 \
  --stats_text --hgvs
```

| Run | Ensembl VEP `00-10%` to `90-100%` | vep-rs |
| --- | --- | --- |
| without `--hgvs` | `171 75 79 95 90 129 65 78 113 158` | the same |
| with `--hgvs` | `166 76 79 95 87 126 66 80 112 145` | `171 75 79 95 90 129 65 78 113 158` |

The `Protein_position` column is identical across the four outputs, and the first vector is its distribution over the proteins' lengths.

### 1.7 How the comparator treats this class

No comparator of this repository reads the statistics report of a `--hgvs` run. The format comparator (`tools/vep-formats-compare`) compares the annotation outputs and the warnings files, and the golden harness compares the statistics text of the `GRCh38-statistics` corpus, whose runs carry no `--hgvs`, so the two charts agree there. The class `stats_protein_position_chart` names the section a comparison of a `--hgvs` run's statistics report sets aside, with every other section compared as written.

## 2. An eleventh row for a position past the end of a short peptide

**In short.** On `21 39428463 tiny_cds_end T C`, both engines write the same coding row on `ENST00000451131` (LCA5L, reverse strand, a `cds_end_NF` model of 10 residues): `incomplete_terminal_codon_variant` with `coding_sequence_variant`, `Protein_position` 11. Ensembl VEP's `[Position in protein]` chart prints that row under an eleventh key, `110-120%`, and `0` under `90-100%`; vep-rs counts it under `90-100%`, the only decile a position at or past a protein's end can belong to.

The shape is tier 2 of [Intended divergences from Ensembl VEP](../intended-divergences.md): 20 of the 116,165 transcripts of the committed caches (172 per million) end in a partial codon after at most ten residues, in 7 corpora, and 6 flag sets run `--stats_text`. No comparison of this repository reads such a report, so nothing is set aside.

### 2.1 The record

The record lies in the trailing partial codon of `ENST00000451131`, whose CDS end is not annotated (`cds_end_NF`): CDS position 31 of 32 bases, residue 11 of a 10-residue peptide. Sixteen other coding rows of the record lie at residue 11 of longer proteins and fall in the first four deciles on both engines. Both engines run `--stats_text` on a GRCh38 cache; the `GRCh38-individual` golden corpus's cache (`tests/golden/116/GRCh38-individual/json_cache`) holds the transcript.

### 2.2 What each engine writes

| Row | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `90-100%` | `0` | `1` |
| `110-120%` | `1`, printed between `00-10%` and `10-20%` | no row |

### 2.3 Why Ensembl VEP writes that

The logger keys a coding row by `int(10 * position / length)` (`Bio/EnsEMBL/VEP/Stats.pm`, line 360), 11 here. When the report is built, the key `10`, the last residue, is folded into `9`, the keys `0` to `9` are filled, and every key present is labelled `<k>0-<k+1>0%` (lines 461-467). A key of 11 therefore keeps a row of its own, sorted by its first character beside `10-20%` (lines 769-776).

### 2.4 Why vep-rs's chart is the correct one

The chart's ten rows are the deciles of each protein's length, and the fold of the key `10` states where a position at the protein's end belongs: the top decile. A position past the end, the partial codon of an unfinished model, has no decile of its own and can only fall there too; an eleventh decile of a protein does not exist. vep-rs counts every key of 10 or more under `90-100%`, so its chart has ten rows on every input.
