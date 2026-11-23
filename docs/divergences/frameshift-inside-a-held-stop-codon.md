# A frameshift inside a stop codon that the protein holds

**The defect in one sentence.** Some transcripts have a protein that holds a stop and names the residues after it, as a polymorphic pseudogene's does. An insertion or a deletion can change the length by no multiple of three and have its first base inside the codon of such a stop. Ensembl VEP 116.3 writes no `frameshift_variant` for it, and under `--hgvs` it names the stop unchanged or in a deletion-insertion. vep-rs writes `frameshift_variant` and names the first residue that changes.

**Why this is a defect.** The Sequence Ontology defines `frameshift_variant` as a variant that disrupts the translational reading frame because the number of nucleotides inserted or deleted is not a multiple of three. The protein of such a transcript goes on past the stop, and Ensembl VEP numbers and describes the residues after it on every other row. An allele of that length inside the stop's codon shifts the frame of each of them.

Ensembl VEP writes the frameshift itself when the same change is spelled one base on. Where the base before the stop codon is the same as the codon's first base, a record can delete either, or add a copy beside either, and describe one sequence. Section 1 shows two such pairs: one change and one HGVSc each, with two terms, two impacts and two protein notations from Ensembl VEP.

What Ensembl VEP writes in the frameshift's place is untrue of the sequence. `p.Ter339=` says that a stop is unchanged where its codon has become a glutamate codon. `inframe_deletion` stands on deletions of 4 to 56 bases whose length is no multiple of three.

**What vep-rs does instead.** vep-rs asks whether the first residue of the allele's codons is a stop that the protein holds with a residue named after it. Where it is, the allele's length decides the frameshift, as at any other codon. The cells that follow the term are then a frameshift's: `IMPACT`, the HGVSp, and the cells of the plugins that annotate a frameshift.

Under `--hgvs` the frameshift is named at the first residue that changes. Where the stop becomes a residue, the stop is the first residue that changes, and the notation is the extension both engines write for a frameshift whose first changed residue is a stop. Where the new frame leaves a stop codon in the stop's place, the protein goes on past it as before, and the notation names the first changed residue after it.

The stop codon that closes a coding sequence has no residue after it, and an allele that begins inside it is no frameshift in either engine: its row says whether that stop is lost or kept. The same holds for one that begins at a stop the protein holds as its last residue, with none named after it. A stop codon that the protein reads as a residue, a selenocysteine codon, is no stop of the protein, and both engines write the frameshift there.

One form inside a stop that the protein goes on past keeps Ensembl VEP's row. Where the allele's alternate codons are one whole codon that is a stop codon, and one or two bases more, both engines write `stop_retained_variant`. Section 2 counts those rows.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1541-1586, `frameshift`: line 1554 returns 0 where the reference residues of the allele's codons begin with a stop. Its comment gives the reason: an allele whose first affected base is the stop codon does not affect the reading frame. Lines 1556-1564, the length test, are then not reached.
- The same file, line 815: the in-frame tests get the allele's codons only where `frameshift` is 0. Lines 1126-1127 and line 1178 can then give `inframe_insertion` and `inframe_deletion`, and lines 1616-1638 give `coding_sequence_variant` where no other coding term holds.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2320-2365, `_get_fs_peptides`: lines 2354-2357 end the walk for a frameshift's first changed residue at a stop that both proteins hold, and write it unchanged. The walk's reference protein is the transcript's own, so a stop it holds ends the walk as the stop appended after it at line 2338 does.

**What correcting it would change.** The rows would read `frameshift_variant` with `IMPACT=HIGH`, their HGVSp would name the first changed residue, and the plugins that read the term would annotate them.

**The class it produces.** `frameshift_inside_a_held_stop_codon`: vep-rs's row holds `frameshift_variant` where Ensembl VEP's holds `coding_sequence_variant`, `inframe_deletion`, or `inframe_insertion` with `coding_sequence_variant`. Section 1.7 says how the comparisons recognise such a pair.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape needs a transcript whose protein holds a stop, and chromosome 21 has none, so the five record sets the index calls the cuts hold no such record. Section 2 counts the rows of the whole-genome population datasets.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. One base deleted from the first stop of a pseudogene's protein

**In short.** A record of 1000 Genomes Phase 3 and of the gnomAD v2.1.1 genomes deletes one base at 7:72209528 on GRCh37, the first base of codon 339 of ENST00000435769. Ensembl VEP 116.3 writes `coding_sequence_variant` with `IMPACT=MODIFIER` and `p.Ter339=`. vep-rs writes `frameshift_variant` with `IMPACT=HIGH` and `p.Ter339GluextTer34`, the row Ensembl VEP writes for the same deletion spelled one base on.

### 1.1 What a reader needs to know first

ENST00000435769 is a `polymorphic_pseudogene` on the reverse strand. Its coding sequence holds stop codons before the one that closes it, and the cache's protein of 669 residues holds each as a stop: 23 of them, the first at residue 339. Both engines number the residues after such a stop and describe changes to them.

A VCF record names an insertion or a deletion from the base before it. Where two neighbouring bases are the same, deleting either gives one sequence, and the record may name either place. Both engines keep the place the record names for the consequence terms. The HGVS 3' rule moves the change to its last place for the notations alone.

### 1.2 The records

```
7 72209527 base_lost_inside_the_stop_codon_7_72209527 CA C . . .
7 72209528 same_base_lost_written_before_the_stop_codon AA A . . .
7 72209527 base_added_inside_the_stop_codon C CA . . .
7 72209528 same_base_added_written_before_the_stop_codon A AA . . .
7 72209526 base_lost_and_the_stop_kept TC T . . .
```

The columns are tab-separated in the file. The transcript runs right to left on the chromosome: 72209529 is the last base of codon 338, and 72209528 and 72209527 are the first and second bases of codon 339, the stop codon. The bases at 72209528 and 72209529 are the same.

The first record is the dataset's. It deletes 72209528, the stop codon's first base. The second, written for this page, deletes 72209529 and leaves the same sequence. The third and the fourth add one copy of that base after 72209527 and after 72209528, one sequence again. The fifth deletes 72209527, the stop codon's second base.

Both engines run `--hgvs` on the release 116 GRCh37 cache with the GRCh37 reference.

### 1.3 What each engine writes

The rows of ENST00000435769. Ensembl VEP 116.3:

| Record | Consequence | IMPACT | HGVSp |
| --- | --- | --- | --- |
| `base_lost_inside_the_stop_codon_7_72209527` | `coding_sequence_variant` | `MODIFIER` | `p.Ter339=` |
| `same_base_lost_written_before_the_stop_codon` | `frameshift_variant` | `HIGH` | `p.Ter339GluextTer34` |
| `base_added_inside_the_stop_codon` | `coding_sequence_variant` | `MODIFIER` | `p.Ter339delinsLeuTer` |
| `same_base_added_written_before_the_stop_codon` | `frameshift_variant` | `HIGH` | `p.Ter339LeuextTer42` |
| `base_lost_and_the_stop_kept` | `coding_sequence_variant` | `MODIFIER` | `p.Ter339=` |

vep-rs:

| Record | Consequence | IMPACT | HGVSp |
| --- | --- | --- | --- |
| `base_lost_inside_the_stop_codon_7_72209527` | `frameshift_variant` | `HIGH` | `p.Ter339GluextTer34` |
| `same_base_lost_written_before_the_stop_codon` | `frameshift_variant` | `HIGH` | `p.Ter339GluextTer34` |
| `base_added_inside_the_stop_codon` | `frameshift_variant` | `HIGH` | `p.Ter339LeuextTer42` |
| `same_base_added_written_before_the_stop_codon` | `frameshift_variant` | `HIGH` | `p.Ter339LeuextTer42` |
| `base_lost_and_the_stop_kept` | `frameshift_variant` | `HIGH` | `p.Glu340SerfsTer34` |

Both engines write `c.1015del` for the first two records and `c.1015dup` for the next two, and `c.1016del` for the fifth. The position cells, `Amino_acids` and `Codons` are the same in both engines on every row.

### 1.4 Why Ensembl VEP writes that

The first record's codons are codon 339 alone. Its reference residue is the stop, so line 1554 of `frameshift` returns 0 and the length test is not made. The alternate codon has two bases, which Ensembl VEP reads as one unknown residue, `X`. For an allele whose alternate residues hold `X`, `stop_lost` and `stop_retained_variant` are asked of the stop codon that closes the coding sequence alone. `coding_sequence_variant` is the term left.

For the HGVSp the frameshift is asked again and is 0, so the notation comes from the two residues of the codon. The stop and the unknown residue are both written `Ter`, which reads as no change: `p.Ter339=`.

The second record's first base lies in codon 338, whose reference residue is no stop. Line 1554 does not apply, the length test gives a frameshift, and the walk finds residue 339 as the first that changes.

The fifth record gets `coding_sequence_variant` and `p.Ter339=` by the same lines as the first. Were its frameshift asked for, the walk for the first changed residue would begin at residue 339, find a stop in both proteins there and write `p.Ter339=` again, by lines 2354-2357.

### 1.5 Why vep-rs's row is the correct one

From the record, the reference and the exons alone: the deleted base is coding base 1015, the first of codon 339. Without it the next three bases spell a glutamate codon, and every codon after it is read one base on. The new frame meets its first stop at residue 373, 34 residues after the lost stop.

So the stop is the first residue that changes, the frame of the 330 residues after it shifts, and the change is the frameshift that `p.Ter339GluextTer34` names. The second record is the same sequence, and Ensembl VEP writes exactly that row for it. For the first it writes a term that names no change of the protein and a notation that says the stop is unchanged.

The fifth record leaves a stop codon at residue 339: the codon's two other bases and the next base spell one. The protein goes on past it as it does in the reference. Residue 340, a glutamate, is the first that changes, to a serine, and the new frame's first stop after it is residue 373, the 34th residue counted from the changed one: `p.Glu340SerfsTer34`. The stop of residue 339 stands before the changed residue, so the count is the one [The stop count of a frameshift after a stop codon inside the coding sequence](hgvsp-frameshift-count-past-a-stop-read-as-a-residue.md) describes.

### 1.6 Try it yourself

With the records of 1.2 saved as `records.vcf`, tab-separated, the command as vep-rs takes it, on a converted release 116 GRCh37 cache:

```
vep -i records.vcf -o out.txt --offline --hgvs \
  --assembly GRCh37 --json_cache <json_cache> \
  --fasta Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
grep ENST00000435769 out.txt
```

Ensembl VEP 116.3 takes the same records and `--hgvs` with the cache options of [Reproducing these examples](reproducing.md). The entry's reproduction, `scripts/divergences/reproductions/frameshift_inside_a_held_stop_codon.toml`, cuts the region file that holds the transcript from the GRCh37 cache, converts it and runs both engines on the five records.

### 1.7 How the comparisons treat this class

The three comparison tools read the two rows of a pair and nothing else. A pair is of the class where four things hold. vep-rs's terms hold `frameshift_variant` and Ensembl VEP's do not. Ensembl VEP's terms hold in its place `coding_sequence_variant`, `inframe_deletion`, or `inframe_insertion` with `coding_sequence_variant`, and every other term stands on both rows. The `Amino_acids` cells are equal, open with `*` and hold `X` after the `/`. The `Codons` cells are equal, and their two parts differ in length by no multiple of three.

The tuple comparison (`tools/vep-outputs-compare`) scores consequence terms and sets the pair's terms aside. The format comparator (`tools/vep-formats-compare`) sets aside the `Consequence`, an `IMPACT` that is the one each side's own terms give, and the HGVSp. It sets aside a `DownstreamProtein`, `ProteinLengthChange`, `NMD` or LoFTEE cell that vep-rs alone writes as well, since those plugins annotate a frameshift. The plugin tuple scorer (`tools/vep-plugin-pairs`) takes a row by the same test and sets aside the same plugin cells, under its class `plugin_value_follows_stop_codon_term_pair`.

A pair on which Ensembl VEP writes `stop_retained_variant`, `stop_lost` or `frameshift_variant` is not of the class, and any other cell that differs stays a difference charged to vep-rs.

## 2. Where real data holds such rows

**In short.** The four whole-genome population datasets hold 55 such rows, on 53 records and 18 transcripts; the genome-wide ClinVar files hold none. A record is counted in each dataset that holds it, and one record lies on three of the transcripts. On each, the row vep-rs writes is the one derived from the record, the reference and the exons.

| Dataset | Rows |
| --- | --- |
| 1000 Genomes Phase 3 (GRCh37) | 1 |
| gnomAD v2.1.1 genomes (GRCh37) | 2 |
| 1000 Genomes high coverage (GRCh38) | 1 |
| gnomAD v4.1 genomes (GRCh38) | 51 |

The transcripts are `polymorphic_pseudogene` on GRCh37 (3 rows) and `protein_coding_LoF` (40) or immunoglobulin D segments (12) on GRCh38. The alleles are 31 deletions of 1 to 56 bases and 24 insertions of 1 to 22.

| What the change does | Rows | Ensembl VEP 116.3 |
| --- | --- | --- |
| The stop becomes a residue | 42 | `coding_sequence_variant` |
| The stop becomes a residue | 6 | `inframe_deletion` |
| The stop codon is kept | 4 | `inframe_insertion,coding_sequence_variant` |
| The stop codon is kept | 2 | `inframe_deletion` |
| The stop codon is kept | 1 | `coding_sequence_variant` |

Beside either in-frame term Ensembl VEP's `IMPACT` is `MODERATE`, and its HGVSp is a deletion-insertion that ends in a stop.

vep-rs writes `frameshift_variant` on all 55. On the 48 rows of the first two lines its HGVSp is the extension named at the stop, with the count to the first stop the new frame meets through the 3' UTR. Eleven of them meet none there and keep `extTer?`. On the 7 rows whose stop codon is kept it names the first changed residue after the stop: in the frameshift form on 6, two of them with the unknown count, and on one as a residue that becomes a stop at once.

Where a stop that the protein holds stands before the first changed residue, the count is the one [The stop count of a frameshift after a stop codon inside the coding sequence](hgvsp-frameshift-count-past-a-stop-read-as-a-residue.md) describes. On a row whose stop codon is kept, that stop itself stands before it.

Eight more pairs of the gnomAD v4.1 genomes begin inside such a stop codon with an allele of such a length. Their alternate codons are one whole codon that is a stop codon, and one or two bases more. Ensembl VEP reads the alternate residues of such an allele as the stop alone, with no `X` (`TranscriptVariationAllele.pm`, line 774). Its `stop_retained_variant` test then holds, and `frameshift` returns 0 for it at line 1549, before line 1554.

Ensembl VEP writes `stop_retained_variant` for them and vep-rs writes the same. On the 7 rows above whose stop codon is kept, the alternate residues hold `X`, and that test is not passed. This entry makes no claim about those eight rows.

## How to reproduce it

`scripts/divergences/reproductions/frameshift_inside_a_held_stop_codon.toml` runs both engines under `--hgvs` on the five GRCh37 records of section 1. It compares the row of ENST00000435769 for each: the term, the impact, the coding and protein positions, the HGVSc and the HGVSp. The two records spelled before the stop codon are its controls, on which the engines agree.
