# A start codon left in place that UTRAnnotator reports gained

**The defect in one sentence.** Ensembl VEP's UTRAnnotator plugin reports an upstream start codon gained when an insertion written inside the codon leaves it where it stood, although that codon is the reference's own and the inserted bases make none. vep-rs reports no gained start for that codon: the row carries whichever of the other analyses fires, or `-` when none does.

**Why this is a defect.** Bases inserted after the `A` of an upstream `ATG` that open with `TG`, or after its `T` that open with `G`, leave an `ATG` at the codon's own place: the changed 5' UTR is the reference with the inserted bases directly after the untouched codon. The plugin's gained-start analysis is meant to report a start codon the variant creates. Here the start codon the plugin names stood at that place before the insertion, and inserted bases that add no start codon leave the UTR with as many as the reference holds.

A VCF writes an insertion that repeats neighbouring bases at its leftmost place on the genome. On a forward-strand transcript, bases put directly after the codon that end with `G`, or with `TG`, are therefore written one or two bases earlier, inside the codon.

**What vep-rs does instead.** vep-rs reports no gained start for that codon. The row carries whichever of the other analyses fires, under the fixed order described on [the precedence page](utrannotator-consequence-precedence.md), or `-` when none does. The frameshift analysis reads the insertion where the inserted bases stand, directly after the codon. An insertion whose length is no multiple of three therefore shifts the codon's own ORF, and the frameshift analysis lists that ORF unless a stop codon, or the UTR's end, follows the codon directly. The other analyses read the insertion where the record writes it.

When the inserted bases also hold a start codon of their own, vep-rs reports that codon gained, not the reference's. Inserted bases that end with the codon's bases before them only move the codon on, the subject of [the lost-start page](utrannotator-lost-start-after-insertion.md). Where they also open with its remaining bases, the changed sequence repeats the codon, and both engines name the copy at the codon's place gained. Inserted bases that do neither part the codon, and the lost-start analysis of both engines fires. The plugin's parameters, fields and analyses are described under [UTRAnnotator](../plugins/transcript-model.md#utrannotator).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `UTRAnnotator.pm`, lines 1285-1313: an insertion's position in the 5' UTR is the index of the base after it, in transcript order on both strands.
- Lines 339 and 345: the gained-start window is cut from the changed sequence, from two bases before that index. For an insertion it holds four bases.
- Lines 346-353: the first `ATG` in the window is reported gained. Nothing holds it against the reference sequence, so a codon that begins one or two bases before the insertion and still reads `ATG` is reported.
- Line 921: the frameshift analysis counts an ORF only when the variant's index lies at least three bases past the ORF's start codon, so an insertion written inside the codon is not counted for the codon's own ORF.

**What correcting it would change.** The plugin would hold the window's start codon against the reference sequence, report no gain for a codon that stood there already, and read the inserted bases after the codon they follow.

**The class it produces.** `utrannotator_gained_start_remade_in_place`: the `5UTR_consequence` and `5UTR_annotation` of an insertion's row on which Ensembl VEP reports a start codon gained at the place of a reference start codon the insertion is written inside, where vep-rs reports no gained start at that place. It is counted once per row. The concordance figures score consequence sets and never see the fields, so nothing is set aside in them.

**How often the shape occurs.** The index, [Intended divergences from Ensembl VEP](../intended-divergences.md), places this entry at its tier 2: ordinary data carries the shape, and the difference needs a plugin. The shape a record alone shows is an insertion. Of the 15,646,481 records of the five chromosome 21 datasets the index counts shapes on, 996,643 are insertions (63,698 per million), and 132 of the golden corpora hold one. The 5' UTR sequence decides whether it stands inside a start codon and leaves it in place. Four committed flag sets run the plugin, two in each UTRAnnotator corpus. On the genome-wide ClinVar files Ensembl VEP 116.3 reports such a codon gained on 59 rows (GRCh37) and 147 rows (GRCh38).

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. One base inserted after the second base of an upstream start codon

**In short.** The constructed record `21 10524584 T TG` inserts a `G` between cDNA positions 147 and 148 of `ENST00000427445`, after the `T` of the upstream ORF's `ATG` at cDNA position 146. The codon reads `ATG` as before, and its own `G` follows it: the changed sequence is the reference with a `G` directly after the codon. Ensembl VEP 116.3 reports that start codon gained; vep-rs reports the frameshift of the codon's own ORF.

### 1.1 What a reader needs to know first

UTRAnnotator runs five analyses over a 5' UTR variant: a start codon gained (`uAUG_gained`) or lost (`uAUG_lost`), a uORF's stop gained (`uSTOP_gained`) or lost (`uSTOP_lost`), and a frameshift in a uORF (`uFrameshift`). The row carries one flag in `5UTR_consequence` and that flag's effects in `5UTR_annotation`. A gained start codon is reported with `CapDistanceToStart`, the number of 5' UTR bases before the codon, which is one less than the cDNA position of its `A`. A frameshift lists each ORF it shifts, with the ORF's length before and after.

### 1.2 The record

A constructed GRCh38 insertion in the 5' UTR of `ENST00000427445` (TPTE, forward strand) on chromosome 21, against the Ensembl 116 GRCh38 cache and the published GRCh38 uORF table:

```
21	10524584	utr_ins_kept_start	T	TG	.	.	.
```

Both engines run `--plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt` with `--fasta`. The record is constructed because no ClinVar record of the shape lies on chromosome 21, the one chromosome the committed reference sequences cover.

### 1.3 What each engine writes

The row on `ENST00000427445`, one of nine transcripts whose rows differ the same way. Ensembl VEP:

```
utr_ins_kept_start	21:10524584-10524585	G	ENSG00000274391	ENST00000427445	Transcript	5_prime_UTR_variant	147-148	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;5UTR_annotation=KozakContext=AGAATGG:CapDistanceToStart=145:type=uORF:KozakStrength=Strong:DistanceToStop=42:DistanceToCDS=208;5UTR_consequence=5_prime_UTR_premature_start_codon_gain_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=6
```

vep-rs:

```
utr_ins_kept_start	21:10524584-10524585	G	ENSG00000274391	ENST00000427445	Transcript	5_prime_UTR_variant	147-148	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;5UTR_annotation=Evidence=False:KozakContext=AGAATGT:KozakStrength=Moderate:alt_type=uORF:alt_type_length=45:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51;5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=6
```

The `Existing_` counts agree.

### 1.4 Why Ensembl VEP writes that

The transcript's 5' UTR reads `CGAGAGAATGTTGGA` around the site, with the `A` of the start codon 145 bases from the cap. The record places a `G` between that codon's `T` and its `G`, giving `CGAGAGAATGGTTGGA`.

`UTRAnnotator.pm`, lines 1285-1313, give the insertion the index of the base after it, 147, the index of the codon's third base. Line 345 cuts the window from two bases before it, `ATGG`, and line 346 finds `ATG` at its first base, index 145. Lines 347-353 report a start codon gained there, with the Kozak context of the changed sequence, `AGAATGG`. The codon at 145 was there before the insertion; the window does not ask.

### 1.5 Why vep-rs's row is the correct one

The changed sequence decides. It is the reference with one `G` after the codon at 145: the codon's three bases are the reference's `A` and `T` and the inserted `G`, and the reference's own `G` follows. No start codon was made, and the changed 5' UTR holds six, as the reference does (`Existing_uORFs=6`).

What the insertion changes is the reading frame after that codon. One base shifts the codon's ORF, whose length falls from 51 bases to 45, and vep-rs reports that frameshift. Both engines report the same frameshift for the record `21 10524585 G GG`, which writes the same changed sequence with the base after the codon.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh38 --plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt
```

Ensembl VEP writes `5UTR_consequence=5_prime_UTR_premature_start_codon_gain_variant` with `CapDistanceToStart=145` on `ENST00000427445`; vep-rs writes `5UTR_consequence=5_prime_UTR_uORF_frameshift_variant`.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) applies `utrannotator_gained_start_remade_in_place` when a comparison case names it. It sets a row's `5UTR_consequence` and `5UTR_annotation` aside, counted once per row, when three things hold. The row's two cDNA positions are an insertion's two flanking bases. Ensembl VEP reports a gained start whose one `CapDistanceToStart` lies one below the lower of them on a row whose allele opens with `TG`, or two below it on a row whose allele opens with `G`. vep-rs's row carries the plugin's counts and reports no gained start at that index.

The class reads the row. The lower cDNA position is the 5' UTR index of the base after the insertion, so a codon reported one or two below it holds the insertion after its first or second base, and the allele's first bases say whether it still reads `ATG`. The allele is read on the transcript's strand, reverse-complemented where the row's `STRAND` is `-1`; a row that states no strand is refused. Where a row's record states its alleles, a reference allele that is not empty refuses the row.

Ensembl VEP writes one term on a row. Where a frameshift fires as well, its hash order may keep the frameshift, which lists the ORFs it counts from inside the codon and so leaves the codon's own ORF out. The class takes such a row when the comparison is given Ensembl VEP's outputs of the same input under other Perl hash seeds and one of them reports the gained start at that index.

UTRAnnotator's classes are tried in a fixed order and a differing field is counted under the first that accepts it; `utrannotator_frameshift_at_stop_boundary` and `utrannotator_lost_stop_written_before_stop_codon` follow this class. When vep-rs reports a frameshift that the seeded outputs also write on the row, the `5UTR_consequence` is counted under `utrannotator_consequence_precedence`. The `5UTR_annotation` is counted there too where vep-rs lists the ORFs those outputs list, and under this class where it lists the codon's own ORF beside them. When vep-rs reports a start codon gained inside the inserted bases and Ensembl VEP another term, the row is counted under `utrannotator_gained_start_window`, provided the comparison reads the record's alleles.

The `GRCh38-plugin-utr-annotator` and `GRCh37-plugin-utr-annotator` golden corpora carry no row of the shape, so the class applies to none of their rows.
