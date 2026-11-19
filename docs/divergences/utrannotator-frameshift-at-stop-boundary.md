# A shifted upstream ORF that UTRAnnotator leaves off its frameshift list

**The defect in one sentence.** Ensembl VEP's UTRAnnotator omits a shifted upstream ORF from its frameshift list for an insertion written inside a start codon. vep-rs lists that ORF whichever way the record writes the insertion, so one changed sequence gets one row.

**Why this is a defect.** An upstream ORF can end inside a start codon. Where an `A` follows an upstream `ATG`, the codon's `T` and `G` open a `TGA`, and an ORF in that frame ends there. Two ways of writing an insertion that only moves the start codon on matter here: directly before the codon, or after its `A` when the inserted bases end with `A`. [The lost-start page](utrannotator-lost-start-after-insertion.md) describes them.

Both ways give one changed sequence, in which the inserted bases stand inside that ORF, before its stop codon. When their length is no multiple of three and the three bases at the stop codon's own index spell no stop codon, the ORF runs on to a later one: it is shifted, and its length changes.

The plugin's frameshift analysis lists the ORF for the first way of writing and not for the second, although the changed sequence, and with it the ORF's length before and after, is the same. No other analysis reports the ORF for the second way. A comment in the source leaves an insertion at the boundary of a stop codon to the stop codon's own analysis, and the lost-stop analysis skips an insertion written at a stop codon's first base.

A VCF writes an insertion that repeats neighbouring bases at its leftmost place on the genome, which on a reverse-strand transcript is its last place on the transcript. There such an insertion is written the second way, after the codon's `A`.

**What vep-rs does instead.** vep-rs lists the ORF for either way of writing, so one changed sequence gets one row. Its frameshift analysis takes an insertion that only moves a start codon on as standing directly before that codon. For the ORF that ended inside the codon it looks at the changed sequence: where the three bases at the stop codon's own index spell no stop codon, the ORF is shifted, and the analysis lists it beside every other ORF the insertion shifts. The row carries the frameshift unless an analysis earlier in the fixed order of [the precedence page](utrannotator-consequence-precedence.md) fires. The plugin's parameters, fields and analyses are described under [UTRAnnotator](../plugins/transcript-model.md#utrannotator).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `UTRAnnotator.pm`, lines 1285-1313: an insertion's position in the 5' UTR is the index of the base after it, in transcript order on both strands.
- Lines 911-913: the frameshift analysis puts an upstream ORF's check point one index before the first base of its stop codon.
- Line 921: an ORF is counted only when the variant's index lies at or before that check point. The line holds the record's own index against it and looks at no base of the changed sequence. An insertion written after the first base of the start codon has the index of the stop codon's first base, one past the check point.
- Lines 919-920: the comment says cases at the boundary of a start or stop codon are ignored "since the effect on start/stop would be evaluated anyway".
- Line 631: the lost-stop analysis skips an insertion whose index does not lie past the stop codon's first base, so it does not evaluate this one.

**What correcting it would change.** The plugin would list the ORF whichever way a record writes the inserted bases. An insertion that only moves a start codon on would be counted for the ORF that ends inside that codon, as the same bases written directly before the codon are.

**The class it produces.** `utrannotator_frameshift_at_stop_boundary`: the `5UTR_annotation` of an insertion's row on which both engines write the frameshift and vep-rs lists, beside every ORF Ensembl VEP lists, an upstream ORF whose stop codon begins directly after the insertion. The row's `5UTR_consequence` is equal and is not set aside. The class is counted once per row. The concordance figures score consequence sets and never see the field, so nothing is set aside in them.

**How often the shape occurs.** The index, [Intended divergences from Ensembl VEP](../intended-divergences.md), places this entry at its tier 2: ordinary data carries the shape, and the difference needs a plugin. The shape a record alone shows is an insertion. Of the 15,646,481 records of the five chromosome 21 datasets the index counts shapes on, 996,643 are insertions (63,698 per million), and 93 of the golden corpora hold one. The 5' UTR sequence decides whether it is written inside a start codon at which an upstream ORF ends. Four committed flag sets run the plugin, two in each UTRAnnotator corpus.

gnomAD v4.1 carries the shape at `2:68048431-68048432` on `ENST00000871248`. Ensembl VEP lists one ORF on that row, the one that begins 1,155 bases before the coding sequence. vep-rs lists two, that one and the one that begins 1,151 bases before the coding sequence.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. One base inserted after the first base of a start codon inside which an upstream ORF ends

**In short.** The constructed record `21 14991217 A AT` inserts one base between cDNA positions 474 and 475 of `ENST00001008934`, a reverse-strand transcript, after the `A` of an upstream `ATG` at cDNA position 474. On the transcript's strand the inserted base is `A`, so the changed sequence is the reference with an `A` directly before that codon. One upstream ORF ends on the `TGA` that begins at the codon's `T`, and a second runs through the codon. The inserted base shifts both. Ensembl VEP 116.2 lists the second alone; vep-rs lists both.

### 1.1 What a reader needs to know first

UTRAnnotator runs five analyses over a 5' UTR variant: a start codon gained (`uAUG_gained`) or lost (`uAUG_lost`), a uORF's stop gained (`uSTOP_gained`) or lost (`uSTOP_lost`), and a frameshift in a uORF (`uFrameshift`). The row carries one flag in `5UTR_consequence` and that flag's effects in `5UTR_annotation`. A frameshift lists each ORF it shifts as one group of `key=value` pairs, the groups parted by commas.

In a group, `ref_StartDistanceToCDS` is the number of 5' UTR bases from the `A` of the ORF's start codon to the coding sequence. `ref_type_length` is the ORF's length with its stop codon, and `alt_type_length` its length in the changed sequence. A row's `Allele` is written on the genome's forward strand: on a reverse-strand transcript the transcript carries its complement.

### 1.2 The record

A constructed GRCh38 insertion in the 5' UTR of `ENST00001008934` (NRIP1, reverse strand) on chromosome 21, against the Ensembl 116 GRCh38 cache and the published GRCh38 uORF table:

```
21	14991217	utr_ins_at_stop_boundary	A	AT	.	.	.
```

Both engines run `--plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt` with `--fasta` and `--hgvs`, which prints the `HGVSc` a comparison takes the insertion point from. The record is constructed because the committed reference sequences cover chromosome 21 alone, and the gnomAD record named above lies on chromosome 2.

### 1.3 What each engine writes

The row on `ENST00001008934`; the record's row on `ENST00001008935` differs the same way. Ensembl VEP:

```
utr_ins_at_stop_boundary	21:14991217-14991218	T	ENSG00000180530	ENST00001008934	Transcript	5_prime_UTR_variant	474-475	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00001008934.1:c.-352dup;5UTR_annotation=alt_type_length=12:KozakContext=GGCATGT:KozakStrength=Moderate:ref_type_length=33:ref_type=uORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=359;5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=8
```

vep-rs:

```
utr_ins_at_stop_boundary	21:14991217-14991218	T	ENSG00000180530	ENST00001008934	Transcript	5_prime_UTR_variant	474-475	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00001008934.1:c.-352dup;5UTR_annotation=Evidence=False:KozakContext=TGCATGG:KozakStrength=Moderate:alt_type=uORF:alt_type_length=102:ref_StartDistanceToCDS=414:ref_type=uORF:ref_type_length=66,Evidence=False:KozakContext=GGCATGT:KozakStrength=Moderate:alt_type=uORF:alt_type_length=12:ref_StartDistanceToCDS=359:ref_type=uORF:ref_type_length=33;5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=8
```

The `5UTR_consequence`, the `Existing_` counts and every other column agree, and the group both rows carry holds the same pairs. The order of the pairs on Ensembl VEP's row is its hash order, which the comparison tools take as a set.

### 1.4 Why Ensembl VEP writes that

The transcript's 5' UTR is `ATGTTACATGATACG` around the site. Its second `ATG` is the start codon the record is written inside, with its `A` 473 bases from the cap. The `TGA` that begins at that codon's `T` closes an upstream ORF whose start codon lies 62 bases before the `A`, 414 bases before the coding sequence. The excerpt's first `ATG`, seven bases before the `A`, opens a second upstream ORF, which runs through the codon in another frame to a stop codon that begins 23 bases after the `A`.

The record's `T` is `A` on the transcript's strand and lands after the codon's `A`, giving `ATGTTACAATGATACG`.

`UTRAnnotator.pm`, lines 1285-1313, give the insertion the index of the base after it, 474: the index of the codon's `T` and of the first ORF's stop codon. Lines 911-913 put an ORF's check point one index before its stop codon, at 473 for the first ORF and 495 for the second. Line 921 counts an ORF only when the insertion's index lies at or before its check point. Index 474 lies past 473, so the first ORF is not counted; the second is.

The comment above that line says cases at the boundary of a start or stop codon are ignored "since the effect on start/stop would be evaluated anyway". For this insertion no analysis evaluates the stop codon: line 631 of the lost-stop analysis skips an insertion whose index does not lie past the stop codon's first base.

Ensembl VEP's lost-start analysis fires on the record as well, since the three bases at index 473 of the changed sequence do not spell `ATG`; that is the subject of [the lost-start page](utrannotator-lost-start-after-insertion.md). The plugin writes one term on a row, the one whose key its hash lists last, as [the precedence page](utrannotator-consequence-precedence.md) describes. The row of 1.3 is the one written under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`, where that term is the frameshift.

### 1.5 Why vep-rs's row is the correct one

The changed sequence decides. It is the reference with one `A` before the start codon at index 473, and it is the same sequence whichever way the record writes that base: `21 14991218 T TT` writes it directly before the codon.

Translated from its own start codon, the first ORF meets the inserted base before its stop codon. The three bases at index 474, where that stop codon stood, are `ATG`, the start codon moved one base on. The ORF runs through them and ends 36 bases later: its length rises from 66 bases to 102. The second ORF meets the moved `TGA` in its own frame, and its length falls from 33 bases to 12.

One inserted base shifts both ORFs, and vep-rs lists both. For the record written directly before the codon vep-rs writes the same two groups, and so does Ensembl VEP: there the insertion's index is 473, which lies at the first ORF's check point.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh38 --hgvs \
  --plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt
```

Both engines write `5UTR_consequence=5_prime_UTR_uORF_frameshift_variant` on `ENST00001008934`, Ensembl VEP under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`. Its `5UTR_annotation` holds one group, with `ref_StartDistanceToCDS=359`; vep-rs's holds that group and a second, with `ref_StartDistanceToCDS=414`.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies `utrannotator_frameshift_at_stop_boundary` when a comparison case names it. It sets a row's `5UTR_annotation` aside, counted once per row, when four things hold. Both rows write the frameshift term. Every group of Ensembl VEP's row is among vep-rs's, and vep-rs's row holds at least one more. Each group vep-rs adds is an ORF that ends in the 5' UTR (`ref_type=uORF`), whose two lengths differ and whose stop codon begins directly after the insertion point. No group both rows carry is such an ORF.

The class takes everything from the row. A group's `ref_StartDistanceToCDS` less its `ref_type_length`, plus three, is the distance of its stop codon's first base before the coding sequence.

The insertion point is read from the row's `HGVSc`: an `ins` between two neighbouring 5' UTR positions, or a `dup`, whose copy follows the duplicated bases. Ensembl VEP writes that notation at the last place on the transcript that gives the same sequence, and prints as `HGVS_OFFSET` the number of bases it moved the insertion to get there, with a minus sign on a reverse-strand transcript. The class moves the point that many bases back towards the cap, to the place where the record writes the insertion; a row without an `HGVS_OFFSET` is taken as its notation stands. A row whose `HGVSc` is neither notation is refused.

On this row the `HGVSc` is `c.-352dup` and no `HGVS_OFFSET` is printed, so the inserted base follows position -352. The group vep-rs adds has its stop codon's first base 414 - 66 + 3 = 351 bases before the coding sequence, directly after it.

The class leaves two things open. A row whose two lists differ in any other way stays a difference. So does a row on which either engine writes another term: whether the insertion also loses the ORF's stop codon is not decided by this class. A run of Ensembl VEP whose hash order writes the lost start on the row gives the two rows different terms, and this class refuses it.

UTRAnnotator's classes are tried in a fixed order and a differing field is counted under the first that accepts it; `utrannotator_lost_stop_written_before_stop_codon` alone follows this class.

One of the others, `utrannotator_gained_start_remade_in_place`, can take a row of this shape. It applies where the inserted bases follow the first or second base of a start codon and open with that codon's remaining bases, so that the codon still stands at its own index. Ensembl VEP then also flags a start codon gained there, and its outputs of the same input under other hash seeds, which a comparison can be given, write that term on the row. With those outputs the row's `5UTR_annotation` is counted under that class; without them it is counted under this one.

The `GRCh38-plugin-utr-annotator` and `GRCh37-plugin-utr-annotator` golden corpora carry no row of the shape, so the class applies to none of their rows.
