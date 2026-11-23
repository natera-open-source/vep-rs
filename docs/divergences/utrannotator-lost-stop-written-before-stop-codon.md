# A lost stop codon that UTRAnnotator reports as a frameshift

**The defect in one sentence.** Ensembl VEP's UTRAnnotator decides whether an insertion removes an upstream ORF's stop codon by the index the record gives. It looks only when the record writes the insertion inside the stop codon, and then reads only the three bases at the stop codon's index, so one changed sequence gets two answers. vep-rs decides from the changed sequence, so one changed sequence gets one row.

**Why this is a defect.** An upstream ORF can end inside a start codon. Where an `A` follows an upstream `ATG`, the codon's `T` and `G` open a `TGA`, and an ORF in that frame ends there. Two ways of writing one insertion matter here: directly before the codon, when the inserted bases open with `AT`, or after the codon's `T`, one base inside that stop codon. [The lost-start page](utrannotator-lost-start-after-insertion.md) describes them.

Both ways give one changed sequence, in which the start codon stands whole, moved on by the inserted bases. When the insertion's length is no multiple of three and the three bases at the stop codon's own index spell no stop codon, the ORF's stop codon is lost. The ORF runs on to its next stop codon, or to none.

The plugin describes its lost-stop analysis as one that annotates a variant that "removes a stop codon of an existing uORF". Its analyses flag the lost stop for the second way of writing and, for the first, a frameshift of the same ORF, although the changed sequence is the same. Its lost-stop analysis reads the changed sequence at a stop codon's index only for an insertion written inside the stop codon, so the first way of writing never reaches that read.

A VCF writes an insertion that repeats neighbouring bases at its leftmost place on the genome, which on a forward-strand transcript is its first place on the transcript. There a record writes such an insertion before the codon, never inside the stop codon.

An insertion written at the stop codon's first base is skipped as well, so its lost stop goes unreported ([section 2](#2-one-base-inserted-at-the-first-base-of-the-stop-codon)). And three inserted bases, or any multiple of three, can leave the stop codon whole directly after them. The ORF then still ends on it, in frame, yet the plugin reads the three bases at the old index and reports the stop lost ([section 3](#3-three-bases-inserted-inside-a-stop-codon-that-stands-whole-after-them)).

**What vep-rs does instead.** vep-rs reads an insertion from all its ways of writing only where one of them lies at a start codon's three bases or directly after the codon; any other insertion it reads at the record's own index, as Ensembl VEP does. Beside a start codon, then, one changed sequence gets one row. An upstream ORF's stop codon is lost when three things hold. One way of writing lies inside the stop codon. The changed sequence reads no stop codon at the stop codon's own index. Every way of writing lies inside the stop codon, or the insertion's length is no multiple of three.

Where one way of writing lies at or before the stop codon's first base, the stop codon stands whole after the inserted bases, and at a multiple of three the ORF still ends on it. The row's group names the ORF's next stop codon, or says that the changed 5' UTR holds none.

A lost stop stands before a lost start and a frameshift in the fixed order of [the precedence page](utrannotator-consequence-precedence.md). The plugin's parameters, fields and analyses are described under [UTRAnnotator](../plugins/transcript-model.md#utrannotator).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `UTRAnnotator.pm`, lines 1285-1313: an insertion's position in the 5' UTR is the index of the base after it, in transcript order on both strands.
- Line 578: the lost-stop analysis is described as one that annotates a variant that "removes a stop codon of an existing uORF".
- Line 628: the lost-stop analysis skips a variant whose index lies more than two bases past the first base of the ORF's stop codon.
- Line 631: it skips an insertion whose index lies at or before that first base. The test at line 633, written for insertions, never holds: it asks for an undefined reference allele, and an insertion's is defined and empty (lines 1160-1161).
- Lines 639-642: for an insertion one or two bases past the stop codon's first base, the three bases at the stop codon's index are read from the changed sequence, and the stop codon is reported lost unless they spell one. An insertion written directly before the start codon has its index one before the stop codon's first base, so the lines are never reached for it. They read no base past those three, so a stop codon that stands whole after the inserted bases is not seen.
- Lines 911-917: the frameshift analysis puts an upstream ORF's check point one index before the first base of its stop codon. Line 921 counts the ORF when the variant's index lies at or before that check point, which the index of such an insertion does.
- Lines 262-267: the row carries one term, the flagged one whose key the hash lists last.

**What correcting it would change.** The plugin would report the lost stop codon whichever way a record writes the inserted bases. An insertion written directly before a start codon, or at the stop codon's first base, would be read as the same bases written inside the stop codon are. Inserted codons that leave the stop codon whole and in frame would report no lost stop.

**The class it produces.** `utrannotator_lost_stop_written_before_stop_codon`: the `5UTR_consequence` and `5UTR_annotation` of an insertion's row on which one engine alone writes the lost stop of an upstream ORF. vep-rs alone writes it where Ensembl VEP writes the frameshift of that ORF (section 1), or no term, the lost start, the gained start or a frameshift of other ORFs (section 2). Ensembl VEP alone writes it for inserted codons that leave the stop codon whole (section 3). The two fields are set aside together and counted once per row. The concordance figures score consequence sets and never see the fields, so nothing is set aside in them.

**How often the shape occurs.** The index, [Intended divergences from Ensembl VEP](../intended-divergences.md), places this entry at its tier 2: ordinary data carries the shape, and the difference needs a plugin. The shape a record alone shows is an insertion. Of the 15,646,481 records of the five chromosome 21 datasets the index counts shapes on, 996,643 are insertions (63,698 per million), and 132 of the golden corpora hold one. The 5' UTR sequence decides whether it is written at a start codon inside which an upstream ORF ends. Four committed flag sets run the plugin, two in each UTRAnnotator corpus.

gnomAD v4.1 carries the shape at `7:80633290-80633291` on `ENST00000956936`, an insertion of two bases. Ensembl VEP lists three ORFs as shifted on that row, the ones that begin 295, 274 and 263 bases before the coding sequence. vep-rs reports the lost stop of the one that begins 263 bases before it, whose next stop codon begins 220 bases before the coding sequence.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Two bases inserted directly before a start codon inside which an upstream ORF ends

**In short.** The constructed record `21 10535884 C CAT` inserts two bases between cDNA positions 360 and 361 of `ENST00001091162`, directly before an upstream `ATG` at cDNA position 361. The inserted bases repeat the codon's `A` and `T`, so the changed sequence is the reference with `AT` written after the codon's `T`, inside the `TGA` that closes an upstream ORF. In the changed sequence that ORF reads `TAT` where its stop codon stood and ends nine bases later. Ensembl VEP 116.3 reports the ORF's frameshift; vep-rs reports its lost stop codon.

### 1.1 What a reader needs to know first

UTRAnnotator runs five analyses over a 5' UTR variant: a start codon gained (`uAUG_gained`) or lost (`uAUG_lost`), a uORF's stop gained (`uSTOP_gained`) or lost (`uSTOP_lost`), and a frameshift in a uORF (`uFrameshift`). The row carries one flag in `5UTR_consequence` and that flag's effects in `5UTR_annotation`, one group of `key=value` pairs for each ORF the flag names.

A lost stop's group says where the ORF ends in the changed sequence. `AltStop` is `True` when the ORF meets another stop codon inside the changed 5' UTR, and `AltStopDistanceToCDS` is the number of 5' UTR bases from that stop codon's first base to the coding sequence.

A frameshift's group gives the ORF's place and lengths. `ref_StartDistanceToCDS` is the number of 5' UTR bases from the `A` of the ORF's start codon to the coding sequence. `ref_type_length` is the ORF's length with its stop codon, and `alt_type_length` its length in the changed sequence. Both kinds of group carry the `KozakContext` and `KozakStrength` of the ORF's start codon.

### 1.2 The record

A constructed GRCh38 insertion in the 5' UTR of `ENST00001091162` (TPTE, forward strand) on chromosome 21, against the Ensembl 116 GRCh38 cache and the published GRCh38 uORF table:

```
21	10535884	utr_ins_before_stop_codon	C	CAT	.	.	.
```

Both engines run `--plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt` with `--fasta` and `--hgvs`. The `HGVSc` that `--hgvs` prints names the inserted bases as a duplication of the start codon's first two bases, which is the second way of writing them. The record is constructed because the committed reference sequences cover chromosome 21 alone, and the gnomAD record named above lies on chromosome 7.

### 1.3 What each engine writes

The row on `ENST00001091162`, the one row of the record's nine on which the engines differ. Ensembl VEP, under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`:

```
utr_ins_before_stop_codon	21:10535884-10535885	AT	ENSG00000274391	ENST00001091162	Transcript	5_prime_UTR_variant,NMD_transcript_variant	360-361	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00001091162.1:c.-101_-100dup;HGVS_OFFSET=2;5UTR_annotation=alt_type_length=45:KozakContext=TACATGG:KozakStrength=Moderate:ref_type_length=36:ref_type=uORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=133;5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;Existing_InFrame_oORFs=1;Existing_OutOfFrame_oORFs=0;Existing_uORFs=9
```

vep-rs:

```
utr_ins_before_stop_codon	21:10535884-10535885	AT	ENSG00000274391	ENST00001091162	Transcript	5_prime_UTR_variant,NMD_transcript_variant	360-361	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00001091162.1:c.-101_-100dup;HGVS_OFFSET=2;5UTR_annotation=AltStop=True:AltStopDistanceToCDS=93:Evidence=False:FrameWithCDS=inFrame:KozakContext=TACATGG:KozakStrength=Moderate;5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant;Existing_InFrame_oORFs=1;Existing_OutOfFrame_oORFs=0;Existing_uORFs=9
```

The `Existing_` counts, the `HGVSc` and every other column agree. The two groups carry the same `KozakContext` and `KozakStrength`: both describe the ORF that begins 133 bases before the coding sequence. The order of the pairs on Ensembl VEP's row is its hash order, which the comparison tools take as a set.

### 1.4 Why Ensembl VEP writes that

The transcript's 5' UTR is `GACAGCCATGATGGG` around the site. Its first `ATG` is the start codon the record is written before, with its `A` 360 bases from the cap. The `TGA` that begins at that codon's `T` closes an upstream ORF whose start codon lies 32 bases before the `A`, 133 bases before the coding sequence. The `A` that ends the `TGA` opens the excerpt's second `ATG`, which begins after the insertion and plays no part here.

The record's `AT` lands directly before the first codon, giving `GACAGCCATATGATGGG`.

`UTRAnnotator.pm`, lines 1285-1313, give the insertion the index of the base after it, 360: the index of the codon's `A`, one before the stop codon's first index, 361. The lost-stop analysis passes line 628, since 360 does not lie more than two past 361. Line 631 then skips the ORF, since the insertion's index lies before the stop codon's first base. Lines 639-642, which read the three bases at index 361 of the changed sequence, are never reached, and no stop codon is reported lost.

The frameshift analysis puts the ORF's check point at index 360, one before its stop codon (lines 911-917). Line 921 counts an ORF when the variant's index lies at or before its check point and at least three past the first base of its start codon. Index 360 lies at the check point, and two inserted bases shift the frame, so the ORF is listed with its two lengths, 36 and 45 bases.

Ensembl VEP's lost-start analysis fires on the record as well, since the three bases at index 360 of the changed sequence do not spell `ATG`; that is the subject of [the lost-start page](utrannotator-lost-start-after-insertion.md). The plugin writes one term on a row, the one whose key its hash lists last, as [the precedence page](utrannotator-consequence-precedence.md) describes. The row of 1.3 is the one written under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`, where that term is the frameshift.

### 1.5 Why vep-rs's row is the correct one

The changed sequence decides. It is the reference with `AT` before the start codon at index 360, and it is the same sequence whichever way the record writes those bases: `21 10535886 T TAT` writes them after the codon's `T`, one base inside the stop codon.

Translated from its own start codon, the ORF reads every base before index 360 unchanged. The three bases at index 361, where its stop codon stood, are `TAT`, which is no stop codon. The ORF reads on to the stop codon that begins at index 370 of the changed sequence, 93 bases before the coding sequence, and its length rises from 36 bases to 45. Its stop codon is lost and the changed 5' UTR gives it another, which vep-rs writes as `AltStop=True` with `AltStopDistanceToCDS=93`.

Ensembl VEP's own row for the other way of writing says the same. For `21 10535886 T TAT` the insertion's index is 362, one past the stop codon's first base. The lost-stop analysis reads index 361, finds no stop codon and reports the lost stop with the group vep-rs writes, `AltStop=True` and `AltStopDistanceToCDS=93`. Its lost-start analysis fires on that record as well. The row was written under `PERL_HASH_SEED=5` with `PERL_PERTURB_KEYS=0`, where the lost stop is the term the plugin keeps. vep-rs writes the lost stop with that group for both records.

Both terms are true of the ORF: its stop codon is lost, and two inserted bases shift it. A row carries one term, and vep-rs's fixed order puts a lost stop codon before a frameshift. Ensembl VEP's two rows give one changed sequence two terms, by the way the record is written.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh38 --hgvs \
  --plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt
```

On `ENST00001091162` Ensembl VEP writes `5UTR_consequence=5_prime_UTR_uORF_frameshift_variant` under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`, and its one group carries `ref_StartDistanceToCDS=133`. vep-rs writes `5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant` with `AltStop=True` and `AltStopDistanceToCDS=93`. For `21 10535886 T TAT` both engines write the lost stop, Ensembl VEP under `PERL_HASH_SEED=5` with `PERL_PERTURB_KEYS=0`.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) applies `utrannotator_lost_stop_written_before_stop_codon` when a comparison case names it. It sets a row's `5UTR_consequence` and `5UTR_annotation` aside, counted once per row, when four things hold. Ensembl VEP's term is the frameshift and vep-rs's the lost stop. The row's two cDNA positions are an insertion's two flanking bases and, where the row's record states its alleles, its reference allele is empty once the first base its alleles share is dropped. Each group of vep-rs's annotation matches exactly one upstream ORF (`ref_type=uORF`) of Ensembl VEP's list, a different one for each group. vep-rs's row carries UTRAnnotator's three `Existing_` counts, which the plugin writes on every row it annotates; their values are not compared by this class.

The class reads nothing but the two engines' rows and, where the row's record states them, the record's alleles. A group matches an ORF when it carries that ORF's `KozakContext` and `KozakStrength` and the alternative stop codon the ORF's printed lengths give. Where the ORF's `alt_type` is `uORF`, the group reads `AltStop=True`, and its `AltStopDistanceToCDS` equals the ORF's `ref_StartDistanceToCDS` plus the allele's length, less its `alt_type_length`, plus three. Otherwise the group reads `AltStop=False` with `AltStopDistanceToCDS=NA`. On this row the listed ORF has `ref_StartDistanceToCDS=133` and `alt_type_length=45`, and the allele is two bases long: 133 + 2 - 45 + 3 = 93, the value of vep-rs's group.

The class leaves two things open. A row whose lost-stop group matches no ORF of Ensembl VEP's list, or more than one, stays a difference. So does a row whose alternative stop differs from the one the listed ORF's lengths give.

UTRAnnotator's classes are tried in a fixed order and a differing field is counted under the first that accepts it; this class comes last, after `utrannotator_frameshift_at_stop_boundary`. Two of the others can take a row of this shape. A row whose two terms Ensembl VEP's outputs under other hash seeds both write is counted under `utrannotator_consequence_precedence`, when the comparison is given those outputs. A row on which Ensembl VEP's hash order writes the lost start is counted under `utrannotator_lost_start_after_insertion`.

The `GRCh38-plugin-utr-annotator` and `GRCh37-plugin-utr-annotator` golden corpora carry no row of the shape, so the class applies to none of their rows.

## 2. One base inserted at the first base of the stop codon

**In short.** The constructed record `21 10535885 A AT` inserts one `T` between cDNA positions 361 and 362 of `ENST00001091162`: after the `A` of the upstream `ATG` of section 1, at the first base of the `TGA` that closes the upstream ORF. The same sequence is the reference with a `T` written one base further on, inside the stop codon. The ORF reads `TTG` where its stop codon stood and ends 84 bases later. Ensembl VEP 116.3 reports the start codon lost and no lost stop; vep-rs reports the lost stop codon.

### 2.1 What a reader needs to know first

The analyses, the two fields and a lost stop's group are described in [1.1](#11-what-a-reader-needs-to-know-first). A lost start's group names its codon by `CapDistanceToStart`, the number of 5' UTR bases before the codon's `A`.

### 2.2 The record

The transcript and the run are those of [1.2](#12-the-record):

```
21	10535885	utr_ins_at_stop_first_base	A	AT	.	.	.
```

The `HGVSc` names the inserted base as a duplication of the codon's `T`, and `HGVS_OFFSET=1` says that place lies one base past the record's own.

### 2.3 What each engine writes

The row on `ENST00001091162`, the one row of the record's nine on which the engines differ. Ensembl VEP, under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`:

```
utr_ins_at_stop_first_base	21:10535885-10535886	T	ENSG00000274391	ENST00001091162	Transcript	5_prime_UTR_variant,NMD_transcript_variant	361-362	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00001091162.1:c.-100dup;HGVS_OFFSET=1;5UTR_annotation=KozakContext=GCCATGA:type=uORF:KozakStrength=Moderate:DistanceToStop=84:DistanceToCDS=101:CapDistanceToStart=360:Evidence=False;5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;Existing_InFrame_oORFs=1;Existing_OutOfFrame_oORFs=0;Existing_uORFs=9
```

vep-rs:

```
utr_ins_at_stop_first_base	21:10535885-10535886	T	ENSG00000274391	ENST00001091162	Transcript	5_prime_UTR_variant,NMD_transcript_variant	361-362	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00001091162.1:c.-100dup;HGVS_OFFSET=1;5UTR_annotation=AltStop=True:AltStopDistanceToCDS=17:Evidence=False:FrameWithCDS=outOfFrame:KozakContext=TACATGG:KozakStrength=Moderate;5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant;Existing_InFrame_oORFs=1;Existing_OutOfFrame_oORFs=0;Existing_uORFs=9
```

Every other column agrees. Ensembl VEP's group names the start codon 360 bases from the cap, vep-rs's the ORF of section 1, 133 bases before the coding sequence.

### 2.4 Why Ensembl VEP writes that

Around the site the 5' UTR is `GACAGCCATGATGGG`, as [1.4](#14-why-ensembl-vep-writes-that) describes it: the first `ATG` begins at index 360 and the ORF's `TGA` at index 361. The record's `T` lands after the codon's `A`, giving `GACAGCCATTGATGGG`.

Lines 1285-1313 give the insertion the index of the base after it, 361, the stop codon's first base. Line 631 skips an insertion whose index does not lie past that base, so lines 639-642 never read the changed sequence and no stop codon is reported lost. The frameshift analysis does not list the ORF either: its check point is index 360 (lines 911-917), and index 361 lies past it (line 921).

The lost-start analysis alone fires: the insertion's index lies inside the start codon (lines 764-769), and the three bases at index 360 of the changed sequence are `ATT` (lines 776-778). The row carries the lost start under any hash order.

### 2.5 Why vep-rs's row is the correct one

The changed sequence decides, and it is the same sequence when the record writes the base one place further on: `21 10535886 T TT` writes it after the codon's `T`, one base inside the stop codon.

Translated from its own start codon, the ORF reads every base before index 361 unchanged. The three bases at index 361 are `TTG`, no stop codon, and one inserted base shifts every codon after it. The ORF reads on to the stop codon at index 445 of the changed 5' UTR, which is 462 bases long: 17 bases before the coding sequence. From its first base at index 328 the ORF grows from 36 bases to 120.

Ensembl VEP's own row for the other way of writing says the same. For `21 10535886 T TT` the insertion's index is 362, one past the stop codon's first base. Lines 639-642 read index 361, find no stop codon and report the lost stop with the group vep-rs writes, under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`.

The start codon is lost as well, since both ways of writing part it. A row carries one term, and vep-rs's fixed order puts a lost stop before a lost start.

### 2.6 Try it yourself

With the command of [1.6](#16-try-it-yourself), Ensembl VEP writes `5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant` on `ENST00001091162` and vep-rs `5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant` with `AltStopDistanceToCDS=17`. For `21 10535886 T TT` both write that lost stop, Ensembl VEP under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`.

### 2.7 How the comparator treats this class

The class sets both fields aside under a second rule, counted once per row, where the first has refused the row. The row is an insertion's as [1.7](#17-how-the-comparator-treats-this-class) tells one, vep-rs's row carries the plugin's counts, the row's `Allele` is a sequence, and each of these holds:

| The rule reads | and finds |
| --- | --- |
| vep-rs's term | the lost stop, beside an annotation whose every group carries a `KozakContext` and a `KozakStrength` |
| Ensembl VEP's two fields | both empty, or the lost start, the gained start or the frameshift; a frameshift's list holds no group with the `KozakContext` and `KozakStrength` of a group of vep-rs's |
| the allele's length | no multiple of three |
| `HGVS_OFFSET` | a whole number other than zero: the insertion can be written further on, inside the stop codon |
| Ensembl VEP's outputs under other hash seeds, where the comparison is given them | none writes the lost stop on the row |

On this row vep-rs's one group carries `KozakContext=TACATGG` and `KozakStrength=Moderate`, and Ensembl VEP's term is the lost start. The allele is one base long, and the row prints `HGVS_OFFSET=1`.

## 3. Three bases inserted inside a stop codon that stands whole after them

**In short.** The constructed record `21 10535887 G GCTG` inserts three bases between cDNA positions 363 and 364 of `ENST00001091162`, after the `G` of the `TGA` that closes the upstream ORF of section 1. The same sequence is the reference with `TGC` written at the stop codon's first base, before the whole `TGA`. The ORF reads `TGC` where its stop codon stood and ends on that `TGA`, three bases on and in frame. Ensembl VEP 116.3 reports the stop codon lost; vep-rs writes no term.

### 3.1 What a reader needs to know first

A lost stop's group is described in [1.1](#11-what-a-reader-needs-to-know-first).

### 3.2 The record

The transcript and the run are those of [1.2](#12-the-record):

```
21	10535887	utr_codon_inside_stop_codon	G	GCTG	.	.	.
```

The `HGVSc` is an `ins` with no `HGVS_OFFSET`: the record writes the bases at the last of the three places that give this sequence.

### 3.3 What each engine writes

The row on `ENST00001091162`, the one row of the record's nine on which the engines differ. Ensembl VEP, under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`:

```
utr_codon_inside_stop_codon	21:10535887-10535888	CTG	ENSG00000274391	ENST00001091162	Transcript	5_prime_UTR_variant,NMD_transcript_variant	363-364	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00001091162.1:c.-99_-98insCTG;5UTR_annotation=KozakContext=TACATGG:KozakStrength=Moderate:FrameWithCDS=outOfFrame:AltStopDistanceToCDS=100:AltStop=True:Evidence=False;5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant;Existing_InFrame_oORFs=1;Existing_OutOfFrame_oORFs=0;Existing_uORFs=9
```

vep-rs:

```
utr_codon_inside_stop_codon	21:10535887-10535888	CTG	ENSG00000274391	ENST00001091162	Transcript	5_prime_UTR_variant,NMD_transcript_variant	363-364	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00001091162.1:c.-99_-98insCTG;5UTR_annotation=-;5UTR_consequence=-;Existing_InFrame_oORFs=1;Existing_OutOfFrame_oORFs=0;Existing_uORFs=9
```

Every other column agrees. Ensembl VEP's group names the ORF of section 1 by its `KozakContext`, `TACATGG`.

### 3.4 Why Ensembl VEP writes that

The record's `CTG` lands after the `G` of the first `ATG` of `GACAGCCATGATGGG`, giving `GACAGCCATGCTGATGGG`. Lines 1285-1313 give the insertion the index 363, two past the stop codon's first base at index 361, so it passes lines 628 and 631. Lines 639-642 read the three bases at index 361 of the changed sequence, `TGC`. They spell no stop codon, and the stop codon is reported lost. No base past those three is read.

The lost-start analysis fires as well, for the excerpt's second `ATG`, which the insertion only moves on; that is the subject of [the lost-start page](utrannotator-lost-start-after-insertion.md). The row of 3.3 is the one written under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`, where the plugin keeps the lost stop.

### 3.5 Why vep-rs's row is the correct one

The changed sequence decides, and it is the same sequence when the record writes the bases two places earlier: `21 10535885 A ATGC` writes `TGC` after the start codon's `A`, at the stop codon's first base. Written that way, the `TGA` follows the inserted bases untouched.

Translated from its own start codon, the ORF reads every base before index 361 unchanged, then the codon `TGC`, then `TGA` at index 364. Three inserted bases keep the frame, so that `TGA` is the ORF's own stop codon, moved three bases on. The ORF gains one codon, from 36 bases to 39, and loses no stop codon.

Ensembl VEP's group names that very codon as the alternative stop: `AltStopDistanceToCDS=100` is index 364 of the changed 5' UTR, 464 bases long. For `21 10535885 A ATGC` line 631 skips the ORF, and Ensembl VEP reports no lost stop. Under any hash order its row carries a gained start there, the start codon that stands in place, where vep-rs writes no term; the entry `utrannotator_gained_start_remade_in_place` reads that pair, on [the page on a start codon left in place](utrannotator-gained-start-remade-in-place.md).

Both start codons stand whole, the first at its own index, and three bases shift no ORF. vep-rs writes no term for either record.

### 3.6 Try it yourself

With the command of [1.6](#16-try-it-yourself), Ensembl VEP writes `5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant` with `AltStopDistanceToCDS=100` on `ENST00001091162`, under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`. vep-rs writes `5UTR_consequence=-`.

### 3.7 How the comparator treats this class

The class sets both fields aside under a third rule, counted once per row, on a row of the kind [2.7](#27-how-the-comparator-treats-this-class) asks for: an insertion's, with the plugin's counts on vep-rs's row and a sequence as its `Allele`. Ensembl VEP's term is the lost stop, and vep-rs's is another of the five terms, or vep-rs fills neither field. The allele's length is a multiple of three. Each group of Ensembl VEP's annotation reads `AltStop=True`, and its `AltStopDistanceToCDS` exceeds the point by a number below the allele's length that is no multiple of three:

| It exceeds the point by | The stop codon the ORF meets next |
| --- | --- |
| one or two | its own, moved on by the allele's length |
| a greater such number | one the inserted bases themselves put in the ORF's frame, before its own |

The point is the number of bases before the coding sequence of the 5' UTR base directly after the insertion, as the record writes it. A row whose `HGVSc` states no insertion point is refused. On this row the `HGVSc` is `c.-99_-98insCTG` with no `HGVS_OFFSET`, so the point is 98, and the group's `AltStopDistanceToCDS=100` exceeds it by two.

`utrannotator_consequence_precedence` and `utrannotator_gained_start_window` are tried before this class and take such a row where they apply. The two golden corpora carry no row of the shapes of sections 2 and 3.
