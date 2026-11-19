# A lost stop codon that UTRAnnotator reports as a frameshift

**The defect in one sentence.** Ensembl VEP's UTRAnnotator misses a lost stop codon when the insertion that removes it is written before the stop codon. For an insertion written that way it flags a frameshift of the same ORF. vep-rs reports the lost stop for either way of writing, so one changed sequence gets one row.

**Why this is a defect.** An upstream ORF can end inside a start codon. Where an `A` follows an upstream `ATG`, the codon's `T` and `G` open a `TGA`, and an ORF in that frame ends there. Two ways of writing one insertion matter here: directly before the codon, when the inserted bases open with `AT`, or after the codon's `T`, one base inside that stop codon. [The lost-start page](utrannotator-lost-start-after-insertion.md) describes them.

Both ways give one changed sequence, in which the start codon stands whole, moved on by the inserted bases. When the insertion's length is no multiple of three and the three bases at the stop codon's own index spell no stop codon, the ORF's stop codon is lost. The ORF runs on to its next stop codon, or to none.

The plugin describes its lost-stop analysis as one that annotates a variant that "removes a stop codon of an existing uORF". Its analyses flag the lost stop for the second way of writing and, for the first, a frameshift of the same ORF, although the changed sequence is the same. Its lost-stop analysis reads the changed sequence at a stop codon's index only for an insertion written inside the stop codon, so the first way of writing never reaches that read.

A VCF writes an insertion that repeats neighbouring bases at its leftmost place on the genome, which on a forward-strand transcript is its first place on the transcript. There a record writes such an insertion before the codon, never inside the stop codon.

**What vep-rs does instead.** vep-rs reports the lost stop for either way of writing, so one changed sequence gets one row. Its stop-codon analyses take an insertion that only moves a start codon on as standing directly before that codon. For the ORF that ended inside the codon vep-rs looks at the changed sequence. The stop codon is lost when three things hold: the insertion's length is no multiple of three, the same sequence can be written with the insertion inside the stop codon, and the three bases at the stop codon's own index spell no stop codon. The row's group names the ORF's next stop codon, or says that the changed 5' UTR holds none.

The ORF is shifted as well, and the lost stop stands before the frameshift in the fixed order of [the precedence page](utrannotator-consequence-precedence.md). The plugin's parameters, fields and analyses are described under [UTRAnnotator](../plugins/transcript-model.md#utrannotator).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `UTRAnnotator.pm`, lines 1285-1313: an insertion's position in the 5' UTR is the index of the base after it, in transcript order on both strands.
- Line 578: the lost-stop analysis is described as one that annotates a variant that "removes a stop codon of an existing uORF".
- Line 628: the lost-stop analysis skips a variant whose index lies more than two bases past the first base of the ORF's stop codon.
- Line 631: it skips an insertion whose index lies at or before that first base. The test at line 633, written for insertions, never holds: it asks for an undefined reference allele, and an insertion's is defined and empty (lines 1160-1161).
- Lines 639-642: for an insertion one or two bases past the stop codon's first base, the three bases at the stop codon's index are read from the changed sequence, and the stop codon is reported lost unless they spell one. An insertion written directly before the start codon has its index one before the stop codon's first base, so the lines are never reached for it.
- Lines 911-917: the frameshift analysis puts an upstream ORF's check point one index before the first base of its stop codon. Line 921 counts the ORF when the variant's index lies at or before that check point, which the index of such an insertion does.
- Lines 262-267: the row carries one term, the flagged one whose key the hash lists last.

**What correcting it would change.** The plugin would report the lost stop codon whichever way a record writes the inserted bases. An insertion written directly before a start codon would be read as the same bases written inside the stop codon are.

**The class it produces.** `utrannotator_lost_stop_written_before_stop_codon`: the `5UTR_consequence` and `5UTR_annotation` of an insertion's row on which Ensembl VEP writes the frameshift and vep-rs the lost stop of an upstream ORF that Ensembl VEP lists as shifted. The two fields are set aside together and counted once per row. The concordance figures score consequence sets and never see the fields, so nothing is set aside in them.

**How often the shape occurs.** The index, [Intended divergences from Ensembl VEP](../intended-divergences.md), places this entry at its tier 2: ordinary data carries the shape, and the difference needs a plugin. The shape a record alone shows is an insertion. Of the 15,646,481 records of the five chromosome 21 datasets the index counts shapes on, 996,643 are insertions (63,698 per million), and 93 of the golden corpora hold one. The 5' UTR sequence decides whether it is written before a start codon inside which an upstream ORF ends. Four committed flag sets run the plugin, two in each UTRAnnotator corpus.

gnomAD v4.1 carries the shape at `7:80633290-80633291` on `ENST00000956936`, an insertion of two bases. Ensembl VEP lists three ORFs as shifted on that row, the ones that begin 295, 274 and 263 bases before the coding sequence. vep-rs reports the lost stop of the one that begins 263 bases before it, whose next stop codon begins 220 bases before the coding sequence.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Two bases inserted directly before a start codon inside which an upstream ORF ends

**In short.** The constructed record `21 10535884 C CAT` inserts two bases between cDNA positions 360 and 361 of `ENST00001091162`, directly before an upstream `ATG` at cDNA position 361. The inserted bases repeat the codon's `A` and `T`, so the changed sequence is the reference with `AT` written after the codon's `T`, inside the `TGA` that closes an upstream ORF. In the changed sequence that ORF reads `TAT` where its stop codon stood and ends nine bases later. Ensembl VEP 116.2 reports the ORF's frameshift; vep-rs reports its lost stop codon.

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

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies `utrannotator_lost_stop_written_before_stop_codon` when a comparison case names it. It sets a row's `5UTR_consequence` and `5UTR_annotation` aside, counted once per row, when four things hold. Ensembl VEP's term is the frameshift and vep-rs's the lost stop. The row's two cDNA positions are an insertion's two flanking bases and, where the row's record states its alleles, its reference allele is empty once the first base its alleles share is dropped. Each group of vep-rs's annotation matches exactly one upstream ORF (`ref_type=uORF`) of Ensembl VEP's list, a different one for each group. vep-rs's row carries UTRAnnotator's three `Existing_` counts, which the plugin writes on every row it annotates; their values are not compared by this class.

The class reads nothing but the two engines' rows and, where the row's record states them, the record's alleles. A group matches an ORF when it carries that ORF's `KozakContext` and `KozakStrength` and the alternative stop codon the ORF's printed lengths give. Where the ORF's `alt_type` is `uORF`, the group reads `AltStop=True`, and its `AltStopDistanceToCDS` equals the ORF's `ref_StartDistanceToCDS` plus the allele's length, less its `alt_type_length`, plus three. Otherwise the group reads `AltStop=False` with `AltStopDistanceToCDS=NA`. On this row the listed ORF has `ref_StartDistanceToCDS=133` and `alt_type_length=45`, and the allele is two bases long: 133 + 2 - 45 + 3 = 93, the value of vep-rs's group.

The class leaves two things open. A row whose lost-stop group matches no ORF of Ensembl VEP's list, or more than one, stays a difference. So does a row whose alternative stop differs from the one the listed ORF's lengths give.

UTRAnnotator's classes are tried in a fixed order and a differing field is counted under the first that accepts it; this class comes last, after `utrannotator_frameshift_at_stop_boundary`. Two of the others can take a row of this shape. A row whose two terms Ensembl VEP's outputs under other hash seeds both write is counted under `utrannotator_consequence_precedence`, when the comparison is given those outputs. A row on which Ensembl VEP's hash order writes the lost start is counted under `utrannotator_lost_start_after_insertion`.

The `GRCh38-plugin-utr-annotator` and `GRCh37-plugin-utr-annotator` golden corpora carry no row of the shape, so the class applies to none of their rows.
