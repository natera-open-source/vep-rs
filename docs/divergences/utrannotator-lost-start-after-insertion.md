# A start codon left whole that UTRAnnotator reports lost

**The defect in one sentence.** Ensembl VEP's UTRAnnotator plugin reports an upstream start codon lost when an insertion only moves the codon on: the codon stands whole after the inserted bases, with every base after it unchanged. vep-rs reports no lost start codon for such an insertion, whichever way the record writes it, and does not report the moved codon gained.

**Why this is a defect.** An insertion lies between two bases. When the base after it is the `A` of an upstream `ATG`, the insertion touches no base of the codon: in the changed 5' UTR the codon stands whole, as many bases further from the cap as the insertion is long. The plugin's lost-start analysis is meant to report a start codon the variant removes, and nothing was removed.

A record can write the same changed sequence a second way. Bases inserted after the codon's `A` that end with `A`, or after its `T` that end with `AT`, give the sequence of the first way, base for base. Read from the cap, the inserted bases and the codon's first base or two trade places, and the codon again stands whole further on. A VCF writes a repeat at its leftmost place on the genome, which on a reverse-strand transcript is its last place on the transcript, so there such an insertion is written inside the codon.

**What vep-rs does instead.** vep-rs reports no lost start codon for an insertion that only moves the codon on, whichever way the record writes it, and does not report the moved codon gained. The row then carries whichever of the other analyses fires, under the fixed order described on [the precedence page](utrannotator-consequence-precedence.md), or `-` when none does; those analyses read the insertion directly before the codon. After any other insertion that follows the codon's first or second base, both engines test the codon by the three bases at its old index. vep-rs's row carries the lost start when that test fires and no analysis earlier in the fixed order does. The plugin's parameters, fields and analyses are described under [UTRAnnotator](../plugins/transcript-model.md#utrannotator).

An upstream ORF whose stop codon begins at the moved codon's second base is read from the changed sequence at that index. Where a stop codon still stands there, the ORF ends where it ended and the frameshift leaves it out. Where none does, the insertion's length is no multiple of three and the inserted bases can be written inside that stop codon, the stop codon is lost, as [the page on that row](utrannotator-lost-stop-written-before-stop-codon.md) shows. With no stop codon there, a length that is a multiple of three shifts no frame, and vep-rs reports neither a lost stop nor a frameshift for that ORF; any other length shifts it, as on [the page on Ensembl VEP's frameshift list](utrannotator-frameshift-at-stop-boundary.md).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `UTRAnnotator.pm`, lines 1285-1313: an insertion's position in the 5' UTR is the index of the base after it, in transcript order on both strands.
- Lines 764-769: a start codon is tested when that index lies at the codon's first, second or third base. Line 769 skips an allele that starts before the codon's first base.
- Lines 776-778: for an allele that is not a deletion, the three bases at the codon's old index are read from the changed sequence, and the codon is reported lost unless they spell `ATG`. After an insertion that only moves the codon on, those three bases hold inserted bases, and the codon is never looked for where it stands.

**What correcting it would change.** The plugin would read the codon where the insertion moved it, and would report the loss only for a codon whose bases the variant parts or replaces.

**The class it produces.** `utrannotator_lost_start_after_insertion`: a `5UTR_consequence` that is the lost-start term on Ensembl VEP's row alone, where the row's `5UTR_annotation` names one codon and the insertion only moves that codon on. It is set aside with the row's `5UTR_annotation` and counted once per row. The concordance figures score consequence sets and never see the field, so nothing is set aside in them.

**How often the shape occurs.** The index, [Intended divergences from Ensembl VEP](../intended-divergences.md), places this entry at its tier 2: ordinary data carries the shape, and the difference needs a plugin. The shape a record alone shows is an insertion. Of the 15,646,481 records of the five chromosome 21 datasets the index counts shapes on, 996,643 are insertions (63,698 per million), and 132 of the golden corpora hold one. The 5' UTR sequence decides whether it only moves a start codon on. Four committed flag sets of the golden corpora run the plugin, two in each UTRAnnotator corpus. On the genome-wide ClinVar files Ensembl VEP 116.3 reports such a codon lost on 89 rows (GRCh37) and 210 rows (GRCh38).

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. Three bases inserted directly before an upstream start codon

**In short.** The constructed record `21 10524582 A ACCC` inserts three bases between cDNA positions 145 and 146 of `ENST00000427445`. The upstream ORF's `ATG` begins at cDNA position 146 and is untouched. Ensembl VEP 116.3 reports that start codon lost; vep-rs reports nothing, since the codon stands whole three bases on and no other analysis fires.

### 1.1 What a reader needs to know first

UTRAnnotator runs five analyses over a 5' UTR variant: a start codon gained (`uAUG_gained`) or lost (`uAUG_lost`), a uORF's stop gained (`uSTOP_gained`) or lost (`uSTOP_lost`), and a frameshift in a uORF (`uFrameshift`). The row carries one flag in `5UTR_consequence` and that flag's effects in `5UTR_annotation`. A lost start codon is reported with `CapDistanceToStart`, the number of 5' UTR bases before the codon, which is one less than the cDNA position of its `A`.

### 1.2 The record

A constructed GRCh38 insertion in the 5' UTR of `ENST00000427445` (TPTE, forward strand) on chromosome 21, against the Ensembl 116 GRCh38 cache and the published GRCh38 uORF table:

```
21	10524582	utr_ins_before_start	A	ACCC	.	.	.
```

Both engines run `--plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt` with `--fasta`. The record is constructed because no ClinVar record on chromosome 21, the one chromosome the committed reference sequences cover, writes an insertion directly before an upstream start codon.

### 1.3 What each engine writes

The row on `ENST00000427445`, one of nine transcripts whose rows differ the same way. Ensembl VEP:

```
utr_ins_before_start	21:10524582-10524583	CCC	ENSG00000274391	ENST00000427445	Transcript	5_prime_UTR_variant	145-146	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;5UTR_annotation=KozakContext=AGAATGT:type=uORF:KozakStrength=Moderate:DistanceToStop=48:DistanceToCDS=207:CapDistanceToStart=145:Evidence=False;5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=6
```

vep-rs:

```
utr_ins_before_start	21:10524582-10524583	CCC	ENSG00000274391	ENST00000427445	Transcript	5_prime_UTR_variant	145-146	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;5UTR_annotation=-;5UTR_consequence=-;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=6
```

The `Existing_` counts agree; the order of the `5UTR_annotation` items is Ensembl VEP's hash order, which the comparison tools read as a set.

### 1.4 Why Ensembl VEP writes that

The transcript's 5' UTR reads `CGAGAGAATGTTGGA` around the site, with the `A` of the start codon 145 bases from the cap. The record places three bases between the `A` at cDNA position 145 and that `A`, giving `CGAGAGACCCATGTTGGA`.

`UTRAnnotator.pm`, lines 1285-1313, give the insertion the index of the base after it, 145, which is the codon's own first index. Line 769 skips an allele only when its index lies before the codon, so this one is tested. Line 776 then reads the three bases at index 145 of the changed sequence, `CCC`, the inserted bases, and line 778 reports the codon lost because they do not spell `ATG`.

### 1.5 Why vep-rs's row is the correct one

The changed sequence decides. Its bases 148 to 150 read `ATG`: the codon the reference has at 145, moved on by the three inserted bases, followed by the same base. No base of the codon was removed or replaced, so no start codon is lost.

The insertion changes what precedes the codon, and the other analyses report what follows from that. Three bases shift no frame, the inserted bases hold no start codon and no stop codon is touched, so nothing fires on this row and vep-rs writes `-`. Where such an insertion lies inside another upstream ORF and its length is no multiple of three, vep-rs reports that ORF's frameshift, or its lost stop where the inserted bases can be written inside its stop codon. An ORF that still ends on a stop codon at its own index is not shifted and is not reported.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh38 --plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt
```

Ensembl VEP writes `5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant` with `CapDistanceToStart=145` on `ENST00000427445`; vep-rs writes `5UTR_consequence=-`.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) applies `utrannotator_lost_start_after_insertion` when a comparison case names it. It sets a row's `5UTR_consequence` and `5UTR_annotation` aside, counted once per row, when three things hold. Ensembl VEP's value is the lost-start term and its annotation names one `CapDistanceToStart`. The row's two cDNA positions are an insertion's two flanking bases, and that number is the lower of them, or lies one below it on a row whose allele ends with `A`, or two below it on a row whose allele ends with `AT`. vep-rs's row carries the plugin's counts and does not report a lost start.

The class reads the row alone. The lower cDNA position of an insertion's row is the 5' UTR index of the base after the insertion, so a codon reported at that index begins directly after it. A codon reported one or two below it holds the insertion after its first or second base, and the allele's last bases say whether it stands whole further on. The allele is read on the transcript's strand, reverse-complemented where the row's `STRAND` is `-1`; a row that states no strand is read at the codon's own index alone.

A deletion's or a substitution's row reports an index below its first position (`UTRAnnotator.pm`, lines 764-769). Where a row's record states its alleles, a reference allele that is not empty refuses the row.

UTRAnnotator's classes are tried in a fixed order and a differing field is counted under the first that accepts it; this class follows the gained-start window and the precedence classes, and precedes `utrannotator_gained_start_remade_in_place`. The first two can take a row of this shape. When vep-rs reports a start codon gained inside the inserted bases, the row is counted under `utrannotator_gained_start_window`, provided the comparison is given the run's input file and reads the record's alleles from it. When vep-rs reports a frameshift, the row is counted under `utrannotator_consequence_precedence`, provided the comparison is given Ensembl VEP's outputs of the same input under other Perl hash seeds and one of them reports the frameshift on the row.

The `GRCh38-plugin-utr-annotator` and `GRCh37-plugin-utr-annotator` golden corpora carry no row of the shape, so the class applies to none of their rows.

## 2. Three bases inserted after the first base of an upstream start codon

**In short.** The constructed record `21 14210906 A ATGG` inserts three bases between cDNA positions 49 and 50 of `ENST00000400211`, a reverse-strand transcript, after the `A` of the upstream ORF's `ATG` at cDNA position 49. Read on the transcript's strand the inserted bases end with `A`, so the changed sequence is the reference with three bases directly before the codon. Ensembl VEP 116.3 reports that start codon lost; vep-rs reports nothing.

### 2.1 What a reader needs to know first

The analyses and the two fields are described in 1.1. A row's `Allele` is written on the genome's forward strand: on a reverse-strand transcript the transcript reads its reverse complement.

### 2.2 The record

A constructed GRCh38 insertion in the 5' UTR of `ENST00000400211` (reverse strand) on chromosome 21, against the cache and the uORF table of 1.2:

```
21	14210906	utr_ins_inside_start	A	ATGG	.	.	.
```

Both engines run the command of 1.2. The record is constructed: the one ClinVar record on chromosome 21 whose insertion only moves a start codon on from inside it is written after the codon's second base, and this section shows the insertion after the first.

### 2.3 What each engine writes

The row on `ENST00000400211`, one of six transcripts whose rows differ the same way. Ensembl VEP:

```
utr_ins_inside_start	21:14210906-14210907	TGG	ENSG00000188992	ENST00000400211	Transcript	5_prime_UTR_variant	49-50	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;5UTR_annotation=KozakContext=TCAATGA:type=uORF:KozakStrength=Weak:DistanceToStop=30:DistanceToCDS=62:CapDistanceToStart=48:Evidence=False;5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;Existing_InFrame_oORFs=1;Existing_OutOfFrame_oORFs=0;Existing_uORFs=1
```

vep-rs:

```
utr_ins_inside_start	21:14210906-14210907	TGG	ENSG00000188992	ENST00000400211	Transcript	5_prime_UTR_variant	49-50	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;5UTR_annotation=-;5UTR_consequence=-;Existing_InFrame_oORFs=1;Existing_OutOfFrame_oORFs=0;Existing_uORFs=1
```

### 2.4 Why Ensembl VEP writes that

The transcript's 5' UTR reads `AACCTCAATGAGTTA` around the site, with the `A` of the start codon 48 bases from the cap. The record's `TGG` reads `CCA` on the transcript's strand and lands after that `A`, giving `AACCTCAACCATGAGTTA`.

`UTRAnnotator.pm`, lines 1285-1313, give the insertion the index of the base after it, 49, the index of the codon's second base, and lines 764-769 let it through. Line 776 then reads the three bases at index 48 of the changed sequence, `ACC`, and line 778 reports the codon lost because they do not spell `ATG`.

### 2.5 Why vep-rs's row is the correct one

The changed sequence decides. Its bases 51 to 53 read `ATG`, followed by every base that follows the codon in the reference: it is the sequence that `ACC` inserted directly before the codon gives, and a forward-strand record would write it that way. The codon was moved on, not removed, and the inserted bases make no start codon: the changed 5' UTR holds as many as the reference.

Three bases shift no frame and no stop codon is touched, so nothing fires on this row and vep-rs writes `-`, as it does for the record written the other way.

### 2.6 Try it yourself

The command of 1.6 on this record. Ensembl VEP writes `5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant` with `CapDistanceToStart=48` on `ENST00000400211`; vep-rs writes `5UTR_consequence=-`.

### 2.7 How the comparator treats this record

Ensembl VEP's annotation names `CapDistanceToStart=48`, one below the row's lower cDNA position, 49. The row's `Allele`, `TGG`, reads `CCA` on the transcript's strand and ends with `A`, so the class sets the row's two fields aside. Had the inserted bases neither ended with `A` nor opened with `TG`, the codon would be parted, and the lost-start analysis of both engines would fire.
