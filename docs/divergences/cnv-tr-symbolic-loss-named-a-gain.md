# A `<CNV:TR>` allele read symbolically that states a loss

This page shows the direction term, `feature_elongation` or `feature_truncation`, of a `<CNV:TR>` allele read symbolically that states a loss or its run's length.

**The defect in one sentence.** Ensembl VEP 116.3 reads every `<CNV:TR>` allele it keeps symbolic as a copy-number gain, so its direction term is `feature_elongation` for an allele whose record states fewer bases than the run it replaces, and for one whose record states exactly as many.

**Why this is a defect.** The two terms are opposites. `feature_elongation` (SO:0001907) is a variant "that causes the extension of a genomic feature, with regard to the reference sequence", and `feature_truncation` (SO:0001906) one "that causes the reduction" of it. A record that states six bases where the reference has nine shortens the exon that holds the run, and Ensembl VEP's row says the exon grows. Its own reading with the reference sequence writes a deletion for the same record. The term comes from the allele's class alone: the unit count, which the same parser reads when it has the sequence, is never asked.

A record that states nine bases where the reference has nine leaves the exon as long as it is, and Ensembl VEP's row says the exon grows, with the impact HIGH. With the reference sequence Ensembl VEP's own reading finds that allele equal to the run and writes no consequence for it ([2.3](#23-what-each-engine-writes)).

**What vep-rs does instead.** vep-rs reads two lengths from the record: the bases the allele states for itself and the bases of the run. Where the allele is the shorter, its row names `feature_truncation` in the place of `feature_elongation`. Both terms carry the impact HIGH, so IMPACT stands. Where the two lengths are equal and no field of the record states another, its row names neither direction term, and its IMPACT is that of the first term left by rank. Every other cell stands in both: the Location, the Allele, the region terms and the positions among them. Every other allele keeps the row Ensembl VEP writes: [1.5](#15-why-vep-rss-row-is-the-correct-one) and [2.5](#25-why-vep-rss-row-is-the-one-the-record-supports) say which.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`; `ensembl-variation` at release/116, `Bio/EnsEMBL/Variation/`):

- `Parser/VCF.pm`, `create_StructuralVariationFeatures`: line 502 opens the tandem-repeat branch. Lines 516-520 fetch the run's reference sequence and, with it, build the allele string in `_expand_tandem_repeat_allele_string` (lines 382-431, which read the repeat fields at lines 385-419). Lines 521-524 warn under `--offline` without `--fasta`, and no allele string is built.
- Lines 560-574 of the same file: without an allele string the record becomes a structural feature of class `tandem_repeat`. `Parser.pm` lines 699-701 read `CNV:TR` as the type `TREP`, whose term that is. Its allele string is the ALT column. The feature keeps the record's line (`_line`, line 573), and no repeat field is read from it.
- `BaseVariationFeatureOverlapAllele.pm`, lines 429-436: an allele whose class matches `tandem` gets the keys `insertion` and `increase_length`. `decrease_length` goes to a class matching `deletion` or `loss`, to a breakpoint class (lines 437-440) and to a copy-number class with a supporting deletion or loss (line 425). Lines 317-333, `_skip_oc`: a consequence whose keys the allele lacks is not evaluated, and `feature_truncation` asks for `decrease_length` (`Utils/Constants.pm`, lines 821-834).
- `Utils/VariationEffect.pm`, lines 246-256, `copy_number_gain`: true for every `tandem_repeat` (lines 282-309), with no test of a count. Lines 258-266, `copy_number_loss`: true only where the class, or the class of a supporting variant, matches `loss`.
- Lines 337-348, `feature_elongation`: true for a gain or an insertion that lies wholly inside the feature and for which `within_cdna` holds. Lines 350-373, `feature_truncation`: asks for a loss or a deletion.
- Lines 1616-1649, `coding_unknown`, the test of `coding_sequence_variant`: for a structural allele, true where the run lies in the coding sequence and none of `inframe_insertion`, `inframe_deletion` and `frameshift` holds (line 1643). It asks no gain.
- `Utils/Constants.pm`, lines 803-816 and 821-834: the two terms, both of impact HIGH, with ranks 9 and 10. Lines 1061-1078: `coding_sequence_variant`, of impact MODIFIER, with rank 23.
- `OutputFactory.pm` of `ensembl-vep`, lines 2043-2050: a row's terms are sorted by rank, and its IMPACT is the first term's.

**What correcting it would change.** Each symbolic `<CNV:TR>` row whose allele states a loss would name `feature_truncation`, and each whose allele states exactly the run's length would name no direction term, as vep-rs's rows do. In JSON a record's `most_severe_consequence` would follow where that row is the record's most severe. No other row would change.

**The class it produces.** `cnv_tr_symbolic_loss_named_a_gain`: a `Transcript` row whose Allele is `tandem_repeat` in both engines and whose Consequence holds `feature_elongation` in Ensembl VEP's row. In the first shape vep-rs's row reads `feature_truncation` in that term's place. In the second vep-rs's row holds no term there, and its IMPACT is that of the terms left. [1.7](#17-how-the-comparison-tools-treat-this-class) and [2.7](#27-how-the-comparison-tools-treat-this-shape) say how each comparison tool reads them.

**How often the shape occurs.** Tier 3 of [the index's order](../intended-divergences.md#differences-outside-the-classes). One committed corpus carries both shapes, `GRCh38-vcf-tandem-repeats`, whose fifteen records are hand-written. Three of its alleles state a loss: `tr_loss` (line 4 of its input), the second allele of `tr_multi` (line 8) and the second allele of `tr_two_svlen` (line 9). One states exactly its run's length: `tr_same` (line 5). Each has a row on two coding transcripts.

The corpus's own run reads no reference sequence and carries eight rows of the class in each output format: six of a stated loss and two of `tr_same`. Its flag set `fasta-max10` carries four: under `--max_sv_size 10` the two records of two alleles stay symbolic, their first allele being twelve bases, and `tr_same` is expanded. Its flag set `fasta` carries none. The release's structural-variant runs give both engines the FASTA, where the record is expanded and the engines agree.

## 1. A loss of one unit, read without the reference sequence

**In short.** The record `21 13620283 tr_loss G <CNV:TR>` with `SVLEN=9;RN=1;RUS=TCA;RUC=2` states a run of nine bases, `21:13620284-13620292`, and an alternate allele of two units of three bases. Read without a FASTA, Ensembl VEP 116.3 writes `feature_elongation,coding_sequence_variant` on the two coding transcripts that hold the run, and vep-rs writes `feature_truncation,coding_sequence_variant`. With the FASTA both engines write an `inframe_deletion` at `21:13620290-13620292`.

### 1.1 What a reader needs to know first

A `<CNV:TR>` record (VCF 4.4) states a change in the copy number of a tandem repeat without spelling the alternate allele out. POS is the base before the run and `SVLEN` the run's length. An alternate allele is `RN` repeat sequences. `RUS` is the unit of a repeat sequence, `RUC` its number of units and `RB` its length in bases. [The `<CNV:TR>` page](cnv-tr.md#121-what-a-reader-needs-to-know-first) explains the fields in full.

Both engines read such a record in one of two ways ([class 12](cnv-tr.md#12-cnvtr-read-symbolically-not-expanded)). With the run's reference sequence they rebuild the allele and annotate its difference from the run as an insertion or deletion of bases. Without it the allele stays symbolic: the row's Location is the run, its Allele is `tandem_repeat`, and its Consequence holds the region terms of the run. A direction term stands beside them where the run lies wholly inside the transcript and overlaps its cDNA.

A record stays symbolic on two routes in vep-rs. The first is `--offline` without `--fasta`, where no reference sequence is read. The second is `--max_sv_size`: with a FASTA, a record one of whose repeat sequences is longer than that many bases is not expanded, and every allele of it stays symbolic. Where a FASTA holds no sequence for the run, both engines expand the record against its REF column (`Parser/VCF.pm`, lines 516-520). Ensembl VEP has a third route, a database that gives it no sequence for the chromosome (lines 525-527); vep-rs connects to none. The page's statements hold on every route.

### 1.2 The record

The record is line 4 of `tests/golden/116/GRCh38-vcf-tandem-repeats/variants.vcf`, hand-written on a repeat the reference carries:

```
21	13620283	tr_loss	G	<CNV:TR>	.	.	SVLEN=9;RN=1;RUS=TCA;RUC=2
```

Its fields give two lengths. The run is POS + 1 to POS + `SVLEN`, `21:13620284-13620292`, nine bases. The alternate allele is one repeat sequence (`RN=1`) of two units (`RUC=2`) of the three-base unit `TCA`, six bases. So the record states six bases where the reference has nine.

The reference sequence agrees. In `tests/golden/116/GRCh38-identifiers/reference.fa.gz` the base at POS is the record's `G`, the nine bases of the run are three copies of `TCA`, and no further copy stands in phase on either side. The allele is the run less one unit.

The run lies in a coding exon of ENST00000299443 and ENST00000620442, both on the forward strand: cDNA positions 1496 to 1504, coding positions 1044 to 1052, residues 348 to 351. A third transcript, ENST00000384745, begins 1,285 bases after the run's last base.

### 1.3 What each engine writes

Both engines ran the record with the corpus's flags, `--uploaded_allele --allele_number --variant_class`, and no FASTA. Ensembl VEP 116.3 writes this row on ENST00000299443, stored in the corpus's `expected/default.txt.gz`:

```
tr_loss	21:13620284-13620292	tandem_repeat	ENSG00000166351	ENST00000299443	Transcript	feature_elongation,coding_sequence_variant	1496-1504	1044-1052	348-351	-	-	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;VARIANT_CLASS=tandem_repeat;OverlapBP=9;OverlapPC=0.02
```

vep-rs writes:

```
tr_loss	21:13620284-13620292	tandem_repeat	ENSG00000166351	ENST00000299443	Transcript	feature_truncation,coding_sequence_variant	1496-1504	1044-1052	348-351	-	-	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;VARIANT_CLASS=tandem_repeat;OverlapBP=9;OverlapPC=0.02
```

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| Consequence | `feature_elongation,coding_sequence_variant` | `feature_truncation,coding_sequence_variant` |
| `IMPACT` | `HIGH` | `HIGH` |

No other cell differs. The row on ENST00000620442 differs in the same term, and the row on ENST00000384745 reads `upstream_gene_variant` in both engines. In JSON the record's `most_severe_consequence` is `feature_elongation` in Ensembl VEP's output and `feature_truncation` in vep-rs's.

Ensembl VEP's three rows for `tr_gain`, line 3 of the same file, whose `RUC=4` states twelve bases, are its three rows for `tr_loss` but for the name. With the FASTA, the corpus's flag set `fasta`, both engines write `tr_loss` as a deletion: Location `21:13620290-13620292`, Allele `-`, `inframe_deletion`, `IMPACT=MODERATE`.

### 1.4 Why Ensembl VEP writes that

The parser builds a literal allele only with the run's reference sequence to hand. Without it the record becomes a structural feature whose class is `tandem_repeat` and whose allele is the ALT text. No repeat field is read from that feature.

The consequence terms are then decided from the class. An allele of a class matching `tandem` is marked as one that increases length, so `feature_truncation`, which asks for a decrease, is not evaluated at all. `feature_elongation` is, and its test of the direction is `copy_number_gain`, which is true of every `tandem_repeat`. So a gain of one unit and a loss of one unit get the same row.

A term of an earlier tier ends the evaluation (`BaseVariationFeatureOverlapAllele.pm`, lines 267 and 276-279). So `feature_elongation`, of the third tier, is not written where the run covers the whole transcript, from its first base to its last, or overlaps the mature sequence of a miRNA. Such a run names `transcript_amplification`, of the first tier (`Utils/VariationEffect.pm`, lines 330-335), or `mature_miRNA_variant`, of the second (`Utils/Constants.pm`, lines 785-801 and 1079-1098).

Without a FASTA Ensembl VEP warns `could not fetch sequence for tandem repeats; consequence calculation for tandem repeats is less precise when using --offline without --fasta`. On the `--max_sv_size` route it writes no warning.

Source (`ensembl-vep` at release/116.3 and `ensembl-variation` at release/116):

- `Parser/VCF.pm`, lines 516-524: the allele string built with the reference sequence, and the warning without it; lines 560-574: the symbolic feature.
- `Parser/VCF.pm`, lines 403 and 417-430: a repeat sequence longer than `--max_sv_size` ends the expansion with no allele string and no warning.
- `BaseVariationFeatureOverlapAllele.pm`, lines 429-436 and 317-333: the keys of a `tandem` class, and the consequences they leave out.
- `Utils/VariationEffect.pm`, lines 246-256 and 337-348: `copy_number_gain` and `feature_elongation`.

### 1.5 Why vep-rs's row is the correct one

The record's own fields decide it. Section 5.7 of the VCF 4.4 specification says: "The SVLEN of the `<CNV:TR>` is the length of the reference allele. It is not the length of the `<CNV:TR>` allele." Section 3 says: "Each `<CNV:TR>` allele consists of RN repeat sequences each containing RUC repeat units with sequence RUS." So the reference allele is nine bases and the alternate allele six. Section 5.7 calls the change a "repeat expansion/contraction", and its own example record carries an allele of 29 bases on a run of 30.

The reference sequence confirms it: the run is three copies of the unit and the allele two. The exon that holds the run is three bases shorter under this allele, which is the reduction of a feature that `feature_truncation` names. Ensembl VEP's own reading with the FASTA is a deletion of the run's last three bases. Its symbolic row states the opposite direction, which is not a less precise form of the same reading.

vep-rs changes that one term. The row keeps the run as its Location and `coding_sequence_variant` as its region term. It names no in-frame or frameshift effect from the lost length and places the loss on no base of the run.

vep-rs reads a loss where the record states the allele's length in whole units of a written unit:

- every ALT of the record is `<CNV:TR>`;
- the record carries `RUS` with `RUC` or `RB`, and neither `CIRUC` nor `CIRB`, the confidence intervals around them; a value of `0` or of no characters reads as an absent field (as at `Parser/VCF.pm` line 385);
- `RN` is absent, which gives each allele one repeat sequence, or holds one whole number for each allele;
- the allele has at least one repeat sequence, and the `RUS` list and the count list, read in the alleles' order, hold a unit and a count for each of its sequences;
- each unit is not `.`, and each count is a whole-number `RUC` or, in a record without `RUC`, an `RB` that is a whole multiple of the unit's length.

The allele's length is the sum over its sequences. The run is the row's Location: POS + 1 to the end that `SVLEN` gives, the largest where a record carries one for each allele, or `END` without `SVLEN`. A loss is a stated length under the run's length.

An allele that states exactly the run's length is read as [section 2](#2-an-allele-as-long-as-its-run-read-without-the-reference-sequence) shows. Every other allele keeps Ensembl VEP's row: one that states more bases than the run, and one whose length the record does not state so. Some of the last kind still suggest a loss, such as a fractional `RUC` whose product is under the run, or an interval around the count. Leaving them as Ensembl VEP writes them is vep-rs's choice of what it takes as proven.

The terms of a run that covers a whole feature stay too: `transcript_amplification`, `regulatory_region_amplification` and `TFBS_amplification`. No committed record that states a loss reaches one of them.

### 1.6 Try it yourself

The test file `crates/vep-cli/tests/tandem_repeat_loss_read_symbolically.rs` runs the corpus through the `vep` binary, without a FASTA and with one under `--max_sv_size 10`. It holds every row to Ensembl VEP's stored row, apart from the twelve rows of the class, ten of a stated loss and the two of `tr_same`:

```
cargo test -p vep-cli --test tandem_repeat_loss_read_symbolically
```

To run both engines, write the record of 1.2 to `record.vcf` under the two header lines of the corpus's input. The commands are those of [Reproducing these examples](reproducing.md#the-commands) without `--fasta` and with the corpus's three flags. vep-rs, from the repository's root, on the corpus's committed cache:

```
vep -i record.vcf -o vep-rs.txt --offline \
  --json_cache tests/golden/116/GRCh38-vcf-tandem-repeats/json_cache \
  --species homo_sapiens --assembly GRCh38 \
  --uploaded_allele --allele_number --variant_class \
  --buffer_size 5000 --force --no_stats
```

Ensembl VEP 116.3, on its GRCh38 cache of release 116:

```
docker run --rm -v "$PWD":/work -v /path/to/vep_cache:/cache \
  ensemblorg/ensembl-vep:release_116.3 \
  vep -i /work/record.vcf -o /work/vep.txt --offline --cache \
  --dir_cache /cache --cache_version 116 \
  --species homo_sapiens --assembly GRCh38 --format vcf \
  --uploaded_allele --allele_number --variant_class \
  --buffer_size 5000 --fork 4 --no_stats --force_overwrite
```

`grep ENST00000299443` selects the row of 1.3 in each output. With `--fasta` naming the FASTA of the corpus `GRCh38-identifiers`, both engines write the deletion, as the rows of the corpus's flag set `fasta` show.

### 1.7 How the comparison tools treat this class

The reproduction `scripts/divergences/reproductions/cnv_tr_symbolic_loss_named_a_gain.toml` runs the records of 1.2 and [2.2](#22-the-record) on both engines over the corpus's cache without a FASTA. It records the three rows of each engine for each record. A further run of the reproduction takes line 8 of the corpus's input with the FASTA under `--max_sv_size 10` and records the rows of its two alleles on ENST00000299443.

The format comparator, `tools/vep-formats-compare`, has the class built in and recognises it from the rows alone. The row is a `Transcript` row whose Allele is `tandem_repeat` on both sides. Its two lists of terms are equal place for place but for one place, where Ensembl VEP writes `feature_elongation` and vep-rs `feature_truncation`. The comparator sets the difference of the Consequence aside and compares every other cell; in JSON it leaves the `most_severe_consequence` of such a record uncompared. A row with another Allele, a pair that differs in any other term, and the two terms the other way round stay differences.

The corpus tool, `tools/vep-golden-corpus`, reads the same pair when it classifies a corpus's rows and writes an entry under the class. The golden harness applies the entries of the corpus's manifest, which names the class in its `masks`. The manifest holds four entries: one for each of the two coding transcripts, in the corpus's own run and in `fasta-max10`.

The structural-variant comparator, `tools/vep-sv-compare`, has no class of this name and reads the rows under its `<CNV:TR>` pairing ([12.7](cnv-tr.md#127-how-the-comparator-treats-this-class)). Where Ensembl VEP's run read a FASTA and vep-rs's did not, vep-rs's row of a stated loss pairs with the reference's deletion, and `feature_truncation` is the direction term set aside. Where neither run read a FASTA, the two symbolic rows pair, since they differ in the direction term alone.

## 2. An allele as long as its run, read without the reference sequence

**In short.** The record `21 13620283 tr_same G <CNV:TR>` with `SVLEN=9;RN=1;RUS=TCA;RUC=3` states a run of nine bases, `21:13620284-13620292`, and an alternate allele of three units of three bases. Read without a FASTA, Ensembl VEP 116.3 writes `feature_elongation,coding_sequence_variant` with `IMPACT=HIGH` on the two coding transcripts that hold the run, and vep-rs `coding_sequence_variant` with `IMPACT=MODIFIER`. With the FASTA both engines find the allele equal to the run and write no default or tab row.

### 2.1 What a reader needs to know first

The fields and the two readings are those of [1.1](#11-what-a-reader-needs-to-know-first).

A row's IMPACT is that of its first term by rank (`OutputFactory.pm`, lines 2043-2050). `feature_elongation` has rank 9 and the impact HIGH, `coding_sequence_variant` rank 23 and the impact MODIFIER (`Utils/Constants.pm`, lines 803-816 and 1061-1078).

With the reference sequence, an expanded allele that equals the run is no allele to annotate (`VariationFeatureOverlap.pm`, line 481). [Class 19](defect-l.md#19-an-empty-or-unchanged-cnvtr-expansion-dropped-without-a-row) sets aside the symbolic rows vep-rs writes for such a record without a FASTA where Ensembl VEP, with one, wrote none.

### 2.2 The record

The record is line 5 of `tests/golden/116/GRCh38-vcf-tandem-repeats/variants.vcf`, hand-written on the run of [1.2](#12-the-record):

```
21	13620283	tr_same	G	<CNV:TR>	.	.	SVLEN=9;RN=1;RUS=TCA;RUC=3
```

Its fields give two lengths. The run is POS + 1 to POS + `SVLEN`, `21:13620284-13620292`, nine bases. The alternate allele is one repeat sequence (`RN=1`) of three units (`RUC=3`) of the three-base unit `TCA`, nine bases. So the record states nine bases where the reference has nine. It writes no `END`, `RB`, `RUL`, `RUB` or `CN`, so no field of it states another length.

The reference sequence agrees. In `tests/golden/116/GRCh38-identifiers/reference.fa.gz` the nine bases of the run are three copies of `TCA`, so this allele is the run itself. The transcripts and the positions are those of 1.2.

### 2.3 What each engine writes

Both engines ran the record with the corpus's flags and no FASTA. Ensembl VEP 116.3 writes this row on ENST00000299443, stored in the corpus's `expected/default.txt.gz`:

```
tr_same	21:13620284-13620292	tandem_repeat	ENSG00000166351	ENST00000299443	Transcript	feature_elongation,coding_sequence_variant	1496-1504	1044-1052	348-351	-	-	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;VARIANT_CLASS=tandem_repeat;OverlapBP=9;OverlapPC=0.02
```

vep-rs writes:

```
tr_same	21:13620284-13620292	tandem_repeat	ENSG00000166351	ENST00000299443	Transcript	coding_sequence_variant	1496-1504	1044-1052	348-351	-	-	-	ALLELE_NUM=1;IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=tandem_repeat;OverlapBP=9;OverlapPC=0.02
```

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| Consequence | `feature_elongation,coding_sequence_variant` | `coding_sequence_variant` |
| `IMPACT` | `HIGH` | `MODIFIER` |

No other cell differs. The row on ENST00000620442 differs in the same two cells, and the row on ENST00000384745 reads `upstream_gene_variant` in both engines. In JSON the record's `most_severe_consequence` is `feature_elongation` in Ensembl VEP's output and `coding_sequence_variant` in vep-rs's.

Ensembl VEP's three rows for `tr_same` are its three rows for `tr_gain` and for `tr_loss` but for the name. With the FASTA, in the corpus's flag sets `fasta` and `fasta-max10`, neither engine writes a default or tab row for the record. Its VCF line stands without consequences, and its JSON object with the `most_severe_consequence` `?`.

### 2.4 Why Ensembl VEP writes that

The mechanism is that of [1.4](#14-why-ensembl-vep-writes-that): the count is never asked, so three units on a run of three get the row of a gain.

Ensembl VEP's test of a tandem repeat on a sequence allele is narrower. There `tandem_repeat` holds only for an alternate longer than its reference and made of exact copies of it (`Utils/VariationEffect.pm`, lines 287-298). On a structural feature it asks the class alone (lines 301-307).

`coding_sequence_variant` stands because its test, `coding_unknown`, asks no gain. For a structural allele it holds where the run lies in the coding sequence and no in-frame or frameshift term does (line 1643). `inframe_deletion` and `frameshift` ask for a deletion or a loss (lines 1192 and 1573-1574), and `inframe_insertion` is false of every structural allele (lines 1132-1135).

### 2.5 Why vep-rs's row is the one the record supports

The record's own fields decide it. By the sentences of the specification quoted in [1.5](#15-why-vep-rss-row-is-the-correct-one), the reference allele is the nine bases of the run and the alternate allele three units of three bases, nine. The exon that holds the run is as long under this allele as in the reference. `feature_elongation` names an extension of a feature and `feature_truncation` a reduction, so neither holds.

Section 5.7 says the copy number of a `<CNV:TR>` allele "must correspond to the sample allelic length divided by the reference allelic length", here nine over nine. Section 5.6 says of copy-number alleles: "GT=0 is equivalent to a `<CNV>` ALT allele with INFO CN of 1 and should be treated identically." For copy number, an allele as long as its run is the reference allele.

That could be a reason to write no row, as both engines do with the FASTA. This page does not claim it where the sequence is not read. With the sequence the engines find the allele equal to the run, base for base. Without it the record gives no base of the run, so vep-rs cannot tell the allele from one of other bases, and a record the input states is annotated ([19.5](defect-l.md#195-why-vep-rss-rows-are-the-correct-ones)). So vep-rs writes the row the record alone supports: it names where the run lies, by its Location and the region term `coding_sequence_variant`, and names no change of length.

That the allele adds no base is a claim about the record's fields, not about a sample: section 5.7 lets a `<CNV:TR>` record that stands beside sequence records of the same repeat leave out small changes.

vep-rs takes the equal length as proven where all of these hold:

- the record states the allele's length by the conditions of [1.5](#15-why-vep-rss-row-is-the-correct-one), and that length is the run's;
- the record writes `SVLEN` or `END`, and the run is longer than one base: for a repeat that the reference lacks, section 5.7 sets `SVLEN` to 1;
- every `SVLEN` value the record writes is the same number, and an `END` beside `SVLEN` names the run's end: `SVLEN` is one value for each allele (section 3), and Ensembl VEP takes the largest (`Parser/VCF.pm`, lines 504-511);
- no other field may state another length for the allele. An `RB` beside `RUC` must equal the unit's length times the count for each of its sequences. An `RUL`, the unit's length, must equal the length of its `RUS`. The record must write no `RUB`, the bases of each single unit, and an INFO `CN`, where written, must be 1 for the allele.

An allele that fails any of them keeps Ensembl VEP's row. Only the INFO fields are read: no sample column changes the reading.

Lines 732, 750 and 768 of the corpus `GRCh37`'s `variants.vcf` write an `RB` equal to `SVLEN` and still state no length by the conditions of 1.5: their `RUC` is fractional, and their `RN` counts more repeat sequences than their lists hold. They keep Ensembl VEP's rows.

The amplification terms of a run that covers a whole feature (`transcript_amplification`, `regulatory_region_amplification`, `TFBS_amplification`) are left as Ensembl VEP writes them, as a choice of scope: no committed record of this shape reaches one.

### 2.6 Try it yourself

The test file `crates/vep-cli/tests/tandem_repeat_as_long_as_its_run_read_symbolically.rs` runs the corpus through the `vep` binary, without a FASTA and with one. It holds the rows and the JSON object of `tr_same` to the values of 2.3:

```
cargo test -p vep-cli --test tandem_repeat_as_long_as_its_run_read_symbolically
```

To run both engines, write the record of 2.2 to `record.vcf` under the two header lines of the corpus's input and run the two commands of [1.6](#16-try-it-yourself). `grep ENST00000299443` selects the row of 2.3 in each output. With `--fasta` naming the FASTA of the corpus `GRCh38-identifiers`, neither output holds a row for the record.

### 2.7 How the comparison tools treat this shape

The format comparator reads the shape under the same class: the two lists of terms are equal but for Ensembl VEP's `feature_elongation`, where vep-rs writes no term. It sets aside the Consequence and the IMPACT each side's terms give, and leaves the JSON record's `most_severe_consequence` uncompared.

The corpus tool writes the pair under the class: the manifest's two entries of the corpus's own run name `tr_same` beside the stated losses.

The structural-variant comparator needs no rule of its own. Where Ensembl VEP read a FASTA and vep-rs did not, Ensembl VEP wrote nothing for the record, and vep-rs's rows at the run are set aside under [class 19](defect-l.md#197-how-the-comparator-treats-this-class). Where neither read one, the two symbolic rows pair: they differ in `feature_elongation` alone, since this comparator reads no IMPACT.
