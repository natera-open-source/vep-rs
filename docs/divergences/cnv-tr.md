# A representation difference, not a defect: `<CNV:TR>`

**In one sentence.** A `<CNV:TR>` record describes a change in the copy number of a tandem repeat without spelling the alternate allele out, and there are two equally faithful ways to read it. Each engine gives the same reading on the same input; the class pairs the readings of two runs that read the record differently.

**The two readings.** With a FASTA, both engines materialise the allele: they repeat the unit the stated number of times, read the run's reference bases from the FASTA and annotate the difference as an ordinary insertion or deletion of bases. Without a FASTA, both keep the allele symbolic, as a copy-number change of the run, and annotate the run's span as a gain whatever the stated count. Both give the symbolic reading when the expanded allele would exceed `--max_sv_size`.

**Why neither is wrong.** The expanded row is a faithful annotation of a real molecule: a run of 121 units where the reference has nine is, base for base, an insertion. The symbolic row keeps the record in the terms its author chose: a copy-number change of a repeat, located on the repeat, with the region term and the direction term. Nothing in either row is false, so the comparator pairs the two readings rather than ranking them when one run read a FASTA and the other did not.

**Where the readings part** (`ensembl-vep`): `Bio/EnsEMBL/VEP/Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`, lines 382-431, called from `create_StructuralVariationFeatures` at line 520.

**The class.**
- [Class 12](#12-cnvtr-read-symbolically-not-expanded): the synthetic tandem-repeat set on both assemblies, 1,368 and 3,751 pairs. Set aside: both members of each pair.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 12. `<CNV:TR>` read symbolically, not expanded

**In short.** A `<CNV:TR>` record states a copy-number change of a tandem repeat without spelling the alternate allele out. A run that reads a FASTA rebuilds that allele and annotates its difference from the reference run as an insertion or deletion of bases. A run without a FASTA keeps the allele symbolic: its Location is the run, its Allele `tandem_repeat`, its consequence the region the run lies in, plus `feature_elongation` where the run touches exonic sequence. Both engines read the record both ways on the same terms, and both give the symbolic reading when the expanded allele would exceed `--max_sv_size`. Nothing in either row is false, so the comparator pairs the readings rather than ranking them.

The rows below are the expanded reading on Ensembl VEP's side and the symbolic reading on vep-rs's, the two runs of the paper's measurement: the reference ran with the FASTA and vep-rs without one.

### 12.1 What a reader needs to know first

A `<CNV:TR>` record (VCF 4.4) describes a change in the copy number of a tandem repeat without spelling the alternate allele out. Its fields are:

- POS is the base before the repeat run;
- `INFO/END` is the run's last base;
- `SVLEN` is its length;
- `RUS` is the repeat unit;
- `RUL` is the unit's length;
- `RUC` is the number of units in the alternate allele;
- `RB` is the alternate allele's length in bases;
- `RN=1` says the alternate allele consists of one repeat sequence.

The record `T <CNV:TR>` with `END=25896436;SVLEN=18;RUS=CA;RUC=121;RB=242` says: the run of nine CA units at 25896419-25896436 carries 121 units in the alternate allele, a gain of 112 units (224 bases).

There are two ways to annotate such a record. One is to reconstruct the alternate sequence, compare it with the reference run and annotate the difference as an ordinary insertion or deletion of bases. The other is to keep the allele symbolic, as a copy-number change of the run itself, and annotate the run's span for its direction. Both engines do the first with a FASTA to read the run from and the second without one (`Parser/VCF.pm` lines 516-527); the rows below set a run with a FASTA against one without.

The terms involved are:

- `intron_variant` (SO:0001627).
- `feature_elongation` (SO:0001907): "a sequence variant that causes the extension of a genomic feature, with regard to the reference sequence".
- `coding_sequence_variant` (SO:0001580): "a sequence variant that changes the coding sequence".
- `inframe_insertion` (SO:0001821): "an inframe non synonymous variant that inserts bases into in the coding sequence".

Both engines follow Ensembl's rule for a span allele (a deletion, duplication, insertion or repeat): `feature_elongation` and `feature_truncation` are written only when the span touches the transcript's exonic sequence. The rule is in `ensembl-variation`, `Utils/VariationEffect.pm` lines 337-348 and 350-373, the `within_cdna` tests at lines 344 and 362. So a run inside an intron gets the region term alone, and a run inside an exon gets the region term plus the direction term. A breakend's `feature_truncation` is a different test, the breakpoint falling anywhere inside the transcript, intron included (lines 357-358; [class 11](defect-e.md#11-breakend-mate-rows-read-from-the-other-chromosomes-coordinate)).

### 12.2 The record

Both records are from the synthetic tandem-repeat set `tests/sv_validation/grch38/07_cnv_repeat.vcf.gz` (GRCh38; the GRCh37 set is `tests/sv_validation/grch37/07_cnv_repeat.vcf.gz`). Every `<CNV:TR>` record of those sets sits on a tandem repeat the chromosome 21 reference actually carries. The first is a gain in an intron:

```
21	25896418	synth_cnv_tr_0001	T	<CNV:TR>	.	PASS	SVTYPE=CNV;END=25896436;SVLEN=18;RN=1;RUS=CA;RUL=2;RUC=121;RB=242
```

The reference at 21:25896410-25896445 reads `AGTAACACT` `CACACACACACACACACA` `AAACAAAAC` (9 + 18 + 9 bases). The run is the 18 bases 25896419-25896436 (nine CA units), and the alternate has 121 units, a gain of 112 units or 224 bases.

The transcript is ENST00000346798 (APP, reverse strand, canonical, 25880550-26170770); the run lies in its intron 16 (25891869-25897572). The record has 41 rows on each side, one per APP transcript within reach of the run.

The second is a gain inside a coding exon:

```
21	26021865	synth_cnv_tr_0023	T	<CNV:TR>	.	PASS	SVTYPE=CNV;END=26021886;SVLEN=21;RN=1;RUS=GTG;RUL=3;RUC=46;RB=138
```

The run is the 21 bases 26021866-26021886 (seven GTG units, reference `GACTCTGTGGTGGTGGTGGTGGTGGTGGCAATGCTG` at 26021860-26021895). The alternate has 46 units, a gain of 39 units or 117 bases, a multiple of three. On ENST00000346798 the run lies in exon 6 (26021840-26022042), coding positions 819-839.

### 12.3 What each engine writes

Ensembl VEP, record 1, and vep-rs, record 1, one row each, the differing cells in bold:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `synth_cnv_tr_0001` | `synth_cnv_tr_0001` |
| Location | **`21:25896436-25896437`** | **`21:25896419-25896436`** |
| Allele | **`CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA`** | **`tandem_repeat`** |
| Gene | `ENSG00000142192` | `ENSG00000142192` |
| Feature | `ENST00000346798` | `ENST00000346798` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `intron_variant` | `intron_variant` |
| cDNA_position | `-` | `-` |
| CDS_position | `-` | `-` |
| Protein_position | `-` | `-` |
| Amino_acids | `-` | `-` |
| Codons | `-` | `-` |
| Existing_variation | `-` | `-` |
| Extra | **`IMPACT=MODIFIER;STRAND=-1`** | **`IMPACT=MODIFIER;STRAND=-1;OverlapBP=18;OverlapPC=0.01`** |

Ensembl VEP, record 2, and vep-rs, record 2, one row each, the differing cells in bold:

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `synth_cnv_tr_0023` | `synth_cnv_tr_0023` |
| Location | **`21:26021886-26021887`** | **`21:26021866-26021886`** |
| Allele | **`GTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTG`** | **`tandem_repeat`** |
| Gene | `ENSG00000142192` | `ENSG00000142192` |
| Feature | `ENST00000346798` | `ENST00000346798` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | **`inframe_insertion`** | **`feature_elongation,coding_sequence_variant`** |
| cDNA_position | **`968-969`** | **`969-989`** |
| CDS_position | **`818-819`** | **`819-839`** |
| Protein_position | **`273`** | **`273-280`** |
| Amino_acids | **`A/ATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT`** | **`-`** |
| Codons | **`gcc/gcCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACc`** | **`-`** |
| Existing_variation | `-` | `-` |
| Extra | **`IMPACT=MODERATE;STRAND=-1`** | **`IMPACT=HIGH;STRAND=-1;OverlapBP=21;OverlapPC=0.01`** |

On the first record the two engines agree on the consequence: the run is intronic, so neither writes a direction term. The whole difference is the Location and the Allele. Ensembl VEP's Location is an insertion between the run's last base and the next, and its Allele is the 224 inserted bases; vep-rs's Location is the run, and its Allele is the symbolic class `tandem_repeat`. The comparator keys on both columns, so it sees two tuples.

On the second record the two readings also name different terms: an in-frame insertion of 39 codons, against a coding-sequence change that lengthens the feature.

The `OverlapBP` and `OverlapPC` fields of the symbolic rows are described under [FLAGS, OverlapBP and OverlapPC](../output/vep-and-tab.md#flags-overlapbp-and-overlappc).

### 12.4 Why Ensembl VEP writes that

Ensembl VEP builds the alternate allele as a literal sequence (it materialises the allele). It reads the reference bases of the run and builds the alternate as `RUS` repeated `RUC` times (or `RB` divided by the unit's length, in whole units, when `RUC` is absent). It then hands the pair to its sequence-variant machinery, which trims the common prefix. What is left is a gain as an insertion of whole units after the run's last base, or a loss as a deletion of whole units ending at the run's last base.

The materialisation is skipped for the whole record when any repeat sequence of any of its alleles is longer than `--max_sv_size` bases (the same 10,000,000 default; under `-1` it is always skipped), and every allele of the record is then annotated symbolically. vep-rs expands under the same limit, so above it both engines give the symbolic gain.

The materialisation also needs the run's reference bases. A run without `--fasta` (`--offline` with no FASTA, the protocol of the golden corpora; the measured structural-variant sets ran with the FASTA, so their reference rows are the expanded alleles) fetches no slice and builds no literal allele. Every `<CNV:TR>` record is then annotated as a symbolic `tandem_repeat` over the run, which `copy_number_gain` reads as a gain, so a loss is labelled `feature_elongation` with the region terms (lines 516-527 of the same file; `Utils/VariationEffect.pm` lines 246-256).

vep-rs reads the record the same way. With a FASTA it expands the allele from the run's bases, the record then a sequence record; without one it keeps the symbolic `tandem_repeat` over the run, a gain whatever `RUC` says, and writes the same warning once per run.

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`, lines 382-431: builds the literal alternate allele; called from `create_StructuralVariationFeatures` at line 520 once the reference sequence of the run has been fetched.

### 12.5 Why the symbolic reading stands, and why Ensembl VEP's reading is not a defect

Ensembl VEP's reading is a faithful annotation of a real molecule. A run of 121 CA units where the reference has nine is, base for base, an insertion of 224 bases after the run, and an in-frame gain of 39 GTG units inside a coding exon does insert 39 codons. Nothing in its rows is false, and this class is not a defect.

When `RUC` is exact the literal reading also carries what the symbolic one does not: the reading frame and the inserted residues, and an IMPACT of MODERATE for an in-frame gain where the symbolic reading's `feature_elongation` is HIGH.

The symbolic reading keeps the record in the terms the VCF author chose: the variant is a copy-number change of a repeat, its location is the repeat, its allele is the repeat class, and its consequence is the region the run lies in. Where the run touches exonic sequence, the consequence adds `feature_elongation`, the term `copy_number_gain` gives a tandem repeat.

That reading is stable under the information the record carries and does not carry:

- it does not depend on the reference sequence of the run being fetched;
- it does not turn a copy-number call into a literal sequence the caller never asserted (a `<CNV:TR>` record may carry confidence intervals on the unit count, `CIRUC` and `CIRB`, which Ensembl VEP's expansion ignores with a warning);
- it is the reading both engines give the record when neither reads a FASTA, and Ensembl VEP gives it with one when it does not expand the allele ([12.6](#126-try-it-yourself)).

Neither reading is wrong; the comparator pairs them rather than ranking them.

### 12.6 Try it yourself

Ensembl VEP 116.2, GRCh38 cache and FASTA, [the command of the closing section](reproducing.md#reproducing-these-examples). Lower the expansion cap below the alternate allele's length and Ensembl VEP annotates the record symbolically. Two settings do that for the two records:

- `--max_sv_size 200` on the first record (alternate 242 bases);
- `--max_sv_size 100` on the second record (alternate 138 bases).

With `--max_sv_size 200` the Ensembl VEP row for ENST00000346798 on the first record becomes

```
synth_cnv_tr_0001	21:25896419-25896436	tandem_repeat	ENSG00000142192	ENST00000346798	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;OverlapBP=18;OverlapPC=0.01
```

With `--max_sv_size 100` the row on the second record becomes

```
synth_cnv_tr_0023	21:26021866-26021886	tandem_repeat	ENSG00000142192	ENST00000346798	Transcript	feature_elongation,coding_sequence_variant	969-989	819-839	273-280	-	-	-	IMPACT=HIGH;STRAND=-1;OverlapBP=21;OverlapPC=0.01
```

Both rows are byte for byte the vep-rs rows of [12.3](#123-what-each-engine-writes). Running Ensembl VEP without `--fasta` gives the same two symbolic rows, with the warning `could not fetch sequence for tandem repeats; consequence calculation for tandem repeats is less precise when using --offline without --fasta`. At the default `--max_sv_size` of 10,000,000 and with the FASTA every `<CNV:TR>` record of the synthetic sets is expanded.

### 12.7 How the comparator treats this class

`filter_cnv_tr_expansion_divergences` in `scripts/validation/compare_sv_concordance.py` pairs a vep-rs tuple whose Allele is `tandem_repeat` with the Ensembl VEP tuple on the same chromosome, Feature and Feature_type that is its literal reading (`_cnv_tr_pair_kind`):

- for a gain, an insertion (Allele all bases) located between the run's last base and the next;
- for a loss, a deletion (Allele `-`) ending at the run's last base;
- for a literal reading the expansion does not trim to a unit indel, any other sequence allele (`N` bases included) whose span touches the run or the base after it. Those readings are a multi-allelic `<CNV:TR>` record, whose allele string is never prefix-trimmed, a unit that is not the reference's repeat, a `.` unit read as `N`, or an `SVLEN` that disagrees with the unit count. Each is a substitution across the run (`AGCAGCAGCAGC protein_altering_variant`, `AGCAGC inframe_deletion`, `NNNN frameshift_variant`);
- for an unexpanded reading, the reference's own `tandem_repeat` tuple at the run, written without a FASTA ([12.4](#124-why-ensembl-vep-writes-that)).

Both members of the pair are set aside when the two consequence sets agree after the symbolic reading's direction term is disregarded. That term, where present, is `feature_elongation` on the symbolic reading, a gain whatever the count, and `feature_elongation` or `feature_truncation` on a literal reading.

Inside a coding sequence, the literal reading's coding effect is also treated as matching the symbolic reading's `coding_sequence_variant`. That effect is `inframe_insertion`, `inframe_deletion` or `frameshift_variant` for a unit indel, and any coding term, `protein_altering_variant` or `missense_variant` included, for a literal substitution. An unexpanded pair may differ in the direction and coding terms alone (`_is_cnv_tr_representation_pair`). Any other difference, a start or splice term one side carries and the other lacks, stays charged.

The pairing reads every shape the expansion can take, the ones below included. The run at 30425588-30425596 is three copies of `AGC`, the run at 30425756-30425763 four copies of `AG`, the run at 30430410-30430417 four copies of `CT`, each with `--fasta` unless stated:

- a unit that is not DNA, which Ensembl VEP repeats verbatim across the run and vep-rs keeps symbolic. `RUS=123;RUC=4` gives `21:30425588-30425596 123123123123 ENSG00000240432 ENST00000390690 Transcript coding_sequence_variant`. `RUS=R;RUC=10` gives `RRRRRRRRRR` with `frameshift_variant`, `RUS=*;RUC=9` gives `*********` with `coding_sequence_variant`, and `RUS=X;RUC=10` gives `XXXXXXXXXX` with `frameshift_variant`. Each pairs with vep-rs's `tandem_repeat feature_elongation,coding_sequence_variant` at the run;
- a loss whose literal reading is a substitution with a coding term other than a deletion's. `RN=3;RUS=AG,AG,C;RUC=1,1,1` gives `21:30425760-30425763 C ENSG00000240432 ENST00000390690 Transcript protein_altering_variant`. `RUS=.;RB=5` gives `21:30425756-30425763 NNNNN ENSG00000240432 ENST00000390690 Transcript coding_sequence_variant`. `SVLEN=30_425_949;RUS=AGC;RUC=4` (an `SVLEN` Perl reads as 30) gives `21:30425597-30425617 AGC ENSG00000240432 ENST00000390690 Transcript protein_altering_variant`. `SVLEN=4000;RUS=CT;RUC=5` on the run at 30430410 gives `21:30430418-30434408 C ENSG00000186971 ENST00000334068 Transcript stop_retained_variant,3_prime_UTR_variant` against the symbolic reading `21:30430410-30434409 tandem_repeat` over the whole run, the literal reading's stop term against the span's;
- a multi-allelic record whose reference rows collapse to one tuple, the reference run without a FASTA. `<CNV:TR>,<CNV:TR>,<CNV:TR>` with `RN=1,1,1;RUS=AGC,AGC,AGC;RUC=4,2,3` is one tuple in the reference, `21:30425588-30425596 tandem_repeat ENSG00000240432 ENST00000390690 Transcript feature_elongation,coding_sequence_variant` for all three alleles (every allele a gain without a FASTA, [12.4](#124-why-ensembl-vep-writes-that)). In vep-rs with a FASTA it is a gain tuple and a loss tuple (`feature_truncation,inframe_deletion`), both paired with it; without one it is the same one tuple.

An expansion the reference drops without a row is not a pair but [class 19](defect-l.md#19-an-empty-or-unchanged-cnvtr-expansion-dropped-without-a-row) ([19.7](defect-l.md#197-how-the-comparator-treats-this-class)): an unchanged expansion over several repeat sequences (`RN=3;RUS=AGC,AGC,AGC;RUC=1,1,1`, `RN=2;RUS=AGC,AG;RUC=3,0`), a unit of no bases under `RUC` (`RUS=,AGC;RUC=4,4`) and a count in fullwidth digits (`RUC=４`, which Perl numifies to 0).

The class is fully masked:

- GRCh37 set: 1,368 of 1,368 pairs;
- GRCh38 set: 3,751 of 3,751 pairs.

The counters are `excluded_cnvtr_rust` and `excluded_cnvtr_perl`, against the class sizes `cnvtr_swap_pairs_total_rust` and `cnvtr_swap_pairs_total_perl`.
