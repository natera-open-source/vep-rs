# The effect UTRAnnotator reports when several analyses fire

**The defect in one sentence.** When several of UTRAnnotator's five analyses fire on one row, Ensembl VEP's plugin keeps the one whose key a hash lists last, so `5UTR_consequence` and `5UTR_annotation` change between runs of the same input. vep-rs reports one effect under a fixed precedence, with the chosen effect's `5UTR_annotation` group, so every run writes the same row.

**Why this is a defect.** A substitution inside an upstream ORF can at once destroy its start codon, create a new start codon and lose or gain a stop codon, and the plugin computes all five effects. It then walks the five flags as a hash and overwrites the output with every truthy one, so the row reports whichever effect the hash happened to list last. Perl randomises hash order per process, and per hash unless `PERL_PERTURB_KEYS=0`, so the effect a user sees is not a property of the variant.

**What vep-rs does instead.** vep-rs computes the same five effects and reports one under a fixed precedence: a gained start codon, then a gained stop codon, then a lost stop codon, then a lost start codon, then a frameshift. The `5UTR_annotation` group is the chosen effect's. The plugin's parameters, fields and analyses are described under [UTRAnnotator](../plugins/transcript-model.md#utrannotator).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `UTRAnnotator.pm`, lines 234-238: the five analyses run and each returns its flag and effect; lines 242-256: the flags and effects are keyed by analysis name in two hashes.
- Lines 262-267: `foreach my $flag (keys %five_prime_flag)` assigns `$output_five_prime_flag` and `$output_five_prime_annotation` for every truthy flag, so the last key in hash order wins; lines 282-284 print that one effect's pairs, again in hash order.

**What correcting it would change.** The plugin would report every firing effect, or one chosen by a stated rule, and two runs over one input would agree.

**The class it produces.** `utrannotator_consequence_precedence`: a `5UTR_consequence` whose two values are both among the five `5_prime_UTR_` terms and differ, on a row where Ensembl VEP's own runs of the input under other Perl hash seeds write both terms. It is set aside with the row's `5UTR_annotation` where each engine's annotation is one those runs write beside its term, and counted once per row. A value of `-` on either side stays a difference, as does a row the runs do not hold. The committed UTRAnnotator corpora carry no such row. The concordance figures score consequence sets and never see the field, so nothing is set aside in them.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape the record alone shows is a sequence substitution, insertion or deletion, 15,645,693 of the cuts' 15,646,481 records (999,950 per million), in 109 corpora. The transcript decides whether the variant lies in a 5' UTR and whether more than one analysis fires there. Four committed flag sets of the golden corpora run the plugin, two in each UTRAnnotator corpus, and the UTRAnnotator comparison cases run it too.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A substitution that both destroys an upstream start codon and removes an upstream stop codon

**In short.** The record `21 35748195 rs1396983766 G T` on `ENST00000399295`, cDNA position 293, lies in a uORF whose `ATG` and a stop codon it alters at once. Ensembl VEP 116.2 writes the start-codon loss in one run and the stop-codon loss in another; vep-rs writes the stop-codon loss in every run. On the chromosome 21 gnomAD v2.1.1 GRCh37 dataset, 36 rows change their `5UTR_consequence` across six runs of Ensembl VEP 116.2 (seeds 0, 42 and 12345, and three unseeded), and vep-rs's value is one of the outcomes on every one.

### 1.1 What a reader needs to know first

UTRAnnotator runs five analyses over a 5' UTR variant: a start codon gained (`uAUG_gained`) or lost (`uAUG_lost`), a uORF's stop gained (`uSTOP_gained`) or lost (`uSTOP_lost`), and a frameshift in a uORF (`uFrameshift`). The row carries one flag in `5UTR_consequence` and that flag's effects in `5UTR_annotation`. When two analyses fire, one of them has to be chosen, and the choice is the subject of this page.

### 1.2 The record

A gnomAD v2.1.1 substitution in the 5' UTR of a forward-strand transcript on chromosome 21, against the Ensembl 116 GRCh37 cache and the published GRCh37 uORF table:

```
21	35748195	rs1396983766	G	T	.	.	.
```

Both engines run `--plugin UTRAnnotator,file=uORF_5UTR_GRCh37_PUBLIC.txt` with `--fasta`.

### 1.3 What each engine writes

The row on `ENST00000399295` in the default format. Ensembl VEP, unseeded:

```
rs1396983766	21:35748195	T	ENSG00000205670	ENST00000399295	Transcript	splice_region_variant,5_prime_UTR_variant	293	-	-	-	-	-	IMPACT=LOW;STRAND=1;5UTR_annotation=CapDistanceToStart=290:KozakStrength=Moderate:type=uORF:Evidence=False:KozakContext=GGAATGA:DistanceToCDS=80:DistanceToStop=63;5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=4
```

Ensembl VEP under `PERL_HASH_SEED=0` with `PERL_PERTURB_KEYS=0`, and in another unseeded run, reports the other effect:

```
rs1396983766	21:35748195	T	ENSG00000205670	ENST00000399295	Transcript	splice_region_variant,5_prime_UTR_variant	293	-	-	-	-	-	IMPACT=LOW;STRAND=1;5UTR_annotation=KozakContext=GGGATGC:KozakStrength=Moderate:FrameWithCDS=outOfFrame:AltStopDistanceToCDS=NA:AltStop=False:Evidence=False;5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=4
```

vep-rs, in every run:

```
rs1396983766	21:35748195	T	ENSG00000205670	ENST00000399295	Transcript	splice_region_variant,5_prime_UTR_variant	293	-	-	-	-	-	IMPACT=LOW;STRAND=1;5UTR_annotation=AltStop=False:AltStopDistanceToCDS=NA:Evidence=False:FrameWithCDS=outOfFrame:KozakContext=GGGATGC:KozakStrength=Moderate;5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=4
```

The consequence terms, the `Existing_` counts and every other column agree. The order of the `5UTR_annotation` pairs is Ensembl VEP's hash order, which can change from one run to the next and which the comparison tools read as a set.

### 1.4 Why Ensembl VEP writes that

The substitution lands on the third base of the `ATG` at cDNA position 291 of a uORF (`CapDistanceToStart=290`) and on the stop codon of an overlapping ORF, so `uAUG_lost` and `uSTOP_lost` both fire. `UTRAnnotator.pm` keys the five flags by analysis name in a hash and, at lines 262-267, assigns the output flag and annotation for every truthy key in turn, so the key the hash lists last wins.

Perl orders a hash's keys by a per-process seed, and perturbs the order per hash unless `PERL_PERTURB_KEYS=0`. Under seeds 0 to 7 Perl 5.34 lists the five keys `uAUG_lost`, `uSTOP_lost`, `uAUG_gained`, `uFrameshift`, `uSTOP_gained`, and under seeds 42 and 12345 `uFrameshift`, `uSTOP_gained`, `uAUG_lost`, `uSTOP_lost`, `uAUG_gained`; in both orders `uSTOP_lost` follows `uAUG_lost`, so the stop-codon loss is written. An unseeded process draws its own order, and the one behind the first row of 1.3 listed `uAUG_lost` after `uSTOP_lost`.

### 1.5 Why vep-rs's row is the correct one

Both effects are facts of the record: the uORF's start codon is destroyed and the overlapping ORF's stop codon is removed. The plugin's output carries one effect per row, so an engine must choose, and the choice has to be a rule of the variant rather than of the process. Ensembl VEP's two rows are each a true effect; which one appears is drawn at random.

vep-rs's rule is the fixed precedence stated under "What vep-rs does instead" at the top of this page: a gained start codon, then a gained stop codon, then a lost stop codon, then a lost start codon, then a frameshift. A lost stop codon stands before a lost start codon in it, so the stop-codon loss is reported here, and every run writes the same row with the same `5UTR_annotation` group. The effect it reports is one Ensembl VEP's own runs report as well, in the second row of 1.3.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh37 --plugin UTRAnnotator,file=uORF_5UTR_GRCh37_PUBLIC.txt
```

vep-rs writes `5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant` on `ENST00000399295` in every run. Ensembl VEP 116.2 writes it under `PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0` and `5_prime_UTR_premature_start_codon_loss_variant` in some unseeded runs; running it several times over the one record shows both.

### 1.7 How the comparator treats this class

The format comparator (`scripts/concordance/compare_vep_formats.py`) applies `utrannotator_consequence_precedence` when a comparison case names it and hands it Ensembl VEP's outputs of the same input under fixed Perl hash seeds (`--utrannotator-sweep`, one default-format output per seed). Each of the plugin's five flags is walked last under one of the seeds, so the terms those runs write on a row are the analyses that fire on it. A differing `5UTR_consequence` is set aside when both values are terms the runs write on the row. The row's `5UTR_annotation` is set aside with it when each side's value is one the runs write beside its term. The row is counted once.

Without the runs the class sets nothing aside, and a VCF or JSON entry, which names no Location, matches no row of them. A `-` on either side, a value outside the five terms, a term the runs never write on the row and an annotation they do not write beside its term stay differences. A row the class `utrannotator_gained_start_window` takes ([UTRAnnotator: the window of the gained-start test](../plugins/transcript-model.md#1-utrannotator-the-window-of-the-gained-start-test)) is counted under that class alone when a case names both. The `GRCh38-plugin-utr-annotator` and `GRCh37-plugin-utr-annotator` golden corpora carry no row of the shape, so the class applies to none of their rows.
