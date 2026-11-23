# The order of the pairs in UTRAnnotator's `5UTR_annotation`

**The defect in one sentence.** Ensembl VEP's UTRAnnotator plugin keeps an effect's details in hashes and writes `5UTR_annotation` by walking their keys, so the `key=value` pairs of a group, and the groups of a value, come in Perl's hash order, which changes from one run to the next. vep-rs writes the pairs Ensembl VEP writes, each group's in byte order of their keys, so every run of an input writes the same value.

**Why this is a defect.** The variant, the transcript and the uORF table fix the pairs: the distances, the Kozak context and strength, the ORF type and the evidence flag of each upstream ORF the variant touches. Their order in Ensembl VEP's output comes from none of them. The same input, flags and cache write a row's pairs in one order on one run and in another on the next, so two runs of one command differ byte for byte on almost every row that carries the field.

**What vep-rs does instead.** vep-rs writes each group's pairs in byte order of their keys, capitalised keys first, and the groups in the order its analysis finds them, so every run of an input writes the same value. The pairs themselves are the ones Ensembl VEP writes. The plugin's parameters, fields and analyses are described under [UTRAnnotator](../plugins/transcript-model.md#utrannotator).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `UTRAnnotator.pm`, lines 828-840, and the same pattern in each of the five analyses: an affected ORF's details are a hash keyed by name, filed in the analysis's hash `%output_effects` under a running number. Lines 250-267 keep the reported analysis's hash as `$output_five_prime_annotation`.
- Line 282: `foreach my $value (keys %{$output_five_prime_annotation})` walks that hash's groups.
- Line 283: inside each group, `keys %{$output_five_prime_annotation->{$value}}` walks the pairs, which are joined with `:`; line 288 joins the groups with `,`.

**What correcting it would change.** Sorting the keys at both walks would give one text for one variant on every run.

**The class it produces.** `utrannotator_annotation_pair_order`: a `5UTR_annotation` value whose groups, in any order, carry vep-rs's pairs group for group, each group's pairs in an order Perl's hash order decided for that run. vep-rs's pairs stand in byte order of their keys. A value with a pair changed, missing or added is not the class and stays a difference. The JSON format writes the field as an object that holds each group under a number, each group an object of its pairs; an object has no order, so the class occurs in the default, tab and VCF formats alone.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md). The trigger is a value of two or more pairs, which the transcript decides and the record alone does not, so the rate is counted from the output rows of `--plugin UTRAnnotator,file=<the assembly's uORF table>` over the cuts. 8,984 of the cuts' 15,646,481 records have a row that carries the field, 574.2 per million, on 41,583 rows, and every such value holds six or more pairs.

Four committed flag sets of the golden corpora run the plugin, two in each UTRAnnotator corpus, and the committed reference outputs of those two corpora hold such values, 87 and 40 under their default runs.

On the chromosome 21 records of the ClinVar GRCh38 dataset, Ensembl VEP 116.3 writes the field on 598 rows of 161 records. Over ten runs every one of them changed its text between two runs: 572 in the order of their pairs alone, and 26 because the plugin also reported another effect ([its own page](utrannotator-consequence-precedence.md)). On the 63 rows whose value has two or more groups in every run, the groups changed their order as well.

vep-rs writes the field on 596 of those rows, on 160 records. The other two are the rows of one deletion, on which Ensembl VEP alone reports a gained start ([the window of the gained-start test](../plugins/transcript-model.md#1-utrannotator-the-window-of-the-gained-start-test)).

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A substitution that destroys an upstream start codon

**In short.** ClinVar variation 4159501 is a substitution at `21:14516984`, cDNA position 410 of `ENST00000619120`, the last base of an upstream ORF's start codon. Both engines report that one lost start codon with the same seven pairs. Ensembl VEP 116.3 wrote the pairs in ten different orders over ten runs of one command; vep-rs writes them in byte order of their keys on every run.

### 1.1 What a reader needs to know first

An upstream ORF is an open reading frame that starts in a transcript's 5' UTR. UTRAnnotator runs five analyses over a 5' UTR variant: a start codon gained or lost, an upstream ORF's stop codon gained or lost, and a frameshift in an upstream ORF. It reports one effect in `5UTR_consequence` and that effect's details in `5UTR_annotation`.

The details of each affected ORF form a group of `key=value` pairs joined by `:`, and a value of several groups joins them by `,`. For a lost start codon the keys are `CapDistanceToStart` (the bases before the codon), `DistanceToCDS`, `DistanceToStop`, `KozakContext` (the bases around the codon), `KozakStrength` (their rating), `type` and `Evidence` (whether the uORF table given with `file=` lists the ORF).

A Perl hash holds its keys in no order of its own. `keys` returns them in an order set by the hash function and its seed, and Perl draws a new seed for each process unless the environment variable `PERL_HASH_SEED` fixes one. Perl also perturbs the order hash by hash as keys are stored, unless `PERL_PERTURB_KEYS=0` switches that off, so the order is fixed only with both variables set.

### 1.2 The record

ClinVar variation 4159501 on GRCh38, a substitution in the 5' UTR of a reverse-strand transcript on chromosome 21, against the Ensembl 116 GRCh38 cache and the published GRCh38 uORF table:

```
21	14516984	4159501	C	G	.	.	.
```

Both engines run `--plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt` with `--fasta`, in the default format. The coding sequence of `ENST00000619120` starts at cDNA position 431, so its 5' UTR is 430 bases long.

### 1.3 What each engine writes

The row on `ENST00000619120`. Ensembl VEP 116.3, run with `PERL_HASH_SEED=0` and `PERL_PERTURB_KEYS=0`:

```
4159501	21:14516984	G	ENSG00000155307	ENST00000619120	Transcript	5_prime_UTR_variant	410	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;5UTR_annotation=KozakContext=ATAATGG:type=OutOfFrame_oORF:KozakStrength=Strong:DistanceToStop=24:DistanceToCDS=23:CapDistanceToStart=407:Evidence=False;5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=4;Existing_uORFs=3
```

Ensembl VEP 116.3 with the two variables unset, one run of ten:

```
4159501	21:14516984	G	ENSG00000155307	ENST00000619120	Transcript	5_prime_UTR_variant	410	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;5UTR_annotation=KozakStrength=Strong:Evidence=False:KozakContext=ATAATGG:type=OutOfFrame_oORF:CapDistanceToStart=407:DistanceToStop=24:DistanceToCDS=23;5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=4;Existing_uORFs=3
```

vep-rs, in every run:

```
4159501	21:14516984	G	ENSG00000155307	ENST00000619120	Transcript	5_prime_UTR_variant	410	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;5UTR_annotation=CapDistanceToStart=407:DistanceToCDS=23:DistanceToStop=24:Evidence=False:KozakContext=ATAATGG:KozakStrength=Strong:type=OutOfFrame_oORF;5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=4;Existing_uORFs=3
```

The three rows carry the same seven pairs, one group, and agree in every other column. Under the fixed seed the keys read `KozakContext`, `type`, `KozakStrength`, `DistanceToStop`, `DistanceToCDS`, `CapDistanceToStart`, `Evidence`. In the unseeded run they read `DistanceToStop`, `CapDistanceToStart`, `KozakContext`, `KozakStrength`, `DistanceToCDS`, `Evidence`, `type`, and each of the nine other unseeded runs wrote an order of its own.

vep-rs writes `CapDistanceToStart`, `DistanceToCDS`, `DistanceToStop`, `Evidence`, `KozakContext`, `KozakStrength`, `type`.

### 1.4 Why Ensembl VEP writes that

The substitution changes the last base of a start codon at cDNA positions 408 to 410, so the plugin's lost-start analysis fires and files the ORF's details in a hash keyed by name. When the row is written, the plugin walks `keys` of the reported analysis's hash to find its groups. Inside each group it walks `keys` of the group's hash and prints one `key=value` pair per key, in the order `keys` returns them.

That order is the hash's own. It depends on where the seeded hash function places each key and on the perturbation Perl applies as keys are stored. The seed is drawn per process, so the order moves from run to run. Nothing sorts the keys.

Source (`VEP_plugins` at release/116):

- `UTRAnnotator.pm`, lines 828-840: the lost-start analysis files the seven pairs in a hash.
- Line 282: the walk over the groups.
- Line 283: the walk over a group's pairs.

### 1.5 Why vep-rs's row is the correct one

The pairs follow from the record, the transcript and the uORF table. The start codon is cDNA positions 408 to 410, so 407 bases precede it (`CapDistanceToStart=407`) and the coding sequence at 431 starts 23 bases after its first base (`DistanceToCDS=23`). Its stop codon starts 24 bases after its first base (`DistanceToStop=24`), one base past the coding sequence's first base and out of its frame (`type=OutOfFrame_oORF`). The bases around the codon rate `Strong`, and the uORF table does not list the ORF (`Evidence=False`). Both engines print those pairs and the context's own bases, seven pairs and no other.

Only the order differs, and Ensembl VEP's own runs of this one input disagree with each other about it, so none of its orders can be read from the variant or the transcript. vep-rs's order is a rule of the text itself, the byte order of the keys. It is the same on every run, so two runs of one input can be compared byte for byte.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache <corpus>/json_cache --fasta <corpus>/reference.fa.gz \
  --species homo_sapiens --assembly GRCh38 \
  --plugin UTRAnnotator,file=<corpus>/data/uORF_5UTR_GRCh38_PUBLIC.chr21.txt
```

`<corpus>` is `tests/golden/116/GRCh38-plugin-utr-annotator`: its cache holds the transcript, and its table is the chromosome 21 rows of the published one. vep-rs writes the seven pairs in key order on the row of `ENST00000619120`. Ensembl VEP 116.3 takes the same plugin flag with the cache options of [Reproducing these examples](reproducing.md). With `PERL_HASH_SEED=0` and `PERL_PERTURB_KEYS=0` in its environment it writes the first row of 1.3 on every run; without them, repeated runs write other orders.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) applies `utrannotator_annotation_pair_order` when a comparison case, one run of both engines over one input, names it. A `5UTR_annotation` is set aside, and counted once per row under the class, when its two values carry the same groups of the same pairs and differ in the order of the groups or of the pairs inside a group. A value with a pair changed, missing or added, a value one engine alone writes, and the row's other fields stay differences.

The class also asks that vep-rs's pairs stand in byte order of their keys inside every group. It documents Ensembl VEP's changing order against vep-rs's fixed one, so a vep-rs value out of key order is a difference of vep-rs's own. It is reported, never set aside, even when its pairs are Ensembl VEP's.

UTRAnnotator has eight classes, this one, `utrannotator_gained_start_window`, `utrannotator_consequence_precedence`, `utrannotator_structural_allele_text`, `utrannotator_lost_start_after_insertion`, `utrannotator_gained_start_remade_in_place`, `utrannotator_frameshift_at_stop_boundary` and `utrannotator_lost_stop_written_before_stop_codon`, and a differing field is counted under the first that accepts it. This class comes first: the window class reads every UTRAnnotator field of an insertion's or a deletion's row, and would otherwise take such a row's reordered value. A value whose pairs differ because another effect was reported belongs to the precedence class.

A case may name the comparator's built-in class `order_only` instead. That class reads `5UTR_annotation` as groups of pairs in any order and counts the same rows under its own name; the golden harness compares the two UTRAnnotator corpora under it. Where a case names both, `order_only` takes the row.
