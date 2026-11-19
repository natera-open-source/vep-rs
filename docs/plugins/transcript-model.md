# Transcript-model plugins: Blosum62, Downstream, HGVSIntronOffset, NearestExonJB, TSSDistance, UTRAnnotator, MaxEntScan

Five Ensembl VEP plugins compute their fields from the transcript model and the
variant alone, with no data file. vep-rs provides each as a built-in with the
plugin's name and the parameters stated below, its output matching Ensembl
VEP's exactly except for the documented divergences. The five share
one corpus comparing every row with Ensembl VEP's in all four output formats.
Two more, UTRAnnotator and MaxEntScan, read the transcript model and the
reference sequence, each with a data directory or table of its own.

All seven write on transcript rows alone. Five follow Ensembl VEP's default and
leave a structural variant's rows empty; TSSDistance and UTRAnnotator ask for
structural variants, as their Ensembl VEP plugins do.
TSSDistance annotates them; UTRAnnotator annotates some
([its section](#utrannotator)).

## Blosum62 {#blosum62}

| Property      | Value                                                        |
| ------------- | ------------------------------------------------------------ |
| Builds        | GRCh37, GRCh38                                               |
| Output fields | `BLOSUM62`                                                   |
| Parameters    | none; a positional one is not read, a `key=value` one refused |

```bash
--plugin Blosum62
```

The value is the BLOSUM62 substitution score of a row whose amino-acid change
is one capital letter over another. It is read from the matrix the plugin
embeds, at the row of the reference residue and the column of the alternate
one. A
synonymous row (one letter), a stop (`*`), an unknown residue (`X`) or a letter
outside the twenty residues of the matrix gets no value.
Under `--hgvs`, Ensembl VEP alone can score a later allele's inserted residue
([a documented class](../divergences/multiallelic-hgvs-shift.md#3-a-record-whose-first-allele-gives-the-second-a-blosum62-score)).

## Downstream {#downstream}

| Property      | Value                                                        |
| ------------- | ------------------------------------------------------------ |
| Builds        | GRCh37, GRCh38                                               |
| Output fields | `DownstreamProtein`, `ProteinLengthChange`                   |
| Parameters    | none; a positional one is not read, a `key=value` one refused |

```bash
--plugin Downstream
```

The plugin writes on a row with `frameshift_variant` and no splice-site term.
The allele's bases on the transcript strand replace its CDS span in the
translateable sequence, and the 3' UTR is appended. The translation starts at
the first codon the change touches, under the transcript's codon table.
Ensembl VEP starts one codon later for a deletion from a codon's third base,
and out of frame when it ends the coding sequence.

`DownstreamProtein` is that
translation up to its first stop, empty when the first codon is a stop.
`ProteinLengthChange` is the whole changed translation's length up to its
first stop minus the reference protein's, negative for a truncation.

Ensembl VEP's plugin needs a reference FASTA offline to spell the 3' UTR.
vep-rs reads the UTR from the cache, or from `--fasta` when the cache carries
none; with neither, the coding sequence translates alone.

vep-rs reads the record's own CDS span, the pair the `CDS_position` column
prints, and the allele as the row prints it, with and without `--hgvs`. Both
runs therefore write the same two fields. Under `--hgvs` Ensembl VEP's plugin reads
the span at the 3'-shifted position its protein notation left on the row, with
the allele's bases rotated for that position. It can write another value:
[the Downstream and NMD plugins at the HGVS-shifted coordinates](../divergences/plugin-hgvs-shifted-cds-coordinates.md).
An insertion whose point lies on a seam of an edited RefSeq model, where the
model carries bases the genome lacks, replaces those bases in both engines.

## HGVSIntronOffset {#hgvsintronoffset}

| Property      | Value                                                        |
| ------------- | ------------------------------------------------------------ |
| Builds        | GRCh37, GRCh38                                               |
| Output fields | `HGVS_IntronStartOffset`, `HGVS_IntronEndOffset`             |
| Parameters    | none; a positional one is not read, a `key=value` one refused |

```bash
--hgvs --plugin HGVSIntronOffset
```

The two values are the intron offsets of the HGVSc start and end coordinates:
`45` for `c.123+45`, `-10` for `c.124-10`, `0` for an exonic coordinate. They
are read before the coordinates are ordered for printing. A row without an HGVSc
gets no value. The plugin reads the coordinates `--hgvs` computes, so a run
without the flag is refused at load. Ensembl VEP's plugin computes the
notation itself and writes offsets whether or not HGVSc is printed.

## NearestExonJB {#nearestexonjb}

| Property      | Value                                                        |
| ------------- | ------------------------------------------------------------ |
| Builds        | GRCh37, GRCh38                                               |
| Output fields | `NearestExonJB`                                              |
| Parameters    | `max_range=<bases>` (default 10000), `intronic=0` or `1`; each is read as Ensembl VEP's numeric comparison reads it: the leading number, else 0, so `max_range=len` reports no exon where Ensembl VEP reports every one ([a documented divergence](../divergences/nearest-exon-jb-len-key-collision.md)). `parent_only=0` or `1`, and a `gff3` value Perl reads as false (`0` or empty), change nothing, as in Ensembl VEP |

```bash
--plugin NearestExonJB
--plugin NearestExonJB,max_range=50000,intronic=1
```

Each item is `ExonID|distance|start/end|length`:

| Part        | What it holds                                                                                          |
| ----------- | ------------------------------------------------------------------------------------------------------ |
| `ExonID`    | the exon's stable id, empty for an exon without one (a RefSeq aligned transcript's under `--all_refseq`) |
| `distance`  | the distance to its nearer junction                                                                    |
| `start/end` | that junction (`start`, `end`, `start_end` on a tie)                                                   |
| `length`    | its length                                                                                             |

Items are joined by `,`, their parts by
`+` under `--vcf`. Two id-less exons of a row keep their own items; Ensembl VEP
folds them into one ([a documented divergence](../divergences/nearest-exon-jb-id-less-exon-fold.md)).

For every exon the variant overlaps, its start is measured to the exon's start
and its end to the exon's end. The shorter is the exon's distance. The
exons at the smallest are reported when it lies within `max_range`. An exon is
stretched by 12 bases on a transcript with an intron of at most 13 bases.

vep-rs applies the stretch on every row. Ensembl VEP 116.2 applies it once an
earlier row has built the transcript's intron index. A row within 12 bases
outside the first or last exon therefore answers by the records before it ([section 3 below](#3-nearestexonjb-a-value-on-an-upstream-or-downstream-row-that-depends-on-the-records-before-it), tier 2).

Under `intronic=1` a row with no overlapped exon whose terms include
`intron_variant` reports the two exons flanking its intron, each at its own
distance, without `max_range`. A variant overlapping exactly two exons reports
both likewise. The header line carries `(coding sequence variants only)`
only without the parameter.

Under `--minimal`, a split allele reads its own trimmed span for its exon set
and distances; Ensembl VEP's rejoin hands the plugin the record's span
([a documented class](../divergences/minimal-plugin-sees-the-rejoined-record.md)).
Each allele answers for its own span, and under `intronic=1` its own terms.
Ensembl VEP's cache hands one row's value to the record's other rows
([a documented divergence](../divergences/nearest-exon-jb-shared-span-cache.md)).
Under `--shift_3prime 1` the exon set is the record's span's, like the
distances; Ensembl VEP reads the shifted span's exons
([a documented divergence](../divergences/nearest-exon-jb-shifted-span.md)).

Ensembl VEP lists several items in [Perl hash order](../divergences/nearest-exon-jb-item-order.md), which changes between runs.
vep-rs lists them in transcript order, and the corpus compares them as sets.
`gff3=<file>` is refused: vep-rs reads the cache's transcript model alone.

## TSSDistance {#tssdistance}

| Property      | Value                                                        |
| ------------- | ------------------------------------------------------------ |
| Builds        | GRCh37, GRCh38                                               |
| Output fields | `TSSDistance`                                                |
| Parameters    | `both_direction=1`; a token without `=` is not read, and is noted in the warnings file |

```bash
--plugin TSSDistance
--plugin TSSDistance,both_direction=1
```

On a forward-strand transcript the distance is the transcript's start minus the
variant's end, on a reverse-strand one the variant's start minus the
transcript's end. It is positive upstream of the transcription start site. An
upstream row carries it; under `both_direction=1` every transcript row does,
zero at the start site and negative inside or downstream of the transcript.
Under `--minimal`, a split allele measures from its own trimmed span, the one its
row's alleles describe; Ensembl VEP's rejoin hands the plugin the record's span
([a documented class](../divergences/minimal-plugin-sees-the-rejoined-record.md)).

The parameter is read by its truth, as Ensembl VEP reads it: `both_direction=0`
and an empty value are off, any other value on. One without `=` leaves Ensembl
VEP's parameter hash empty, the flag off, in vep-rs too.

## UTRAnnotator {#utrannotator}

UTRAnnotator reports the effect of a 5' UTR variant on upstream open reading
frames, as Ensembl VEP's `UTRAnnotator.pm` computes it.

| Property      | Value                                                                                                                                     |
| ------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| Builds        | GRCh37, GRCh38                                                                                                                            |
| Output fields | `5UTR_consequence`, `5UTR_annotation`, `Existing_uORFs`, `Existing_OutOfFrame_oORFs`, `Existing_InFrame_oORFs`, on rows carrying `5_prime_UTR_variant` |
| Data source   | `uORF_5UTR_GRCh37_PUBLIC.txt` and `uORF_5UTR_GRCh38_PUBLIC.txt` of the Ensembl UTRannotator repository, the uORF evidence table             |
| Needs         | `--fasta`: the 5' UTR and coding sequences are read from the reference through the exon model; without it no row is annotated and the run warns once. Ensembl VEP 116.2 warns once too and goes on, its exon slices reading as `N`, so it prints `-` for the two consequence fields and three `0` counts on every 5' UTR row |

| Parameter               | Meaning                                                                                                     |
| ----------------------- | ----------------------------------------------------------------------------------------------------------- |
| `file=<path>`           | The uORF evidence table (`chr pos gene strand type stop_pos`); without it every `Evidence` value is `NA`. A path that names no file is refused when the plugin loads (`Failed to instantiate plugin UTRAnnotator: UTRAnnotator: uORF file <path> not found`, a warning, or the run's error under `--safe`) and the run goes on without the plugin. Ensembl VEP's plugin tests the path, prints a warning on standard output and annotates with every `Evidence` value `NA`, as if `file=` were absent, although its own documentation marks the parameter `(Required)` (`UTRAnnotator.pm`, lines 40, 92 and 109-112; VEP_plugins at release/116). No committed run names a missing file (tier 3 of [Intended divergences from Ensembl VEP](../intended-divergences.md); 4 flag sets run the plugin). |
| `max_overlap=<percent>` | Leave out a UTR exon the variant covers at or above this share; the default is 100.                          |
| a token without `=`     | Not read, as Ensembl VEP's plugin reads it (its parameter hash is empty without `=`), and noted in the warnings file; the plugin runs without the file. |

```bash
--fasta reference.fa \
--plugin UTRAnnotator,file=/path/to/uORF_5UTR_GRCh38_PUBLIC.txt
```

Five analyses run over the 5' UTR with the allele spliced in:

| Analysis       | What it reports                                        |
| -------------- | ------------------------------------------------------ |
| `uAUG_gained`  | a start codon gained                                   |
| `uAUG_lost`    | a start codon lost                                     |
| `uSTOP_gained` | a uORF's stop gained, substitutions alone              |
| `uSTOP_lost`   | a uORF's stop lost                                     |
| `uFrameshift`  | a frameshift in a uORF, insertions and deletions alone |

`5UTR_consequence` names the analysis that fired and `5UTR_annotation` its
effects: the Kozak context and strength around the start codon, its distance to
the cap and to the coding sequence, and the evidence from the uORF table.

A symbolic deletion inside one 5' UTR exon is analysed as its stated span's
deletion. An `<INS>` whose `END` is absent, `POS` or `POS` plus `SVLEN` gets
the frameshift its length decides and no key its bases would decide. Other
symbolic alleles keep `-` in `5UTR_consequence` and `5UTR_annotation`
([the rule](../output/structural-variants.md#utrannotator-on-a-symbolic-allele-the-text-del-scored-as-a-five-base-deletion)).

When more than one analysis fires, the flag is the first of `uAUG_gained`,
`uSTOP_gained`, `uSTOP_lost`, `uAUG_lost`, `uFrameshift`. Ensembl VEP keeps
the last in its hash order
([its own page](../divergences/utrannotator-consequence-precedence.md)). The
`key=value` pairs of `5UTR_annotation` are written in byte order, capitalised
keys first, where Ensembl VEP writes [hash order](../divergences/utrannotator-annotation-pair-order.md). The comparison tools read them
as a set.

## MaxEntScan {#maxentscan}

MaxEntScan scores the splice-site strength of the reference and the alternate
sequence around a variant with the maximum-entropy models in the matrix
files the Ensembl VEP plugin opens (`MaxEntScan.pm` lines 104-105 and 824-831).

| Property      | Value                                                                                                                                      |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| Builds        | GRCh37, GRCh38                                                                                                                             |
| Output fields | `MaxEntScan_ref`, `MaxEntScan_alt`, `MaxEntScan_diff`; with `SWA` the `MES-SWA_donor_` and `MES-SWA_acceptor_` fields; with `NCSS` the `MES-NCSS_` fields |
| Data source   | the directory of matrix files the Ensembl VEP plugin opens (`MaxEntScan.pm` lines 104-105 and 824-831): `me2x5`, `splicemodels/splice5sequences` and `splicemodels/me2x3acc1` to `me2x3acc9` |
| Matrix lines  | Where Ensembl VEP's plugin reads on, vep-rs refuses a matrix file by its line: an empty line, text that is no number, a number followed by other text, or a word for an infinity or a NaN other than `inf`, `infinity` or `nan` in any letter case. |
| Needs         | `--fasta`: every window is read from the reference                                                                                          |

| Parameter   | Meaning                                                                                                                         |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `<path>`    | The matrix directory, as the first parameter.                                                                                   |
| `SWA`       | Also run the sliding window: the highest-scoring donor 9-mer and acceptor 23-mer of the reference and of the alternate context, a reference comparison score and the difference. |
| `NCSS`      | Also score the nearest canonical splice sites upstream and downstream of the variant's exon or intron.                           |
| `verbose`   | Also write the sequences, k-mer positions and contexts behind the scores.                                                        |
| other       | A token that is none of the three is not read, as Ensembl VEP's plugin reads it, and is noted in the warnings file.              |

```bash
--fasta reference.fa --plugin MaxEntScan,/path/to/matrices
--fasta reference.fa --plugin MaxEntScan,/path/to/matrices,SWA,NCSS,verbose
```

The donor model scores a 9-mer, the last 3 bases of an exon and the first 6 of
the intron after it. The acceptor model scores a 23-mer, the last 20 bases of an
intron and the first 3 of the exon after it. `MaxEntScan_ref`, `MaxEntScan_alt`
and `MaxEntScan_diff` (reference minus alternate) are written for a single-base
substitution inside one of those canonical windows of one of the transcript's
introns. The window is scored on the reference and on the sequence with the
allele substituted. Every score prints with three decimals.

Within the transcript's span, the introns tested are the ones the substitution
lies in or within three bases of, the list Ensembl VEP's plugin walks. An intron
of 17 bases or more has both canonical windows inside that reach; a shorter
intron's window bases further out are left blank. A substitution outside the
transcript's span is tested against every intron.

## Where the transcript-model plugins differ from Ensembl VEP

Four behaviours of Ensembl VEP 116.2's plugins are not reproduced here. The
first two are classes the comparison tools set aside. The next two, on
NearestExonJB and Downstream, are documented classes no tool codes, since no
committed row falls under them. Each section names the rows a comparison would
set aside.

Six more UTRAnnotator behaviours have their own pages:

- [several analyses firing](../divergences/utrannotator-consequence-precedence.md),
  the class `utrannotator_consequence_precedence`
- [a start codon lost after an insertion](../divergences/utrannotator-lost-start-after-insertion.md)
- [one gained that stood in place](../divergences/utrannotator-gained-start-remade-in-place.md)
- [an ORF off the frameshift list](../divergences/utrannotator-frameshift-at-stop-boundary.md)
- [a lost stop reported as a frameshift](../divergences/utrannotator-lost-stop-written-before-stop-codon.md)
- [a symbolic allele's reading](../output/structural-variants.md#utrannotator-on-a-symbolic-allele-the-text-del-scored-as-a-five-base-deletion)

[Intended divergences from Ensembl VEP](../intended-divergences.md) lists them
all. Its
[How to read the counts](../intended-divergences.md#how-to-read-the-counts)
defines the format comparator and the comparison case the sections below name.
The concordance figures score consequence sets and never read a plugin field, so
nothing is set aside in them.

Sections 3 and 4 carry their tier under that page's rule. Tier 1 is a shape the
chromosome 21 population cuts carry with no flag. Tier 2 is one that needs a
flag or plugin on ordinary data. Tier 3 is one ordinary data does not carry. The counts
are over the cuts' 15,646,481 records, the 111 golden corpora and their 392 flag
sets.

### 1. UTRAnnotator: the window of the gained-start test

**In short.** Ensembl VEP's UTRAnnotator plugin tests whether a 5' UTR variant
creates a start codon over a window sized by the reference allele, not the
alternate. An insertion's reference allele is empty, so its window is four
bases. A start codon the inserted bases make whose third base lies past
the window is missed. vep-rs sizes the window by the alternate allele and
reports the gained start.

#### 1.1 What a reader needs to know first

UTRAnnotator's five analyses and the two fields that report them,
`5UTR_consequence` and `5UTR_annotation`, are described under
[UTRAnnotator](#utrannotator) above.

#### 1.2 The record

A constructed GRCh38 insertion in the 5' UTR of a reverse-strand gene on
chromosome 21, against the Ensembl 116 GRCh38 cache and the published GRCh38
uORF table:

```
21	25734954	utr_ins	A	ACAT	.	.	.
```

Both engines run `--plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt` with
`--fasta`.

#### 1.3 What each engine writes

The row on `ENST00000924761`, one of eleven reverse-strand transcripts whose
rows differ the same way. Ensembl VEP:

```
utr_ins	21:25734954-25734955	CAT	ENSG00000154723	ENST00000924761	Transcript	5_prime_UTR_variant	64-65	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;5UTR_annotation=KozakContext=ATAATGC:type=uORF:KozakStrength=Moderate:DistanceToStop=18:DistanceToCDS=110:CapDistanceToStart=63:Evidence=False;5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=1;Existing_uORFs=1
```

vep-rs:

```
utr_ins	21:25734954-25734955	CAT	ENSG00000154723	ENST00000924761	Transcript	5_prime_UTR_variant	64-65	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1;5UTR_annotation=CapDistanceToStart=64:DistanceToCDS=112:DistanceToStop=207:KozakContext=TAAATGT:KozakStrength=Weak:type=OutOfFrame_oORF;5UTR_consequence=5_prime_UTR_premature_start_codon_gain_variant;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=1;Existing_uORFs=1
```

The `Existing_` counts agree; the `5UTR_annotation` items are in Ensembl VEP's
hash order, read as a set.

#### 1.4 Why Ensembl VEP writes that

`UTRAnnotator.pm`, line 339, sets the variable that sizes the gained-start
window, `$alt_coding`, to the reference allele. Line 345
cuts the window from two bases before the variant to the reference allele's
length plus one base after it. An insertion's window therefore holds four
bases, and a codon the inserted bases begin reaches one base past it.

The transcript's 5' UTR reads `CCCATAATGCATTGC` around the site, with the
start codon 63 bases from the cap. The inserted bases read `ATG` on the
transcript strand and land after the first base of that codon, giving
`CCCATAAATGTGCATTGC`. The codon at 63 is parted and lost, since those three
bases read `AAT`, and the inserted bases stand as a start codon at 64. The
window of four bases ends at the `T` of the gained codon, so Ensembl VEP never
sees it and reports the loss alone.

For a deletion the window is longer than the alternate allele by the deleted
bases, so it reaches an `ATG` that stood after them and reports it gained.
Opened before the UTR's first base, it reads from the sequence's end, where
Perl counts negative indices.

#### 1.5 Why vep-rs's row is the correct one

The mutated UTR sequence decides. An `ATG` stands at 64 in the mutated
sequence where the reference has `TGC`, made by the inserted bases. The
plugin's gained-start test exists to report exactly that. When a start is lost
and another gained on one row, vep-rs reports the gain. A start codon that an
insertion only moves on is neither: vep-rs reports it
[neither gained nor lost](../divergences/utrannotator-lost-start-after-insertion.md).

#### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh38 --plugin UTRAnnotator,file=uORF_5UTR_GRCh38_PUBLIC.txt
```

`record.vcf` holds the record of 1.2, `<cache>` is the Ensembl 116 GRCh38
cache and `<fasta>` the GRCh38 reference FASTA.

Ensembl VEP writes `5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant`
on `ENST00000924761` and its siblings; vep-rs writes
`5_prime_UTR_premature_start_codon_gain_variant` with the Kozak context
`TAAATGT`.

#### 1.7 How the comparator treats this class

What a comparison sets aside for this class is stated on
[its own page](../divergences/utrannotator-gained-start-window-comparison.md).

### 2. MaxEntScan: two splice-site windows over one substitution

**In short.** A substitution can lie inside two of MaxEntScan's canonical
windows at once. In an exon of at most five bases, the acceptor window of the
intron before it and the donor window of the intron after it both cover it. The
plugin scores one window per row.

Ensembl VEP scores the window of whichever of
the two introns its interval tree returns first, the donor before the acceptor
of that intron. vep-rs scores the first window in transcript order, the
upstream intron's acceptor. The engines agree on the transcripts whose tree
returns the upstream intron first and differ on the others. Each MaxEntScan
fixture corpus carries one record of the shape.

#### 2.1 What a reader needs to know first

On the forward strand an intron from `start`
to `end` has its donor window at `start - 3` to `start + 5` and its acceptor
window at `end - 19` to `end + 3`. The reverse strand mirrors both.

The acceptor window of one intron and the donor window of the next therefore
share bases when the exon between them has at most 5 bases. One intron's own
two windows share bases when `start + 5 >= end - 19`, an intron of at most 25
bases. There both engines score the donor, since both test an intron's donor
window before its acceptor window, so that shape is not a divergence.

#### 2.2 The record

A hand-written GRCh38 substitution in the three-base exon `1:201374980-201374982`
of `ENST00000714314` (TNNT2, reverse strand), between the intron
`1:201374983-201377622`, whose acceptor window is `1:201374980-201375002`, and
the intron `1:201373269-201374979`, whose donor window is `1:201374974-201374982`:

```
1	201374981	tnnt2_both_windows	T	C	.	.	.
```

Both engines run `--fasta <fasta> --plugin MaxEntScan,<synthetic_matrices>,verbose`
against the Ensembl 116 GRCh38 cache. The matrix directory is the synthetic one
the `GRCh38-plugin-maxentscan-fixture` golden corpus generates
(`scripts/fixtures/maxentscan_fixture_matrices.py --seed 1`), so the scores are
the fixture's. `verbose` adds the sequence each engine scored.

#### 2.3 What each engine writes

On the synthetic matrices, Ensembl VEP scores the downstream intron's donor
window, a 9-mer:

```
tnnt2_both_windows	1:201374981	C	ENSG00000118194	ENST00000714314	Transcript	splice_region_variant,5_prime_UTR_variant	2153	-	-	-	-	-	IMPACT=LOW;STRAND=-1;MaxEntScan_alt=-0.284;MaxEntScan_alt_seq=GGGGTGGGC;MaxEntScan_diff=4.852;MaxEntScan_ref=4.568;MaxEntScan_ref_seq=GAGGTGGGC
```

On the same matrices, vep-rs scores the upstream intron's acceptor window, a
23-mer:

```
tnnt2_both_windows	1:201374981	C	ENSG00000118194	ENST00000714314	Transcript	splice_region_variant,5_prime_UTR_variant	2153	-	-	-	-	-	IMPACT=LOW;STRAND=-1;MaxEntScan_alt=-4.441;MaxEntScan_alt_seq=GGCACGCGTGAGTGGGGCAGGGG;MaxEntScan_diff=1.000;MaxEntScan_ref=-3.441;MaxEntScan_ref_seq=GGCACGCGTGAGTGGGGCAGGAG
```

Every other column agrees, and every score here is the synthetic matrices'. A
substitution at `1:201374979`, inside the donor
window alone, scores the same 9-mer in both engines (`MaxEntScan_ref=4.568`).
One at `1:201374983`, inside the acceptor window alone, scores the same 23-mer
(`MaxEntScan_ref=-3.441`). On the three-base exon `1:45616701-45616703` of
`ENST00000531612` (NASP, forward strand) the substitution `1 45616702 A G` lies
in two windows too. Both engines score the upstream intron's acceptor
window, `TTTCCCACCTTGTCAGAAAGAAG`, `MaxEntScan_ref=6.526`.

#### 2.4 Why Ensembl VEP writes that

`MaxEntScan.pm`, line 270, walks the introns `_overlapped_introns` returns and,
for each in turn, returns the donor window's scores when the substitution lies
in it (lines 294-310) and otherwise the acceptor window's (lines 312-328). The
first window met wins, and nothing is written twice. The introns come from the
`Set::IntervalTree` of `ensembl-variation`'s `BaseTranscriptVariation.pm`
(line 813; each intron inserted over its span plus three bases on each side,
line 1027). Its traversal order depends on the transcript's whole set of
introns and is neither genomic nor transcript order.

For `ENST00000714314` the tree returns the downstream intron first, so its donor
window is tested first and scored. For `ENST00000531612` it returns the upstream
intron first, whose donor window misses and whose acceptor window is scored.
Four TNNI1 transcripts share the four-base exon `1:201417116-201417119`. On one
the tree returns the upstream intron first and on three the downstream one. The
same substitution is therefore scored in the acceptor window on one transcript
and the donor window on the others.

#### 2.5 Why vep-rs's row is the correct one

The substitution lies in both windows, and the plugin writes one. The record
decides nothing between them, so the choice must be a rule of the transcript
model. vep-rs's rule is the first window in transcript order, the upstream
intron's acceptor, which the exon and intron coordinates alone determine. The
same substitution in the same exon is therefore scored in the same window on
every transcript that shares the exon.

Ensembl VEP's choice is the traversal order of
a search tree built from the transcript's intron set. The TNNI1 transcripts show
it scoring the same site in different windows for no reason in the record.

#### 2.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --assembly GRCh38 --plugin MaxEntScan,<synthetic_matrices>,verbose
```

`record.vcf` holds one record of the table, `<cache>` and `<fasta>` are those
of 1.6, and `<synthetic_matrices>` is the matrix directory of 2.2.

| Record | Ensembl VEP `MaxEntScan_ref_seq` | vep-rs `MaxEntScan_ref_seq` |
| --- | --- | --- |
| `1 201374981 T C` (TNNT2, `ENST00000714314`) | `GAGGTGGGC`, the donor 9-mer | `GGCACGCGTGAGTGGGGCAGGAG`, the acceptor 23-mer |
| `1 45616702 A G` (NASP, `ENST00000531612`) | `TTTCCCACCTTGTCAGAAAGAAG` | the same |
| `1 201374979 C T` (donor window alone) | `GAGGTGGGC` | the same |

Both engines' scores are the synthetic matrix directory's: a matrix set decides
the numbers and never the window chosen.

#### 2.7 How the comparator treats this class

The format comparator, the plugin tuple scorer and the corpus builder apply
`maxentscan_overlapping_windows` only on a row whose one-base record lies in
both windows of its transcript, read from the candidate's JSON cache
(`--json-cache`). Without it nothing is set aside. There a `MaxEntScan_` number
both engines fill differently is set aside, with the two window texts under
`verbose`, counted once per row. A field one engine alone fills stays a
difference.

Each MaxEntScan fixture corpus, scored on its synthetic matrix set, holds one
such record, the GRCh38 one that of 2.2:

| Corpus | Record, transcripts | vep-rs's acceptor window, Ensembl VEP's donor window | `MaxEntScan_ref` |
| --- | --- | --- | --- |
| `GRCh38-plugin-maxentscan-fixture` | `1 201374981 T C`, `ENST00000714314` | `1:201374980-201375002`, `1:201374974-201374982` | Ensembl VEP `4.568`, vep-rs `-3.441` |
| `GRCh37-plugin-maxentscan-fixture` | `21 30973548 A G`, `ENST00000309434` and `ENST00000535441` | `21:30973547-30973569`, `21:30973542-30973550` | Ensembl VEP `9.975`, vep-rs `9.468` |

### 3. NearestExonJB: a value on an upstream or downstream row that depends on the records before it

**The defect in one sentence.** Ensembl VEP's NearestExonJB plugin widens a
frameshift-intron transcript's exons by 12 bases only once a record inside the
transcript has flagged it, so a record within 12 bases outside its first or last
exon gets a value or nothing by the records before it.

**Why this is a defect.** The widening is a property of the transcript model:
its own comment says the exons are stretched when the transcript has a
frameshift intron. The flag is set while the transcript's intron index is built,
which only a record overlapping the transcript does. A record outside the
transcript reads the flag as it finds it, so one record gets two answers in two
files.

**What vep-rs does instead.** vep-rs decides the widening from the transcript's
introns on every row. A record within 12 bases outside such a transcript's
first or last exon therefore reports that exon whatever the records around it.

**Where it is in Ensembl VEP's source** (`VEP_plugins` and `ensembl-variation`
at release/116):

- `NearestExonJB.pm`, line 184, `run`: the exons overlapping the variant are
  read from `_overlapped_exons` with no argument. The plugin's own call
  therefore fills the memo for a record outside the transcript.
- `Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm`, lines 858-865,
  `_overlapped_exons`: the stretch is 12 when the transcript's
  `_has_frameshift_intron` flag is set, else 0. Lines 1032-1034,
  `_create_intron_trees`: the flag is set on the transcript object while its
  intron index is built.
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 480-500,
  `_bvfo_preds`: the index is built only for a variant within the feature.

**What correcting it would change.** Every row of such a transcript would read
the stretch from the model, as vep-rs's rows do.

**The class it produces.** `nearest_exon_jb_stretch_priming`: a `NearestExonJB`
(JSON `nearestexonjb`) the candidate alone fills on a row whose `DISTANCE` is at
most 12. The `GRCh38-plugin-computed` corpus's cache holds no transcript with an
intron of 13 bases or fewer, so no row of its two comparison cases is set aside.

**How often the shape occurs.** Tier 2: the shape is a transcript with an
intron of 13 bases or fewer, 78 of the 77,601 transcripts of the committed
caches (1,005 per million), in 12 corpora. A record within 12 bases of such an
exon boundary shows it once an earlier record primed the stretch; 2 flag sets
run the plugin.

#### 3.1 A substitution nine bases before the first exon of a lncRNA with two 2-base introns

**In short.** `ENST00000624847` (GRCh38 chromosome 21, forward strand) has
2-base introns at 6358575-6358576 and 6358739-6358740; its first exon,
`ENSE00003757672`, spans 6354264-6354522. For `21 6354255 . G T`, an
`upstream_gene_variant` row with `DISTANCE=9`, vep-rs writes
`NearestExonJB=ENSE00003757672|9|start|259`, the exon 9 bases away, 259 bases
long. Ensembl VEP 116.2 writes it only when a record inside the transcript,
`21 6358575 . G T` in the first short intron, precedes it in the same buffer, so
the upstream row of a position-sorted file gets no value. The
intronic record's rows, `ENSE00003758864|1|end|20`, agree in every file.

Nine bases past the last exon, `21 6359776 . A T` gets
`ENSE00003755945|9|end|1027` from vep-rs whatever the file holds and from Ensembl
VEP only after the intronic record.

### 4. Downstream: the rotation of a shifted allele on a reverse-strand transcript

**The defect in one sentence.** Under `--hgvs`, on a reverse-strand transcript
Ensembl VEP rotates a 3'-shifted allele's bases by the allele's length less the
shift, so a longer shift rotates nothing and its Downstream plugin splices
unrotated bases, not the variant's unless the shift is a multiple of that
length.

**Why this is a defect.** An inserted or deleted sequence moved 3' along the
transcript keeps its bases and changes only which comes first. `TG` moved one
place along a `TGTG` run reads `GT`. Counting `length - shift` places on a
reverse-strand transcript gives a negative count whenever the shift exceeds the
length, and a negative count rotates nothing. The plugin therefore reads an
unrotated spelling at the shifted position.

**What vep-rs does instead.** vep-rs's Downstream plugin reads the record's own
span and the allele as written, with and without `--hgvs`, so it rotates
nothing. Where vep-rs rotates a shifted allele, for HGVSp and the
`Allele` column, it turns the bases as read on the transcript's strand by the
shift modulo the allele's length. They are therefore the variant's own however
long the shift.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1330-1345,
  `shift_feature_seqs`: `shift_length` single-base rotations on a forward-strand
  transcript and `seq_length - shift_length` on a reverse-strand one (lines
  1337-1339), none when that count is zero or negative. Line 1778,
  `hgvs_protein`, calls it, and the Downstream plugin splices the rotated
  `feature_seq` into the coding sequence (`Downstream.pm`, line 110).

**What correcting it would change.** On a transcript of either strand the bases,
as read on its strand, would turn by the shift modulo the allele's length. The
HGVSp would describe the variant's own bases. The Downstream plugin, while
it reads the shifted span, would splice the
variant's own bases there: `SP` on the example below, the last two residues of
vep-rs's `ALSP`.

#### The class for a rotated allele and how often it occurs

**The class it produces.** `downstream_reverse_strand_rotation`: a
`DownstreamProtein` or `ProteinLengthChange` (JSON `downstreamprotein`,
`proteinlengthchange`) that differs on a frameshift row of a reverse-strand
transcript whose `HGVS_OFFSET`, read without its sign, exceeds the inserted or
deleted length and is not a multiple of it, under `--hgvs`. Every row the class
covers is also a row of
[`plugin_hgvs_shifted_cds_coordinates`](../divergences/plugin-hgvs-shifted-cds-coordinates.md):
Ensembl VEP's plugin reads the shifted span, here with unrotated bases, and
vep-rs's the record's own. The genome-wide ClinVar datasets count such rows
under that class and none under this one. Ensembl VEP's rotation reaches HGVSp
and the `Allele` column too
([Allele on a coding row after a shift](../output/vep-and-tab.md#allele-on-a-coding-row-after-a-shift)).

**How often the shape occurs.** Tier 2, a flag and a plugin being needed. The
count is of the records whose allele a shift can rotate, not of the class's
rows: an insertion or deletion of two or more bases not all one letter. They are
1,203,766 of the 15,646,481 records of the five chromosome 21 population cuts
(76,935 per million), in 89 golden corpora.

The transcript's strand and the shift's length decide
whether the rotation goes wrong. The `GRCh38-plugin-computed` corpus runs the
plugin under `--hgvs` on 13 such records, and none of its 139 reverse-strand
frameshift rows shifts at all.

#### 4.1 A two-base insertion at the end of an `AG` run in CLDN17

**In short.** `21 30166352 cldn17_ins_AG G GAG` inserts `AG` after the run
`GAGAGAG` at 21:30166346-30166352, in the single exon of CLDN17
(`ENST00000286808`, reverse strand). On the transcript the allele is `CT`,
shifted seven bases to `c.271_272dup`, `HGVS_OFFSET=-7` on both engines. Ensembl
VEP 116.2 counts `2 - 7` rotations, none, and splices `CT` unturned at the
shifted point: `TCC TCT TGA`, `p.Leu92SerfsTer2`, `DownstreamProtein=SS`.
Turned one place (7 modulo 2) the bases are `TC`, codons `TCT CCT TGA`, and
vep-rs writes `HGVSp=ENSP00000286808.3:p.Leu92ProfsTer2`.

vep-rs's Downstream plugin reads the record's own point, between coding bases
265 and 266, and writes `DownstreamProtein=ALSP` from codon 89, with and
without `--hgvs`. `HGVSc`, `Codons` and `ProteinLengthChange=-132` agree.

Written one base on, `21 30166351 A AGA` shifts six bases: both engines write
`ProfsTer2`, vep-rs `ALSP` again and Ensembl VEP `SP`, its last two residues.
Both shifts stay inside one exon. Shifted six, a multiple of the allele's
length, the unturned bases are the variant's own and both readings spell one
changed sequence. Shifted seven they spell two, the reverse-strand case
[the shifted-coordinates page](../divergences/plugin-hgvs-shifted-cds-coordinates.md)
names.

## Where next

- [Plugins](../plugins.md): every plugin's support status and how one is
  enabled.
- [Intended divergences from Ensembl VEP](../intended-divergences.md): the
  index of every documented difference between the engines.
- [Reproducing these examples](../divergences/reproducing.md): the engines,
  inputs and commands behind the examples.
