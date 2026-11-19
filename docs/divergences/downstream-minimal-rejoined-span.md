# The Downstream protein of a `--minimal` split allele

**The defect in one sentence.** Under `--minimal`, Ensembl VEP's default and tab formats run the Downstream plugin on the rejoined record's CDS coordinates, so `DownstreamProtein` and `ProteinLengthChange` on a split insertion's row describe the allele's bases written over the record's span, a change the allele does not make. vep-rs writes the two fields from each split allele's own CDS coordinates, as on a record of its own, in every output format.

**Why this is a defect.** The plugin predicts the protein a frameshift leaves: it writes the allele into the coding sequence at the row's CDS coordinates and translates the changed sequence from there to the first stop. `--minimal` annotates each allele of a multi-allelic record on its own trimmed span, and the row's terms, codons and amino acids are that span's. The plugin runs after the rejoin and reads the record's span, which the other alleles set. An insertion written beside a deletion is spliced in place of the deleted base, and the same insertion alone on a record gets another value. Ensembl VEP's VCF and JSON formats write the allele's own value for the same record and flags.

**What vep-rs does instead.** vep-rs reads each split allele's own CDS coordinates in every output format. The two fields equal the ones the allele gets on a record of its own, whatever alleles stand beside it. The plugin's fields are described under [Downstream](../plugins/transcript-model.md#downstream).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-vep` at release/116.2, `ensembl-variation` at release/116):

- `Downstream.pm`, line 108: `my ($start, $end) = ($tv->cds_start, $tv->cds_end);` reads the row's CDS coordinates, and line 110 replaces that span of the coding sequence with the allele's bases. Lines 112-113 take the lower of the two coordinates, round it down to a multiple of three and start the translation after that many bases.
- `OutputFactory.pm`, lines 2284-2316: the rejoin annotates each split feature on its own span (lines 2289 and 2303), then writes the record's span back onto the first feature (lines 2292-2293) and points every merged feature's TranscriptVariations at it (line 2314). The output pass that runs the plugin (line 2226) follows.
- `OutputFactory.pm`, lines 294-298: the default and tab formats rejoin, then clear every cached CDS position (`clear_shifting_variables`, `TranscriptVariationAllele.pm` lines 509-511), so the output pass maps `$tv->cds_start` and `$tv->cds_end` from the record's span (`BaseTranscriptVariation.pm` lines 511-513). `OutputFactory/VCF.pm` lines 305-308 and `OutputFactory/JSON.pm` lines 217-220 clear first and rejoin after, so those formats keep the coordinates cached on the allele's own span.

**What correcting it would change.** The plugin would read each allele's own CDS coordinates, and a split allele's two fields would equal the ones it gets on a record of its own, in every format.

**The class it produces.** `downstream_minimal_rejoined_span`: a `DownstreamProtein` or `ProteinLengthChange` that differs between the engines in the default or tab format, on a row of a multi-allelic record whose alleles trim to different spans under `--minimal`. No comparison tool sets a row of the shape aside: no committed flag set runs the plugin with `--minimal`, so a disagreement of this shape would be reported as a difference and charged to vep-rs. The concordance figures score consequence sets and never see the fields.

**How often the shape occurs.** Tier 2, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape the record alone shows is a multi-allelic record whose sequence alleles trim to different spans, 1,894 of the cuts' 15,646,481 records (121.0 per million), in 24 corpora. The transcript decides whether a split allele is an insertion or a deletion on a frameshift row the plugin writes, and whether its own CDS coordinates and the record's differ. Of the committed flag sets, none runs the plugin with `--minimal`.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A one-base insertion written beside the deletion of the base before it

**In short.** `21 36896882 ins_beside_del TG TGA,T` on GRCh38 lies in the coding sequence of HLCS on ENST00000336648. Under `--minimal`, `TGA` trims to a one-base insertion after 36896883 and `T` to the deletion of the base at 36896883, coding base 1428. In the default and tab formats Ensembl VEP 116.2 writes a `DownstreamProtein` of 250 residues and `ProteinLengthChange=0` on the insertion's frameshift row, the tail of the unchanged protein. vep-rs writes 100 residues and `-151` in every format, as both engines do for the insertion alone.

### 1.1 What a reader needs to know first

The Downstream plugin writes two fields on a row whose terms include `frameshift_variant` and no splice term. It builds the changed coding sequence by replacing the bases from the row's CDS start to its CDS end with the allele's bases. For an insertion the start lies one base past the end, so nothing is replaced. The translation starts after the lower coordinate rounded down to a multiple of three, runs on into the 3' UTR and ends at the first stop; `DownstreamProtein` is its residues. `ProteinLengthChange` is the changed protein's length less the unchanged protein's.

Under `--minimal` Ensembl VEP trims each allele of a multi-allelic record, annotates the trimmed alleles one at a time and rejoins their rows under the record's span. The `Location` column prints the record's span, and in the default and tab formats the `CDS_position` column does too.

### 1.2 The record

A multi-allelic record against the Ensembl 116 GRCh38 cache. Its second alternate allele is ClinVar variation 1069237, a record of the `GRCh38-plugin-computed` golden corpus; the first is written beside it for this example.

```
21	36896882	ins_beside_del	TG	TGA,T	.	.	.
```

Both engines drop the shared first base, so the record spans the one base 36896883. On ENST00000336648, a reverse-strand transcript, that base is coding base 1428, the third base of codon 476. The first alternate allele trims to the insertion of one base between 36896883 and 36896884, that is between coding bases 1427 and 1428. The second is the deletion of coding base 1428. Both are frameshifts. The transcript's coding sequence is 2,181 bases long, 726 codons and the stop codon, so its protein has 726 residues.

Both engines run `--minimal --plugin Downstream` with `--fasta`, and `--tab`, `--vcf` or `--json` for the other formats.

### 1.3 What each engine writes

The insertion's row on ENST00000336648 in the default format. `<250 residues>` and `<100 residues>` stand for the residues of `DownstreamProtein`, which this page gives by their count. Ensembl VEP 116.2:

```
ins_beside_del	21:36896883	A	ENSG00000159267	ENST00000336648	Transcript	frameshift_variant	1899	1428	476	P/PX	ccc/ccTc	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1;MINIMISED=1;DownstreamProtein=<250 residues>;ProteinLengthChange=0
```

vep-rs:

```
ins_beside_del	21:36896883	A	ENSG00000159267	ENST00000336648	Transcript	frameshift_variant	1899	1428	476	P/PX	ccc/ccTc	-	ALLELE_NUM=1;IMPACT=HIGH;STRAND=-1;MINIMISED=1;DownstreamProtein=<100 residues>;ProteinLengthChange=-151
```

The record writes 34 rows in each engine, 17 per allele, equal in every column but these two fields. The insertion's values on ENST00000336648, by format:

| Format | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- |
| default and tab | 250 residues, `0` | 100 residues, `-151` |
| VCF and JSON | 100 residues, `-151` | 100 residues, `-151` |

The insertion has ten frameshift rows, and in the default and tab formats all ten differ:

| Transcript | Ensembl VEP 116.2 | vep-rs |
| --- | --- | --- |
| ENST00000336648 and five others | 250 residues, `0` | 100 residues, `-151` |
| ENST00000675057, ENST00001083022 | 128 residues, `0` | 100 residues, `-29` |
| ENST00001017196 | 157 residues, `0` | 100 residues, `-58` |
| ENST00001111936 | 44 residues, `0` | 61 residues, `16` |

The five others are ENST00000399120, ENST00000612277, ENST00000674895, ENST00000675307 and ENST00000876941. The deletion's ten frameshift rows carry the same `ProteinLengthChange` in both engines in every format, `-242` on ENST00000336648: the record's span is the deletion's own.

### 1.4 Why Ensembl VEP writes that

The rejoin writes the record's span back onto the first allele's feature and points the other allele's TranscriptVariations at it. The default and tab formats then clear the cached CDS coordinates, so the plugin's read maps them from the record's span: `cds_start` 1428 and `cds_end` 1428 on the insertion's row, the position its `CDS_position` column prints. Line 110 replaces coding base 1428 with the allele's one base. One base in place of one leaves the sequence's length and its frame as they were.

Line 113 starts the translation after 1,428 bases, 476 whole codons, since 1,428 is a multiple of three. The 753 bases left are 250 codons and the stop codon, so the value is residues 477 to 726 of the unchanged protein, and 476 + 250 = 726 gives a change of 0. In the VCF and JSON formats the coordinates cached on the allele's own span survive, `cds_start` 1428 and `cds_end` 1427. Line 110 then replaces nothing and inserts the base, and the plugin writes 100 residues and `-151`.

### 1.5 Why vep-rs's row is the correct one

The row names its allele: one base inserted between coding bases 1427 and 1428, with the codons `ccc/ccTc`, the amino acids `P/PX` and the term `frameshift_variant` in both engines. An insertion keeps every reference base, so the changed sequence is one base longer and every codon from 476 on is read in a new frame. The lower coordinate, 1427, rounds down to 1425, so the translation starts after 475 whole codons and meets a stop after 100 residues: 475 + 100 = 575 = 726 - 151.

A changed protein as long as the unchanged one, in the unchanged frame, belongs to a substitution of coding base 1428, which the record does not hold. That base is in the record's span only because the deletion beside the insertion removes it. Ensembl VEP agrees when the insertion is given alone: `21 36896883 ins_alone G GA` under the same flags writes 100 residues and `-151` on the ENST00000336648 row in both engines, in every format.

### 1.6 Try it yourself

With the record of 1.2 saved as `record.vcf`, the command as vep-rs takes it:

```
vep -i record.vcf -o out.txt --offline \
  --json_cache <corpus>/json_cache --fasta <corpus>/reference.fa.gz \
  --species homo_sapiens --assembly GRCh38 \
  --minimal --plugin Downstream
```

`<corpus>` is `tests/golden/116/GRCh38-plugin-computed`. vep-rs writes `ProteinLengthChange=-151` on the insertion's ENST00000336648 row. Ensembl VEP 116.2 takes the same flags with the cache options of [Reproducing these examples](reproducing.md). It writes `ProteinLengthChange=0` on that row, and `-151` under `--vcf` or `--json`. Writing the insertion on a record of its own gives `-151` in both engines.

### 1.7 How the comparator treats this class

Nothing is set aside. The format comparator (`scripts/concordance/compare_vep_formats.py`) has no class of this shape. The golden harness compares the two fields on every row of the corpora that run the plugin, and none of their flag sets passes `--minimal`. A disagreement of the shape would be reported as a difference and charged to vep-rs.
