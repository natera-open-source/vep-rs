# The Downstream protein of a deletion that runs to the last coding base

**The defect in one sentence.** When a deletion starts on a codon's third base and runs to the last coding base, Ensembl VEP's Downstream plugin translates the 3' UTR from its first base, one base out of the changed transcript's frame. Its `DownstreamProtein` then spells residues the changed transcript does not encode, while `ProteinLengthChange` on the same row counts the residues of the right frame. vep-rs starts the translation at the first codon the deletion touches and reads on in the changed transcript's frame.

**Why this is a defect.** Such a deletion leaves the first two bases of the codon it starts in, and the first base of the 3' UTR completes that codon. Every later codon of the changed transcript therefore starts on a UTR base numbered 3n + 2, and a translation that starts on the UTR's first base reads another frame. `17 59760656 530351 ATTAC A` leaves `AA` before a UTR that opens `TAATACTTAA`. The changed transcript reads `AAT AAT ACT TAA`, three residues and a stop, as the row's own `HGVSp`, `p.Lys1249AsnfsTer4`, and its `ProteinLengthChange=2` say. Ensembl VEP reads `TAA` first and writes an empty `DownstreamProtein`.

**What vep-rs does instead.** vep-rs starts the downstream translation at the first codon the deletion touches, the two bases it leaves completed by the UTR's first base, and reads on in that frame: `DownstreamProtein=NNT` for the record above. On a transcript whose 3' UTR neither the cache nor the FASTA spells, the two bases code nothing and both engines write an empty value.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `Downstream.pm`):

- Line 110 replaces the allele's CDS span, here with nothing, so the changed coding sequence ends two bases into the codon the deletion starts in.
- Lines 112-113 take the offset of the downstream translation from the first deleted base's own coordinate, `$low_pos - ($low_pos % 3)`. On a codon's third base the remainder is zero and the offset is the coordinate itself, one more than the length of the changed coding sequence.
- Line 115 asks `substr` for the changed sequence from that offset on. The offset lies outside the string, so Perl warns `substr outside of string` and returns nothing, and lines 118-119 join that nothing to the 3' UTR.
- Line 125 builds the whole changed sequence from the bases before that offset and the string of lines 118-119. Asked for the bases before an offset past the end, `substr` returns all of the changed coding sequence, so the `ProteinLengthChange` of line 148 is counted in the right frame.

**What correcting it would change.** The offset would be the number of bases before the first deleted base rounded down to a codon, `($start - 1) - (($start - 1) % 3)`. The downstream translation would start at the codon the deletion starts in, in the changed transcript's frame, and agree with the row's `ProteinLengthChange` and `HGVSp`. The same offset corrects [the start one codon late](downstream-start-past-own-codon.md), the wider class this one is cut from.

**The class it produces.** `downstream_utr_frame_past_coding_end`: a `DownstreamProtein` (JSON `downstreamprotein`) that differs on a frameshift row without a splice term whose allele is a deletion, both engines writing one `ProteinLengthChange`. The class is built into the format comparator, which decides it from the candidate cache named by `--json-cache`. It recomputes both values from the transcript's coding sequence, 3' UTR and codon table, under the `CDS_position` column as printed. It accepts the row only when the deletion starts on a codon's third base and ends on the last coding base, and each engine's value equals its recomputation. A field either side omits, any other value, an insertion and a deletion that ends earlier stay differences.

A row on which `--hgvs` changes Ensembl VEP's value belongs to [`plugin_hgvs_shifted_cds_coordinates`](plugin-hgvs-shifted-cds-coordinates.md).

**How often the shape occurs.** Tier 2, the plugin being needed. Under `--hgvs --plugin Downstream`, vep-rs writes the plugin's fields on 557,264 rows of the genome-wide ClinVar GRCh37 file, and 49 of them have the shape. The 23 on a transcript with a 3' UTR, from 11 records, differ between the engines; the 26 on a transcript without one are equal. On the GRCh38 file 108 rows, from 33 records, differ. Two committed flag sets reach the class, both of the `GRCh37-plugin-downstream-coding-end` corpus: four rows under `--hgvs --plugin Downstream` and four under `--plugin Downstream`.

Tier 2 is the second of the three tiers of [Intended divergences from Ensembl VEP](../intended-divergences.md): ordinary data carries the shape, and a plugin or a flag is needed to reach it. The corpus's `manifest.json` lists the rows of its two flag sets by location and transcript under the class.

## 1. A four-base deletion over the stop codon of BRIP1

**In short.** `17 59760656 530351 ATTAC A`, a ClinVar record, deletes the last base of codon 1249 and the stop codon of ENST00000259008 (BRIP1, reverse strand). Both engines write `frameshift_variant,stop_lost`, `p.Lys1249AsnfsTer4` and `ProteinLengthChange=2`: three new residues and a stop. vep-rs writes those three residues, `DownstreamProtein=NNT`. Ensembl VEP 116.2 writes `DownstreamProtein=`, an empty value, because it reads the 3' UTR from its first codon, `TAA`.

### 1.1 What a reader needs to know first

The Downstream plugin writes two fields on a frameshift row. `DownstreamProtein` is the translation of the changed transcript from the variant's codon on, 3' UTR included, up to the first stop; the plugin's own description says the sequence includes "any amino acids overlapped by the variant itself". `ProteinLengthChange` is the length of the whole changed protein less the reference protein's. Both describe one changed sequence.

The coding sequence of ENST00000259008 is 3,750 bases, 1,249 codons and the stop codon `TAA`; its 3' UTR is 2,030 bases. Codon 1249 is `AAG`, lysine.

### 1.2 The record

```
17	59760656	530351	ATTAC	A	.	.	.
```

On the forward strand the record deletes `TTAC`, 17:59760657-59760660. On the transcript that is `GTAA`: coding bases 3747 to 3750, the `G` of codon 1249 and the three bases of the stop codon. Both engines run `--hgvs --plugin Downstream` on GRCh37 with a reference FASTA.

### 1.3 What each engine writes

Ensembl VEP 116.2, default format:

```
530351	17:59760657-59760660	-	ENSG00000136492	ENST00000259008	Transcript	frameshift_variant,stop_lost	4015-4018	3747-3750	1249-1250	K*/X	aaGTAA/aa	-	IMPACT=HIGH;STRAND=-1;HGVSc=ENST00000259008.2:c.3747_3750del;HGVSp=ENSP00000259008.2:p.Lys1249AsnfsTer4;DownstreamProtein=;ProteinLengthChange=2
```

vep-rs:

```
530351	17:59760657-59760660	-	ENSG00000136492	ENST00000259008	Transcript	frameshift_variant,stop_lost	4015-4018	3747-3750	1249-1250	K*/X	aaGTAA/aa	-	IMPACT=HIGH;STRAND=-1;HGVSc=ENST00000259008.2:c.3747_3750del;HGVSp=ENSP00000259008.2:p.Lys1249AsnfsTer4;DownstreamProtein=NNT;ProteinLengthChange=2
```

Every other column is equal. The tab, VCF and JSON outputs carry the same two values, the JSON under `downstreamprotein`. Ensembl VEP's run also warns `substr outside of string` from line 115 of `Downstream.pm` for the row.

### 1.4 Why Ensembl VEP writes that

Line 110 removes the four bases and leaves 3,746. Line 113 computes `3747 - (3747 % 3)`, which is 3747, since 3747 is three times 1,249. Line 115 asks for the changed sequence from offset 3747, one past its end, and gets nothing. Line 119 joins that to the UTR, which opens `TAA TAC TTA ACT`, and the translation of line 139 ends at its first codon. Line 125 builds the 3,746 bases and the UTR whole, which translate to 1,251 residues, so line 148 writes 2.

### 1.5 Why vep-rs's row is the correct one

The record removes coding bases 3747 to 3750, so 3,746 coding bases remain: 1,248 whole codons and `AA`. A frame set by the start codon takes those two bases with the UTR's first base, `T`, and reads `AAT AAT ACT TAA`: asparagine, asparagine, threonine and a stop. That is the protein the row's own `HGVSp` describes, lysine 1249 replaced by asparagine with the stop four codons on, and the 1,251 residues its `ProteinLengthChange` counts. A translation that opens on the UTR's first base starts one base early and belongs to no reading of the changed transcript.

### 1.6 Try it yourself

The corpus `tests/golden/116/GRCh37-plugin-downstream-coding-end/` holds the record with its cache slice and masked reference:

```
vep -i variants.vcf -o out.txt --offline --json_cache json_cache \
  --assembly GRCh37 --fasta reference.fa.gz --hgvs --plugin Downstream
```

vep-rs writes `DownstreamProtein=NNT` on the row. Ensembl VEP 116.2, run on the same record against its own GRCh37 cache and the reference FASTA, writes `DownstreamProtein=`; the corpus's `expected/` directory holds its four outputs. Dropping `--hgvs` changes neither value, the record having no 3' shift.

### 1.7 How the comparator treats this class

The format comparator accepts the row under `downstream_utr_frame_past_coding_end` when it is given the candidate cache: the deletion starts on base 3747, a codon's third, and ends on base 3,750, the last coding base, and the two values are the two readings of the cache's sequences. The corpus manifest documents the row's `DownstreamProtein` under the class, and the golden harness reads vep-rs's value from it. The same record with a wrong value on either side, or with the field omitted on either side, stays a difference. Ensembl VEP's empty value is the value its reading gives, and it is accepted.

## 2. A twenty-two-base deletion that ends a globin's coding sequence

**In short.** `16 223577 4685676 <23 bases> T` deletes coding bases 408 to 429 of ENST00000251595 (HBA2, forward strand): the last base of codon 136, the six codons after it and the stop codon. `<23 bases>` stands for the reference allele, `T` and the 22 deleted bases after it, which the corpus's `variants.vcf` holds whole. The bases left, `GT`, and the UTR read `GTG CTG GAG CCT CGG TAG`. Both engines write `ProteinLengthChange=-2`, a protein of 140 residues against 142. vep-rs writes `DownstreamProtein=VLEPR`, the five residues from codon 136 to that stop. Ensembl VEP writes 30 residues that open `AGASV`, the UTR read from its first base, where its own length change leaves room for five.

The record's 3' shift runs past the coding sequence (`c.412_*4del`, `HGVS_OFFSET=4`), so Ensembl VEP computes no protein notation and the plugin reads the column as printed, with and without `--hgvs`. The row on ENST00000397806, a second transcript of the gene, carries the same pair.
