# The first codon of the Downstream protein for a deletion from a codon's third base

**The defect in one sentence.** When a deletion, with or without inserted bases, starts on the third base of a codon, Ensembl VEP's Downstream plugin starts its downstream translation at the next codon and leaves out the codon the deletion changes. Its `DownstreamProtein` then lacks the first changed residue, lacks an unchanged residue the variant lies in, or carries residues past a stop the changed codon spells. vep-rs starts the translation at the first codon the change touches, and writes an empty value where that codon is a stop.

**Why this is a defect.** A deletion from base 3k of the coding sequence leaves bases 3k-2 and 3k-1 in place and puts the base after the removed span at position 3k. Codon k is therefore the first codon the change touches, and the changed protein is the reference's first k-1 residues, then whatever codon k encodes, then the rest. The plugin's description says the sequence includes "any amino acids overlapped by the variant itself", and the row's own `ProteinLengthChange` counts codon k. `21 37512050 2081367 AC A` turns codon 604 of ENST00000398960 from `AAC` into `AAA`, `p.Asn604LysfsTer52`. Ensembl VEP's `DownstreamProtein` opens at codon 605 and lacks that lysine.

**What vep-rs does instead.** vep-rs starts the downstream translation at the first codon the change touches. For the record above it writes the lysine and the 50 residues after it. Where the changed codon is a stop it writes an empty value, as both engines do for any frameshift whose first downstream codon is a stop. For an insertion between two adjacent bases vep-rs writes what Ensembl VEP writes, the translation from the first codon that takes an inserted base: an insertion at a codon boundary overlaps no residue before it.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `Downstream.pm`):

- Line 108 reads the allele's CDS start and end, and line 112 takes the lower of the two as `$low_pos`. For an insertion that is the base before the insertion point. For a deletion, with or without inserted bases, it is the first changed base itself.
- Line 113 computes the offset of the downstream translation, `$low_pos - ($low_pos % 3)`. The offset wanted is the number of bases before the first changed base, rounded down to a codon. `$low_pos` is that number for an insertion and one more than it for a deletion. The two round alike except on a codon's third base, where the remainder is zero and the offset is the codon after the changed one.
- Line 115 reads the changed sequence from that offset, and line 139 translates it.
- Line 125 builds the whole changed sequence, the skipped codon included, so the `ProteinLengthChange` of line 148 counts the residue that `DownstreamProtein` leaves out.

**What correcting it would change.** The offset would be `($start - 1) - (($start - 1) % 3)` wherever `$start` is not past `$end`. `DownstreamProtein` would open at the codon the deletion changes on every such row, and would be empty where that codon is a stop, as the row's `HGVSp` and `ProteinLengthChange` say.

**The class it produces.** `downstream_start_past_own_codon`: a `DownstreamProtein` (JSON `downstreamprotein`) that differs on a frameshift row without a splice term, both engines writing one `ProteinLengthChange`. The format comparator decides it from the candidate cache named by `--json-cache`: it recomputes both values from the transcript's coding sequence, 3' UTR and codon table, under the `CDS_position` column as printed. It accepts the row when the first changed base is a codon's third, Ensembl VEP's value is the translation from the next codon and vep-rs's from the changed one. An insertion between two adjacent bases, a field either side omits and any other value stay differences.

A deletion that also runs to the last coding base belongs to [`downstream_utr_frame_past_coding_end`](downstream-utr-frame-past-coding-end.md), and a row on which `--hgvs` changes Ensembl VEP's value to [`plugin_hgvs_shifted_cds_coordinates`](plugin-hgvs-shifted-cds-coordinates.md).

**How often the shape occurs.** Tier 2, the plugin being needed. Under `--hgvs --plugin Downstream`, vep-rs writes the plugin's fields on 557,264 rows of the genome-wide ClinVar GRCh37 file, and 116,617 of them, from 32,235 records, differ between the engines by this class. The changed codon keeps its residue on 87,842 of them, changes it on 25,483 and is a stop on 3,292. On the GRCh38 file 515,443 of 2,495,728 rows differ, from 34,853 records. Five committed flag sets reach the class: 120 rows in each of the two runs of the `GRCh38-plugin-computed` corpus, three rows and nine in the two runs of `GRCh37-plugin-downstream-own-span`, and one row of the `without-hgvs` flag set of `GRCh37-plugin-downstream-coding-end`.

Where a deletion lies in a repeat, and the repeat reaches an exon's last base, the spliced sequence lets the deletion sit further 3' than the coordinates the plugin holds. This page keeps the value that follows from the coordinates the plugin holds, vep-rs's, because the plugin's description says the sequence includes "any amino acids overlapped by the variant itself". The other reading is the translation from the spliced sequence's first changed codon, on such a row a later codon than the first one the change touches at those coordinates.

Of the GRCh37 and GRCh38 rows of this class and `downstream_utr_frame_past_coding_end`, 51 and 246 have that shape. The format comparator accepts all of them, 7 and 37 under `downstream_utr_frame_past_coding_end` and the rest under this class. The other reading's value is Ensembl VEP's on 14 and 71 of them, neither engine's on 37 and 175 and vep-rs's on none.

## 1. A one-base deletion that changes the codon it starts in

**In short.** `21 37512050 2081367 AC A`, a ClinVar record, deletes the last base of codon 604 of ENST00000398960 (DYRK1A, forward strand). The next base completes the codon as `AAA`: asparagine becomes lysine and the frame shifts, `p.Asn604LysfsTer52` on both engines. vep-rs writes the 51 residues from that lysine to the new stop. Ensembl VEP 116.2 writes the last 50 of them: its `DownstreamProtein` opens at codon 605 and lacks the first residue the variant changes.

### 1.1 What a reader needs to know first

The Downstream plugin writes two fields on a frameshift row. `DownstreamProtein` is the translation of the changed transcript from the variant's codon on, 3' UTR included, up to the first stop. `ProteinLengthChange` is the length of the whole changed protein less the reference protein's. Both describe one changed sequence.

The coding sequence of ENST00000398960 is 2,292 bases and its protein 763 residues. Codons 603 to 607 read `GGT AAC AGT TCC CAT`: glycine, asparagine, serine, serine, histidine.

### 1.2 The record

```
21	37512050	2081367	AC	A	.	.	.
```

The record deletes the `C` at 21:37512051, coding base 1812, the third base of codon 604. It is a record of the `GRCh38-plugin-computed` corpus, whose flags both engines run: `--hgvs` and five plugins, Downstream among them.

### 1.3 What each engine writes

`<50 residues>` and `<51 residues>` stand for the residues of `DownstreamProtein`, which this page gives by their count. vep-rs's 51 are `K`, the lysine of codon 604, and then the 50 that Ensembl VEP writes. The command of 1.6 writes vep-rs's row whole, and the corpus's `expected/` directory holds Ensembl VEP's.

Ensembl VEP 116.2, default format:

```
2081367	21:37512051	-	ENSG00000157540	ENST00000398960	Transcript	frameshift_variant	1912	1812	604	N/X	aaC/aa	-	IMPACT=HIGH;STRAND=1;HGVSc=ENST00000398960.8:c.1812del;HGVSp=ENSP00000381932.2:p.Asn604LysfsTer52;DownstreamProtein=<50 residues>;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0;NearestExonJB=ENSE00003827888|140|start|14448;ProteinLengthChange=-109
```

vep-rs:

```
2081367	21:37512051	-	ENSG00000157540	ENST00000398960	Transcript	frameshift_variant	1912	1812	604	N/X	aaC/aa	-	IMPACT=HIGH;STRAND=1;HGVSc=ENST00000398960.8:c.1812del;HGVSp=ENSP00000381932.2:p.Asn604LysfsTer52;DownstreamProtein=<51 residues>;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0;NearestExonJB=ENSE00003827888|140|start|14448;ProteinLengthChange=-109
```

Every other column is equal. The tab, VCF and JSON outputs carry the same two values, the JSON under `downstreamprotein`. The record's rows on 14 more transcripts of the gene differ the same way.

### 1.4 Why Ensembl VEP writes that

Line 113 computes `1812 - (1812 % 3)`, which is 1,812, since 1,812 is three times 604. Offset 1,812 is the first base of codon 605. From codon 603 the changed sequence reads `GGT AAA GTT CCC ATC`, and the translation of line 139 starts at `GTT`, valine. Line 125 builds the whole changed sequence, which translates to 654 residues, so line 148 writes -109.

### 1.5 Why vep-rs's row is the correct one

The record removes coding base 1812 and no other. Bases 1810 and 1811, `AA`, stay where they are, and base 1813, `A`, takes the place of the removed `C`, so codon 604 reads `AAA`. That lysine is the first residue the variant changes. The row's `HGVSp` names it and counts 52 codons from it to the new stop, the stop included: 51 residues. The changed protein is 603 unchanged residues and those 51, the 654 that `ProteinLengthChange=-109` describes against 763. A downstream sequence of 50 residues leaves one of them out.

### 1.6 Try it yourself

From `tests/golden/116/GRCh38-plugin-computed/`:

```
vep -i variants.vcf -o out.txt --offline --json_cache json_cache \
  --assembly GRCh38 --fasta reference.fa.gz --hgvs --plugin Downstream
```

vep-rs writes a `DownstreamProtein` of 51 residues that opens `KVPITITT` on the row. The corpus's `expected/` directory holds Ensembl VEP 116.2's four outputs of the same records under the corpus's flags, with the 50 residues that open `VPITITT`.

### 1.7 How the comparator treats this class

The format comparator accepts the row under `downstream_start_past_own_codon` when it is given the candidate cache: base 1812 is a codon's third, Ensembl VEP's value is the translation of the changed sequence from codon 605 and vep-rs's the translation from codon 604. The corpus manifest documents the `DownstreamProtein` of 120 rows of each run under the class, and the golden harness reads vep-rs's value from it. The same row with a wrong value on either side, or with the field omitted on either side, stays a difference. An empty value is a value and is compared.

## 2. A deletion whose changed codon is a stop

**In short.** `2 179425400 1066899 TG T`, a ClinVar record, deletes the last base of codon 19613 of ENST00000342175 (TTN, reverse strand), `TAC`: coding base 58839, `c.58839del`, with no 3' shift. The next base completes the codon as `TAA`, a stop. Both engines write `p.Tyr19613Ter` and `ProteinLengthChange=-7506`, a protein that ends after residue 19612. vep-rs writes an empty `DownstreamProtein`. Ensembl VEP writes one residue, `S`.

From codon 19612 the changed sequence reads `TAT TAA TCG TAG AAA`. Line 113 gives offset 58839, the first base of codon 19614, so the plugin translates `TCG` and ends at `TAG`: a codon after the stop the variant makes, read to the stop after it. No reading of the changed transcript encodes it, and the row's own notation and length change say the protein ends at codon 19613. The record's rows on five more transcripts of the gene carry the same pair of values.

## 3. A deletion whose changed codon keeps its residue

**In short.** `21 43058850 1162233 CG C` deletes the last base of codon 447 of ENST00000352178 (CBS, reverse strand), `CCC`. The next base completes the codon as `CCG`, proline again, so the first changed residue is 448: `p.Val448TrpfsTer7` and `ProteinLengthChange=-98` on both engines. vep-rs writes `DownstreamProtein=PWWMRRG`, the residue the deleted base lies in and the six new ones. Ensembl VEP writes `WWMRRG`.

From codon 446 the changed sequence reads `GCG CCG TGG TGG ATG`. The deleted base is a base of codon 447, so that codon is one the variant overlaps, and the plugin's description counts it into the sequence. Line 113 starts at the variant's codon when the deleted base is the codon's first or second, and one codon later when it is the third.
