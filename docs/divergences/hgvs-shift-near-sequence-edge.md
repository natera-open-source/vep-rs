# The 3' shift of an indel within 1,000 bases of a sequence region's edge

**The defect in one sentence.** Ensembl VEP reads the flanking sequence it shifts an insertion or deletion along out of a window cut at fixed 1,000-base offsets from the ends of a slice it has already clipped to the sequence region. Within 1,000 bases of a chromosome's or an organelle genome's edge the flank it reads is therefore not the flank beside the variant, and the HGVS shift it reports is not the shift the sequence admits. It stops short of that shift in each example on this page, and can run past it where the misread bases happen to match. vep-rs reads the flank the sequence region holds beside the variant and reports the shift the sequence admits.

**Why this is a defect.** HGVS notation places an insertion or a deletion at the most 3' position of the transcript at which the altered sequence is the same, and both engines shift the allele for that purpose. The shift is a property of the bases around the variant alone. Ensembl VEP's window logic is right everywhere except where the window would run past the start or end of the sequence region. There the clipping shortens the sequence it reads while the offsets it cuts at stay fixed, and the shift is computed against bases that are not adjacent to the variant.

**What vep-rs does instead.** vep-rs reads the bases on each side of the variant from the reference FASTA, as many as the sequence region holds, and shifts the allele along them until the altered sequence changes. At the edge of a sequence region the available flank is shorter and the shift is computed against it. A deletion that contains a transcript whole before and after the shift keeps the notation of the whole transcript and carries the `HGVS_OFFSET` of the shift, as the allele's other rows on that strand do.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 427-437, `_genomic_shift`: the slice around the variant is expanded to 1,000 bases on each side and then `constrain_to_seq_region` clips it. `pre_seq` is the first 1,000 bases of the clipped slice and `post_seq` its last 1,000, which are the variant's flanks only when nothing was clipped.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, line 462: `perform_shift` is called with the transcript's strand, and lines 315 and 320 walk `pre_seq` for a reverse-strand transcript and `post_seq` for a forward-strand one, so only the clipped side's flank misleads the walk.
- `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm` (`ensembl-vep` at release/116.2), line 122: the variant takes the slice of the first transcript annotated against it, as the cache stores it, so the clipping ends at the region length the cache records. A FASTA record of another length does not move the edge.

**What correcting it would change.** This class vanishes: the HGVSc, the HGVSp and `HGVS_OFFSET` of an indel near a sequence edge carry the shift the sequence supports, as they do everywhere else. A row that another documented defect of Ensembl VEP reaches wherever the flank is read whole, such as the whole-transcript rows of 1.1, then differs for that reason and falls to its class.

**The class it produces.** `hgvs_shift_sequence_edge_window`, documented, where a corpus carries such a record, in its manifest for a row whose `HGVSc`, `HGVSp` or `HGVS_OFFSET` differs for this reason, and applied by the golden harness and the format comparator from the manifest. The concordance figures score consequence sets and never see the fields, so nothing is set aside in them.

Under `--shift_3prime 1 --shift_length` the same walk is printed as `SHIFT_LENGTH`. The deletion `MT:6468-15586` of the human mitochondrial genome ends 983 bases before the region's end and admits a shift of two. On its 27 forward-strand rows, on GRCh38 and on GRCh37 alike, Ensembl VEP writes `SHIFT_LENGTH=0` and vep-rs `SHIFT_LENGTH=2`, with the same `Consequence`; on its nine reverse-strand rows both write 0. On the two transcripts it overlaps in part, `cDNA_position` and `CDS_position` stand two places further on vep-rs's row and `Protein_position` one place further on one of them; under `--shift_3prime 1` alone the same fields differ on those rows. The class covers none of these fields, so these rows stay outside it, though their cause is this defect.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion 291 bases into the Arabidopsis mitochondrial genome

**In short.** A three-base deletion at `Mt:291-293` of the Arabidopsis thaliana TAIR10 assembly lies in a repeat that admits a shift of four bases towards the transcript's 3' end. Ensembl VEP shifts it two bases, `c.444_446del` with `HGVS_OFFSET=-2`; vep-rs shifts it four, `c.446_448del` with `HGVS_OFFSET=-4`. The same deletion anywhere further than 1,000 bases from the sequence's start is shifted alike by both engines.

### 1.1 What a reader needs to know first

For HGVS, an insertion or deletion inside a repeated sequence is written at the 3'-most position along the transcript at which deleting or inserting the bases gives the same result, and `HGVS_OFFSET` reports how far the written position moved from the input's. On a reverse-strand transcript the 3' direction is towards lower genomic coordinates.

A sequence region is one FASTA record: a chromosome, a plastid or a mitochondrial genome. A variant within 1,000 bases of its first or last base has less than 1,000 bases of flank on that side.

The 3' walk reads one flank: the downstream flank for a forward-strand transcript and the upstream flank for a reverse-strand one (`perform_shift`, lines 315 and 320, chosen by the strand at line 462). The shift is therefore wrong for a reverse-strand transcript near a region's start and for a forward-strand transcript near its end, and right for the other two combinations, whose clipped flank the walk never reads. Walked against the wrong bases it can also run further than the true flank admits, so the class covers a shift that differs for this reason in either direction.

The forward-strand case shows on the human mitochondrial genome, 16,569 bases long. A deletion of the one base at `MT:15934`, 635 bases before the region's end, lies in the forward-strand transcript `ENST00000387460`, on the first base of a run of four identical bases, so the sequence admits a shift of three. Ensembl VEP's clipped slice runs from `MT:14934` to `MT:16569`, and `post_seq`, its last 1,000 bases, starts at `MT:15570` where the flank starts at `MT:15935`. The first comparison fails: Ensembl VEP writes `n.47del` with no `HGVS_OFFSET`, vep-rs `n.50del` with `HGVS_OFFSET=3`.

A deletion that contains a transcript whole carries the same `HGVSc` on both engines, the deletion of the transcript's first to last base, and differs in `HGVS_OFFSET` alone. The deletion `MT:6468-15586` ends 983 bases before the region's end and admits a shift of two: vep-rs writes `HGVS_OFFSET=2` on the fifteen forward-strand transcripts inside it, and Ensembl VEP, whose walk stops at zero, writes none.

The zero offset is also why Ensembl VEP writes a notation on those rows at all: it passes the test at line 1472 of `TranscriptVariationAllele.pm`, which refuses a clamped end carried past the transcript's last base. Where the flank is read whole, as for a deletion more than 1,000 bases from an edge, the walk finds the shift the sequence admits, and where that shift is not zero the test refuses the notation. Ensembl VEP then writes no `HGVSc` on the rows of a transcript the deletion contains whole, a difference of another class, [`hgvs_shifted_indel_written_span_gate`](hgvs-shifted-indel-written-span-gate.md).

### 1.2 The record

The record is a hand-written deletion on the Arabidopsis mitochondrial genome, annotated against the TAIR10 cache of Ensembl Plants release 63 and its reference FASTA (the `TAIR10-arabidopsis` golden corpus, `tests/golden/116/TAIR10-arabidopsis/`, carries the cache; this record is not among its records):

```
Mt	290	mt_del3_290_TCTT_T	TCTT	T	.	.	.
```

The deleted bases are `Mt:291-293`, `CTT`. The bases `Mt:286-296` read `CTCTTCTTGCA`. The transcript `ATMG00010.1` lies on the reverse strand and its codon 148 covers the deleted bases. Both engines run with `--hgvs`.

### 1.3 What each engine writes

The transcript's row in the default format. Ensembl VEP:

```
mt_del3_290_TCTT_T	Mt:291-293	-	ATMG00010	ATMG00010.1	Transcript	inframe_deletion	442-444	442-444	148	K/-	AAG/-	-	IMPACT=MODERATE;STRAND=-1;HGVSc=ATMG00010.1:c.444_446del;HGVSp=ATMG00010.1:p.Lys149del;HGVS_OFFSET=-2
```

vep-rs:

```
mt_del3_290_TCTT_T	Mt:291-293	-	ATMG00010	ATMG00010.1	Transcript	inframe_deletion	442-444	442-444	148	K/-	AAG/-	-	IMPACT=MODERATE;STRAND=-1;HGVSc=ATMG00010.1:c.446_448del;HGVSp=ATMG00010.1:p.Lys149del;HGVS_OFFSET=-4
```

The consequence, the positions, the amino acids, the codons and the HGVSp agree. The HGVSc and `HGVS_OFFSET` differ by the two bases Ensembl VEP did not shift.

### 1.4 Why Ensembl VEP writes that

`_genomic_shift` expands the variant's slice by 1,000 bases on each side, then clips the expanded slice to the sequence region. For this record the clipped slice runs from `Mt:1` to `Mt:1293`, 1,293 bases, in place of the 2,003 the expansion asked for. The routine then takes `pre_seq` as the first 1,000 bases of the slice's sequence and `post_seq` as its last 1,000, offsets that assume the variant sits 1,000 bases from each end.

Here `post_seq` is right, `Mt:294-1293`, but `pre_seq` is `Mt:1-1000`, which ends 707 bases past the deletion instead of at the base before it. The walk towards the transcript's 3' end, lower coordinates on this strand, compares the deleted bases, last base first, with the end of `pre_seq`. That end reads `AGATT` over `Mt:996-1000`: the `T` and `T` match, the `C` against `A` does not, and the shift stops at two. Walked against the true flank, `Mt:1-290`, the same test matches four times.

Source (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 427-437: the expansion, the clipping and the two fixed cuts.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, line 462: the walk over the two sequences.

### 1.5 Why vep-rs's row is the correct one

The sequence decides. Over `Mt:286-296`, `CTCTTCTTGCA`, deleting any three consecutive bases starting at 287, 288, 289, 290 or 291 leaves the same sequence, `CTCTTGCA`; deleting from 286 does not. The 3' end of this reverse-strand transcript lies towards lower coordinates, so the 3'-most placement is 287-289, four bases from the written 291-293: `HGVS_OFFSET=-4`, and the notation stands two coding bases further along the transcript than Ensembl VEP's, `c.446_448` against `c.444_446`. Ensembl VEP's placement at 289-291 is one of the equivalent placements but not the 3'-most, which is what the notation requires.

### 1.6 Try it yourself

Annotate the record of 1.2 with `--hgvs` against the TAIR10 cache and FASTA. vep-rs's command is below; Ensembl VEP takes the same record with its own Ensembl Plants cache and the FASTA under `--offline --cache`:

```
vep -i record.vcf -o out.txt --offline --json_cache <cache> --fasta <fasta> \
  --species arabidopsis_thaliana --assembly TAIR10 --hgvs
```

| Field | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `HGVSc` | `ATMG00010.1:c.444_446del` | `ATMG00010.1:c.446_448del` |
| `HGVS_OFFSET` | `-2` | `-4` |

The bases can be read from the FASTA: `samtools faidx <fasta> Mt:286-296` prints `CTCTTCTTGCA`, and deleting any three consecutive bases starting at 287, 288, 289, 290 or 291 leaves `CTCTTGCA` each time. A deletion of the same bases placed more than 1,000 bases into a sequence region is shifted the same way by both engines.

### 1.7 How the comparator treats this class

The golden harness (`crates/vep-cli/tests/common/mod.rs`) and the format comparator (`scripts/concordance/compare_vep_formats.py`) apply `hgvs_shift_sequence_edge_window` from a corpus manifest. A `field_divergences` entry of the class names the row and the field (`HGVSc`, `HGVSp` or `HGVS_OFFSET`) with the value vep-rs prints, and the comparison checks vep-rs's value against the documented one in place of Ensembl VEP's. No predicate decides the class from the two outputs, since neither output states the sequence region's length.

No committed golden corpus carries a row of the class. Where one holds an insertion or deletion within 1,000 bases of a sequence edge, its rows with a notation lie on the strand whose walk reads the unclipped flank, or its flag sets write no HGVS. The population cases of the comparison oracle, the harness that runs both engines over the same inputs, carry no such record under `--hgvs` either. The mitochondrial records of the ClinVar GRCh37 and GRCh38 files, 3,122 on each, carry four: under `--hgvs` they give 20 rows of the class on each assembly, `HGVS_OFFSET` on all 20 and `HGVSc` on five.
