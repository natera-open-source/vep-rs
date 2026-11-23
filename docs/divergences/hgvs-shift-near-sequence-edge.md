# The 3' shift of an indel within 1,000 bases of a sequence region's edge

**The defect in one sentence.** Ensembl VEP reads the flanking sequence it shifts an insertion or deletion along out of a window cut at fixed 1,000-base offsets from the ends of a slice around the variant, already clipped to the sequence region. Within 1,000 bases of a chromosome's or an organelle genome's edge the flank it reads is not the flank beside the variant, and the HGVS shift it reports is not the shift the sequence admits. It stops short of that shift in this page's first example and runs past it in the second, where the misread bases happen to match. vep-rs reads the flank the sequence region holds beside the variant and reports the shift the sequence admits.

**Why this is a defect.** HGVS notation places an insertion or a deletion at the most 3' position of the transcript at which the altered sequence is the same, and both engines shift the allele for that purpose. The shift, found by a walk along the flank one base at a time, is a property of the bases around the variant alone. Ensembl VEP's window logic is right everywhere except where the window would run past the start or end of the sequence region. There the clipping shortens the sequence it reads while the offsets it cuts at stay fixed, and the shift is computed against bases that are not adjacent to the variant.

**What vep-rs does instead.** vep-rs reads the bases on each side of the variant from the reference FASTA, as many as the sequence region holds, and shifts the allele along them until the altered sequence changes. At the edge of a sequence region the available flank is shorter and the shift is computed against it. A deletion that contains a transcript whole before and after the shift keeps the notation of the whole transcript and carries the `HGVS_OFFSET` of the shift, as the allele's other rows with a notation on that strand do.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 427-437, `_genomic_shift`: the slice around the variant is expanded to 1,000 bases on each side and then `constrain_to_seq_region` clips it. `pre_seq` is the first 1,000 bases of the clipped slice and `post_seq` its last 1,000, which are the variant's flanks only when nothing was clipped.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, line 462: `perform_shift` is called with the transcript's strand, and lines 315 and 320 walk `pre_seq` for a reverse-strand transcript and `post_seq` for a forward-strand one, so only the clipped side's flank misleads the walk.
- `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm` (`ensembl-vep` at release/116.3), line 122: the variant takes the slice of the first transcript annotated against it, as the cache stores it, so the clipping ends at the region length the cache records. A FASTA record of another length does not move the edge.

**What correcting it would change.** This class vanishes: the HGVSc, the HGVSp and `HGVS_OFFSET` of an indel near a sequence edge carry the shift the sequence supports, as they do everywhere else. Some of its rows then differ for another reason. Wherever the flank is read whole, Ensembl VEP writes no `HGVSc` for a deletion that contains a transcript whole and admits a shift towards that transcript's 3' end, the documented difference [`hgvs_shifted_indel_written_span_gate`](hgvs-shifted-indel-written-span-gate.md). The rows of such a deletion near an edge, which 1.1 describes, then fall to that class.

**The class it produces.** `hgvs_shift_sequence_edge_window`, documented, where a corpus carries such a record, in its manifest for a row whose `HGVSc`, `HGVSp` or `HGVS_OFFSET` differs for this reason, and applied by the golden harness and the format comparator from the manifest. A corpus is a committed set of records with Ensembl VEP's outputs, and its manifest lists the rows documented to differ (1.7). The concordance figures score consequence sets and never see the fields, so nothing is set aside in them.

Under `--shift_3prime 1 --shift_length` the same walk is printed as `SHIFT_LENGTH`. The deletion `MT:6468-15586` of the human mitochondrial genome ends 983 bases before the region's end and admits a shift of two towards that end. On its 27 forward-strand rows, on GRCh38 and on GRCh37 alike, Ensembl VEP writes `SHIFT_LENGTH=0` and vep-rs `SHIFT_LENGTH=2`; on its nine reverse-strand rows both write 0. The 27 are the rows of the fifteen forward-strand transcripts the deletion contains whole, of the two it overlaps in part and of ten more within 5,000 bases of it, which get an upstream or a downstream row.

The two engines write the same `Consequence` on 26 of the 27 rows. On the row of `ENST00000361624`, one of the two transcripts the deletion overlaps in part, Ensembl VEP writes `stop_retained_variant` and vep-rs `stop_lost`, each beside `3_prime_UTR_variant`, with these flags and without them. That row's terms are another documented difference, [`stop_codon_removed_without_stop_codon_test`](stop-codon-removed-without-stop-codon-test.md), and no part of this one.

On the two transcripts the deletion overlaps in part, `cDNA_position` and `CDS_position` stand two places further on vep-rs's row and `Protein_position` one place further on one of them, `ENST00000361789`; under `--shift_3prime 1` alone the same fields differ on those rows. Each engine annotates the deletion at the span its own walk moved it to, so every field that gives a place on the transcript follows that walk.

The class covers `SHIFT_LENGTH` and these three position fields too. Given the reference FASTA and vep-rs's JSON cache, the format comparator sets aside a `SHIFT_LENGTH` where Ensembl VEP's value is the length of the walk along its clipped slice and vep-rs's the length the sequence admits. It sets aside a `cDNA_position`, a `CDS_position` or a `Protein_position` where each engine's value is the one the transcript gives the span that engine's walk moved the allele to. A position printed with a total under `--total_length` stays a difference. So do the codons, the amino acids and the consequence terms of a row, which this class never takes.

Within 1,000 bases of an end of a region the genome-wide ClinVar file of each assembly holds 14 insertions or deletions, and the two walks differ on eight of them. Under `--shift_3prime 1 --shift_length`, eight of the 14 insertions or deletions that the genome-wide ClinVar file of each assembly holds within 1,000 bases of an end of a region have rows that differ in these fields, all on the mitochondrial genome. They are the four of 1.7, two more within 1,000 bases of the region's end and two within 1,000 bases of its start.

The four others give no row of the class under `--hgvs`: on the strand whose walk is misled they have downstream rows alone, which carry no notation and still print `SHIFT_LENGTH`. Near the start Ensembl VEP's walk on the reverse strand runs one base for one record and two for the other, where the sequence admits none. On each assembly the eight give 69 rows that differ in `SHIFT_LENGTH`; of those, five differ in `cDNA_position`, two in `CDS_position` and one in `Protein_position`. The comparator sets all of them aside.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A deletion 291 bases into the Arabidopsis mitochondrial genome

**In short.** A three-base deletion at `Mt:291-293` of the Arabidopsis thaliana TAIR10 assembly lies in a repeat that admits a shift of four bases towards the transcript's 3' end. Ensembl VEP shifts it two bases, `c.444_446del` with `HGVS_OFFSET=-2`; vep-rs shifts it four, `c.446_448del` with `HGVS_OFFSET=-4`. The same deletion anywhere further than 1,000 bases from the sequence's start is shifted alike by both engines.

### 1.1 What a reader needs to know first

For HGVS, an insertion or deletion inside a repeated sequence is written at the 3'-most position along the transcript at which deleting or inserting the bases gives the same result. `HGVS_OFFSET` reports how many bases the notation's position moved from the record's. On a reverse-strand transcript the 3' direction is towards lower genomic coordinates and the offset is negative; on a forward-strand transcript it is positive.

A sequence region is one FASTA record: a chromosome, a plastid or a mitochondrial genome. A variant within 1,000 bases of its first or last base has less than 1,000 bases of flank on that side.

The 3' walk reads one flank: the downstream flank for a forward-strand transcript and the upstream flank for a reverse-strand one (`perform_shift`, line 320 for the first and line 315 for the second, by the strand it is given at line 462). The shift is therefore wrong for a reverse-strand transcript near a region's start and for a forward-strand transcript near its end, and right for the other two combinations, whose clipped flank the walk never reads. Walked against the wrong bases it can also run further than the true flank admits, so the class covers a shift that differs for this reason in either direction.

The forward-strand case shows on the human mitochondrial genome, 16,569 bases long. A deletion of the one base at `MT:15934`, 635 bases before the region's end, lies in the forward-strand transcript `ENST00000387460`, on the first base of a run of four identical bases, so the sequence admits a shift of three. Ensembl VEP's clipped slice runs from `MT:14934` to `MT:16569`, and `post_seq`, its last 1,000 bases, starts at `MT:15570` where the flank starts at `MT:15935`. The first comparison fails: Ensembl VEP writes `n.47del` with no `HGVS_OFFSET`, vep-rs `n.50del` with `HGVS_OFFSET=3`.

The defect also reaches a deletion that contains a transcript whole. Where Ensembl VEP's walk stops at zero, both engines write the same `HGVSc`, the deletion of the transcript's first to last base, and the rows differ in `HGVS_OFFSET` alone. The deletion `MT:6468-15586` ends 983 bases before the region's end and admits a shift of two: vep-rs writes `HGVS_OFFSET=2` on the fifteen forward-strand transcripts it contains whole, and Ensembl VEP, whose walk stops at zero, writes none.

The zero offset is also why Ensembl VEP writes a notation on those rows at all. Where a deletion runs past a transcript's last base, Ensembl VEP holds the end of the notation at that base, a clamped end (line 2733 of `TranscriptVariationAllele.pm`). The test at line 1472 of the same file refuses a notation whose clamped end the offset would carry past that base, and only a zero offset passes it.

Where the flank is read whole, as for a deletion more than 1,000 bases from an edge, the walk finds the shift the sequence admits, and where that shift is not zero the test refuses a notation with a clamped end. Ensembl VEP then writes no `HGVSc` on the rows of a transcript the deletion contains whole, a difference of another class, [`hgvs_shifted_indel_written_span_gate`](hgvs-shifted-indel-written-span-gate.md).

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

The HGVSp names residue 149 where `Protein_position` reads 148. Residues 148 and 149 are both lysine, so the protein that loses either is the same, and the notation names the last of the two, the HGVS rule for a run of equal residues. Without `--shift_3prime 1`, `Protein_position` keeps the record's place, as the other position columns do. The HGVSp is equal on the two rows because every placement of the deletion that 1.5 lists leaves the same coding sequence.

### 1.4 Why Ensembl VEP writes that

`_genomic_shift` expands the variant's slice by 1,000 bases on each side, then clips the expanded slice to the sequence region. For this record the clipped slice runs from `Mt:1` to `Mt:1293`, 1,293 bases, in place of the 2,003 the expansion asked for. The routine then takes `pre_seq` as the first 1,000 bases of the slice's sequence and `post_seq` as its last 1,000, offsets that assume the variant sits 1,000 bases from each end.

Here `post_seq` is right, `Mt:294-1293`, but `pre_seq` is `Mt:1-1000`, which ends 707 bases past the deletion instead of at the base before it. The walk towards the transcript's 3' end, lower coordinates on this strand, compares the deleted bases, last base first, with the end of `pre_seq`. After each match it turns the deleted bases by one, last base to the front, so a walk can run further than their length. That end reads `AGATT` over `Mt:996-1000`: the `T` and `T` match, the `C` against `A` does not, and the shift stops at two. Walked against the true flank, `Mt:1-290`, the same test matches four times.

Source (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 427-437: the expansion, the clipping and the two fixed cuts.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 306-349, the loop of `perform_shift`, which line 462 calls: the walk. For a reverse-strand transcript it reads `pre_seq` from its last base backwards, at line 315, and takes a step for each match.

### 1.5 Why vep-rs's row is the correct one

The sequence decides. Over `Mt:286-296`, `CTCTTCTTGCA`, deleting any three consecutive bases starting at 287, 288, 289, 290 or 291 leaves the same sequence, `CTCTTGCA`; deleting from 286 does not. The 3' end of this reverse-strand transcript lies towards lower coordinates, so the 3'-most placement is 287-289, four bases from the record's 291-293: `HGVS_OFFSET=-4`, and the notation stands two coding bases further along the transcript than Ensembl VEP's, `c.446_448` against `c.444_446`. Ensembl VEP's placement at 289-291 is one of the equivalent placements but not the 3'-most, which is what the notation requires.

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

The golden harness (`crates/vep-cli/tests/common/mod.rs`) and the format comparator (`tools/vep-formats-compare`) apply `hgvs_shift_sequence_edge_window` from a corpus manifest. A `field_divergences` entry of the class names the row and the field (`HGVSc`, `HGVSp` or `HGVS_OFFSET`) with the value vep-rs prints, and the comparison checks vep-rs's value against the documented one in place of Ensembl VEP's. No predicate decides the class from the two outputs alone, since neither output states the sequence region's length; beside the reference FASTA and vep-rs's JSON cache the format comparator decides it (2.7).

Two committed golden corpora carry rows of the class, `ARS-UCD2.0-cattle` and `GRCz11-zebrafish` (2.7). Where another holds an insertion or deletion within 1,000 bases of a sequence edge, its rows with a notation lie on the strand whose walk reads the unclipped flank, or its flag sets, each the flags of one committed run, write no HGVS.

The comparison oracle is the harness that runs both engines over the same inputs, and its population cases are the ones that annotate the datasets the concordance figures are measured on. The chromosome 21 files they annotate under `--hgvs` hold no record within 1,000 bases of that chromosome's ends. The population cases that annotate the genome-wide ClinVar files run under `--everything`, which turns `--hgvs` on: they read the four records of the next paragraph, and its counts hold for them.

The mitochondrial records of the ClinVar GRCh37 and GRCh38 files, 3,122 on each, hold four insertions or deletions with rows of the class: the two human deletions of 1.1 (`MT:15934` and `MT:6468-15586`), an insertion between `MT:15939` and `MT:15940` and a deletion of `MT:15940`. Under `--hgvs` the four give 20 rows of the class on each assembly, which differ in `HGVS_OFFSET` on all 20 and in `HGVSc` on five. No corpus manifest documents these rows: the format comparator decides all 20 beside the reference FASTA and vep-rs's JSON cache (2.7).

The deletion `MT:6468-15586` gives 17 of the 20: the fifteen transcripts it contains whole, where `HGVS_OFFSET` alone differs, and the two it overlaps in part. Its ten forward-strand upstream and downstream rows carry neither a notation nor an offset in either engine. Each of the other three records gives one row, on the transcript `ENST00000387460`. `HGVSc` differs on those three rows and on the rows of the two transcripts the long deletion overlaps in part, five in all.

## 2. An insertion 107 bases before the end of a zebrafish scaffold

**In short.** A two-base insertion between `KN150614.1:8912` and `KN150614.1:8913` of the zebrafish GRCz11 assembly admits no shift: the base after it is not the first inserted base. Ensembl VEP shifts it one base, `c.333_334insAG` with `HGVS_OFFSET=1`, and derives the protein change of the insertion it moved, `p.Asn112ArgfsTer4`. vep-rs leaves it where the record puts it, `c.332_333insGA` with no offset and `p.Asn112ThrfsTer4`. Here the walk runs past the shift the sequence admits, and the HGVSp differs with the HGVSc.

### 2.1 What a reader needs to know first

A step of the 3' walk is admitted where the base after an insertion equals the first inserted base. The insertion can then be written one base further with its bases turned by one, and the changed sequence is the same. Where that base differs, a placement one base further describes another sequence.

A frameshift's HGVSp names the first residue that changes, the residue that takes its place and the count of residues to the new stop, the changed residue counted as the first and the stop as the last. All three are read from the changed coding sequence. An insertion written one base from its place, where the sequence admits no step, can therefore name another residue.

The scaffold `KN150614.1` is 9,019 bases long. The forward-strand transcript `ENSDART00000186300` ends at the scaffold's base 9,018. Its coding sequence of 438 bases runs to the transcript's last base, with no 3' UTR after it. A forward-strand transcript's walk reads the flank after the variant (1.1), the side this scaffold's end clips.

### 2.2 The record

The record is line 364 of the input of the `GRCz11-zebrafish` golden corpus (`tests/golden/116/GRCz11-zebrafish/variants.vcf`), annotated against that corpus's cache of Ensembl release 116 and its reference FASTA:

```
KN150614.1	8912	frameshift_insertion_12	T	TGA	.	.	.
```

The inserted bases are `GA`, between bases 8912 and 8913, the coding bases 332 and 333 of the transcript. The bases `KN150614.1:8911-8916` read `GTCAAC`, and the bases `KN150614.1:8020-8021` read `GG`. Both engines run with `--hgvs --sift b`, the corpus's flag set `hgvs-sift`.

### 2.3 What each engine writes

The transcript's row in the default format. Ensembl VEP:

```
frameshift_insertion_12	KN150614.1:8912-8913	GA	ENSDARG00000111916	ENSDART00000186300	Transcript	frameshift_variant	332-333	332-333	111	V/VX	gtc/gtGAc	-	IMPACT=HIGH;STRAND=1;HGVSc=ENSDART00000186300.1:c.333_334insAG;HGVSp=ENSDARP00000147465.1:p.Asn112ArgfsTer4;HGVS_OFFSET=1
```

vep-rs:

```
frameshift_insertion_12	KN150614.1:8912-8913	GA	ENSDARG00000111916	ENSDART00000186300	Transcript	frameshift_variant	332-333	332-333	111	V/VX	gtc/gtGAc	-	IMPACT=HIGH;STRAND=1;HGVSc=ENSDART00000186300.1:c.332_333insGA;HGVSp=ENSDARP00000147465.1:p.Asn112ThrfsTer4
```

The consequence, the positions, the amino acids and the codons agree. The HGVSc, the HGVSp and `HGVS_OFFSET` differ.

### 2.4 Why Ensembl VEP writes that

Ensembl VEP numbers an insertion's start one past its end, here 8913 and 8912. `_genomic_shift` asks for the slice `KN150614.1:7913-9912` and clips it to the scaffold, `7913-9019`: 1,107 bases where the expansion asked for 2,000. `pre_seq`, the first 1,000 bases, is `7913-8912` and ends at the base before the insertion, as it should. `post_seq`, the last 1,000, is `8020-9019`. The flank that follows the insertion starts at 8913, 893 bases further on.

The walk compares the first inserted base, `G`, with the first base of `post_seq`. Base 8020 is `G`, so the walk takes a step and turns the inserted bases to `AG`. Its next comparison, `A` against base 8021, `G`, fails, and the shift stops at one. Walked against the true flank, the first comparison is `G` against base 8913, `C`, and no step is taken.

The HGVSc is written at the moved place with the turned bases. For the HGVSp, Ensembl VEP writes the turned bases into the coding sequence at the moved coordinates, after coding base 333. It translates the result and names the first residue that differs from the protein.

Source (`ensembl-variation` at release/116), all in `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`:

- lines 427-437: the expansion, the clipping and the two fixed cuts.
- lines 319-344, inside `perform_shift`, which line 462 calls: the forward-strand walk along `post_seq`, a step and a turn of the allele for each match.
- lines 1774-1778 and 1330-1348: with a shift, the protein notation is built from the turned allele.
- lines 844-858 and 2404-2421: the changed coding sequence takes the allele at the shifted coordinates.
- lines 2320-2365 and 1995-2015: the first residue that differs, and the count to the new stop.

### 2.5 Why vep-rs's row is the correct one

The sequence decides. Base 8913 is `C`, not `G`, so the insertion written one base further gives another sequence. Over `8911-8916` the record's insertion leaves `GTGACAAC`, and the moved one `GTCAGAAC`. HGVS takes the 3'-most of the placements that give the same sequence, and this insertion has one placement: between coding bases 332 and 333, `c.332_333insGA`, with no offset.

The protein follows. Codon 112 of the transcript reads `AAC`, asparagine. The record's insertion stands after the second base of codon 111, which keeps its residue, and makes codon 112 `ACA`, threonine, and the changed reading frame stops at codon 115, the fourth counted from 112: `p.Asn112ThrfsTer4`. Ensembl VEP's moved insertion makes codon 112 `AGA`, arginine. Its `p.Asn112ArgfsTer4` is the protein change of an insertion the record does not hold.

### 2.6 Try it yourself

Annotate the record of 2.2 with `--hgvs`, which writes the three fields, against the corpus's cache and its `reference.fa.gz`. vep-rs's command is below; Ensembl VEP takes the same record with its own cache of the species and the FASTA under `--offline --cache`:

```
vep -i record.vcf -o out.txt --offline --hgvs \
  --json_cache tests/golden/116/GRCz11-zebrafish/json_cache \
  --fasta tests/golden/116/GRCz11-zebrafish/reference.fa.gz \
  --species danio_rerio --assembly GRCz11
```

| Field | Ensembl VEP | vep-rs |
| --- | --- | --- |
| `HGVSc` | `ENSDART00000186300.1:c.333_334insAG` | `ENSDART00000186300.1:c.332_333insGA` |
| `HGVSp` | `ENSDARP00000147465.1:p.Asn112ArgfsTer4` | `ENSDARP00000147465.1:p.Asn112ThrfsTer4` |
| `HGVS_OFFSET` | `1` | none |

The bases can be read from the FASTA: `samtools faidx <fasta> KN150614.1:8913-8913` prints `C`, and `KN150614.1:8020-8021` prints `GG`, the two bases Ensembl VEP's walk reads. The corpus's `expected_hgvs-sift/default.txt.gz` holds Ensembl VEP's row as its run with `--hgvs --sift b` wrote it.

### 2.7 How the comparator treats this class

The `GRCz11-zebrafish` corpus documents the three fields of this row under the class in each of its flag sets `everything` and `hgvs-sift`. The `ARS-UCD2.0-cattle` corpus documents one `HGVSc` and two `HGVS_OFFSET` values under it in each of the same two flag sets.

Beside the reference FASTA and vep-rs's JSON cache (`--reference-fasta`, `--json-cache`), the format comparator decides the class without a manifest. It repeats both walks on the reference, and accepts a row's HGVSc and `HGVS_OFFSET` where each engine's value is the one its walk leaves on the transcript. It accepts a frameshift's HGVSp where each engine's value is the frameshift its walk gives on the cache's coding sequence. An HGVSp that differs in another way stays a difference unless a manifest documents it. `vep-golden-corpus classify --reference-fasta` asks the comparator and writes the entries.

Under `--shift_3prime 1` the comparator accepts a `SHIFT_LENGTH` and a position field in the same way, each engine's value the one its walk gives.
