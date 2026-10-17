# Defect E: a breakend's mate-side row evaluated at the record's own coordinate

**The defect in one sentence.** For a breakend, Ensembl VEP chooses which transcripts get a mate-side row by the mate's chromosome and coordinate, then evaluates almost every consequence test on that row with the record's own coordinate, which may lie on another chromosome.

**Why this is a defect.** The row is written for the mate breakpoint: its transcript was chosen because the mate lands in or near it. Asking whether the record's chromosome 1 position lies inside a chromosome 21 transcript's intron answers a question about the wrong breakpoint. Only `feature_truncation` looks at the mate, so the row comes out as bare `feature_truncation` or, when even that does not fire, falls through to `intergenic_variant` on a row that names a transcript.

**What vep-rs does instead.** vep-rs evaluates the mate breakpoint where it is: inside the transcript it writes `feature_truncation` plus the exon, intron, UTR or coding region the breakpoint lands in; within 5,000 bases of the transcript, `upstream_gene_variant` or `downstream_gene_variant`.

**Where it is in Ensembl VEP's source** (`ensembl-variation`):
- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 257 and 273: the record's own feature is passed to every predicate.
- `Utils/VariationEffect.pm`, lines 350-358: `feature_truncation` alone tests `within_feature` on the breakend; lines 444-458: `upstream` and `downstream` read the record; lines 629-634: `within_intron` reads the record's differing region, which `BaseTranscriptVariationAllele.pm` line 118 derives from the record's start.
- `StructuralVariationOverlap.pm`, lines 86-102 and 144-160: the mate allele is created from the mate's coordinate and chromosome.

**What correcting it would change.** If every predicate were evaluated at the breakpoint the row is for, [class 11](#11-breakend-mate-rows-read-from-the-other-chromosomes-coordinate) would vanish.

**The classes it produces.**
- [Class 11](#11-breakend-mate-rows-read-from-the-other-chromosomes-coordinate): mate-side rows of point breakends whose mate lies on another chromosome or far from the record (the synthetic `09_breakends` and `13_vcf45_features` sets and gnomAD-SV v4.1 chr21: 1,407 and 10,518 pairs). Set aside: both members of each pair.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 11. Breakend mate rows read from the other chromosome's coordinate

**In short.** The synthetic breakend `synth_bnd_inter_0066_B` sits on chromosome 1 and names its mate at 21:33034355, inside intron 1 of SOD1 and 2,542 bases upstream of a lincRNA. Ensembl VEP picks the transcripts for the mate-allele rows by the mate's coordinate but evaluates the rows with the record's own chromosome 1 coordinate, so it writes `feature_truncation` alone for SOD1 and `intergenic_variant` for the lincRNA. vep-rs writes `feature_truncation,intron_variant` and `upstream_gene_variant` with `DISTANCE=2542`. vep-rs's rows state what the mate breakpoint does where it is, which is what Ensembl VEP itself writes for the same breakpoint when it is a record's own coordinate, so they are the correct ones.

### 11.1 What a reader needs to know first

A breakend record has two coordinates: its own (POS on CHROM) and its mate's (the chromosome and position inside the allele's brackets).

Ensembl VEP represents a breakend in one of two ways:

- A `<BND>` whose mate is on the same chromosome and given in `INFO/CHR2` and `END2` becomes a span from POS+1 to the mate ([class 5](defect-c.md#5-a-record-ensembl-vep-drops-before-annotation) shows one). Every allele of a span record is evaluated over the span.
- A breakend written in bracket notation, or a `<BND>` whose mate is on another chromosome, becomes a point record whose Location is the single base after POS.

This class is about the mate-allele rows of point records. For a transcript within 5,000 bases of the mate, Ensembl VEP writes a mate-allele row (`StructuralVariationOverlap.pm` lines 86-102 and 144-160). What that row should say is what the mate breakpoint does to that transcript:

- when the breakpoint falls inside it: `feature_truncation` (SO:0001906), plus the region the breakpoint lands in (`intron_variant`, SO:0001627, or an exon term);
- when it falls near it: `upstream_gene_variant` (SO:0001631, "a sequence variant located 5' of a gene") and `downstream_gene_variant` (SO:0001632, "a sequence variant located 3' of a gene") with a DISTANCE.

`intergenic_variant` (SO:0001628) is the row for a position near no transcript. On a row that names a transcript it is a contradiction in terms.

### 11.2 The record

**An inter-chromosomal breakend** from the synthetic breakend set `tests/sv_validation/grch37/09_breakends.vcf.gz` (GRCh37; the GRCh38 set `tests/sv_validation/grch38/09_breakends.vcf.gz` has the same construction at other coordinates). It is a pair of records that are each other's mates:

```
1	14374333	synth_bnd_inter_0066_B	A	]21:33034355]A	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0066_A
21	33034355	synth_bnd_inter_0066_A	G	]1:14374333]G	.	PASS	SVTYPE=BND;MATEID=synth_bnd_inter_0066_B
```

The rows in question are the first record's mate-allele rows: the transcripts near 21:33034355.

- ENST00000270142 (SOD1, forward strand, canonical, GRCh37 33031935-33041244) contains the mate breakpoint in its intron 1 (33032155-33036102).
- ENST00000449339 (a lincRNA on the reverse strand at 33025906-33031813) ends 2,542 bases before it, which for a reverse-strand transcript is upstream.

The record's own coordinate, 1:14374334, is on another chromosome altogether.

**The coincidence** `gnomAD-SV_v2.1_BND_21_56716`, gnomAD-SV v2.1 chr21 (GRCh37, `gnomad_sv_v2.1_chr21.vcf.gz`):

```
21	15437970	gnomAD-SV_v2.1_BND_21_56716	N	<BND>	892	UNRESOLVED;UNSTABLE_AF_PCRMINUS;LOW_CALL_RATE	END=15437971;SVTYPE=BND;SVLEN=-1;CHR2=Y;POS2=15437970;END2=15437971
```

Its mate, Y:15437971, lies in intron 17 (Y:15436587-15438100) of ENST00000329134 (UTY, reverse strand, Y:15434914-15592550). Its own coordinate, 21:15437971, is the same number on chromosome 21.

### 11.3 What each engine writes

Ensembl VEP, record `synth_bnd_inter_0066_B`, mate allele (2 of 7 rows), and vep-rs, the same rows, one table per transcript, the differing cells in bold. First ENST00000270142 (SOD1):

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `synth_bnd_inter_0066_B` | `synth_bnd_inter_0066_B` |
| Location | `1:14374334` | `1:14374334` |
| Allele | `]21:33034355]A` | `]21:33034355]A` |
| Gene | `ENSG00000142168` | `ENSG00000142168` |
| Feature | `ENST00000270142` | `ENST00000270142` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | **`feature_truncation`** | **`feature_truncation,intron_variant`** |
| cDNA_position | `-` | `-` |
| CDS_position | `-` | `-` |
| Protein_position | `-` | `-` |
| Amino_acids | `-` | `-` |
| Codons | `-` | `-` |
| Existing_variation | `-` | `-` |
| Extra | `IMPACT=HIGH;STRAND=1` | `IMPACT=HIGH;STRAND=1` |

Record `synth_bnd_inter_0066_B`, mate allele, on ENST00000449339 (the lincRNA):

| Column | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Uploaded_variation | `synth_bnd_inter_0066_B` | `synth_bnd_inter_0066_B` |
| Location | `1:14374334` | `1:14374334` |
| Allele | `]21:33034355]A` | `]21:33034355]A` |
| Gene | `ENSG00000234509` | `ENSG00000234509` |
| Feature | `ENST00000449339` | `ENST00000449339` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | **`intergenic_variant`** | **`upstream_gene_variant`** |
| cDNA_position | `-` | `-` |
| CDS_position | `-` | `-` |
| Protein_position | `-` | `-` |
| Amino_acids | `-` | `-` |
| Codons | `-` | `-` |
| Existing_variation | `-` | `-` |
| Extra | **`IMPACT=MODIFIER;STRAND=-1`** | **`IMPACT=MODIFIER;DISTANCE=2542;STRAND=-1`** |

The record has seven mate-allele rows and every one differs the same way:

- The four transcripts that contain the breakpoint get `feature_truncation` alone from Ensembl VEP and `feature_truncation` plus the region from vep-rs.
- The three nearby ones get `intergenic_variant` from Ensembl VEP and `upstream_gene_variant` or `downstream_gene_variant` with a distance from vep-rs.

Ensembl VEP's own rows for the second record of the pair are what vep-rs writes for the first. The second record is the same breakpoint written as the record's own line. Its rows below are under its local allele, the allele for its own breakpoint, printed `.G`.

Ensembl VEP, record `synth_bnd_inter_0066_A`, local allele (2 of 7 rows), one column per transcript:

| Column | Ensembl VEP, ENST00000270142 | Ensembl VEP, ENST00000449339 |
| --- | --- | --- |
| Uploaded_variation | `synth_bnd_inter_0066_A` | `synth_bnd_inter_0066_A` |
| Location | `21:33034356` | `21:33034356` |
| Allele | `.G` | `.G` |
| Gene | `ENSG00000142168` | `ENSG00000234509` |
| Feature | `ENST00000270142` | `ENST00000449339` |
| Feature_type | `Transcript` | `Transcript` |
| Consequence | `feature_truncation,intron_variant` | `upstream_gene_variant` |
| cDNA_position | `-` | `-` |
| CDS_position | `-` | `-` |
| Protein_position | `-` | `-` |
| Amino_acids | `-` | `-` |
| Codons | `-` | `-` |
| Existing_variation | `-` | `-` |
| Extra | `IMPACT=HIGH;STRAND=1` | `IMPACT=MODIFIER;DISTANCE=2543;STRAND=-1` |

The distance differs by one because the two rows are the two sides of one junction. Ensembl VEP places a record's own breakpoint on the first base after POS (33034356 for the second record) and a mate on the base the bracket names (33034355 in the first record's allele). The two are one base apart by convention, in both engines.

On the coincidence record both engines write the same row.

Ensembl VEP and vep-rs, record `gnomAD-SV_v2.1_BND_21_56716`, the one row both write:

| Column | Ensembl VEP and vep-rs |
| --- | --- |
| Uploaded_variation | `gnomAD-SV_v2.1_BND_21_56716` |
| Location | `21:15437971` |
| Allele | `N[Y:15437971[` |
| Gene | `ENSG00000183878` |
| Feature | `ENST00000329134` |
| Feature_type | `Transcript` |
| Consequence | `feature_truncation,intron_variant` |
| cDNA_position | `-` |
| CDS_position | `-` |
| Protein_position | `-` |
| Amino_acids | `-` |
| Codons | `-` |
| Existing_variation | `-` |
| Extra | `IMPACT=HIGH;STRAND=-1` |

### 11.4 Why Ensembl VEP writes that

Ensembl VEP decides which transcripts get a mate-allele row by the mate's coordinate and chromosome. It then evaluates almost every consequence predicate (the test that decides one consequence term) on that row with the record's own coordinate. Only `feature_truncation` looks at the mate breakpoint; `within_intron`, the exon tests, `upstream` and `downstream` all read the record's POS, on chromosome 1, against a transcript on chromosome 21.

For SOD1 the chromosome 1 number 14374334 is nowhere near the transcript's chromosome 21 coordinates, so the intron test finds nothing and only `feature_truncation` fires. For the lincRNA nothing fires at all, and the row falls through to `intergenic_variant`. When the record's own number happens to land inside or near the transcript, as on the coincidence record, the row comes out right.

Source (`ensembl-variation`):

- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 257 and 273: `$bvf` is the record's own feature and is passed to every predicate.
- `Utils/VariationEffect.pm`, lines 350-358: `feature_truncation` alone tests `within_feature` on `$bvfoa->breakend`; lines 444-458: `upstream` and `downstream` read `$bvf`; lines 629-634: `within_intron` reads the record's differing region, which `BaseTranscriptVariationAllele.pm` line 118 derives from `$vf->{start}`.
- `StructuralVariationOverlap.pm`, lines 86-102 and 144-160: the mate allele is created from the mate's coordinate and chromosome.

### 11.5 Why vep-rs's row is the correct one

The mate breakpoint 21:33034355 lies in intron 1 of SOD1 (33032155-33036102) and 2,542 bases upstream of the lincRNA's end at 33031813. Those are facts about the coordinates the record states, and they are what vep-rs writes.

Ensembl VEP's row for SOD1 says the transcript is truncated but names no region. Its row for the lincRNA calls a position 2,542 bases from the transcript intergenic, on a line that names the transcript.

Ensembl VEP itself annotates the same breakpoint correctly when it is the record's own coordinate ([11.3](#113-what-each-engine-writes), the second record of the pair). Its mate row changes when the record's chromosome 1 position is moved while the mate stays fixed ([11.6](#116-try-it-yourself)), which no correct reading of the mate row could do.

### 11.6 Try it yourself

Annotate with Ensembl VEP 116.2, GRCh37 cache and FASTA, using [the command of the closing section](reproducing.md#reproducing-these-examples). Keep the mate at `21:33034355` and move the record's own position:

| Record | Ensembl VEP, mate-allele rows on ENST00000270142 and ENST00000449339 |
| --- | --- |
| `1 14374333 . A ]21:33034355]A` (as in the file) | `feature_truncation`; `intergenic_variant` |
| `1 33034300 . A ]21:33034355]A` | `feature_truncation,intron_variant`; `upstream_gene_variant` with `DISTANCE=2488` |

The second input places the chromosome 1 number inside SOD1's chromosome 21 span and 2,488 bases from the lincRNA, and the rows follow the number. The mate rows describe wherever the record's own coordinate happens to point. vep-rs writes `feature_truncation,intron_variant` and `upstream_gene_variant` with `DISTANCE=2542` for both inputs.

### 11.7 How the comparator treats this class

`filter_breakend_mate_local_read_pairs` in `scripts/validation/compare_sv_concordance.py`, with the derivation in `scripts/validation/breakend_mate_context.py`, is computed, not pattern-matched. For a pair whose row is a mate-allele row of a point breakend ([11.1](#111-what-a-reader-needs-to-know-first)) it derives both of these sets from the transcript cache:

- the set the mate coordinate supports for the transcript: inside it, `feature_truncation` plus the exon, intron, UTR or coding region the breakpoint lands in, with the non-coding and NMD context terms; within 5,000 bases of it, `upstream_gene_variant` or `downstream_gene_variant`;
- the set the record's own coordinate yields by Ensembl VEP's reading: `feature_truncation` alone when the mate is inside, `intergenic_variant` when it is not.

It sets a pair aside only when the two derivations differ, Ensembl VEP's set equals the one from the record's own coordinate and vep-rs's set equals the one the mate coordinate supports. A vep-rs over- or under-annotation matches neither and stays charged.

Counts, reported as `excluded_breakend_mate_local_read_rust` and `excluded_breakend_mate_local_read_perl`:

- GRCh37: 1,407 pairs (`09_breakends.vcf.gz` 1,369; `13_vcf45_features.vcf.gz` 38).
- GRCh38: 10,518 pairs (`09_breakends.vcf.gz` 6,333; `13_vcf45_features.vcf.gz` 27; `gnomad_sv_v4.1_chr21.vcf.gz` 4,158).

Rows like the coincidence record, where the two derivations agree, are concordant and not in the class. The class is fully masked, both sides.

vep-rs releases before 0.3.0 evaluated no position on the mate side and so reproduced Ensembl VEP's fall-through rows; from 0.3.0 vep-rs evaluates the mate breakpoint where it is, and these pairs are what that costs the raw figure, which sets nothing aside.
