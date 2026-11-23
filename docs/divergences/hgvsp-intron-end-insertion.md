# The HGVSp of an insertion at either end of an intron

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 writes a definite `HGVSp` for an insertion that lies in the gap between an intron and a coding exon, as if the inserted bases were bases of the exon; vep-rs writes no `HGVSp` there.

**Why this is a defect.** The gap is a splice junction, and the inserted bases can be placed on either side of it, with the intron or with the exon. Ensembl VEP's `HGVSp` states the outcome of one placement as the consequence: the junction used where it stands, the inserted bases translated as the exon's.

The HGVS Nomenclature (`hgvs-nomenclature.org/stable/`) predicts no protein outcome at the junction. Its Splicing page, `recommendations/RNA/splicing/`, under Discussion, gives the unknown-effect form `p.?` as the protein description for a change of an intron's last two positions. Its protein page for substitutions, `recommendations/protein/substitution/`, gives the same form in its splicing example. The nomenclature makes a protein description optional, so an empty field conforms to it and one definite consequence does not.

A copy of the intron's last base or bases is described by intron positions alone. An insertion that copies none has no sentence of its own in the nomenclature. Either placement of its bases falls in one of the Splicing page's two lists, the changes expected to affect splicing and the changes that might. The step from the second list to `p.?` rests on that one protein example.

After an exon, an insertion that copies no base before the gap is the same case with the sides exchanged. Its bases can go with the exon's end or with the intron's start, and the value that conforms is again the `HGVSc` with no protein description. This page cites no further sentence of the nomenclature for that side.

HGVS describes an insertion at its most 3' position, so for a record that writes the insertion alone the `HGVSc` is the same wherever the record writes the change. Before an exon, Ensembl VEP's `HGVSp` is not. On a transcript with no intron of 13 bases or fewer, it is written for the record that writes the insertion in the gap and not for a record that writes the same change inside the intron. Section 1 shows such a pair.

**What vep-rs does instead.** vep-rs writes no `HGVSp` for an insertion in either gap. Every other column of the row is the one both engines write.

Before an exon the rule has three conditions. The record writes the insertion between the intron's last base and the exon's first coding base. The 3' rule leaves it there. The intron is longer than 13 bases. The inserted bases may be any. The exon may open anywhere in the coding sequence, the stop codon included; the ClinVar rows counted below hold no gap before an exon that opens inside the stop codon. On both engines an unknown reference base matches no inserted base: at the exon's first position it stops the 3' rule, among the intron's last bases it leaves the `HGVSc` an insertion, and vep-rs writes no `HGVSp` in that gap.

After an exon the rule has two conditions. The 3' rule leaves the insertion between the exon's last coding base and the intron's first base, whether the record writes it there or inside the exon. The inserted bases are no copy of the reference bases before the gap. The intron may be of any length. On both engines an unknown reference base matches no inserted base: at the intron's first position it stops the 3' rule, and before the gap the inserted bases are no copy of it, so the `HGVSc` is an insertion and vep-rs writes no `HGVSp`.

A record holds the insertion alone, trimmed to the inserted bases between two positions, or beside another allele, with a base the alleles share. Both engines clip the shared base when they write the `HGVSc`. vep-rs reads the allele as that notation clips it. Where the clipped insertion lies in one of the two gaps, vep-rs asks of the reference sequence whether the 3' rule would leave it there. If it would, the conditions above decide, as they do for the same insertion written alone in that gap. Section 6 shows a record of the second kind, and section 7 the neighbouring shapes that keep their value.

**Where it is in Ensembl VEP's source** (`ensembl-variation` and `ensembl` at release/116, the modules the Ensembl VEP 116.3 image holds; each file lies under `modules/` in its repository). The routines are the same for both gaps:

- `Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm`, lines 454-571, `_bvfo_preds`, flags the allele as the record writes it. Lines 340, 363 and 394 of that file compute the flags once and store them, so the 3' shift of the next bullet leaves them as they are. Lines 462-463 take the allele's span, for an insertion the two bases beside its gap. Lines 511-516 set `exon` when an exon of the transcript overlaps the span (`BaseTranscriptVariation.pm`, lines 863-865, `_overlapped_exons`). Without an exon, lines 535-537 flag the allele `non_coding`. With one, line 540 fetches the CDS coordinates of the span, a list of mapped pieces and gaps. Lines 543-545 set `coding` where that list is one mapped piece, as it is for an insertion in either gap.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, `hgvs_protein`: line 1711 shifts the allele 3'. Lines 1730-1738 then return nothing unless `coding` is set (line 1731) and the shifted position has a first and a last residue (lines 1732-1733).
- Before an exon, two records of one change part at line 1731. The record written in the gap has the exon's first base beside it, so `exon` and `coding` are set and the test passes. The record written inside the intron has two intron bases beside it, is flagged `non_coding`, and `hgvs_protein` returns at line 1737. A transcript with an intron of 13 bases or fewer is the exception (`BaseTranscriptVariation.pm`, lines 1032-1035). There line 861 sets a stretch of 12 bases, by which line 865 widens the span against every exon of the transcript. `_bvfo_preds` then sets `coding`, at its lines 555-557, for a record written inside the intron within that reach, and the test passes for it too.
- After an exon, two records of one change do not part. For a transcript with an Ensembl identifier the shift of line 1711 is made on the genome (`_return_3prime`, lines 109-259, at lines 140-145; `_genomic_shift`, lines 411-466). `perform_shift`, lines 291-352, carries the insertion while the inserted base at its leading end equals the next base of the genome. An insertion whose first base differs from the intron's first base stays in the gap. One written inside the exon has two exon bases beside it, so `coding` is set, and is carried to the gap, where both records pass the test. An insertion carried on into the intron has no residue there, and `hgvs_protein` returns at line 1737.
- `Bio/EnsEMBL/Mapper.pm`, lines 485-580, `map_insert`, reached from `cds_coords` (`BaseTranscriptVariation.pm`, line 513), places the insertion. Lines 489-491 map the two bases beside it as one span, which comes back as a gap for the intron base and a coordinate for the exon base. Lines 504-558 drop the gap and set the insertion against the exon base. In the gap before an exon it is put before that base, at line 541 on a forward-strand transcript and line 525 on a reverse-strand one. In the gap after an exon it is put after that base, at line 527 on a forward-strand transcript and line 539 on a reverse-strand one. No line on the way reads the inserted bases or the intron's.
- `Bio/EnsEMBL/TranscriptMapper.pm`, lines 447-448 and 467-468, `genomic2cds`, number that place in the coding sequence, and lines 528-529, `genomic2pep`, give its residues. The coding sequence holds no intron, so the insertion stands there between two coding bases.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm` then makes the notation. Where the shift moved the record, lines 1774-1778 turn the inserted bases to their order at the gap (`shift_feature_seqs`, lines 1330-1348). Lines 2409-2410, in `_get_alternate_cds`, cut the reference coding sequence at the insertion's coordinate, and line 2421 joins it around the inserted bases. Lines 1802-1805 clip the residues the two peptides share. Lines 2053-2054 set a frameshift, written at lines 2012-2014. Otherwise lines 2063-2065 set an insertion, written at lines 1987-1991, or lines 2072-2078 a replacement, written at lines 1982-1983.

These lines were followed on GRCh37 records that write the insertion alone, in runs without `--shift_3prime 1`. Before an exon they were followed for three forms: a copy of the intron's last base, a copy of two or more bases ending on it, and bases that copy none and do not end as an intron ends. Bases that end as an intron ends were counted and not followed. A record that holds the insertion beside another allele was run on both engines and not followed. A transcript without an Ensembl identifier, whose shift goes on to lines 211-256 of `_return_3prime`, was not followed either.

**What correcting it would change.** The rows would carry an empty `HGVSp`, and no other column would change. Before an exon, that is what Ensembl VEP writes for the same change written inside the intron, on a transcript with no intron of 13 bases or fewer. On a transcript that has one, Ensembl VEP fills the field for the change written inside the intron as well, and vep-rs writes none there. That difference has an entry of its own, [`hgvsp_coding_test_widened_by_a_short_intron`](hgvsp-coding-test-widened-by-a-short-intron.md). Its two records stand 1 and 3 bases before the intron's end; the reach of 12 bases is read from line 861 of `BaseTranscriptVariation.pm`.

**The classes it produces.** Two classes of the format comparator, `tools/vep-formats-compare`. Each sets aside an `HGVSp` that Ensembl VEP fills and vep-rs leaves empty, on a row where both engines write one `HGVSc`, a notation of the row's transcript.

`hgvsp_intron_end_insertion` is the gap before an exon. The `HGVSc` is a duplication of the intron position before an exon's first coding base, a duplication of two or more intron positions ending on that position, or an insertion between that position and the coding base. No `HGVS_OFFSET` stands on either side. The record is written where an intron of more than 13 bases of that transcript ends before a coding exon base.

`hgvsp_intron_start_insertion` is the gap after an exon. The `HGVSc` is an insertion between an exon's last coding base and the intron's first position. The two sides carry no `HGVS_OFFSET`, or one and the same offset, the bases the 3' rule moved the record. The gap the insertion is left in is one where an exon of that transcript ends inside the coding region before an intron of any length. A run under `--shift_genomic 1` prints the moved Location, and the class then takes the gap that Location names.

Both classes read the transcript's exons and coding bounds from the JSON cache of the vep-rs run under comparison, named by `--json-cache`, and set nothing aside without it. A record that holds the insertion with a shared base is placed by the gap of its allele once clipped to the bases the `HGVSc` inserts. Its row in the default or tab format gives the record's own Location, so the class works the gap out from that Location and those bases. A VCF or JSON entry is placed from its line of the input file.

These stay differences: an `HGVSp` on both sides, an `HGVSc` that differs between the engines or has another shape, offsets that differ, and a transcript the cache does not hold. The class of the gap before an exon also sets nothing aside on a row that carries an `HGVS_OFFSET`, or at the end of an intron of 13 bases or fewer. An `HGVSp` on one side alone is reported there as a difference. The tests of both classes are in `tools/vep-formats-compare/tests/intron_ends.rs`. The concordance figures compare consequence sets and do not read this field.

On the rows Ensembl VEP 116.3 writes for the genome-wide ClinVar files under `--hgvs`, vep-rs leaves empty an `HGVSp` that Ensembl VEP fills on these rows:

| The insertion | GRCh37 rows | GRCh38 rows |
| --- | --- | --- |
| Before an exon, a copy of the intron's last base | 34 | 112 |
| Before an exon, a copy of two or more bases ending on the intron's last base | 33 | 162 |
| Before an exon, no copy, not ending as an intron ends | 89 | 410 |
| Before an exon, no copy, ending as an intron ends | 8 | 79 |
| After an exon | 231 | 1,133 |

Of the rows after an exon, 199 on GRCh37 and 986 on GRCh38 are written at the gap, and 32 and 147 are carried to it from inside the exon. The counts are of ClinVar alone; no other dataset was counted.

## 1. A copy of the intron's last base, before an exon of ENST00000312957

**In short.** On GRCh37, a record at 21:27070988 with one reference base and two alternate bases inserts one base between 27070988, the last base of intron 4 of the forward-strand transcript ENST00000312957, and 27070989, the first base of exon 5 and coding base 395. The inserted base equals the base at 27070988, and the 3' rule leaves it in the gap. Both engines write `frameshift_variant,splice_region_variant`, `CDS_position` 394-395, `Protein_position` 132 and one `HGVSc`, a duplication of the intron position before coding base 395, with no `HGVS_OFFSET`. Ensembl VEP writes an `HGVSp`, a frameshift at residue 132; vep-rs writes none.

A record at 21:27070987, with alleles of the same two lengths, writes the same base one position earlier. The two records give one changed sequence, and the 3' rule carries this one 1 base, to the same gap. Both engines write `splice_acceptor_variant`, no CDS or protein position, the same `HGVSc` with `HGVS_OFFSET=1`, and no `HGVSp`. So Ensembl VEP's `HGVSp` depends on which of the two records holds the change.

## 2. A copy of the last six bases of an intron of ENST00000297268

**In short.** On GRCh37, a record at 7:94049852 with one reference base and seven alternate bases inserts six bases between 94049852, the last base of an intron of the forward-strand transcript ENST00000297268, and 94049853, an exon's first base and coding base 2188. The six bases equal the six reference bases that end the intron. Both engines write `inframe_insertion,splice_region_variant` and one `HGVSc`, a duplication of the six intron positions before coding base 2188, with no `HGVS_OFFSET`. Ensembl VEP writes an `HGVSp`, an insertion of two residues between residues 729 and 730; vep-rs writes none.

A record at 7:94049846 writes the same change six positions earlier, and the 3' rule carries it 6 bases, to the same gap. Both engines write `splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant`, the same `HGVSc` with `HGVS_OFFSET=6`, and no `HGVSp`.

## 3. Three bases that copy no intron base, before an exon of ENST00000321219

**In short.** On GRCh37, a record at 21:38877585 with one reference base and four alternate bases inserts three bases between 38877585, the last base of intron 8 of the forward-strand transcript ENST00000321219, and 38877586, the first base of exon 9 and coding base 1240. The three bases are no copy of the bases before them and do not end in the two bases that end an intron. Both engines write `stop_gained,inframe_insertion,splice_region_variant`, `CDS_position` 1239-1240, `Protein_position` 413-414 and one `HGVSc`, an insertion between the intron's last position and coding base 1240, with no `HGVS_OFFSET`. Ensembl VEP writes an `HGVSp`, the insertion of a stop between residues 413 and 414; vep-rs writes none.

This change can be written in one gap alone, so no second record shows it without an `HGVSp`.

On a reverse-strand transcript the positions mirror. On GRCh37, a record at 15:45657058 inserts one base between 45657059, the last base of an intron of ENST00000396659, and 45657058, an exon's first base and coding base 979. Ensembl VEP writes an `HGVSp`, a frameshift at residue 327; vep-rs writes none.

## 4. Bases that end as an intron ends, before an exon of ENST00000353479

**In short.** On GRCh37, a record at 10:105800876 inserts 16 bases into the reverse-strand transcript ENST00000353479. They lie between 105800877, the last base of an intron of 184 bases, and 105800876, the first base of exon 39 and coding base 2648. The 16 bases are no copy of the bases before them, and their last two are the two bases that end an intron. Ensembl VEP writes an `HGVSp`, a frameshift at residue 883; vep-rs writes none.

With this form the rule before an exon has no condition on the inserted bases. Where the intron is longer than 13 bases, every insertion that the record writes between the intron's last base and the exon's first coding base, and that the 3' rule leaves there, has no `HGVSp` on vep-rs's row.

## 5. Three bases after the last base of an exon of ENST00000231790

**In short.** On GRCh37, a record at 3:37045965 with one reference base and four alternate bases inserts three bases after exon 4 of the forward-strand transcript ENST00000231790. They lie between 37045965, the exon's last base and coding base 380, and 37045966, the intron's first base. The three bases are no copy of the reference bases before the gap, and the 3' rule leaves them there. Both engines write `protein_altering_variant,splice_region_variant`, `CDS_position` 380-381, `Protein_position` 127 and one `HGVSc`, an insertion between coding base 380 and the intron's first position, with no `HGVS_OFFSET`. Ensembl VEP writes an `HGVSp`, a replacement of residue 127 by two residues; vep-rs writes none.

A record written inside the exon reaches the same kind of gap. On GRCh37, a record at 6:24528392 inserts three bases between 24528392 and 24528393. Both are bases of an exon of the forward-strand transcript ENST00000348925, and the exon ends at 24528394, coding base 1382. The 3' rule carries the insertion 2 bases, to the gap between 24528394 and 24528395, the intron's first base. Both engines write the `HGVSc` of an insertion in that gap, with `HGVS_OFFSET=2`. Ensembl VEP writes an `HGVSp`, an insertion of one residue between residues 460 and 461; vep-rs writes none.

The rule after an exon has no bound on the intron's length. Of the ClinVar rows after an exon counted above, 9 on GRCh37 and 1 on GRCh38 stand before an intron of 13 bases or fewer. Of the same rows after an exon, 2 on GRCh37 and 7 on GRCh38 follow an exon that ends inside the stop codon.

## 6. The same insertion held beside another allele

**In short.** Give the first record of section 5 a second alternate allele, a substitution of the base at 37045965. The record then has one reference base and two alternates, of four bases and of one, and cannot be trimmed. The insertion's row has the Location 3:37045965, a four-base Allele that opens with the shared base, and `CDS_position` 380. Both engines write the terms, the `Protein_position` and the `HGVSc` of section 5, and that `HGVSc` names the three inserted bases without the shared base. Ensembl VEP writes the `HGVSp` of section 5; vep-rs writes none. On the substitution's row the two engines agree.

Neither engine shifts an allele that a record holds this way, so the row's `HGVSc` is the insertion where the record writes it. ClinVar holds no record of this kind. The rows were counted on 70 GRCh37 records, each a ClinVar insertion at a location where a row of these classes stands, written again with the base before it as shared base and a substitution of that base as second allele.

On 119 of their rows Ensembl VEP writes an `HGVSp` and vep-rs none, beside the same `HGVSc`. On 111 the shared base is an exon's last base, on a forward-strand transcript. On 8 it is an exon's first base, on a reverse-strand transcript, and Ensembl VEP's `HGVSp` is a frameshift on all 8. Where the shared base lies in the intron, neither engine writes an `HGVSp`.

## 7. Neighbouring shapes that keep their value on both engines

**In short.** The rule and the classes stop at the shapes above. The rule changes nothing on the shapes below. On each record named below both engines write the same `HGVSp`, or both write none.

Before an exon:

- An insertion that the 3' rule carries from the gap on into the exon keeps the `HGVSp` it has. So does a substitution or a deletion of the exon's first base, which is described in the coding sequence.
- An insertion at the end of an intron of 13 bases or fewer keeps its `HGVSp`. The comment at lines 858-860 of `BaseTranscriptVariation.pm` says that a variant in such an intron is to be called exonic. On GRCh37, a record at 11:67765163 inserts two bases at the end of an intron of 2 bases of ENST00000227471, and both engines write the same frameshift.
- The rule reads where the record writes the insertion, so it changes nothing on a record written anywhere but the gap, wherever the 3' rule carries the insertion. On GRCh37, a record at 7:100371474 inserts one base inside an exon of ENST00000538115, before an intron of 2 bases. The 3' rule carries it 3 bases, to the gap before the next exon, and both engines write the same frameshift.

After an exon:

- A copy of the reference bases before the gap keeps its `HGVSp` on both engines. A copy of the exon's last bases is a duplication of coding positions. A copy longer than its exon is written by the `HGVSc` as an insertion in the gap, and keeps its `HGVSp` too. On GRCh37, a record at 22:29120964 inserts 391 bases after an exon of ENST00000439346, such a copy, and both engines write the same `HGVSp` on its row.
- An insertion that the record writes alone, and whose first base equals the intron's first base, is carried into the intron by the 3' rule. It has no residue there, and neither engine writes an `HGVSp`.
- Held beside another allele, that insertion gets one `HGVSc` and one `HGVSp` from both engines: an insertion between the exon's last coding base and the intron's first position, with its protein notation. Two of the 70 records of section 6 are such insertions, and they give 12 such rows, on six transcripts. Of those rows this page states only that the engines agree.

## 8. Running the records

**In short.** The entry's reproduction, `scripts/divergences/reproductions/hgvsp_intron_end_insertion.toml`, runs five records on both engines under `--hgvs` and holds six rows to their values. They are both records of section 1, the first record of section 3, the first record of section 5 and the record of section 6 with its two rows. It spells the notations this page gives in words. It cuts three region files from an Ensembl VEP GRCh37 cache and converts them with `tools/vep-storable-to-json`.

The reproduction does not hold the records of sections 2, 4 and 7, the second records of sections 3 and 5, or the 70 records of section 6. Of the records this page names, those at 21:27070987, 7:94049846, 21:38877585, 15:45657058, 10:105800876, 3:37045965 and 6:24528392, and the three of section 7, are records of ClinVar's GRCh37 file. The records at 21:27070988 and 7:94049852 write the first two of those changes at the gap. Strands, exon numbers, coding bases and reference bases are read from the GRCh37 cache and the FASTA that the reproduction names.
