# The HGVSp of an indel moved onto the start codon

How the two engines write the protein notation of an insertion or deletion that `--shift_3prime 1` moves onto a transcript's start codon.

**The defect in one sentence.** Under `--shift_3prime 1`, Ensembl VEP 116.3 makes the second start codon test of `start_lost` at a place that lies the shift length past the allele. Its output step reads the allele's cDNA span at the shifted place, and the test adds the shift length to that span again. An insertion or deletion that the shift moves onto the start codon, and that leaves the whole coding sequence standing, can pass that test on the bases further on.

Its row then carries `start_retained_variant` and no `start_lost`, and its HGVSp is read from the coding sequence edited without the 5' UTR. That is a frameshift at residue 1 with the unknown count, or an insertion between a residue 0 and residue 1. vep-rs tests the start codon at the moved allele and writes the first residue unchanged, `p.Met1=`.

**Why this is a defect.** A protein notation describes the protein of the edited transcript. On these rows the edited transcript is the reference transcript with a shorter or a longer 5' UTR. The start codon stands, formed again by the last bases of the 5' UTR or by the added bases, and every codon follows it. Translation gives the reference protein. Ensembl VEP's own row says so in its term, `start_retained_variant`.

Its notation says otherwise. `p.Met1TrpfsTer?` states that the first residue changes and that the frame shifts from it. `p.Asn0_Asn1insLys` names a residue 0, which no protein has, and gives it and residue 1 the name of the protein's last residue. Both are read from the coding sequence with the allele written into it and no 5' UTR base before it. The edited transcript does not hold that sequence from its start codon on.

The notation is reached only because `start_lost` is not found, and it is not found because the test is made at another place than the allele's. Made at the allele, the test answers that the three bases at the old coding start are no start codon. Ensembl VEP then writes `start_lost` beside `start_retained_variant`, with its lost-start notation. That row is the subject of [the page of `start_cooccurrence_swap`](defect-b.md#38-the-protein-notation-under---hgvs), where vep-rs writes `start_retained_variant` and `p.Met1=`.

So one allele gets three cells from Ensembl VEP, by flag and by where its record writes it. Written in the 5' UTR and annotated under `--hgvs` alone, it gets none: the entry [`hgvsp_shifted_indel_written_span_gate`](hgvsp-shifted-indel-written-span-gate.md). Under `--shift_3prime 1` it gets the frameshift or the insertion at a residue 0. Written at its 3' place it gets the lost start. vep-rs writes `p.Met1=` in all three.

**What vep-rs does instead.** Under `--shift_3prime 1` vep-rs reads every term and every notation on the moved allele, and makes each start codon test at the span of that allele. The notation is the first residue unchanged where two things hold. The edited 5' UTR and coding sequence end in the whole coding sequence, and the first codon does not stand at the old coding start behind an unchanged 5' UTR. An allele that leaves the first codon in place is described from that codon. The unchanged first residue is the cell vep-rs writes for the same allele written at its 3' place without the flag.

The row's terms are those Ensembl VEP writes: `start_retained_variant`, with `frameshift_variant` where the allele's length is no multiple of three. vep-rs writes no `start_lost` beside a retained start, whatever the second test answers, as the page of `start_cooccurrence_swap` describes. The default, tab and VCF outputs write the `=` of a notation as `%3D` in both engines.

**Where it is in Ensembl VEP's source** (`ensembl-vep` and `ensembl-variation` at release/116):

- `Bio/EnsEMBL/VEP/OutputFactory.pm`, lines 297-298 and 306-323: before any row is made, `reset_shifted_positions` calls `clear_shifting_variables` on every allele.
- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 502-527, `clear_shifting_variables`: the cDNA, CDS and translation coordinates of the transcript variation are deleted. Where the allele carries a shift they are read again with the shift length, and serve the moved span from then on.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1028-1075, `_ins_del_start_altered`: line 1047 reads that span, and line 1058 edits the 5' UTR and coding sequence there. The edited sequence ends in the whole coding sequence (line 1071), so `start_retained_variant` holds (lines 947-963).
- The same file, lines 906-945, `_inv_start_altered`: line 920 reads the same span and lines 928-930 add the shift length to it. Line 938 makes the edit there, and lines 939-941 ask whether the three bases at the old coding start are `ATG`.
- The same file, lines 851-904, `start_lost`: line 870 takes the first test's answer and line 871 the second's. Lines 872-882 then compare the peptides. Line 875 returns 0 where the alternate peptide is the one unknown residue `X`, as for a deletion that leaves a partial codon. Lines 878-882 give 0 for an alternate peptide that ends in the reference residue.
- `TranscriptVariationAllele.pm`, line 2161, `_get_hgvs_peptides`: the lost start is written only where `start_lost` holds. Lines 2394-2439, `_get_alternate_cds`, write the allele into the coding sequence alone.
- The same file, lines 2320-2365 with lines 2003 and 2483: a frameshift names the first residue that differs and asks its count from the residue before it. At residue 1 that is residue 0, for which the count routine returns nothing.
- The same file, lines 2231-2238, 2134-2144 and 2383: an inserted peptide that ends in the reference residue is clipped to an insertion that ends at residue 0. Its flanks are read at offset -1 of the translation, which is its last residue.

**What correcting it would change.** With the second test made at the span the output step serves, Ensembl VEP would write `start_lost` beside `start_retained_variant` on these rows, and its lost-start notation as their HGVSp. That is the row of `start_cooccurrence_swap`, beside which vep-rs's row stands as it is here.

**The class it produces.** `hgvsp_start_codon_test_past_the_shifted_allele` sets aside the `HGVSp` cell of a pair of rows, Ensembl VEP's and vep-rs's, on an insertion's or deletion's row. Both rows carry one equal, positive `SHIFT_LENGTH`: the row was annotated at its shifted place. Their `Consequence` cells are equal, hold `start_retained_variant` and hold no `start_lost`, and neither end of their CDS position is unknown. Both cells are protein notations of one protein identifier. vep-rs's is the first residue unchanged, `p.Met1=` or `p.Met1%3D`, and Ensembl VEP's is one of two:

- a frameshift at residue 1 with the unknown count that names the first residue vep-rs's cell names, such as `p.Met1TrpfsTer?`;
- an insertion between a residue 0 and residue 1 whose two flanks carry one name, such as `p.Asn0_Asn1insLys`.

Every other pair stays a difference: an empty cell on either side, a count in Ensembl VEP's frameshift, another residue number, two flank names that differ, a cell of Ensembl VEP's in any other form, or a row whose terms hold `start_lost`. A row with no `SHIFT_LENGTH` field is not of the class. A comparison of two runs made under `--shift_3prime 1` without `--shift_length` therefore reports these pairs.

The class is built into the format comparator, `tools/vep-formats-compare`, which applies it once per row and reads the `=` of a notation in either spelling. The concordance figures score consequence sets and never see the field.

**How often the shape occurs.** Tier 2, `--shift_3prime 1` and a flag that writes the HGVSp being needed. Under `--hgvs --shift_3prime 1 --shift_length`, the chromosome 21 records of gnomAD v4.1 on GRCh38 hold 51 rows of the class against the Ensembl cache of release 116. They stand on two records: a one-base deletion in a 5' UTR gives 46 rows, one a transcript, and a duplication of fifteen bases gives 5. The chromosome 21 records of the 1000 Genomes high-coverage dataset hold none. The allele needs a repeat that carries it from the 5' UTR onto the start codon.

An insertion of a frameshift length that is moved the same way is not of the class. Its alternate peptide is no single `X` and does not end in the reference residue, so Ensembl VEP's peptide test (lines 878-882) finds `start_lost` for it. Its row is the pair of `start_cooccurrence_swap`.

Two worked examples follow, one for each kind of record. Each record is made up, and stands on the cache and the reference of the committed corpus `tests/golden/116/GRCh38-custom`.

## 1. A deletion moved onto coding base 1 of ENST00000481605

**In short.** On GRCh38, the record at chromosome 21, position 37124710, has a two-base reference allele and a one-base alternate. It deletes the base at 37124711, two bases before the first coding base of the forward-strand transcript ENST00000481605, which is 37124713 and base 393 of its cDNA. The deleted base, the base after it and the first coding base are equal, so the 3' rule carries the deletion two bases on, onto coding base 1. Under `--hgvs --shift_3prime 1 --shift_length` the two rows read:

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `frameshift_variant,start_retained_variant` | the same |
| `CDS_position` | `1` | `1` |
| `HGVSc` | `ENST00000481605.6:c.1del` | the same |
| `HGVSp` | `ENSP00000512892.1:p.Met1TrpfsTer?` | `ENSP00000512892.1:p.Met1%3D` |
| `SHIFT_LENGTH` | `2` | `2` |

The edited transcript has a 5' UTR one base shorter. The last base of that UTR and the second and third bases of the old start codon are the start codon, and the coding sequence follows it unbroken. Translation gives the reference protein.

The output step serves the deletion's cDNA span as position 393. `_ins_del_start_altered` removes that base and finds the whole coding sequence at the end of the edited sequence. `_inv_start_altered` adds the shift length once more and removes position 395, the third base of the start codon. The first two bases of the start codon and coding base 4 then stand at the old coding start. Coding base 4 equals the removed base, so the test reads a start codon there. Line 875 ends the peptide test, and line 2161 does not fire.

A record at 37124712 with a two-base reference allele and a one-base alternate gives the same deletion at its 3' place. Under `--hgvs` alone Ensembl VEP annotates it `frameshift_variant,start_lost,start_retained_variant` with `HGVSp=ENSP00000512892.1:p.Met1?`. vep-rs annotates it `frameshift_variant,start_retained_variant` with `HGVSp=ENSP00000512892.1:p.Met1%3D`, the cell it writes under the flag.

## 2. A duplication moved to behind coding base 1

**In short.** On GRCh38, the record at chromosome 21, position 37124710, with a one-base reference allele and a four-base alternate adds three bases after 37124710. They equal the last two bases of the 5' UTR of ENST00000481605 and its first coding base. The 3' rule carries them three bases on, to between coding bases 1 and 2, where they are a duplication that ends on coding base 1. A record of the same form at position 15730391 does the same on the forward-strand transcript ENST00000400183, whose first coding base is 15730394. Under `--hgvs --shift_3prime 1 --shift_length` the rows read:

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `start_retained_variant` | the same |
| `CDS_position` | `1-2` | `1-2` |
| `HGVSc` | `ENST00000481605.6:c.-2_1dup`, `ENST00000400183.7:c.-2_1dup` | the same |
| `HGVSp`, ENST00000481605 | `ENSP00000512892.1:p.Asn0_Asn1insLys` | `ENSP00000512892.1:p.Met1%3D` |
| `HGVSp`, ENST00000400183 | `ENSP00000383044.2:p.Arg0_Arg1insThr` | `ENSP00000383044.2:p.Met1%3D` |
| `SHIFT_LENGTH` | `3` | `3` |

The edited transcript has a 5' UTR three bases longer. The added bases end in the first coding base, so the start codon stands after them and the coding sequence follows it unbroken. Translation gives the reference protein, and no residue is added.

The output step serves the insertion between coding bases 1 and 2. `_ins_del_start_altered` writes the bases there and finds the whole coding sequence at the end of the edited sequence. `_inv_start_altered` writes them three bases further on, between coding bases 4 and 5, and finds the start codon untouched. The alternate peptide ends in the reference residue, so lines 880-881 give 0 and line 2161 does not fire.

The coding sequence with the bases written into it reads one added residue and then the start codon's residue. The clipping of lines 2231-2238 leaves an insertion that ends at residue 0. `Arg` is the last residue of the protein of ENST00000400183. `Asn` is the last residue the cache holds for ENST00000481605, whose coding sequence is cut short at its end.

For the same three bases written at their 3' place, Ensembl VEP writes `start_lost,start_retained_variant` under `--hgvs` alone. Its HGVSp is its lost-start form over a span that ends at residue 0, `ENSP00000512892.1:p.Asn1_?0` and `ENSP00000383044.2:p.Arg1_?0`. vep-rs writes `start_retained_variant` and `p.Met1%3D` there, as under the flag.

The reproduction `scripts/divergences/reproductions/hgvsp_start_codon_test_past_the_shifted_allele.toml` holds the three records and runs them on both engines on GRCh38 under `--hgvs --shift_3prime 1 --shift_length`. It records the consequence, `HGVSc`, `HGVSp` and `SHIFT_LENGTH` of each engine's rows on ENST00000481605 and ENST00000400183. For vep-rs it names the JSON cache and the reference of the committed corpus `tests/golden/116/GRCh38-custom`.
