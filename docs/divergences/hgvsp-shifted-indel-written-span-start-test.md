# The HGVSp of an indel shifted 3' into the first codon

How the two engines write the protein notation of an insertion or deletion that the 3' shift moves into a transcript's first codon.

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 reads the frameshift, the deleted residues and the peptides of a shifted insertion or deletion on the allele moved 3', and makes its start codon test on the allele as the record writes it. An allele can pass Ensembl VEP's coding test as written, have no start codon test there, and reach the first codon only at its 3' place. Ensembl VEP writes such an allele as a frameshift at residue 1 with the unknown count, or as a deletion of residues from the first. vep-rs makes the start codon test on the moved allele, and writes the notation it writes for the same allele given at its 3' place.

**Why this is a defect.** A protein notation describes the allele that the HGVSc beside it describes. HGVS describes an insertion or deletion at its most 3' position, and both engines move the allele there before they write a notation. The row's HGVSc describes the change at that place, and its `HGVS_OFFSET` gives the number of bases the shift moved it. The protein notation is the predicted consequence of the change that the HGVSc describes. On these rows Ensembl VEP builds it from two readings: the residues of the moved allele, and the start codon verdict of the allele as written.

For the same allele given at its 3' place Ensembl VEP writes another cell, its lost-start form. It writes that form for the three records of this page under `--shift_3prime 1` as well, a flag that does not change where a notation is written. So one allele gets two protein notations by where its record writes it.

The verdict that the notation reads is not always the result of a test. A deletion that begins before the first base of a transcript has no cDNA start as written, so Ensembl VEP's start codon test returns before it is made, although the written deletion itself removes the start codon. An insertion written at the last base before the coding start does not reach the start codon's cDNA positions as written, and lies inside the first codon one or two bases on.

The sequence supports neither of Ensembl VEP's notations. Where the whole coding sequence stands after an insertion, every codon is the reference's, and no reading frame is shifted at residue 1. Where a deletion removes the start codon and no start codon stands after it, a frameshift at residue 1 and a deletion of the first residues each describe a translated protein, and the edited transcript holds no start codon to begin one.

**What vep-rs does instead.** vep-rs reads the start codon on the moved allele, where both engines read the frameshift and the residues. Where the whole coding sequence stands after the allele, the first residue is unchanged and vep-rs writes `p.Met1=`. That is the notation [the page of `start_cooccurrence_swap`](defect-b.md#38-the-protein-notation-under---hgvs) gives a retained start whose coding sequence stands whole. Where no start codon stands after the allele, vep-rs writes the lost start over the residues of the moved span: `p.Met1_?3` for a deletion of a frameshift length that ends in residue 3, and `p.Met1?` where the span ends inside the first codon.

For a deletion of an in-frame length that does not lie on codon boundaries at its 3' place, vep-rs names every residue that Ensembl VEP's deletion spans, then `1_?` and the last residue's number: `p.MetValPro1_?3`. Each of these lost-start forms is the cell vep-rs writes for the same deletion given at its 3' place. Ensembl VEP writes it there too, except for a deletion of whole codons (section 4).

vep-rs moves the start codon test only for an allele that has none where it is written: an allele that overlaps the start codon as written keeps that test's verdict in both engines. Nothing else in the row changes. The consequence terms and the position columns stay those of the allele as written, as in Ensembl VEP's row. An allele written wholly inside the 5' UTR fails the coding test, which reads the span as written (`TranscriptVariationAllele.pm`, line 1731), and gets no HGVSp from either engine. The default, tab and VCF outputs write the `=` of a notation as `%3D` in both engines.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1709-1711, `hgvs_protein`: the allele is shifted 3' for the notation alone. Line 1796 deletes the shift from the allele once the peptides of the moved allele are read.
- The same file, line 2053, `_get_hgvs_protein_type`: `frameshift` reads the coding positions the shift left behind, those of the moved allele.
- The same file, line 2161, `_get_hgvs_peptides`: `start_lost` is asked after the shift is deleted, and where it holds the notation becomes the lost-start form.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 855-862, `start_lost`: the verdict is served from the allele's cache, where the row's consequence terms left the verdict of the allele as written.
- The same file, lines 977-981, `_overlaps_start_codon`: the test reads the unshifted cDNA span plus the shift the allele carries, which line 1796 has deleted. Line 981 returns 0 where the cDNA start is undefined, before any test.
- `TranscriptVariationAllele.pm`, line 2003 with line 2483: a frameshift's count to the new stop is asked from the residue before the first changed one. At residue 1 that is residue 0, for which `_stop_loss_extra_AA` returns nothing: the unknown count.

**What correcting it would change.** With the start codon test made on the moved allele, Ensembl VEP would write on these rows the cell it writes for the same allele given at its 3' place. For a deletion that is vep-rs's cell. For an insertion that leaves the coding sequence whole it is `p.Met1?`, beside which vep-rs writes `p.Met1=`, as the page of `start_cooccurrence_swap` describes. No other column would change.

**The class it produces.** `hgvsp_shifted_indel_written_span_start_test` sets aside the `HGVSp` cell of a pair of rows, Ensembl VEP's and vep-rs's, on an insertion's or deletion's row. Both rows carry one equal, non-zero `HGVS_OFFSET`. Their `Consequence` cells are equal, hold `5_prime_UTR_variant` and hold neither `start_lost` nor `start_retained_variant`: the allele as written has no start codon test. Both cells are protein notations of one protein identifier, and the pair is one of three:

- Ensembl VEP's cell is a frameshift at residue 1 with the unknown count, such as `p.Met1AsnfsTer?`, and vep-rs's is the same first residue unchanged, `p.Met1=` or `p.Met1%3D`. The row's CDS position is empty: the allele as written lies before the coding start.
- Ensembl VEP's cell is a frameshift at residue 1 with the unknown count, and vep-rs's is the lost start named at the same residue: `p.Met1?`, or `p.Met1_?` and a residue number of 2 or more. The row's `Consequence` holds `coding_sequence_variant`, and its CDS position reads `?-` and a number, since the written span has no CDS start. Where vep-rs's HGVSc is a deletion, the residue number is that of the deletion's last coding base, and the cell is `p.Met1?` where that residue is 1.
- Ensembl VEP's cell is a deletion of residues from the first, such as `p.Met1_Pro3del`. vep-rs's cell names as many residues as that deletion spans, the first and the last equal to its two, then `1_?` and the last residue's number. The row is as in the second pair.

Every other pair stays a difference: an empty cell on either side, another residue number, a count in Ensembl VEP's frameshift, a cell of vep-rs's in another form, or a row whose terms hold a start term. The three records of this page annotated under `--shift_3prime 1` carry a start term in both engines' rows, so none of those rows is of the class.

The class is built into the format comparator, `tools/vep-formats-compare`, which applies it once per row and reads the `=` of a notation in either spelling. The concordance figures score consequence sets and never see the field. On the rows of sections 2 and 3 the HGVSc differs as well, under [`hgvs_shifted_deletion_clamped_span`](hgvs-shifted-deletion-clamped-span.md).

**How often the shape occurs.** Tier 2, `--hgvs` being needed. Under `--hgvs`, the genome-wide ClinVar files of both assemblies hold 95 rows of the class against the Ensembl and RefSeq caches of release 116. Of these, 83 are insertions written at the last base before the coding start: 38 and 22 on GRCh38 and 9 and 14 on GRCh37, the Ensembl cache first. Deletions of a frameshift length that begin before a transcript's first base give 5 rows, 2 on GRCh38 and 3 on GRCh37. Deletions of an in-frame length there give 7 rows, 3 on GRCh38 and 4 on GRCh37; one of the 3 removes whole codons. All twelve deletion rows stand on the Ensembl caches.

The committed corpus `tests/golden/116/BDGP6.54-fly` carries one row of the second kind. Its record at MT:9832 deletes four bases across the first base of the reverse-strand transcript FBtr0100880. The stored Ensembl VEP 116.3 row reads `p.Met1LeufsTer?` where vep-rs writes `p.Met1_?2`.

Four worked examples follow, one for each kind of record and one for a deletion of whole codons. Each record is made up, and stands on the cache and the reference of the committed corpus `tests/golden/116/GRCh38-custom`.

## 1. An insertion at the last base before the coding start of ENST00000400183

**In short.** On GRCh38, the record at chromosome 21, position 15730393, has a one-base reference allele and a two-base alternate that adds an `A` between 15730393 and 15730394. Position 15730394 is the first coding base of the forward-strand transcript ENST00000400183, base 413 of its cDNA, and the added base equals it. The 3' rule therefore carries the insertion one base on, between coding bases 1 and 2. Under `--hgvs` each engine's row holds the cells of the allele as written and the HGVSc of the moved one:

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `5_prime_UTR_variant` | `5_prime_UTR_variant` |
| `cDNA_position` | `412-413` | `412-413` |
| `HGVSc` | `ENST00000400183.7:c.1dup` | `ENST00000400183.7:c.1dup` |
| `HGVSp` | `ENSP00000383044.2:p.Met1AsnfsTer?` | `ENSP00000383044.2:p.Met1%3D` |
| `HGVS_OFFSET` | `1` | `1` |

A record at 15730394 with a one-base reference allele and a two-base alternate gives the same insertion at its 3' place. Both engines write `HGVSc=ENST00000400183.7:c.1dup` for it, at CDS position `1-2`. Ensembl VEP annotates it `frameshift_variant,start_lost,start_retained_variant` with `HGVSp=ENSP00000383044.2:p.Met1?`. vep-rs annotates it `frameshift_variant,start_retained_variant` with `HGVSp=ENSP00000383044.2:p.Met1%3D`, the cell it writes for the record as written. The two terms and the two notations of that row are the subject of the page of `start_cooccurrence_swap`.

The three reference bases after the insertion point are the start codon. In the edited sequence the added base stands before them, and the reference bases from the first coding base on follow it unbroken. The 5' UTR is one base longer and every codon of the coding sequence stands, so translation from the start codon gives the reference protein. No residue changes, and no frame is shifted.

As written, the insertion lies between cDNA bases 412 and 413 and does not reach bases 413 to 415, the start codon. The test of lines 977-981 answers 0 for the row's terms, and the answer is cached. Line 2053 reads a frameshift between coding bases 1 and 2, and line 2161 reads the cached answer beside it.

## 2. A deletion of a frameshift length across the first base of ENST00000652508

**In short.** On GRCh38, the record at chromosome 21, position 46326286, has an eight-base reference allele and a one-base alternate. It deletes the seven bases from 46326287 to 46326293: the base before the first base of the forward-strand transcript ENST00000652508, which has no 5' UTR, and its coding bases 1 to 6. The base after them equals the first of them, so the 3' rule carries the deletion one base on, onto coding bases 1 to 7. Under `--hgvs` the two rows read:

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `coding_sequence_variant,5_prime_UTR_variant` | the same |
| `cDNA_position` | `?-6` | `?-6` |
| `HGVSc` | `ENST00000652508.1:c.2_7del` | `ENST00000652508.1:c.1_7del` |
| `HGVSp` | `ENSP00000498568.1:p.Met1ProfsTer?` | `ENSP00000498568.1:p.Met1_?3` |
| `HGVS_OFFSET` | `1` | `1` |

A record at 46326287 with an eight-base reference allele and a one-base alternate gives the same deletion at its 3' place. Both engines annotate it `frameshift_variant,start_lost` at cDNA position `1-7`, with `HGVSc=ENST00000652508.1:c.1_7del` and `HGVSp=ENSP00000498568.1:p.Met1_?3`. That is the cell vep-rs writes for the record as written.

The deleted bases hold the whole start codon, and seven is no multiple of three. After the deletion no start codon stands where the transcript opens, and the coding sequence is not left whole. Coding base 7 lies in the codon of residue 3, the last residue the notation names. The term `5_prime_UTR_variant` names the deleted base before the coding start; the transcript itself holds no 5' UTR base.

The written deletion begins before the transcript's first base, so it has no cDNA start. Line 981 returns 0 before any test, and the row's terms carry no `start_lost`. Line 2161 reads that verdict beside the frameshift of the moved deletion, and line 2003 asks the count from residue 0.

## 3. A deletion of an in-frame length across the first base of ENST00000652508

**In short.** On GRCh38, the record at chromosome 21, position 46326286, has a ten-base reference allele and a one-base alternate. It deletes the nine bases from 46326287 to 46326295: the base before the first base of ENST00000652508 and its coding bases 1 to 8. The 3' rule carries the deletion two bases on, onto coding bases 2 to 10. Under `--hgvs` the two rows read:

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `coding_sequence_variant,5_prime_UTR_variant` | the same |
| `cDNA_position` | `?-8` | `?-8` |
| `HGVSc` | `ENST00000652508.1:c.3_10del` | `ENST00000652508.1:c.2_10del` |
| `HGVSp` | `ENSP00000498568.1:p.Met1_Pro3del` | `ENSP00000498568.1:p.MetValPro1_?3` |
| `HGVS_OFFSET` | `2` | `2` |

A record at 46326288 with a ten-base reference allele and a one-base alternate gives the same deletion at its 3' place. Both engines annotate it `start_lost,inframe_deletion` at cDNA position `2-10`, with `HGVSc=ENST00000652508.1:c.2_10del` and `HGVSp=ENSP00000498568.1:p.MetValPro1_?3`. That is the cell vep-rs writes for the record as written.

The moved deletion takes the second and third bases of the start codon and ends on the first base of the codon of residue 4. The transcript then opens on the start codon's first base and the last two bases of that codon. Together they spell the codon of residue 4 and no start codon. Ensembl VEP's cell says that residues 1 to 3 are deleted, so that the protein begins at residue 4, where no start codon stands to begin it.

As in section 2, the written deletion has no cDNA start and line 981 returns 0 before any test. The notation's type is the deletion of the moved allele's residues, and line 2161 leaves it so.

## 4. A deletion of whole codons across the first base of ENST00000652508

**In short.** On GRCh38, the record at chromosome 21, position 46326286, has a seven-base reference allele and a one-base alternate. It deletes the six bases from 46326287 to 46326292: the base before the first base of ENST00000652508 and its coding bases 1 to 5. The 3' rule carries the deletion one base on, onto coding bases 1 to 6, the first two codons. Under `--hgvs` the two rows read:

| Column | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `coding_sequence_variant,5_prime_UTR_variant` | the same |
| `cDNA_position` | `?-5` | `?-5` |
| `HGVSc` | `ENST00000652508.1:c.2_6del` | `ENST00000652508.1:c.1_6del` |
| `HGVSp` | `ENSP00000498568.1:p.Met1_Val2del` | `ENSP00000498568.1:p.MetVal1_?2` |
| `HGVS_OFFSET` | `1` | `1` |

For the same deletion given at its 3' place the engines differ as well. Ensembl VEP annotates it `inframe_deletion` with `HGVSp=ENSP00000498568.1:p.Met1_Val2del`, and vep-rs `start_lost,inframe_deletion` with `HGVSp=ENSP00000498568.1:p.MetVal1_?2`. That difference is the entry [`start_codon_removed_without_five_prime_utr`](start-codon-removed-without-five-prime-utr.md): Ensembl VEP's start codon test of a removal of whole codons needs a 5' UTR.

vep-rs writes for the record as written the cell it writes for the deletion at its 3' place, and the class takes the pair as its third kind. After the deletion the transcript opens on the codon of residue 3, which is no start codon.

The reproduction `scripts/divergences/reproductions/hgvsp_shifted_indel_written_span_start_test.toml` holds the three records as written and runs them on both engines on GRCh38 under `--hgvs`. It records the consequence, `HGVSc`, `HGVSp` and `HGVS_OFFSET` of each engine's rows on ENST00000400183 and ENST00000652508. For vep-rs it names the JSON cache and the reference of the committed corpus `tests/golden/116/GRCh38-custom`.
