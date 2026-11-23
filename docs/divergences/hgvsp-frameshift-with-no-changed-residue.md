# The HGVSp of a frameshift that changes no residue

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 names a frameshift at a residue that does not change. Where the alternate protein ends before any residue differs from the reference protein, it writes the last residue it compared, at the number after it, as changed to itself: `p.Leu11LeufsTer?` on a protein of ten residues. vep-rs writes no HGVSp there.

**Why this is a defect.** The HGVS nomenclature's Frameshift page describes a frameshift from "the first amino acid changed": that residue, its position and the residue that replaces it. It gives `p.His150HisfsTer10`, a frameshift that names one residue twice, as a description that "is not correct". On these rows the residue named does not change, the residue written in its place is the same one, and the number is not its own but the one after it. At that number the reference protein holds its stop, another residue, or nothing at all where the coding sequence ends before its stop codon.

**What vep-rs does instead.** vep-rs reads the alternate protein as Ensembl VEP does, over the whole codons of the changed coding sequence and the 3' UTR. Where it ends before a residue differs from the reference protein, the transcript's bases give the frameshift no changed residue to be named by, and vep-rs writes no HGVSp. The row's consequence terms, its positions and its HGVSc follow from the bases, and both engines write them alike.

Two bases can stand after the last whole codon of the changed sequence, and they can spell one residue whatever base comes third. Neither engine reads a residue from an incomplete codon, of a transcript or of a changed sequence. Ensembl VEP translates whole codons, cuts an incomplete last codon off a transcript's protein so that it is not guessed, and writes an incomplete codon as the unknown residue in the row's own `Amino_acids` cell. vep-rs keeps that reading and names no residue from such a pair, as section 4 shows.

Where Ensembl VEP has no description of a change it writes no HGVSp, and neither engine writes a bare `p.?`. vep-rs leaves the field empty in the same way.

A frameshift whose first changed residue the transcript's bases give keeps its notation, on the same transcripts. On the rows of this entry only the HGVSp differs from Ensembl VEP's row.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2320-2365, `_get_fs_peptides`: line 2333 translates the alternate coding sequence, with the 3' UTR where the transcript has one. Line 2348 walks the two translations from the row's residue while the position lies inside the alternate translation. Lines 2351-2352 read the residue of each at the position, line 2360 leaves the walk where the two differ, and line 2361 moves one position on.
- Line 2364 is the defect: it returns the position and the two residues with no test of how the walk ended. Where it ended by its bound, they are the pair read last, two equal residues, and the position after them.
- Lines 2342-2346 test for a position past the alternate translation once, before the walk, and write a deletion there. Such a row is no row of this entry.
- Lines 1898-2044, `_get_hgvs_protein_format`: line 1911 withholds the form of an unchanged residue, `=`, from a frameshift. Line 2013 writes the residue, the position, the second residue, `fsTer` and the count, which line 2006 leaves `?`.
- `Bio/Tools/CodonTable.pm` of the BioPerl release in the Ensembl VEP 116.3 image, line 448: two bases after the last whole codon are translated only under the option `complete_codons`, which line 2333 does not set.
- `Bio/EnsEMBL/Transcript.pm` of that image, lines 2118-2133, `translate`: a transcript's protein is translated with an incomplete last codon cut off, so that it is not guessed. Lines 774-776 of `TranscriptVariationAllele.pm` write the unknown residue for an incomplete codon of the row's own peptides.

**What correcting it would change.** The rows would carry no HGVSp. No other column would change.

**The class it produces.** `hgvsp_frameshift_with_no_changed_residue`: Ensembl VEP's description is a frameshift with the unknown count and one residue name on both sides of its number, vep-rs writes none, and the row's terms hold `frameshift_variant`.

The format comparator, `tools/vep-formats-compare`, tells the form from the row alone. Line 2013 is the one line that writes `fsTer`, and the walk leaves by a changed residue only on two different residues, so nothing else writes one name twice. Two different names, a number for the count in Ensembl VEP's description, a notation in vep-rs and the two engines' values exchanged stay differences. The corpus tool, `vep-golden-corpus classify`, records such a pair of a committed corpus under the same class.

## 1. A deletion that keeps the last residue where no stop codon closes the coding sequence

**In short.** On GRCh38, ClinVar variation 4081419, a record at chromosome 3, position 158646779, deletes one base of ENST00000464732, a forward-strand transcript of GFM1 flagged `cds_end_NF`. Ensembl VEP writes `<protein>:p.Leu61LeufsTer?`; vep-rs writes no HGVSp. `<protein>` stands for the protein's identifier, which both engines write alike.

The coding sequence has 181 bases, 60 whole codons and one base of a 61st, and no 3' UTR follows it. The protein has 60 residues. The record deletes coding base 180, the third base of codon 60. The codon closes on base 181 and still codes leucine, so the changed sequence holds 60 whole codons, the reference's 60 residues. The frame of whatever follows shifts, and the transcript holds no base of it.

Ensembl VEP's walk starts at residue 60, finds a leucine on both sides and moves to 61. The alternate translation has 60 residues, so the walk ends there, and the pair it read at 60 is written under the number 61. The protein has no residue 61.

## 2. A deletion that keeps every residue and breaks the stop codon

**In short.** On GRCh38, ClinVar variation 11411, a record at chromosome X, position 49251335, deletes two bases of FOXP3, a gene on the reverse strand. On ENST00000518685 Ensembl VEP writes `<protein>:p.Pro405ProfsTer?`, and on ENST00000557224 `<protein>:p.Pro457ProfsTer?`. vep-rs writes no HGVSp on either.

No flag marks either transcript incomplete. A stop codon closes each coding sequence, and each transcript ends on its stop codon's last base, with no 3' UTR. On ENST00000518685 the coding sequence has 404 codons and the stop codon. The record deletes coding bases 1212 and 1213, the third base of codon 404 and the first base of the stop codon. Codon 404 closes on the stop codon's second base and still codes proline, and one base is left after it.

Every residue of the protein is kept and the stop codon is broken, so the row's terms hold `stop_lost` beside `frameshift_variant` in both engines. The protein goes on past residue 404 with bases the transcript does not hold. Ensembl VEP names the stop's position, 405, by the proline before it. A lost stop is described by the residue that takes its place, and the transcript holds no codon to read one from.

The same holds on ENST00000557224, whose coding sequence has 456 codons and the stop codon: the record deletes coding bases 1368 and 1369.

## 3. A deletion that removes the next residue's codon

**In short.** On GRCh38, ClinVar variation 3367156, a record at chromosome 16, position 1460840, deletes four bases of ENST00000564568, a reverse-strand transcript of CLCN7 flagged `cds_end_NF`. Ensembl VEP writes `<protein>:p.Ala118AlafsTer?`; vep-rs writes no HGVSp.

The coding sequence has 356 bases, 118 whole codons and two more bases. The protein has 118 residues: residue 117 is an alanine and residue 118 a glycine. The record deletes coding bases 351 to 354, the third base of codon 117 and the whole of codon 118. Codon 117 closes on base 355 and still codes alanine, and one base is left after it.

The changed sequence holds 117 whole codons, the reference's first 117 residues. Ensembl VEP's walk ends at 118 and writes the alanine of residue 117 there, where the reference protein holds a glycine. The codon of residue 118 is deleted and the one base left spells no residue, so what replaces the glycine is not in the transcript.

## 4. Two bases left after the last whole codon

**In short.** On GRCh38, ClinVar variation 2972051, a record at chromosome 16, position 57651588, deletes one base of ENST00000566187, a forward-strand transcript of ADGRG1 flagged `cds_end_NF`. Ensembl VEP writes `<protein>:p.Pro153ProfsTer?`; vep-rs writes no HGVSp.

The coding sequence has 459 bases, 153 whole codons, and the protein has 153 residues: residue 152 is a proline and residue 153 a serine. The record deletes one base of codon 152, whose three bases are equal, so both engines describe it at the codon's last base, coding base 456, and print `HGVS_OFFSET=2`. Codon 152 closes on base 457 and still codes proline, and two bases are left after it.

Those two bases spell a valine codon whatever base comes third, where the reference protein holds a serine. Neither engine reads a residue from an incomplete codon. Ensembl VEP's translation stops at the last whole codon, and its walk ends at 153 with the proline of residue 152. vep-rs reads the alternate protein over whole codons as well and writes no notation. Of the 300 rows counted below, the reference protein holds a residue at the number written on 113, and on 23 of those such a pair spells another residue.

## 5. An immunoglobulin segment whose protein holds stops

**In short.** On GRCh38, the gnomAD v4.1 record at chromosome 14, position 105886031, deletes one base of ENST00000390571, the reverse-strand D segment IGHD3-22, flagged `cds_start_NF` and `cds_end_NF`. Ensembl VEP writes `<protein>:p.Leu11LeufsTer?`; vep-rs writes no HGVSp.

The segment is one exon of 31 bases, all of them coding: 10 whole codons and one base. Its protein has 10 residues, three of them stops, at residues 4 to 6, and is read on past them; the form does not turn on them. The record deletes coding base 30, the third base of codon 10, which closes on base 31 and still codes leucine. The walk ends at 11 with the leucine of residue 10.

Three further records of gnomAD v4.1 give the same form on such segments. One inserts a base at the same codon of this transcript. One deletes a base of ENST00000431870, a protein of five residues, where Ensembl VEP writes `<protein>:p.Leu6LeufsTer?`. One inserts a base into ENST00000390591, a protein of ten, where it writes `<protein>:p.Cys11CysfsTer?`.

## Where real data holds such rows

The form needs a transcript whose bases end soon after the change: a coding sequence with no 3' UTR behind it, or with one too short to hold a whole codon of the new frame. Deletions and insertions both give it. The Ensembl VEP 116.3 outputs of the six whole-genome datasets under `--hgvs` hold it on 300 rows. On each, the record's alleles, the reference sequence and the transcript's exons show that no whole codon of the changed sequence differs from the reference protein. vep-rs writes no HGVSp on all 300, and every other column of those rows equals Ensembl VEP's.

| Dataset | Rows |
|---|---|
| ClinVar, GRCh38 | 56 |
| ClinVar, GRCh37 | 78 |
| 1000 Genomes Phase 3, GRCh37 | 1 |
| gnomAD v2.1.1 genomes, GRCh37 | 36 |
| 1000 Genomes high coverage, GRCh38 | 4 |
| gnomAD v4.1 genomes, GRCh38 | 125 |

The ClinVar rows stand on 55 of 4,398,127 records of the GRCh38 file and on 73 of 4,388,172 of the GRCh37 file. The rows are counted in each dataset, and five variants stand in two.

At the number Ensembl VEP writes, the reference protein holds nothing on 154 rows, where no stop codon closes the coding sequence, as in section 1. It holds its stop on 33, as in section 2, and a residue on 113, as in sections 3 and 4. The transcript has no 3' UTR on 298 rows and a 3' UTR of two bases on 2.

The flags do not decide it. Of the 300 rows, 247 stand on transcripts flagged `cds_end_NF`, with `cds_start_NF` or without it, 7 on transcripts flagged `cds_start_NF` alone and 46 on transcripts that no flag marks incomplete.

None of the 300 rows is on chromosome 21, so no record of [the chromosome 21 cuts](../intended-divergences.md#differences-outside-the-classes) carries the form.

One committed corpus holds a row of the form: line 80 of `tests/golden/116/GRCz11-zebrafish/variants.vcf` deletes coding bases 2544 and 2545 of ENSDART00000155375, a reverse-strand zebrafish transcript flagged `cds_end_NF` whose coding sequence has 849 whole codons and one more base. Codon 848 still codes cysteine and two bases are left after it. Ensembl VEP writes `<protein>:p.Cys849CysfsTer?`, the cysteine of residue 848 at 849, where the reference protein holds a proline. The corpus's manifest records the row under the class in its two flag sets that write HGVSp.

## How to reproduce it

`scripts/divergences/reproductions/hgvsp_frameshift_with_no_changed_residue.toml` runs both engines on the five records above under `--hgvs`. It cuts the five region files that hold the transcripts from the GRCh38 cache and compares the six rows of those records on the six transcripts of sections 1 to 5. [Reproducing the examples](reproducing.md) gives the commands.
