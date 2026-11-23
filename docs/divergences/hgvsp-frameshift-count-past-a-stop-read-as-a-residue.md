# The stop count of a frameshift after a stop codon inside the coding sequence

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 writes the unknown count, `fsTer?`, for a frameshift whose first changed residue lies after a stop codon that the transcript's protein goes on past, also where the new reading frame meets a stop codon. The protein reads such a codon as a residue, such as a selenocysteine codon, or holds it as a stop and names the residues after it, as a polymorphic pseudogene's does. vep-rs writes the count to the new frame's stop where its three bases are not those of any stop codon the protein reads as a residue.

A frameshift whose first changed residue is a stop is written by both engines as an extension named at that stop, `Ter<n><Xaa>extTer<count>`. Where a stop that the protein holds stands before it, Ensembl VEP writes `extTer?` the same way. vep-rs writes the count from the lost stop to the new frame's stop.

**Why this is a defect.** The HGVS nomenclature's Frameshift page counts the position of the stop in the new reading frame from the first changed amino acid, "ending at the first stop codon encountered". It keeps the unknown count for a frame that "does not encounter a new translation termination (stop) codon".

A stop codon before the change is no stop of the variant protein. No residue before the first changed one differs, and the reference protein reads that codon as a residue and goes on past it. On the rows where vep-rs writes a count the new frame meets a stop codon whose bases no sequence edit of the transcript reads as a residue, so the unknown count says what the sequence does not.

A protein that holds a stop and names the residues after it goes on past that stop as well. Ensembl VEP numbers and describes those residues: a substitution after such a stop gets a protein notation, and the frameshift itself is named at its first changed residue. The count speaks of the new reading frame from that residue on, and a stop the protein holds before it is one of the unchanged residues. Where no such stop stands before the change, Ensembl VEP writes the count on the same transcript, as section 3 shows.

The Extension page gives an extension's length as "The number of amino acids added beyond the original Ter". It keeps `extTer?` for a frame in which "no new termination codon is encountered". An extension that names a frameshift at a stop the protein holds speaks of the new frame from that stop on. Where no stop stands before it, Ensembl VEP writes that length itself, as section 4 shows.

**What vep-rs does instead.** Where Ensembl VEP's count is not positive because of such a codon, vep-rs counts from the first changed residue to the first stop codon the new frame meets, through the 3' UTR. It writes that count where the three bases of that stop are not those of any stop codon the protein reads as a residue.

Where they are the same three bases, a sequence edit of the transcript reads a codon of those bases as a residue in the reference frame, and the transcript's model holds no edit for the new frame. The sequence alone does not prove that they end the protein there. No defect is charged to Ensembl VEP on such a row, and vep-rs keeps the unknown count, as Ensembl VEP writes it.

A protein that holds every stop codon of its coding sequence as a stop reads none as a residue. On such a transcript every stop codon of the new frame is a stop, and vep-rs counts to the first.

The stop is sought through the alternate coding sequence and the transcript's 3' UTR, the bases Ensembl VEP's count reads, and no farther. vep-rs keeps the unknown count where the new frame meets no stop there. It keeps it as well where the transcript's peptide differs from the translation of its coding sequence at a codon that is no stop codon, since the peptide is then not the protein of those bases. A first residue that the peptide reads as a methionine over another codon is no such difference. On such a transcript vep-rs writes what Ensembl VEP writes. On the rows that get a count, only the HGVSp changes.

For the extension, where Ensembl VEP's count is not positive, vep-rs counts from the lost stop to the first stop codon the new frame meets. The count is the new stop's residue number less the lost stop's: the residue in the lost stop's place is counted and the new stop is not. The lost stop is a stop the protein holds, as in section 4, or the stop codon that closes the coding sequence behind one. The stop is sought through the alternate coding sequence and the 3' UTR and no farther, and the unknown count stays where the new frame meets none there.

vep-rs counts an extension so only on a protein that holds a stop, reads no stop codon as a residue and is the translation of its coding sequence at every codon. On any other protein such an extension keeps the count Ensembl VEP writes. The test of the new stop's bases, above, belongs to the frameshift's own count and is not made for an extension. A first residue that the peptide reads as a methionine over another codon is no difference from the translation, as for the frameshift's count. This entry makes no claim about the count that stays on another protein.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2477-2526, `_stop_loss_extra_AA`: line 2492 translates the alternate coding sequence and 3' UTR with no sequence edit of the translation applied. Line 2503 takes the first stop of that whole translation, and line 2508 its distance from the residue before the first changed one. Lines 2519-2525 return the count only where it is positive.
- Lines 1898-2044, `_get_hgvs_protein_format`: line 2003 asks for the count of a frameshift, line 2006 writes `?` where none came back, and line 2013 writes `fsTer` and the count. Where the first changed residue is a stop, line 2007 takes the extension, line 2009 takes one off the count, which leaves the new stop out, and line 2010 writes `extTer` and the count.
- Lines 2320-2365, `_get_fs_peptides`: line 2336 reads the reference protein with the sequence edits applied, so the first changed residue is found past a codon the count stops at.

**What correcting it would change.** The rows would carry the count in place of `?`. No other column would change.

**The class it produces.** `hgvsp_frameshift_count_past_a_stop_read_as_a_residue`: Ensembl VEP's description is a frameshift with the unknown count, and vep-rs's is the same frameshift, at the same residue and with the same new residue, with a count.

The format comparator, `tools/vep-formats-compare`, sets such a pair aside only in a comparison that gives it vep-rs's cache with `--json-cache`. It rebuilds the alternate sequence from the transcript's coding sequence and 3' UTR and from the row's own allele and coding position. The pair is set aside where a stop codon that the peptide reads as a residue or holds as a stop stands before the first changed residue, and the rebuilt sequence gives vep-rs's count. A count the sequence does not give stays a difference. So does a count to a stop that has the three bases of a stop codon the peptide reads as a residue.

The class reads a pair of extensions too, on a row whose terms hold `frameshift_variant`: Ensembl VEP's of unknown count and vep-rs's, named at the same stop with the same new residue, with a count. The pair is set aside where the peptide holds a stop and is the translation of its coding sequence at every codon, and a stop stands before the lost stop, which the peptide holds or which closes the coding sequence. The rebuilt sequence must also give vep-rs's count from the lost stop.

## 1. One base inserted after the selenocysteine codon of a selenoprotein

**In short.** On GRCh37, ClinVar variation 4700090, a record at chromosome 19, position 1105661, inserts one base into the coding sequence of ENST00000354171, a forward-strand transcript of the selenoprotein gene GPX4. Ensembl VEP writes `<protein>:p.Gly111ArgfsTer?`; vep-rs writes `<protein>:p.Gly111ArgfsTer3`. `<protein>` stands for the protein's identifier, which both engines write alike.

The coding sequence has 198 codons. It holds two stop codons: the codon of residue 73, which the protein reads as a selenocysteine, and codon 198, which closes it. The record inserts its base between coding bases 329 and 330. Coding base 330 is that same base, so both engines describe the insertion after base 330 and print `HGVS_OFFSET=1`. The first residue that changes is 111, a glycine that becomes an arginine. The new frame meets its first stop codon at residue 113, three residues counted from the first changed one. Its three bases are not those of the codon of residue 73.

Ensembl VEP's translation of the alternate sequence has its first stop at residue 73. The distance from residue 110 is 73 - 110, which is not positive, and the count is left unknown. The same record gives the same two values, at each transcript's own residue number, on ENST00000589115 and ENST00000587648, which read an earlier stop codon as a selenocysteine too. On ENST00000585480, a transcript of the gene whose coding sequence opens after that codon, the record changes residue 22, and both engines write `<protein>:p.Gly22ArgfsTer3`.

## 2. Two bases deleted: a count on one transcript, the unknown count on two

**In short.** On GRCh37, ClinVar variation 871898, a record at chromosome 19, position 1106405, deletes two coding bases of the same gene. On ENST00000589115 Ensembl VEP writes `<protein>:p.Ser162AspfsTer?` and vep-rs `<protein>:p.Ser162AspfsTer28`. On ENST00000354171 both engines write `<protein>:p.Ile170ArgfsTer?`.

| Transcript | First stop of the new frame | Ensembl VEP | vep-rs |
|---|---|---|---|
| ENST00000589115 | residue 189, other bases | `fsTer?` | `fsTer28` |
| ENST00000354171 | residue 220, the same bases | `fsTer?` | `fsTer?` |
| ENST00000587648 | none | `fsTer?` | `fsTer?` |
| ENST00000585480 | residue 120 | `fsTer40` | `fsTer40` |

The four transcripts differ in their exons, so the new frame reads other bases on each. The first three read the codon of an earlier residue as a selenocysteine: residue 73 on the first two and residue 33 on the third. ENST00000589115 has 181 codons: the new frame runs into the 3' UTR and ends on a stop codon of other bases than that codon's. On ENST00000354171, of 198 codons, it also runs into the 3' UTR, and its first stop there has the three bases of the selenocysteine codon. ENST00000587648 is flagged `cds_end_NF`, and its sequence ends before the new frame meets a stop. ENST00000585480 holds no stop codon before the change, and Ensembl VEP's own count is positive.

## 3. Thirty-five bases deleted after the stops a polymorphic pseudogene's protein holds

**In short.** On GRCh37, ClinVar variation 769334, a record at chromosome 7, position 100385561, deletes 35 coding bases of ZAN, whose six transcripts of biotype `polymorphic_pseudogene` each have a protein that holds stops. On ENST00000348028 Ensembl VEP writes `<protein>:p.Ala2344ValfsTer?`; vep-rs writes `<protein>:p.Ala2344ValfsTer367`.

The coding sequence of ENST00000348028 has 2,721 whole codons and two bases over, with no sequence edit. Its protein has 2,721 residues, 20 of them stops, the first at residue 1940, and it is the translation of those codons at every residue. The record deletes coding bases 7030 to 7064 and leaves 2,710 whole codons. The first residue that changes is 2344, an alanine that becomes a valine, after 11 of those stops. From that residue on the new frame meets no stop before residue 2710, the last codon of the alternate coding sequence: 367 residues counted from the first changed one.

Ensembl VEP's translation of the alternate sequence has its first stop at residue 1940. The distance from residue 2343 is 1940 - 2343, which is not positive, and the count is left unknown. On each of the six transcripts the first changed residue is 2344, and the new frame's first stop is the last codon of its alternate coding sequence.

| Transcript | First stop of the new frame | Ensembl VEP | vep-rs |
|---|---|---|---|
| ENST00000348028 | residue 2710 | `fsTer?` | `fsTer367` |
| ENST00000349350 | residue 2801 | `fsTer?` | `fsTer458` |
| ENST00000421100 | residue 2590 | `fsTer?` | `fsTer247` |
| ENST00000427578 | residue 2713 | `fsTer?` | `fsTer370` |
| ENST00000443370 | residue 2613 | `fsTer?` | `fsTer270` |
| ENST00000449052 | residue 2678 | `fsTer?` | `fsTer335` |

A second record of the same file changes a residue before the first of those stops. ClinVar variation 403617, at chromosome 7, position 100371474, inserts one coding base, and the first residue that changes is 1924, a glycine that becomes an arginine. No stop stands before it, so Ensembl VEP's own count is positive: on ENST00000348028 both engines write `<protein>:p.Gly1924ArgfsTer799`. That frame's first stop is residue 2722, the last codon of its alternate coding sequence. It is the codon that closes the deletion's frame at residue 2710, past all 20 stops of the reference frame.

## 4. Two bases deleted before a later stop that the same pseudogene's protein holds

**In short.** On GRCh37, a record of the gnomAD v2.1.1 genomes at chromosome 7, position 100388594, deletes two coding bases of ZAN and makes a stop of its protein the first changed residue. On ENST00000348028 Ensembl VEP writes `<protein>:p.Ter2464SerextTer?`; vep-rs writes `<protein>:p.Ter2464SerextTer257`.

The record deletes coding bases 7388 and 7389 of ENST00000348028. Residue 2463 stays as it is, and the first residue that changes is 2464, the sixteenth of the 20 stops its protein holds, which becomes a serine. Both engines write the change as an extension named at that stop. The new frame meets its first stop at residue 2721, the last codon of the alternate coding sequence. The count is 2721 - 2464, 257: residue 2464 and the 256 after it.

Ensembl VEP's translation of the alternate sequence has its first stop at residue 1940. The distance from residue 2463 is 1940 - 2463, which is not positive, and the count is left unknown. On each of the six transcripts the lost stop is residue 2464, and the new frame's first stop is the last codon of its alternate coding sequence.

| Transcript | First stop of the new frame | Ensembl VEP | vep-rs |
|---|---|---|---|
| ENST00000348028 | residue 2721 | `extTer?` | `extTer257` |
| ENST00000349350 | residue 2812 | `extTer?` | `extTer348` |
| ENST00000421100 | residue 2601 | `extTer?` | `extTer137` |
| ENST00000427578 | residue 2724 | `extTer?` | `extTer260` |
| ENST00000443370 | residue 2624 | `extTer?` | `extTer160` |
| ENST00000449052 | residue 2689 | `extTer?` | `extTer225` |

A record of the same dataset changes the one stop that another pseudogene's protein holds. At chromosome 17, position 79614932, it deletes four coding bases of TSPAN10, and on ENST00000572675 the first residue that changes is 229, that stop, which becomes a cysteine. No stop stands before it, so Ensembl VEP's own count is positive: both engines write `<protein>:p.Ter229CysextTer127`. The new frame's first stop is residue 356, the last codon of its alternate coding sequence, 127 residues from the lost stop.

## Where real data holds such rows

The genome-wide ClinVar files hold frameshifts after such a codon on selenoprotein genes. On GRCh38, 9 of 4,398,127 records carry a row on which vep-rs writes a count where Ensembl VEP writes none, 17 rows in all. On GRCh37, 6 of 4,388,172 records carry 10 such rows. Over both files the counts run from 3 to 118. Each is the count that the record's alleles, the reference sequence and the transcript's exons give, and every other column of those rows equals Ensembl VEP's.

On 58 further rows of GRCh38 and 41 of GRCh37 the new stop has the bases of the codon read as a residue, and on one further row of each assembly the new frame meets no stop. Both engines write the unknown count on all of them.

A protein holds a stop on 91 transcripts of the Ensembl VEP release 116 cache of GRCh37 and on 156 of GRCh38. On GRCh37, 70 are polymorphic pseudogenes and 21 are immunoglobulin or T-cell receptor segments. On GRCh38, 124 are `protein_coding_LoF` transcripts, 27 are such segments, three are LRG transcripts and two are `nonsense_mediated_decay` transcripts. On each, the protein is the translation of the coding sequence at every residue. The RefSeq caches of release 116 hold such proteins too, on 220 transcripts of GRCh37 and 615 of GRCh38.

In the ClinVar files the record of section 3 is the one frameshift after a stop that a protein holds: one further record of the GRCh37 file, with six rows, and none of the GRCh38 file.

The four whole-genome population datasets hold more of them. The table counts the rows on which Ensembl VEP 116.3 writes the unknown count after a stop the protein holds, and how many of them get a count in vep-rs. Each count is the one that the record's alleles, the reference sequence and the transcript's exons give. Of the other rows, four stand on short immunoglobulin segments, where the alternate sequence ends before a residue differs and vep-rs writes [no HGVSp](hgvsp-frameshift-with-no-changed-residue.md). On the rest the new frame meets no stop through the 3' UTR, and both engines write the unknown count.

| Dataset | Rows with `fsTer?` in Ensembl VEP | With a count in vep-rs |
|---|---|---|
| 1000 Genomes Phase 3, GRCh37 | 28 | 28 |
| gnomAD v2.1.1 genomes, GRCh37 | 174 | 159 |
| 1000 Genomes high coverage, GRCh38 | 27 | 25 |
| gnomAD v4.1 genomes, GRCh38 | 936 | 888 |

The same four datasets hold 43 rows on which a frameshift's first changed residue is a stop of such a protein and both engines write the extension; the ClinVar files hold none. On 7 of them both engines write the count: no stop stands before the changed one, and the new frame meets a stop. The table counts the other 36: the rows on which Ensembl VEP 116.3 writes `extTer?`, and the 28 of them that get a count in vep-rs. On each of those 28 a stop the protein holds stands before the lost one. On the other 8, rows of short immunoglobulin segments, the new frame meets no stop through the 3' UTR, and both engines write `extTer?`.

On one of the 28 the lost stop is the stop codon that closes the coding sequence, behind the one stop its protein holds. On GRCh38, a record of the gnomAD v4.1 genomes at chromosome 11, position 104886798, inserts four bases before the closing stop of ENST00000446862. Ensembl VEP writes `<protein>:p.Ter278HisextTer?`; vep-rs writes `<protein>:p.Ter278HisextTer7`.

| Dataset | Rows with `extTer?` in Ensembl VEP | With a count in vep-rs |
|---|---|---|
| 1000 Genomes Phase 3, GRCh37 | 0 | 0 |
| gnomAD v2.1.1 genomes, GRCh37 | 8 | 7 |
| 1000 Genomes high coverage, GRCh38 | 3 | 3 |
| gnomAD v4.1 genomes, GRCh38 | 25 | 18 |

No chromosome 21 transcript of the Ensembl VEP release 116 caches of GRCh37 and GRCh38 reads a stop codon as a residue, and none has a protein that holds a stop, in the RefSeq caches either. So no record of [the chromosome 21 cuts](../intended-divergences.md#differences-outside-the-classes) carries the shape.

## How to reproduce it

`scripts/divergences/reproductions/hgvsp_frameshift_count_past_a_stop_read_as_a_residue.toml` runs both engines under `--hgvs` on the six GRCh37 records of sections 1 to 4; the GRCh38 record at a closing stop is not among them. It cuts the three region files that hold the genes from the GRCh37 cache. It compares the eight rows of the first two records on the four transcripts of the second example and the twelve rows of the next two on the six transcripts of the third. Of the two records of the fourth example it compares the six rows of the first on those six transcripts and the one row of the second. [Reproducing the examples](reproducing.md) gives the commands.
