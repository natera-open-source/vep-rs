# The HGVSp of a change at the stop codon

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 describes a change at the stop codon, and an insertion before it that brings a stop to the stop's residue, from the peptides of the allele's codon window. A change that keeps the whole protein is written as a deletion-insertion of the stop, as an insertion before the stop or an earlier residue, as the stop unchanged where the stop is lost, or not at all. A change that reaches from before the stop codon through it is written as a deletion-insertion of the range ending at the stop. vep-rs writes the first as the stop unchanged or an extension named at the stop, and the second as a frameshift.

**Why this is a defect.** The HGVS nomenclature derives a protein description "from comparing the variant protein sequence with the reference protein sequence" and adds that knowledge of the change on the DNA level should not be used. A change that keeps every residue and adds residues after the last one extends "the reference amino acid sequence at the N- or C-terminal end", which is the nomenclature's definition of an extension. Its Insertion, Deletion and Deletion-insertion pages each send a change that extends the full-length protein at its C-terminal end to the Extension page. A deletion-insertion that replaces a residue by the same residue is, in the nomenclature's words, one "which effectively is no change".

The Frameshift page describes deletions "starting 5' of and including the translation termination (stop) codon" as a frameshift, and the Deletion-Insertion page says the same of deletion-insertions. A frameshift is, in the Frameshift page's words, a change where translation "shifts to another reading frame".

**What vep-rs does instead.** vep-rs compares the alternate protein with the reference protein at the stop codon, and for an insertion before the last residue's codon that brings a stop codon of its own to the stop's residue. Where the two are equal and a stop stands at the same residue, it writes the stop at its residue number as unchanged. Where the whole reference protein is kept and residues follow it, it writes an extension named at the stop. The extension carries the stop's residue number, the first new residue and the count of residues added up to the new stop, or `?` where the transcript holds none. Only the HGVSp changes.

For an insertion inside the last residue's codon that leaves the protein as it was, vep-rs writes the stop as unchanged at any length (section 10).

Where a deletion or deletion-insertion reaches from before the stop codon into it, and runs to its last base or is not a whole number of codons, vep-rs writes a frameshift. The change is read less the bases its two alleles share at either end. It names the first changed residue, the residue in its place and the position of the new stop counted from it, or `?` where the transcript holds none.

Where a transcript's peptide goes on past a stop codon inside the coding sequence, as a residue or a held stop, Ensembl VEP leaves a lost stop's extension count unknown or short. It takes one and the peptide's length from the first stop's place in the whole alternate translation. vep-rs counts from the lost stop to the next stop (sections 7 and 11).

Where an in-frame deletion takes the stop codon and keeps every residue, Ensembl VEP writes the stop as deleted, then the extension's count, with no residue in the stop's place. vep-rs writes the extension named at the stop with the residue that stands in the stop's place and the same count.

On a coding sequence that no stop codon closes, vep-rs writes no HGVSp where Ensembl VEP's description names the stop alone, as unchanged, deleted, replaced by one residue or extended, at a residue where the transcript's peptide carries no stop. At an incomplete first codon it writes none where Ensembl VEP writes the stop as unchanged at residue 1. Both are a difference of another class, [`hgvsp_stop_named_on_incomplete_coding_sequence`](hgvsp-stop-named-on-incomplete-coding-sequence.md). Every other description on such rows keeps the codon window's value, as vep-rs reads that window.

A lost stop with no whole codon after it keeps the codon window's value. So do an allele at an exon's edge, where the spliced sequence is not known, and a stop written as unchanged on a row whose own terms do not hold `stop_lost`. A change into the stop codon after which the protein keeps its length keeps that value too.

A change of a whole number of codons that ends before the stop codon's last base keeps the codon window's value where it changes or removes a residue before the stop. A deletion of that reach that keeps every residue and loses the stop is the in-frame deletion above, written as the extension. The frameshift is written only where a residue before the stop changes.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1898-2044, `_get_hgvs_protein_format`: one chain of tests over a type and the reference and alternate peptides of the codon window. Line 1911 writes the stop as unchanged when the two peptides are equal. Line 1916 writes an extension only when `stop_lost` holds and the type is a deletion or a substitution. Lines 1958-1993 write a deletion-insertion or an insertion, and line 1961 cuts every residue after the alternate's first stop.
- Lines 2047-2111, `_get_hgvs_protein_type`: lines 2072-2077 type two peptides of unequal length, neither of them empty, a deletion-insertion. The stop alone against a longer alternate peptide is so written as a deletion-insertion of the stop. A range that ends at the stop takes its type at the same lines.
- Lines 2063-2064 type a codon window an insertion where its reference peptide is empty once the residues both peptides share are clipped. Where inserted bases read as the last residue and a stop, line 1961 keeps both, and line 1991 writes them as inserted between the two residues that flank them. Lines 1808-1821 call the type, the peptides and the format; they compare no proteins unless the type is a frameshift.
- Lines 2188-2295, `_clip_alleles`: lines 2231-2242 chop the stop that ends both peptides. What is left of an alternate that ends in a re-formed stop is typed an insertion at lines 2271-2285, and lines 1985-1992 write it between the last residue and the stop.
- Lines 692-787, `peptide`: lines 774-776 give the one or two bases left of the stop codon the unknown-residue symbol, and line 2060 rewrites the reference's stop to the same symbol. The two peptides are then equal, and line 1911 writes the stop as unchanged before line 1916 is asked.
- Lines 1916-1919 take the extension's new residue from the alternate peptide, which lines 2167-2169 of `_get_hgvs_peptides` set to the word for a deletion: the stop is written as deleted before the count. The count is that of `_stop_loss_extra_AA`, lines 2477-2526: the place of the first stop in the translated alternate sequence, less one, less the length of the reference peptide. It is unknown when no stop is met, the figure is not above zero, or the lost stop is residue 1 (line 2483).

**What correcting it would change.** The rows would carry the stop written as unchanged, the extension with its first new residue and its count, or the frameshift. No other column would change: the consequence terms, the coding notation, the protein position, the amino acids and the codons are the same in both engines' rows.

**The class it produces.** `hgvsp_stop_codon_protein_comparison`: eleven shapes, all but the tenth of two `HGVSp` values that name one protein. In the first five both name the stop at one residue number above 1. In the first, Ensembl VEP's value is the stop replaced by a stop and vep-rs's the stop unchanged.

In the next four vep-rs's value is an extension named at the stop. Ensembl VEP's is then the stop alone replaced by the first residues that follow the protein, with or without a stop after them that the protein does not hold. Or it is the stop replaced by every added residue and the stop the inserted bases bring. Or it is residues inserted between the last residue and the stop, or the stop unchanged on a row whose terms hold `stop_lost`.

In the sixth shape Ensembl VEP's value is a range from a residue to the stop replaced by residues, none of them a stop, on a row whose terms hold `stop_lost`. vep-rs's is a frameshift named at that residue with the first of them. Its count is unknown, or above their number and not the one that leaves the new stop at the stop's own residue.

In the seventh shape both values are an extension at the same stop with the same first residue, on a row whose terms hold `stop_lost`: Ensembl VEP's with an unknown count, vep-rs's with the count of residues added.

In the eighth shape Ensembl VEP's value is the stop written as deleted before an extension's count, on a row whose terms hold `stop_lost`, and vep-rs's the extension at the same residue with the same count and the first new residue.

In the ninth shape Ensembl VEP's value is residues and a stop inserted between two neighbouring residues. It names one residue, the second neighbour's, on a row whose terms hold `stop_retained_variant`, or several, the protein's from that neighbour to its last. vep-rs's value is the stop unchanged at the number after the protein's last residue's. Section 9 shows both, the second in a table.

The tenth shape, Ensembl VEP's empty value against the stop unchanged, is read as section 10 lists. The eleventh, an extension on a transcript whose peptide holds a stop, is read from vep-rs's cache as section 11 lists.

The class reads the two values and the row alone, but for the cache of the ninth and eleventh shapes; a row whose `FLAGS` name `cds_end_NF` is a row of the eleventh alone. The terms it reads are those of Ensembl VEP's row.

It holds the extension's first residue to the first residue Ensembl VEP names, and its count to the number of residues Ensembl VEP names wherever that value names them all. It cannot tell from the row that the whole reference protein is kept, nor check a count where Ensembl VEP names only the first residues. The same values with the sides exchanged, an empty value on Ensembl VEP's side off the tenth shape, and the stop written as deleted with no extension's count after it stay differences.

For the sixth shape the class cannot tell from the row that the change runs to the stop codon's last base or shifts the reading frame. A range that ends at the stop against a frameshift at another residue or with another first residue stays a difference, and so does one on a row whose terms hold no `stop_lost`.

For the ninth shape's one residue the class reads the alternate protein from Ensembl VEP's value: the reference protein up to the second neighbour's residue, and then a stop. That value does not number the reference protein's stop. The class takes the number from vep-rs's value and holds it to the residue after the second neighbour.

Several residues before the inserted stop are read from vep-rs's cache; without it the pair stays a difference. So do an inserted residue other than the second neighbour's, an insertion between two residues that are not neighbours, the same values with the sides exchanged, and one residue on a row whose terms hold no `stop_retained_variant`.

The class reads the ninth shape's pair inside prediction parentheses, as `--hgvsp_use_prediction` writes both values, the way it reads it outside them. One value alone inside parentheses stays a difference. So does any value of Ensembl VEP's where vep-rs writes no `HGVSp`: the ninth shape's, and the stop unchanged.

## 1. Three bases inserted in the stop codon that leave a stop at its residue

**In short.** On GRCh37, a record at chromosome 1, position 90470803, with a one-base reference allele and a four-base alternate, inserts three bases after coding base 209 of ENST00000361911, a forward-strand transcript whose stop codon is coding bases 211 to 213, residue 71. The three bases copy the three reference bases that follow, so the 3' rule carries the insertion three bases on, to between the stop codon's second and third base. The first inserted base closes a stop codon with the stop codon's first two bases. Both engines write `inframe_insertion,stop_retained_variant`; Ensembl VEP writes `HGVSp=<the stop at residue 71 replaced by a stop, as a deletion-insertion>` and vep-rs `HGVSp=<the stop at residue 71, unchanged>`.

The alternate protein is the reference protein, residue for residue, and its stop stands at residue 71. A substitution that turns one stop codon into another is written the same way by both engines.

## 2. One base inserted in the stop codon, and a new stop in the 3' UTR

**In short.** On GRCh37, a record at chromosome 1, position 17345377, with a one-base reference allele and a two-base alternate, inserts one base after the first base of the stop codon of ENST00000375499, a reverse-strand transcript whose stop is residue 281. The stop codon's place is read as a residue, and the new frame runs through the 3' UTR to a stop codon 22 residues on. Both engines write `stop_lost`; Ensembl VEP writes `HGVSp=<the stop at residue 281 replaced by the first added residue and a stop, as a deletion-insertion>` and vep-rs `HGVSp=<an extension at residue 281: the first new residue, 22 residues added>`.

Ensembl VEP's value names one of the 22 added residues and a stop the protein does not hold there. On a transcript whose 3' UTR holds no stop codon in the new frame, such as ENST00000421986 for the record at chromosome 1, position 151318416, vep-rs writes the extension with `?` for the count.

## 3. Inserted bases that hold the new stop

**In short.** On GRCh38, a record at chromosome 3, position 81591146, with a one-base reference allele and a 42-base alternate, inserts 41 bases after the second base of the stop codon of ENST00001019268, a reverse-strand transcript whose stop is residue 377. The inserted bases bring ten residues and then a stop codon of their own. Both engines write `stop_lost,inframe_insertion`; Ensembl VEP writes `HGVSp=<the stop at residue 377 replaced by the ten added residues and a stop, as a deletion-insertion>` and vep-rs `HGVSp=<an extension at residue 377: the first new residue, 10 residues added>`.

The protein is the whole reference protein and ten residues more. The nomenclature has a second sentence that fits a stop codon encoded in an inserted sequence, by which the description would be an insertion of the ten residues and a stop between the last residue and the stop. No sentence ranks the two, and the Extension page calls a change that extends the protein at its C-terminal end "by definition an extension".

## 4. Residues inserted before a stop re-formed from the inserted bases

**In short.** On GRCh37, a record at chromosome 16, position 79245647, with a one-base reference allele and a 46-base alternate, inserts 45 bases after coding base 659 of ENST00000406884, a forward-strand transcript whose stop is residue 235. The bases copy the 45 reference bases that follow, so the 3' rule carries the insertion 45 bases on, to between the stop codon's second and third base. Fifteen residues follow the last residue, and the last two inserted bases re-form a stop codon with the stop codon's last base. Both engines write `inframe_insertion`; Ensembl VEP writes `HGVSp=<fifteen residues inserted between residue 234 and the stop at residue 235>` and vep-rs `HGVSp=<an extension at residue 235: the first new residue, 15 residues added>`.

The nomenclature's Insertion page describes an insertion that extends the full-length protein at its C-terminal end as an extension.

## 5. One base deleted in the stop codon

**In short.** On GRCh37, a record at chromosome 1, position 2238201, with a two-base reference allele and a one-base alternate, deletes the first base of the stop codon of ENST00000378536, a forward-strand transcript whose stop is residue 729. The stop is lost, and the new frame runs through the 3' UTR to a stop codon 47 residues on. Both engines write `stop_lost`; Ensembl VEP writes `HGVSp=<the stop at residue 729, unchanged>` and vep-rs `HGVSp=<an extension at residue 729: the first new residue, 47 residues added>`.

Ensembl VEP's own consequence says the stop is lost, and its HGVSp says the stop is unchanged. The two bases left of the stop codon are read as a partial codon and written as the stop.

## 6. The last residue's codon and the stop codon replaced together

**In short.** On GRCh38, a record at chromosome 3, position 128908990, with a six-base reference allele and a three-base alternate, replaces the last residue's codon and the stop codon of ENST00000681585 by three bases. The transcript lies on the forward strand and its stop is residue 477. Residue 476 changes, the stop is lost, and the alternate protein runs on through the 3' UTR to a stop at position 24 counted from the changed residue. Ensembl VEP writes `<protein>:p.<Aaa>476_Ter477delins<Bbb>`, a deletion-insertion of the range that ends at the stop. vep-rs writes `<protein>:p.<Aaa>476<Bbb>fsTer24`.

The record at chromosome 12, position 112450359, on ENST00001096176 is the same shape with a four-base reference allele and a one-base alternate. The last base of the last residue's codon and the whole stop codon are replaced by one base, and the new stop stands at position 4 from residue 58.

Two rows of the same records' kind keep the codon window's value, because the alternate protein keeps its length and a stop stands at the stop's own residue: chromosome 9, position 127818308, on ENST00001020284, and chromosome X, position 49076647, on ENST00001057692. A change of a whole number of codons that takes one or two bases of the stop codon keeps it as well where it changes or removes a residue before the stop. A deletion of that reach that keeps every residue and loses the stop does not keep it: it is the eighth shape, and section 8's last sentence gives what each engine writes.

## 7. A lost stop on a transcript whose peptide reads a stop codon as a residue

**In short.** On GRCh37, a record at chromosome 1, position 26142208, a substitution of one base, changes the second base of the stop codon of ENST00000374315, a forward-strand transcript whose stop is residue 557. Its coding sequence holds a stop codon at residue 428 that the transcript's peptide reads as a residue. The stop is lost, the whole protein is kept, and 35 residues follow up to a stop of the 3' UTR that is another codon than the one inside. Ensembl VEP writes `<protein>:p.Ter557<Xaa>extTer?`; vep-rs writes `<protein>:p.Ter557<Xaa>extTer35`.

| Record, on GRCh37 | The first stop codon of the 3' UTR | Ensembl VEP's count | vep-rs's count |
|---|---|---|---|
| Chromosome 1, position 26142208, ENST00000374315 | A codon of other bases than the one the peptide reads as a residue | unknown | 35 |
| Chromosome 22, position 19864628, ENST00000400519 | A codon of the same bases as the one the peptide reads as a residue: the model does not say whether it ends the protein | unknown | unknown |

## 8. The stop codon deleted, with residues of the 3' UTR after the protein

**In short.** On GRCh37, a record at chromosome 7, position 151292534, with a four-base reference allele and a one-base alternate, deletes the three bases of the stop codon of ENST00000488258, a reverse-strand transcript whose stop is residue 230. The whole protein is kept and 11 residues follow it up to a stop of the 3' UTR. Ensembl VEP writes `<protein>:p.Ter230delextTer11`, the stop deleted and then the count; vep-rs writes `<protein>:p.Ter230<Xaa>extTer11`.

The Deletion page describes deletions "extending the amino acid sequence at the C-terminal end with one or more amino acids" as an extension, and the extension's syntax names "A single amino acid" after the stop's position. A deletion that reaches from the last codon into the stop codon and leaves the last residue as it is, as at chromosome 20, position 43251710, on ENST00000539235, is written the same way.

## 9. Nine bases inserted before the last residue's codon that hold that residue and a stop

**In short.** On GRCh38, a record written for this example at chromosome 21, position 30282175, with a one-base reference allele and a ten-base alternate, inserts nine bases after coding base 757 of ENST00000340345. The transcript lies on the reverse strand, its protein has 254 residues and its stop is residue 255. Coding base 757 is the first base of the codon of residue 253. The first two inserted bases complete that codon to the same residue. The next six read as a codon of the residue at 254 and a stop codon, which stands at residue 255. The ninth inserted base and the two bases that remain of the codon of residue 253 follow that stop and are not translated.

The rows are those of a run with `--hgvs --hgvsg --shift_hgvs 0`. `--shift_hgvs 0` keeps the description at the record's own place, and `--hgvsg` adds to the `Extra` column an `HGVSg` item that plays no part here. Both engines write `inframe_insertion,stop_retained_variant`. Ensembl VEP writes `<protein>:p.<Aaa>253_<Bbb>254ins<Bbb>Ter`, the residue at 254 and a stop inserted between residues 253 and 254. vep-rs writes `<protein>:p.Ter255=`, the stop at residue 255 unchanged. Without `--shift_hgvs 0` both engines carry the description eight bases on, past the stop codon, and neither writes an `HGVSp` for the record.

The alternate protein is the reference protein, residue for residue: 254 residues and a stop at residue 255. The protein gains neither the residue nor the stop that Ensembl VEP's value names. With the ninth base left out, the eight that remain shift the frame, and the protein is the same. Both engines write `inframe_insertion,stop_retained_variant` for them as well, so Ensembl VEP reads the description's type from the codon window's two peptides, as for the nine. Neither engine's `frameshift` test holds where `stop_retained` does (`Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, line 1549), and `inframe_insertion` holds because the alternate peptide, cut after its first stop, opens with the reference peptide (lines 1126-1127). Each engine writes the `HGVSp` it writes for the nine bases.

**Where Ensembl VEP writes it.** The codon window's reference peptide is residue 253, and its alternate peptide opens with the same residue. Once the shared residue is clipped the reference side is empty, and lines 2063-2064 of `TranscriptVariationAllele.pm` type the window an insertion. Line 1961 cuts the alternate peptide after its first stop, which leaves the residue at 254 and the stop. Line 1991 writes them between residues 253 and 254, the two residues that flank the insertion. Between the type and the format, lines 1808-1821 compare no proteins unless the type is a frameshift.

**What the nomenclature says.** Its Insertion page defines the type in one sentence. "Insertion: a sequence change between the translation initiation (start) and termination (stop) codon where, compared to the reference sequence, one or more amino acids are inserted, which is not a frameshift and where the insertion is not a copy of a sequence immediately N-terminal (5')." Its General page gives the sign `=` (equals) "to indicate a sequence was tested but found unchanged; p.(Arg234=)."

**Why vep-rs's value is the one to keep.** Compared to the reference protein, this protein has no inserted residue, so the Insertion page's definition does not fit it. A reader of Ensembl VEP's value takes the protein for changed. vep-rs's value names the stop at its residue number with the sign for unchanged. The nomenclature's pages on protein changes show that sign on a residue, and none of them shows it on a stop.

The Insertion page has a second sentence, on variants that "encode a translation stop codon in the inserted sequence". It describes them "as an insertion of this sequence, not as a deletion-insertion removing the entire C-terminal amino acid sequence". Read on the inserted bases, that sentence fits this record as it fits the record of section 3, and Ensembl VEP's value is the description it gives. The nomenclature derives a protein description from the comparison of the two proteins, and says that knowledge of the change on the DNA level should not be used. The comparison finds no change here, and vep-rs keeps to it.

**Where an open dataset holds such a row.** gnomAD v4.1 holds one on chromosome 21 of GRCh38, between positions 44809112 and 44809113, on ENST00000997828. The transcript lies on the reverse strand, its protein has 52 residues and its stop is residue 53. The record inserts 42 bases between the second and third base of the codon of residue 51. The protein keeps its 52 residues and its stop at residue 53. Under `--hgvs`, Ensembl VEP writes `<protein>:p.<Aaa>51_<Bbb>52ins<Bbb>Ter` and vep-rs `<protein>:p.Ter53=`. Under `--hgvsp_use_prediction` each engine writes the same description inside parentheses.

| More than one residue before the inserted stop | What stands there |
|---|---|
| The numbers | M numbers the first of the two neighbours in the reference protein, and N that protein's stop. |
| `HGVSp`, Ensembl VEP | `<protein>:p.<Aaa>M_<Bbb>M+1ins<residues>Ter`: the residues from M+1 to the last one, then a stop. |
| `HGVSp`, vep-rs | `<protein>:p.Ter<N>=`, N the number after the last residue's. |
| What the row leaves open | The reference residues after M+1. The class reads them from the transcript's coding sequence in the cache vep-rs ran on, which the format comparator (`tools/vep-formats-compare`) is given as `--json-cache`. |
| What the class holds | Under the transcript's own genetic code the coding sequence spells residue M, the named residues from M+1 to the last, and a stop codon after them at N. The row's terms are not read: the coding sequence shows that the protein is kept, which the term stands for where the row alone is read. |
| What stays a difference | The pair without the cache, under a genetic code the comparator holds no table for, on a coding sequence no stop codon closes, and with a named residue that is not the protein's at its number. |
| A record of the open datasets | None carries the pair: a made-up record of the tests holds it. |
| A record beside it | ClinVar variation 100660, at 2:197705405 on GRCh38 and 2:198570129 on GRCh37, on ENST00000282276: 1,782 inserted bases, a copy of the whole coding sequence, written before coding base 1. |
| The rows | Both assemblies, under `--hgvs`. The transcript lies on the forward strand, and both engines write `5_prime_UTR_variant`, the term of the insertion where the record writes it. |
| Where Ensembl VEP reads it | Where the walk of its 3' shift ends, 1,000 positions on ([`hgvs_insertion_left_where_the_shift_stops`](hgvs-insertion-left-where-the-shift-stops.md)): between coding bases 1000 and 1001, where it writes residues 335 to 593 and a stop as inserted between residues 334 and 335. |
| What vep-rs writes | No `HGVSp`. The copy can stand at every point up to one base past the stop codon, in the 3' UTR, and vep-rs reads both notations at that last place. |
| Which entry the row is | The `HGVSc`, the `HGVSp` and the `HGVS_OFFSET` differ under the entry linked above. Written at a point whose walk ends on a base that differs, such an insertion gives this section's pair. |
| Why Ensembl VEP's residues are the protein's own | The 1,782 bases are 594 codons, and coding base 1000 is the first base of the codon of residue 334. Moved 1,000 positions on, the copy's bases are turned by as many, so it opens with coding base 1001. Read from there, it completes that codon to the same residue and spells residues 335 to 593 and the stop codon at residue 594. What follows that stop is not translated: the alternate protein is the reference protein. |

## 10. Six bases inserted inside the last residue's codon that keep the protein

**In short.** Ensembl VEP writes no `HGVSp` for this record, and vep-rs the stop unchanged.

| The record | On GRCh38, under `--hgvs --shift_hgvs 0` |
|---|---|
| Input | Chromosome 21, position 17565688, a one-base reference allele and a seven-base alternate, written for this example. `crates/vep-cli/tests/hgvsp_insertion_in_the_last_codon.rs` holds it as `fwd_6_pure`. |
| Transcript | ENST00000284878, on the forward strand. Its protein has 365 residues and its stop is residue 366. |
| Place | After coding base 1094, between the second and third base of the codon of residue 365. |
| Alternate protein | The reference protein. The first inserted base completes the codon to the same residue, and the next three are a stop codon at residue 366. |
| Terms, both engines | `inframe_insertion,stop_retained_variant` |
| `HGVSp`, Ensembl VEP | none |
| `HGVSp`, vep-rs | `<protein>:p.Ter366=` |
| Without `--shift_hgvs 0` | The same rows: the base after the insertion is not its first base, so nothing is carried on. |
| Why vep-rs's value is kept | The two proteins are equal, residue for residue. Section 9 gives the nomenclature's sign for that. |

| Where Ensembl VEP writes none | Lines | On this record |
|---|---|---|
| `stop_retained` holds for an alternate peptide of whole codons, so `frameshift` returns 0 | `Utils/VariationEffect.pm`, 1313-1315 and 1549 | The alternate peptide is residue 365, a stop and one more residue. |
| The clipped codon window is typed an insertion | `TranscriptVariationAllele.pm`, 2063-2064 | The reference side is empty, the start 366 and the end 365. |
| No residue follows the last one to flank the insertion | 2140-2144, 2368-2380 | The peptide has 365 residues, and the flank is asked for from residue 365. |
| The description is given up | 2146, 1815-1818 | No `HGVSp`. |

| Other records of the kind | `HGVSp`, Ensembl VEP | `HGVSp`, vep-rs |
|---|---|---|
| Nine bases at the same place | none | the stop unchanged |
| Five, six, eight or nine bases on a reverse-strand transcript, ENST00000284971, whose stop is residue 109 | none | the stop unchanged |
| Five or eight bases at the same place, where Ensembl VEP's terms are `stop_gained,frameshift_variant` ([`stop_codon_overlap_insertion_range`](stop-codon-overlap-insertion-range.md)) and it walks the two translations to the stop | the stop unchanged | the stop unchanged |
| Six or nine bases held with a base the two alleles share, in a record of two alternate alleles | none | the stop unchanged |

| What the class reads for the tenth shape | On the row |
|---|---|
| Ensembl VEP's value | none |
| vep-rs's value | the stop unchanged, alone or inside prediction parentheses |
| The terms | they hold `stop_retained_variant` |
| `Protein_position` | one codon, numbered one below the stop in vep-rs's value |
| `CDS_position` | inside that codon |
| `Amino_acids` | that codon's residue against the same residue and then a stop |
| What it cannot tell | The row does not number the reference protein's stop. The class takes the number from vep-rs's value and holds it to the residue after the codon. |
| What stays a difference | The same pair with the sides exchanged, and every other empty value on Ensembl VEP's side. |

## 11. A stop lost on a transcript whose peptide holds a stop

**In short.** Ensembl VEP writes these extensions with the unknown count, or with a count from the peptide's last residue; vep-rs counts from the lost stop to the next stop of the frame.

| Record, a one-base substitution | Transcript, the stop lost and the next stop of the frame | `HGVSp`, Ensembl VEP | `HGVSp`, vep-rs |
|---|---|---|---|
| GRCh37, chromosome 4, position 100268190 | ENST00000515683, a polymorphic pseudogene whose peptide holds one stop, residue 78 of 375. The next stop is the stop codon that closes the coding sequence, residue 376. | `<protein>:p.Ter78GlyextTer?` | `<protein>:p.Ter78GlyextTer298` |
| GRCh37, chromosome 12, position 6940463 | ENST00000606935, whose peptide of 736 residues holds 18 stops, the first at residue 167. The stop of residue 356 is lost; the next is the stop the peptide holds at residue 435. | `<protein>:p.Ter356TrpextTer?` | `<protein>:p.Ter356TrpextTer79` |
| GRCh38, chromosome 11, position 104886747 | ENST00000703420, a `protein_coding_LoF` transcript whose peptide of 263 residues holds four stops, the first at residue 125, and whose coding sequence no stop codon closes. The last of the four is lost, residue 253; the next stop is residue 280, in the 3' UTR. | `<protein>:p.Ter253SerextTer?` | `<protein>:p.Ter253SerextTer27` |
| GRCh37, chromosome 18, position 61390631 | ENST00000382749, whose peptide of 392 residues holds a stop at residue 90. The stop codon that closes the coding sequence is lost, residue 393; the next stop is residue 436, in the 3' UTR. | `<protein>:p.Ter393GluextTer?` | `<protein>:p.Ter393GluextTer43` |
| GRCh38, chromosome 19, position 39247459 | ENST00000634680, whose peptide holds a stop at residue 52 of 108 and whose coding sequence no stop codon closes. The next stop is residue 144, in the 3' UTR. | `<protein>:p.Ter52CysextTer35` | `<protein>:p.Ter52CysextTer92` |
| The same record | ENST00000606380, whose peptide of 123 residues holds no stop. The same codon closes its coding sequence, residue 124; the next stop is residue 216, the bases of residue 144 of the transcript above. | `<protein>:p.Ter124CysextTer92` | `<protein>:p.Ter124CysextTer92` |

| Record, a deletion of whole residues | Transcript, the stop lost and the first stop after the change | `HGVSp`, Ensembl VEP | `HGVSp`, vep-rs |
|---|---|---|---|
| GRCh38, chromosome 2, position 227608903, 36 bases deleted | ENST00000641246, whose peptide of 105 residues holds one stop, residue 48. Residues 41 to 52 go, that stop among them. The stop codon that closes the coding sequence, residue 106, then stands at residue 94. | `<protein>:p.Gly41_Gly52delextTer?` | `<protein>:p.Gly41_Gly52delextTer46` |
| GRCh38, chromosome 12, position 40485606, 90 bases deleted | ENST00000454784, whose peptide of 8,222 residues holds 85 stops, the first at residue 4198. Residues 4223 to 4252 go, four stops among them, the first at residue 4224. The stop of residue 4278 then stands at residue 4248. | `<protein>:p.Arg4223_Asp4252delextTer?` | `<protein>:p.Arg4223_Asp4252delextTer24` |
| GRCh38, chromosome 12, position 40493980, 18 bases deleted | The same transcript. The deletion opens at the stop of residue 6131. The stop of residue 6171 then stands at residue 6165. | `<protein>:p.Ter6131_Ser6136delextTer?` | `<protein>:p.Ter6131_Ser6136delextTer34` |
| GRCh37, chromosome 4, position 100268187, 6 bases deleted, written for this example | ENST00000515683 of the first table. Residue 77 and the stop of residue 78 go. The stop codon that closes the coding sequence then stands at residue 374. | `<protein>:p.Val77_Ter78delextTer?` | `<protein>:p.Val77_Ter78delextTer296` |
| GRCh37, chromosome 12, position 6940460, 6 bases deleted, written for this example | ENST00000606935 of the first table. The stop of residue 356 and residue 357 go. The stop of residue 435 then stands at residue 433. | `<protein>:p.Ter356_Thr357delextTer?` | `<protein>:p.Ter356_Thr357delextTer77` |
| GRCh37, chromosome 18, position 61390627, 6 bases deleted, written for this example | ENST00000382749 of the first table. Residue 392 and the stop codon that closes the coding sequence, residue 393, go. The stop of the 3' UTR then stands at residue 434. | `<protein>:p.Pro392_Ter393delextTer?` | `<protein>:p.Pro392_Ter393delextTer41` |
| GRCm39, chromosome 19, position 13259544, 12 bases deleted, a record of the committed mouse corpora | ENSMUST00000208760, whose peptide of 307 residues holds six stops, the first at residue 210. Residues 289 to 292 go, the stop of residue 290 among them. The stop of residue 305 then stands at residue 301. | `<protein>:p.Ala289_Thr292delextTer?` | `<protein>:p.Ala289_Thr292delextTer11` |
| GRCh38, chromosome 12, position 40486056, 30 bases deleted | ENST00000454784. Residues 4375 to 4384 go, the stop of residue 4378 among them. The stop of residue 4388 then stands at residue 4378, the lost stop's own. | `<protein>:p.Asn4375_Trp4384delextTer?` | `<protein>:p.Asn4375_Trp4384delextTer?` |
| GRCh38, chromosome 12, position 40488980, 12 bases deleted | ENST00000454784. Residues 5344 to 5347 go, the last of them the stop that is lost. The stop of residue 5373 then stands at residue 5369. | `<protein>:p.Thr5344_Ter5347delextTer?` | `<protein>:p.Thr5344_Ter5347delextTer22` |
| GRCh38, chromosome 2, position 241873267, 6 bases deleted | ENST00000343216, whose peptide of 572 residues holds no stop. Residue 572 and the stop codon that closes the coding sequence, residue 573, go. The stop of the 3' UTR then stands at residue 588. The codon window types the change a deletion, so both engines write the deletion and the count; vep-rs writes the frameshift of this page's opening part only for a range that the codon window types a deletion-insertion. | `<protein>:p.Val572_Ter573delextTer15` | `<protein>:p.Val572_Ter573delextTer15` |

| Where Ensembl VEP's count comes from | `TranscriptVariationAllele.pm` | On these records |
|---|---|---|
| The first stop of the whole alternate translation is taken, wherever it stands | Line 2503 | On the second, third and fourth it is a stop the peptide holds nearer its start than the lost one: residues 167, 125 and 90. |
| One and the peptide's length are taken from that stop's residue number | Line 2513 | 376 - 1 - 375 = 0 on the first. 167 - 1 - 736, 125 - 1 - 263 and 90 - 1 - 392 are below zero on the next three. 144 - 1 - 108 = 35 on the fifth: the residues after the peptide's last one. |
| A figure that is not above zero gives no count | Lines 2519-2525 and 1924 | The unknown count. |
| No count is made where the lost stop is residue 1 | Line 2483, called at line 1922 | A peptide that opens with a stop, as a short immunoglobulin segment's can. |
| A deletion's count is asked for the same way | Lines 1916 and 1922 | The same arithmetic: 94 - 1 - 105 on the first deletion, 4198 - 1 - 8222 on the second, and 588 - 1 - 572 = 15 on the last, whose peptide holds no stop. |

| What vep-rs writes | |
|---|---|
| The count | The next stop's residue number less the lost stop's. The next stop is sought from the lost stop through the coding sequence and the 3' UTR, and no farther. |
| No stop there | The unknown count, as Ensembl VEP writes it. |
| The stop that is lost | A stop the peptide holds, the first residue's included, or the stop codon that closes the coding sequence of a transcript whose peptide holds a stop before it. |
| A peptide that is not the translation of its coding sequence at every residue, a first residue written as a methionine over another first codon set aside | The value both engines write. A peptide that reads a stop codon as a residue is of this kind; its lost closing stop is counted as section 7 says. |
| A deletion of whole residues that takes a stop the peptide holds, or the closing stop behind one | The form both engines write, the residues deleted and then `delextTer`, with a count. The count is the residue number that the first stop at or after the first residue that differs has in the protein after the deletion, less the number the lost stop had before it. On a peptide that holds no stop the same subtraction is the figure Ensembl VEP writes, so nothing changes there. |
| The stop a deletion's count runs from | The first stop among the residues deleted. The description names it only where its range opens or ends at that stop. |
| A stop the peptide holds, or the closing stop behind one, deleted alone | Both engines write the stop as deleted and then an extension's count: Ensembl VEP the unknown count, vep-rs the count derived as above. Section 8's form, the residue in the stop's place, is written only where vep-rs compares the two proteins: on a peptide that holds no stop, whose coding sequence one stop codon closes. |
| A deletion of either of the two rows above whose count so derived is not above nought | The unknown count both engines write. Ensembl VEP's own figure is not above nought there either, and lines 2519-2525 give no count for such a figure on any transcript. |
| A frameshift whose first changed residue is a stop the peptide holds | An extension named at that stop in both engines. Its count is another entry's: [`hgvsp_frameshift_count_past_a_stop_read_as_a_residue`](hgvsp-frameshift-count-past-a-stop-read-as-a-residue.md). |
| Why vep-rs's count is the right one | The nomenclature's Extension page writes the unknown count where the frame "does not contain a new stop codon". On these records the transcript's own bases hold one. For the lost codon of the fifth record both engines count 92 on the transcript whose coding sequence it closes. |

| What the class reads for the eleventh shape | On the row |
|---|---|
| Both values | An extension at one residue number with one first residue, or one deletion with an extension's count after it, alone or inside prediction parentheses. |
| The counts | vep-rs's is a number. Ensembl VEP's is the one its lines give on the cache's sequences: unknown, or the residues after the peptide's last one. |
| The terms | They hold `stop_lost`. |
| vep-rs's cache | The row's transcript has a peptide that holds a stop and is the translation of its coding sequence. For an extension named at a stop, the residue named is a stop it holds, or the closing stop, and the next stop of the frame stands as many residues after it as vep-rs's count. |
| `CDS_position`, for a deletion | The coding bases the row names are three for each residue the deletion names. With them out of the cache's coding sequence, and the 3' UTR read after it, the first residue that differs is the deletion's first. The stop lost is the first stop of the reference at or after that residue: a stop the peptide holds, or the closing stop. The first stop of the alternate at or after that residue stands vep-rs's count of places after the lost stop's. |
| `FLAGS` | The column is not consulted: a row whose `FLAGS` name `cds_end_NF` is taken like any other. |
| What stays a difference | For this shape: another count on either side, the sides exchanged, and a peptide that holds no stop. Without the cache the seventh shape, read from the row alone, still takes an unknown count against a count above residue 1 on a row whose `FLAGS` name no `cds_end_NF`; every other pair stays a difference. |

| Where real data holds such rows, under `--hgvs` | Rows on which vep-rs's count is not Ensembl VEP's value | Rows on which both write the unknown count: no stop follows through the 3' UTR |
|---|---|---|
| ClinVar, GRCh37 | 4 | 0 |
| ClinVar, GRCh38 | 14 | 0 |
| 1000 Genomes Phase 3, GRCh37 | 47 | 19 |
| gnomAD v2.1.1, GRCh37 | 74 | 10 |
| 1000 Genomes high coverage, GRCh38 | 50 | 2 |
| gnomAD v4.1, GRCh38 | 324: Ensembl VEP writes the unknown count on 320 and a count from the peptide's last residue on 4 | 48 |

| Where real data holds such deletions, under `--hgvs`: rows the table above does not count | Rows on which vep-rs writes a count and Ensembl VEP the unknown count | Rows on which both write the unknown count: the count is not above nought, or no stop follows |
|---|---|---|
| ClinVar, GRCh37 and GRCh38 | 0 | 0 |
| The committed mouse corpora (`tests/golden/116/GRCm39-mouse` and `GRCm39-mouse-merged`), one record of GRCm39 | 1 | 0 |
| 1000 Genomes Phase 3, GRCh37 | 0 | 1 |
| gnomAD v2.1.1, GRCh37 | 0 | 2 |
| 1000 Genomes high coverage, GRCh38 | 5 | 3 |
| gnomAD v4.1, GRCh38 | 97 | 37 |
