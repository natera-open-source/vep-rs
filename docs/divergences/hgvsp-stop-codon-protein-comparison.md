# The HGVSp of a change at the stop codon

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.2 describes a change at the stop codon, and an insertion before it that brings a stop to the stop's residue, from the peptides of the allele's codon window. A change that keeps the whole protein is so written as a deletion-insertion of the stop, as an insertion before the stop or the last residue, or as the stop unchanged where the stop is lost. A change that reaches from before the stop codon through it is written as a deletion-insertion of the range that ends at the stop. vep-rs writes the first as the stop unchanged or as an extension named at the stop, and the second as a frameshift.

**Why this is a defect.** The HGVS nomenclature derives a protein description "from comparing the variant protein sequence with the reference protein sequence" and adds that knowledge of the change on the DNA level should not be used. A change that keeps every residue and adds residues after the last one extends "the reference amino acid sequence at the N- or C-terminal end", which is the nomenclature's definition of an extension. Its Insertion, Deletion and Deletion-insertion pages each send a change that extends the full-length protein at its C-terminal end to the Extension page. A deletion-insertion that replaces a residue by the same residue is, in the nomenclature's words, one "which effectively is no change".

The Frameshift page describes deletions "starting 5' of and including the translation termination (stop) codon" as a frameshift, and the Deletion-Insertion page says the same of deletion-insertions. A frameshift is, in the Frameshift page's words, a change where translation "shifts to another reading frame".

**What vep-rs does instead.** vep-rs compares the alternate protein with the reference protein at the stop codon, and for an insertion before the last residue's codon that brings a stop codon of its own to the stop's residue. Where the two are equal and a stop stands at the same residue, it writes the stop at its residue number as unchanged. Where the whole reference protein is kept and residues follow it, it writes an extension named at the stop. The extension carries the stop's residue number, the first new residue and the count of residues added up to the new stop, or `?` where the transcript holds none. Only the HGVSp changes.

Where a deletion or deletion-insertion reaches from before the stop codon into it, and runs to its last base or is not a whole number of codons, vep-rs writes a frameshift. The change is read less the bases its two alleles share at either end. It names the first changed residue, the residue in its place and the position of the new stop counted from it, or `?` where the transcript holds none.

On a transcript whose peptide reads a stop codon inside the coding sequence as a residue, a substitution that loses the closing stop is an extension whose count Ensembl VEP leaves unknown, because its count ends at the first stop codon of the translated bases. vep-rs takes the peptide as the reference protein and writes the count of residues added, unless the first stop codon of the 3' UTR is a codon the peptide reads as a residue.

Where an in-frame deletion takes the stop codon and keeps every residue, Ensembl VEP writes the stop as deleted, then the extension's count, with no residue in the stop's place. vep-rs writes the extension named at the stop with the residue that stands in the stop's place and the same count.

On a coding sequence that no stop codon closes, vep-rs writes no HGVSp where Ensembl VEP's description names the stop alone, as unchanged, deleted, replaced by one residue or extended, at a residue where the transcript's peptide carries no stop. At an incomplete first codon it writes none where Ensembl VEP writes the stop as unchanged at residue 1. Both are a difference of another class, [`hgvsp_stop_named_on_incomplete_coding_sequence`](hgvsp-stop-named-on-incomplete-coding-sequence.md). Every other description on such rows keeps the value both engines write.

A lost stop with no whole codon after it keeps the value both engines write. So do an allele at an exon's edge, where the spliced sequence is not known, and a stop written as unchanged on a row whose own terms do not hold `stop_lost`. A change into the stop codon after which the protein keeps its length keeps that value too.

A change of a whole number of codons that ends before the stop codon's last base keeps the value both engines write where it changes or removes a residue before the stop. A deletion of that reach that keeps every residue and loses the stop is the in-frame deletion above, written as the extension. The frameshift is written only where a residue before the stop changes.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 1898-2044, `_get_hgvs_protein_format`: one chain of tests over a type and the reference and alternate peptides of the codon window. Line 1911 writes the stop as unchanged when the two peptides are equal. Line 1916 writes an extension only when `stop_lost` holds and the type is a deletion or a substitution. Lines 1958-1993 write a deletion-insertion or an insertion, and line 1961 cuts every residue after the alternate's first stop.
- Lines 2047-2111, `_get_hgvs_protein_type`: lines 2072-2077 type two peptides of unequal length, neither of them empty, a deletion-insertion. The stop alone against a longer alternate peptide is so written as a deletion-insertion of the stop. A range that ends at the stop takes its type at the same lines.
- Lines 2063-2064 type a codon window an insertion where its reference peptide is empty once the residues both peptides share are clipped. Where inserted bases read as the last residue and a stop, line 1961 keeps both, and line 1991 writes them as inserted between the two residues that flank them. Lines 1808-1821 call the type, the peptides and the format; they compare no proteins unless the type is a frameshift.
- Lines 2188-2295, `_clip_alleles`: lines 2231-2242 chop the stop that ends both peptides. What is left of an alternate that ends in a re-formed stop is typed an insertion at lines 2271-2285, and lines 1985-1992 write it between the last residue and the stop.
- Lines 692-787, `peptide`: lines 774-776 give the one or two bases left of the stop codon the unknown-residue symbol, and line 2060 rewrites the reference's stop to the same symbol. The two peptides are then equal, and line 1911 writes the stop as unchanged before line 1916 is asked.
- Lines 1916-1919 take the extension's new residue from the alternate peptide, which lines 2167-2169 of `_get_hgvs_peptides` set to the word for a deletion: the stop is written as deleted before the count. The count is that of `_stop_loss_extra_AA`, lines 2477-2526: the place of the first stop in the translated alternate sequence, less one, less the length of the reference peptide. It is unknown when no stop is met or the figure is not above zero.

**What correcting it would change.** The rows would carry the stop written as unchanged, the extension with its first new residue and its count, or the frameshift. No other column would change: the consequence terms, the coding notation, the protein position, the amino acids and the codons are the same in both engines' rows.

**The class it produces.** `hgvsp_stop_codon_protein_comparison`: two `HGVSp` values that name one protein, in nine shapes. In the first five both name the stop at one residue number above 1. In the first, Ensembl VEP's value is the stop replaced by a stop and vep-rs's the stop unchanged.

In the other four vep-rs's value is an extension named at the stop. Ensembl VEP's is then the stop alone replaced by the first residues that follow the protein, with or without a stop after them that the protein does not hold. Or it is the stop replaced by every added residue and the stop the inserted bases bring. Or it is residues inserted between the last residue and the stop, or the stop unchanged on a row whose terms hold `stop_lost`.

In the sixth shape Ensembl VEP's value is a range from a residue to the stop replaced by residues, none of them a stop, on a row whose terms hold `stop_lost`. vep-rs's is a frameshift named at that residue with the first of them. Its count is unknown, or above their number and not the one that leaves the new stop at the stop's own residue.

In the seventh shape both values are an extension at the same stop with the same first residue, on a row whose terms hold `stop_lost`: Ensembl VEP's with an unknown count, vep-rs's with the count of residues added.

In the eighth shape Ensembl VEP's value is the stop written as deleted before an extension's count, on a row whose terms hold `stop_lost`, and vep-rs's the extension at the same residue with the same count and the first new residue.

In the ninth shape Ensembl VEP's value is one residue and a stop inserted between two neighbouring residues, on a row whose terms hold `stop_retained_variant`. The inserted residue is the second of the two neighbours. vep-rs's value is the stop unchanged at the residue number after that neighbour's. Section 9 shows the pair.

The class reads the two values and the row alone; a row whose `FLAGS` name `cds_end_NF` is no row of it. The terms it reads are those of Ensembl VEP's row. It holds the extension's first residue to the first residue Ensembl VEP names, and its count to the number of residues Ensembl VEP names wherever that value names them all. It cannot tell from the row that the whole reference protein is kept, nor check a count where Ensembl VEP names only the first residues. The same values with the sides exchanged, an empty value on Ensembl VEP's side, and the stop written as deleted with no extension's count after it stay differences.

For the sixth shape the class cannot tell from the row that the change runs to the stop codon's last base or shifts the reading frame. A range that ends at the stop against a frameshift at another residue or with another first residue stays a difference, and so does one on a row whose terms hold no `stop_lost`.

For the ninth shape the class reads the alternate protein from Ensembl VEP's value: the reference protein up to the second neighbour's residue, and then a stop. That value does not number the reference protein's stop. The class takes the number from vep-rs's value and holds it to the residue after the second neighbour.

With more than one residue before the inserted stop, the row does not name the reference residues past the second neighbour, and the class cannot tell that the protein is kept. Such a pair stays a difference. So do an inserted residue other than the second neighbour's, an insertion between two residues that are not neighbours, the same values with the sides exchanged, and a row whose terms hold no `stop_retained_variant`.

The class reads the ninth shape's pair inside prediction parentheses, as `--hgvsp_use_prediction` writes both values, the way it reads it outside them. One value alone inside parentheses stays a difference. So does Ensembl VEP's value against a row on which vep-rs writes no `HGVSp`, and Ensembl VEP's stop unchanged against such a row.

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

Two rows of the same records' kind keep the value both engines write, because the alternate protein keeps its length and a stop stands at the stop's own residue: chromosome 9, position 127818308, on ENST00001020284, and chromosome X, position 49076647, on ENST00001057692. A change of a whole number of codons that takes one or two bases of the stop codon keeps it as well where it changes or removes a residue before the stop. A deletion of that reach that keeps every residue and loses the stop does not keep it: it is the eighth shape, and section 8's last sentence gives what each engine writes.

## 7. A lost stop on a transcript whose peptide reads a stop codon as a residue

**In short.** On GRCh37, a record at chromosome 1, position 26142208, a substitution of one base, changes the second base of the stop codon of ENST00000374315, a forward-strand transcript whose stop is residue 557. Its coding sequence holds a stop codon at residue 428 that the transcript's peptide reads as a residue. The stop is lost, the whole protein is kept, and 35 residues follow up to a stop of the 3' UTR that is another codon than the one inside. Ensembl VEP writes `<protein>:p.Ter557<Xaa>extTer?`; vep-rs writes `<protein>:p.Ter557<Xaa>extTer35`.

Where the first stop codon of the 3' UTR is the same codon as the one the peptide reads as a residue, as for the record at chromosome 22, position 19864628, on ENST00000400519, the model does not say whether it ends the protein. Both engines write the unknown count there.

## 8. The stop codon deleted, with residues of the 3' UTR after the protein

**In short.** On GRCh37, a record at chromosome 7, position 151292534, with a four-base reference allele and a one-base alternate, deletes the three bases of the stop codon of ENST00000488258, a reverse-strand transcript whose stop is residue 230. The whole protein is kept and 11 residues follow it up to a stop of the 3' UTR. Ensembl VEP writes `<protein>:p.Ter230delextTer11`, the stop deleted and then the count; vep-rs writes `<protein>:p.Ter230<Xaa>extTer11`.

The Deletion page describes deletions "extending the amino acid sequence at the C-terminal end with one or more amino acids" as an extension, and the extension's syntax names "A single amino acid" after the stop's position. A deletion that reaches from the last codon into the stop codon and leaves the last residue as it is, as at chromosome 20, position 43251710, on ENST00000539235, is written the same way.

## 9. Nine bases inserted before the last residue's codon that hold that residue and a stop

**In short.** On GRCh38, a record written for this example at chromosome 21, position 30282175, with a one-base reference allele and a ten-base alternate, inserts nine bases after coding base 757 of ENST00000340345. The transcript lies on the reverse strand, its protein has 254 residues and its stop is residue 255. Coding base 757 is the first base of the codon of residue 253. The first two inserted bases complete that codon to the same residue. The next six read as a codon of the residue at 254 and a stop codon, which stands at residue 255. The ninth inserted base and the two bases that remain of the codon of residue 253 follow that stop and are not translated.

The rows are those of a run with `--hgvs --hgvsg --shift_hgvs 0`. `--shift_hgvs 0` keeps the description at the record's own place, and `--hgvsg` adds to the `Extra` column an `HGVSg` item that plays no part here. Both engines write `inframe_insertion,stop_retained_variant`. Ensembl VEP writes `<protein>:p.<Aaa>253_<Bbb>254ins<Bbb>Ter`, the residue at 254 and a stop inserted between residues 253 and 254. vep-rs writes `<protein>:p.Ter255=`, the stop at residue 255 unchanged. Without `--shift_hgvs 0` both engines carry the description eight bases on, past the stop codon, and neither writes an `HGVSp` for the record.

The alternate protein is the reference protein, residue for residue: 254 residues and a stop at residue 255. The protein gains neither the residue nor the stop that Ensembl VEP's value names. With the ninth base left out, the eight that remain shift the frame, and the protein is the same. Both engines write `inframe_insertion,stop_retained_variant` for them as well, so Ensembl VEP reads the type from the peptides as it does for the nine. Each engine writes the `HGVSp` it writes for the nine bases.

**Where Ensembl VEP writes it.** The codon window's reference peptide is residue 253, and its alternate peptide opens with the same residue. Once the shared residue is clipped the reference side is empty, and lines 2063-2064 of `TranscriptVariationAllele.pm` type the window an insertion. Line 1961 cuts the alternate peptide after its first stop, which leaves the residue at 254 and the stop. Line 1991 writes them between residues 253 and 254, the two residues that flank the insertion. Between the type and the format, lines 1808-1821 compare no proteins unless the type is a frameshift.

**What the nomenclature says.** Its Insertion page defines the type in one sentence. "Insertion: a sequence change between the translation initiation (start) and termination (stop) codon where, compared to the reference sequence, one or more amino acids are inserted, which is not a frameshift and where the insertion is not a copy of a sequence immediately N-terminal (5')." Its General page gives the sign `=` (equals) "to indicate a sequence was tested but found unchanged; p.(Arg234=)."

**Why vep-rs's value is the one to keep.** Compared to the reference protein, this protein has no inserted residue, so the Insertion page's definition does not fit it. A reader of Ensembl VEP's value takes the protein for changed. vep-rs's value names the stop at its residue number with the sign for unchanged. The nomenclature's pages on protein changes show that sign on a residue, and none of them shows it on a stop.

The Insertion page has a second sentence, on variants that "encode a translation stop codon in the inserted sequence". It describes them "as an insertion of this sequence, not as a deletion-insertion removing the entire C-terminal amino acid sequence". Read on the inserted bases, that sentence fits this record as it fits the record of section 3, and Ensembl VEP's value is the description it gives. The nomenclature derives a protein description from the comparison of the two proteins, and says that knowledge of the change on the DNA level should not be used. The comparison finds no change here, and vep-rs keeps to it.

**Where an open dataset holds such a row.** gnomAD v4.1 holds one on chromosome 21 of GRCh38, between positions 44809112 and 44809113, on ENST00000997828. The transcript lies on the reverse strand, its protein has 52 residues and its stop is residue 53. The record inserts 42 bases between the second and third base of the codon of residue 51. The protein keeps its 52 residues and its stop at residue 53. Under `--hgvs`, Ensembl VEP writes `<protein>:p.<Aaa>51_<Bbb>52ins<Bbb>Ter` and vep-rs `<protein>:p.Ter53=`. Under `--hgvsp_use_prediction` each engine writes the same description inside parentheses.
