# The HGVSp that names a stop on an incomplete coding sequence

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.3 writes the unknown residue of a partial codon, and of a first codon the coding sequence opens inside, as the stop. It so names a stop at a residue where the reference protein holds none, unless the record shows polyadenylation completing the partial codon to a stop codon. It names one too where a frameshift reaches past the last whole codon, since it closes a frameshift's reference protein with a stop whatever the coding sequence's end. vep-rs writes no HGVSp where the description would name such a stop alone.

**Why this is a defect.** The HGVS nomenclature's Reference Sequences page takes as a reference only a transcript with "a complete CDS", whose "translation termination codon must be clearly annotated", and names the Ensembl transcripts flagged `cds_start_NF` or `cds_end_NF` as incomplete. A description prefixes a position with the reference residue found there; where the coding sequence ends inside a codon its record does not show completed to a stop codon, or opens inside one, the reference holds no stop to name.

**What vep-rs does instead.** vep-rs writes no HGVSp where the description would name the stop alone, as unchanged, deleted, replaced by one residue or extended, on a coding sequence that no stop codon closes; nor for the stop written as unchanged at residue 1. A coding sequence that ends inside a codon is closed by no stop codon unless its record shows polyadenylation completing that codon to one, whatever its last three bases spell: they are no codon of the frame its protein is read in. A stop that the transcript's own peptide carries at the residue named keeps its notation, and so does every other description. Only the HGVSp changes.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 692-787, `peptide`: lines 774-776 append the unknown-residue symbol to a translation whose codon string ends in an incomplete codon.
- Lines 2114-2185, `_get_hgvs_peptides`: lines 2181-2182 rewrite every unknown-residue name of the reference and the alternate peptide to the stop name, so the residue of an incomplete codon is written as a stop.
- Lines 2320-2365, `_get_fs_peptides`: line 2338 closes the reference translation of a frameshift with a stop whatever the coding sequence's end, so the residue after the last whole codon is read as a stop.

**What correcting it would change.** The rows would carry no HGVSp. No other column would change.

**The class it produces.** `hgvsp_stop_named_on_incomplete_coding_sequence`: Ensembl VEP's description names the stop alone in one of the four ways and vep-rs writes none, on a row that shows the incomplete coding sequence. The row shows it by `FLAGS` naming `cds_end_NF`, by the term `incomplete_terminal_codon_variant`, or by the stop unchanged at residue 1 under `FLAGS` naming `cds_start_NF`. A row on a transcript that no flag marks incomplete, whose terms hold no such term, is shown by the candidate's cache: the format comparator, given `--json-cache`, and the corpus builder read there that the transcript's coding sequence is no whole number of codons. A comparison that names no cache leaves such a row a difference.

Given the candidate's cache, the class takes no row on a transcript whose record shows polyadenylation completing its stop codon: the reference protein holds that stop.

## 1. A substitution in a partial last codon

**In short.** On GRCh37, a record at chromosome 1, position 3548777, a substitution of one base, changes the last base of the coding sequence of ENST00000424367, a reverse-strand transcript flagged `cds_end_NF` whose coding sequence ends on the first base of a codon. Ensembl VEP writes `<protein>:p.Ter305=`, a stop unchanged at a residue the protein does not have; vep-rs writes no HGVSp.

## 2. A substitution in an incomplete first codon

**In short.** On GRCh37, a record at chromosome 1, position 865692, a substitution of one base, changes the first coding base of ENST00000341065, a forward-strand transcript flagged `cds_start_NF` whose coding sequence opens on the second base of its first codon and is closed by a stop codon. Ensembl VEP writes `<protein>:p.Ter1=`; the reference protein holds an unknown residue there, and vep-rs writes no HGVSp.

**Where the unknown residue arises.** The transcript's start exon has phase 1. `Bio/EnsEMBL/Transcript.pm` of Ensembl core at release/116, lines 917-920, `translateable_seq`, puts one unknown base before such a coding sequence, and Ensembl VEP's cache stores the result (`Bio/EnsEMBL/VEP/AnnotationSource/Database/Transcript.pm` of `ensembl-vep` at release/116.3, line 351). `TranscriptVariationAllele.pm`, lines 868 and 877, `codon`, cuts the first codon from it: the unknown base and the first two coding bases. Lines 730-736, `peptide`, translate that codon, and a codon with an unknown first base names the unknown residue, on the reference and on the changed allele. Lines 1791-1793, `hgvs_protein`, take the two peptides as the notation's, and lines 2181-2182, listed above, write both as the stop.

## 3. A partial last codon whose base completes a stop codon in another frame

**In short.** On GRCh37, ClinVar holds the three substitutions of chromosome 22, position 23235998, the last base of the coding sequence of ENST00000390320, a forward-strand transcript flagged `cds_start_NF`. Its protein is read from the first of its 127 coding bases: 42 codons and the first base of a 43rd. The last three bases spell a stop codon in another frame; one substitution makes another stop codon of them and two a residue's codon. Ensembl VEP writes `<protein>:p.Ter43=` on all three rows, whose terms hold `incomplete_terminal_codon_variant`; vep-rs writes no HGVSp on any.

## 4. A frameshift carried past the last whole codon

**In short.** On GRCh38, ClinVar holds at chromosome 3, position 37050649, an insertion of one base in ENST00000413740, a forward-strand transcript that no flag marks incomplete: its coding sequence is 611 codons and two more bases. The 3' rule carries the insertion two bases on, to between the last whole codon and those two bases, with which the inserted base spells a stop codon. Ensembl VEP closes the reference translation with a stop and writes `<protein>:p.Ter612=`, on a row whose terms are `frameshift_variant` alone; the reference protein holds 611 residues and no stop, and vep-rs writes no HGVSp.
