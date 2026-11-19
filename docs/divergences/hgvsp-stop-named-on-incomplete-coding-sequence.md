# The HGVSp that names a stop on an incomplete coding sequence

**The defect in one sentence.** Under `--hgvs`, Ensembl VEP 116.2 writes the unknown residue of a partial codon, and of a first codon the coding sequence opens inside, as the stop, and so names a stop at a residue where the reference protein holds none. vep-rs writes no HGVSp where the description would name such a stop alone.

**Why this is a defect.** The HGVS nomenclature's Reference Sequences page takes as a reference only a transcript with "a complete CDS", whose "translation termination codon must be clearly annotated", and names the Ensembl transcripts flagged `cds_start_NF` or `cds_end_NF` as incomplete. A description prefixes a position with the reference residue found there; where the coding sequence ends inside a codon or opens inside one, the reference holds no stop to name.

**What vep-rs does instead.** vep-rs writes no HGVSp where the description would name the stop alone, as unchanged, deleted, replaced by one residue or extended, on a coding sequence that no stop codon closes; nor for the stop written as unchanged at residue 1. A stop that the transcript's own peptide carries at the residue named keeps its notation, and so does every other description. Only the HGVSp changes.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 692-787, `peptide`: lines 774-776 append the unknown-residue symbol to a translation whose codon string ends in an incomplete codon.
- Lines 2114-2185, `_get_hgvs_peptides`: lines 2181-2182 rewrite every unknown-residue name of the reference and the alternate peptide to the stop name, so the residue of an incomplete codon is written as a stop.

**What correcting it would change.** The rows would carry no HGVSp. No other column would change.

**The class it produces.** `hgvsp_stop_named_on_incomplete_coding_sequence`: Ensembl VEP's description names the stop alone in one of the four ways and vep-rs writes none, on a row that shows the incomplete coding sequence. The row shows it by `FLAGS` naming `cds_end_NF`, by the term `incomplete_terminal_codon_variant`, or by the stop unchanged at residue 1 under `FLAGS` naming `cds_start_NF`. A row on a transcript that no flag marks incomplete, whose terms hold no such term, stays a difference in the format comparator; the corpus builder also reads the corpus cache.

## 1. A substitution in a partial last codon

**In short.** On GRCh37, a record at chromosome 1, position 3548777, a substitution of one base, changes the last base of the coding sequence of ENST00000424367, a reverse-strand transcript flagged `cds_end_NF` whose coding sequence ends on the first base of a codon. Ensembl VEP writes `<protein>:p.Ter305=`, a stop unchanged at a residue the protein does not have; vep-rs writes no HGVSp.

## 2. A substitution in an incomplete first codon

**In short.** On GRCh37, a record at chromosome 1, position 865692, a substitution of one base, changes the first coding base of ENST00000341065, a forward-strand transcript whose coding sequence opens on the second base of its first codon and is closed by a stop codon. Ensembl VEP writes `<protein>:p.Ter1=`; the reference protein holds an unknown residue there, and vep-rs writes no HGVSp.
