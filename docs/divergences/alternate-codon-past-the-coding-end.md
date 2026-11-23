# The alternate codon that holds bases of the 3' UTR

**The defect in one sentence.** Ensembl VEP 116.3 cuts a row's alternate codon from the coding sequence with the 3' UTR joined on, so where the coding sequence ends inside the codon, or a deletion leaves under three coding bases, that codon holds UTR bases. vep-rs cuts both codons from coding bases alone.

**Why this is a defect.** `Codons` sets the reference codon beside the alternate one, and `Amino_acids` the residues read from each. Ensembl VEP reads the reference codon from the coding sequence, which ends where the transcript's coding region ends. It reads the alternate codon from a longer string. A substitution of one base then turns a codon of two bases into one of three, whose third base the record does not change and the coding sequence does not hold. The row names a residue for the alternate allele at a position where the reference allele has none: `X/K`, an unknown residue replaced by a lysine.

Read from the same longer string, the reference codon would be three bases and name a residue too. No one reading of the transcript gives both halves of the row.

**What vep-rs does instead.** vep-rs cuts the alternate codon from the alternate coding sequence: the coding bases before the record, the allele, and the coding bases after it. A codon window ends where that sequence ends. In a last codon the coding sequence holds part of, both codons are partial and both residues are `X`, which the row writes once. Where a deletion leaves one or two coding bases, the alternate codon is those bases. vep-rs writes its other columns as before.

**Where it is in Ensembl VEP's source** (`ensembl-variation` at release/116):

- `Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm`, lines 2394-2439, `_get_alternate_cds`: lines 2409-2421 join the coding bases before the record, the allele and the coding bases after it. Line 2422 hands the result to `_trim_incomplete_codon`, and lines 2427-2433 join the 3' UTR onto what comes back.
- Lines 2529-2540, `_trim_incomplete_codon`: line 2537 assigns the trimmed length where a comparison was meant. A sequence of three bases or more comes back whole, as in the first shape below, and a shorter one comes back empty, as in the second.
- Lines 798-888, `codon`: line 868 takes the reference allele's codons from the coding sequence alone and the alternate allele's from the joined string, and line 877 cuts each window. The reference window stops at the coding end; the alternate window runs on into the UTR.
- Lines 692-787, `peptide`: lines 724-737 translate the window's whole codons, and lines 774-776 append `X` for the bases left over.
- `Bio/EnsEMBL/Variation/Utils/VariationEffect.pm`, lines 1211-1232, `stop_gained`: line 1228 writes the term where the alternate residues hold a stop and the reference residues hold none.

**What correcting it would change.** The alternate half of `Codons` would lose its UTR bases on every such row, and `Amino_acids` would read `X` wherever it names a residue for the alternate allele alone. A row whose UTR base completes a stop codon would lose `stop_gained`. Positions and every other term would stay.

**The class it produces.** `alternate_codon_past_the_coding_end` reads `Codons` and `Amino_acids` on a `Transcript` row, in two shapes. In the first, the two engines' reference codons are equal and no whole number of codons, and Ensembl VEP's alternate codon is vep-rs's with one or two lower-case bases behind it, no more than the last codon lacks. Each engine's alternate residue is then the translation of that engine's alternate codon under the standard genetic code.

In the second shape, the row's protein positions open at residue 1, the reference codons are whole, and vep-rs's alternate codon is their one or two lower-case bases, where Ensembl VEP's is as many lower-case bases that spell something else. The format comparator applies the class once per row and compares every other field.

Where a transcript of the second shape has no 3' UTR, Ensembl VEP writes a dash for the alternate codon and residue, and vep-rs the bases left and `X`. The class does not read that row, which stays a difference, and no record of the datasets below carries it.

The class does not read the terms. Where a UTR base completes a stop codon, Ensembl VEP's row carries `stop_gained` and vep-rs's does not, and that row stays a difference. No record of the datasets below carries such an alternate; section 1 shows a made-up one.

Under `--hgvs`, Ensembl VEP builds the HGVSp of such a substitution from the same alternate residue: it writes a stop at the last codon's position as replaced by that residue. That notation belongs to another entry, [`hgvsp_stop_named_on_incomplete_coding_sequence`](hgvsp-stop-named-on-incomplete-coding-sequence.md), under which vep-rs writes no HGVSp.

The six whole-genome datasets hold 51 rows of the class. Genome-wide ClinVar holds 5 on GRCh37 and 20 on GRCh38: 11 of the first shape, and 14 on the deletion of section 2. The four population datasets hold 26, all of the first shape: 1000 Genomes Phase 3 4, gnomAD v2.1.1 7, 1000 Genomes high coverage 3 and gnomAD v4.1 12. The chromosome 21 cuts hold none. No chromosome 21 transcript of the release 116 caches has a coding sequence that ends inside a codon before a 3' UTR, and no record of the cuts leaves one or two bases of a coding sequence. The golden corpora `GRCh37` and `GRCh38` carry 3 and 15 rows, documented in their manifests.

## 1. A substitution in a last codon of two bases

**In short.** On GRCh37, the coding sequence of ENST00000434274, a forward-strand transcript, holds 62 bases, cDNA 52 to 113: twenty codons and two bases of a twenty-first. Its 3' UTR follows in the same exon. ClinVar variation 4047389, a substitution at `2:109276125`, changes CDS base 61, the first base of that codon. Ensembl VEP writes an alternate codon of three bases, the third the UTR's first base at `2:109276127`, and a lysine read from it. vep-rs writes the two coding bases and no residue.

| Engine | `Amino_acids` | `Codons` |
| --- | --- | --- |
| Ensembl VEP 116.3 | `X/K` | `Ga/Aag` |
| vep-rs | `X` | `Ga/Aa` |

**The two other substitutions of the base.** Neither stands in a measured dataset; the reproduction writes both. One gives `X/Q` and `Ga/Cag` in Ensembl VEP, and `X` and `Ga/Ca` in vep-rs. The other makes a stop codon of the alternate base, the second coding base and the UTR base. Ensembl VEP writes `X/*` and `Ga/Tag` and adds `stop_gained` before `incomplete_terminal_codon_variant`. vep-rs writes `X` and `Ga/Ta` with the terms of the other two rows.

Both engines write `incomplete_terminal_codon_variant,coding_sequence_variant,NMD_transcript_variant`, CDS position 61 and protein position 21 on the three rows, but for that one term.

## 2. A deletion that leaves one coding base

**In short.** On GRCh37, ClinVar variation 4533193 deletes the 116 bases `3:37035039-37035154`. Three transcripts of one gene, ENST00000442249, ENST00000454028 and ENST00000457004, each hold a coding sequence of 117 bases: those 116 in one exon and the last in the next. The deletion leaves that last base, so the alternate coding sequence is one base long. Ensembl VEP empties it and writes the first base of each transcript's 3' UTR as the alternate codon. vep-rs writes the base that is left.

| Transcript | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| ENST00000442249 | `a` | `a` |
| ENST00000454028 | `t` | `a` |
| ENST00000457004 | `t` | `a` |

The table gives the alternate half of `Codons`. The reference half, the 117 bases with the deleted ones in upper case, is equal on both engines, and so are `Amino_acids`, CDS position 1-116, protein position 1-39 and the terms. The UTR of ENST00000442249 opens with the base the deletion leaves, so its row agrees. On GRCh38 the same record reaches thirteen such transcripts and twelve rows differ.

The entry's reproduction runs the two records and the two made-up substitutions on the cache of the `GRCh37` corpus. The commands of [Reproducing these examples](reproducing.md#reproducing-these-examples) print the same rows; `grep -E 'ENST00000434274|ENST000004(42249|54028|57004)'` selects them.
