// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Standard genetic code (NCBI table 1) codon translation and amino acid utilities.
//!
//! Provides [`translate_codon`] to convert a 3-base DNA codon to a single-letter
//! amino acid, [`complement_base`] / [`reverse_complement`] for strand conversion,
//! and [`aa_three_letter`] for IUPAC three-letter amino acid names.

use std::collections::HashMap;
use std::sync::LazyLock;

/// Standard genetic code (NCBI translation table 1).
/// Maps DNA codon triplets to single-letter amino acid codes.
/// '*' represents a stop codon.
static CODON_TABLE: LazyLock<HashMap<&'static [u8; 3], u8>> = LazyLock::new(|| {
    let mut m = HashMap::new();
    m.insert(b"TTT", b'F');
    m.insert(b"TTC", b'F');
    m.insert(b"TTA", b'L');
    m.insert(b"TTG", b'L');
    m.insert(b"CTT", b'L');
    m.insert(b"CTC", b'L');
    m.insert(b"CTA", b'L');
    m.insert(b"CTG", b'L');
    m.insert(b"ATT", b'I');
    m.insert(b"ATC", b'I');
    m.insert(b"ATA", b'I');
    m.insert(b"ATG", b'M');
    m.insert(b"GTT", b'V');
    m.insert(b"GTC", b'V');
    m.insert(b"GTA", b'V');
    m.insert(b"GTG", b'V');
    m.insert(b"TCT", b'S');
    m.insert(b"TCC", b'S');
    m.insert(b"TCA", b'S');
    m.insert(b"TCG", b'S');
    m.insert(b"AGT", b'S');
    m.insert(b"AGC", b'S');
    m.insert(b"CCT", b'P');
    m.insert(b"CCC", b'P');
    m.insert(b"CCA", b'P');
    m.insert(b"CCG", b'P');
    m.insert(b"ACT", b'T');
    m.insert(b"ACC", b'T');
    m.insert(b"ACA", b'T');
    m.insert(b"ACG", b'T');
    m.insert(b"GCT", b'A');
    m.insert(b"GCC", b'A');
    m.insert(b"GCA", b'A');
    m.insert(b"GCG", b'A');
    m.insert(b"TAT", b'Y');
    m.insert(b"TAC", b'Y');
    m.insert(b"TAA", b'*');
    m.insert(b"TAG", b'*');
    m.insert(b"TGA", b'*');
    m.insert(b"CAT", b'H');
    m.insert(b"CAC", b'H');
    m.insert(b"CAA", b'Q');
    m.insert(b"CAG", b'Q');
    m.insert(b"AAT", b'N');
    m.insert(b"AAC", b'N');
    m.insert(b"AAA", b'K');
    m.insert(b"AAG", b'K');
    m.insert(b"GAT", b'D');
    m.insert(b"GAC", b'D');
    m.insert(b"GAA", b'E');
    m.insert(b"GAG", b'E');
    m.insert(b"TGT", b'C');
    m.insert(b"TGC", b'C');
    m.insert(b"TGG", b'W');
    m.insert(b"CGT", b'R');
    m.insert(b"CGC", b'R');
    m.insert(b"CGA", b'R');
    m.insert(b"CGG", b'R');
    m.insert(b"AGA", b'R');
    m.insert(b"AGG", b'R');
    m.insert(b"GGT", b'G');
    m.insert(b"GGC", b'G');
    m.insert(b"GGA", b'G');
    m.insert(b"GGG", b'G');
    m
});

/// Three-letter amino acid codes.
static AA_THREE_LETTER: LazyLock<HashMap<u8, &'static str>> = LazyLock::new(|| {
    let mut m = HashMap::new();
    m.insert(b'A', "Ala");
    m.insert(b'R', "Arg");
    m.insert(b'N', "Asn");
    m.insert(b'D', "Asp");
    m.insert(b'C', "Cys");
    m.insert(b'E', "Glu");
    m.insert(b'Q', "Gln");
    m.insert(b'G', "Gly");
    m.insert(b'H', "His");
    m.insert(b'I', "Ile");
    m.insert(b'L', "Leu");
    m.insert(b'K', "Lys");
    m.insert(b'M', "Met");
    m.insert(b'F', "Phe");
    m.insert(b'P', "Pro");
    m.insert(b'S', "Ser");
    m.insert(b'T', "Thr");
    m.insert(b'W', "Trp");
    m.insert(b'Y', "Tyr");
    m.insert(b'V', "Val");
    m.insert(b'*', "Ter");
    m.insert(b'X', "Xaa");
    m
});

/// Translate a DNA codon (3 uppercase bases) to a single-letter amino acid using
/// the standard genetic code (NCBI translation table 1).
///
/// Returns 'X' for unknown codons (e.g., containing N).
///
/// Transcripts on a non-standard genetic code must go through
/// [`translate_codon_with_table`] instead; see that function for why.
pub fn translate_codon(codon: &[u8]) -> u8 {
    translate_codon_with_table(codon, 1)
}

/// Translate a DNA codon under the given NCBI translation table.
///
/// Perl VEP selects the table per transcript from the slice's `codon_table`
/// attribute (`BaseTranscriptVariation::_codon_table`, defaulting to 1) and
/// threads it into every translation. Human vertebrate-mitochondrial transcripts
/// carry table 2, which differs from table 1 in exactly four codons: `ATA` codes
/// Met rather than Ile, `TGA` codes Trp rather than a stop, and `AGA` / `AGG`
/// are stops rather than Arg. Translating MT transcripts under table 1 silently
/// mis-calls the consequence for variants touching those codons.
///
/// Any table other than 2 falls back to table 1, matching Perl's default.
pub fn translate_codon_with_table(codon: &[u8], codon_table: u8) -> u8 {
    if codon.len() != 3 {
        return b'X';
    }
    let upper: [u8; 3] = [
        codon[0].to_ascii_uppercase(),
        codon[1].to_ascii_uppercase(),
        codon[2].to_ascii_uppercase(),
    ];
    if codon_table == 2 {
        match &upper {
            b"ATA" => return b'M',
            b"TGA" => return b'W',
            b"AGA" | b"AGG" => return b'*',
            _ => {}
        }
    }
    CODON_TABLE.get(&upper).copied().unwrap_or(b'X')
}

/// Get three-letter amino acid code from single-letter code.
pub fn aa_three_letter(aa: u8) -> &'static str {
    AA_THREE_LETTER.get(&aa).unwrap_or(&"Xaa")
}

/// Reverse complement a DNA sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// Complement a single DNA base.
pub fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'a' => b't',
        b'T' => b'A',
        b't' => b'a',
        b'G' => b'C',
        b'g' => b'c',
        b'C' => b'G',
        b'c' => b'g',
        b'N' => b'N',
        b'n' => b'n',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_codon() {
        assert_eq!(translate_codon(b"ATG"), b'M');
        assert_eq!(translate_codon(b"TAA"), b'*');
        assert_eq!(translate_codon(b"TAG"), b'*');
        assert_eq!(translate_codon(b"TGA"), b'*');
        assert_eq!(translate_codon(b"GCT"), b'A');
        assert_eq!(translate_codon(b"TGG"), b'W');
        assert_eq!(translate_codon(b"NNN"), b'X');
    }

    #[test]
    fn test_translate_codon_case_insensitive() {
        assert_eq!(translate_codon(b"atg"), b'M');
        assert_eq!(translate_codon(b"Atg"), b'M');
    }

    #[test]
    fn test_translate_codon_table2_diverging_codons() {
        // Each divergent codon is checked against table 1 in the same pair so a
        // change in either direction is visible.
        assert_eq!(translate_codon_with_table(b"ATA", 1), b'I');
        assert_eq!(translate_codon_with_table(b"ATA", 2), b'M');

        assert_eq!(translate_codon_with_table(b"TGA", 1), b'*');
        assert_eq!(translate_codon_with_table(b"TGA", 2), b'W');

        assert_eq!(translate_codon_with_table(b"AGA", 1), b'R');
        assert_eq!(translate_codon_with_table(b"AGA", 2), b'*');

        assert_eq!(translate_codon_with_table(b"AGG", 1), b'R');
        assert_eq!(translate_codon_with_table(b"AGG", 2), b'*');
    }

    #[test]
    fn test_translate_codon_table2_agrees_with_table1_elsewhere() {
        // Every codon outside the four divergent ones must translate identically
        // under both tables, so enabling table 2 cannot perturb unrelated calls.
        let bases = b"ACGT";
        let divergent: [&[u8; 3]; 4] = [b"ATA", b"TGA", b"AGA", b"AGG"];
        let mut compared = 0;
        for &a in bases {
            for &b in bases {
                for &c in bases {
                    let codon = [a, b, c];
                    if divergent.iter().any(|d| d.as_slice() == codon) {
                        continue;
                    }
                    assert_eq!(
                        translate_codon_with_table(&codon, 1),
                        translate_codon_with_table(&codon, 2),
                        "table 1 and 2 must agree on {}",
                        std::str::from_utf8(&codon).unwrap()
                    );
                    compared += 1;
                }
            }
        }
        assert_eq!(
            compared, 60,
            "expected 64 codons minus the 4 divergent ones"
        );
    }

    #[test]
    fn test_translate_codon_unknown_table_falls_back_to_table1() {
        // Perl defaults to table 1 for any transcript without a codon_table
        // attribute; an unrecognised table id must not silently change calls.
        for table in [0u8, 1, 3, 11, 255] {
            assert_eq!(translate_codon_with_table(b"ATA", table), b'I');
            assert_eq!(translate_codon_with_table(b"TGA", table), b'*');
        }
    }

    #[test]
    fn test_translate_codon_defaults_to_table1() {
        assert_eq!(translate_codon(b"ATA"), b'I');
        assert_eq!(translate_codon(b"TGA"), b'*');
        assert_eq!(translate_codon(b"AGA"), b'R');
        assert_eq!(translate_codon(b"AGG"), b'R');
    }

    #[test]
    fn test_aa_three_letter() {
        assert_eq!(aa_three_letter(b'M'), "Met");
        assert_eq!(aa_three_letter(b'*'), "Ter");
        assert_eq!(aa_three_letter(b'X'), "Xaa");
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_all_64_codons_translate() {
        let bases = *b"ATGC";
        for &b1 in &bases {
            for &b2 in &bases {
                for &b3 in &bases {
                    let codon = [b1, b2, b3];
                    let aa = translate_codon(&codon);
                    assert_ne!(
                        aa,
                        b'X',
                        "Standard codon {:?} should translate",
                        std::str::from_utf8(&codon).unwrap()
                    );
                }
            }
        }
    }
}
