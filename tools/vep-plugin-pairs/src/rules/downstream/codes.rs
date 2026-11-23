// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The genetic codes the Downstream plugin translates with and the matrix the Blosum62
//! plugin scores with.

/// The genetic codes Bio::Tools::CodonTable translates with, by table id, as
/// `crates/vep-core/src/codon.rs` lists them: each codon's amino acid, the codons in the
/// order TTT TTC TTA TTG TCT ... GGG. Ids 7, 8 and 17 to 20 hold no table there and ids from
/// 25 on are unknown to it.
const CODES: [(u8, &[u8; 64]); 18] = [
    (
        1,
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        2,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    ),
    (
        3,
        b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        4,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        5,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        6,
        b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        9,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        10,
        b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        11,
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        12,
        b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        13,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    ),
    (
        14,
        b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        15,
        b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        16,
        b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        21,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        22,
        b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        23,
        b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        24,
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
];

/// A genetic code: the amino acid of each codon, a stop `*`.
pub type GeneticCode = [u8; 64];

/// The genetic code of the table `table`, `None` for an id that holds no table.
pub fn genetic_code(table: u8) -> Option<&'static GeneticCode> {
    CODES
        .iter()
        .find(|(id, _)| *id == table)
        .map(|(_, code)| *code)
}

/// The place of a base in the order T, C, A, G of the codes' codons.
fn base(letter: char) -> Option<usize> {
    match letter {
        'T' => Some(0),
        'C' => Some(1),
        'A' => Some(2),
        'G' => Some(3),
        _ => None,
    }
}

/// The bases a nucleotide letter stands for, none for a letter that is no nucleotide's.
fn readings(letter: char) -> &'static [usize] {
    const T: usize = 0;
    const C: usize = 1;
    const A: usize = 2;
    const G: usize = 3;
    match letter {
        'A' => &[A],
        'C' => &[C],
        'G' => &[G],
        'T' | 'U' => &[T],
        'R' => &[A, G],
        'Y' => &[C, T],
        'S' => &[C, G],
        'W' => &[A, T],
        'K' => &[G, T],
        'M' => &[A, C],
        'B' => &[C, G, T],
        'D' => &[A, G, T],
        'H' => &[A, C, T],
        'V' => &[A, C, G],
        'N' | 'X' => &[A, C, G, T],
        _ => &[],
    }
}

/// The residue Bio::Tools::CodonTable writes for a codon: its amino acid; for a codon with
/// an ambiguous letter the amino acids its readings code, one amino acid as itself,
/// aspartate with asparagine `B`, glutamate with glutamine `Z` and anything else `X`, as is
/// a codon with a letter that is no nucleotide's.
fn residue(code: &GeneticCode, codon: [char; 3]) -> u8 {
    if let (Some(first), Some(second), Some(third)) =
        (base(codon[0]), base(codon[1]), base(codon[2]))
    {
        return code[16 * first + 4 * second + third];
    }
    let mut coded: Vec<u8> = Vec::new();
    for first in readings(codon[0]) {
        for second in readings(codon[1]) {
            for third in readings(codon[2]) {
                let amino_acid = code[16 * first + 4 * second + third];
                if !coded.contains(&amino_acid) {
                    coded.push(amino_acid);
                }
            }
        }
    }
    coded.sort_unstable();
    match coded[..] {
        [only] => only,
        [b'D', b'N'] => b'B',
        [b'E', b'Q'] => b'Z',
        _ => b'X',
    }
}

/// Hands `take` each residue of the translation of `sequence` under `code` up to its first
/// stop, as `Bio::Tools::CodonTable::translate` reads it for Downstream.pm:139-142: the
/// whole codons in order, a trailing one or two bases dropped.
fn residues(sequence: impl Iterator<Item = char>, code: &GeneticCode, mut take: impl FnMut(u8)) {
    let mut letters = sequence;
    while let (Some(first), Some(second), Some(third)) =
        (letters.next(), letters.next(), letters.next())
    {
        match residue(code, [first, second, third]) {
            b'*' => break,
            amino_acid => take(amino_acid),
        }
    }
}

/// The translation of `sequence` under `code` up to its first stop (`residues`).
pub fn translate_until_stop(sequence: impl Iterator<Item = char>, code: &GeneticCode) -> String {
    let mut translation = String::new();
    residues(sequence, code, |amino_acid| {
        translation.push(char::from(amino_acid))
    });
    translation
}

/// The number of residues `translate_until_stop` gives.
pub fn translated_length(sequence: impl Iterator<Item = char>, code: &GeneticCode) -> usize {
    let mut length = 0;
    residues(sequence, code, |_| length += 1);
    length
}

/// The residue `code` reads from the three letters `codon`, `None` for a codon that is not
/// three of `T`, `C`, `A` and `G`.
pub fn codon_residue(code: &GeneticCode, codon: &str) -> Option<char> {
    let mut letters = codon.chars();
    let (first, second, third) = (
        base(letters.next()?)?,
        base(letters.next()?)?,
        base(letters.next()?)?,
    );
    if letters.next().is_some() {
        return None;
    }
    Some(char::from(code[16 * first + 4 * second + third]))
}

/// The residue `code` reads from a codon of nucleotide letters (`residue`): the amino acid
/// every reading of the codon shares, `B` for aspartate or asparagine, `Z` for glutamate or
/// glutamine and `X` otherwise. `None` for a codon with a letter that names no base.
pub fn nucleotide_codon_residue(code: &GeneticCode, codon: [char; 3]) -> Option<char> {
    codon
        .iter()
        .all(|letter| !readings(*letter).is_empty())
        .then(|| char::from(residue(code, codon)))
}

const BLOSUM62_RESIDUES: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";
/// The BLOSUM62 matrix of Blosum62.pm:50-71, its rows and columns in the order of
/// `BLOSUM62_RESIDUES`.
const BLOSUM62: [[i8; 20]; 20] = [
    [
        4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0,
    ],
    [
        -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3,
    ],
    [
        -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3,
    ],
    [
        -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3,
    ],
    [
        0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
    ],
    [
        -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2,
    ],
    [
        -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2,
    ],
    [
        0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3,
    ],
    [
        -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3,
    ],
    [
        -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3,
    ],
    [
        -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1,
    ],
    [
        -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2,
    ],
    [
        -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1,
    ],
    [
        -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1,
    ],
    [
        -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2,
    ],
    [
        1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2,
    ],
    [
        0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0,
    ],
    [
        -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3,
    ],
    [
        -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1,
    ],
    [
        0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4,
    ],
];

/// The matrix cell of a residue pair, `None` for a residue the matrix does not hold.
pub fn blosum62(first: char, second: char) -> Option<i8> {
    let place = |residue: char| {
        BLOSUM62_RESIDUES
            .iter()
            .position(|held| char::from(*held) == residue)
    };
    Some(BLOSUM62[place(first)?][place(second)?])
}
