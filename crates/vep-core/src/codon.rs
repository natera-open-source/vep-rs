// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Codon translation under the genetic codes of BioPerl's `Bio::Tools::CodonTable`,
//! the module Ensembl VEP translates with, and under NCBI's genetic code 26,
//! which that module lacks, and amino acid utilities.
//!
//! Provides [`translate_codon`] (the standard code) and
//! [`translate_codon_with_table`] (a transcript's own code) to convert a 3-base
//! DNA codon to a single-letter amino acid, the start and stop predicates
//! [`is_start_codon`] and [`is_ter_codon`] with [`known_table`], the set of codes
//! read here, and [`module_table`], the table the module itself reads an id
//! under, [`complement_base`] / [`reverse_complement`] for strand
//! conversion, and [`aa_three_letter`] for IUPAC three-letter amino acid names.

use std::collections::HashMap;
use std::sync::LazyLock;

/// One genetic code as `Bio::Tools::CodonTable` holds one: the amino acid of
/// each codon and the start marks (`M` at a start codon, `-` elsewhere), both
/// in the module's codon order `TTT TTC TTA TTG TCT ... GGG`, the first base
/// the most significant and `T`, `C`, `A`, `G` counted 0 to 3.
#[derive(Clone, Copy)]
struct GeneticCode {
    amino_acids: &'static [u8; 64],
    starts: &'static [u8; 64],
}

const fn code(amino_acids: &'static [u8; 64], starts: &'static [u8; 64]) -> Option<GeneticCode> {
    Some(GeneticCode {
        amino_acids,
        starts,
    })
}

/// The genetic codes of BioPerl 1.006924 by id, with the module's translation
/// and start strings. Ids 7, 8 and 17 to 20 hold no table there (the module's
/// entries are the two-character text `''`, read as an apostrophe for `TTT` and
/// `TTC` and as an empty text for every other codon, marking no start or stop)
/// and ids from 25 on are unknown to it (it warns and translates under its last
/// table): both kinds are `None` here. Table 24 is the module's own `Strict`
/// code, the standard translation with `ATG` as its one start codon, not the
/// NCBI table of that number.
const CODES: [Option<GeneticCode>; 25] = [
    None,
    // 1 Standard
    code(
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"---M---------------M---------------M----------------------------",
    ),
    // 2 Vertebrate Mitochondrial
    code(
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        b"--------------------------------MMMM---------------M------------",
    ),
    // 3 Yeast Mitochondrial
    code(
        b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"----------------------------------MM----------------------------",
    ),
    // 4 Mold, Protozoan and Coelenterate Mitochondrial; Mycoplasma and Spiroplasma
    code(
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"--MM---------------M------------MMMM---------------M------------",
    ),
    // 5 Invertebrate Mitochondrial
    code(
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        b"---M----------------------------MMMM---------------M------------",
    ),
    // 6 Ciliate, Dasycladacean and Hexamita Nuclear
    code(
        b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"-----------------------------------M----------------------------",
    ),
    None,
    None,
    // 9 Echinoderm Mitochondrial
    code(
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        b"-----------------------------------M----------------------------",
    ),
    // 10 Euplotid Nuclear
    code(
        b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"-----------------------------------M----------------------------",
    ),
    // 11 Bacterial
    code(
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"---M---------------M------------MMMM---------------M------------",
    ),
    // 12 Alternative Yeast Nuclear
    code(
        b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"-------------------M---------------M----------------------------",
    ),
    // 13 Ascidian Mitochondrial
    code(
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        b"-----------------------------------M----------------------------",
    ),
    // 14 Flatworm Mitochondrial
    code(
        b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        b"-----------------------------------M----------------------------",
    ),
    // 15 Blepharisma Nuclear
    code(
        b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"-----------------------------------M----------------------------",
    ),
    // 16 Chlorophycean Mitochondrial
    code(
        b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"-----------------------------------M----------------------------",
    ),
    None,
    None,
    None,
    None,
    // 21 Trematode Mitochondrial
    code(
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        b"-----------------------------------M---------------M------------",
    ),
    // 22 Scenedesmus obliquus Mitochondrial
    code(
        b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"-----------------------------------M----------------------------",
    ),
    // 23 Thraustochytrium Mitochondrial
    code(
        b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"--------------------------------M--M---------------M------------",
    ),
    // 24 Strict
    code(
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"-----------------------------------M----------------------------",
    ),
];

/// The id of the module's last table, its `Strict` code.
const MODULE_LAST_TABLE: u8 = 24;

/// The id of the nuclear genetic code of Pachysolen tannophilus, which the
/// module holds no table for.
const PACHYSOLEN_NUCLEAR: u8 = 26;

/// Genetic code 26 as NCBI's genetic code table gives it: `CTG` codes alanine
/// where the standard code reads leucine, and `CTG` and `ATG` are its start
/// codons. It stands outside `CODES`, which holds the module's tables alone.
const PACHYSOLEN_NUCLEAR_CODE: GeneticCode = GeneticCode {
    amino_acids: b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    starts: b"-------------------M---------------M----------------------------",
};

/// The genetic code of `table`: a table of the module, or code 26;
/// `Bio::Tools::CodonTable->new(-id => 0)` is table 1, as is a transcript
/// whose sequence region carries no `codon_table` attribute
/// (`BaseTranscriptVariation::_codon_table`).
#[inline]
fn genetic_code(table: u8) -> Option<GeneticCode> {
    if table == PACHYSOLEN_NUCLEAR {
        return Some(PACHYSOLEN_NUCLEAR_CODE);
    }
    let id = if table == 0 { 1 } else { usize::from(table) };
    CODES.get(id).copied().flatten()
}

/// Whether `table` is a genetic code read here: one BioPerl 1.006924, and so
/// Ensembl VEP, translates with (0 to 24 except 7, 8 and 17 to 20), or 26.
pub fn known_table(table: u8) -> bool {
    genetic_code(table).is_some()
}

/// The table `Bio::Tools::CodonTable` reads `table` under: an id from 25 on
/// holds no table there, so `id` warns and falls to id 0
/// (`CodonTable.pm:328-338`), under which `translate` and `_codon_is` read the
/// module's last table (`CodonTable.pm:412-460` and `:735-749`); every other
/// id is read as itself. Ensembl VEP translates a transcript on such an id
/// under the table returned, and a dumper on that module writes the
/// transcript's cached peptide under it.
pub fn module_table(table: u8) -> u8 {
    if table > MODULE_LAST_TABLE {
        MODULE_LAST_TABLE
    } else {
        table
    }
}

/// The module's rank of a base, in either case, `U` read as `T` as its
/// `tr/u/t/` reads it; `None` for any other byte.
#[inline]
fn base_rank(base: u8) -> Option<usize> {
    match base {
        b'T' | b't' | b'U' | b'u' => Some(0),
        b'C' | b'c' => Some(1),
        b'A' | b'a' => Some(2),
        b'G' | b'g' => Some(3),
        _ => None,
    }
}

/// The codon's position in the module's 64-codon strings, `None` for a codon
/// not three bases long or carrying a byte that is no base.
#[inline]
fn codon_index(codon: &[u8]) -> Option<usize> {
    let [first, second, third] = codon else {
        return None;
    };
    Some((base_rank(*first)? << 4) | (base_rank(*second)? << 2) | base_rank(*third)?)
}

/// Translate a DNA codon (3 bases) to a single-letter amino acid using the
/// standard genetic code (NCBI translation table 1).
///
/// Returns 'X' for unknown codons (e.g., containing N).
///
/// Transcripts on a non-standard genetic code must go through
/// [`translate_codon_with_table`] instead; see that function for why.
pub fn translate_codon(codon: &[u8]) -> u8 {
    translate_codon_with_table(codon, 1)
}

/// Translate a DNA codon under the genetic code of the given `Bio::Tools::CodonTable`
/// id, as `Bio::Tools::CodonTable::translate` translates one unambiguous codon.
///
/// Ensembl VEP selects the table per transcript from the slice's `codon_table`
/// attribute (`BaseTranscriptVariation::_codon_table`, defaulting to 1) and
/// threads it into every translation (`TranscriptVariationAllele::peptide`,
/// `Transcript::translate`, the stop-codon predicates of `VariationEffect`).
/// Human vertebrate-mitochondrial transcripts carry table 2, which differs from
/// table 1 in exactly four codons: `ATA` codes Met rather than Ile, `TGA`
/// codes Trp rather than a stop, and `AGA` / `AGG` are stops rather than Arg;
/// an invertebrate's mitochondrial transcripts carry table 5, a yeast's table 3,
/// a plastid's table 11. Translating such a transcript under table 1 silently
/// mis-calls the consequence for variants touching the codons that differ.
///
/// Returns `X` for a codon not three bases long, for a codon with a byte
/// outside `ACGTU` in either case (an `N`, a gap) and under a table that is
/// not known ([`known_table`]); the cache loader refuses such a table before
/// any transcript is annotated under it.
pub fn translate_codon_with_table(codon: &[u8], codon_table: u8) -> u8 {
    match (genetic_code(codon_table), codon_index(codon)) {
        (Some(code), Some(index)) => code.amino_acids[index],
        _ => b'X',
    }
}

/// `Bio::Tools::CodonTable::translate` under `table`, as `Bio::Seq::translate`
/// calls it for a sequence with no frame and no partial-codon reading: the
/// sequence's whole codons in order, a trailing one or two bases dropped; a
/// codon of unambiguous bases by the table, `---` as `-`, and a codon with an
/// ambiguous letter by the amino acids its readings code
/// (`_translate_ambiguous_codon`): one, that one; aspartate and asparagine
/// `B`; glutamate and glutamine `Z`; anything else `X`, as is a codon with a
/// byte that is no nucleotide letter. A stop stays `*`. Under a table that is
/// not known every codon is `X`.
pub fn translate_sequence(seq: &[u8], table: u8) -> Vec<u8> {
    let Some(code) = genetic_code(table) else {
        return vec![b'X'; seq.len() / 3];
    };
    seq.as_chunks::<3>()
        .0
        .iter()
        .map(|codon| match codon_index(codon) {
            Some(index) => code.amino_acids[index],
            None if codon == b"---" => b'-',
            None => translate_ambiguous_codon(codon, code.amino_acids),
        })
        .collect()
}

/// `Bio::Tools::CodonTable::_translate_ambiguous_codon`: the amino acid of a
/// codon that is no unambiguous codon, from the distinct amino acids
/// `amino_acids` gives its readings; `X` for a slice not three bytes long.
fn translate_ambiguous_codon(codon: &[u8], amino_acids: &[u8; 64]) -> u8 {
    let [first, second, third] = codon else {
        return b'X';
    };
    let mut distinct = [0u8; 64];
    let mut count = 0;
    for a in iupac_ranks(*first) {
        for b in iupac_ranks(*second) {
            for c in iupac_ranks(*third) {
                let amino_acid = amino_acids[(a << 4) | (b << 2) | c];
                if !distinct[..count].contains(&amino_acid) {
                    distinct[count] = amino_acid;
                    count += 1;
                }
            }
        }
    }
    match &distinct[..count] {
        [amino_acid] => *amino_acid,
        [b'D', b'N'] | [b'N', b'D'] => b'B',
        [b'E', b'Q'] | [b'Q', b'E'] => b'Z',
        _ => b'X',
    }
}

/// The bases an IUPAC nucleotide letter stands for, as ranks, in either case;
/// empty for a byte that is no nucleotide letter.
fn iupac_ranks(base: u8) -> &'static [usize] {
    match base.to_ascii_uppercase() {
        b'T' | b'U' => &[0],
        b'C' => &[1],
        b'A' => &[2],
        b'G' => &[3],
        b'Y' => &[0, 1],
        b'W' => &[0, 2],
        b'K' => &[0, 3],
        b'M' => &[1, 2],
        b'S' => &[1, 3],
        b'R' => &[2, 3],
        b'H' => &[0, 1, 2],
        b'B' => &[0, 1, 3],
        b'D' => &[0, 2, 3],
        b'V' => &[1, 2, 3],
        b'N' | b'X' => &[0, 1, 2, 3],
        _ => &[],
    }
}

/// `Bio::Tools::CodonTable::_codon_is`: the codon is three bytes long, each a
/// nucleotide letter, and every unambiguous codon it stands for carries `mark`
/// in `marks`. A byte that is no nucleotide letter (a gap, a `.`) is `false`,
/// where the module's loop over no readings is vacuously true.
fn every_reading_marked(codon: &[u8], marks: &[u8; 64], mark: u8) -> bool {
    let [first, second, third] = codon else {
        return false;
    };
    let (firsts, seconds, thirds) = (
        iupac_ranks(*first),
        iupac_ranks(*second),
        iupac_ranks(*third),
    );
    if firsts.is_empty() || seconds.is_empty() || thirds.is_empty() {
        return false;
    }
    firsts.iter().all(|a| {
        seconds.iter().all(|b| {
            thirds
                .iter()
                .all(|c| marks[(a << 4) | (b << 2) | c] == mark)
        })
    })
}

/// `Bio::Tools::CodonTable::is_start_codon` under `table`: every reading of
/// the codon is a start codon of the table (`TTG`, `CTG` and `ATG` under
/// table 1; `ATT`, `ATC`, `ATA`, `ATG` and `GTG` under table 2), so `NTG` is a
/// start under table 11, whose starts include all four `NTG` codons, and not
/// under table 1. `Transcript::translate` writes `M` over the first residue of
/// a peptide whose first codon satisfies it. `false` under a table that is
/// not known.
pub fn is_start_codon(codon: &[u8], table: u8) -> bool {
    genetic_code(table).is_some_and(|code| every_reading_marked(codon, code.starts, b'M'))
}

/// `Bio::Tools::CodonTable::is_ter_codon` under `table`: every reading of the
/// codon is a stop codon of the table (`TAA`, `TAG` and `TGA` under table 1;
/// `TAA`, `TAG`, `AGA` and `AGG` under table 2), so `TAR` is a stop under
/// table 1 and `TRA` is not under table 2, where `TGA` codes Trp.
/// `Transcript::translate` drops the terminal residue of a peptide whose last
/// codon satisfies it. `false` under a table that is not known.
pub fn is_ter_codon(codon: &[u8], table: u8) -> bool {
    genetic_code(table).is_some_and(|code| every_reading_marked(codon, code.amino_acids, b'*'))
}

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

/// Get three-letter amino acid code from single-letter code.
pub fn aa_three_letter(aa: u8) -> &'static str {
    AA_THREE_LETTER.get(&aa).unwrap_or(&"Xaa")
}

/// Reverse complement a DNA sequence, as Ensembl's `reverse_comp` does
/// (`Bio::EnsEMBL::Utils::Sequence`): every byte through [`complement_base`],
/// in reverse order.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// Whether an allele is DNA by Ensembl's `seq_is_dna`
/// (VariationFeatureOverlapAllele.pm:315-321): one or more nucleotide or IUPAC
/// letters or `-`, in either case. `feature_seq` reverse-complements only such an
/// allele; a symbolic allele (`<DEL>`, a breakend, `*`) is read as written.
pub fn seq_is_dna(allele: &[u8]) -> bool {
    !allele.is_empty()
        && allele
            .iter()
            .all(|b| b"ACGTUMRWSYKVHDBXN-".contains(&b.to_ascii_uppercase()))
}

/// Complement a single DNA base, in either case: the four nucleotides and the
/// IUPAC letters with a complement (`R` and `Y`, `M` and `K`, `B` and `V`, `D`
/// and `H`), as Ensembl's `reverse_comp` complements them; `S`, `W`, `N`, `X`,
/// `-` and every other byte are their own complement.
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
        b'R' => b'Y',
        b'r' => b'y',
        b'Y' => b'R',
        b'y' => b'r',
        b'M' => b'K',
        b'm' => b'k',
        b'K' => b'M',
        b'k' => b'm',
        b'B' => b'V',
        b'b' => b'v',
        b'V' => b'B',
        b'v' => b'b',
        b'D' => b'H',
        b'd' => b'h',
        b'H' => b'D',
        b'h' => b'd',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The 64 codons in the module's order, `TTT` first and `GGG` last.
    fn codons_in_table_order() -> Vec<[u8; 3]> {
        let mut codons = Vec::with_capacity(64);
        for first in b"TCAG" {
            for second in b"TCAG" {
                for third in b"TCAG" {
                    codons.push([*first, *second, *third]);
                }
            }
        }
        codons
    }

    fn translation_string(table: u8) -> String {
        codons_in_table_order()
            .iter()
            .map(|codon| translate_codon_with_table(codon, table) as char)
            .collect()
    }

    fn codons_where(table: u8, predicate: fn(&[u8], u8) -> bool) -> Vec<String> {
        codons_in_table_order()
            .iter()
            .filter(|codon| predicate(codon.as_slice(), table))
            .map(|codon| String::from_utf8(codon.to_vec()).unwrap())
            .collect()
    }

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

    /// `Bio::Tools::CodonTable::translate` lower-cases a codon and reads `u`
    /// as `t`, so an RNA spelling translates; `translate_codon` returns `X` for
    /// a wrong length and for a byte that is no base (an ambiguity letter, a gap).
    #[test]
    fn rna_bases_translate_and_other_bytes_are_x() {
        assert_eq!(translate_codon(b"AUG"), b'M');
        assert_eq!(translate_codon(b"uaa"), b'*');
        for refused in [
            &b"AT"[..],
            b"ATGA",
            b"",
            b"---",
            b"A-G",
            b"AT.",
            b"NTG",
            b"RAY",
        ] {
            assert_eq!(
                translate_codon(refused),
                b'X',
                "{:?}",
                String::from_utf8_lossy(refused)
            );
        }
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

    /// Every genetic code BioPerl 1.006924 ships, the `Bio::Tools::CodonTable`
    /// of Ensembl VEP release 116.3: the 64-codon translation in the module's
    /// `TCAG` order, its start codons (`is_start_codon`) and its stop codons
    /// (`is_ter_codon`). Table 11 translates as table 1 and differs in its
    /// starts alone; table 24 is the module's `Strict` code, table 1 with `ATG`
    /// as its one start.
    #[test]
    fn every_bioperl_table_translates_and_starts_as_bioperl_1_6_924() {
        let standard = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        let expected: [(u8, &str, &str, &str); 18] = [
            (1, standard, "TTG CTG ATG", "TAA TAG TGA"),
            (
                2,
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
                "ATT ATC ATA ATG GTG",
                "TAA TAG AGA AGG",
            ),
            (
                3,
                "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "ATA ATG",
                "TAA TAG",
            ),
            (
                4,
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "TTA TTG CTG ATT ATC ATA ATG GTG",
                "TAA TAG",
            ),
            (
                5,
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
                "TTG ATT ATC ATA ATG GTG",
                "TAA TAG",
            ),
            (
                6,
                "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "ATG",
                "TGA",
            ),
            (
                9,
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "ATG",
                "TAA TAG",
            ),
            (
                10,
                "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "ATG",
                "TAA TAG",
            ),
            (11, standard, "TTG CTG ATT ATC ATA ATG GTG", "TAA TAG TGA"),
            (
                12,
                "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "CTG ATG",
                "TAA TAG TGA",
            ),
            (
                13,
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
                "ATG",
                "TAA TAG",
            ),
            (
                14,
                "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "ATG",
                "TAG",
            ),
            (
                15,
                "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "ATG",
                "TAA TGA",
            ),
            (
                16,
                "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "ATG",
                "TAA TGA",
            ),
            (
                21,
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                "ATG GTG",
                "TAA TAG",
            ),
            (
                22,
                "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "ATG",
                "TCA TAA TGA",
            ),
            (
                23,
                "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                "ATT ATG GTG",
                "TTA TAA TAG TGA",
            ),
            (24, standard, "ATG", "TAA TAG TGA"),
        ];
        for (table, translation, starts, stops) in expected {
            assert!(known_table(table), "table {table}");
            assert_eq!(translation_string(table), translation, "table {table}");
            assert_eq!(
                codons_where(table, is_start_codon).join(" "),
                starts,
                "starts of table {table}"
            );
            assert_eq!(
                codons_where(table, is_ter_codon).join(" "),
                stops,
                "stops of table {table}"
            );
        }
        let differing = |table: u8| -> Vec<String> {
            codons_in_table_order()
                .iter()
                .filter(|c| translate_codon_with_table(*c, table) != translate_codon(*c))
                .map(|c| {
                    format!(
                        "{}>{}",
                        std::str::from_utf8(c).unwrap(),
                        translate_codon_with_table(c, table) as char
                    )
                })
                .collect()
        };
        assert_eq!(differing(2), ["TGA>W", "ATA>M", "AGA>*", "AGG>*"]);
        assert_eq!(
            differing(3),
            ["TGA>W", "CTT>T", "CTC>T", "CTA>T", "CTG>T", "ATA>M"]
        );
        assert_eq!(differing(5), ["TGA>W", "ATA>M", "AGA>S", "AGG>S"]);
        assert!(differing(11).is_empty());
        // `Bio::Tools::CodonTable->new(-id => 0)` is table 1.
        assert_eq!(translation_string(0), standard);
        assert_eq!(codons_where(0, is_start_codon), ["TTG", "CTG", "ATG"]);
    }

    /// Ids 7, 8 and 17 to 20 hold no table in BioPerl 1.006924 (a codon reads
    /// as an apostrophe or as nothing there) and 25 and above are unknown to it
    /// (it warns and translates under its last table): none is a table Ensembl
    /// VEP translates with, and none but 26 is known. Under an unknown one
    /// every codon is `X` and no codon is a start or a stop. Table 24, the last
    /// one the module ships, is known, and so is 26, NCBI's genetic code of
    /// that number, under which `CTG` codes alanine and starts a translation.
    #[test]
    fn an_unknown_table_is_refused_not_translated_empty() {
        for table in [7u8, 8, 17, 18, 19, 20, 25, 27, 100, 255] {
            assert!(!known_table(table), "table {table}");
            assert_eq!(translate_codon_with_table(b"ATG", table), b'X', "{table}");
            assert_eq!(translate_codon_with_table(b"TAA", table), b'X', "{table}");
            assert!(!is_start_codon(b"ATG", table), "{table}");
            assert!(!is_ter_codon(b"TAA", table), "{table}");
        }
        for table in (0u8..=6).chain(9..=16).chain(21..=24) {
            assert!(known_table(table), "table {table}");
            assert_eq!(module_table(table), table, "table {table}");
        }
        assert!(known_table(24));
        assert!(!known_table(25));
        assert!(known_table(26));
        assert!(!known_table(27));
        assert_eq!(module_table(26), 24);
        assert_eq!(translate_codon_with_table(b"CTG", 26), b'A');
        assert_eq!(translate_codon_with_table(b"CTG", module_table(26)), b'L');
        assert_eq!(codons_where(26, is_start_codon), ["CTG", "ATG"]);
        assert_eq!(codons_where(26, is_ter_codon), ["TAA", "TAG", "TGA"]);
    }

    /// `_codon_is` reads an ambiguous codon through every codon it stands for
    /// and answers true only when all of them carry the mark: `TAR` is a stop
    /// under table 1 (`TAA`, `TAG`), `TRA` is not under table 2 (`TGA` is Trp),
    /// `NTG` is a start under table 11 (all four `NTG` codons start there) and
    /// not under table 1, `ATN` and `ATH` are starts under table 2; the codon
    /// is read in either case and with `U` as `T`; a wrong length is false.
    /// A byte that is no nucleotide letter is false, where the module's loop
    /// over no readings answers true for `---`.
    #[test]
    fn start_and_stop_predicates_read_ambiguous_codons_as_bioperl() {
        assert!(is_ter_codon(b"TAR", 1));
        assert!(is_ter_codon(b"TRA", 1));
        assert!(!is_ter_codon(b"TRA", 2));
        assert!(!is_ter_codon(b"NNN", 1));
        assert!(!is_ter_codon(b"TAX", 1));
        assert!(is_ter_codon(b"tag", 1));
        assert!(is_ter_codon(b"UAG", 1));
        assert!(!is_ter_codon(b"TA", 1));
        assert!(!is_ter_codon(b"TAGA", 1));
        assert!(is_start_codon(b"NTG", 11));
        assert!(!is_start_codon(b"NTG", 1));
        assert!(is_start_codon(b"HTG", 11));
        assert!(is_start_codon(b"VTG", 11));
        assert!(is_start_codon(b"ATN", 2));
        assert!(is_start_codon(b"ATH", 2));
        assert!(!is_start_codon(b"ATN", 1));
        assert!(is_start_codon(b"atg", 1));
        assert!(is_start_codon(b"AUG", 1));
        assert!(!is_start_codon(b"TAG", 1));
        for gap in [&b"---"[..], b"T-A", b"TA.", b"..G"] {
            assert!(
                !is_start_codon(gap, 1),
                "{:?}",
                String::from_utf8_lossy(gap)
            );
            assert!(!is_ter_codon(gap, 1), "{:?}", String::from_utf8_lossy(gap));
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

    /// Ensembl's `reverse_comp` table (`Bio::EnsEMBL::Utils::Sequence`):
    /// `tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/`,
    /// every other character kept.
    #[test]
    fn complement_follows_ensembl_reverse_comp_for_iupac_letters() {
        assert_eq!(reverse_complement(b"RYMKBVDH"), b"DHBVMKRY");
        assert_eq!(reverse_complement(b"rymkbvdh"), b"dhbvmkry");
        assert_eq!(reverse_complement(b"SWNX-*<>"), b"><*-XNWS");
    }

    /// Ensembl's `$ALL_NUCLEOTIDES` (VariationFeatureOverlapAllele.pm:76,
    /// `qr/^[ACGTUMRWSYKVHDBXN-]+$/i`): one or more of those letters in either
    /// case or `-`; the other letters, the empty allele and every symbolic
    /// allele are refused.
    #[test]
    fn seq_is_dna_matches_perls_all_nucleotides_class() {
        assert!(seq_is_dna(b"ACGTUMRWSYKVHDBXN-"));
        assert!(seq_is_dna(b"acgtumrwsykvhdbxn-"));
        for letter in b"EFIJLOPQZ" {
            assert!(!seq_is_dna(&[*letter]), "{}", *letter as char);
            assert!(
                !seq_is_dna(&[letter.to_ascii_lowercase()]),
                "{}",
                letter.to_ascii_lowercase() as char
            );
        }
        for refused in [
            &b""[..],
            b"*",
            b"AC*",
            b"<DEL>",
            b"N[21:200[",
            b"A.",
            b"(3 BP INSERTION)",
            b"A/G",
            b" A",
        ] {
            assert!(
                !seq_is_dna(refused),
                "{:?}",
                String::from_utf8_lossy(refused)
            );
        }
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

    /// `Bio::Tools::CodonTable::translate("ATGGCNTGAAA")` is `MA*`: whole codons
    /// only, `GCN` resolved to the one amino acid its readings share, the stop
    /// kept; `RAY` codes aspartate or asparagine (`B`), `SAR` glutamate or
    /// glutamine (`Z`), `NNN` anything (`X`), and `---` is a gap.
    #[test]
    fn translate_sequence_follows_the_module_for_whole_and_ambiguous_codons() {
        assert_eq!(translate_sequence(b"ATGGCNTGAAA", 1), b"MA*");
        assert_eq!(translate_sequence(b"atgtga", 1), b"M*");
        assert_eq!(translate_sequence(b"AUGUGA", 1), b"M*");
        assert_eq!(translate_sequence(b"RAYSARNNN---A-C", 1), b"BZX-X");
        assert_eq!(translate_sequence(b"TGAAGAATA", 2), b"W*M");
        assert_eq!(translate_sequence(b"AT", 1), b"");
        assert_eq!(translate_sequence(b"ATGTGA", 7), b"XX");
    }
}
