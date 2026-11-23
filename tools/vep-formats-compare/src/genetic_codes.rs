// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The genetic codes the comparison translates under: the tables of
//! `Bio::Tools::CodonTable`, read from the repository's table of the 64 codons under every
//! code of that module, and genetic code 26, which vep-rs reads and the module lacks; each
//! engine's reading of a code, the codes that read every codon as the standard code does,
//! and the refusal of a code no table is held for.

use std::collections::BTreeMap;
use std::sync::LazyLock;

/// The 64 codons under every genetic code of `Bio::Tools::CodonTable`, under `#` header
/// lines one line a codon: the table id, the codon, its amino acid, its start mark and its
/// stop mark, tab-separated. The rows of id 0 repeat table 1's.
const MODULE_CODES_TABLE: &str =
    include_str!("../../../crates/vep-core/tests/data/genetic_codes.tsv");

/// Genetic code 26 as NCBI's genetic code table gives it, the amino acid of each codon in the
/// order `TTT TTC TTA TTG TCT ... GGG`: vep-rs reads it (`crates/vep-core/src/codon.rs`,
/// `PACHYSOLEN_NUCLEAR_CODE`) and `Bio::Tools::CodonTable` holds no table for it.
pub const CODES_THE_MODULE_LACKS: &[(i64, &[u8; 64])] = &[(
    26,
    b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
)];

/// The id of the standard code.
pub const STANDARD_CODE: i64 = 1;

/// One genetic code: the amino acid of each codon, a stop `*`, the codons in the order
/// `TTT TTC TTA TTG TCT ... GGG`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeneticCode([u8; 64]);

/// The place of a codon of the letters `T`, `C`, `A` and `G` in that order; `None` for any
/// other text.
fn codon_place(codon: &[u8]) -> Option<usize> {
    let &[first, second, third] = codon else {
        return None;
    };
    let place = |base: u8| match base {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    Some(16 * place(first)? + 4 * place(second)? + place(third)?)
}

impl GeneticCode {
    /// The amino acid of a codon of the letters `T`, `C`, `A` and `G`, a stop `*`; `None`
    /// for any other text.
    pub fn amino_acid(&self, codon: &[u8]) -> Option<u8> {
        Some(self.0[codon_place(codon)?])
    }

    /// The amino acid of each codon, in the order of the codons.
    pub fn amino_acids(&self) -> &[u8; 64] {
        &self.0
    }
}

/// The amino acids of the 64 codons under each table id of a table in the layout of
/// `MODULE_CODES_TABLE`; a line that is not one, or a table short of a codon, is a panic,
/// since the text is part of the build.
fn codes_by_id(table: &str) -> BTreeMap<i64, [u8; 64]> {
    let mut codes: BTreeMap<i64, [u8; 64]> = BTreeMap::new();
    for line in table.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut cells = line.split('\t');
        let id = cells.next().and_then(|id| id.parse::<i64>().ok());
        let place = cells.next().and_then(|codon| codon_place(codon.as_bytes()));
        let residue = cells.next().map(str::as_bytes);
        let (Some(id), Some(place), Some(&[residue])) = (id, place, residue) else {
            panic!("the genetic-code table holds a line that names no codon of a table: {line}")
        };
        codes.entry(id).or_insert([0; 64])[place] = residue;
    }
    for (id, amino_acids) in &codes {
        assert!(
            amino_acids.iter().all(|residue| *residue != 0),
            "the genetic-code table holds fewer than 64 codons under id {id}"
        );
    }
    codes
}

/// The genetic codes of `Bio::Tools::CodonTable` by table id. Id 0, which the module answers
/// as table 1, is left out: a cache's `codon_table` of 0 is read as table 1 before any
/// lookup here.
pub static MODULE_CODES: LazyLock<BTreeMap<i64, GeneticCode>> = LazyLock::new(|| {
    codes_by_id(MODULE_CODES_TABLE)
        .into_iter()
        .filter(|(id, _)| *id != 0)
        .map(|(id, amino_acids)| (id, GeneticCode(amino_acids)))
        .collect()
});

/// The engine whose reading of a genetic code a value is derived under.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reading {
    /// Ensembl VEP translates with `Bio::Tools::CodonTable`, which reads an id it holds no
    /// table for under its last table (CodonTable.pm:328-338, :412-460).
    EnsemblVep,
    /// vep-rs reads every code it holds as that code.
    VepRs,
}

/// The genetic code of `table` as `reading` reads it; `None` for an id no table is held for.
pub fn genetic_code(table: i64, reading: Reading) -> Option<GeneticCode> {
    if let Some(code) = MODULE_CODES.get(&table) {
        return Some(*code);
    }
    let (_, amino_acids) = CODES_THE_MODULE_LACKS.iter().find(|(id, _)| *id == table)?;
    match reading {
        Reading::VepRs => Some(GeneticCode(**amino_acids)),
        Reading::EnsemblVep => MODULE_CODES.last_key_value().map(|(_, code)| *code),
    }
}

/// The ids a table is held for, ascending: the module's and the codes it lacks.
pub fn held_tables() -> Vec<i64> {
    let mut ids: Vec<i64> = MODULE_CODES
        .keys()
        .copied()
        .chain(CODES_THE_MODULE_LACKS.iter().map(|(id, _)| *id))
        .collect();
    ids.sort_unstable();
    ids
}

/// Whether Ensembl VEP reads every codon under `table` as it reads it under the standard
/// code; `None` for an id no table is held for.
pub fn reads_as_the_standard_code(table: i64) -> Option<bool> {
    let standard = genetic_code(STANDARD_CODE, Reading::EnsemblVep)?;
    Some(genetic_code(table, Reading::EnsemblVep)? == standard)
}

/// The genetic code a cache transcript names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CodonTable {
    /// A code a table is held for, by id.
    Held(i64),
    /// A `codon_table` that names no code a table is held for, as the cache spells it.
    Lacking(String),
}

impl CodonTable {
    /// The id of the code, or the refusal of a code no table is held for, which names
    /// `transcript`: a row of such a transcript is never read under another code.
    pub fn held(&self, transcript: &str) -> Result<i64, String> {
        match self {
            CodonTable::Held(table) => Ok(*table),
            CodonTable::Lacking(value) => Err(lacking_table(transcript, value)),
        }
    }
}

/// The refusal of the transcript `transcript`, whose cache record's `codon_table` is `value`.
pub fn lacking_table(transcript: &str, value: &str) -> String {
    let held: Vec<String> = held_tables().iter().map(i64::to_string).collect();
    format!(
        "transcript {transcript}: the cache's codon_table {value} names no genetic code the \
         comparison holds a table for (it holds {})",
        held.join(", ")
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_is_read_whole_and_its_id_0_repeats_table_1() {
        let read = codes_by_id(MODULE_CODES_TABLE);
        assert_eq!(read.get(&0), read.get(&STANDARD_CODE));
        assert_eq!(read.len(), MODULE_CODES.len() + 1);
        for (id, code) in MODULE_CODES.iter() {
            assert_eq!(Some(code.amino_acids()), read.get(id), "table {id}");
        }
        assert!(!MODULE_CODES.contains_key(&0));
    }

    #[test]
    fn a_code_the_module_lacks_is_read_by_each_engine_in_its_own_way() {
        let last = *MODULE_CODES.last_key_value().expect("a table").1;
        for (id, amino_acids) in CODES_THE_MODULE_LACKS {
            assert!(!MODULE_CODES.contains_key(id), "code {id}");
            assert_eq!(genetic_code(*id, Reading::EnsemblVep), Some(last));
            assert_eq!(
                genetic_code(*id, Reading::VepRs)
                    .as_ref()
                    .map(GeneticCode::amino_acids),
                Some(*amino_acids)
            );
        }
        for (id, code) in MODULE_CODES.iter() {
            for reading in [Reading::EnsemblVep, Reading::VepRs] {
                assert_eq!(genetic_code(*id, reading), Some(*code), "table {id}");
            }
        }
        for id in [-1, 0, 7, 8, 17, 20, 25, 27, i64::MAX] {
            for reading in [Reading::EnsemblVep, Reading::VepRs] {
                assert_eq!(genetic_code(id, reading), None, "id {id}");
            }
            assert_eq!(reads_as_the_standard_code(id), None, "id {id}");
        }
    }

    /// The codes that read as the standard one are found by comparing the tables: each one
    /// listed here holds the standard code's amino acid under every codon in the table the
    /// module reads it under (`genetic_codes.tsv`, the 64 rows of each id), and no other
    /// held code does.
    #[test]
    fn the_codes_that_read_as_the_standard_code_are_derived_from_the_table() {
        let standard: Vec<i64> = held_tables()
            .into_iter()
            .filter(|id| reads_as_the_standard_code(*id) == Some(true))
            .collect();
        assert_eq!(standard, [1, 11, 24, 26]);
        let read = codes_by_id(MODULE_CODES_TABLE);
        for id in held_tables() {
            let last = *read.keys().next_back().expect("a table");
            let module_id = if read.contains_key(&id) { id } else { last };
            assert_eq!(
                reads_as_the_standard_code(id),
                Some(read[&module_id] == read[&STANDARD_CODE]),
                "code {id}"
            );
        }
    }

    #[test]
    fn a_code_no_table_is_held_for_is_refused_by_the_transcript_and_the_value() {
        assert_eq!(CodonTable::Held(4).held("T1"), Ok(4));
        let refused = CodonTable::Lacking("7".to_string())
            .held("ENST00000000106")
            .unwrap_err();
        assert_eq!(
            refused,
            "transcript ENST00000000106: the cache's codon_table 7 names no genetic code the \
             comparison holds a table for (it holds 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, \
             15, 16, 21, 22, 23, 24, 26)"
        );
    }
}
