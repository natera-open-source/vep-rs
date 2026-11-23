// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Blosum62 plugin: the BLOSUM62 substitution score of a missense change, as
//! Ensembl VEP's `Blosum62.pm` writes it.
//!
//! **Usage:** `--plugin Blosum62` (the plugin takes no parameter; a positional one
//! is not read, and is logged; a `key=value` one is refused)
//!
//! Output field: `BLOSUM62`, the matrix score for the reference and alternate
//! amino acids of a row whose amino-acid change is one capital letter over
//! another (`Blosum62.pm:110`), looked up in the plugin's embedded 20 by 20
//! matrix (`:50-73`); a pair outside the matrix (a stop, an unknown residue)
//! and every other row carry nothing (`:114-121`). A structural variant's rows
//! carry nothing: Ensembl VEP hands a plugin only the variant classes its
//! `variant_feature_types` names (`OutputFactory.pm:2229`), and `Blosum62.pm`
//! keeps the default, sequence variants alone (`BaseVepPlugin.pm:100`).
//!
//! The amino-acid change is read from the row's `Amino_acids` column, which
//! prints `pep_allele_string` as the plugin reads it. On a row of the
//! documented class `alternate_codon_past_the_coding_end` the column's
//! alternate half is read from vep-rs's own alternate codon, and neither
//! engine writes a score there: in both, the reference codons or the alternate
//! ones end in a partial codon (`TranscriptVariationAllele.pm:774-776`).
//! Under any other documented class whose rows differ between the engines in
//! that column, each engine scores the residues it reads itself, so a score
//! can differ or stand in one engine alone: under
//! `hgvs_shift_shared_across_alleles` Ensembl VEP alone can write one, on a
//! later allele's row, which the class
//! `blosum62_hgvs_shift_shared_across_alleles` reads.
//!
//! Per-transcript plugin: returns annotation from [`crate::BuiltinPlugin::run`]
//! via the default [`crate::BuiltinPlugin::run_batch`] dispatcher.

use indexmap::IndexMap;

use crate::params::parse_none;
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// Output field name.
const OUTPUT_FIELD: &str = "BLOSUM62";

/// The residues in the matrix's row and column order (`Blosum62.pm:73`).
const RESIDUES: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

/// The BLOSUM62 matrix, row by row in [`RESIDUES`] order (`Blosum62.pm:50-71`).
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

/// BLOSUM62 substitution score plugin.
#[derive(Debug, Default)]
pub struct Blosum62Plugin;

impl Blosum62Plugin {
    /// Create the plugin; `--plugin Blosum62` carries no parameter to initialize.
    pub fn new() -> Self {
        Self
    }
}

impl crate::traits::sealed::Sealed for Blosum62Plugin {}

impl BuiltinPlugin for Blosum62Plugin {
    fn name(&self) -> &str {
        "Blosum62"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            OUTPUT_FIELD.to_string(),
            "BLOSUM62 substitution score for the reference and alternative amino acids".to_string(),
        )]
    }

    /// A positional parameter is not read, as the Ensembl VEP module reads
    /// none (`Blosum62.pm:75-88`); it is logged. A `key=value` list is refused
    /// as every built-in refuses an undeclared key.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        parse_none(self.name(), params)
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();
        if variant.is_structural {
            return Ok(result);
        }
        if let Some(score) = consequence.amino_acids.as_deref().and_then(score) {
            result.insert(OUTPUT_FIELD.to_string(), score.to_string());
        }
        Ok(result)
    }
}

/// The score of an amino-acid change written `X/Y` with one capital letter on
/// each side (`Blosum62.pm:110-118`); `None` for any other column text and for
/// a letter outside the matrix.
fn score(amino_acids: &str) -> Option<i8> {
    let [reference, b'/', alternate] = amino_acids.as_bytes() else {
        return None;
    };
    if !(reference.is_ascii_uppercase() && alternate.is_ascii_uppercase()) {
        return None;
    }
    let index = |residue: &u8| RESIDUES.iter().position(|r| r == residue);
    Some(BLOSUM62[index(reference)?][index(alternate)?])
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::consequence::FeatureType;

    fn row(amino_acids: Option<&str>) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.amino_acids = amino_acids.map(str::to_string);
        consequence.feature_type = FeatureType::Transcript;
        consequence
    }

    fn snv() -> InputVariant {
        InputVariant::new("21".into(), 100, 100, b"C".to_vec(), b"T".to_vec())
    }

    fn blosum62(consequence: &TranscriptConsequence, variant: &InputVariant) -> Option<String> {
        Blosum62Plugin::new()
            .run(consequence, variant)
            .unwrap()
            .get(OUTPUT_FIELD)
            .cloned()
    }

    #[test]
    fn the_matrix_is_symmetric_with_the_published_diagonal() {
        for (i, row) in BLOSUM62.iter().enumerate() {
            for (j, &score) in row.iter().enumerate() {
                assert_eq!(score, BLOSUM62[j][i], "{i} {j}");
            }
        }
        let diagonal: Vec<i8> = BLOSUM62.iter().enumerate().map(|(i, row)| row[i]).collect();
        assert_eq!(
            diagonal,
            [4, 5, 6, 6, 9, 5, 5, 6, 8, 4, 4, 5, 5, 6, 7, 4, 5, 11, 7, 4]
        );
    }

    /// Ensembl VEP writes `BLOSUM62=0` for `A/V`, `-1` for `A/R` and `11` for
    /// `W/W`: the matrix cell of the reference row and alternate column.
    #[test]
    fn a_missense_change_scores_its_matrix_cell() {
        assert_eq!(blosum62(&row(Some("A/V")), &snv()).as_deref(), Some("0"));
        assert_eq!(blosum62(&row(Some("A/R")), &snv()).as_deref(), Some("-1"));
        assert_eq!(blosum62(&row(Some("W/W")), &snv()).as_deref(), Some("11"));
        assert_eq!(blosum62(&row(Some("D/K")), &snv()).as_deref(), Some("-1"));
        assert_eq!(blosum62(&row(Some("C/W")), &snv()).as_deref(), Some("-2"));
    }

    /// A synonymous row prints one letter, a stop-gained row `Q/*`, a
    /// frameshift row `X`, an unknown residue a letter outside the matrix:
    /// none matches `^([A-Z])\/([A-Z])$` with both letters in the matrix.
    #[test]
    fn other_amino_acid_columns_carry_nothing() {
        for column in [
            None,
            Some("A"),
            Some("Q/*"),
            Some("*/Q"),
            Some("X"),
            Some("A/X"),
        ] {
            assert_eq!(blosum62(&row(column), &snv()), None, "{column:?}");
        }
        assert_eq!(blosum62(&row(Some("A/U")), &snv()), None);
        assert_eq!(blosum62(&row(Some("a/v")), &snv()), None);
        assert_eq!(blosum62(&row(Some("AA/V")), &snv()), None);
    }

    #[test]
    fn a_structural_variant_row_carries_nothing() {
        let mut deletion = InputVariant::new("21".into(), 100, 500, b"N".to_vec(), b"-".to_vec());
        deletion.is_structural = true;
        assert_eq!(blosum62(&row(Some("A/V")), &deletion), None);
    }

    /// The Ensembl VEP module reads no parameter (`Blosum62.pm:75-88`), so a
    /// positional one is not read and the plugin loads; a `key=value` list is
    /// refused as the undeclared key every built-in refuses.
    #[test]
    fn a_parameter_is_not_read() {
        let mut plugin = Blosum62Plugin::new();
        assert!(plugin.init(&[]).is_ok());
        assert!(plugin.init(&["x".to_string()]).is_ok());
        assert!(plugin.init(&["x".to_string(), "y".to_string()]).is_ok());
        assert_eq!(
            plugin
                .init(&["file=x".to_string()])
                .unwrap_err()
                .to_string(),
            "plugin initialization failed: Blosum62: unknown parameter 'file'; Blosum62 takes no \
             parameter"
        );
        assert_eq!(
            plugin
                .run(&row(Some("A/V")), &snv())
                .unwrap()
                .get(OUTPUT_FIELD)
                .map(String::as_str),
            Some("0")
        );
    }
}
