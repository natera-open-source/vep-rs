// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The names a consequence type is counted under in `[Consequences (all)]`.
//!
//! `data/consequence_terms.tsv` holds one row for each entry of `%OVERLAP_CONSEQUENCES`
//! (ensembl-variation `Utils/Constants.pm:653-1388` at release/116, the entry
//! `intergenic_variant` being `$DEFAULT_OVERLAP_CONSEQUENCE` of `:636-650`): its `SO_term`,
//! its `display_term` and its `NCBI_term`, an absent key written as `-`. A run counts each
//! consequence under the name `--terms` chooses, the SO term where the vocabulary has none
//! (`OutputFactory.pm:1340-1341`), so a key of the chart is a consequence type exactly
//! where it stands in some cell of the table.

const TABLE: &str = include_str!("../data/consequence_terms.tsv");

/// The cell of a term the vocabulary has no name for.
const ABSENT: &str = "-";

/// Whether `key` names a consequence type in the `SO`, the `display` or the `NCBI`
/// vocabulary.
pub fn is_consequence_term(key: &str) -> bool {
    key != ABSENT
        && TABLE
            .lines()
            .skip(1)
            .flat_map(|row| row.split('\t'))
            .any(|cell| cell == key)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The table is three columns wide under its header, 41 rows long, and names every
    /// consequence once in its first column.
    #[test]
    fn the_table_holds_one_row_a_consequence() {
        let mut rows = TABLE.lines();
        assert_eq!(rows.next(), Some("SO_term\tdisplay_term\tNCBI_term"));
        let rows: Vec<Vec<&str>> = rows.map(|row| row.split('\t').collect()).collect();
        assert_eq!(rows.len(), 41);
        assert!(rows.iter().all(|row| row.len() == 3), "{rows:?}");
        let mut so_terms: Vec<&str> = rows.iter().map(|row| row[0]).collect();
        so_terms.sort_unstable();
        so_terms.dedup();
        assert_eq!(so_terms.len(), 41);
        assert!(!so_terms.contains(&ABSENT));
    }

    /// A term of each vocabulary is one, and a name of another kind is none: the dash of
    /// an absent cell, a column's title, a term in another letter case, the identifier of
    /// an annotation record.
    #[test]
    fn a_key_is_a_term_where_some_vocabulary_names_it() {
        for term in [
            "intergenic_variant",
            "5_prime_UTR_variant",
            "TFBS_ablation",
            "INTERGENIC",
            "WITHIN_MATURE_miRNA",
            "splice-3",
            "cds-synon",
        ] {
            assert!(is_consequence_term(term), "{term}");
        }
        for other in [
            "-",
            "",
            "SO_term",
            "display_term",
            "Intergenic_variant",
            "intergenic",
            "exact_2",
            "rs1444031311",
            "intergenic_variant\tINTERGENIC",
        ] {
            assert!(!is_consequence_term(other), "{other:?}");
        }
    }
}
