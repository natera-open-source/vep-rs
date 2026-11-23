// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The source text the `implementation` tables of the fixture records name. Nothing
//! compiles it: `--check` reads it for the items the tables' rows name.

/// The row the still-present fixture records for vep-rs.
pub fn downstream_row() -> &'static str {
    "21:100 | A | ENST1 | downstream_gene_variant | 5000"
}

/// The value the closed fixture records for vep-rs.
pub const PLAIN_EQUALS_SIGN: &str = "ENSP2:p.Leu4=";

/// Whether a pair of values is the closed fixture's difference.
pub fn is_escaped_equals_sign(ensembl_vep: &str, vep_rs: &str) -> bool {
    ensembl_vep.replace("%3D", "=") == vep_rs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_downstream_row_is_written() {
        assert!(downstream_row().ends_with("5000"));
    }

    #[test]
    fn the_equals_sign_is_written_plain() {
        assert!(PLAIN_EQUALS_SIGN.ends_with('='));
    }
}
