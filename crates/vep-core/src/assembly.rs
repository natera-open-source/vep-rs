// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reference genome assembly identity.
//!
//! Several plugins read a single file that carries coordinates for both GRCh37
//! and GRCh38 in different columns (REVEL, and dbNSFP as `dbNSFP.pm:316-334`
//! reads it), or read structurally different bundles per assembly (LoFTEE's GERP
//! data is a tabix
//! TSV on GRCh37 and a bigWig on GRCh38). Such plugins cannot pick the right
//! column or the right reader without knowing which assembly the run targets,
//! so the resolved assembly is threaded from `--assembly` down into plugin
//! initialization.
//!
//! Ensembl VEP does the same: `REVEL.pm` (Ensembl/VEP_plugins release/115) selects
//! `hg19_pos` vs `grch38_pos` from `$self->{config}->{assembly}` and dies when the
//! required column is absent.

use std::fmt;

/// A reference genome assembly.
///
/// Parsed from the `--assembly` CLI value by its aliases, each whole and in any
/// case: `GRCh37`/`hg19`/`37` and `GRCh38`/`hg38`/`38`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Assembly {
    /// GRCh37, also known as hg19.
    Grch37,
    /// GRCh38, also known as hg38.
    Grch38,
}

impl Assembly {
    /// Parse an assembly from a user-supplied string.
    ///
    /// Returns `None` for anything unrecognized, so callers can distinguish
    /// "no assembly given" from "an unrecognized assembly" and decide
    /// whether that is fatal for their purposes. Only a whole alias matches, so
    /// another species' assembly whose name carries a human assembly's digits
    /// (`GRCm38`) is unrecognized rather than read as the human one.
    pub fn parse(value: &str) -> Option<Self> {
        match value.trim().to_ascii_lowercase().as_str() {
            "grch37" | "hg19" | "37" => Some(Assembly::Grch37),
            "grch38" | "hg38" | "38" => Some(Assembly::Grch38),
            _ => None,
        }
    }

    /// The canonical Ensembl name (`GRCh37` or `GRCh38`).
    pub fn as_str(&self) -> &'static str {
        match self {
            Assembly::Grch37 => "GRCh37",
            Assembly::Grch38 => "GRCh38",
        }
    }

    /// The UCSC-style name (`hg19` or `hg38`).
    ///
    /// Some data files are named and self-describe with this spelling
    /// (AlphaMissense's `genome` column, the SpliceAI file names of the plugin's
    /// usage example, SpliceAI.pm:117-129), so plugins that
    /// validate a file against the requested assembly need both forms.
    pub fn ucsc_name(&self) -> &'static str {
        match self {
            Assembly::Grch37 => "hg19",
            Assembly::Grch38 => "hg38",
        }
    }
}

impl fmt::Display for Assembly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_canonical_ensembl_names() {
        assert_eq!(Assembly::parse("GRCh37"), Some(Assembly::Grch37));
        assert_eq!(Assembly::parse("GRCh38"), Some(Assembly::Grch38));
    }

    #[test]
    fn parsing_is_case_insensitive_and_trims() {
        assert_eq!(Assembly::parse("  grch38  "), Some(Assembly::Grch38));
        assert_eq!(Assembly::parse("GRCH37"), Some(Assembly::Grch37));
    }

    #[test]
    fn parses_ucsc_and_bare_numeric_aliases() {
        assert_eq!(Assembly::parse("hg19"), Some(Assembly::Grch37));
        assert_eq!(Assembly::parse("hg38"), Some(Assembly::Grch38));
        assert_eq!(Assembly::parse("37"), Some(Assembly::Grch37));
        assert_eq!(Assembly::parse("38"), Some(Assembly::Grch38));
    }

    /// Another species' assembly is unrecognized even when its name carries a
    /// human assembly's digits, as the mouse `GRCm38` does, and so is a human
    /// name with a patch suffix: a plugin reading the hint gets none rather
    /// than the wrong human one.
    #[test]
    fn unrecognized_and_empty_values_are_none() {
        assert_eq!(Assembly::parse(""), None);
        assert_eq!(Assembly::parse("   "), None);
        assert_eq!(Assembly::parse("GRCm39"), None);
        assert_eq!(Assembly::parse("GRCm38"), None);
        assert_eq!(Assembly::parse("BDGP6.54"), None);
        assert_eq!(Assembly::parse("GRCh38.p14"), None);
        assert_eq!(Assembly::parse("T2T-CHM13"), None);
    }

    #[test]
    fn exposes_both_naming_conventions() {
        assert_eq!(Assembly::Grch37.as_str(), "GRCh37");
        assert_eq!(Assembly::Grch37.ucsc_name(), "hg19");
        assert_eq!(Assembly::Grch38.as_str(), "GRCh38");
        assert_eq!(Assembly::Grch38.ucsc_name(), "hg38");
    }

    #[test]
    fn display_uses_canonical_name() {
        assert_eq!(Assembly::Grch38.to_string(), "GRCh38");
    }
}
