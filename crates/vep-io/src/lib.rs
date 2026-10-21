// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Input parsers and output formatters for VEP.
//!
//! This crate provides:
//! - Input parsers for the Ensembl default, HGVS (genomic, transcript and
//!   protein notation), SPDI, region and variant identifier formats, and the
//!   detection of a format from the first data line; VCF is read by the
//!   command-line crate's line parser
//! - Default, tab, VCF, JSON and Parquet-intermediate output formatters

pub mod error;
pub mod input;
pub mod output;

pub use error::IoError;
pub use input::InputParser;
pub use output::OutputFormatter;

#[cfg(test)]
mod tests {
    use std::path::Path;

    /// The one VCF reader is the command-line crate's line parser: this crate
    /// declares no VCF input module, so no second reader can drift from it.
    #[test]
    fn no_second_vcf_parser() {
        let input = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/input");
        assert!(!input.join("vcf.rs").exists());
        assert!(!input.join("vcf").exists());
        let declared = std::fs::read_to_string(input.join("mod.rs")).unwrap();
        assert!(!declared.contains("mod vcf"));
    }
}
