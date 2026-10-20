// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Output formatters for VEP results.

pub mod fields;
pub mod frequency;
pub mod json;
pub mod parquet;
pub mod tab;
pub mod vcf_output;

use std::io::Write;

use vep_core::variant::InputVariant;

use crate::error::IoError;

/// Trait for formatting annotated variants into output lines.
pub trait OutputFormatter: Send {
    /// Write the header section to the output.
    fn write_header(&mut self, writer: &mut dyn Write) -> Result<(), IoError>;

    /// Format a single annotated variant into output lines.
    /// Returns one line per transcript consequence (or one line for intergenic).
    fn format_variant(&self, variant: &InputVariant) -> Result<Vec<String>, IoError>;

    /// Format every allele split from one input record together, so a format
    /// that writes one unit per record (VCF, JSON) can do so and the others can
    /// order a record's rows the way VEP does. The default formats each allele
    /// on its own.
    fn format_record(&self, record: &[&InputVariant]) -> Result<Vec<String>, IoError> {
        let mut lines = Vec::new();
        for variant in record {
            lines.extend(self.format_variant(variant)?);
        }
        Ok(lines)
    }

    /// [`Self::format_record`] as bytes, for a format that writes input bytes
    /// back and so must carry a line that is not UTF-8 as it was read. The
    /// default is the text lines' bytes.
    fn format_record_bytes(&self, record: &[&InputVariant]) -> Result<Vec<Vec<u8>>, IoError> {
        Ok(self
            .format_record(record)?
            .into_iter()
            .map(String::into_bytes)
            .collect())
    }

    /// Write any footer/closing content to the output.
    fn finish(&mut self, writer: &mut dyn Write) -> Result<(), IoError>;
}
