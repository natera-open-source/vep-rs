// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Ensembl VEP default input format parser.
//!
//! Parses lines of the form:
//! ```text
//! chr  start  end  allele_string  strand  [identifier]
//! 21   25585733  25585733  A/G  1  rs699
//! 21   25585656  25585660  ACGTG/-  1  .
//! ```
//!
//! Tab or whitespace delimited, 5-6 columns. Multi-allelic allele strings
//! (e.g., "A/G/T") produce one `InputVariant` per ALT allele.

use std::io::BufRead;
use std::sync::Arc;

use vep_core::coordinate::{parse_coordinate, Strand};
use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::InputParser;

/// Parser for Ensembl VEP default input format.
pub struct EnsemblParser<R: BufRead> {
    reader: R,
    line_buf: String,
    /// Buffer for multi-allelic variants from a single line.
    buffer: Vec<InputVariant>,
}

impl<R: BufRead> EnsemblParser<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: String::new(),
            buffer: Vec::new(),
        }
    }

    /// Parse a single line into one or more `InputVariant`s.
    fn parse_line(&self, line: &str) -> Result<Vec<InputVariant>, IoError> {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return Ok(Vec::new());
        }

        let parts: Vec<&str> = if line.contains('\t') {
            line.split('\t').collect()
        } else {
            line.split_whitespace().collect()
        };

        if parts.len() < 5 {
            return Err(IoError::VcfParse(format!(
                "Ensembl format requires at least 5 columns, got {}",
                parts.len()
            )));
        }

        let chr = parts[0].to_string();
        let start = parse_coordinate(parts[1])
            .ok_or_else(|| IoError::VcfParse(format!("Invalid start position: '{}'", parts[1])))?;
        let end = parse_coordinate(parts[2])
            .ok_or_else(|| IoError::VcfParse(format!("Invalid end position: '{}'", parts[2])))?;
        let allele_string = parts[3];
        let strand_val: i8 = parts[4]
            .parse()
            .map_err(|_| IoError::VcfParse(format!("Invalid strand: '{}'", parts[4])))?;

        let strand = Strand::from_i8(strand_val).unwrap_or(Strand::Forward);

        // The ID column whole, `.` included, as Parser/VEP_input.pm:127 takes it; a
        // `;` is part of the name here, where a VCF ID field is split on it.
        let id = parts
            .get(5)
            .filter(|id| !id.is_empty())
            .map(|id| id.to_string());

        let alleles: Vec<&str> = allele_string.split('/').collect();
        if alleles.len() < 2 {
            return Err(IoError::VcfParse(format!(
                "Invalid allele string: '{}' (expected REF/ALT)",
                allele_string
            )));
        }

        let ref_allele = alleles[0];
        let alts = &alleles[1..];
        let multi_allelic = alts.len() > 1;

        let mut variants = Vec::with_capacity(alts.len());
        // OutputFactory.pm `VariationFeature_to_output_hash`, the ensembl `_line`
        // branch: the lower coordinate, so an insertion written `100 99` is named
        // at 99.
        let line_name: Arc<str> = format!("{chr}_{}_{allele_string}", start.min(end)).into();

        for (i, alt) in alts.iter().enumerate() {
            let mut variant = InputVariant::new(
                chr.clone(),
                start,
                end,
                ref_allele.as_bytes().to_vec(),
                alt.as_bytes().to_vec(),
            );
            variant.strand = strand;
            variant.id = id.clone();
            variant.allele_index = i;
            variant.minimised = multi_allelic;
            variant.input_line_name = Some(line_name.clone());

            variants.push(variant);
        }

        Ok(variants)
    }
}

impl<R: BufRead + Send> InputParser for EnsemblParser<R> {
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>> {
        if let Some(variant) = self.buffer.pop() {
            return Some(Ok(variant));
        }

        loop {
            self.line_buf.clear();
            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => return None,
                Ok(_) => {
                    match self.parse_line(&self.line_buf) {
                        Ok(variants) => {
                            if variants.is_empty() {
                                continue; // Skip empty/comment lines
                            }
                            let mut variants = variants;
                            variants.reverse();
                            let first = variants.pop().unwrap();
                            self.buffer = variants;
                            return Some(Ok(first));
                        }
                        Err(e) => return Some(Err(e)),
                    }
                }
                Err(e) => return Some(Err(IoError::Io(e))),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use vep_core::coordinate::MAX_COORD;

    #[test]
    fn test_parse_snv() {
        let data = "21\t25585733\t25585733\tA/G\t1\trs699\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = EnsemblParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585733);
        assert_eq!(v.end, 25585733);
        assert_eq!(v.ref_allele, b"A");
        assert_eq!(v.alt_allele(), b"G");
        assert_eq!(v.strand, Strand::Forward);
        assert_eq!(v.id, Some("rs699".to_string()));

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_parse_deletion() {
        let data = "21\t25585656\t25585660\tACGTG/-\t1\t.\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = EnsemblParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585656);
        assert_eq!(v.end, 25585660);
        assert_eq!(v.ref_allele, b"ACGTG");
        assert_eq!(v.alt_allele(), b"-");
        assert_eq!(v.id, Some(".".to_string()));

        assert!(parser.next_variant().is_none());
    }

    /// OutputFactory.pm `VariationFeature_to_output_hash`, the ensembl `_line`
    /// branch: an ID-less line is named `chr_start_alleles` with the lower of its
    /// two coordinates (99 for the insertion `100 99`) and the allele column as
    /// written (`A/G/T` on every allele of a multi-allelic line); an ID names it.
    #[test]
    fn id_less_line_is_named_by_its_lower_coordinate_and_allele_column() {
        let data = "21\t100\t99\t-/A\t1\n21\t200\t200\tA/G/T\t1\t.\n21\t300\t300\tC/T\t1\trs1\n";
        let mut parser = EnsemblParser::new(Cursor::new(data.as_bytes()));
        let mut names = Vec::new();
        while let Some(v) = parser.next_variant() {
            names.push(v.unwrap().uploaded_variation());
        }
        assert_eq!(names, ["21_99_-/A", "21_200_A/G/T", "21_200_A/G/T", "rs1"]);
    }

    /// Parser/VEP_input.pm:127 takes the ID column whole, so `a;b` names the record
    /// `a;b` in `Uploaded_variation` and in the JSON `id` (a VCF ID field is split on
    /// `;`). A `.` stays the JSON `id` (Perl reads it as true) while the record is
    /// named by its line; `0` is replaced by `validate_vf`'s generated name, here
    /// `21_100_-/A` (start 100, allele string `-/A`), which OutputFactory.pm rejects,
    /// so the record is named by its line at the lower coordinate, `21_99_-/A`.
    #[test]
    fn an_ensembl_format_id_column_is_kept_whole() {
        let data = "21\t300\t300\tC/T\t1\ta;b\n21\t200\t200\tA/G\t1\t.\n21\t100\t99\t-/A\t1\t0\n";
        let mut parser = EnsemblParser::new(Cursor::new(data.as_bytes()));
        let mut names = Vec::new();
        while let Some(v) = parser.next_variant() {
            let v = v.unwrap();
            names.push((v.uploaded_variation(), v.variation_name()));
        }
        assert_eq!(
            names,
            [
                ("a;b".to_string(), "a;b".to_string()),
                ("21_200_A/G".to_string(), ".".to_string()),
                ("21_99_-/A".to_string(), "21_100_-/A".to_string()),
            ]
        );
    }

    #[test]
    fn test_parse_multi_allelic() {
        let data = "21\t25585733\t25585733\tA/G/T\t1\trs699\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = EnsemblParser::new(reader);

        let v1 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v1.alt_allele(), b"G");
        assert_eq!(v1.allele_index, 0);
        assert!(v1.minimised);

        let v2 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v2.alt_allele(), b"T");
        assert_eq!(v2.allele_index, 1);
        assert!(v2.minimised);

        assert!(parser.next_variant().is_none());
    }

    /// A coordinate is read up to `MAX_COORD`; a start or end past it fails the line.
    #[test]
    fn coordinates_past_the_ceiling_fail_the_line() {
        let past = MAX_COORD + 1;
        let data = format!(
            "21\t{MAX_COORD}\t{MAX_COORD}\tA/G\t1\n\
             21\t{past}\t{past}\tA/G\t1\n\
             21\t100\t{}\tA/G\t1\n",
            i64::MAX
        );
        let mut parser = EnsemblParser::new(Cursor::new(data.as_bytes()));
        let kept = parser.next_variant().unwrap().unwrap();
        assert_eq!((kept.start, kept.end), (MAX_COORD, MAX_COORD));
        let start_err = parser.next_variant().unwrap().unwrap_err().to_string();
        assert!(start_err.contains("Invalid start position"), "{start_err}");
        let end_err = parser.next_variant().unwrap().unwrap_err().to_string();
        assert!(end_err.contains("Invalid end position"), "{end_err}");
    }

    #[test]
    fn test_missing_identifier() {
        let data = "21\t25585733\t25585733\tA/G\t1\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = EnsemblParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.id, None);
    }

    #[test]
    fn test_whitespace_delimited() {
        let data = "21  25585733  25585733  A/G  1  rs699\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = EnsemblParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585733);
        assert_eq!(v.id, Some("rs699".to_string()));
    }

    #[test]
    fn test_reverse_strand() {
        let data = "21\t25585733\t25585733\tA/G\t-1\trs699\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = EnsemblParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.strand, Strand::Reverse);
    }

    #[test]
    fn test_skip_comments_and_empty() {
        let data = "# comment\n\n21\t25585733\t25585733\tA/G\t1\trs699\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = EnsemblParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_multiple_lines() {
        let data = "21\t25585733\t25585733\tA/G\t1\trs699\n\
                     21\t25587758\t25587758\tG/A\t1\trs4149056\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = EnsemblParser::new(reader);

        let v1 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v1.start, 25585733);
        assert_eq!(v1.id, Some("rs699".to_string()));

        let v2 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v2.start, 25587758);
        assert_eq!(v2.id, Some("rs4149056".to_string()));

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_too_few_columns() {
        let data = "21\t25585733\t25585733\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = EnsemblParser::new(reader);

        let result = parser.next_variant().unwrap();
        assert!(result.is_err());
    }
}
