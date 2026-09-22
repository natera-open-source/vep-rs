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

use vep_core::coordinate::Strand;
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
        let start: u64 = parts[1]
            .parse()
            .map_err(|_| IoError::VcfParse(format!("Invalid start position: '{}'", parts[1])))?;
        let end: u64 = parts[2]
            .parse()
            .map_err(|_| IoError::VcfParse(format!("Invalid end position: '{}'", parts[2])))?;
        let allele_string = parts[3];
        let strand_val: i8 = parts[4]
            .parse()
            .map_err(|_| IoError::VcfParse(format!("Invalid strand: '{}'", parts[4])))?;

        let strand = Strand::from_i8(strand_val).unwrap_or(Strand::Forward);

        let id = if parts.len() >= 6 {
            let id_str = parts[5];
            if id_str == "." || id_str.is_empty() {
                None
            } else {
                Some(id_str.to_string())
            }
        } else {
            None
        };

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
        assert_eq!(v.id, None);

        assert!(parser.next_variant().is_none());
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
