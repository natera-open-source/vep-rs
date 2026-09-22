// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Region format parser.
//!
//! Parses lines of the form:
//! ```text
//! chr:start-end:strand/allele
//! chr:start-end/allele
//! ```
//!
//! Examples:
//! ```text
//! 21:25585733-25585733:1/G
//! 21:25585656-25585660:-1/ACGTG
//! X:100000-100000/T
//! ```

use std::io::BufRead;

use vep_core::coordinate::Strand;
use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::InputParser;

/// Parser for region format input.
pub struct RegionParser<R: BufRead> {
    reader: R,
    line_buf: String,
}

impl<R: BufRead> RegionParser<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: String::new(),
        }
    }

    /// Parse a single line into an `InputVariant`.
    fn parse_line(&self, line: &str) -> Result<Option<InputVariant>, IoError> {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return Ok(None);
        }

        let slash_pos = line
            .rfind('/')
            .ok_or_else(|| IoError::VcfParse(format!("No '/' found in region: '{}'", line)))?;

        let coord_part = &line[..slash_pos];
        let allele = &line[slash_pos + 1..];

        if allele.is_empty() {
            return Err(IoError::VcfParse(format!(
                "Empty allele in region: '{}'",
                line
            )));
        }

        let parts: Vec<&str> = coord_part.split(':').collect();

        if parts.len() < 2 {
            return Err(IoError::VcfParse(format!(
                "Invalid region format: '{}' (expected chr:start-end)",
                line
            )));
        }

        let chr = parts[0].to_string();

        let range_part = parts[1];
        let (start, end) = if let Some((s, e)) = range_part.split_once('-') {
            let start: u64 = s
                .parse()
                .map_err(|_| IoError::VcfParse(format!("Invalid start position: '{}'", s)))?;
            let end: u64 = e
                .parse()
                .map_err(|_| IoError::VcfParse(format!("Invalid end position: '{}'", e)))?;
            (start, end)
        } else {
            let pos: u64 = range_part
                .parse()
                .map_err(|_| IoError::VcfParse(format!("Invalid position: '{}'", range_part)))?;
            (pos, pos)
        };

        let strand = if parts.len() >= 3 {
            let strand_val: i8 = parts[2]
                .parse()
                .map_err(|_| IoError::VcfParse(format!("Invalid strand: '{}'", parts[2])))?;
            Strand::from_i8(strand_val).unwrap_or(Strand::Forward)
        } else {
            Strand::Forward
        };

        // The region format carries no reference allele, so ref is `-`.
        let mut variant =
            InputVariant::new(chr, start, end, b"-".to_vec(), allele.as_bytes().to_vec());
        variant.strand = strand;

        Ok(Some(variant))
    }
}

impl<R: BufRead + Send> InputParser for RegionParser<R> {
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>> {
        loop {
            self.line_buf.clear();
            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => return None,
                Ok(_) => match self.parse_line(&self.line_buf) {
                    Ok(Some(variant)) => return Some(Ok(variant)),
                    Ok(None) => continue,
                    Err(e) => return Some(Err(e)),
                },
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
    fn test_parse_with_strand() {
        let data = "21:25585733-25585733:1/G\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = RegionParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585733);
        assert_eq!(v.end, 25585733);
        assert_eq!(v.alt_allele(), b"G");
        assert_eq!(v.strand, Strand::Forward);
    }

    #[test]
    fn test_parse_without_strand() {
        let data = "X:100000-100000/T\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = RegionParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "X");
        assert_eq!(v.start, 100000);
        assert_eq!(v.end, 100000);
        assert_eq!(v.alt_allele(), b"T");
        assert_eq!(v.strand, Strand::Forward);
    }

    #[test]
    fn test_parse_reverse_strand() {
        let data = "21:25585656-25585660:-1/ACGTG\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = RegionParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585656);
        assert_eq!(v.end, 25585660);
        assert_eq!(v.alt_allele(), b"ACGTG");
        assert_eq!(v.strand, Strand::Reverse);
    }

    #[test]
    fn test_ref_is_dash() {
        let data = "21:25585733-25585733:1/G\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = RegionParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.ref_allele, b"-");
    }

    #[test]
    fn test_skip_comments() {
        let data = "# comment\n21:100-200:1/A\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = RegionParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
    }

    #[test]
    fn test_no_slash_error() {
        let data = "21:100-200\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = RegionParser::new(reader);

        let result = parser.next_variant().unwrap();
        assert!(result.is_err());
    }

    #[test]
    fn test_multiple_lines() {
        let data = "21:100-100:1/A\n22:200-200:1/G\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = RegionParser::new(reader);

        let v1 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v1.chr, "21");

        let v2 = parser.next_variant().unwrap().unwrap();
        assert_eq!(v2.chr, "22");

        assert!(parser.next_variant().is_none());
    }
}
