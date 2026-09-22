// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Simplified genomic HGVS notation parser.
//!
//! Supports basic genomic HGVS expressions:
//! ```text
//! 21:g.25585733A>G          (substitution)
//! 21:g.25585656_25585660del (deletion)
//! 21:g.25585733delAinsGT    (delins)
//! 21:g.25585734_25585735insACG (insertion)
//! ```
//!
//! Transcript-level HGVS (e.g., `ENST00000366667:c.803C>T`) is unsupported
//! and returns an error.

use std::io::BufRead;

use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::InputParser;

/// Parser for genomic HGVS notation.
pub struct HgvsParser<R: BufRead> {
    reader: R,
    line_buf: String,
}

impl<R: BufRead> HgvsParser<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: String::new(),
        }
    }

    /// Parse a single HGVS expression into an `InputVariant`.
    fn parse_line(&self, line: &str) -> Result<Option<InputVariant>, IoError> {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return Ok(None);
        }

        if line.contains(":c.")
            || line.contains(":p.")
            || line.contains(":n.")
            || line.contains(":r.")
        {
            return Err(IoError::NotImplemented(
                "transcript HGVS is not supported".into(),
            ));
        }

        let (chr, rest) = line
            .split_once(":g.")
            .ok_or_else(|| IoError::VcfParse(format!("Invalid HGVS notation: '{}'", line)))?;

        let chr = chr.to_string();

        if rest.contains('>') {
            self.parse_substitution(&chr, rest)
        } else if rest.contains("delins") || (rest.contains("del") && rest.contains("ins")) {
            self.parse_delins(&chr, rest)
        } else if rest.contains("ins") {
            self.parse_insertion(&chr, rest)
        } else if rest.contains("del") {
            self.parse_deletion(&chr, rest)
        } else if rest.contains("dup") {
            self.parse_duplication(&chr, rest)
        } else {
            Err(IoError::VcfParse(format!(
                "Unrecognized HGVS type in: '{}'",
                line
            )))
        }
    }

    /// Parse a substitution like "25585733A>G".
    fn parse_substitution(&self, chr: &str, rest: &str) -> Result<Option<InputVariant>, IoError> {
        let gt_pos = rest
            .find('>')
            .ok_or_else(|| IoError::VcfParse(format!("No '>' in substitution: '{}'", rest)))?;

        let before_gt = &rest[..gt_pos];
        let alt = &rest[gt_pos + 1..];

        let pos_end = before_gt
            .find(|c: char| c.is_ascii_alphabetic())
            .unwrap_or(before_gt.len());

        let pos_str = &before_gt[..pos_end];
        let ref_allele = &before_gt[pos_end..];

        let pos: u64 = pos_str
            .parse()
            .map_err(|_| IoError::VcfParse(format!("Invalid position: '{}'", pos_str)))?;

        if ref_allele.is_empty() || alt.is_empty() {
            return Err(IoError::VcfParse(format!(
                "Missing ref or alt allele in substitution: '{}'",
                rest
            )));
        }

        let variant = InputVariant::new(
            chr.to_string(),
            pos,
            pos,
            ref_allele.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );

        Ok(Some(variant))
    }

    /// Parse a deletion like "25585656_25585660del" or "25585656_25585660delACGTG".
    fn parse_deletion(&self, chr: &str, rest: &str) -> Result<Option<InputVariant>, IoError> {
        let del_pos = rest
            .find("del")
            .ok_or_else(|| IoError::VcfParse(format!("No 'del' in deletion: '{}'", rest)))?;

        let pos_part = &rest[..del_pos];
        let deleted_seq = &rest[del_pos + 3..]; // After "del"

        let (start, end) = self.parse_position_range(pos_part)?;

        let ref_allele = if deleted_seq.is_empty() {
            b"-".to_vec()
        } else {
            deleted_seq.as_bytes().to_vec()
        };

        let variant = InputVariant::new(chr.to_string(), start, end, ref_allele, b"-".to_vec());

        Ok(Some(variant))
    }

    /// Parse an insertion like "25585734_25585735insACG".
    fn parse_insertion(&self, chr: &str, rest: &str) -> Result<Option<InputVariant>, IoError> {
        let ins_pos = rest
            .find("ins")
            .ok_or_else(|| IoError::VcfParse(format!("No 'ins' in insertion: '{}'", rest)))?;

        let pos_part = &rest[..ins_pos];
        let inserted_seq = &rest[ins_pos + 3..]; // After "ins"

        if inserted_seq.is_empty() {
            return Err(IoError::VcfParse(format!(
                "Missing inserted sequence in: '{}'",
                rest
            )));
        }

        let (start, end) = self.parse_position_range(pos_part)?;

        let variant = InputVariant::new(
            chr.to_string(),
            start,
            end,
            b"-".to_vec(),
            inserted_seq.as_bytes().to_vec(),
        );

        Ok(Some(variant))
    }

    /// Parse a delins like "25585733delAinsGT".
    fn parse_delins(&self, chr: &str, rest: &str) -> Result<Option<InputVariant>, IoError> {
        let del_pos = rest
            .find("del")
            .ok_or_else(|| IoError::VcfParse(format!("No 'del' in delins: '{}'", rest)))?;
        let ins_pos = rest
            .find("ins")
            .ok_or_else(|| IoError::VcfParse(format!("No 'ins' in delins: '{}'", rest)))?;

        let pos_part = &rest[..del_pos];
        let deleted_seq = &rest[del_pos + 3..ins_pos];
        let inserted_seq = &rest[ins_pos + 3..];

        if inserted_seq.is_empty() {
            return Err(IoError::VcfParse(format!(
                "Missing inserted sequence in delins: '{}'",
                rest
            )));
        }

        let (start, end) = self.parse_position_range(pos_part)?;

        let ref_allele = if deleted_seq.is_empty() {
            b"-".to_vec()
        } else {
            deleted_seq.as_bytes().to_vec()
        };

        let variant = InputVariant::new(
            chr.to_string(),
            start,
            end,
            ref_allele,
            inserted_seq.as_bytes().to_vec(),
        );

        Ok(Some(variant))
    }

    /// Parse a duplication like "25585733dupA".
    fn parse_duplication(&self, chr: &str, rest: &str) -> Result<Option<InputVariant>, IoError> {
        let dup_pos = rest
            .find("dup")
            .ok_or_else(|| IoError::VcfParse(format!("No 'dup' in duplication: '{}'", rest)))?;

        let pos_part = &rest[..dup_pos];
        let dup_seq = &rest[dup_pos + 3..]; // After "dup"

        let (start, end) = self.parse_position_range(pos_part)?;

        // A duplication is effectively an insertion of the duplicated sequence
        let alt = if dup_seq.is_empty() {
            b"-".to_vec()
        } else {
            dup_seq.as_bytes().to_vec()
        };

        let variant = InputVariant::new(chr.to_string(), start, end, b"-".to_vec(), alt);

        Ok(Some(variant))
    }

    /// Parse a position or position range from a string.
    /// Handles both "123" and "123_456".
    fn parse_position_range(&self, s: &str) -> Result<(u64, u64), IoError> {
        if let Some((start_str, end_str)) = s.split_once('_') {
            let start: u64 = start_str.parse().map_err(|_| {
                IoError::VcfParse(format!("Invalid start position: '{}'", start_str))
            })?;
            let end: u64 = end_str
                .parse()
                .map_err(|_| IoError::VcfParse(format!("Invalid end position: '{}'", end_str)))?;
            Ok((start, end))
        } else {
            let pos: u64 = s
                .parse()
                .map_err(|_| IoError::VcfParse(format!("Invalid position: '{}'", s)))?;
            Ok((pos, pos))
        }
    }
}

impl<R: BufRead + Send> InputParser for HgvsParser<R> {
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
    fn test_parse_substitution() {
        let data = "21:g.25585733A>G\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585733);
        assert_eq!(v.end, 25585733);
        assert_eq!(v.ref_allele, b"A");
        assert_eq!(v.alt_allele(), b"G");
    }

    #[test]
    fn test_parse_deletion() {
        let data = "21:g.25585656_25585660del\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585656);
        assert_eq!(v.end, 25585660);
        assert_eq!(v.alt_allele(), b"-");
    }

    #[test]
    fn test_parse_deletion_with_seq() {
        let data = "21:g.25585656_25585660delACGTG\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.ref_allele, b"ACGTG");
        assert_eq!(v.alt_allele(), b"-");
    }

    #[test]
    fn test_parse_delins() {
        let data = "21:g.25585733delAinsGT\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585733);
        assert_eq!(v.ref_allele, b"A");
        assert_eq!(v.alt_allele(), b"GT");
    }

    #[test]
    fn test_parse_insertion() {
        let data = "21:g.25585734_25585735insACG\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585734);
        assert_eq!(v.end, 25585735);
        assert_eq!(v.ref_allele, b"-");
        assert_eq!(v.alt_allele(), b"ACG");
    }

    #[test]
    fn test_transcript_hgvs_error() {
        let data = "ENST00000366667:c.803C>T\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let result = parser.next_variant().unwrap();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("transcript HGVS is not supported"));
    }

    #[test]
    fn test_skip_comments_and_empty() {
        let data = "# comment\n\n21:g.25585733A>G\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");

        assert!(parser.next_variant().is_none());
    }

    #[test]
    fn test_parse_duplication() {
        let data = "21:g.25585733dupA\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "21");
        assert_eq!(v.start, 25585733);
        assert_eq!(v.alt_allele(), b"A");
    }

    #[test]
    fn test_invalid_notation() {
        let data = "21:25585733\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let result = parser.next_variant().unwrap();
        assert!(result.is_err());
    }

    #[test]
    fn test_x_chromosome() {
        let data = "X:g.100000A>T\n";
        let reader = Cursor::new(data.as_bytes());
        let mut parser = HgvsParser::new(reader);

        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.chr, "X");
        assert_eq!(v.start, 100000);
    }
}
