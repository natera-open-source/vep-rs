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
//!
//! The format carries no reference allele. Ensembl VEP 116.2 annotates a line
//! covering bases (`start <= end`) against the genome's bases over the span, so
//! a parser given the reference FASTA ([`RegionParser::with_reference`]) reads
//! them as the reference allele, on the line's strand; an insertion (`start > end`) has `-`. Without
//! the bases, a deletion keeps one `N` per deleted base, so it still reads as a
//! deletion of its span's length, and any other allele keeps `-`. The allele
//! string every output names the record by stays `-/allele` either way.

use std::io::BufRead;
use std::sync::Arc;

use vep_core::coordinate::{parse_coordinate, Strand};
use vep_core::variant::InputVariant;
use vep_fasta::IndexedFasta;

use crate::error::IoError;

use super::InputParser;

/// The longest span a deletion without reference bases is read over: its `N`
/// placeholder is allocated base by base, so a line spanning more than 4.3 Gb,
/// over ten times the longest human chromosome, fails rather than allocating it.
const MAX_PLACEHOLDER_SPAN: u64 = u32::MAX as u64;

/// Parser for region format input.
pub struct RegionParser<R: BufRead> {
    reader: R,
    line_buf: String,
    reference: Option<Arc<IndexedFasta>>,
}

impl<R: BufRead> RegionParser<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: String::new(),
            reference: None,
        }
    }

    /// Read the reference allele of every line covering bases from `reference`.
    pub fn with_reference(mut self, reference: Option<Arc<IndexedFasta>>) -> Self {
        self.reference = reference;
        self
    }

    /// The reference allele of a line over `start..=end` with `allele`: `-` for
    /// an insertion; the genome's bases over the span when the FASTA holds them;
    /// else one `N` per base for a deletion and `-` for any other allele.
    fn reference_allele(
        &self,
        chr: &str,
        start: u64,
        end: u64,
        allele: &[u8],
    ) -> Result<Vec<u8>, IoError> {
        if start > end {
            return Ok(b"-".to_vec());
        }
        if let Some(bases) = self
            .reference
            .as_deref()
            .and_then(|fasta| fasta.sequence(chr, start, end))
        {
            return Ok(bases);
        }
        if allele != b"-" {
            return Ok(b"-".to_vec());
        }
        let span = end - start + 1;
        if span > MAX_PLACEHOLDER_SPAN {
            return Err(IoError::VcfParse(format!(
                "Deletion {chr}:{start}-{end} spans {span} bases, more than {MAX_PLACEHOLDER_SPAN} \
                 with no reference bases to read"
            )));
        }
        Ok(vec![b'N'; span as usize])
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
            let start = parse_coordinate(s)
                .ok_or_else(|| IoError::VcfParse(format!("Invalid start position: '{}'", s)))?;
            let end = parse_coordinate(e)
                .ok_or_else(|| IoError::VcfParse(format!("Invalid end position: '{}'", e)))?;
            (start, end)
        } else {
            let pos = parse_coordinate(range_part)
                .ok_or_else(|| IoError::VcfParse(format!("Invalid position: '{}'", range_part)))?;
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

        let mut ref_allele = self.reference_allele(&chr, start, end, allele.as_bytes())?;
        // VariationFeatureOverlap.pm reads the bases on the line's own strand
        // (`sub_Slice($start, $end, $vf->strand)`), the strand its allele is on.
        if strand == Strand::Reverse && ref_allele != b"-" {
            ref_allele = vep_core::codon::reverse_complement(&ref_allele);
        }
        let mut variant =
            InputVariant::new(chr, start, end, ref_allele, allele.as_bytes().to_vec());
        // Outputs name the record by the line's own alleles, which carry no reference.
        variant.allele_string = format!("-/{allele}");
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
    use vep_core::coordinate::MAX_COORD;
    use vep_core::variant::VariantClass;

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

    /// A coordinate is read up to `MAX_COORD`; a start, end or position past it
    /// fails the line.
    #[test]
    fn coordinates_past_the_ceiling_fail_the_line() {
        let data = "21:4611686018427387903-4611686018427387903:1/G\n\
                    21:4611686018427387904-4611686018427387904:1/G\n\
                    21:100-18446744073709551615:1/G\n\
                    21:9223372036854775807:1/G\n";
        let mut parser = RegionParser::new(Cursor::new(data.as_bytes()));
        let kept = parser.next_variant().unwrap().unwrap();
        assert_eq!((kept.start, kept.end), (MAX_COORD, MAX_COORD));
        for expected in [
            "Invalid start position",
            "Invalid end position",
            "Invalid position",
        ] {
            let err = parser.next_variant().unwrap().unwrap_err().to_string();
            assert!(err.contains(expected), "{err}");
        }
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

    /// Contig `21` reading `ACGTACGTAC` from position 1, and the handle over it.
    fn reference_fixture() -> (tempfile::TempDir, Arc<IndexedFasta>) {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, ">21\nACGTACGTAC\n").unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "21\t10\t4\t10\t11\n").unwrap();
        let fasta = Arc::new(IndexedFasta::from_path(&fa).unwrap());
        (dir, fasta)
    }

    /// (ref_allele, variant_class, allele_string, Uploaded_variation) of each line.
    fn read_all(
        data: &str,
        reference: Option<Arc<IndexedFasta>>,
    ) -> Vec<(Vec<u8>, VariantClass, String, String)> {
        let mut parser = RegionParser::new(Cursor::new(data.as_bytes())).with_reference(reference);
        std::iter::from_fn(|| parser.next_variant())
            .map(|v| {
                let v = v.unwrap();
                let name = v.uploaded_variation();
                (v.ref_allele, v.variant_class, v.allele_string, name)
            })
            .collect()
    }

    /// With the reference, a line covering bases reads the genome's bases over
    /// its span as the reference allele, so `-` is a deletion of them; an
    /// insertion keeps `-`. The record is still named by the line's `-/allele`.
    #[test]
    fn reads_the_genome_bases_over_the_span_as_reference() {
        let (_dir, fasta) = reference_fixture();
        let data = "21:3-4:1/-\n21:2-2:1/T\n21:5-4:1/AA\n21:2-4/A\nchr21:10:1/-\n";
        let rows = read_all(data, Some(fasta));
        let expected: [(&[u8], VariantClass, &str, &str); 5] = [
            (b"GT", VariantClass::Deletion, "-/-", "21_3_-/-"),
            (b"C", VariantClass::Snv, "-/T", "21_2_-/T"),
            (b"-", VariantClass::Insertion, "-/AA", "21_5_-/AA"),
            (b"CGT", VariantClass::Indel, "-/A", "21_2_-/A"),
            (b"C", VariantClass::Deletion, "-/-", "chr21_10_-/-"),
        ];
        for (row, (ref_allele, class, allele_string, name)) in rows.iter().zip(expected) {
            assert_eq!(
                (row.0.as_slice(), row.1, row.2.as_str(), row.3.as_str()),
                (ref_allele, class, allele_string, name)
            );
        }
        assert_eq!(rows.len(), 5);
    }

    /// A `-1` line reads the bases on its own strand, as Ensembl VEP's
    /// `sub_Slice($start, $end, $vf->strand)` does, so the reference sits on the
    /// same strand as the line's allele.
    #[test]
    fn a_reverse_strand_line_reads_the_reverse_complement() {
        let (_dir, fasta) = reference_fixture();
        let rows = read_all("21:3-5:-1/-\n21:2-2:-1/A\n21:5-4:-1/AA\n", Some(fasta));
        let refs: Vec<(&[u8], VariantClass)> = rows.iter().map(|r| (r.0.as_slice(), r.1)).collect();
        assert_eq!(
            refs,
            [
                (&b"TAC"[..], VariantClass::Deletion),
                (b"G", VariantClass::Snv),
                (b"-", VariantClass::Insertion),
            ]
        );
        let without = read_all("21:3-4:-1/-\n", None);
        assert_eq!(without[0].0, b"NN");
    }

    /// Without the bases (no FASTA, or a contig or span it lacks) a deletion
    /// reads one `N` per base and stays a deletion; any other allele keeps `-`.
    /// A deletion too long to hold its placeholder fails the line.
    #[test]
    fn a_deletion_without_reference_bases_reads_n_per_base() {
        let data = "21:3-4:1/-\n21:2-2:1/T\n";
        let without = read_all(data, None);
        assert_eq!(
            (without[0].0.as_slice(), without[0].1),
            (&b"NN"[..], VariantClass::Deletion)
        );
        assert_eq!(without[1].0, b"-");

        let (_dir, fasta) = reference_fixture();
        let off_reference = read_all("22:3-4:1/-\n21:10-12:1/-\n21:11-11:1/T\n", Some(fasta));
        let refs: Vec<&[u8]> = off_reference.iter().map(|r| r.0.as_slice()).collect();
        assert_eq!(refs, [&b"NN"[..], b"NNN", b"-"]);

        let mut parser = RegionParser::new(Cursor::new("21:1-4294967296:1/-\n".as_bytes()));
        let err = parser.next_variant().unwrap().unwrap_err().to_string();
        assert!(err.contains("spans 4294967296 bases"), "{err}");
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
