// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parse Ensembl protein FASTA files to extract transcript -> protein sequence mapping.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

use crate::core_dumps::BoundedLines;

/// The longest protein FASTA line read, 16 MiB. Ensembl wraps sequences at 60
/// residues, and the longest human protein written on one line runs to some
/// 35,000, so the cap refuses a file with no line end in it long before the
/// line buffer grows to its size.
const MAX_LINE_BYTES: usize = 16 << 20;

/// Parse a protein FASTA file (optionally gzipped) into a map of ENST -> sequence,
/// with the length of a line capped: a line longer than `MAX_LINE_BYTES` fails
/// the parse, naming the file and the line.
///
/// Ensembl protein FASTA headers look like:
/// ```text
/// >ENSP00000123456.2 pep chromosome:GRCh38:1:100:500:1 gene:ENSG00000123 transcript:ENST00000456 ...
/// MKVLWAALLLLAAMYTISVQTPMGKSRCMKDRHGAFEEHG...
/// ```
///
/// The transcript ID (ENST...) is extracted from the header and mapped to the
/// amino acid sequence.
pub fn parse_protein_fasta(path: &Path) -> Result<HashMap<String, String>> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;

    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut sequences: HashMap<String, String> = HashMap::new();
    let mut current_transcript: Option<String> = None;
    let mut current_seq = String::new();

    for line in BoundedLines::new(reader, MAX_LINE_BYTES, path.display().to_string()) {
        let line = String::from_utf8(line?).context("Failed to read FASTA line")?;
        let line = line.trim();

        if line.starts_with('>') {
            if let Some(ref tid) = current_transcript {
                if !current_seq.is_empty() {
                    sequences.insert(tid.clone(), current_seq.clone());
                }
            }

            current_transcript = extract_transcript_id(line);
            current_seq.clear();
        } else if !line.is_empty() {
            current_seq.push_str(line);
        }
    }

    if let Some(ref tid) = current_transcript {
        if !current_seq.is_empty() {
            sequences.insert(tid.clone(), current_seq);
        }
    }

    Ok(sequences)
}

/// Extract the ENST transcript ID from a FASTA header line.
///
/// Looks for a `transcript:ENST<digits>` token in the header. Returns the
/// ENST ID without version suffix.
fn extract_transcript_id(header: &str) -> Option<String> {
    for part in header.split_whitespace() {
        if let Some(tid) = part.strip_prefix("transcript:") {
            // Strip version suffix (e.g., "ENST00000456.3" -> "ENST00000456").
            let base_id = tid.split('.').next().unwrap_or(tid);
            return Some(base_id.to_string());
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_transcript_id() {
        let header = ">ENSP00000123456.2 pep chromosome:GRCh38:1:100:500:1 gene:ENSG00000123 transcript:ENST00000456.3 gene_biotype:protein_coding";
        assert_eq!(
            extract_transcript_id(header),
            Some("ENST00000456".to_string())
        );
    }

    #[test]
    fn test_extract_transcript_id_no_version() {
        let header = ">ENSP00000123456 pep transcript:ENST00000789";
        assert_eq!(
            extract_transcript_id(header),
            Some("ENST00000789".to_string())
        );
    }

    #[test]
    fn test_extract_transcript_id_missing() {
        let header = ">ENSP00000123456 pep gene:ENSG00000123";
        assert_eq!(extract_transcript_id(header), None);
    }

    /// A sequence line one byte over the cap fails the parse, naming the file
    /// and the line, 2; a file within the cap parses.
    #[test]
    fn line_past_the_cap_fails_the_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pep.fa");
        let header = ">ENSP00000123456.2 pep transcript:ENST00000456.3\n";
        std::fs::write(&path, format!("{header}MKV\n")).unwrap();
        assert_eq!(parse_protein_fasta(&path).unwrap()["ENST00000456"], "MKV");
        let mut text = header.as_bytes().to_vec();
        text.extend(std::iter::repeat_n(b'M', MAX_LINE_BYTES + 1));
        text.push(b'\n');
        std::fs::write(&path, text).unwrap();
        let err = parse_protein_fasta(&path).unwrap_err();
        assert!(
            format!("{err:#}").ends_with(&format!(
                "pep.fa line 2 is longer than the {MAX_LINE_BYTES}-byte limit for one line"
            )),
            "{err:#}"
        );
    }
}
