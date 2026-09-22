// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! SIFT/PolyPhen prediction matrix decoder.
//!
//! Ensembl VEP stores pre-computed SIFT and PolyPhen predictions as compressed
//! binary matrices per transcript. Each matrix encodes predictions for every
//! possible amino acid substitution at every protein position.
//!
//! ## Binary Format (from `ProteinFunctionPredictionMatrix.pm` in ensembl-variation release/115)
//!
//! The matrix is gzip-compressed. After decompression:
//!
//! - **Header**: 3 bytes `"VEP"`
//! - **Per position**: 20 amino acids × 2 bytes = 40 bytes
//! - **Cell**: `u16` little-endian (Perl `pack 'v'`)
//!   - Top 2 bits (15-14): prediction category index. Ensembl reserves
//!     `NUM_PRED_BITS = ceil(log2(4)) = 2` bits for the four categories of the
//!     largest tool and shifts by `16 - 2`; bits 13-10 are unused.
//!   - Bottom 10 bits (9-0): score × 1000 (0-1000 → 0.000-1.000)
//!   - `0xFFFF`: no prediction available
//! - **Amino acid order**: A C D E F G H I K L M N P Q R S T V W Y
//!
//! ## Data Source
//!
//! Prediction matrices come from Ensembl's variation database:
//! - Public MySQL: `ensembldb.ensembl.org:3306` (anonymous, read-only)
//! - Tables: `translation_md5`, `protein_function_predictions`, `attrib`
//! - Keyed by MD5 hash of the protein sequence

use std::io::Read;

use flate2::read::GzDecoder;

/// Header bytes at the start of every prediction matrix.
const HEADER: &[u8] = b"VEP";

/// The 20 standard amino acids in Ensembl matrix column order.
const AMINO_ACIDS: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Number of amino acids per position.
const AA_COUNT: usize = 20;

/// Bytes per amino acid cell (u16 little-endian).
const BYTES_PER_CELL: usize = 2;

/// Sentinel value indicating no prediction is available.
const NO_PREDICTION: u16 = 0xFFFF;

/// Mask for the score (bottom 10 bits), Perl's `$val & (2**10 - 1)`.
const SCORE_MASK: u16 = 0x03FF;

/// Number of bits to shift right to get the category index: Perl's
/// `$val >> (16 - $NUM_PRED_BITS)` with `NUM_PRED_BITS = 2`.
const CATEGORY_SHIFT: u32 = 14;

/// Divisor to convert raw score to 0.000-1.000 range.
const SCORE_DIVISOR: f64 = 1000.0;

/// Prediction analysis type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalysisType {
    Sift,
    PolyPhenHumDiv,
    PolyPhenHumVar,
}

/// SIFT prediction categories (indexed by the top 2 bits of the cell).
const SIFT_CATEGORIES: &[&str] = &[
    "tolerated",
    "deleterious",
    "tolerated_low_confidence",
    "deleterious_low_confidence",
];

/// PolyPhen prediction categories (indexed by the top 2 bits of the cell).
const POLYPHEN_CATEGORIES: &[&str] = &[
    "probably_damaging",
    "possibly_damaging",
    "benign",
    "unknown",
];

/// A decoded prediction: category label and numeric score.
#[derive(Debug, Clone, PartialEq)]
pub struct Prediction {
    /// Categorical prediction (e.g., "tolerated", "probably_damaging").
    pub label: String,
    /// Numeric score (0.0 to 1.0, 3 decimal places).
    pub score: f64,
}

/// Decompress a gzip-compressed prediction matrix.
///
/// Returns the raw bytes after decompression (header + position data).
pub fn gunzip_matrix(compressed: &[u8]) -> Result<Vec<u8>, std::io::Error> {
    let mut decoder = GzDecoder::new(compressed);
    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed)?;
    Ok(decompressed)
}

/// Decode a prediction from a decompressed binary matrix.
///
/// # Arguments
/// - `data`: Decompressed matrix bytes (must start with "VEP" header).
/// - `position`: Protein position (1-based).
/// - `amino_acid`: Single-letter amino acid code (e.g., `b'V'`).
/// - `analysis`: Which prediction type to decode categories for.
///
/// # Returns
/// `Some(Prediction)` if a prediction exists, `None` if the position/AA is out of
/// range or no prediction is available (0xFFFF sentinel).
pub fn decode_prediction(
    data: &[u8],
    position: usize,
    amino_acid: u8,
    analysis: AnalysisType,
) -> Option<Prediction> {
    if data.len() < HEADER.len() || &data[..HEADER.len()] != HEADER {
        return None;
    }

    let aa_idx = AMINO_ACIDS.iter().position(|&b| b == amino_acid)?;

    let offset = HEADER.len() + ((position - 1) * AA_COUNT + aa_idx) * BYTES_PER_CELL;
    if offset + 1 >= data.len() {
        return None;
    }

    let raw = u16::from_le_bytes([data[offset], data[offset + 1]]);
    if raw == NO_PREDICTION {
        return None;
    }

    let cat_idx = (raw >> CATEGORY_SHIFT) as usize;
    let score = f64::from(raw & SCORE_MASK) / SCORE_DIVISOR;

    let categories = match analysis {
        AnalysisType::Sift => SIFT_CATEGORIES,
        AnalysisType::PolyPhenHumDiv | AnalysisType::PolyPhenHumVar => POLYPHEN_CATEGORIES,
    };

    let label = categories.get(cat_idx)?.to_string();
    Some(Prediction { label, score })
}

/// Format a prediction for VEP output.
///
/// The score is printed the way Perl stringifies `$prob / 1000`: no fixed number of
/// decimals, so 10 prints as `0.01`, 1000 as `1` and 796 as `0.796`.
///
/// # Mode values
/// - `"p"`: prediction label only (e.g., `"tolerated"`)
/// - `"s"`: score only (e.g., `"0.456"`)
/// - `"b"` or anything else: both (e.g., `"tolerated(0.456)"`)
pub fn format_prediction(prediction: &Prediction, mode: &str) -> String {
    match mode {
        "p" => prediction.label.clone(),
        "s" => format_score(prediction.score),
        _ => format!("{}({})", prediction.label, format_score(prediction.score)),
    }
}

/// A score as Perl prints a number: shortest round-trip decimal, no trailing zeros.
fn format_score(score: f64) -> String {
    format!("{score}")
}

/// Map an analysis name string to the corresponding [`AnalysisType`].
pub fn analysis_from_name(name: &str) -> Option<AnalysisType> {
    match name {
        "sift" => Some(AnalysisType::Sift),
        "polyphen_humdiv" => Some(AnalysisType::PolyPhenHumDiv),
        "polyphen_humvar" => Some(AnalysisType::PolyPhenHumVar),
        _ => None,
    }
}

/// Build a synthetic prediction matrix for testing.
///
/// Creates a valid gzip-compressed matrix with the given predictions.
/// Each entry is `(position, amino_acid, category_index, score_x1000)`.
#[cfg(test)]
fn build_test_matrix(peptide_length: usize, entries: &[(usize, u8, u16, u16)]) -> Vec<u8> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    let data_len = HEADER.len() + peptide_length * AA_COUNT * BYTES_PER_CELL;
    let mut raw = vec![0xFF; data_len]; // fill with 0xFF (no prediction)

    raw[..HEADER.len()].copy_from_slice(HEADER);

    for &(pos, aa, cat, score) in entries {
        if let Some(aa_idx) = AMINO_ACIDS.iter().position(|&b| b == aa) {
            let offset = HEADER.len() + ((pos - 1) * AA_COUNT + aa_idx) * BYTES_PER_CELL;
            let value = (cat << CATEGORY_SHIFT) | (score & SCORE_MASK);
            let bytes = value.to_le_bytes();
            raw[offset] = bytes[0];
            raw[offset + 1] = bytes[1];
        }
    }

    // Gzip compress
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(&raw).unwrap();
    encoder.finish().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_sift_tolerated() {
        let compressed = build_test_matrix(
            5,
            &[
                // position 3, amino acid V (index 17), category 0 (tolerated), score 456
                (3, b'V', 0, 456),
            ],
        );
        let data = gunzip_matrix(&compressed).unwrap();
        let pred = decode_prediction(&data, 3, b'V', AnalysisType::Sift).unwrap();
        assert_eq!(pred.label, "tolerated");
        assert!((pred.score - 0.456).abs() < 1e-9);
    }

    #[test]
    fn test_decode_sift_deleterious() {
        let compressed = build_test_matrix(
            5,
            &[
                // position 1, amino acid A (index 0), category 1 (deleterious), score 10
                (1, b'A', 1, 10),
            ],
        );
        let data = gunzip_matrix(&compressed).unwrap();
        let pred = decode_prediction(&data, 1, b'A', AnalysisType::Sift).unwrap();
        assert_eq!(pred.label, "deleterious");
        assert!((pred.score - 0.010).abs() < 1e-9);
    }

    #[test]
    fn test_decode_polyphen_probably_damaging() {
        let compressed = build_test_matrix(
            10,
            &[
                // position 7, amino acid R (index 14), category 0 (probably_damaging), score 950
                (7, b'R', 0, 950),
            ],
        );
        let data = gunzip_matrix(&compressed).unwrap();
        let pred = decode_prediction(&data, 7, b'R', AnalysisType::PolyPhenHumVar).unwrap();
        assert_eq!(pred.label, "probably_damaging");
        assert!((pred.score - 0.950).abs() < 1e-9);
    }

    #[test]
    fn test_decode_polyphen_benign() {
        let compressed = build_test_matrix(
            5,
            &[
                (2, b'D', 2, 50), // category 2 (benign), score 0.050
            ],
        );
        let data = gunzip_matrix(&compressed).unwrap();
        let pred = decode_prediction(&data, 2, b'D', AnalysisType::PolyPhenHumDiv).unwrap();
        assert_eq!(pred.label, "benign");
        assert!((pred.score - 0.050).abs() < 1e-9);
    }

    #[test]
    fn test_no_prediction_ffff() {
        let compressed = build_test_matrix(5, &[]); // all cells are 0xFFFF
        let data = gunzip_matrix(&compressed).unwrap();
        let result = decode_prediction(&data, 1, b'A', AnalysisType::Sift);
        assert!(result.is_none());
    }

    #[test]
    fn test_out_of_bounds_position() {
        let compressed = build_test_matrix(3, &[(1, b'A', 0, 500)]);
        let data = gunzip_matrix(&compressed).unwrap();
        // Position 10 is beyond a 3-residue protein
        let result = decode_prediction(&data, 10, b'A', AnalysisType::Sift);
        assert!(result.is_none());
    }

    #[test]
    fn test_invalid_amino_acid() {
        let compressed = build_test_matrix(5, &[(1, b'A', 0, 500)]);
        let data = gunzip_matrix(&compressed).unwrap();
        // 'X' is not in the standard 20 amino acids
        let result = decode_prediction(&data, 1, b'X', AnalysisType::Sift);
        assert!(result.is_none());
    }

    #[test]
    fn test_invalid_header() {
        let result = decode_prediction(b"BAD", 1, b'A', AnalysisType::Sift);
        assert!(result.is_none());
    }

    #[test]
    fn test_empty_data() {
        let result = decode_prediction(&[], 1, b'A', AnalysisType::Sift);
        assert!(result.is_none());
    }

    #[test]
    fn test_format_prediction_both() {
        let pred = Prediction {
            label: "tolerated".to_string(),
            score: 0.456,
        };
        assert_eq!(format_prediction(&pred, "b"), "tolerated(0.456)");
    }

    #[test]
    fn test_format_prediction_label_only() {
        let pred = Prediction {
            label: "deleterious".to_string(),
            score: 0.01,
        };
        assert_eq!(format_prediction(&pred, "p"), "deleterious");
    }

    #[test]
    fn test_format_prediction_score_only() {
        let pred = Prediction {
            label: "probably_damaging".to_string(),
            score: 0.95,
        };
        assert_eq!(format_prediction(&pred, "s"), "0.95");
    }

    /// Perl prints `$prob / 1000` as a bare number, so a whole score and a two-place
    /// score carry no padding zeros.
    #[test]
    fn test_format_prediction_matches_perl_number_stringification() {
        let one = Prediction {
            label: "tolerated".to_string(),
            score: 1.0,
        };
        assert_eq!(format_prediction(&one, "b"), "tolerated(1)");
        let hundredth = Prediction {
            label: "deleterious".to_string(),
            score: 0.01,
        };
        assert_eq!(format_prediction(&hundredth, "b"), "deleterious(0.01)");
        let zero = Prediction {
            label: "deleterious".to_string(),
            score: 0.0,
        };
        assert_eq!(format_prediction(&zero, "s"), "0");
    }

    /// A cell packed by Ensembl's own rule (`$val = $prob * 1000; $val |= $pred_val << 14`)
    /// decodes to the category and score it was packed from. The two examples are the
    /// ones ProteinFunctionPredictionMatrix.pm's synopsis adds: probably damaging 0.967
    /// at position 1 A, benign 0.09 at position 2 C.
    #[test]
    fn test_decode_cells_packed_by_perl_rule() {
        let mut raw = vec![0xFFu8; HEADER.len() + 2 * AA_COUNT * BYTES_PER_CELL];
        raw[..HEADER.len()].copy_from_slice(HEADER);
        let cell_a: u16 = 967; // probably damaging = 0 << 14
        let cell_c: u16 = 90 | (2 << 14); // benign = 2
        raw[HEADER.len()..HEADER.len() + 2].copy_from_slice(&cell_a.to_le_bytes());
        let off_c = HEADER.len() + (AA_COUNT + 1) * BYTES_PER_CELL;
        raw[off_c..off_c + 2].copy_from_slice(&cell_c.to_le_bytes());
        let a = decode_prediction(&raw, 1, b'A', AnalysisType::PolyPhenHumVar).unwrap();
        assert_eq!(a.label, "probably_damaging");
        assert!((a.score - 0.967).abs() < 1e-9);
        let c = decode_prediction(&raw, 2, b'C', AnalysisType::PolyPhenHumVar).unwrap();
        assert_eq!(c.label, "benign");
        assert!((c.score - 0.09).abs() < 1e-9);
        assert_eq!(format_prediction(&c, "b"), "benign(0.09)");
    }

    #[test]
    fn test_gunzip_roundtrip() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let original = b"VEP\x00\x01\x02\x03\x04\x05";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(original).unwrap();
        let compressed = encoder.finish().unwrap();

        let decompressed = gunzip_matrix(&compressed).unwrap();
        assert_eq!(&decompressed, original);
    }

    #[test]
    fn test_analysis_from_name() {
        assert_eq!(analysis_from_name("sift"), Some(AnalysisType::Sift));
        assert_eq!(
            analysis_from_name("polyphen_humvar"),
            Some(AnalysisType::PolyPhenHumVar)
        );
        assert_eq!(
            analysis_from_name("polyphen_humdiv"),
            Some(AnalysisType::PolyPhenHumDiv)
        );
        assert_eq!(analysis_from_name("unknown"), None);
    }

    #[test]
    fn test_all_amino_acids_indexed() {
        // Verify every standard AA has a valid index
        for &aa in b"ACDEFGHIKLMNPQRSTVWY" {
            assert!(
                AMINO_ACIDS.iter().position(|&b| b == aa).is_some(),
                "AA '{}' not found in index",
                aa as char
            );
        }
    }

    #[test]
    fn test_multiple_positions_and_aas() {
        let compressed = build_test_matrix(
            5,
            &[
                (1, b'A', 0, 100), // pos 1, Ala, tolerated(0.100)
                (1, b'V', 1, 5),   // pos 1, Val, deleterious(0.005)
                (3, b'G', 0, 999), // pos 3, Gly, tolerated(0.999)
                (5, b'Y', 3, 0),   // pos 5, Tyr, deleterious_low_confidence(0.000)
            ],
        );
        let data = gunzip_matrix(&compressed).unwrap();

        let p1 = decode_prediction(&data, 1, b'A', AnalysisType::Sift).unwrap();
        assert_eq!(p1.label, "tolerated");
        assert!((p1.score - 0.100).abs() < 1e-9);

        let p2 = decode_prediction(&data, 1, b'V', AnalysisType::Sift).unwrap();
        assert_eq!(p2.label, "deleterious");
        assert!((p2.score - 0.005).abs() < 1e-9);

        let p3 = decode_prediction(&data, 3, b'G', AnalysisType::Sift).unwrap();
        assert_eq!(p3.label, "tolerated");
        assert!((p3.score - 0.999).abs() < 1e-9);

        let p4 = decode_prediction(&data, 5, b'Y', AnalysisType::Sift).unwrap();
        assert_eq!(p4.label, "deleterious_low_confidence");
        assert!((p4.score - 0.0).abs() < 1e-9);

        // Unpopulated cell should return None
        let p5 = decode_prediction(&data, 2, b'C', AnalysisType::Sift);
        assert!(p5.is_none());
    }
}
