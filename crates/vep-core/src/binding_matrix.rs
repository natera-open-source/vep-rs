// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A transcription factor binding matrix as the regulatory cache stores it
//! (Ensembl funcgen `BindingMatrix`): nucleotide counts by position, the
//! matrix's descriptive fields and the complexes that bind it, with the
//! scoring tables funcgen derives from the counts
//! (`BindingMatrix/Converter.pm`) and the two scores Ensembl VEP reads from
//! them (`BindingMatrix.pm`): a sequence's similarity to the matrix relative to
//! its best and worst sites, and whether a position is informative.
//!
//! The cache keeps every matrix once in `regulatory/binding_matrices.json`,
//! keyed by stable id, and each motif feature names its matrix by that id, so
//! the tables are derived once per matrix and shared by every motif.

use std::sync::{Arc, OnceLock};

/// The nucleotide order of every position's counts: A, C, G, T.
pub const NUCLEOTIDES: [u8; 4] = *b"ACGT";

/// The information content, in bits, a position must exceed to be informative
/// (the default threshold of `is_position_informative`, `BindingMatrix.pm:631`).
pub const INFORMATIVE_BITS: f64 = 1.5;

/// The pseudocount added to every count before a probability is taken
/// (`Converter.pm:116`).
const PSEUDOCOUNT: f64 = 0.1;

/// The background frequency every weight is taken against (`Converter.pm:215`).
const EXPECTED_FREQUENCY: f64 = 0.25;

/// A binding matrix of the regulatory cache. `frequencies` holds one row of
/// A, C, G, T counts per position, positions 1 to `length` in order, and is
/// empty for a matrix known only by its id and length (a motif source without
/// counts, whose motif rows then carry no score arithmetic); `length` is the
/// number of positions either way.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct BindingMatrix {
    /// The matrix stable id (`ENSPFM0364`), the table's key.
    #[serde(skip)]
    pub stable_id: Arc<str>,
    pub name: Option<String>,
    pub source: Option<String>,
    pub threshold: Option<f64>,
    pub length: u64,
    #[serde(default)]
    pub frequencies: Vec<[f64; 4]>,
    /// The display names of the transcription factor complexes, in the
    /// cache's order; TRANSCRIPTION_FACTORS prints them.
    #[serde(default)]
    pub transcription_factors: Arc<[Arc<str>]>,
    /// The tables the scores read, derived from `frequencies` on first use.
    #[serde(skip)]
    scoring: OnceLock<Scoring>,
}

impl PartialEq for BindingMatrix {
    /// Two matrices are equal by their stored fields; the derived tables
    /// follow from the counts and are left out.
    fn eq(&self, other: &Self) -> bool {
        self.stable_id == other.stable_id
            && self.name == other.name
            && self.source == other.source
            && self.threshold == other.threshold
            && self.length == other.length
            && self.frequencies == other.frequencies
            && self.transcription_factors == other.transcription_factors
    }
}

impl BindingMatrix {
    /// A matrix from its parts, for a source that builds its matrices rather
    /// than reading the cache's table: no name, source or threshold, and
    /// `frequencies` empty for a matrix known by its length alone. Refuses
    /// counts whose row count disagrees with `length`.
    pub fn new(
        stable_id: &str,
        length: u64,
        frequencies: Vec<[f64; 4]>,
        transcription_factors: Arc<[Arc<str>]>,
    ) -> Result<Self, String> {
        Self::from_table_entry(
            stable_id,
            Self {
                stable_id: Arc::from(""),
                name: None,
                source: None,
                threshold: None,
                length,
                frequencies,
                transcription_factors,
                scoring: OnceLock::new(),
            },
        )
    }

    /// Reads one table entry under its stable id, refusing counts whose row
    /// count disagrees with `length`.
    pub fn from_table_entry(stable_id: &str, mut entry: BindingMatrix) -> Result<Self, String> {
        entry.stable_id = Arc::from(stable_id);
        if !entry.frequencies.is_empty() && entry.frequencies.len() as u64 != entry.length {
            return Err(format!(
                "binding matrix {stable_id} declares length {} but holds {} positions",
                entry.length,
                entry.frequencies.len()
            ));
        }
        Ok(entry)
    }

    /// Whether the matrix carries counts, the condition for the motif score
    /// arithmetic (`OutputFactory.pm:1935` writes HIGH_INF_POS and
    /// MOTIF_SCORE_CHANGE only for a matrix with elements).
    pub fn has_frequencies(&self) -> bool {
        !self.frequencies.is_empty()
    }

    /// The derived tables, `None` for a matrix without counts.
    fn scoring(&self) -> Option<&Scoring> {
        self.has_frequencies()
            .then(|| self.scoring.get_or_init(|| Scoring::of(&self.frequencies)))
    }

    /// The similarity of `sequence` to the matrix relative to the matrix's best
    /// and worst sites (`relative_sequence_similarity_score`,
    /// `BindingMatrix.pm:577-601`): the sum of the weights of the sequence's
    /// bases, placed between the sums of the smallest and of the largest weight
    /// of every position (`:491-509`), so the best site scores 1 and the worst
    /// 0. Case does not matter. `None` for a matrix without counts, a sequence
    /// of another length than the matrix or with a base other than A, C, G or
    /// T (funcgen throws on both, `:534-549`), and a matrix whose best and
    /// worst sites score alike.
    pub fn relative_similarity(&self, sequence: &[u8]) -> Option<f64> {
        let scoring = self.scoring()?;
        if sequence.len() != scoring.weights.len() {
            return None;
        }
        let range = scoring.max_score - scoring.min_score;
        if range == 0.0 {
            return None;
        }
        let mut score = 0.0;
        for (weights, base) in scoring.weights.iter().zip(sequence) {
            score += weights[nucleotide_index(*base)?];
        }
        Some((score - scoring.min_score) / range)
    }

    /// Whether `position` (1-based) of the matrix carries more than
    /// [`INFORMATIVE_BITS`] of information (`is_position_informative`,
    /// `BindingMatrix.pm:616-655`: the position's probabilities each times
    /// `2 - H`, `H` the position's Shannon entropy in bits, summed and compared
    /// strictly against the threshold). `false` for a matrix without counts
    /// and for a position outside the matrix (funcgen throws, `:624-629`).
    pub fn is_position_informative(&self, position: u32) -> bool {
        let Some(scoring) = self.scoring() else {
            return false;
        };
        position
            .checked_sub(1)
            .and_then(|index| scoring.bits.get(index as usize))
            .is_some_and(|bits| *bits > INFORMATIVE_BITS)
    }
}

/// The tables `Converter.pm` derives from a matrix's counts, held once per
/// matrix.
#[derive(Debug, Clone, PartialEq)]
struct Scoring {
    /// Each position's weights, A, C, G, T (`from_frequencies_to_weights`).
    weights: Vec<[f64; 4]>,
    /// Each position's information content in bits, the sum of its four
    /// `bits` elements (`from_frequencies_to_bits`).
    bits: Vec<f64>,
    /// The sum over positions of the smallest weight: the worst site's score.
    min_score: f64,
    /// The sum over positions of the largest weight: the best site's score.
    max_score: f64,
}

impl Scoring {
    fn of(frequencies: &[[f64; 4]]) -> Self {
        let mut weights = Vec::with_capacity(frequencies.len());
        let mut bits = Vec::with_capacity(frequencies.len());
        let (mut min_score, mut max_score) = (0.0, 0.0);
        for counts in frequencies {
            let probabilities = probabilities(counts);
            let position_weights = weights_of(&probabilities);
            min_score += position_weights
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min);
            max_score += position_weights
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            weights.push(position_weights);
            bits.push(bits_of(&probabilities).iter().sum());
        }
        Self {
            weights,
            bits,
            min_score,
            max_score,
        }
    }
}

/// A position's probabilities from its counts (`from_frequencies_to_probabilities`,
/// `Converter.pm:116-144`): each count plus the pseudocount, over the
/// position's total plus four pseudocounts.
fn probabilities(counts: &[f64; 4]) -> [f64; 4] {
    let total: f64 = counts.iter().sum();
    counts.map(|count| (count + PSEUDOCOUNT) / (total + 4.0 * PSEUDOCOUNT))
}

/// A position's weights from its probabilities (`from_probabilities_to_weights`,
/// `Converter.pm:215-239`): the log2 of each probability over the background
/// frequency, taken as funcgen takes it, a natural logarithm over `ln 2`.
fn weights_of(probabilities: &[f64; 4]) -> [f64; 4] {
    probabilities.map(|probability| log2(probability / EXPECTED_FREQUENCY))
}

/// A position's bits from its probabilities (`from_probabilities_to_bits`,
/// `Converter.pm:266-300`): each probability times the position's information
/// content, `2 - H` with `H` the Shannon entropy of the four probabilities.
fn bits_of(probabilities: &[f64; 4]) -> [f64; 4] {
    let mut entropy = 0.0;
    for probability in probabilities {
        entropy -= probability * log2(*probability);
    }
    let information = 2.0 - entropy;
    probabilities.map(|probability| probability * information)
}

/// `_log2` (`Converter.pm:426-431`): `log(n) / log(2)`.
fn log2(n: f64) -> f64 {
    n.ln() / std::f64::consts::LN_2
}

/// The column of `base` in a position's A, C, G, T row, either case; `None`
/// for any other byte.
fn nucleotide_index(base: u8) -> Option<usize> {
    NUCLEOTIDES
        .iter()
        .position(|nucleotide| *nucleotide == base.to_ascii_uppercase())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(json: &str) -> BindingMatrix {
        serde_json::from_str(json).unwrap()
    }

    /// The first position of ENSPFM0364 as the GRCh38 cache stores it.
    const ENSPFM0364_POSITION_1: [f64; 4] = [4881.0, 1573.0, 1789.0, 1847.0];

    fn close(actual: f64, expected: f64) -> bool {
        (actual - expected).abs() < 5e-7
    }

    #[test]
    fn table_entry_reads_its_fields_and_takes_the_key_as_id() {
        let matrix = BindingMatrix::from_table_entry(
            "ENSPFM0364",
            entry(
                r#"{"name":"MAFK","source":"SELEX","threshold":4.4,"length":2,
                    "frequencies":[[4881,1573,1789,1847],[1,2,3,4]],
                    "transcription_factors":["MAFK"]}"#,
            ),
        )
        .unwrap();
        assert_eq!(matrix.stable_id.as_ref(), "ENSPFM0364");
        assert_eq!(matrix.name.as_deref(), Some("MAFK"));
        assert_eq!(matrix.source.as_deref(), Some("SELEX"));
        assert_eq!(matrix.threshold, Some(4.4));
        assert_eq!(matrix.length, 2);
        assert_eq!(matrix.frequencies[0], ENSPFM0364_POSITION_1);
        assert_eq!(matrix.transcription_factors.len(), 1);
        assert!(matrix.has_frequencies());
        assert_eq!(NUCLEOTIDES, *b"ACGT");
        // Equality reads the stored fields: scoring one copy derives its tables
        // and leaves the two equal.
        let scored = matrix.clone();
        assert!(scored.relative_similarity(b"AC").is_some());
        assert_eq!(scored, matrix);
    }

    #[test]
    fn a_matrix_known_by_length_alone_has_no_frequencies() {
        let matrix = BindingMatrix::from_table_entry("STUB", entry(r#"{"length":6}"#)).unwrap();
        assert!(!matrix.has_frequencies());
        assert_eq!(matrix.length, 6);
        assert!(matrix.transcription_factors.is_empty());
        assert!(matrix.name.is_none() && matrix.source.is_none() && matrix.threshold.is_none());
        assert_eq!(matrix.relative_similarity(b"ACGTAC"), None);
        assert!(!matrix.is_position_informative(1));
    }

    #[test]
    fn a_row_count_disagreeing_with_the_length_is_refused() {
        let error = BindingMatrix::from_table_entry(
            "ENSPFM0001",
            entry(r#"{"length":3,"frequencies":[[1,2,3,4]]}"#),
        )
        .unwrap_err();
        assert_eq!(
            error,
            "binding matrix ENSPFM0001 declares length 3 but holds 1 positions"
        );
        assert!(
            serde_json::from_str::<BindingMatrix>(r#"{"length":1,"frequencies":[[1,2,3]]}"#)
                .is_err()
        );
    }

    /// Position 1 of ENSPFM0364 (counts 4881, 1573, 1789, 1847, total 10090):
    /// p(A) = 4881.1 / 10090.4 (`Converter.pm:142-144`), w(A) = log2(p(A) / 0.25)
    /// (`:238-239`), and the position's bits sum to `(2 - H) * sum(p)` with `H` its
    /// entropy (`:288-296`, `:433-453`); a uniform position has 0 bits and a
    /// position of one base close to 2.
    #[test]
    fn probabilities_weights_and_bits_follow_the_converter() {
        let p = probabilities(&ENSPFM0364_POSITION_1);
        assert!(close(p[0], 0.483_737), "{}", p[0]);
        assert!(close(p[0], 4881.1 / 10090.4));
        assert!(close(p[1], 1573.1 / 10090.4));
        assert!(close(p[2], 1789.1 / 10090.4));
        assert!(close(p[3], 1847.1 / 10090.4));
        assert!(close(p.iter().sum::<f64>(), 1.0));
        let w = weights_of(&p);
        assert!(close(w[0], 0.952_295), "{}", w[0]);
        assert!(close(w[0], (0.483_737_f64 / 0.25).log2()));
        assert!(w[1] < 0.0 && w[2] < 0.0 && w[3] < 0.0);
        let bits = bits_of(&p);
        let mut entropy = 0.0;
        for q in p {
            entropy -= q * q.log2();
        }
        assert!(close(bits.iter().sum::<f64>(), 2.0 - entropy));
        assert!(close(
            bits_of(&probabilities(&[7.0; 4])).iter().sum::<f64>(),
            0.0
        ));
        let single = bits_of(&probabilities(&[10_000.0, 0.0, 0.0, 0.0]));
        let sum: f64 = single.iter().sum();
        assert!(sum > 1.99 && sum < 2.0, "{sum}");
    }

    /// `relative_sequence_similarity_score` (`BindingMatrix.pm:594-597`) places
    /// a sequence's weight sum between the matrix's worst and best sites: the
    /// sequence of each position's most frequent base scores 1, of its least
    /// frequent 0, and a sequence of another length, with a base outside ACGT,
    /// or against a matrix without counts has no score; case does not matter.
    #[test]
    fn relative_similarity_of_the_best_site_is_one_and_of_the_worst_is_zero() {
        let matrix = BindingMatrix::from_table_entry(
            "ENSPFM0042",
            entry(
                r#"{"length":4,"frequencies":[[5461,545,1745,4231],[2613,450,11439,632],[0,8101,0,668],[64,0,12567,86]]}"#,
            ),
        )
        .unwrap();
        assert!(close(matrix.relative_similarity(b"AGCG").unwrap(), 1.0));
        assert!(close(matrix.relative_similarity(b"agcg").unwrap(), 1.0));
        assert!(close(matrix.relative_similarity(b"CCAC").unwrap(), 0.0));
        let middle = matrix.relative_similarity(b"AGTG").unwrap();
        assert!(middle > 0.0 && middle < 1.0, "{middle}");
        // The weight sums, as funcgen sums them position by position.
        let scoring = matrix.scoring().unwrap();
        let expected_max: f64 = scoring
            .weights
            .iter()
            .map(|w| w.iter().copied().fold(f64::NEG_INFINITY, f64::max))
            .sum();
        assert_eq!(scoring.max_score, expected_max);
        assert_eq!(matrix.relative_similarity(b"AGC"), None);
        assert_eq!(matrix.relative_similarity(b"AGCGA"), None);
        assert_eq!(matrix.relative_similarity(b"AGNG"), None);
        assert_eq!(matrix.relative_similarity(b"AG-G"), None);
        let flat = BindingMatrix::from_table_entry(
            "FLAT",
            entry(r#"{"length":2,"frequencies":[[5,5,5,5],[5,5,5,5]]}"#),
        )
        .unwrap();
        assert_eq!(flat.relative_similarity(b"AC"), None);
    }

    /// `is_position_informative` (`BindingMatrix.pm:616-655`) compares the
    /// position's bits strictly against 1.5 (`:650`): ENSPFM0042's position 3
    /// (1.611 bits) is informative and its position 2 (0.915) is not; a position
    /// whose bits sum to exactly 1.5 is not; a position outside the matrix is not
    /// (funcgen throws there, `:624-629`).
    #[test]
    fn position_informative_strictly_over_one_point_five() {
        let matrix = BindingMatrix::from_table_entry(
            "ENSPFM0042",
            entry(
                r#"{"length":3,"frequencies":[[5461,545,1745,4231],[2613,450,11439,632],[0,8101,0,668]]}"#,
            ),
        )
        .unwrap();
        let bits = &matrix.scoring().unwrap().bits;
        assert!((bits[2] - 1.611).abs() < 5e-4, "{}", bits[2]);
        assert!((bits[1] - 0.915).abs() < 5e-4, "{}", bits[1]);
        assert!(!matrix.is_position_informative(1));
        assert!(!matrix.is_position_informative(2));
        assert!(matrix.is_position_informative(3));
        assert!(!matrix.is_position_informative(0));
        assert!(!matrix.is_position_informative(4));
        let mut exact = matrix.clone();
        exact.scoring = OnceLock::new();
        exact
            .scoring
            .set(Scoring {
                weights: vec![[0.0; 4]; 3],
                bits: vec![INFORMATIVE_BITS, INFORMATIVE_BITS + f64::EPSILON, 2.0],
                min_score: 0.0,
                max_score: 1.0,
            })
            .unwrap();
        assert!(!exact.is_position_informative(1));
        assert!(exact.is_position_informative(2));
        assert!(exact.is_position_informative(3));
    }
}
