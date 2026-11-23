// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The made-up values. Each is drawn from the stream that the seed, the name of the
//! value's table and the row's key name, so it depends on nothing else: a value does not
//! move when the cache gains or loses a transcript, and the three tables draw from
//! streams of their own. A value has at most three decimals, so Ensembl VEP, which prints
//! a number with up to fifteen significant digits, and vep-rs print it as the same text.

use vep_seeded_variants::random::{keyed_stream, SplitMix64};

/// The stream name of the table of base scores.
const BASE_SCORES: &[u8] = b"GERP_scores";

/// The stream name of the `gerp_exons` table.
const EXON_TOTALS: &[u8] = b"gerp_exons";

/// The stream name of the `phylocsf_data` table.
const PHYLOCSF_SCORES: &[u8] = b"phylocsf_data";

/// The thousandths a score is a multiple of: three, two, one or no decimals.
const STEPS: [i64; 4] = [1, 10, 100, 1000];

/// The lowest and the highest score of a base, in thousandths.
pub const BASE_SCORE_RANGE: (i64, i64) = (-6_000, 9_000);

/// The lowest and the highest total of an exon, in hundredths.
pub const EXON_TOTAL_RANGE: (i64, i64) = (-30_000, 90_000);

/// The largest annotated PhyloCSF score that is not negative, in thousandths.
pub const MOST_ANNOTATED: i64 = 400_000;

/// The size of the lowest negative annotated PhyloCSF score, in thousandths.
pub const MOST_NEGATIVE: i64 = 300_000;

/// The largest highest PhyloCSF score, in thousandths.
pub const MOST_HIGHEST: i64 = 400_000;

/// One annotated score in this many that are not negative is an exact 0.
const ANNOTATED_ZERO_ONE_IN: u64 = 8;

/// One highest score in this many that are not positive is an exact 0.
const HIGHEST_ZERO_ONE_IN: u64 = 3;

/// The decimal digits of `number`, written at the end of `buffer`.
pub fn decimal(mut number: u64, buffer: &mut [u8; 20]) -> &[u8] {
    let mut at = buffer.len();
    loop {
        at -= 1;
        buffer[at] = b'0' + (number % 10) as u8;
        number /= 10;
        if number == 0 {
            break;
        }
    }
    &buffer[at..]
}

/// Appends the shortest decimal of `thousandths` to `text`: no exponent, no trailing
/// zero, and `0` for zero.
pub fn push_thousandths(text: &mut Vec<u8>, thousandths: i64) {
    if thousandths < 0 {
        text.push(b'-');
    }
    let size = thousandths.unsigned_abs();
    let mut buffer = [0; 20];
    text.extend_from_slice(decimal(size / 1000, &mut buffer));
    let fraction = size % 1000;
    if fraction == 0 {
        return;
    }
    let digits = [
        b'0' + (fraction / 100) as u8,
        b'0' + (fraction / 10 % 10) as u8,
        b'0' + (fraction % 10) as u8,
    ];
    let kept = match digits {
        [_, b'0', b'0'] => 1,
        [_, _, b'0'] => 2,
        _ => 3,
    };
    text.push(b'.');
    text.extend_from_slice(&digits[..kept]);
}

/// The shortest decimal of `thousandths`.
pub fn thousandths_text(thousandths: i64) -> String {
    let mut text = Vec::new();
    push_thousandths(&mut text, thousandths);
    String::from_utf8_lossy(&text).into_owned()
}

/// The double nearest to `thousandths` as a number: both operands of the division are
/// exact, so its result is the double the decimal text reads as.
pub fn thousandths_real(thousandths: i64) -> f64 {
    thousandths as f64 / 1000.0
}

/// The double nearest to `hundredths` as a number.
pub fn hundredths_real(hundredths: i64) -> f64 {
    hundredths as f64 / 100.0
}

fn step(stream: &mut SplitMix64) -> i64 {
    STEPS[stream.below(STEPS.len() as u64) as usize]
}

/// A multiple of `step` from `step` to `most`, every multiple equally likely.
fn multiple(stream: &mut SplitMix64, step: i64, most: i64) -> i64 {
    step * stream.between(1, (most / step) as u64) as i64
}

/// The score of the base at `position` of `sequence`, in thousandths: a step is drawn
/// first, then a multiple of it inside [`BASE_SCORE_RANGE`], both ends included. The
/// scores average 1.5, so the bases of a run of about 120 add up to the 180 that
/// LoFTEE's `LoF.pm:268` holds the walk's total against.
pub fn base_score(seed: u64, sequence: &str, position: u64) -> i64 {
    let mut buffer = [0; 20];
    base_score_at(seed, sequence, decimal(position, &mut buffer))
}

/// [`base_score`] for a position given as its decimal digits, the form the stream's key
/// holds it in.
pub fn base_score_at(seed: u64, sequence: &str, digits: &[u8]) -> i64 {
    let mut stream = keyed_stream(seed, &[BASE_SCORES, sequence.as_bytes(), digits]);
    let step = step(&mut stream);
    let (lowest, highest) = BASE_SCORE_RANGE;
    lowest + step * stream.between(0, ((highest - lowest) / step) as u64) as i64
}

/// The total of exon `exon` of `transcript`, in hundredths inside [`EXON_TOTAL_RANGE`],
/// both ends included. The key is the stable id and the exon's number: no base score
/// enters the draw, so a total read from the table's rows differs from this one.
pub fn exon_total(seed: u64, transcript: &str, exon: u64) -> i64 {
    let mut buffer = [0; 20];
    let number = decimal(exon, &mut buffer);
    let mut stream = keyed_stream(seed, &[EXON_TOTALS, transcript.as_bytes(), number]);
    let (lowest, highest) = EXON_TOTAL_RANGE;
    lowest + stream.between(0, (highest - lowest) as u64) as i64
}

/// The two PhyloCSF scores of an exon, in thousandths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Phylocsf {
    /// The score of the annotated reading frame, `corresponding_orf_score`.
    pub annotated: i64,
    /// The highest score of the exon, `max_score`.
    pub highest: i64,
}

/// What LoFTEE's `LoF.pm:295-296` makes of an exon's two scores.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Flag {
    /// The annotated score is not negative: no flag.
    None,
    /// The annotated score is negative and the highest is above 0:
    /// `PHYLOCSF_UNLIKELY_ORF`.
    UnlikelyOrf,
    /// The annotated score is negative and the highest is not above 0: `PHYLOCSF_WEAK`.
    Weak,
}

impl Phylocsf {
    /// The flag the two scores give.
    pub fn flag(&self) -> Flag {
        match (self.annotated < 0, self.highest > 0) {
            (false, _) => Flag::None,
            (true, true) => Flag::UnlikelyOrf,
            (true, false) => Flag::Weak,
        }
    }
}

/// The PhyloCSF scores of exon `exon` of `transcript`. One of the three flags is drawn
/// first, each as likely as the others, then the step both scores are multiples of, then
/// the scores: an annotated score that is not negative lies within [`MOST_ANNOTATED`], a
/// negative one within [`MOST_NEGATIVE`] and a positive highest score within
/// [`MOST_HIGHEST`]. The highest score is never below the annotated one, and exact
/// zeros stand among the annotated scores and among the highest ones.
pub fn phylocsf_scores(seed: u64, transcript: &str, exon: u64) -> Phylocsf {
    let mut buffer = [0; 20];
    let number = decimal(exon, &mut buffer);
    let mut stream = keyed_stream(seed, &[PHYLOCSF_SCORES, transcript.as_bytes(), number]);
    let flag = [Flag::None, Flag::UnlikelyOrf, Flag::Weak][stream.below(3) as usize];
    let step = step(&mut stream);
    match flag {
        Flag::None => {
            let annotated = if stream.below(ANNOTATED_ZERO_ONE_IN) == 0 {
                0
            } else {
                multiple(&mut stream, step, MOST_ANNOTATED)
            };
            let above = ((MOST_HIGHEST - annotated) / step) as u64;
            Phylocsf {
                annotated,
                highest: annotated + step * stream.between(0, above) as i64,
            }
        }
        Flag::UnlikelyOrf => Phylocsf {
            annotated: -multiple(&mut stream, step, MOST_NEGATIVE),
            highest: multiple(&mut stream, step, MOST_HIGHEST),
        },
        Flag::Weak => {
            let annotated = -multiple(&mut stream, step, MOST_NEGATIVE);
            let highest = if stream.below(HIGHEST_ZERO_ONE_IN) == 0 {
                0
            } else {
                -multiple(&mut stream, step, -annotated)
            };
            Phylocsf { annotated, highest }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_score_reads_as_its_shortest_decimal() {
        for (thousandths, text) in [
            (0, "0"),
            (-500, "-0.5"),
            (3_000, "3"),
            (2_125, "2.125"),
            (-6_000, "-6"),
            (9_000, "9"),
            (10, "0.01"),
            (120, "0.12"),
            (-1, "-0.001"),
            (-299_990, "-299.99"),
            (400_000, "400"),
            (1_050, "1.05"),
        ] {
            assert_eq!(thousandths_text(thousandths), text);
        }
        let mut buffer = [0; 20];
        assert_eq!(decimal(0, &mut buffer), b"0");
        assert_eq!(decimal(536_870_912, &mut buffer), b"536870912");
        assert_eq!(decimal(u64::MAX, &mut buffer), b"18446744073709551615");
    }

    /// The digits of `text` from its first digit above zero to its last.
    fn significant(text: &str) -> String {
        let digits: String = text.chars().filter(char::is_ascii_digit).collect();
        digits.trim_matches('0').to_owned()
    }

    /// The fifteen significant digits of `real`, without the zeros that close them.
    fn fifteen_digits(real: f64) -> String {
        let rounded = format!("{:.14e}", real.abs());
        let (digits, _) = rounded.split_once('e').unwrap();
        significant(digits)
    }

    /// The double a value is stored as prints as the value's own decimal, in Rust's
    /// shortest form and in the fifteen significant digits Ensembl VEP prints a number
    /// with.
    #[test]
    fn a_stored_value_prints_as_its_decimal() {
        for thousandths in (-400_000..=400_000).step_by(7) {
            let text = thousandths_text(thousandths);
            let real = thousandths_real(thousandths);
            assert_eq!(real, text.parse::<f64>().unwrap(), "{text}");
            assert_eq!(format!("{real}"), text);
            if thousandths != 0 {
                assert_eq!(fifteen_digits(real), significant(&text), "{text}");
            }
        }
        for hundredths in (EXON_TOTAL_RANGE.0..=EXON_TOTAL_RANGE.1).step_by(3) {
            let real = hundredths_real(hundredths);
            assert_eq!(format!("{real}"), thousandths_text(hundredths * 10));
        }
        assert_eq!(format!("{}", thousandths_real(0)), "0");
        assert!(thousandths_real(0).is_sign_positive());
    }

    #[test]
    fn a_base_score_is_a_multiple_of_its_step_inside_the_range() {
        let (lowest, highest) = BASE_SCORE_RANGE;
        let scores: Vec<i64> = (1..=40_000)
            .map(|position| base_score(1, "made_up", position))
            .collect();
        assert!(scores
            .iter()
            .all(|score| (lowest..=highest).contains(score)));
        assert!(scores.contains(&lowest) && scores.contains(&highest));
        assert!(scores.contains(&0));
        assert!(scores.iter().any(|score| *score < 0));
        assert!(scores.iter().any(|score| *score > 0));
        for step in STEPS {
            assert!(
                scores
                    .iter()
                    .any(|score| score % step == 0 && score % (step * 10) != 0),
                "no score ends on the step of {step} thousandths"
            );
        }
        let mean = scores.iter().sum::<i64>() / scores.len() as i64;
        assert!((1_400..=1_600).contains(&mean), "the mean is {mean}");
    }

    /// Runs of a hundred bases add up to totals on both sides of the 180 of
    /// `LoF.pm:268`.
    #[test]
    fn a_run_of_a_hundred_bases_crosses_180() {
        let totals: Vec<i64> = (0..200)
            .map(|run| {
                (1..=100)
                    .map(|base| base_score(1, "made_up", run * 100 + base))
                    .sum()
            })
            .collect();
        assert!(totals.iter().any(|total| *total > 180_000));
        assert!(totals.iter().any(|total| *total <= 180_000));
    }

    #[test]
    fn an_exon_total_stays_inside_its_range_on_both_sides_of_180() {
        let (lowest, highest) = EXON_TOTAL_RANGE;
        let totals: Vec<i64> = (1..=20_000)
            .map(|exon| exon_total(1, "MADE_t", exon))
            .collect();
        assert!(totals
            .iter()
            .all(|total| (lowest..=highest).contains(total)));
        assert!(totals.iter().any(|total| *total > 18_000));
        assert!(totals.iter().any(|total| *total <= 18_000));
        assert!(totals.iter().any(|total| *total < 0));
        assert!(totals.iter().any(|total| total % 100 != 0));
    }

    #[test]
    fn the_three_phylocsf_flags_are_drawn_evenly_inside_their_bounds() {
        let rows: Vec<Phylocsf> = (1..=30_000)
            .map(|exon| phylocsf_scores(1, "MADE_t", exon))
            .collect();
        for flag in [Flag::None, Flag::UnlikelyOrf, Flag::Weak] {
            let share = rows.iter().filter(|row| row.flag() == flag).count();
            assert!((9_500..=10_500).contains(&share), "{flag:?}: {share}");
        }
        for row in &rows {
            assert!(row.highest >= row.annotated, "{row:?}");
            assert!((-MOST_NEGATIVE..=MOST_ANNOTATED).contains(&row.annotated));
            assert!((-MOST_NEGATIVE..=MOST_HIGHEST).contains(&row.highest));
        }
        let zero_annotated = rows.iter().filter(|row| row.annotated == 0).count();
        assert!(
            (1_000..=1_500).contains(&zero_annotated),
            "{zero_annotated}"
        );
        let weak_zero = rows
            .iter()
            .filter(|row| row.flag() == Flag::Weak && row.highest == 0)
            .count();
        assert!((3_000..=3_700).contains(&weak_zero), "{weak_zero}");
        assert!(rows.iter().any(|row| row.highest < 0));
        assert!(rows.iter().any(|row| row.annotated % 10 != 0));
        assert!(rows
            .iter()
            .any(|row| row.annotated % 1_000 == 0 && row.annotated != 0));
    }

    #[test]
    fn a_value_is_a_function_of_the_seed_and_its_key() {
        assert_eq!(base_score(7, "1", 1_234), base_score(7, "1", 1_234));
        assert_eq!(exon_total(7, "MADE_t", 3), exon_total(7, "MADE_t", 3));
        assert_eq!(
            phylocsf_scores(7, "MADE_t", 3),
            phylocsf_scores(7, "MADE_t", 3)
        );
        let moved = |other: &dyn Fn(u64) -> bool| (1..=200).filter(|key| other(*key)).count();
        assert!(moved(&|key| base_score(7, "1", key) != base_score(8, "1", key)) > 150);
        assert!(moved(&|key| base_score(7, "1", key) != base_score(7, "2", key)) > 150);
        assert!(moved(&|key| exon_total(7, "MADE_t", key) != exon_total(8, "MADE_t", key)) > 150);
        assert!(moved(&|key| exon_total(7, "MADE_t", key) != exon_total(7, "MADE_u", key)) > 150);
        assert!(
            moved(&|key| phylocsf_scores(7, "MADE_t", key) != phylocsf_scores(8, "MADE_t", key))
                > 150
        );
    }
}
