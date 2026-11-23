// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Locations, spans and breakend alleles as the rules read them.

use super::py_int;
use regex::Regex;
use std::cmp::Ordering;
use std::fmt;
use std::sync::LazyLock;

/// Ensembl VEP's `--max_sv_size` default (ensembl-vep Config.pm:310); a span above it is
/// marked skipped (Parser.pm:493-497).
const MAX_SV_SIZE: i128 = 10_000_000;

/// The symbolic ALT types without a Sequence Ontology term, marked skipped at
/// Parser/VCF.pm:477-481; the Allele column carries the type string (Parser/VCF.pm:469-475).
const VEP_SKIP_UNSUPPORTED_ALLELES: [&str; 4] = ["CPX", "CTX", "<CPX>", "<CTX>"];

/// The smallest coordinate at which Ensembl VEP 116.2 annotates every record of the input
/// buffer as intergenic_variant: the coordinate plus the default `--distance` flank reaches
/// 2^31.
pub const BND_MATE_32_BIT_LIMIT: i128 = (1 << 31) - 5000;

/// A coordinate with fewer digits than the bound cannot reach it.
pub const BND_MATE_32_BIT_LIMIT_DIGITS: usize = 10;

/// The largest coordinate vep-rs holds; a POS past it fails its line in every input format.
const VEP_RS_COORDINATE_CEILING: i128 = (1 << 62) - 1;

/// The class term for a breakend's own end; also the allele of a bare breakpoint record.
pub const BARE_BREAKPOINT_ALLELE: &str = "chromosome_breakpoint";

/// The one symbolic record both engines keep but write at different Locations.
pub const TANDEM_REPEAT_ALLELE: &str = "tandem_repeat";

pub const INTERGENIC_FEATURE: &str = "-";

static BRACKET_MATE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]").expect("regex"));
static OWN_BREAKEND: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^(?:\.[ACGTN]+|[ACGTN]+\.|chromosome_breakpoint)$").expect("regex")
});
static LOCATION_SPAN: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^([^:]+):(-?[0-9]+)(?:-(-?[0-9]+))?$").expect("regex"));
static INTEGER_COORDS: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^-?[0-9]+(?:--?[0-9]+)?$").expect("regex"));
static VEP_RS_POS_NUMBER: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^\+?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$").expect("regex")
});

/// A coordinate of any length, read as Python's `int` reads it: a value within `i128` as
/// that word, a longer one as its sign and its digits without leading zeros, so one value has
/// one form.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Coordinate {
    Word(i128),
    Long { negative: bool, digits: String },
}

impl Coordinate {
    /// Reads an optional `-` before the digits `0` to `9`; None for any other text.
    fn read(text: &str) -> Option<Coordinate> {
        let digits = text.strip_prefix('-').unwrap_or(text);
        if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
            return None;
        }
        Some(match text.parse() {
            Ok(value) => Coordinate::Word(value),
            Err(_) => Coordinate::Long {
                negative: digits.len() < text.len(),
                digits: digits.trim_start_matches('0').to_string(),
            },
        })
    }

    /// Whether the coordinate is `bound` or past it.
    pub fn at_least(&self, bound: i128) -> bool {
        match self {
            Coordinate::Word(value) => *value >= bound,
            Coordinate::Long { negative, .. } => !negative,
        }
    }

    pub fn is_negative(&self) -> bool {
        match self {
            Coordinate::Word(value) => *value < 0,
            Coordinate::Long { negative, .. } => *negative,
        }
    }

    /// The coordinate after this one.
    pub fn next(&self) -> Coordinate {
        match self {
            Coordinate::Word(value) => match value.checked_add(1) {
                Some(next) => Coordinate::Word(next),
                None => Coordinate::Long {
                    negative: false,
                    digits: step_digits(&value.to_string(), true),
                },
            },
            Coordinate::Long { negative, digits } => {
                let sign = if *negative { "-" } else { "" };
                let stepped = step_digits(digits, !negative);
                Coordinate::read(&format!("{sign}{stepped}")).expect("an integer")
            }
        }
    }
}

/// The digits of a positive integer one above (`up`) or one below the integer `digits` spell.
fn step_digits(digits: &str, up: bool) -> String {
    let (turned, carried) = if up { (b'9', b'0') } else { (b'0', b'9') };
    let mut stepped = digits.as_bytes().to_vec();
    let mut at = stepped.len();
    while at > 0 && stepped[at - 1] == turned {
        stepped[at - 1] = carried;
        at -= 1;
    }
    if at == 0 {
        stepped.insert(0, b'1');
    } else if up {
        stepped[at - 1] += 1;
    } else {
        stepped[at - 1] -= 1;
    }
    String::from_utf8(stepped).expect("ASCII digits")
}

impl Ord for Coordinate {
    fn cmp(&self, other: &Coordinate) -> Ordering {
        use Coordinate::{Long, Word};
        match (self, other) {
            (Word(a), Word(b)) => a.cmp(b),
            (Long { negative, .. }, Word(_)) if *negative => Ordering::Less,
            (Long { .. }, Word(_)) => Ordering::Greater,
            (Word(_), Long { .. }) => other.cmp(self).reverse(),
            (
                Long {
                    negative: a_negative,
                    digits: a,
                },
                Long {
                    negative: b_negative,
                    digits: b,
                },
            ) => {
                let magnitude = a.len().cmp(&b.len()).then_with(|| a.cmp(b));
                match (a_negative, b_negative) {
                    (false, false) => magnitude,
                    (true, true) => magnitude.reverse(),
                    (false, true) => Ordering::Greater,
                    (true, false) => Ordering::Less,
                }
            }
        }
    }
}

impl PartialOrd for Coordinate {
    fn partial_cmp(&self, other: &Coordinate) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Coordinate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Coordinate::Word(value) => write!(f, "{value}"),
            Coordinate::Long { negative, digits } => {
                write!(f, "{}{digits}", if *negative { "-" } else { "" })
            }
        }
    }
}

/// `chrom:start-end` or `chrom:pos` as (chrom, start, end), a coordinate the reference
/// numified below zero included and a coordinate of any length read whole; None when
/// unparseable.
pub fn location_span(loc: &str) -> Option<(&str, Coordinate, Coordinate)> {
    let m = LOCATION_SPAN.captures(loc)?;
    let start = Coordinate::read(&m[2])?;
    let end = match m.get(3) {
        Some(end) => Coordinate::read(end.as_str())?,
        None => start.clone(),
    };
    Some((m.get(1).expect("chrom group").as_str(), start, end))
}

/// Whether the Location's span exceeds Ensembl VEP's `--max_sv_size` default; a point
/// Location never does, and the magnitude is taken since an insertion's end precedes its
/// start.
fn span_exceeds_max_sv_size(loc: &str) -> bool {
    let Some((_, coords)) = loc.split_once(':') else {
        return false;
    };
    let Some((start, end)) = coords.split_once('-') else {
        return false;
    };
    match (py_int(end), py_int(start)) {
        // The width of two 128-bit coordinates can pass 128 bits.
        (Some(end), Some(start)) => end.abs_diff(start) > MAX_SV_SIZE.unsigned_abs(),
        _ => false,
    }
}

/// Whether Ensembl VEP annotates the record only against its batch's loaded regions: a span
/// above `--max_sv_size` (Parser.pm:493-497) or a symbolic type without a Sequence Ontology
/// term (Parser/VCF.pm:477-481), both through the `vep_skip` flag.
pub fn batch_dependent_annotation_scope(loc: &str, allele: &str) -> bool {
    span_exceeds_max_sv_size(loc) || VEP_SKIP_UNSUPPORTED_ALLELES.contains(&allele)
}

/// The mate (chromosome, position) a bracket-notation breakend allele names, the chromosome
/// with any `chr` prefix removed as the reference removes it
/// (AnnotationType/Transcript.pm:134), the position read whole at any length; None for any
/// other allele, a position in digits other than `0` to `9` among them.
pub fn bracket_mate(allele: &str) -> Option<(&str, Coordinate)> {
    let m = BRACKET_MATE.captures(allele)?;
    let chrom = m.get(1).expect("chrom group").as_str();
    let position = Coordinate::read(&m[2])?;
    Some((strip_chr(chrom), position))
}

/// Whether the allele column names a breakend record's own end rather than a mate.
pub fn is_own_breakend_allele(allele: &str) -> bool {
    OWN_BREAKEND.is_match(allele)
}

/// Whether the allele column names a breakend's own end or its mate.
pub fn is_breakend_allele(allele: &str) -> bool {
    is_own_breakend_allele(allele) || bracket_mate(allele).is_some()
}

/// The Locations at which vep-rs can write the point of a dot-notation breakend without
/// SVTYPE whose reference row spans (chromosome, low, high): the span's start and the base
/// after it, the reference taking the sequence rules for the coordinates
/// (BaseVCF4.pm:265-274, 297-317).
pub fn dot_breakend_points(span: &(&str, Coordinate, Coordinate)) -> [String; 2] {
    let (chrom, start, _end) = span;
    [
        format!("{chrom}:{start}"),
        format!("{chrom}:{}", start.next()),
    ]
}

/// The mate chromosome of a bracket-notation breakend allele as the reference resolves it:
/// a lower-case `chr` prefix removed and any other spelling kept, so `N[chr21:...[` names
/// `21` and `N[Chr21:...[` names `Chr21` (StructuralVariationOverlap.pm:86-103).
pub fn reference_mate_chromosome(allele: &str) -> Option<&str> {
    let m = BRACKET_MATE.captures(allele)?;
    let chrom = m.get(1).expect("chrom group").as_str();
    Some(chrom.strip_prefix("chr").unwrap_or(chrom))
}

/// The reference's normalisation of a breakend's chromosome name at the slice test
/// (AnnotationType/Transcript.pm:134): a `chr` prefix of any case removed, nothing else.
pub fn strip_chr(chrom: &str) -> &str {
    if chrom.len() >= 3 && chrom.is_char_boundary(3) && chrom[..3].eq_ignore_ascii_case("chr") {
        &chrom[3..]
    } else {
        chrom
    }
}

/// Whether a coordinate of the Location is not an integer (`21:30425651-30425949.5`,
/// `21:30425949.0`); a coordinate numified below zero is an integer.
pub fn non_integer_end(location: &str) -> bool {
    match location.split_once(':') {
        Some((_, coords)) => !INTEGER_COORDS.is_match(coords),
        None => false,
    }
}

/// The POS column as vep-rs's parsers read it: an unpadded decimal integer up to the
/// coordinate ceiling exactly, or another spelling of an exact non-negative integer below
/// the ceiling that Rust's `f64` parser accepts, blanks around it ignored; None where vep-rs
/// fails the line.
pub fn vep_rs_pos(text: &str) -> Option<i128> {
    let stripped = text.trim_matches([' ', '\t', '\n', '\r', '\x0c']);
    if !VEP_RS_POS_NUMBER.is_match(stripped) {
        return None;
    }
    let digits = text.strip_prefix('+').unwrap_or(text);
    if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
        let value = py_int(text)?;
        return (value <= VEP_RS_COORDINATE_CEILING).then_some(value);
    }
    let number: f64 = stripped.parse().ok()?;
    (number.is_finite() && number.fract() == 0.0 && number < 2f64.powi(62))
        .then_some(number as i128)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn coordinate(text: &str) -> Coordinate {
        Coordinate::read(text).expect("an integer")
    }

    #[test]
    fn a_coordinate_longer_than_a_word_keeps_its_value() {
        let long = "9".repeat(46);
        let location = format!("21:5000-{long}");
        let (chrom, start, end) = location_span(&location).expect("a span");
        assert_eq!((chrom, start), ("21", Coordinate::Word(5000)));
        assert_eq!(end.to_string(), long);
        assert!(end.at_least(BND_MATE_32_BIT_LIMIT));
        assert!(!coordinate(&format!("-{long}")).at_least(BND_MATE_32_BIT_LIMIT));
        assert!(coordinate(&format!("-{long}")).is_negative());
        let allele = format!("N[chr5:{long}[");
        let (chrom, position) = bracket_mate(&allele).expect("a mate");
        assert_eq!((chrom, position), ("5", end));
    }

    #[test]
    fn one_value_has_one_form_whatever_its_leading_zeros() {
        let word = i128::MAX.to_string();
        assert_eq!(
            coordinate(&format!("000{word}")),
            Coordinate::Word(i128::MAX)
        );
        assert_eq!(coordinate("-0"), Coordinate::Word(0));
        let long = format!("{word}0");
        assert_eq!(coordinate(&format!("00{long}")), coordinate(&long));
        assert_eq!(coordinate(&format!("00{long}")).to_string(), long);
    }

    #[test]
    fn coordinates_order_by_value_across_the_word() {
        let long = format!("{}0", i128::MAX);
        let ascending = [
            format!("-{long}9"),
            format!("-{long}1"),
            format!("-{long}"),
            i128::MIN.to_string(),
            "-1".to_string(),
            "0".to_string(),
            i128::MAX.to_string(),
            long.clone(),
            format!("{long}1"),
            format!("{long}9"),
        ]
        .map(|text| coordinate(&text));
        for (at, low) in ascending.iter().enumerate() {
            assert_eq!(low.cmp(low), Ordering::Equal, "{low}");
            for high in &ascending[at + 1..] {
                assert!(low < high && high > low, "{low} against {high}");
            }
        }
    }

    #[test]
    fn the_coordinate_after_one_crosses_the_word_and_carries() {
        let past_the_word = "170141183460469231731687303715884105728";
        assert_eq!(Coordinate::Word(41).next(), Coordinate::Word(42));
        assert_eq!(
            Coordinate::Word(i128::MAX).next().to_string(),
            past_the_word
        );
        assert_eq!(
            coordinate(&"9".repeat(40)).next().to_string(),
            format!("1{}", "0".repeat(40))
        );
        let below_the_word = coordinate(&format!("-{past_the_word}").replace("728", "729"));
        assert_eq!(below_the_word.next(), Coordinate::Word(i128::MIN));
        assert_eq!(
            coordinate(&format!("-1{}", "0".repeat(40)))
                .next()
                .to_string(),
            format!("-{}", "9".repeat(40))
        );
        let location = format!("21:{}-{}", "9".repeat(40), "9".repeat(41));
        let span = location_span(&location).expect("a span");
        assert_eq!(
            dot_breakend_points(&span),
            [
                format!("21:{}", "9".repeat(40)),
                format!("21:1{}", "0".repeat(40))
            ]
        );
    }

    #[test]
    fn digits_outside_zero_to_nine_are_no_coordinate() {
        assert_eq!(Coordinate::read("١٢"), None);
        assert_eq!(Coordinate::read("-"), None);
        assert_eq!(Coordinate::read("+5"), None);
        assert_eq!(Coordinate::read(""), None);
        assert_eq!(bracket_mate("N[5:١٢["), None);
        assert_eq!(reference_mate_chromosome("N[Chr21:١٢["), None);
        assert!(!is_breakend_allele("N[5:١٢["));
    }

    #[test]
    fn a_span_whose_width_passes_128_bits_exceeds_the_size_limit() {
        let bound = i128::MAX;
        // The width is 2^128 - 2, then 2^127 + 1, then 2^127 exactly.
        for location in [
            format!("21:{bound}--{bound}"),
            format!("21:2--{bound}"),
            format!("21:1--{bound}"),
        ] {
            assert!(span_exceeds_max_sv_size(&location), "{location}");
            assert!(
                batch_dependent_annotation_scope(&location, "deletion"),
                "{location}"
            );
        }
        let narrow = format!("21:{bound}-{}", bound - MAX_SV_SIZE);
        assert!(!span_exceeds_max_sv_size(&narrow));
        assert!(span_exceeds_max_sv_size(&format!(
            "21:{bound}-{}",
            bound - MAX_SV_SIZE - 1
        )));
    }
}

#[cfg(test)]
#[path = "breakend_ten_digits_tests.rs"]
mod ten_digits_tests;

#[cfg(test)]
#[path = "breakend_coordinate_readers_tests.rs"]
mod coordinate_readers_tests;
