// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The span Ensembl VEP gives a symbolic structural-variant record from its INFO
//! fields, shared by every VCF reader.
//!
//! A symbolic record starts at POS + 1 (ensembl-io `BaseVCF4.pm` `get_start`, lines
//! 255-264: the base before the event is the anchor) and ends where `get_end` (lines
//! 297-317) says: at POS + |SVLEN| when SVLEN is defined, else at INFO/END, else at the
//! start. A breakpoint record never reads END for its end (`Parser/VCF.pm` lines
//! 487-489 delete it), a tandem repeat or tandem duplication reads the largest SVLEN
//! (lines 502-512), and `validate_vf` (`Parser.pm` lines 475-488) drops a record whose
//! end is not a number or lies more than one base before its start.

use std::fmt;

use crate::coordinate::{parse_coordinate, MAX_COORD};
use crate::variant::VariantClass;

/// Perl's numeric reading of an INFO value, as `abs($svlen)` and `max(@svlen)` apply
/// it: leading whitespace, an optional sign, digits with an optional fraction and
/// exponent, and nothing else counts, so `300abc` is 300, ` 300` and `3e2` are 300 and
/// `abc`, `.` or an empty value is 0. The integer part is returned, saturating at the
/// machine word (`1e400` is `Inf` to Perl); a fractional SVLEN ends Ensembl VEP's span
/// at a fraction no coordinate can carry.
pub(crate) fn perl_numeric_value(text: &str) -> i64 {
    perl_numeric_f64(text).trunc() as i64
}

/// Perl's numeric reading of a value before any truncation, the number `x` and a range
/// read before they truncate it (`2abc` is 2.0, `1.9` is 1.9, `abc` is 0.0).
pub(crate) fn perl_numeric_f64(text: &str) -> f64 {
    let rest = text.trim_start();
    let bytes = rest.as_bytes();
    let mut i = 0;
    if matches!(bytes.first(), Some(b'+' | b'-')) {
        i += 1;
    }
    let digits_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    let mut mantissa_end = i;
    if i < bytes.len() && bytes[i] == b'.' {
        let mut j = i + 1;
        while j < bytes.len() && bytes[j].is_ascii_digit() {
            j += 1;
        }
        if j > i + 1 || i > digits_start {
            mantissa_end = j;
            i = j;
        }
    }
    if mantissa_end == digits_start
        || (mantissa_end == digits_start + 1 && bytes[digits_start] == b'.')
    {
        return 0.0;
    }
    let mut end = mantissa_end;
    if i < bytes.len() && (bytes[i] == b'e' || bytes[i] == b'E') {
        let mut j = i + 1;
        if matches!(bytes.get(j), Some(b'+' | b'-')) {
            j += 1;
        }
        let exponent_start = j;
        while j < bytes.len() && bytes[j].is_ascii_digit() {
            j += 1;
        }
        if j > exponent_start {
            end = j;
        }
    }
    rest[..end].parse::<f64>().unwrap_or(0.0)
}

/// Every comma-separated SVLEN value as Perl reads it (`split(',', $svlen)` then the
/// numeric reading of each); a key written with an empty value is one value, 0.
pub fn svlen_values(svlen: &str) -> Vec<i64> {
    svlen.split(',').map(perl_numeric_value).collect()
}

/// How a record's class reads its span.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpanClass {
    /// `chromosome_breakpoint`: INFO/END is deleted before the end is read.
    Breakpoint,
    /// A tandem repeat or tandem duplication: the largest SVLEN, signed.
    Tandem,
    /// Every other structural class.
    Other,
}

impl SpanClass {
    /// The span rule a variant class reads under.
    pub fn of(class: VariantClass) -> Self {
        match class {
            VariantClass::Translocation => SpanClass::Breakpoint,
            VariantClass::TandemRepeat | VariantClass::TandemDuplication => SpanClass::Tandem,
            _ => SpanClass::Other,
        }
    }
}

/// The first base of a structural record at VCF position `pos`: the base after the
/// anchor for a class that spans bases (BaseVCF4.pm:260-264, `$start++`), `pos` itself
/// for one that does not.
pub fn structural_start(pos: u64, class: VariantClass) -> u64 {
    if class.is_span_type() {
        pos.saturating_add(1)
    } else {
        pos
    }
}

/// Why Ensembl VEP skips a record on its span.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpanError {
    /// INFO/END is not a coordinate (`validate_vf`: "Invalid start or end coordinate").
    EndNotACoordinate(String),
    /// The end lies more than one base before the start (`validate_vf`: `start > end+1`).
    StartPastEnd { start: u64, end: u64 },
    /// A tandem class whose largest SVLEN is negative ends before its start.
    NegativeTandemLength(i64),
    /// The record's first base, `pos + 1`, lies past [`MAX_COORD`].
    StartPastCeiling(u64),
}

impl fmt::Display for SpanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SpanError::EndNotACoordinate(text) => write!(f, "Invalid END: '{text}'"),
            SpanError::StartPastEnd { start, end } => {
                write!(f, "start > end+1 : (START={start}, END={end})")
            }
            SpanError::NegativeTandemLength(len) => write!(
                f,
                "Invalid SVLEN for a tandem repeat or a tandem duplication: '{len}'"
            ),
            SpanError::StartPastCeiling(pos) => write!(f, "Invalid POS: '{pos}'"),
        }
    }
}

impl std::error::Error for SpanError {}

/// The end of the symbolic record at VCF position `pos` (its first base is `pos + 1`)
/// with the INFO values `svlen` and `end` as written (`END=` is `Some("")`), `None` for
/// an absent key or a bare `SVLEN` or `END` flag with no `=`.
///
/// With SVLEN defined the end is `pos + |SVLEN|` from the first value, or from the
/// largest for a tandem class, which refuses a negative one; SVLEN=0 or a value that
/// is not a number ends the record at `pos`, one base before its start, the pair
/// Ensembl VEP carries for an insertion of unknown length. Without SVLEN a breakpoint
/// ends at its start unless END is written empty or `0`, the two values Perl reads as
/// false and so never replaces (the record is then skipped as Ensembl VEP skips it);
/// every other class ends at END, read with the blanks around it dropped, refused when
/// it is not a coordinate or lies more than one base before the start, and at the start
/// when END is absent. The end saturates at [`MAX_COORD`].
pub fn structural_end(
    pos: u64,
    svlen: Option<&str>,
    end: Option<&str>,
    class: SpanClass,
) -> Result<u64, SpanError> {
    let start = pos.saturating_add(1);
    if let Some(svlen) = svlen {
        let values = svlen_values(svlen);
        let len = if class == SpanClass::Tandem {
            let largest = values.iter().copied().max().unwrap_or(0);
            if largest < 0 {
                return Err(SpanError::NegativeTandemLength(largest));
            }
            largest
        } else {
            values.first().copied().unwrap_or(0)
        };
        return Ok(pos.saturating_add(len.unsigned_abs()).min(MAX_COORD));
    }
    match (class, end) {
        (SpanClass::Breakpoint, Some("")) => Err(SpanError::EndNotACoordinate(String::new())),
        (SpanClass::Breakpoint, Some("0")) => Err(SpanError::StartPastEnd { start, end: 0 }),
        (SpanClass::Breakpoint, _) | (_, None) => Ok(start),
        (_, Some(text)) => {
            // `looks_like_number` (Parser.pm:475) admits blanks around the number.
            let end = parse_coordinate(text.trim_ascii())
                .ok_or_else(|| SpanError::EndNotACoordinate(text.to_string()))?;
            if end.saturating_add(1) < start {
                return Err(SpanError::StartPastEnd { start, end });
            }
            Ok(end)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Perl's `abs()` numifies the string: Ensembl VEP 116.2 ends a `<DUP>` at POS 30425650
    /// at 30425950 for `SVLEN=300abc`, `SVLEN= 300` and `SVLEN=3e2`, and at 30425650 (the
    /// insertion pair) for `SVLEN=abc`, `SVLEN=` and `SVLEN=0`.
    #[test]
    fn a_perl_numeric_value_is_the_leading_number_or_zero() {
        for (text, value) in [
            ("300", 300),
            ("-300", -300),
            ("+300", 300),
            ("300abc", 300),
            (" 300", 300),
            ("3e2", 300),
            ("3E+2", 300),
            ("300.7", 300),
            (".5", 0),
            ("abc", 0),
            ("", 0),
            (".", 0),
            ("e5", 0),
            ("-", 0),
            ("1e", 1),
            ("99999999999999999999", i64::MAX),
            ("1e400", i64::MAX),
            ("-1e400", i64::MIN),
        ] {
            assert_eq!(perl_numeric_value(text), value, "{text:?}");
        }
        assert_eq!(perl_numeric_f64("2abc"), 2.0);
        assert_eq!(perl_numeric_f64("1.9"), 1.9);
        assert_eq!(perl_numeric_f64("abc"), 0.0);
        assert_eq!(svlen_values(""), vec![0]);
        assert_eq!(svlen_values(",300"), vec![0, 300]);
        assert_eq!(svlen_values("-300,400"), vec![-300, 400]);
    }

    /// `get_end` (BaseVCF4.pm:302-315): SVLEN first, by its first value, else END, else the
    /// start; Ensembl VEP 116.2 writes `21:30425651-30425950` for `<DEL>` with
    /// `SVLEN=-100;SVLEN=-300` read as its last value by the caller, `21:30425651-30425950`
    /// for `SVLEN=-300,400`, and the pair `21:30425650-30425651` for `SVLEN=0`, `SVLEN=abc`
    /// and `END=30425650`.
    #[test]
    fn svlen_precedes_end_and_a_zero_length_is_the_pair_before_the_start() {
        let pos = 30425650;
        assert_eq!(
            structural_end(pos, Some("-300"), Some("30425700"), SpanClass::Other),
            Ok(30425950)
        );
        assert_eq!(
            structural_end(pos, Some("-300,400"), None, SpanClass::Other),
            Ok(30425950)
        );
        assert_eq!(
            structural_end(pos, Some("-300,400"), None, SpanClass::Tandem),
            Ok(30426050)
        );
        assert_eq!(
            structural_end(pos, Some("abc"), Some("30425899"), SpanClass::Other),
            Ok(pos)
        );
        assert_eq!(
            structural_end(pos, Some("0"), None, SpanClass::Other),
            Ok(pos)
        );
        assert_eq!(
            structural_end(pos, Some(""), None, SpanClass::Other),
            Ok(pos)
        );
        assert_eq!(
            structural_end(pos, None, Some("30425650"), SpanClass::Other),
            Ok(pos)
        );
        assert_eq!(
            structural_end(pos, None, Some("30425651"), SpanClass::Other),
            Ok(pos + 1)
        );
        assert_eq!(
            structural_end(pos, None, Some("+30425899"), SpanClass::Other),
            Ok(30425899)
        );
        assert_eq!(
            structural_end(pos, None, None, SpanClass::Other),
            Ok(pos + 1)
        );
        assert_eq!(
            structural_end(pos, Some("-300"), None, SpanClass::Tandem),
            Err(SpanError::NegativeTandemLength(-300))
        );
    }

    /// `validate_vf` (Parser.pm:475-488) skips the record: Ensembl VEP 116.2 writes
    /// `start > end+1 : (START=30425651, END=0)` for `<DUP>` with `END=0`, the same for
    /// `END=30425000`, and `Invalid start '30425651' or end '' coordinate` for `END=`.
    #[test]
    fn an_end_before_the_start_or_not_a_coordinate_is_refused() {
        let pos = 30425650;
        assert_eq!(
            structural_end(pos, None, Some("0"), SpanClass::Other),
            Err(SpanError::StartPastEnd {
                start: pos + 1,
                end: 0
            })
        );
        assert_eq!(
            structural_end(pos, None, Some("30425000"), SpanClass::Other),
            Err(SpanError::StartPastEnd {
                start: pos + 1,
                end: 30425000
            })
        );
        for text in ["", "abc", "-5", "4611686018427387904"] {
            assert_eq!(
                structural_end(pos, None, Some(text), SpanClass::Other),
                Err(SpanError::EndNotACoordinate(text.to_string())),
                "{text:?}"
            );
        }
        assert_eq!(
            SpanError::StartPastEnd {
                start: 30425651,
                end: 0
            }
            .to_string(),
            "start > end+1 : (START=30425651, END=0)"
        );
    }

    /// A breakpoint's END is deleted before the end is read again (Parser/VCF.pm:487-489),
    /// so Ensembl VEP 116.2 writes `21:30425651` for `<BND>` with `END=30425750` and with
    /// `END=abc`; the deletion is skipped for an END Perl reads as false, so `END=0` keeps
    /// the 0 (`start > end+1`) and `END=` the empty string (invalid coordinate), and the
    /// record is skipped. With SVLEN the span is SVLEN's whatever END says.
    #[test]
    fn a_breakpoint_reads_no_end_unless_it_is_empty_or_zero() {
        let pos = 30425650;
        for end in [None, Some("30425750"), Some("abc"), Some("-5")] {
            assert_eq!(
                structural_end(pos, None, end, SpanClass::Breakpoint),
                Ok(pos + 1),
                "{end:?}"
            );
        }
        assert_eq!(
            structural_end(pos, None, Some("0"), SpanClass::Breakpoint),
            Err(SpanError::StartPastEnd {
                start: pos + 1,
                end: 0
            })
        );
        assert_eq!(
            structural_end(pos, None, Some(""), SpanClass::Breakpoint),
            Err(SpanError::EndNotACoordinate(String::new()))
        );
        assert_eq!(
            structural_end(pos, Some("100"), Some("0"), SpanClass::Breakpoint),
            Ok(pos + 100)
        );
        assert_eq!(
            structural_end(pos, Some("0"), Some("0"), SpanClass::Breakpoint),
            Ok(pos)
        );
    }

    #[test]
    fn a_span_saturates_at_the_ceiling() {
        assert_eq!(
            structural_end(
                MAX_COORD - 10,
                Some("-9223372036854775808"),
                None,
                SpanClass::Other
            ),
            Ok(MAX_COORD)
        );
    }
}
