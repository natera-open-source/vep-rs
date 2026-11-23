// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The canonical form of a `--tab` cell: a DOUBLE or `LIST<DOUBLE>` cell with each number as
//! its shortest round-trip decimal, an INT32 cell as its integer, every other cell as
//! written. A token outside the number grammar of a DOUBLE column and a token of an INT32
//! column that is no INT32 stay as written too. The columns of each type are named here as
//! the writer's cell contract names them.

use std::io::Write;

use crate::pyfmt::decimal_digit;

pub const NULL_TEXT: &[u8] = b"-";

/// The cache columns behind the population frequency fields; the frequency column of each
/// is `<column>_AF`, `AF` for `AF`.
const FREQUENCY_POPULATIONS: [&str; 28] = [
    "AF",
    "AFR",
    "AMR",
    "ASN",
    "EAS",
    "EUR",
    "SAS",
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_MID",
    "gnomADe_NFE",
    "gnomADe_REMAINING",
    "gnomADe_SAS",
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_REMAINING",
    "gnomADg_SAS",
];
const DOUBLE_COLUMNS: [&str; 2] = ["MAX_AF", "MOTIF_SCORE_CHANGE"];
const INT32_COLUMNS: [&str; 7] = [
    "STRAND",
    "DISTANCE",
    "ALLELE_NUM",
    "TSL",
    "HGVS_OFFSET",
    "OverlapBP",
    "MOTIF_POS",
];

/// How the canonical form prints the cells of a column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanonicalKind {
    Verbatim,
    Double,
    DoubleList,
    Int32,
}

fn is_double_list_column(column: &[u8]) -> bool {
    column == b"AF"
        || column.strip_suffix(b"_AF").is_some_and(|population| {
            population != b"AF"
                && FREQUENCY_POPULATIONS
                    .iter()
                    .any(|p| p.as_bytes() == population)
        })
}

/// The canonical kind of the column a header names.
pub fn canonical_kind(column: &[u8]) -> CanonicalKind {
    if DOUBLE_COLUMNS.iter().any(|c| c.as_bytes() == column) {
        CanonicalKind::Double
    } else if is_double_list_column(column) {
        CanonicalKind::DoubleList
    } else if INT32_COLUMNS.iter().any(|c| c.as_bytes() == column) {
        CanonicalKind::Int32
    } else {
        CanonicalKind::Verbatim
    }
}

/// Appends the shortest decimal that reads back to `value`, in plain notation: `0.095`,
/// `0.000098` for 9.8e-05, `1` for 1.0, `-0` for -0.0, `NaN`, `inf` and `-inf`. Of two
/// shortest decimals equally close to the value, the one of larger magnitude.
pub fn push_shortest_decimal(out: &mut Vec<u8>, value: f64) {
    write!(out, "{value}").expect("a write into memory");
}

/// A DOUBLE token read under the cell contract's number grammar.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Number {
    /// Outside the grammar: the token stays as written.
    Outside,
    Value(f64),
    /// Inside the grammar as a case-insensitive pattern reads it and no number: `inf`,
    /// `infinity` or `nan` spelt with a dotted capital or a dotless small `i`.
    Unconvertible,
}

/// `token` under the grammar: an optional sign, then digits with an optional fraction or a
/// fraction alone and an optional exponent, or `inf`, `infinity` or `nan` in any case. A
/// digit is any Unicode decimal digit, and one newline may end the token.
pub fn parse_double(token: &str) -> Number {
    let body = token.strip_suffix('\n').unwrap_or(token);
    let unsigned = body.strip_prefix(['+', '-']).unwrap_or(body);
    let negative = body.starts_with('-');
    if let Some(word) = special_word(unsigned) {
        return match word {
            Word::Unconvertible => Number::Unconvertible,
            Word::Infinity if negative => Number::Value(f64::NEG_INFINITY),
            Word::Infinity => Number::Value(f64::INFINITY),
            Word::NotANumber => Number::Value(f64::NAN),
        };
    }
    let mut ascii = String::with_capacity(body.len());
    if negative {
        ascii.push('-');
    }
    let mut rest = unsigned.chars().peekable();
    let digits = |rest: &mut std::iter::Peekable<std::str::Chars>, ascii: &mut String| {
        let mut count = 0usize;
        while let Some(digit) = rest.peek().and_then(|c| decimal_digit(*c)) {
            ascii.push((b'0' + digit) as char);
            rest.next();
            count += 1;
        }
        count
    };
    let whole = digits(&mut rest, &mut ascii);
    let mut fraction = 0;
    if rest.peek() == Some(&'.') {
        rest.next();
        ascii.push('.');
        fraction = digits(&mut rest, &mut ascii);
    }
    if whole == 0 && fraction == 0 {
        return Number::Outside;
    }
    if matches!(rest.peek(), Some('e' | 'E')) {
        rest.next();
        ascii.push('e');
        if let Some(sign) = rest.next_if(|c| matches!(c, '+' | '-')) {
            ascii.push(sign);
        }
        if digits(&mut rest, &mut ascii) == 0 {
            return Number::Outside;
        }
    }
    if rest.next().is_some() {
        return Number::Outside;
    }
    ascii.parse().map_or(Number::Outside, Number::Value)
}

enum Word {
    Infinity,
    NotANumber,
    Unconvertible,
}

/// `inf`, `infinity` or `nan` as a case-insensitive pattern matches them: the ASCII letters
/// in either case, and for `i` also U+0130 and U+0131, which no number spells.
fn special_word(text: &str) -> Option<Word> {
    let mut folded = String::with_capacity(8);
    let mut ascii = true;
    for c in text.chars() {
        match c {
            '\u{130}' | '\u{131}' => {
                ascii = false;
                folded.push('i');
            }
            c if c.is_ascii_alphabetic() => folded.push(c.to_ascii_lowercase()),
            _ => return None,
        }
        if folded.len() > 8 {
            return None;
        }
    }
    let word = match folded.as_str() {
        "inf" | "infinity" => Word::Infinity,
        "nan" => Word::NotANumber,
        _ => return None,
    };
    Some(if ascii { word } else { Word::Unconvertible })
}

/// A token of a DOUBLE cell that the grammar takes and no number spells; the run ends on
/// it.
#[derive(Debug, PartialEq, Eq)]
pub struct Unconvertible(pub Vec<u8>);

fn push_canonical_number(out: &mut Vec<u8>, token: &[u8]) -> Result<(), Unconvertible> {
    let number = std::str::from_utf8(token).map_or(Number::Outside, parse_double);
    match number {
        Number::Outside => out.extend_from_slice(token),
        Number::Value(value) => push_shortest_decimal(out, value),
        Number::Unconvertible => return Err(Unconvertible(token.to_vec())),
    }
    Ok(())
}

/// An INT32 token: an optional sign and ASCII digits print as their integer (`007` is `7`),
/// as an INT32 column prints the value. A token that is no INT32 stays as written, as the
/// string column a plugin or a `--custom` source declares under the name stores it; an
/// INT32 column stores NULL for it, which prints `-` and so differs.
fn push_canonical_int32(out: &mut Vec<u8>, token: &[u8]) {
    let (negative, digits) = match token.split_first() {
        Some((b'-', digits)) => (true, digits),
        Some((b'+', digits)) => (false, digits),
        _ => (false, token),
    };
    let value = (!digits.is_empty() && digits.iter().all(u8::is_ascii_digit))
        .then(|| {
            let zeros = digits.iter().take_while(|d| **d == b'0').count();
            &digits[zeros..]
        })
        // Ten digits hold every INT32.
        .filter(|significant| significant.len() <= 10)
        .map(|significant| {
            let magnitude = significant
                .iter()
                .fold(0i64, |n, d| n * 10 + i64::from(d - b'0'));
            if negative {
                -magnitude
            } else {
                magnitude
            }
        })
        .and_then(|value| i32::try_from(value).ok());
    match value {
        Some(value) => out.extend_from_slice(itoa::Buffer::new().format(value).as_bytes()),
        None => out.extend_from_slice(token),
    }
}

/// Appends the canonical form of one cell of a column of `kind`.
pub fn push_canonical_cell(
    out: &mut Vec<u8>,
    kind: CanonicalKind,
    cell: &[u8],
) -> Result<(), Unconvertible> {
    if cell == NULL_TEXT {
        out.extend_from_slice(cell);
        return Ok(());
    }
    match kind {
        CanonicalKind::Verbatim => out.extend_from_slice(cell),
        CanonicalKind::Double => push_canonical_number(out, cell)?,
        CanonicalKind::DoubleList => {
            for (i, token) in cell.split(|b| *b == b',').enumerate() {
                if i > 0 {
                    out.push(b',');
                }
                push_canonical_number(out, token)?;
            }
        }
        CanonicalKind::Int32 => push_canonical_int32(out, cell),
    }
    Ok(())
}

/// Appends the canonical form of one data line, without its line end: the cells under the
/// header's columns each in their canonical form, a cell past the last column as written.
pub fn push_canonical_line(
    out: &mut Vec<u8>,
    kinds: &[CanonicalKind],
    line: &[u8],
) -> Result<(), Unconvertible> {
    for (i, cell) in line.split(|b| *b == b'\t').enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        let kind = kinds.get(i).copied().unwrap_or(CanonicalKind::Verbatim);
        push_canonical_cell(out, kind, cell)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shortest(value: f64) -> String {
        let mut out = Vec::new();
        push_shortest_decimal(&mut out, value);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn shortest_decimal_prints_the_writers_text() {
        let tiny = format!("0.{}5", "0".repeat(323));
        let cases: [(f64, &str); 14] = [
            (0.095, "0.095"),
            (1.0, "1"),
            (9.8e-05, "0.000098"),
            (0.0, "0"),
            (-0.0, "-0"),
            (1e16, "10000000000000000"),
            (1e21, "1000000000000000000000"),
            (0.1 + 0.2, "0.30000000000000004"),
            (5e-324, &tiny),
            (f64::NAN, "NaN"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            // Two shortest decimals equally close: the one of larger magnitude.
            (-126430982921026.62, "-126430982921026.63"),
            (126430982921026.62, "126430982921026.63"),
        ];
        for (value, text) in cases {
            assert_eq!(shortest(value), text, "{value:?}");
            if !matches!(text, "NaN" | "inf" | "-inf") {
                let back: f64 = text.parse().unwrap();
                assert_eq!(back.to_bits(), value.to_bits(), "{text}");
            }
        }
    }

    #[test]
    fn parse_double_follows_the_contracts_number_grammar() {
        for token in [
            "0.0950",
            "1",
            "-0.012",
            "9.8e-05",
            "1E5",
            ".5",
            "3.",
            "+1",
            "inf",
            "-Infinity",
            "nan",
        ] {
            assert!(matches!(parse_double(token), Number::Value(_)), "{token}");
        }
        for token in ["n/a", "1_000", " 1", "0x10", "1e", "-", "", "1,2"] {
            assert_eq!(parse_double(token), Number::Outside, "{token}");
        }
    }

    #[test]
    fn parse_double_reads_what_a_python_pattern_and_float_read() {
        assert_eq!(parse_double("1.50\n"), Number::Value(1.5));
        assert_eq!(parse_double("1.50\n\n"), Number::Outside);
        assert_eq!(parse_double("1.50\r"), Number::Outside);
        assert_eq!(parse_double("1.50 "), Number::Outside);
        assert_eq!(parse_double("\n"), Number::Outside);
        assert_eq!(parse_double("."), Number::Outside);
        assert_eq!(parse_double("1.2.3"), Number::Outside);
        assert_eq!(parse_double("1e+"), Number::Outside);
        assert_eq!(parse_double("e5"), Number::Outside);
        assert_eq!(parse_double("infinit"), Number::Outside);
        assert_eq!(parse_double("\u{661}\u{662}.\u{665}"), Number::Value(12.5));
        assert_eq!(parse_double("\u{FF11}e\u{FF12}"), Number::Value(100.0));
        assert_eq!(parse_double("1e999"), Number::Value(f64::INFINITY));
        assert_eq!(parse_double("-1e999"), Number::Value(f64::NEG_INFINITY));
        assert_eq!(parse_double("1e-999"), Number::Value(0.0));
        assert_eq!(parse_double("+iNfInItY"), Number::Value(f64::INFINITY));
        assert_eq!(parse_double("-inf\n"), Number::Value(f64::NEG_INFINITY));
        assert!(matches!(parse_double("-NaN"), Number::Value(v) if v.is_nan()));
        assert_eq!(parse_double("\u{131}nf"), Number::Unconvertible);
        assert_eq!(parse_double("\u{130}NF"), Number::Unconvertible);
        assert_eq!(
            parse_double("-\u{131}nf\u{131}n\u{131}ty\n"),
            Number::Unconvertible
        );
        assert_eq!(parse_double("n\u{e4}n"), Number::Outside);
    }

    fn cell(column: &str, text: &str) -> String {
        let mut out = Vec::new();
        push_canonical_cell(&mut out, canonical_kind(column.as_bytes()), text.as_bytes()).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn canonical_tab_cell_rewrites_a_double_and_an_int32_cell() {
        assert_eq!(cell("MAX_AF", "0.0950"), "0.095");
        assert_eq!(cell("MAX_AF", "1.0"), "1");
        assert_eq!(cell("MAX_AF", "9.8e-05"), "0.000098");
        assert_eq!(cell("MAX_AF", "-"), "-");
        assert_eq!(cell("MAX_AF", "n/a"), "n/a");
        assert_eq!(cell("MOTIF_SCORE_CHANGE", "0.100"), "0.1");
        assert_eq!(cell("gnomADe_AF", "0.0950,0.1,1"), "0.095,0.1,1");
        assert_eq!(cell("Feature", "0.0950"), "0.0950");
        // An INT32 cell prints as its integer, and a text that is no INT32 as written.
        assert_eq!(cell("TSL", "007"), "7");
        assert_eq!(cell("TSL", "00"), "0");
        assert_eq!(cell("STRAND", "-1"), "-1");
        assert_eq!(cell("TSL", "2147483648"), "2147483648");
        assert_eq!(cell("TSL", "1.5"), "1.5");
        assert_eq!(cell("TSL", &format!("{}7", "0".repeat(5000))), "7");
        assert_eq!(cell("TSL", &"9".repeat(5000)), "9".repeat(5000));
        assert_eq!(cell("TSL", " 7"), " 7");
        assert_eq!(cell("TSL", "-"), "-");
        assert_eq!(cell("Feature", "007"), "007");
        let double_lists: Vec<String> = FREQUENCY_POPULATIONS
            .iter()
            .map(|population| match *population {
                "AF" => "AF".to_string(),
                population => format!("{population}_AF"),
            })
            .collect();
        assert_eq!(double_lists.len(), 28);
        for column in &double_lists {
            assert_eq!(
                canonical_kind(column.as_bytes()),
                CanonicalKind::DoubleList,
                "{column}"
            );
        }
        for column in [
            "AF",
            "AFR_AF",
            "gnomADe_AF",
            "gnomADe_REMAINING_AF",
            "gnomADg_AF",
            "gnomADg_SAS_AF",
        ] {
            assert!(double_lists.iter().any(|c| c == column), "{column}");
        }
        assert_eq!(canonical_kind(b"MAX_AF_POPS"), CanonicalKind::Verbatim);
        assert_eq!(canonical_kind(b"AF_AF"), CanonicalKind::Verbatim);
        assert_eq!(canonical_kind(b"_AF"), CanonicalKind::Verbatim);
    }

    #[test]
    fn an_int32_cell_keeps_its_range_and_drops_its_sign_spelling() {
        assert_eq!(cell("TSL", "2147483647"), "2147483647");
        assert_eq!(cell("TSL", "-2147483648"), "-2147483648");
        assert_eq!(cell("TSL", "-2147483649"), "-2147483649");
        assert_eq!(cell("TSL", "+7"), "7");
        assert_eq!(cell("TSL", "-0"), "0");
        assert_eq!(cell("TSL", "-00"), "0");
        assert_eq!(cell("TSL", ""), "");
        assert_eq!(cell("TSL", "+"), "+");
        assert_eq!(cell("TSL", "7\n"), "7\n");
        assert_eq!(cell("TSL", "\u{667}"), "\u{667}");
        assert_eq!(cell("MAX_AF", ""), "");
        assert_eq!(cell("AF", "1.0,,x,2.50"), "1,,x,2.5");
    }

    #[test]
    fn a_cell_of_an_int32_column_that_is_no_int32_stays_as_written() {
        // A list, as a plugin's field of the name holds one, and each of its integers alone.
        assert_eq!(cell("TSL", "1,5,2"), "1,5,2");
        assert_eq!(cell("TSL", "1"), "1");
        assert_eq!(cell("TSL", "05"), "5");
        assert_eq!(cell("DISTANCE", "1;5;2"), "1;5;2");
        assert_eq!(cell("STRAND", "+1,-1"), "+1,-1");
        assert_eq!(cell("HGVS_OFFSET", "n/a"), "n/a");
        assert_eq!(cell("TSL", "--1"), "--1");
        assert_eq!(cell("TSL", "+-1"), "+-1");
        let kinds = [CanonicalKind::Verbatim, CanonicalKind::Int32];
        let mut out = Vec::new();
        push_canonical_line(&mut out, &kinds, b"ENST1\t1,5,2").unwrap();
        assert_eq!(out, b"ENST1\t1,5,2");
    }

    #[test]
    fn a_line_is_rewritten_under_its_header_columns_only() {
        let kinds: Vec<CanonicalKind> = ["MAX_AF", "TSL"]
            .iter()
            .map(|c| canonical_kind(c.as_bytes()))
            .collect();
        let line = |text: &str| {
            let mut out = Vec::new();
            push_canonical_line(&mut out, &kinds, text.as_bytes()).map(|()| out)
        };
        assert_eq!(line("1.50\t007\t1.50\t007").unwrap(), b"1.5\t7\t1.50\t007");
        assert_eq!(line("1.50").unwrap(), b"1.5");
        assert_eq!(line("\t").unwrap(), b"\t");
        assert_eq!(
            line("\u{131}nf\t1"),
            Err(Unconvertible("\u{131}nf".as_bytes().to_vec()))
        );
    }
}
