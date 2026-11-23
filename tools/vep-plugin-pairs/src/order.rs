// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The byte rules shared by every pass: the order `sort` gives tuple lines under C
//! collation, the reading of a Location, and the encodings of the scratch.

use std::cmp::Ordering;

/// A Location longer than this is left to the run's resolver, which keeps a coordinate of
/// more digits than CPython's `int()` converts as written.
const PLAIN_LOCATION_BYTES: usize = 255;

/// The order of two columns that each end at a tab inside a whole line: where one is the
/// start of the other, the tab that ends it is compared with the other's next byte.
pub fn cmp_column(a: &[u8], b: &[u8]) -> Ordering {
    let shared = a.len().min(b.len());
    a[..shared]
        .cmp(&b[..shared])
        .then_with(|| match a.len().cmp(&b.len()) {
            Ordering::Less => b'\t'.cmp(&b[shared]),
            Ordering::Greater => a[shared].cmp(&b'\t'),
            Ordering::Equal => Ordering::Equal,
        })
}

/// `bytes.strip()` with no argument: Python's ASCII whitespace, the vertical tab included.
pub fn strip(mut bytes: &[u8]) -> &[u8] {
    const SPACE: &[u8] = b" \t\n\r\x0b\x0c";
    while let [first, rest @ ..] = bytes {
        if !SPACE.contains(first) {
            break;
        }
        bytes = rest;
    }
    while let [rest @ .., last] = bytes {
        if !SPACE.contains(last) {
            break;
        }
        bytes = rest;
    }
    bytes
}

/// The Location with each coordinate written as its number; `None` for a Location left to
/// the run's resolver: one with a byte outside ASCII or longer than `PLAIN_LOCATION_BYTES`.
pub fn canonical_location(written: &[u8]) -> Option<Vec<u8>> {
    let location = strip(written);
    if location.len() > PLAIN_LOCATION_BYTES || !location.is_ascii() {
        return None;
    }
    let Some(colon) = memchr::memchr(b':', location) else {
        return Some(location.to_vec());
    };
    let mut canonical = location[..=colon].to_vec();
    for (index, coordinate) in location[colon + 1..]
        .split(|&byte| byte == b'-')
        .enumerate()
    {
        if index > 0 {
            canonical.push(b'-');
        }
        let digits = coordinate.strip_prefix(b"+").unwrap_or(coordinate);
        if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
            canonical.extend_from_slice(coordinate);
            continue;
        }
        let significant = digits.iter().position(|&digit| digit != b'0');
        canonical.extend_from_slice(significant.map_or(b"0", |at| &digits[at..]));
    }
    Some(canonical)
}

/// A value with the percent-encoding of `=` and `;` undone, as `decode_extra` of
/// `tools/vep-formats-compare` reads it: `%3D` first, then `%3B`, upper case alone.
pub fn decode_extra(value: &[u8], decoded: &mut Vec<u8>) {
    decoded.clear();
    let mut rest = value;
    while let Some(at) = memchr::memchr(b'%', rest) {
        let escaped = match rest.get(at + 1..at + 3) {
            Some(b"3D") => Some(b'='),
            Some(b"3B") => Some(b';'),
            _ => None,
        };
        decoded.extend_from_slice(&rest[..at]);
        match escaped {
            Some(byte) => {
                decoded.push(byte);
                rest = &rest[at + 3..];
            }
            None => {
                decoded.push(b'%');
                rest = &rest[at + 1..];
            }
        }
    }
    decoded.extend_from_slice(rest);
}

pub fn put_varint(out: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        out.push(value as u8 | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

pub fn take_varint(bytes: &[u8], at: &mut usize) -> u64 {
    let mut value = 0u64;
    let mut shift = 0;
    loop {
        let byte = bytes[*at];
        *at += 1;
        value |= u64::from(byte & 0x7f) << shift;
        if byte < 0x80 {
            return value;
        }
        shift += 7;
    }
}

pub fn hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut text = String::with_capacity(bytes.len() * 2);
    for &byte in bytes {
        text.push(DIGITS[usize::from(byte >> 4)] as char);
        text.push(DIGITS[usize::from(byte & 0xf)] as char);
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    fn canonical(written: &str) -> Option<String> {
        canonical_location(written.as_bytes()).map(|bytes| String::from_utf8(bytes).unwrap())
    }

    #[test]
    fn location_coordinates_are_written_as_numbers() {
        assert_eq!(canonical("21:100-0200").as_deref(), Some("21:100-200"));
        assert_eq!(
            canonical("21:30425651-+30425949").as_deref(),
            Some("21:30425651-30425949")
        );
        assert_eq!(canonical(" 21:000 ").as_deref(), Some("21:0"));
        assert_eq!(canonical("21:+").as_deref(), Some("21:+"));
        assert_eq!(canonical("21:++5").as_deref(), Some("21:++5"));
        assert_eq!(canonical("21:-05").as_deref(), Some("21:-5"));
        assert_eq!(canonical("21:30425949.5").as_deref(), Some("21:30425949.5"));
        assert_eq!(canonical("HLA:1:05").as_deref(), Some("HLA:1:05"));
        assert_eq!(
            canonical("no_coordinates").as_deref(),
            Some("no_coordinates")
        );
        assert_eq!(canonical("").as_deref(), Some(""));
    }

    #[test]
    fn location_outside_ascii_or_over_long_is_left_to_the_resolver() {
        assert_eq!(canonical("21:²"), None);
        assert_eq!(canonical(&format!("21:{}", "9".repeat(300))), None);
    }

    #[test]
    fn strip_takes_pythons_ascii_whitespace() {
        assert_eq!(strip(b" \x0b\x0c\r\tA T\x0b "), b"A T");
        assert_eq!(strip(b"\x1cA\x1c"), b"\x1cA\x1c");
        assert_eq!(strip(b" \t "), b"");
    }

    #[test]
    fn percent_encoding_of_equals_and_semicolon_is_undone() {
        let mut decoded = Vec::new();
        for (value, expected) in [
            ("a%3Bb%3Dc", "a;b=c"),
            ("%3d%3b", "%3d%3b"),
            ("%3%3DB", "%3=B"),
            ("%%3D3B", "%=3B"),
            ("100%", "100%"),
            ("%3", "%3"),
        ] {
            decode_extra(value.as_bytes(), &mut decoded);
            assert_eq!(decoded, expected.as_bytes(), "{value}");
        }
    }

    #[test]
    fn columns_order_as_their_lines_do_under_c_collation() {
        let mut lines = ["a\tz", "a\x01\ty", "ab\tx", "a \tw", "\tv"];
        lines.sort();
        let mut columns: Vec<&[u8]> = vec![b"a", b"a\x01", b"ab", b"a ", b""];
        columns.sort_by(|a, b| cmp_column(a, b));
        let firsts: Vec<&[u8]> = lines
            .iter()
            .map(|line| line.split('\t').next().unwrap().as_bytes())
            .collect();
        assert_eq!(columns, firsts);
    }

    #[test]
    fn varints_round_trip_and_bytes_are_written_as_hex() {
        let mut bytes = Vec::new();
        let values = [0, 1, 127, 128, 300, u64::from(u32::MAX), u64::MAX];
        for value in values {
            put_varint(&mut bytes, value);
        }
        let mut at = 0;
        for value in values {
            assert_eq!(take_varint(&bytes, &mut at), value);
        }
        assert_eq!(at, bytes.len());
        assert_eq!(hex(b"\x00K\xc3\xa9y\xff"), "004bc3a979ff");
    }
}
