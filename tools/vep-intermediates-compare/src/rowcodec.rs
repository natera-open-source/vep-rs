// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row each side hands to `sort`: three key fields, a sequence number, then the value
//! fields, on one line.
//!
//! A field is a marker byte (absent or present), the field's bytes with five of them
//! escaped, and a terminator. The escapes keep the order of the bytes they replace and the
//! terminator is below every byte an escaped field holds, so the bytewise order of two
//! rows is the order of their keys as tuples of strings compared by code point, a shorter
//! string before any string it starts, and then the order of their sequence numbers. No
//! row holds a line feed or a NUL.

use memchr::{memchr, memchr2};
use std::borrow::Cow;

const TERMINATOR: u8 = 0x01;
/// A field the row lacks; it sorts before every present field.
const ABSENT: u8 = b'-';
const PRESENT: u8 = b'=';
/// Introduces the escapes of 0x00, 0x01 and 0x02.
const LOW_ESCAPE: u8 = 0x02;
/// Introduces the escapes of 0x0a and 0x0b.
const LINE_ESCAPE: u8 = 0x0b;
pub const SEQUENCE_DIGITS: usize = 16;

fn needs_escape(byte: u8) -> bool {
    matches!(byte, 0x00..=0x02 | 0x0a | 0x0b)
}

pub fn push_field(out: &mut Vec<u8>, field: Option<&[u8]>) {
    match field {
        None => out.push(ABSENT),
        Some(bytes) => {
            out.push(PRESENT);
            if bytes.iter().any(|b| needs_escape(*b)) {
                for &byte in bytes {
                    match byte {
                        0x00..=0x02 => out.extend_from_slice(&[LOW_ESCAPE, byte + 2]),
                        0x0a | 0x0b => out.extend_from_slice(&[LINE_ESCAPE, byte - 8]),
                        _ => out.push(byte),
                    }
                }
            } else {
                out.extend_from_slice(bytes);
            }
        }
    }
    out.push(TERMINATOR);
}

/// Appends `n` as sixteen hexadecimal digits, so that digit order is numeric order.
pub fn push_sequence(out: &mut Vec<u8>, n: u64) {
    out.extend_from_slice(format!("{n:016x}").as_bytes());
}

/// The length of the three key fields at the start of `row`.
pub fn key_len(row: &[u8]) -> usize {
    let mut end = 0;
    for _ in 0..3 {
        end += memchr(TERMINATOR, &row[end..]).expect("a row starts with three key fields") + 1;
    }
    end
}

/// Takes the next field off the front of `row`.
pub fn take_field<'a>(row: &mut &'a [u8]) -> Option<Cow<'a, [u8]>> {
    let end = memchr(TERMINATOR, row).expect("a field ends with its terminator");
    let (marker, body) = (row[0], &row[1..end]);
    *row = &row[end + 1..];
    if marker == ABSENT {
        return None;
    }
    if memchr2(LOW_ESCAPE, LINE_ESCAPE, body).is_none() {
        return Some(Cow::Borrowed(body));
    }
    let mut bytes = Vec::with_capacity(body.len());
    let mut escaped = body.iter();
    while let Some(&byte) = escaped.next() {
        bytes.push(match byte {
            LOW_ESCAPE => escaped.next().expect("an escape has two bytes") - 2,
            LINE_ESCAPE => escaped.next().expect("an escape has two bytes") + 8,
            _ => byte,
        });
    }
    Some(Cow::Owned(bytes))
}

/// Whether two keys differ first in a field that one of them lacks. Python orders a tuple
/// by its first differing element and has no order between `None` and a string, so two
/// such keys in one comparison end `sorted()` on a `TypeError`.
pub fn unorderable(a: &[u8], b: &[u8]) -> bool {
    let (mut a, mut b) = (a, b);
    for _ in 0..3 {
        let a_end = memchr(TERMINATOR, a).expect("a key field ends with its terminator");
        let b_end = memchr(TERMINATOR, b).expect("a key field ends with its terminator");
        if a[..a_end] != b[..b_end] {
            return a[0] != b[0];
        }
        a = &a[a_end + 1..];
        b = &b[b_end + 1..];
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(key: [Option<&[u8]>; 3], sequence: u64) -> Vec<u8> {
        let mut out = Vec::new();
        for field in key {
            push_field(&mut out, field);
        }
        push_sequence(&mut out, sequence);
        out
    }

    #[test]
    fn fields_round_trip_with_every_escaped_byte() {
        let fields: [Option<&[u8]>; 6] = [
            Some(b""),
            None,
            Some(b"1:100-102"),
            Some(b"\x00\x01\x02\x03\t\n\x0b\x0c\xff"),
            Some(b"-"),
            Some("\u{e9}=".as_bytes()),
        ];
        let mut encoded = Vec::new();
        for field in fields {
            push_field(&mut encoded, field);
        }
        assert!(!encoded.contains(&b'\n') && !encoded.contains(&0));
        let mut rest = encoded.as_slice();
        for field in fields {
            assert_eq!(take_field(&mut rest).as_deref(), field);
        }
        assert!(rest.is_empty());
    }

    #[test]
    fn row_order_is_tuple_order_then_sequence_order() {
        let alphabet = [0x00, 0x01, 0x02, 0x03, b'\t', b'\n', 0x0b, 0x0c, b'a', 0xff];
        let mut strings: Vec<Vec<u8>> = vec![Vec::new()];
        for length in 1..=3 {
            let shorter: Vec<Vec<u8>> = strings
                .iter()
                .filter(|s| s.len() == length - 1)
                .cloned()
                .collect();
            for prefix in shorter {
                for byte in alphabet {
                    strings.push([prefix.as_slice(), &[byte]].concat());
                }
            }
        }
        assert_eq!(strings.len(), 1111);
        let single = |s: &[u8]| {
            let mut out = Vec::new();
            push_field(&mut out, Some(s));
            out
        };
        let encoded: Vec<Vec<u8>> = strings.iter().map(|s| single(s)).collect();
        for (a, ea) in strings.iter().zip(&encoded) {
            for (b, eb) in strings.iter().zip(&encoded) {
                assert_eq!(ea.cmp(eb), a.cmp(b), "{a:?} against {b:?}");
            }
        }
        // A key whose first field starts another's sorts first whatever follows it, and
        // rows of one key sort by sequence.
        let ordered = [
            row([None, Some(b"z"), Some(b"z")], 9),
            row([Some(b"1:100"), None, None], 9),
            row([Some(b"1:100"), Some(b"A"), None], 9),
            row([Some(b"1:100"), Some(b"A"), Some(b"ENST2")], 3),
            row([Some(b"1:100"), Some(b"A"), Some(b"ENST2")], 0x10),
            row([Some(b"1:100\x00"), Some(b""), Some(b"")], 0),
            row([Some(b"1:1000"), Some(b""), Some(b"")], 0),
        ];
        let mut sorted = ordered.to_vec();
        sorted.sort();
        assert_eq!(sorted, ordered);
    }

    #[test]
    fn key_length_stops_after_the_third_field() {
        let mut encoded = row([Some(b"1:100"), None, Some(b"\x01")], 7);
        let key = encoded.len() - SEQUENCE_DIGITS;
        push_field(&mut encoded, Some(b"value"));
        assert_eq!(key_len(&encoded), key);
        let sequence = std::str::from_utf8(&encoded[key..key + SEQUENCE_DIGITS]).unwrap();
        assert_eq!(u64::from_str_radix(sequence, 16), Ok(7));
        assert_eq!(encoded.len(), key + SEQUENCE_DIGITS + b"=value".len() + 1);
    }

    #[test]
    fn keys_are_unorderable_when_they_first_differ_in_a_field_one_lacks() {
        let key = |fields: [Option<&[u8]>; 3]| {
            let encoded = row(fields, 0);
            encoded[..key_len(&encoded)].to_vec()
        };
        let short = key([Some(b"1:100"), None, None]);
        let other_location = key([Some(b"1:200"), Some(b"A"), Some(b"ENST1")]);
        let same_location = key([Some(b"1:100"), Some(b"A"), Some(b"ENST1")]);
        let also_short = key([Some(b"1:100"), None, None]);
        let lacks_transcript = key([Some(b"1:100"), Some(b"A"), None]);
        assert!(!unorderable(&short, &other_location));
        assert!(unorderable(&short, &same_location));
        assert!(unorderable(&same_location, &short));
        assert!(!unorderable(&short, &also_short));
        assert!(unorderable(&lacks_transcript, &same_location));
        assert!(unorderable(&short, &lacks_transcript));
        assert!(!unorderable(&same_location, &other_location));
    }
}
