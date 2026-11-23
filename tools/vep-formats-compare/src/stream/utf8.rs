// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The UTF-8 check of a piece of an output, a failure worded as CPython words a
//! `UnicodeDecodeError` of its `utf-8` codec.

fn continues(byte: u8) -> bool {
    byte & 0xc0 == 0x80
}

/// The bytes at the start of `rest` the codec refuses together, and its reason. `rest`
/// starts at the first byte that is not UTF-8 and runs to the end of the text.
fn refused(rest: &[u8]) -> (usize, &'static str) {
    const START: &str = "invalid start byte";
    const CONTINUATION: &str = "invalid continuation byte";
    let end = (rest.len(), "unexpected end of data");
    let lead = rest[0];
    if !(0xc2..0xf5).contains(&lead) {
        return (1, START);
    }
    let Some(&second) = rest.get(1) else {
        return end;
    };
    // The second byte's range is narrower after four lead bytes: an overlong form, a
    // surrogate and a code point past U+10FFFF are refused at the second byte.
    let second_fits = continues(second)
        && match lead {
            0xe0 => second >= 0xa0,
            0xed => second < 0xa0,
            0xf0 => second >= 0x90,
            0xf4 => second < 0x90,
            _ => true,
        };
    if !second_fits {
        return (1, CONTINUATION);
    }
    let length = match lead {
        0xc2..=0xdf => 2,
        0xe0..=0xef => 3,
        _ => 4,
    };
    for place in 2..length {
        match rest.get(place) {
            None => return end,
            Some(&byte) if !continues(byte) => return (place, CONTINUATION),
            Some(_) => {}
        }
    }
    (1, START)
}

/// `None` for UTF-8 text; else the codec's message for the first place that is not,
/// `offset` being the place of `bytes` in the whole text and `bytes` running to a line end
/// or to the end of the text.
pub(crate) fn failure(bytes: &[u8], offset: u64) -> Option<String> {
    let at = std::str::from_utf8(bytes).err()?.valid_up_to();
    let (length, reason) = refused(&bytes[at..]);
    let start = offset + at as u64;
    Some(if length == 1 {
        format!(
            "'utf-8' codec can't decode byte {:#04x} in position {start}: {reason}",
            bytes[at]
        )
    } else {
        format!(
            "'utf-8' codec can't decode bytes in position {start}-{}: {reason}",
            start + length as u64 - 1
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn message(bytes: &[u8]) -> String {
        failure(bytes, 0).unwrap()
    }

    #[test]
    fn utf8_text_passes() {
        assert_eq!(failure(b"", 0), None);
        assert_eq!(failure("5\u{2032}UTR \u{10348}\n".as_bytes(), 9), None);
        assert_eq!(failure(b"\xef\xbb\xbfa\n", 0), None);
    }

    #[test]
    fn a_refused_byte_is_named_with_its_place_and_reason() {
        for (bytes, expected) in [
            (
                &b"ab\xff\n"[..],
                "byte 0xff in position 2: invalid start byte",
            ),
            (b"ab\x80\n", "byte 0x80 in position 2: invalid start byte"),
            (b"ab\xc0\xaf", "byte 0xc0 in position 2: invalid start byte"),
            (b"ab\xc1\xbf", "byte 0xc1 in position 2: invalid start byte"),
            (b"ab\xf5\x80", "byte 0xf5 in position 2: invalid start byte"),
            (
                b"a\xf8\x88\x80\x80\x80",
                "byte 0xf8 in position 1: invalid start byte",
            ),
            (
                b"ab\xc2\n",
                "byte 0xc2 in position 2: invalid continuation byte",
            ),
            (
                b"ab\xe2\n",
                "byte 0xe2 in position 2: invalid continuation byte",
            ),
            (
                b"ab\xe0\x80\x80",
                "byte 0xe0 in position 2: invalid continuation byte",
            ),
            (
                b"ab\xed\xa0\x80",
                "byte 0xed in position 2: invalid continuation byte",
            ),
            (
                b"ab\xf0\x80\x80\x80",
                "byte 0xf0 in position 2: invalid continuation byte",
            ),
            (
                b"ab\xf4\x90\x80\x80",
                "byte 0xf4 in position 2: invalid continuation byte",
            ),
            (b"ab\xc2", "byte 0xc2 in position 2: unexpected end of data"),
            (b"ab\xf0", "byte 0xf0 in position 2: unexpected end of data"),
            (b"ab\xe1", "byte 0xe1 in position 2: unexpected end of data"),
        ] {
            assert_eq!(
                message(bytes),
                format!("'utf-8' codec can't decode {expected}")
            );
        }
    }

    #[test]
    fn refused_bytes_are_named_as_a_range() {
        for (bytes, expected) in [
            (&b"ab\xe2\x82"[..], "2-3: unexpected end of data"),
            (b"ab\xe2\x82\n", "2-3: invalid continuation byte"),
            (b"ab\xf0\x9f\x98", "2-4: unexpected end of data"),
            (b"ab\xf0\x9f\x98\n", "2-4: invalid continuation byte"),
            (b"ab\xf0\x9f\n", "2-3: invalid continuation byte"),
            (b"\xf0\x9f", "0-1: unexpected end of data"),
        ] {
            assert_eq!(
                message(bytes),
                format!("'utf-8' codec can't decode bytes in position {expected}")
            );
        }
    }

    #[test]
    fn a_place_counts_from_the_start_of_the_whole_text() {
        assert_eq!(
            failure(b"row\t\xff\n", 4_294_967_296).unwrap(),
            "'utf-8' codec can't decode byte 0xff in position 4294967300: invalid start byte"
        );
    }
}
