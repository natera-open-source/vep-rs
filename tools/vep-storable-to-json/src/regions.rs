// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `sequence_regions.tsv`: the length and the circular flag of every chromosome whose
//! transcripts state them.
//!
//! The cache dumper leaves `seq_region_length` and, where it was primed before the dump,
//! `circular` on each transcript's slice (DumpVEP/Dumper/Core.pm:122-123 deletes the
//! adaptors alone; AnnotationSource/Database.pm:114). Every transcript of a chromosome
//! shares one slice, so the first one that states a length is kept.

use crate::json::PerlStr;
use crate::perl::{Die, Perl};
use crate::storable::NodeId;

/// The file's first line.
pub const HEADER: &[u8] = b"name\tlength\tcircular\n";

/// The length and the circular flag one slice states.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRegion {
    /// The length as Perl prints the scalar.
    pub length: Vec<u8>,
    /// Whether `length` is a string under Perl's UTF-8 flag.
    pub utf8: bool,
    /// Perl's truth of `circular`; `None` where it is undef.
    pub circular: Option<bool>,
}

/// What the walk over a shard's transcripts read of their chromosome's sequence, up to
/// the first transcript that states a length.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Noted {
    /// No transcript walked states a length.
    #[default]
    Nothing,
    Region(SequenceRegion),
    /// The pattern match on a length died. The shard's walk ends there, unless a shard
    /// before it in the chromosome stated a length: no slice of this shard is read then.
    Dies(Die),
}

/// The region a transcript states: its `slice` is a reference to a hash whose class, if
/// it has one, is not `0`, and the slice's `seq_region_length` is defined and matches
/// `/^[1-9]\d*$/`.
pub fn of_transcript(perl: Perl<'_>, transcript: NodeId) -> Result<Option<SequenceRegion>, Die> {
    let slice = perl.fetch(transcript, b"slice");
    if !perl.ref_is_true(slice) {
        return Ok(None);
    }
    let Some(slice) = perl.hash_target(slice) else {
        return Ok(None);
    };
    let length = perl.fetch(slice, b"seq_region_length");
    if !perl.defined(length) {
        return Ok(None);
    }
    // The string of a reference is its address, which the pattern never matches.
    let Ok(text) = perl.string(length) else {
        return Ok(None);
    };
    if !matches_length(text.as_perl())? {
        return Ok(None);
    }
    let circular = perl.fetch(slice, b"circular");
    Ok(Some(SequenceRegion {
        length: text.bytes.into_owned(),
        utf8: text.utf8,
        circular: perl.defined(circular).then(|| perl.is_true(circular)),
    }))
}

/// Perl's `$text =~ /^[1-9]\d*$/`: a digit from 1 to 9, then digits to the end of the
/// string or to one newline that ends it. Under the UTF-8 flag `\d` is any decimal
/// digit of Unicode, and the match dies on some strings whose bytes are not UTF-8.
pub fn matches_length(text: PerlStr<'_>) -> Result<bool, Die> {
    if text.utf8 && end_hop_dies(text.bytes) {
        return Err(Die::new("Malformed UTF-8 character (fatal)"));
    }
    let body = text.bytes.strip_suffix(b"\n").unwrap_or(text.bytes);
    let Some((first, mut rest)) = body.split_first() else {
        return Ok(false);
    };
    if !(b'1'..=b'9').contains(first) {
        return Ok(false);
    }
    while let Some(&byte) = rest.first() {
        if byte.is_ascii_digit() {
            rest = &rest[1..];
        } else if !text.utf8 {
            return Ok(false);
        } else {
            match decode(rest) {
                Some((len, character)) if is_decimal_digit(character) => rest = &rest[len..],
                _ => return Ok(false),
            }
        }
    }
    Ok(true)
}

/// Whether the pattern match dies before it tries the string. Perl first finds where
/// `$` can hold, the end of the string or the newline that ends it, and steps back one
/// character from there over any continuation bytes: a byte that starts no character
/// there is fatal (`reghop3` of regexec.c, perl 5.34).
fn end_hop_dies(bytes: &[u8]) -> bool {
    let Some(&first) = bytes.first() else {
        return false;
    };
    // Perl looks for the newline after the string's first character, so a first byte
    // that announces a character as long as the whole string leaves it none to find.
    let newline = bytes.last() == Some(&b'\n') && announced_length(first) != bytes.len();
    let end = bytes.len() - usize::from(newline);
    let Some(mut at) = end.checked_sub(1) else {
        return false;
    };
    if bytes[at] < 0x80 {
        return false;
    }
    while at > 0 && bytes[at] & 0xc0 == 0x80 {
        at -= 1;
    }
    bytes[at] < 0xc2
}

/// The length in bytes of the character a first byte announces: `UTF8SKIP` of perl 5.34,
/// which goes on past the four bytes of UTF-8 to 13 for the byte 0xff.
fn announced_length(first: u8) -> usize {
    match first {
        0x00..=0xbf => 1,
        0xc0..=0xdf => 2,
        0xe0..=0xef => 3,
        0xf0..=0xf7 => 4,
        0xf8..=0xfb => 5,
        0xfc..=0xfd => 6,
        0xfe => 7,
        0xff => 13,
    }
}

/// The character that opens `bytes` in its shortest UTF-8 form of two to four bytes,
/// with the form's length. Every decimal digit past ASCII has such a form, and Perl
/// matches no digit in a longer one.
fn decode(bytes: &[u8]) -> Option<(usize, u32)> {
    let first = *bytes.first()?;
    let (len, min, mut character) = match first {
        0xc2..=0xdf => (2, 0x80, u32::from(first & 0x1f)),
        0xe0..=0xef => (3, 0x800, u32::from(first & 0x0f)),
        0xf0..=0xf7 => (4, 0x1_0000, u32::from(first & 0x07)),
        _ => return None,
    };
    let tail = bytes.get(1..len)?;
    for &byte in tail {
        if byte & 0xc0 != 0x80 {
            return None;
        }
        character = (character << 6) | u32::from(byte & 0x3f);
    }
    (character >= min).then_some((len, character))
}

/// The decimal digits past ASCII, as inclusive ranges: general category Nd of Unicode
/// 13.0, the tables of perl 5.34.
const DECIMAL_DIGITS: [(u32, u32); 60] = [
    (0x0660, 0x0669),
    (0x06f0, 0x06f9),
    (0x07c0, 0x07c9),
    (0x0966, 0x096f),
    (0x09e6, 0x09ef),
    (0x0a66, 0x0a6f),
    (0x0ae6, 0x0aef),
    (0x0b66, 0x0b6f),
    (0x0be6, 0x0bef),
    (0x0c66, 0x0c6f),
    (0x0ce6, 0x0cef),
    (0x0d66, 0x0d6f),
    (0x0de6, 0x0def),
    (0x0e50, 0x0e59),
    (0x0ed0, 0x0ed9),
    (0x0f20, 0x0f29),
    (0x1040, 0x1049),
    (0x1090, 0x1099),
    (0x17e0, 0x17e9),
    (0x1810, 0x1819),
    (0x1946, 0x194f),
    (0x19d0, 0x19d9),
    (0x1a80, 0x1a89),
    (0x1a90, 0x1a99),
    (0x1b50, 0x1b59),
    (0x1bb0, 0x1bb9),
    (0x1c40, 0x1c49),
    (0x1c50, 0x1c59),
    (0xa620, 0xa629),
    (0xa8d0, 0xa8d9),
    (0xa900, 0xa909),
    (0xa9d0, 0xa9d9),
    (0xa9f0, 0xa9f9),
    (0xaa50, 0xaa59),
    (0xabf0, 0xabf9),
    (0xff10, 0xff19),
    (0x1_04a0, 0x1_04a9),
    (0x1_0d30, 0x1_0d39),
    (0x1_1066, 0x1_106f),
    (0x1_10f0, 0x1_10f9),
    (0x1_1136, 0x1_113f),
    (0x1_11d0, 0x1_11d9),
    (0x1_12f0, 0x1_12f9),
    (0x1_1450, 0x1_1459),
    (0x1_14d0, 0x1_14d9),
    (0x1_1650, 0x1_1659),
    (0x1_16c0, 0x1_16c9),
    (0x1_1730, 0x1_1739),
    (0x1_18e0, 0x1_18e9),
    (0x1_1950, 0x1_1959),
    (0x1_1c50, 0x1_1c59),
    (0x1_1d50, 0x1_1d59),
    (0x1_1da0, 0x1_1da9),
    (0x1_6a60, 0x1_6a69),
    (0x1_6b50, 0x1_6b59),
    (0x1_d7ce, 0x1_d7ff),
    (0x1_e140, 0x1_e149),
    (0x1_e2f0, 0x1_e2f9),
    (0x1_e950, 0x1_e959),
    (0x1_fbf0, 0x1_fbf9),
];

fn is_decimal_digit(character: u32) -> bool {
    DECIMAL_DIGITS
        .binary_search_by(|&(first, last)| {
            if character < first {
                std::cmp::Ordering::Greater
            } else if character > last {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Equal
            }
        })
        .is_ok()
}

impl SequenceRegion {
    /// Appends the row of the chromosome directory `name`, the bytes Perl's `print`
    /// writes for `join("\t", $name, $length, $circular), "\n"` on a handle without an
    /// encoding layer, and returns whether Perl warns of a wide character for it. A
    /// length under the UTF-8 flag makes the joined string one of characters: with a
    /// digit past ASCII it is written as UTF-8, a byte above 127 of the name as the
    /// character of that number.
    pub fn write_row(&self, name: &[u8], out: &mut Vec<u8>) -> bool {
        let wide = self.utf8 && !self.length.is_ascii();
        if wide {
            out.extend_from_slice(&PerlStr::bytes(name).characters());
        } else {
            out.extend_from_slice(name);
        }
        out.push(b'\t');
        out.extend_from_slice(&self.length);
        out.push(b'\t');
        out.push(match self.circular {
            Some(true) => b'1',
            Some(false) => b'0',
            None => b'-',
        });
        out.push(b'\n');
        wide
    }
}

#[cfg(test)]
#[path = "regions_tests.rs"]
mod tests;
