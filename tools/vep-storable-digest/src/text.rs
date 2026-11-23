// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Perl strings a program builds and prints: joined with the UTF-8 flag of their
//! parts, written as `print` writes them on a handle without an encoding layer.

use crate::json::{print_bytes, PerlStr};

/// A Perl string: bytes and the UTF-8 flag.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Text {
    pub bytes: Vec<u8>,
    pub utf8: bool,
}

impl Text {
    pub fn bytes(bytes: impl Into<Vec<u8>>) -> Self {
        Text {
            bytes: bytes.into(),
            utf8: false,
        }
    }

    /// Appends a part as Perl's `.` does: a string with the flag makes the whole a
    /// string of characters.
    pub fn push(&mut self, bytes: &[u8], utf8: bool) {
        if self.utf8 == utf8 || bytes.is_ascii() {
            // The flag stays what it was, unless only this part carries it.
            if utf8 && !self.utf8 && !self.bytes.is_ascii() {
                self.bytes = PerlStr::bytes(&self.bytes).characters().into_owned();
            }
            self.utf8 |= utf8;
            self.bytes.extend_from_slice(bytes);
        } else if utf8 {
            self.bytes = PerlStr::bytes(&self.bytes).characters().into_owned();
            self.utf8 = true;
            self.bytes.extend_from_slice(bytes);
        } else {
            self.bytes
                .extend_from_slice(&PerlStr::bytes(bytes).characters());
        }
    }

    /// The string's characters as UTF-8: equal for equal Perl strings, so the key a
    /// Perl hash holds it under.
    pub fn characters(&self) -> Vec<u8> {
        PerlStr {
            bytes: &self.bytes,
            utf8: self.utf8,
        }
        .characters()
        .into_owned()
    }

    /// The bytes `print` writes: the string's own without the flag; with it, one byte a
    /// character when every character is below 256, else the UTF-8.
    pub fn printed(&self) -> Vec<u8> {
        if self.utf8 {
            print_bytes(&self.bytes)
        } else {
            self.bytes.clone()
        }
    }
}

/// `s/\.\d+$//` on an identifier: a final `.<digits>` goes, also before one final
/// newline. The digits are ASCII digits.
pub fn strip_version(text: &mut Vec<u8>) {
    let newline = text.last() == Some(&b'\n');
    let end = text.len() - usize::from(newline);
    let digits = text[..end]
        .iter()
        .rev()
        .take_while(|byte| byte.is_ascii_digit())
        .count();
    if digits == 0 || digits == end || text[end - digits - 1] != b'.' {
        return;
    }
    text.drain(end - digits - 1..end);
}

#[cfg(test)]
mod tests {
    use super::{strip_version, Text};

    fn stripped(text: &str) -> String {
        let mut bytes = text.as_bytes().to_vec();
        strip_version(&mut bytes);
        String::from_utf8(bytes).expect("ASCII")
    }

    #[test]
    fn a_final_version_goes_and_nothing_else() {
        assert_eq!(stripped("ENST00000335137.4"), "ENST00000335137");
        assert_eq!(stripped("ENST00000335137"), "ENST00000335137");
        assert_eq!(stripped("a.1.22"), "a.1");
        assert_eq!(stripped("a.b"), "a.b");
        assert_eq!(stripped("a."), "a.");
        assert_eq!(stripped(".7"), "");
        assert_eq!(stripped("7"), "7");
        assert_eq!(stripped(""), "");
        // Perl's `$` also stands before one final newline.
        assert_eq!(stripped("a.3\n"), "a\n");
        assert_eq!(stripped("a.3\n\n"), "a.3\n\n");
    }

    #[test]
    fn a_part_with_the_flag_makes_the_whole_a_string_of_characters() {
        let mut text = Text::bytes(vec![b'a', 0xe9]);
        text.push(b"\t", false);
        assert_eq!(text.printed(), vec![b'a', 0xe9, b'\t']);
        // U+00E9 under the flag: the whole prints one byte a character.
        text.push(&[0xc3, 0xa9], true);
        assert!(text.utf8);
        assert_eq!(text.printed(), vec![b'a', 0xe9, b'\t', 0xe9]);
        // U+0100 has no byte: the whole prints as UTF-8.
        text.push(&[0xc4, 0x80], true);
        assert_eq!(
            text.printed(),
            vec![b'a', 0xc3, 0xa9, b'\t', 0xc3, 0xa9, 0xc4, 0x80]
        );
        // A byte above 127 joins a string of characters as that character.
        text.push(&[0xe9], false);
        assert_eq!(text.characters().last(), Some(&0xa9));
    }
}
