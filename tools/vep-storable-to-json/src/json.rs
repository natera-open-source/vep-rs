// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! JSON text as JSON::XS 4.04 writes it under the conversion's two encoders:
//! `JSON->new->utf8->canonical` for the shards and the binding matrices, and
//! `JSON->new->canonical->pretty` for `info.json`.

/// A Perl string: bytes and the UTF-8 flag. With the flag the bytes are UTF-8 for the
/// string's characters; without it each byte is one character.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PerlStr<'a> {
    pub bytes: &'a [u8],
    pub utf8: bool,
}

impl<'a> PerlStr<'a> {
    pub fn bytes(bytes: &'a [u8]) -> Self {
        Self { bytes, utf8: false }
    }

    /// The string's characters as UTF-8: the bytes themselves under the flag or when
    /// they are ASCII, else each byte as the character of that number. Two Perl strings
    /// are equal, and order by `cmp`, as these bytes do.
    pub fn characters(&self) -> std::borrow::Cow<'a, [u8]> {
        if self.utf8 || self.bytes.is_ascii() {
            return std::borrow::Cow::Borrowed(self.bytes);
        }
        let mut out = Vec::with_capacity(self.bytes.len() + 8);
        for &byte in self.bytes {
            push_latin1(byte, &mut out);
        }
        std::borrow::Cow::Owned(out)
    }

    /// Perl's `length`: characters.
    pub fn char_len(&self) -> usize {
        if !self.utf8 {
            return self.bytes.len();
        }
        // Every byte that is not a continuation byte starts a character.
        self.bytes
            .iter()
            .filter(|byte| (**byte & 0xc0) != 0x80)
            .count()
    }
}

fn push_latin1(byte: u8, out: &mut Vec<u8>) {
    if byte < 0x80 {
        out.push(byte);
    } else {
        out.push(0xc0 | (byte >> 6));
        out.push(0x80 | (byte & 0x3f));
    }
}

/// A string JSON::XS refuses: its bytes under the UTF-8 flag are not UTF-8, or name a
/// character past U+10FFFF. The text is JSON::XS's own message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BadString(pub Vec<u8>);

/// What an escaped byte below 0x80 is written as; 0 for a byte written as itself.
const fn escape_of(byte: u8) -> u8 {
    match byte {
        b'"' => b'"',
        b'\\' => b'\\',
        0x08 => b'b',
        0x09 => b't',
        0x0a => b'n',
        0x0c => b'f',
        0x0d => b'r',
        // Every other control character is written as \u00XX.
        0x00..=0x1f => b'u',
        _ => 0,
    }
}

const HEX: &[u8; 16] = b"0123456789abcdef";

/// The length of the UTF-8 sequence that opens `bytes` as Perl's lax decoder
/// accepts it (surrogates and non-characters pass; overlong and cut sequences do not),
/// with the character's number.
fn decode_utf8(bytes: &[u8]) -> Option<(usize, u32)> {
    let first = *bytes.first()?;
    let (len, min, mut value) = match first {
        0xc2..=0xdf => (2, 0x80, u32::from(first & 0x1f)),
        0xe0..=0xef => (3, 0x800, u32::from(first & 0x0f)),
        0xf0..=0xf7 => (4, 0x1_0000, u32::from(first & 0x07)),
        _ => return None,
    };
    if bytes.len() < len {
        return None;
    }
    for &byte in &bytes[1..len] {
        if byte & 0xc0 != 0x80 {
            return None;
        }
        value = (value << 6) | u32::from(byte & 0x3f);
    }
    (value >= min).then_some((len, value))
}

/// Writes the characters of a string between the quotes of a JSON string, as the
/// `utf8` encoder does: `"` and `\` escaped, `\b \t \n \f \r` by letter, other control
/// characters as `\u00xx`, and every other character as its UTF-8.
pub fn write_string_body(string: PerlStr<'_>, out: &mut Vec<u8>) -> Result<(), BadString> {
    let bytes = string.bytes;
    out.reserve(bytes.len());
    let mut at = 0;
    while at < bytes.len() {
        // The run of bytes written as themselves.
        let run = bytes[at..]
            .iter()
            .position(|&byte| byte >= 0x80 || escape_of(byte) != 0)
            .unwrap_or(bytes.len() - at);
        out.extend_from_slice(&bytes[at..at + run]);
        at += run;
        let Some(&byte) = bytes.get(at) else {
            break;
        };
        if byte < 0x80 {
            out.push(b'\\');
            match escape_of(byte) {
                b'u' => {
                    out.extend_from_slice(b"u00");
                    out.push(HEX[(byte >> 4) as usize]);
                    out.push(HEX[(byte & 15) as usize]);
                }
                letter => out.push(letter),
            }
            at += 1;
        } else if string.utf8 {
            let Some((len, value)) = decode_utf8(&bytes[at..]) else {
                let shown = &bytes[at..bytes.len().min(at + 11)];
                let mut text = b"malformed or illegal unicode character in string [".to_vec();
                text.extend_from_slice(shown);
                text.extend_from_slice(b"], cannot convert to JSON");
                return Err(BadString(text));
            };
            if value >= 0x11_0000 {
                return Err(BadString(
                    format!(
                        "out of range codepoint (0x{value:x}) encountered, unrepresentable in JSON"
                    )
                    .into_bytes(),
                ));
            }
            out.extend_from_slice(&bytes[at..at + len]);
            at += len;
        } else {
            push_latin1(byte, out);
            at += 1;
        }
    }
    Ok(())
}

/// A whole JSON string: the body between quotes.
pub fn write_string(string: PerlStr<'_>, out: &mut Vec<u8>) -> Result<(), BadString> {
    out.push(b'"');
    write_string_body(string, out)?;
    out.push(b'"');
    Ok(())
}

/// A JSON string of bytes that need no escape: a key or a value the conversion makes.
pub fn write_plain(text: &[u8], out: &mut Vec<u8>) {
    debug_assert!(text.iter().all(|&byte| byte < 0x80 && escape_of(byte) == 0));
    out.push(b'"');
    out.extend_from_slice(text);
    out.push(b'"');
}

/// The members of one JSON object whose keys the conversion writes in sorted order.
/// `key` writes the separator and the key; the caller writes the value after it.
pub struct Object<'a> {
    out: &'a mut Vec<u8>,
    first: bool,
    #[cfg(debug_assertions)]
    last: &'static str,
}

impl<'a> Object<'a> {
    pub fn open(out: &'a mut Vec<u8>) -> Self {
        out.push(b'{');
        Self {
            out,
            first: true,
            #[cfg(debug_assertions)]
            last: "",
        }
    }

    /// Starts the member `name` and returns the buffer for its value. The names of one
    /// object come in byte order, the order `canonical` sorts them into.
    pub fn key(&mut self, name: &'static str) -> &mut Vec<u8> {
        #[cfg(debug_assertions)]
        {
            debug_assert!(
                self.first || self.last < name,
                "{} before {name}",
                self.last
            );
            self.last = name;
        }
        if !self.first {
            self.out.push(b',');
        }
        self.first = false;
        write_plain(name.as_bytes(), self.out);
        self.out.push(b':');
        self.out
    }

    pub fn close(self) {
        self.out.push(b'}');
    }
}

/// A value of the `info.json` document.
#[derive(Debug, Clone, PartialEq)]
pub enum Pretty {
    /// A string held as Perl characters in UTF-8 (see [`PerlStr::characters`]).
    String(Vec<u8>),
    Integer(i128),
    Double(f64),
    True,
    Array(Vec<Pretty>),
    /// Members by key, the key as Perl characters in UTF-8.
    Object(std::collections::BTreeMap<Vec<u8>, Pretty>),
}

/// The document as `JSON->new->canonical->pretty->encode` returns it: three spaces a
/// level, ` : ` between key and value, one member or item a line, an empty container
/// as `[]` or `{}`, and a newline at the end. The result is the encoder's
/// character string in UTF-8.
pub fn write_pretty(value: &Pretty) -> Result<Vec<u8>, BadString> {
    let mut out = Vec::new();
    write_pretty_value(value, 0, &mut out)?;
    out.push(b'\n');
    Ok(out)
}

fn indent(level: usize, out: &mut Vec<u8>) {
    out.resize(out.len() + level * 3, b' ');
}

fn characters(text: &[u8]) -> PerlStr<'_> {
    PerlStr {
        bytes: text,
        utf8: true,
    }
}

fn write_pretty_value(value: &Pretty, level: usize, out: &mut Vec<u8>) -> Result<(), BadString> {
    match value {
        Pretty::String(text) => write_string(characters(text), out)?,
        Pretty::Integer(number) => crate::number::write_int(*number, out),
        Pretty::Double(number) => crate::number::write_g15(*number, out),
        Pretty::True => out.extend_from_slice(b"true"),
        Pretty::Array(items) => {
            out.push(b'[');
            if !items.is_empty() {
                out.push(b'\n');
                for (index, item) in items.iter().enumerate() {
                    indent(level + 1, out);
                    write_pretty_value(item, level + 1, out)?;
                    if index + 1 < items.len() {
                        out.push(b',');
                    }
                    out.push(b'\n');
                }
                indent(level, out);
            }
            out.push(b']');
        }
        Pretty::Object(members) => {
            out.push(b'{');
            if !members.is_empty() {
                out.push(b'\n');
                for (index, (key, member)) in members.iter().enumerate() {
                    indent(level + 1, out);
                    write_string(characters(key), out)?;
                    out.extend_from_slice(b" : ");
                    write_pretty_value(member, level + 1, out)?;
                    if index + 1 < members.len() {
                        out.push(b',');
                    }
                    out.push(b'\n');
                }
                indent(level, out);
            }
            out.push(b'}');
        }
    }
    Ok(())
}

/// The bytes Perl's `print` writes for a character string held as UTF-8, on a handle
/// with no encoding layer: one byte a character when every character is below 256,
/// else the UTF-8 itself.
pub fn print_bytes(characters: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(characters.len());
    let mut at = 0;
    while at < characters.len() {
        let byte = characters[at];
        if byte < 0x80 {
            out.push(byte);
            at += 1;
        } else if (byte == 0xc2 || byte == 0xc3) && at + 1 < characters.len() {
            out.push(((byte & 0x03) << 6) | (characters[at + 1] & 0x3f));
            at += 2;
        } else {
            // A character of 256 or more: Perl writes the whole string as UTF-8.
            return characters.to_vec();
        }
    }
    out
}

/// Base64 with padding and no line breaks: `MIME::Base64::encode_base64($bytes, '')`.
pub fn write_base64(bytes: &[u8], out: &mut Vec<u8>) {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    out.reserve(bytes.len().div_ceil(3) * 4);
    let mut chunks = bytes.chunks_exact(3);
    for chunk in &mut chunks {
        let group = (u32::from(chunk[0]) << 16) | (u32::from(chunk[1]) << 8) | u32::from(chunk[2]);
        out.extend_from_slice(&[
            ALPHABET[(group >> 18) as usize],
            ALPHABET[((group >> 12) & 63) as usize],
            ALPHABET[((group >> 6) & 63) as usize],
            ALPHABET[(group & 63) as usize],
        ]);
    }
    match chunks.remainder() {
        [one] => {
            let group = u32::from(*one) << 16;
            out.extend_from_slice(&[
                ALPHABET[(group >> 18) as usize],
                ALPHABET[((group >> 12) & 63) as usize],
                b'=',
                b'=',
            ]);
        }
        [one, two] => {
            let group = (u32::from(*one) << 16) | (u32::from(*two) << 8);
            out.extend_from_slice(&[
                ALPHABET[(group >> 18) as usize],
                ALPHABET[((group >> 12) & 63) as usize],
                ALPHABET[((group >> 6) & 63) as usize],
                b'=',
            ]);
        }
        _ => {}
    }
}

#[cfg(test)]
#[path = "json_tests.rs"]
mod tests;
