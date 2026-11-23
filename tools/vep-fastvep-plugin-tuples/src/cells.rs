// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One populated supplementary cell: its entries and parts, the value each rule of the
//! column reads, percent-decoded, and the parts no rule reads.

use crate::columns::{Column, Reading};
use crate::map::Half;

/// Appends `text` percent-decoded: `%XX` with two hexadecimal digits is its byte, and any
/// other `%` stays as written, as `urllib.parse.unquote_to_bytes` reads it.
fn unquote_into(out: &mut Vec<u8>, text: &[u8]) {
    let hex = |byte: u8| char::from(byte).to_digit(16).map(|digit| digit as u8);
    let mut at = 0;
    while at < text.len() {
        if text[at] == b'%' && at + 2 < text.len() {
            if let (Some(high), Some(low)) = (hex(text[at + 1]), hex(text[at + 2])) {
                out.push(high << 4 | low);
                at += 3;
                continue;
            }
        }
        out.push(text[at]);
        at += 1;
    }
}

/// Whether a byte may open a character that is written percent-encoded in a value.
fn may_break(byte: u8) -> bool {
    byte < 0x20 || byte == 0xc2 || byte == 0xe2
}

/// Appends `text` with the tab and every character `str.splitlines` ends a line at
/// written percent-encoded over its UTF-8 bytes, as `urllib.parse.quote` writes it.
fn encode_breaks_into(out: &mut Vec<u8>, text: &[u8]) {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    let mut at = 0;
    while at < text.len() {
        let byte = text[at];
        let encoded = match byte {
            0x09..=0x0d | 0x1c..=0x1e => 1,
            0xc2 if text.get(at + 1) == Some(&0x85) => 2,
            0xe2 if text.get(at + 1) == Some(&0x80)
                && matches!(text.get(at + 2), Some(0xa8 | 0xa9)) =>
            {
                3
            }
            _ => 0,
        };
        if encoded == 0 {
            out.push(byte);
            at += 1;
            continue;
        }
        for &byte in &text[at..at + encoded] {
            out.extend_from_slice(&[
                b'%',
                HEX[usize::from(byte >> 4)],
                HEX[usize::from(byte & 0xf)],
            ]);
        }
        at += encoded;
    }
}

/// Appends one value percent-decoded, a tab and every line-breaking character left
/// encoded. The escapes are decoded as bytes among the value's own bytes, so the two read
/// as one UTF-8 text and a byte that is not UTF-8 stays that byte.
pub fn decode_into(out: &mut Vec<u8>, text: &[u8], scratch: &mut Vec<u8>) {
    let escaped = memchr::memchr(b'%', text).is_some();
    if !escaped && !text.iter().copied().any(may_break) {
        out.extend_from_slice(text);
        return;
    }
    if !escaped {
        encode_breaks_into(out, text);
        return;
    }
    scratch.clear();
    unquote_into(scratch, text);
    encode_breaks_into(out, scratch);
}

/// The two halves of a part written `word(score)`: the text before its last `(` and the
/// text between that and the `)` it ends with. None for a part of another shape.
pub fn word_score(part: &[u8]) -> Option<(&[u8], &[u8])> {
    let inner = part.strip_suffix(b")")?;
    let open = memchr::memrchr(b'(', inner)?;
    Some((&inner[..open], &inner[open + 1..]))
}

/// What one populated cell of a mapped column holds.
#[derive(Debug, Default)]
pub struct Decoded {
    /// The values of the kept fields in the order the rows write them: the place of the
    /// field among the kept fields and the end of a tuple line, `field<TAB>value<LF>`.
    pub pairs: Vec<(usize, Vec<u8>)>,
    /// The populated parts and `word(score)` halves no rule reads.
    pub dropped: Vec<(u64, Half)>,
    /// The populated values `--fields` leaves out.
    pub deselected: u64,
    /// The parts a `word(score)` rule could not split.
    pub unshaped: u64,
    /// Whether the cell holds several entries.
    pub several: bool,
}

/// Buffers one worker decodes its cells with.
#[derive(Default)]
pub struct Scratch {
    unquoted: Vec<u8>,
    unshaped: Vec<(usize, u64)>,
}

/// Decodes one populated cell of `column` into `decoded`, which is emptied first.
pub fn decode_cell(cell: &[u8], column: &Column, decoded: &mut Decoded, scratch: &mut Scratch) {
    decoded.pairs.clear();
    decoded.dropped.clear();
    decoded.deselected = 0;
    scratch.unshaped.clear();
    let mut entries = 0usize;
    let mut parts: Vec<&[u8]> = Vec::new();
    for (number, entry) in cell.split(|b| *b == b',').enumerate() {
        entries += 1;
        parts.clear();
        parts.extend(entry.split(|b| *b == b'|'));
        let count = parts.len() as u64;
        for (rule, slot) in column.rules.iter().zip(&column.slots) {
            if rule.first >= count {
                continue;
            }
            let chunk =
                &parts[rule.first as usize..rule.last.saturating_add(1).min(count) as usize];
            if chunk.iter().all(|part| part.is_empty()) {
                continue;
            }
            // The half of a `word(score)` part the rule reads, where it reads one.
            let half = if rule.half == Half::Whole {
                None
            } else {
                let Some((word, score)) = word_score(chunk[0]) else {
                    if !scratch.unshaped.contains(&(number, rule.first)) {
                        scratch.unshaped.push((number, rule.first));
                    }
                    continue;
                };
                let half = if rule.half == Half::Word { word } else { score };
                if half.is_empty() {
                    continue;
                }
                Some(half)
            };
            let Some(slot) = *slot else {
                decoded.deselected += 1;
                continue;
            };
            let mut tail = Vec::with_capacity(rule.field.len() + 2 + chunk[0].len());
            tail.extend_from_slice(rule.field.as_bytes());
            tail.push(b'\t');
            if let Some(half) = half {
                decode_into(&mut tail, half, &mut scratch.unquoted);
            } else if rule.last > rule.first {
                for (at, part) in chunk.iter().enumerate() {
                    if at > 0 {
                        tail.push(b'|');
                    }
                    decode_into(&mut tail, part, &mut scratch.unquoted);
                }
            } else if !rule.separator.is_empty() {
                for (at, item) in chunk[0].split(|b| *b == b'&').enumerate() {
                    if at > 0 {
                        tail.extend_from_slice(rule.separator.as_bytes());
                    }
                    decode_into(&mut tail, item, &mut scratch.unquoted);
                }
            } else {
                decode_into(&mut tail, chunk[0], &mut scratch.unquoted);
            }
            tail.push(b'\n');
            let pair = (slot, tail);
            if !decoded.pairs.contains(&pair) {
                decoded.pairs.push(pair);
            }
        }
        for (index, part) in parts.iter().enumerate().skip(1) {
            if part.is_empty() {
                continue;
            }
            let index = index as u64;
            match column.reading(index) {
                Reading::Whole => {}
                Reading::Unread => decoded.dropped.push((index, Half::Whole)),
                Reading::Halves { word, score } => {
                    if let Some((word_text, score_text)) = word_score(part) {
                        if !word && !word_text.is_empty() {
                            decoded.dropped.push((index, Half::Word));
                        }
                        if !score && !score_text.is_empty() {
                            decoded.dropped.push((index, Half::Score));
                        }
                    }
                }
            }
        }
    }
    decoded.unshaped = scratch.unshaped.len() as u64;
    decoded.several = entries > 1;
}

#[cfg(test)]
#[path = "cells_decoding_tests.rs"]
mod decoding_tests;
