// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Identifiers under restricted redistribution terms in the text a corpus carries: the two
//! shapes the tool reads, and the fabricated identifier of the same shape (the letters, then
//! a zero for each digit) that it writes in the place of each. `stored` reads the files of
//! a corpus for them, and `run` is the `identifiers` subcommand.

mod command;
mod gzip;
pub mod stored;

pub use command::run;

use memchr::memchr3_iter;
use std::borrow::Cow;
use std::collections::BTreeSet;
use std::ops::Range;

/// The letter pairs and digit counts of the second shape of identifier under restricted
/// redistribution terms: every one the variation files of Ensembl VEP's release 116 human
/// caches give a name of that source, on GRCh37 and GRCh38. A row is a first letter, the
/// letters that follow it and the counts of digits a name of such a pair takes.
pub const SECOND_SHAPE_BY_FIRST_LETTER: [(u8, &[u8], &[usize]); 5] = [
    (b'C', b"DGIMRSX", &[6, 7]),
    (b'C', b"P", &[5, 6]),
    (b'B', b"M", &[7]),
    (b'H', b"DIMRS", &[4, 6]),
    (b'H', b"X", &[6]),
];

/// The column of an output table that names the known variants at a record's position.
pub const EXISTING_VARIATION: &str = "Existing_variation";

/// The read of the two shapes of identifier under restricted redistribution terms. The
/// first is `COS`, a type letter and digits, read wherever it stands as a word. The second
/// is a letter pair and digits of `SECOND_SHAPE_BY_FIRST_LETTER`, which a public sequence
/// accession can share and so is read only as a whole member of a list (after `&`, `,`, `|`
/// or `;`, or opening a field before `&`, `,` or `|`) or as a JSON `id`.
/// `lone_cell_identifiers` reads the second shape in the one table cell that holds no
/// accession.
#[derive(Clone, Copy, Debug)]
pub struct RestrictedIdentifier;

/// The reader of both shapes.
pub const RESTRICTED_IDENTIFIER: RestrictedIdentifier = RestrictedIdentifier;

impl RestrictedIdentifier {
    /// The place of every identifier `text` holds in a position the shapes are read in, the
    /// fabricated ones among them, in the order of the text.
    pub fn find_iter(self, text: &str) -> Vec<Range<usize>> {
        let bytes = text.as_bytes();
        let mut found = Vec::new();
        let mut from = 0;
        for at in memchr3_iter(b'C', b'B', b'H', bytes) {
            if at < from {
                continue;
            }
            if let Some(end) = first_shape_end(bytes, at).or_else(|| second_shape_end(bytes, at)) {
                found.push(at..end);
                from = end;
            }
        }
        found
    }
}

/// A letter, a digit or an underscore of ASCII.
fn is_word(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || byte == b'_'
}

/// Whether the two characters before `at` are a tab or a line feed as a JSON string writes
/// it, `\t` or `\n`: the word after them starts at `at`, though a letter stands before it.
fn after_escape(bytes: &[u8], at: usize) -> bool {
    at >= 2 && bytes[at - 2] == b'\\' && matches!(bytes[at - 1], b't' | b'n')
}

/// The end of an identifier of the first shape that starts at `at` as a word.
fn first_shape_end(bytes: &[u8], at: usize) -> Option<usize> {
    let rest = &bytes[at..];
    if !rest.starts_with(b"COS") || !rest.get(3).is_some_and(u8::is_ascii_uppercase) {
        return None;
    }
    let digits = rest[4..]
        .iter()
        .take_while(|byte| byte.is_ascii_digit())
        .count();
    let end = at + 4 + digits;
    let opens = at == 0 || !is_word(bytes[at - 1]) || after_escape(bytes, at);
    let closes = !bytes.get(end).is_some_and(|byte| is_word(*byte));
    (digits > 0 && opens && closes).then_some(end)
}

/// The length of the second shape where `bytes` holds it at `at`: its letters and the whole
/// run of digits after them.
fn second_shape_len(bytes: &[u8], at: usize) -> Option<usize> {
    let rest = &bytes[at..];
    let (first, second) = (*rest.first()?, *rest.get(1)?);
    let pair = |(letter, seconds, _): &&(u8, &[u8], &[usize])| {
        *letter == first && seconds.contains(&second)
    };
    let (_, _, counts) = SECOND_SHAPE_BY_FIRST_LETTER.iter().find(pair)?;
    let digits = rest[2..]
        .iter()
        .take_while(|byte| byte.is_ascii_digit())
        .count();
    counts.contains(&digits).then_some(2 + digits)
}

/// The end of an identifier of the second shape that starts at `at` in one of the three
/// positions `RESTRICTED_IDENTIFIER` reads the shape in.
fn second_shape_end(bytes: &[u8], at: usize) -> Option<usize> {
    let end = at + second_shape_len(bytes, at)?;
    let before = at.checked_sub(1).map(|place| bytes[place]);
    let after = bytes.get(end).copied();
    // A member after a separator, up to the next separator or the end of its field or line.
    let member = before.is_some_and(|byte| b"&,|;".contains(&byte))
        && after.is_none_or(|byte| b"\n&,|;\t\r\"\\".contains(&byte));
    // The member that opens a line or a field, before a separator.
    let opens = (before.is_none_or(|byte| b"\n=\t\"".contains(&byte)) || after_escape(bytes, at))
        && after.is_some_and(|byte| b"&,|".contains(&byte));
    let head = &bytes[..at];
    let id = (head.ends_with(b"\"id\":\"") || head.ends_with(b"\"id\": \"")) && after == Some(b'"');
    (member || opens || id).then_some(end)
}

/// Whether `data` holds what every identifier the bounded reads find contains: the first
/// shape's letters and a digit, or the second shape where no word character and no version
/// follows it. Data that holds neither needs no bounded read.
pub fn holds_identifier_bytes(data: &[u8]) -> bool {
    memchr3_iter(b'C', b'B', b'H', data).any(|at| {
        let rest = &data[at..];
        let first = rest.starts_with(b"COS")
            && rest.get(3).is_some_and(u8::is_ascii_uppercase)
            && rest.get(4).is_some_and(u8::is_ascii_digit);
        first
            || second_shape_len(data, at).is_some_and(|length| {
                let after = &rest[length..];
                let word = after.first().is_some_and(|byte| is_word(*byte));
                let version =
                    after.first() == Some(&b'.') && after.get(1).is_some_and(u8::is_ascii_digit);
                !word && !version
            })
    })
}

/// The count of capital letters an identifier starts with.
fn letters(identifier: &str) -> usize {
    let rest = identifier.trim_start_matches(|c: char| c.is_ascii_uppercase());
    identifier.len() - rest.len()
}

/// Whether an identifier of either restricted shape is the fabricated form: every digit a
/// zero.
pub fn is_fabricated(identifier: &str) -> bool {
    identifier[letters(identifier)..]
        .trim_matches('0')
        .is_empty()
}

/// The second-shape identifiers that are the whole `Existing_variation` cell of a row of an
/// output table in `text` and are not the fabricated form. A table is the rows below a
/// header line, one `#` and tab-separated names, that names the column; a header line that
/// does not name it ends the table. That column names known variants and never a sequence,
/// so the shape alone is read there, where in any other cell it could be a sequence
/// accession's.
pub fn lone_cell_identifiers(text: &str) -> BTreeSet<String> {
    let mut found = BTreeSet::new();
    if !text.contains(EXISTING_VARIATION) {
        return found;
    }
    let mut column = None;
    for line in text.split('\n') {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if let Some(names) = line.strip_prefix('#') {
            if !names.starts_with('#') {
                column = names
                    .split('\t')
                    .position(|name| name == EXISTING_VARIATION);
            }
            continue;
        }
        let Some(cell) = column.and_then(|column| line.split('\t').nth(column)) else {
            continue;
        };
        if second_shape_len(cell.as_bytes(), 0) == Some(cell.len()) && !is_fabricated(cell) {
            found.insert(cell.to_string());
        }
    }
    found
}

/// The identifiers `text` holds in another form than the fabricated one: those
/// `RESTRICTED_IDENTIFIER` reads and those of `lone_cell_identifiers`.
pub fn restricted_identifiers(text: &str) -> BTreeSet<String> {
    let read = RESTRICTED_IDENTIFIER.find_iter(text);
    let read = read.into_iter().map(|place| &text[place]);
    let mut found = lone_cell_identifiers(text);
    found.extend(
        read.filter(|identifier| !is_fabricated(identifier))
            .map(str::to_string),
    );
    found
}

/// The places where one of `identifiers` stands as a word in `text`, a table cell or a JSON
/// string included, though not after a colon or before a version. An identifier read once
/// in a list is the same identifier alone in a cell, where its shape by itself would be a
/// sequence accession's.
fn read_elsewhere(text: &str, identifiers: &BTreeSet<&str>) -> Vec<Range<usize>> {
    let bytes = text.as_bytes();
    let firsts: Vec<u8> = identifiers
        .iter()
        .filter_map(|i| i.bytes().next())
        .collect();
    let mut found = Vec::new();
    let mut from = 0;
    for (at, byte) in bytes.iter().enumerate() {
        if at < from || !firsts.contains(byte) {
            continue;
        }
        let opens = at == 0
            || !(is_word(bytes[at - 1]) || bytes[at - 1] == b':')
            || after_escape(bytes, at);
        if !opens {
            continue;
        }
        let word = bytes[at..]
            .iter()
            .take_while(|byte| is_word(**byte))
            .count();
        let end = at + word;
        let versioned =
            bytes.get(end) == Some(&b'.') && bytes.get(end + 1).is_some_and(u8::is_ascii_digit);
        if !versioned && identifiers.contains(&text[at..end]) {
            found.push(at..end);
            from = end;
        }
    }
    found
}

/// `text` with the identifier at each of `places` in its fabricated form, and how many of
/// them that changed.
fn fabricate<'t>(text: &'t str, places: &[Range<usize>]) -> (Cow<'t, str>, usize) {
    let changed = |place: &&Range<usize>| !is_fabricated(&text[(*place).clone()]);
    let count = places.iter().filter(changed).count();
    if count == 0 {
        return (Cow::Borrowed(text), 0);
    }
    let mut made = String::with_capacity(text.len());
    let mut from = 0;
    for place in places {
        let kept = place.start + letters(&text[place.clone()]);
        made.push_str(&text[from..kept]);
        made.extend(std::iter::repeat_n('0', place.end - kept));
        from = place.end;
    }
    made.push_str(&text[from..]);
    (Cow::Owned(made), count)
}

/// `text` with every identifier under restricted redistribution terms replaced by the
/// fabricated one of the same shape (its letters, then a zero for each of its digits), and
/// how many were replaced. The identifiers are those `RESTRICTED_IDENTIFIER` reads, and each
/// second-shape one that `restricted_identifiers` reads in `text` or that `also` names,
/// wherever that one stands as a word. `also` holds identifiers `restricted_identifiers`
/// gave for another text.
pub fn replace_restricted_identifiers<'t>(
    text: &'t str,
    also: &BTreeSet<String>,
) -> (Cow<'t, str>, usize) {
    let read = RESTRICTED_IDENTIFIER.find_iter(text);
    let alone = lone_cell_identifiers(text);
    let mut named: BTreeSet<&str> = also.iter().chain(&alone).map(String::as_str).collect();
    named.extend(read.iter().map(|place| &text[place.clone()]));
    named.retain(|identifier| !identifier.starts_with("COS") && !is_fabricated(identifier));
    let (listed, count) = fabricate(text, &read);
    if named.is_empty() {
        return (listed, count);
    }
    let elsewhere = read_elsewhere(&listed, &named);
    let (whole, more) = fabricate(&listed, &elsewhere);
    match whole {
        Cow::Owned(whole) => (Cow::Owned(whole), count + more),
        Cow::Borrowed(_) => (listed, count),
    }
}

/// `text` as `replace_restricted_identifiers` rewrites it. Two identifiers of one prefix
/// and length become the same one, and text that holds only fabricated ones comes back
/// unchanged.
pub fn fabricated_identifiers<'t>(text: &'t str, also: &BTreeSet<String>) -> Cow<'t, str> {
    replace_restricted_identifiers(text, also).0
}

/// How many identifiers under restricted redistribution terms `text` holds in another form
/// than the fabricated one, counted as `replace_restricted_identifiers` replaces them.
pub fn restricted_identifier_count(text: &str, also: &BTreeSet<String>) -> usize {
    replace_restricted_identifiers(text, also).1
}

#[cfg(test)]
#[path = "identifiers_tests.rs"]
mod tests;
