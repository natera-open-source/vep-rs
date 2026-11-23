// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A dump of `vep-storable-dump` as the classification reads it: the blocks by
//! transcript id, each a list of lines with number quoting and two keys neutralised,
//! and what a line names.

use std::collections::HashMap;

use crate::py::{is_decimal, universal_newlines};

/// The marker that opens a block.
const MARKER: &str = "### ";

/// `str.isspace()` of one character: the characters `str.strip()` takes off.
fn is_space(c: char) -> bool {
    matches!(
        c,
        '\u{9}'..='\u{d}'
            | '\u{1c}'..='\u{20}'
            | '\u{85}'
            | '\u{a0}'
            | '\u{1680}'
            | '\u{2000}'..='\u{200a}'
            | '\u{2028}'
            | '\u{2029}'
            | '\u{202f}'
            | '\u{205f}'
            | '\u{3000}'
    )
}

/// `re.sub(r"=> '(-?\d+)'(,?)$", r"=> \1\2", line)`: an integer in quotes at the end of
/// a line loses its quotes. `\d` is a decimal digit of any script.
pub fn unquote_integer(line: &str) -> Option<String> {
    let (body, comma) = match line.strip_suffix(',') {
        Some(body) => (body, ","),
        None => (line, ""),
    };
    let inner = body.strip_suffix('\'')?;
    let digits_from = inner
        .char_indices()
        .rev()
        .take_while(|(_, c)| is_decimal(*c))
        .last()
        .map(|(at, _)| at)?;
    let number_from = if inner[..digits_from].ends_with('-') {
        digits_from - 1
    } else {
        digits_from
    };
    let head = inner[..number_from].strip_suffix("=> '")?;
    Some(format!("{head}=> {}{comma}", &inner[number_from..]))
}

/// The lines of one block as the classification compares them.
fn block_lines(body: &str) -> Vec<String> {
    body.split('\n')
        .map(|line| unquote_integer(line).unwrap_or_else(|| line.to_string()))
        .filter(|line| !line.contains("'dbID' =>") && !line.contains("'_gene_version' =>"))
        .collect()
}

/// A block that holds no line feed: the marker's line never ends.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NoLineEnd;

/// The blocks of a dump by transcript id. The text is the file decoded as UTF-8 with
/// `errors="replace"`; a block runs from `### ` to the next `### `, its first line is
/// the id without the white space around it, and a later block of an id replaces an
/// earlier one.
pub fn load(text: &str) -> Result<HashMap<String, Vec<String>>, NoLineEnd> {
    let text = universal_newlines(text);
    let mut blocks = HashMap::new();
    for chunk in text.split(MARKER).skip(1) {
        let (id, body) = chunk.split_once('\n').ok_or(NoLineEnd)?;
        blocks.insert(id.trim_matches(is_space).to_string(), block_lines(body));
    }
    Ok(blocks)
}

/// `re.search(r"'peptide' => '([A-Z*]*)'", line)` on the first line that matches: the
/// stored translation.
pub fn peptide_of(lines: &[String]) -> Option<&str> {
    const OPEN: &str = "'peptide' => '";
    lines.iter().find_map(|line| {
        line.match_indices(OPEN).find_map(|(at, _)| {
            let rest = &line[at + OPEN.len()..];
            let end = rest
                .find(|c: char| !(c.is_ascii_uppercase() || c == '*'))
                .unwrap_or(rest.len());
            rest[end..].starts_with('\'').then(|| &rest[..end])
        })
    })
}

/// `re.search(r"'([A-Za-z_0-9]+)' =>", line)`: the first key a line names.
pub fn key_of(line: &str) -> Option<&str> {
    let bytes = line.as_bytes();
    let mut at = 0;
    while at < bytes.len() {
        if bytes[at] == b'\'' {
            let word = bytes[at + 1..]
                .iter()
                .take_while(|byte| byte.is_ascii_alphanumeric() || **byte == b'_')
                .count();
            if word > 0 && bytes[at + 1 + word..].starts_with(b"' =>") {
                return Some(&line[at + 1..at + 1 + word]);
            }
        }
        at += 1;
    }
    None
}

#[cfg(test)]
#[path = "dumps_tests.rs"]
mod tests;
