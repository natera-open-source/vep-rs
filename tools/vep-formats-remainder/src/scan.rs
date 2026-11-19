// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One pass over a side's data lines on every worker: each line's Location as a hash, the
//! lines whose group the comparator's lists or the line's own bytes hand over, and the
//! shapes under which every line is handed over.

use crate::input::{Line, Side, HEADER};
use crate::layout::{Layout, Whole, KEY_COLUMNS};
use crate::Args;
use memchr::{memchr, memchr2, memchr2_iter, memchr_iter, memmem};
use std::collections::HashSet;
use std::fs;
use std::path::Path;

/// What the comparator hands over at run time: the Locations, texts and column values
/// under which a group is handed over whole.
#[derive(Default)]
pub struct Force {
    pub locations: HashSet<Vec<u8>>,
    pub texts: Vec<Vec<u8>>,
    pub keep: Option<(Vec<u8>, HashSet<Vec<u8>>)>,
}

/// The lines of a list file, an empty line being the empty value.
fn list_lines(path: &Path) -> Result<Vec<Vec<u8>>, String> {
    let bytes = fs::read(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut lines: Vec<Vec<u8>> = bytes.split(|&b| b == b'\n').map(<[u8]>::to_vec).collect();
    if lines.last().is_some_and(Vec::is_empty) {
        lines.pop();
    }
    Ok(lines)
}

impl Force {
    pub fn load(args: &Args) -> Result<Force, String> {
        let mut force = Force::default();
        if let Some(path) = &args.force_locations {
            force.locations = list_lines(path)?.into_iter().collect();
        }
        if let Some(path) = &args.force_text {
            force.texts = list_lines(path)?;
            if force.texts.iter().any(Vec::is_empty) {
                return Err(format!(
                    "{}: an empty text matches every line",
                    path.display()
                ));
            }
        }
        if let Some((column, path)) = &args.keep_values {
            force.keep = Some((column.clone(), list_lines(path)?.into_iter().collect()));
        }
        Ok(force)
    }
}

/// The scan of one side: per line, the hash of its Location (0 on a header line) and
/// whether its group is handed over; and the bit set of `Whole` reasons met.
pub struct Scan {
    pub location: Vec<u64>,
    pub forced: Vec<bool>,
    pub whole: u32,
}

/// FNV-1a over a Location's bytes. Two Locations of one hash form one group, which is kept
/// back only when the lines of both are the same on the two sides.
fn hash(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h = (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

struct LineScanner<'a> {
    layout: &'a Layout<'a>,
    force: &'a Force,
    reference: bool,
    keep_column: Option<usize>,
    text_finders: Vec<memmem::Finder<'a>>,
    last_location: &'a [u8],
    last_forced: bool,
}

impl<'a> LineScanner<'a> {
    fn new(layout: &'a Layout<'a>, force: &'a Force, reference: bool) -> LineScanner<'a> {
        let keep_column = match &force.keep {
            Some((column, _)) if reference => layout.names.iter().position(|n| n == column),
            _ => None,
        };
        let text_finders = if reference {
            force.texts.iter().map(memmem::Finder::new).collect()
        } else {
            Vec::new()
        };
        LineScanner {
            layout,
            force,
            reference,
            keep_column,
            text_finders,
            last_location: b"\n",
            last_forced: false,
        }
    }

    /// One data line: its Location's hash, whether it hands its group over, and the
    /// `Whole` bits it sets.
    fn line(&mut self, text: &'a [u8]) -> (u64, bool, u32) {
        let layout = self.layout;
        let mut whole = 0u32;
        let mut forced = false;
        let mut location: &[u8] = b"";
        let mut kept_value: &[u8] = b"";
        let mut columns = 0usize;
        let mut start = 0usize;
        let finders = &self.text_finders;
        let mut field = |index: usize, value: &'a [u8], whole: &mut u32| {
            if layout.key.contains(&Some(index)) {
                if memchr2(b'|', b'#', value).is_some() {
                    *whole |= 1 << Whole::KeySeparator as u32;
                }
                // A rewrite of a key column's text would move the row to another key.
                if finders.iter().any(|f| f.find(value).is_some()) {
                    *whole |= 1 << Whole::ListedTextInKey as u32;
                }
                if layout.key[0] == Some(index) {
                    // The comparator labels a record by its Location, a space and its
                    // Uploaded_variation; a space in the Location would let one label
                    // name rows of two groups.
                    if memchr(b' ', value).is_some() {
                        *whole |= 1 << Whole::LocationSpace as u32;
                    }
                    location = value;
                }
            }
            if self.keep_column == Some(index) {
                kept_value = value;
            }
            if layout.extra == Some(index) && value != b"-" && overwrites_a_key(value) {
                *whole |= 1 << Whole::KeyOverwritten as u32;
            }
        };
        for tab in memchr_iter(b'\t', text) {
            field(columns, &text[start..tab], &mut whole);
            columns += 1;
            start = tab + 1;
        }
        field(columns, &text[start..], &mut whole);
        columns += 1;
        if columns != layout.names.len() {
            whole |= 1 << Whole::LineWidth as u32;
        }
        if location == self.last_location {
            forced |= self.last_forced;
        } else {
            self.last_forced = self.force.locations.contains(location);
            self.last_location = location;
            forced |= self.last_forced;
        }
        forced |= past_a_reader_s_limit(text);
        if self.reference {
            if let Some((_, kept)) = &self.force.keep {
                // A side without the column reads every row's value as absent, which the
                // comparator reads as the empty value.
                forced |= !kept.contains(kept_value);
            }
            forced |= self.text_finders.iter().any(|f| f.find(text).is_some());
        }
        (hash(location), forced, whole)
    }
}

/// Whether an `Extra` value holds a pair named as a key column: the comparator reads each
/// pair as a field of its own after the columns, so such a pair replaces the column the
/// entry's key is built from.
fn overwrites_a_key(extra: &[u8]) -> bool {
    extra.split(|&b| b == b';').any(|pair| {
        let name = memchr(b'=', pair).map_or(pair, |i| &pair[..i]);
        KEY_COLUMNS.contains(&name)
    })
}

/// The longest run of digits every build of the comparator's interpreter converts to an
/// integer: a longer one may be refused.
const DIGIT_RUN: usize = 640;
/// The opening brackets under which an object's nesting stays well inside the
/// interpreter's call depth.
const OPENING_BRACKETS: usize = 64;

/// Whether the comparator may fail to read a value of the line: it converts some values
/// (a whole number, an object's text) before it compares them, and its interpreter takes a
/// digit outside ASCII for a digit without converting it, refuses a number past a length,
/// and stops at an object nested past its call depth. A line of such a shape hands its
/// group over, so the comparator meets it on the whole outputs and on the remainder alike.
fn past_a_reader_s_limit(text: &[u8]) -> bool {
    !text.is_ascii()
        || memchr2_iter(b'[', b'{', text)
            .nth(OPENING_BRACKETS)
            .is_some()
        || (text.len() > DIGIT_RUN
            && text
                .split(|b| !b.is_ascii_digit())
                .any(|run| run.len() > DIGIT_RUN))
}

fn scan_lines<'a>(
    side: &'a Side,
    lines: &'a [Line],
    layout: &'a Layout<'a>,
    force: &'a Force,
    reference: bool,
) -> Scan {
    let mut scanner = LineScanner::new(layout, force, reference);
    let mut scan = Scan {
        location: Vec::with_capacity(lines.len()),
        forced: Vec::with_capacity(lines.len()),
        whole: 0,
    };
    for line in lines {
        if line.flags & HEADER != 0 {
            scan.location.push(0);
            scan.forced.push(false);
            continue;
        }
        let (location, forced, whole) = scanner.line(side.text(line));
        scan.location.push(location);
        scan.forced.push(forced);
        scan.whole |= whole;
    }
    scan
}

pub fn scan<'a>(
    side: &'a Side,
    layout: &'a Layout<'a>,
    force: &'a Force,
    reference: bool,
    threads: usize,
) -> Scan {
    let per_worker = side.lines.len().div_ceil(threads.max(1)).max(1);
    let parts: Vec<Scan> = std::thread::scope(|scope| {
        let workers: Vec<_> = side
            .lines
            .chunks(per_worker)
            .map(|lines| scope.spawn(move || scan_lines(side, lines, layout, force, reference)))
            .collect();
        workers
            .into_iter()
            .map(|w| w.join().expect("line scanner thread"))
            .collect()
    });
    let mut all = Scan {
        location: Vec::with_capacity(side.lines.len()),
        forced: Vec::with_capacity(side.lines.len()),
        whole: 0,
    };
    for part in parts {
        all.location.extend(part.location);
        all.forced.extend(part.forced);
        all.whole |= part.whole;
    }
    all
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_extra_pair_named_as_a_key_column_is_seen() {
        assert!(overwrites_a_key(b"IMPACT=LOW;Location=1:2"));
        assert!(overwrites_a_key(b"Feature"));
        assert!(overwrites_a_key(b"Allele=;X=1"));
        assert!(!overwrites_a_key(b"IMPACT=LOW;Locations=1:2;X=Location"));
        assert!(!overwrites_a_key(b""));
    }

    #[test]
    fn a_line_past_a_limit_of_the_comparator_s_interpreter_is_seen() {
        assert!(!past_a_reader_s_limit(
            b"rec1\t21:1001\tIMPACT=LOW;N=12345678901234567890"
        ));
        assert!(!past_a_reader_s_limit(b""));
        // A digit outside ASCII, and any other byte outside it.
        assert!(past_a_reader_s_limit("N=\u{b2}".as_bytes()));
        assert!(past_a_reader_s_limit("name=\u{e9}".as_bytes()));
        // A run of digits at the length every build converts, and one digit past it.
        let digits = "7".repeat(DIGIT_RUN);
        assert!(!past_a_reader_s_limit(
            format!("N={digits};M={digits}").as_bytes()
        ));
        assert!(past_a_reader_s_limit(format!("N={digits}7").as_bytes()));
        // Opening brackets at the bound, whatever closes them, and one past it.
        let open = "[{".repeat(OPENING_BRACKETS / 2);
        assert!(!past_a_reader_s_limit(format!("V={open}]]").as_bytes()));
        assert!(past_a_reader_s_limit(format!("V={open};W={{}}").as_bytes()));
    }

    #[test]
    fn a_list_file_keeps_an_empty_line_as_the_empty_value() {
        let path =
            std::env::temp_dir().join(format!("vep-formats-remainder-list-{}", std::process::id()));
        fs::write(&path, "Transcript\n\n-\n").unwrap();
        let lines = list_lines(&path).unwrap();
        fs::remove_file(&path).unwrap();
        assert_eq!(lines, [b"Transcript".to_vec(), b"".to_vec(), b"-".to_vec()]);
    }
}
