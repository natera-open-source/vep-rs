// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The default and tab formats: the column line, the columns `parse_default`
//! (`crate::parse::text`) builds an entry's key from, what hands a line's group over,
//! the shapes under which a line's text does not settle its entry, and the populated
//! fields of a line as `populated_counts` counts them.
//!
//! A record group is every data line of a side that shares a Location.

use super::digest;
use super::grouping::{Grouping, Scanned, Side, Tally};
use super::Plan;
use crate::error::{Error, Result};
use memchr::{memchr, memchr2, memchr2_iter, memchr_iter, memmem};
use std::collections::{HashMap, HashSet};
use std::ops::Range;
use std::sync::Arc;

/// The columns the comparison joins into an entry's key.
const KEY_COLUMNS: [&[u8]; 4] = [b"Location", b"Allele", b"Feature", b"Feature_type"];
/// The column whose `key=value` pairs the comparison reads as fields of their own.
const EXTRA: &[u8] = b"Extra";

const HEADER_AFTER_DATA: u32 = 0;
const NO_COLUMN_LINE: u32 = 1;
const REPEATED_COLUMN: u32 = 2;
const COLUMN_LINES_DIFFER: u32 = 3;
const LINE_WIDTH: u32 = 4;
const KEY_SEPARATOR: u32 = 5;
const LOCATION_SPACE: u32 = 6;
const KEY_OVERWRITTEN: u32 = 7;
const LISTED_TEXT_IN_KEY: u32 = 8;

/// What a column is read for.
const KEY: u8 = 1;
const LOCATION: u8 = 2;
const KEPT: u8 = 4;
const PAIRS: u8 = 8;

/// The Locations, texts and column values under which a group is handed over whole.
struct Force {
    locations: HashSet<Vec<u8>>,
    texts: Vec<Vec<u8>>,
    keep: Option<(Vec<u8>, HashSet<Vec<u8>>)>,
}

pub(crate) struct Columns {
    force: Arc<Force>,
}

fn bytes(values: &[String]) -> impl Iterator<Item = Vec<u8>> + '_ {
    values.iter().map(|value| value.as_bytes().to_vec())
}

impl Columns {
    pub(crate) fn new(plan: &Plan) -> Result<Columns> {
        if plan.force_text.iter().any(String::is_empty) {
            return Err(Error::comparator(
                "the texts that hand a group over: an empty text matches every line",
            ));
        }
        Ok(Columns {
            force: Arc::new(Force {
                locations: bytes(&plan.force_locations).collect(),
                texts: bytes(&plan.force_text).collect(),
                keep: plan
                    .keep_values
                    .as_ref()
                    .map(|(column, values)| (column.as_bytes().to_vec(), bytes(values).collect())),
            }),
        })
    }
}

/// The last column line before a side's first data line, without its `#`.
#[derive(Default)]
pub(crate) struct ColumnLine(Option<Vec<u8>>);

impl Grouping for Columns {
    type Header = ColumnLine;
    type Side = ColumnSide;

    const SHAPES: &'static [&'static str] = &[
        "a header line follows a data line",
        "a data line precedes every column line",
        "the column line names a column twice",
        "the two outputs' column lines differ",
        "a data line's width is not the column line's",
        "a key column holds | or #",
        "a Location holds a space",
        "an Extra pair is named as a key column",
        "a listed text stands in the reference's column line or a key column",
    ];

    fn is_header(line: &[u8]) -> bool {
        line.first() == Some(&b'#')
    }

    fn header_line(header: &mut ColumnLine, line: &[u8]) {
        // A header line that does not start with `##` names the columns.
        if !line.starts_with(b"##") {
            header.0 = Some(line[1..].to_vec());
        }
    }

    fn side(&self, header: ColumnLine, data: bool, reference: bool) -> ColumnSide {
        let names: Vec<Vec<u8>> = match &header.0 {
            Some(line) => line.split(|&b| b == b'\t').map(<[u8]>::to_vec).collect(),
            None => Vec::new(),
        };
        let place = |name: &[u8]| names.iter().position(|n| n == name);
        let mut roles = vec![0u8; names.len()];
        let key = KEY_COLUMNS.map(place);
        for at in key.into_iter().flatten() {
            roles[at] |= KEY;
        }
        if let Some(at) = key[0] {
            roles[at] |= LOCATION;
        }
        if let Some(at) = place(EXTRA) {
            roles[at] |= PAIRS;
        }
        // The lists name reference lines alone.
        let texts: &[Vec<u8>] = if reference { &self.force.texts } else { &[] };
        if let Some((column, _)) = self.force.keep.as_ref().filter(|_| reference) {
            if let Some(at) = place(column) {
                roles[at] |= KEPT;
            }
        }
        ColumnSide {
            force: Arc::clone(&self.force),
            reference,
            finders: texts
                .iter()
                .map(|text| memmem::Finder::new(text).into_owned())
                .collect(),
            repeated: (1..names.len()).any(|i| names[..i].contains(&names[i])),
            columned: header.0.is_some(),
            line: header.0.unwrap_or_default(),
            data,
            names,
            roles,
        }
    }

    fn between(&self, expected: &ColumnSide, actual: &ColumnSide) -> Option<u32> {
        if expected.names != actual.names {
            return Some(COLUMN_LINES_DIFFER);
        }
        // A rewrite of the reference's column line would rename a column on one side alone.
        self.force
            .texts
            .iter()
            .any(|text| memmem::find(&expected.line, text).is_some())
            .then_some(LISTED_TEXT_IN_KEY)
    }
}

pub(crate) struct ColumnSide {
    force: Arc<Force>,
    reference: bool,
    finders: Vec<memmem::Finder<'static>>,
    /// The column line without its `#`; empty for a side without one.
    line: Vec<u8>,
    columned: bool,
    data: bool,
    names: Vec<Vec<u8>>,
    repeated: bool,
    /// Per column, what it is read for.
    roles: Vec<u8>,
}

/// The Location of the line before and whether the lists name it.
#[derive(Default)]
pub(crate) struct LastLocation {
    seen: bool,
    location: Vec<u8>,
    listed: bool,
}

/// FNV-1a over the bytes of a Location or of a record's key.
pub(super) fn hash(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h = (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// Whether an `Extra` value holds a pair named as a key column: the comparison reads each
/// pair as a field of its own after the columns, so such a pair replaces the column the
/// entry's key is built from.
fn overwrites_a_key(extra: &[u8]) -> bool {
    let mut start = 0;
    for end in memchr_iter(b';', extra).chain(std::iter::once(extra.len())) {
        let pair = &extra[start..end];
        // The first letters of the key columns' names.
        if matches!(pair.first(), Some(b'L' | b'A' | b'F')) {
            let name = memchr(b'=', pair).map_or(pair, |i| &pair[..i]);
            if KEY_COLUMNS.contains(&name) {
                return true;
            }
        }
        start = end + 1;
    }
    false
}

/// The longest run of digits every build of the comparison's interpreter converts to an
/// integer: a longer one may be refused.
const DIGIT_RUN: usize = 640;
/// The opening brackets under which an object's nesting stays well inside the
/// interpreter's call depth.
const OPENING_BRACKETS: usize = 64;

/// Whether the comparison may fail to read a value of the line: it converts some values
/// (a whole number, an object's text) before it compares them, takes a digit outside
/// ASCII for a digit without converting it, refuses a number past a length, and stops at
/// an object nested past its call depth. A line of
/// such a shape hands its group over, so the comparison meets it whatever is kept back.
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

impl Side for ColumnSide {
    type Memo = LastLocation;

    fn shape(&self, header_after_data: bool) -> Option<u32> {
        // The comparison reads a column line wherever it stands, so only a side whose
        // header lines all precede its data lines has one layout.
        if header_after_data {
            Some(HEADER_AFTER_DATA)
        } else if self.data && !self.columned {
            Some(NO_COLUMN_LINE)
        } else if self.repeated {
            Some(REPEATED_COLUMN)
        } else {
            None
        }
    }

    fn scan(&self, memo: &mut LastLocation, text: &[u8]) -> Scanned {
        let mut shapes = 0u32;
        let mut location: &[u8] = b"";
        let mut kept_value: &[u8] = b"";
        let mut columns = 0usize;
        let mut start = 0usize;
        let ends = memchr_iter(b'\t', text).chain(std::iter::once(text.len()));
        for end in ends {
            let value = &text[start..end];
            let role = self.roles.get(columns).copied().unwrap_or(0);
            if role & KEY != 0 {
                if memchr2(b'|', b'#', value).is_some() {
                    shapes |= 1 << KEY_SEPARATOR;
                }
                // A rewrite of a key column's text would move the row to another key.
                if self.finders.iter().any(|f| f.find(value).is_some()) {
                    shapes |= 1 << LISTED_TEXT_IN_KEY;
                }
                if role & LOCATION != 0 {
                    // The comparison labels a record by its Location, a space and its
                    // Uploaded_variation; a space in the Location would let one label
                    // name rows of two groups.
                    if memchr(b' ', value).is_some() {
                        shapes |= 1 << LOCATION_SPACE;
                    }
                    location = value;
                }
            }
            if role & KEPT != 0 {
                kept_value = value;
            }
            if role & PAIRS != 0 && value != b"-" && overwrites_a_key(value) {
                shapes |= 1 << KEY_OVERWRITTEN;
            }
            columns += 1;
            start = end + 1;
        }
        if columns != self.names.len() {
            shapes |= 1 << LINE_WIDTH;
        }
        if !memo.seen || memo.location != location {
            memo.seen = true;
            memo.location.clear();
            memo.location.extend_from_slice(location);
            memo.listed = self.force.locations.contains(location);
        }
        let fragile = past_a_reader_s_limit(text);
        let mut forced = memo.listed || fragile;
        if self.reference {
            if let Some((_, kept)) = &self.force.keep {
                // A side without the column reads every row's value as absent, which the
                // comparison reads as the empty value.
                forced |= !kept.contains(kept_value);
            }
            forced |= self.finders.iter().any(|f| f.find(text).is_some());
        }
        Scanned {
            key: hash(location),
            hash: digest::of_line(text),
            forced,
            fragile,
            shapes,
            strict: false,
            carried: 0,
        }
    }

    fn populated(
        &self,
        block: &[u8],
        kept: &[Range<usize>],
        tally: &mut Tally,
        _entries: &mut Vec<u32>,
    ) {
        let mut counter = Counter::new(&self.names);
        let pairs = self.roles.iter().position(|role| role & PAIRS != 0);
        for line in kept {
            counter.line(&block[line.clone()], pairs);
        }
        for (name, &count) in counter.names.iter().zip(&counter.counts) {
            if count > 0 {
                *tally.populated.entry(name.to_vec()).or_insert(0) += count;
            }
        }
    }
}

fn valued(value: &[u8]) -> bool {
    !value.is_empty() && value != b"-"
}

/// Per field name, the lines that carry a value, as `populated_counts` counts entries:
/// every column, and every `key=value` pair of the `Extra` column as a field of its own,
/// the last pair of a name standing for every pair of that name before it and for the
/// column of that name; `-` and the empty text are no value.
struct Counter<'a> {
    ids: HashMap<&'a [u8], u32>,
    names: Vec<&'a [u8]>,
    counts: Vec<u64>,
    /// Per name, the line that last carried it and its place in that line's `fields`.
    seen: Vec<(u64, u32)>,
    fields: Vec<(u32, bool)>,
    column_ids: Vec<u32>,
    /// Per place among a line's pairs, the name of the last pair there.
    pair_ids: Vec<(&'a [u8], u32)>,
    serial: u64,
}

impl<'a> Counter<'a> {
    fn new(columns: &'a [Vec<u8>]) -> Counter<'a> {
        let mut counter = Counter {
            ids: HashMap::new(),
            names: Vec::new(),
            counts: Vec::new(),
            seen: Vec::new(),
            fields: Vec::new(),
            column_ids: Vec::new(),
            pair_ids: Vec::new(),
            serial: 0,
        };
        counter.column_ids = columns.iter().map(|name| counter.id(name)).collect();
        counter
    }

    fn id(&mut self, name: &'a [u8]) -> u32 {
        if let Some(&id) = self.ids.get(name) {
            return id;
        }
        let id = self.names.len() as u32;
        self.ids.insert(name, id);
        self.names.push(name);
        self.counts.push(0);
        self.seen.push((0, 0));
        id
    }

    fn field(&mut self, id: u32, valued: bool) {
        let seen = &mut self.seen[id as usize];
        if seen.0 == self.serial {
            self.fields[seen.1 as usize].1 = valued;
        } else {
            *seen = (self.serial, self.fields.len() as u32);
            self.fields.push((id, valued));
        }
    }

    /// Counts one line, `pairs` being the place of its `Extra` column.
    fn line(&mut self, text: &'a [u8], pairs: Option<usize>) {
        self.serial += 1;
        self.fields.clear();
        let mut extra_value: Option<&'a [u8]> = None;
        let mut start = 0;
        let ends = memchr_iter(b'\t', text).chain(std::iter::once(text.len()));
        for (column, end) in ends.enumerate() {
            let value = &text[start..end];
            if let Some(&id) = self.column_ids.get(column) {
                self.field(id, valued(value));
            }
            if pairs == Some(column) {
                extra_value = Some(value);
            }
            start = end + 1;
        }
        // The comparison leaves an absent `Extra` (`-`) unread and, after the pairs, sets
        // the field `Extra` to the column's own text again.
        if let Some(extra_value) = extra_value.filter(|v| *v != b"-") {
            let mut start = 0;
            let ends = memchr_iter(b';', extra_value).chain(std::iter::once(extra_value.len()));
            for (place, end) in ends.enumerate() {
                let pair = &extra_value[start..end];
                start = end + 1;
                let (name, value) = match memchr(b'=', pair) {
                    Some(i) => (&pair[..i], &pair[i + 1..]),
                    None => (pair, &pair[pair.len()..]),
                };
                if name == EXTRA {
                    continue;
                }
                // A line's pairs mostly carry the names of the line before, in its order.
                let id = match self.pair_ids.get(place) {
                    Some(&(last, id)) if last == name => id,
                    _ => {
                        let id = self.id(name);
                        if place >= self.pair_ids.len() {
                            // A place no pair has named holds `;`, which no pair's name
                            // holds, the empty name being a pair's.
                            self.pair_ids.resize(place + 1, (b";", u32::MAX));
                        }
                        self.pair_ids[place] = (name, id);
                        id
                    }
                };
                self.field(id, valued(value));
            }
        }
        for &(id, valued) in &self.fields {
            if valued {
                self.counts[id as usize] += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stream::lines::lines;
    use crate::stream::Format;

    fn plan() -> Plan {
        Plan {
            format: Format::Default,
            force_locations: Vec::new(),
            force_text: Vec::new(),
            keep_values: None,
            unordered: false,
            threads: 1,
            input: None,
            force_records: Vec::new(),
            floor_fields: Vec::new(),
        }
    }

    /// The side of `text`'s header lines and the text's data lines.
    fn side(columns: &Columns, text: &str, reference: bool) -> (ColumnSide, Vec<String>) {
        let mut header = ColumnLine::default();
        let mut data = Vec::new();
        for line in text.lines() {
            if data.is_empty() && Columns::is_header(line.as_bytes()) {
                Columns::header_line(&mut header, line.as_bytes());
            } else {
                data.push(line.to_string());
            }
        }
        (columns.side(header, !data.is_empty(), reference), data)
    }

    #[test]
    fn the_last_column_line_before_the_data_is_the_layout() {
        let columns = Columns::new(&plan()).unwrap();
        let (side, _) = side(
            &columns,
            "## note\n#old\n#Uploaded_variation\tLocation\tAllele\tFeature\tExtra\nv\tl\ta\tf\t-\n",
            true,
        );
        assert_eq!(side.names.len(), 5);
        assert_eq!(side.roles, [0, KEY | LOCATION, KEY, KEY, PAIRS]);
        assert_eq!(side.shape(false), None);
    }

    #[test]
    fn shapes_without_one_layout_are_named() {
        let columns = Columns::new(&plan()).unwrap();
        let shape = |text: &str, header_after_data: bool| {
            side(&columns, text, true)
                .0
                .shape(header_after_data)
                .map(|shape| Columns::SHAPES[shape as usize])
        };
        assert_eq!(
            shape("#a\tb\nx\ty\n", true),
            Some("a header line follows a data line")
        );
        assert_eq!(
            shape("## note\nx\ty\n", false),
            Some("a data line precedes every column line")
        );
        assert_eq!(
            shape("#a\ta\nx\ty\n", false),
            Some("the column line names a column twice")
        );
        let (only_a_note, _) = side(&columns, "## only a note\n", true);
        assert!(only_a_note.names.is_empty());
        assert_eq!(only_a_note.shape(false), None);
        // A column line of a `#` alone names one column, the empty name.
        assert_eq!(side(&columns, "#\n", true).0.names, [b""]);
    }

    #[test]
    fn an_extra_pair_named_as_a_key_column_is_seen() {
        assert!(overwrites_a_key(b"IMPACT=LOW;Location=1:2"));
        assert!(overwrites_a_key(b"Feature"));
        assert!(overwrites_a_key(b"Allele=;X=1"));
        assert!(!overwrites_a_key(b"IMPACT=LOW;Locations=1:2;X=Location"));
        assert!(!overwrites_a_key(b""));
    }

    #[test]
    fn a_line_past_a_limit_of_the_comparison_s_interpreter_is_seen() {
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
    fn a_line_names_its_location_and_what_hands_its_group_over() {
        let mut plan = plan();
        plan.force_locations = vec!["1:7".to_string(), String::new()];
        plan.force_text = vec!["misspelt".to_string()];
        plan.keep_values = Some((
            "Feature_type".to_string(),
            vec!["Transcript".to_string(), "-".to_string()],
        ));
        let columns = Columns::new(&plan).unwrap();
        let header = "#Location\tAllele\tFeature\tFeature_type\tExtra\n";
        let scan = |reference: bool, row: &str| {
            let (side, _) = side(&columns, header, reference);
            side.scan(&mut LastLocation::default(), row.as_bytes())
        };
        let plain = scan(true, "1:5\tA\ttx\tTranscript\tIMPACT=LOW");
        assert_eq!(plain.key, hash(b"1:5"));
        assert!(!plain.forced && !plain.fragile && plain.shapes == 0);
        assert!(scan(true, "1:7\tA\ttx\tTranscript\t-").forced);
        assert!(scan(false, "1:7\tA\ttx\tTranscript\t-").forced);
        // The empty text among the listed Locations is the empty Location.
        assert!(scan(true, "\tA\ttx\tTranscript\t-").forced);
        assert!(!scan(true, "1:70\tA\ttx\tTranscript\t-").forced);
        // The kept values and the listed texts name reference lines alone.
        assert!(scan(true, "1:5\tA\treg\tRegulatoryFeature\t-").forced);
        assert!(!scan(false, "1:5\tA\treg\tRegulatoryFeature\t-").forced);
        assert!(scan(true, "1:5\tA\ttx\tTranscript\tmisspeltKey=1").forced);
        assert!(!scan(false, "1:5\tA\ttx\tTranscript\tmisspeltKey=1").forced);
        let shapes = |row: &str| scan(true, row).shapes;
        assert_eq!(shapes("1:5\tA\ttx"), 1 << LINE_WIDTH);
        assert_eq!(shapes("1:5\tA|C\ttx\tTranscript\t-"), 1 << KEY_SEPARATOR);
        assert_eq!(shapes("1:5 6\tA\ttx\tTranscript\t-"), 1 << LOCATION_SPACE);
        assert_eq!(
            shapes("1:5\tA\ttx\tTranscript\tAllele=C"),
            1 << KEY_OVERWRITTEN
        );
        assert_eq!(
            shapes("1:5\tA\tmisspelt7\tTranscript\t-"),
            1 << LISTED_TEXT_IN_KEY
        );
        let fragile = scan(false, "1:5\tA\ttx\tTranscript\tN=\u{b2}");
        assert!(fragile.forced && fragile.fragile);
        assert!(Columns::new(&Plan {
            force_text: vec![String::new()],
            ..plan
        })
        .is_err());
    }

    #[test]
    fn a_side_without_the_kept_column_reads_the_empty_value() {
        let mut plan = plan();
        plan.keep_values = Some(("Feature_type".to_string(), vec!["Transcript".to_string()]));
        let columns = Columns::new(&plan).unwrap();
        let (without, _) = side(&columns, "#Location\tAllele\n", true);
        assert!(without.scan(&mut LastLocation::default(), b"1:5\tA").forced);
        plan.keep_values = Some(("Feature_type".to_string(), vec![String::new()]));
        let columns = Columns::new(&plan).unwrap();
        let (without, _) = side(&columns, "#Location\tAllele\n", true);
        assert!(!without.scan(&mut LastLocation::default(), b"1:5\tA").forced);
    }

    #[test]
    fn the_two_column_lines_are_held_against_each_other_and_the_listed_texts() {
        let mut plan = plan();
        plan.force_text = vec!["misspelt".to_string()];
        let columns = Columns::new(&plan).unwrap();
        let of = |text: &str, reference: bool| side(&columns, text, reference).0;
        let shape = |e: &str, a: &str| {
            columns
                .between(&of(e, true), &of(a, false))
                .map(|shape| Columns::SHAPES[shape as usize])
        };
        assert_eq!(shape("#Location\tGene\n", "#Location\tGene\n"), None);
        assert_eq!(
            shape("#Location\tGene\n", "#Location\tGENE\n"),
            Some("the two outputs' column lines differ")
        );
        assert_eq!(
            shape("#\n", "## a note\n"),
            Some("the two outputs' column lines differ")
        );
        assert_eq!(
            shape("#Location\tmisspelt100way\n", "#Location\tmisspelt100way\n"),
            Some("a listed text stands in the reference's column line or a key column")
        );
    }

    fn counts(text: &str) -> Vec<(String, u64)> {
        let columns = Columns::new(&plan()).unwrap();
        let (side, _) = side(&columns, text, true);
        let kept: Vec<Range<usize>> = lines(text.as_bytes())
            .filter(|line| !Columns::is_header(&text.as_bytes()[line.clone()]))
            .collect();
        let mut tally = Tally::default();
        let mut entries = Vec::new();
        let (first, second) = kept.split_at(kept.len() / 2);
        side.populated(text.as_bytes(), first, &mut tally, &mut entries);
        side.populated(text.as_bytes(), second, &mut tally, &mut entries);
        assert!(entries.is_empty());
        tally
            .populated
            .into_iter()
            .map(|(name, n)| (String::from_utf8(name).unwrap(), n))
            .collect()
    }

    fn pairs(expected: &[(&str, u64)]) -> Vec<(String, u64)> {
        expected
            .iter()
            .map(|(n, c)| ((*n).to_string(), *c))
            .collect()
    }

    #[test]
    fn columns_and_extra_pairs_count_once_per_line() {
        let text =
            "#Location\tGene\tExtra\nl1\tg\tIMPACT=LOW;STRAND=1\nl2\t-\tIMPACT=HIGH\nl3\t\t-\n";
        assert_eq!(
            counts(text),
            pairs(&[
                ("Extra", 2),
                ("Gene", 1),
                ("IMPACT", 2),
                ("Location", 3),
                ("STRAND", 1)
            ])
        );
    }

    #[test]
    fn a_dash_an_empty_value_and_a_bare_name_are_no_value() {
        let text = "#Location\tExtra\nl1\tA=-;B=;C;D=x\n";
        assert_eq!(
            counts(text),
            pairs(&[("D", 1), ("Extra", 1), ("Location", 1)])
        );
    }

    #[test]
    fn the_last_pair_of_a_name_stands_for_the_pairs_before_it_and_the_column() {
        let text = "#Location\tGene\tExtra\nl1\tg\tA=x;A=-;Gene=-;B=-;B=y\n";
        assert_eq!(
            counts(text),
            pairs(&[("B", 1), ("Extra", 1), ("Location", 1)])
        );
    }

    #[test]
    fn an_empty_extra_and_a_pair_named_extra_leave_the_column_s_own_text() {
        let text = "#Location\tExtra\nl1\t\nl2\tExtra=-;A=1\n";
        assert_eq!(
            counts(text),
            pairs(&[("A", 1), ("Extra", 1), ("Location", 2)])
        );
    }

    #[test]
    fn a_pair_without_a_name_is_a_field_of_its_own_at_any_place() {
        // The first line's pair named `Extra` leaves its place among the pairs unnamed.
        let text = "#Location\tExtra\nl1\tExtra=-;A=1\nl2\t;A=1\nl3\t=x;A=-\nl4\t-\nl5\t-\nl6\t-\n";
        assert_eq!(
            counts(text),
            pairs(&[("", 1), ("A", 2), ("Extra", 3), ("Location", 6)])
        );
    }
}
