// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The column line of a default- or tab-format output and the places of the columns the
//! comparator builds an entry's key from, as `parse_default` reads them.

use crate::input::{Side, COLUMNS, HEADER};

/// The columns the comparator joins into an entry's key.
pub const KEY_COLUMNS: [&[u8]; 4] = [b"Location", b"Allele", b"Feature", b"Feature_type"];
/// The column whose `key=value` pairs the comparator reads as fields of their own.
pub const EXTRA: &[u8] = b"Extra";

/// Why every line is handed over: a shape of the input under which a line's text does not
/// settle its entry, listed in the order the reason is reported when several hold.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Whole {
    HeaderAfterData = 0,
    NoColumnLine = 1,
    RepeatedColumn = 2,
    ColumnLinesDiffer = 3,
    LineWidth = 4,
    KeySeparator = 5,
    LocationSpace = 6,
    KeyOverwritten = 7,
    ListedTextInKey = 8,
    OrderUnsettled = 9,
}

impl Whole {
    pub const ALL: [Whole; 10] = [
        Whole::HeaderAfterData,
        Whole::NoColumnLine,
        Whole::RepeatedColumn,
        Whole::ColumnLinesDiffer,
        Whole::LineWidth,
        Whole::KeySeparator,
        Whole::LocationSpace,
        Whole::KeyOverwritten,
        Whole::ListedTextInKey,
        Whole::OrderUnsettled,
    ];

    pub fn text(self) -> &'static str {
        match self {
            Whole::HeaderAfterData => "a header line follows a data line",
            Whole::NoColumnLine => "a data line precedes every column line",
            Whole::RepeatedColumn => "the column line names a column twice",
            Whole::ColumnLinesDiffer => "the two outputs' column lines differ",
            Whole::LineWidth => "a data line's width is not the column line's",
            Whole::KeySeparator => "a key column holds | or #",
            Whole::LocationSpace => "a Location holds a space",
            Whole::KeyOverwritten => "an Extra pair is named as a key column",
            Whole::ListedTextInKey => {
                "a listed text stands in the reference's column line or a key column"
            }
            Whole::OrderUnsettled => "the kept-back lines found no common order",
        }
    }
}

pub struct Layout<'a> {
    /// The column line without its `#`; empty for an output without one.
    pub line: &'a [u8],
    pub names: Vec<&'a [u8]>,
    /// The places of `KEY_COLUMNS`, in that order; `None` for a column the line lacks.
    pub key: [Option<usize>; 4],
    pub extra: Option<usize>,
}

impl<'a> Layout<'a> {
    /// The layout in force for every data line of `side`, or the reason there is none: the
    /// comparator reads a column line wherever it stands, so only an output whose header
    /// lines all precede its data lines has one layout.
    pub fn of(side: &'a Side) -> Result<Layout<'a>, Whole> {
        let first_data = side
            .lines
            .iter()
            .position(|l| l.flags & HEADER == 0)
            .unwrap_or(side.lines.len());
        if side.lines[first_data..]
            .iter()
            .any(|l| l.flags & HEADER != 0)
        {
            return Err(Whole::HeaderAfterData);
        }
        let column_line = side.lines[..first_data]
            .iter()
            .rev()
            .find(|l| l.flags & COLUMNS != 0);
        let line: &[u8] = match column_line {
            Some(line) => &side.text(line)[1..],
            None if first_data == side.lines.len() => b"",
            None => return Err(Whole::NoColumnLine),
        };
        let names: Vec<&[u8]> = match column_line {
            Some(_) => line.split(|&b| b == b'\t').collect(),
            None => Vec::new(),
        };
        let place = |name: &[u8]| names.iter().position(|n| *n == name);
        for (i, name) in names.iter().enumerate() {
            if names[..i].contains(name) {
                return Err(Whole::RepeatedColumn);
            }
        }
        Ok(Layout {
            line,
            key: KEY_COLUMNS.map(place),
            extra: place(EXTRA),
            names,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn side(text: &str) -> Side {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-remainder-layout-{}-{:x}",
            std::process::id(),
            text.len() * 31 + text.bytes().map(usize::from).sum::<usize>()
        ));
        std::fs::write(&path, text).unwrap();
        let side = Side::read(&path, 1).unwrap();
        std::fs::remove_file(&path).unwrap();
        side
    }

    #[test]
    fn the_last_column_line_before_the_data_is_the_layout() {
        let s = side(
            "## note\n#old\n#Uploaded_variation\tLocation\tAllele\tFeature\tExtra\nv\tl\ta\tf\t-\n",
        );
        let layout = Layout::of(&s).unwrap();
        assert_eq!(layout.names.len(), 5);
        assert_eq!(layout.key, [Some(1), Some(2), Some(3), None]);
        assert_eq!(layout.extra, Some(4));
    }

    #[test]
    fn shapes_without_one_layout_are_named() {
        assert_eq!(
            Layout::of(&side("#a\tb\nx\ty\n#a\tb\n")).err(),
            Some(Whole::HeaderAfterData)
        );
        assert_eq!(
            Layout::of(&side("## note\nx\ty\n")).err(),
            Some(Whole::NoColumnLine)
        );
        assert_eq!(
            Layout::of(&side("#a\ta\nx\ty\n")).err(),
            Some(Whole::RepeatedColumn)
        );
        assert!(Layout::of(&side("## only a note\n"))
            .unwrap()
            .names
            .is_empty());
    }
}
