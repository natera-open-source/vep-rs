// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The header of a discordant table: which of the two schemas it is, where each named
//! column is read, and the start of a Location.
//!
//! A row is the mapping `csv.DictReader` makes of it: a name that the header holds more
//! than once reads its last column, and a column behind the row's last field holds no
//! value.

use crate::csv::Record;
use memchr::memchr;

/// The name of the column the classes are written to.
pub const CLASS_COLUMN: &[u8] = b"mismatch_type";

/// The name under which a row keeps the start of its Location.
const START_COLUMN: &[u8] = b"_start";

/// Where the value that a column of the written table shows is read.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Shown {
    /// The field of the row at this index, empty where the row ends before it.
    Field(usize),
    /// The start of the row's Location.
    Start,
    /// The row's class.
    Class,
}

#[derive(Debug)]
pub struct Header {
    /// The column names in file order.
    pub names: Vec<Vec<u8>>,
    /// What each column shows in the written table.
    pub shown: Vec<Shown>,
    pub location: usize,
    pub allele: usize,
    pub feature: usize,
    pub consequence: usize,
    /// The column that says which engine wrote the row, and its name.
    pub source: usize,
    pub source_name: &'static str,
    /// The column that names the row's file, and its name.
    pub file: usize,
    pub file_name: &'static str,
    /// The last column named `mismatch_type`, where the table read holds one.
    pub class: Option<usize>,
}

impl Header {
    /// Reads the header record. A table of the structural-variant comparison names its
    /// columns `direction`, `source_file` and `consequence`; any other is read as the tuple
    /// comparator's, with `source`, `file_name` and `consequence_set`.
    pub fn read(record: &Record<'_>) -> Result<Header, String> {
        let names: Vec<Vec<u8>> = record.fields().map(<[u8]>::to_vec).collect();
        let last = |name: &str| names.iter().rposition(|held| held == name.as_bytes());
        let either = |first: &'static str, second: &'static str| match last(first) {
            Some(_) => first,
            None => second,
        };
        let source_name = either("direction", "source");
        let file_name = either("source_file", "file_name");
        let consequence_name = either("consequence", "consequence_set");
        let required = [
            "location",
            "allele",
            "feature",
            source_name,
            file_name,
            consequence_name,
        ];
        let mut missing: Vec<&str> = required
            .into_iter()
            .filter(|name| last(name).is_none())
            .collect();
        missing.sort_unstable();
        let [location, allele, feature, source, file, consequence] = required.map(last);
        let (
            Some(location),
            Some(allele),
            Some(feature),
            Some(source),
            Some(file),
            Some(consequence),
        ) = (location, allele, feature, source, file, consequence)
        else {
            return Err(format!(
                "Unsupported discordant TSV schema; missing required columns: {}",
                missing.join(", ")
            ));
        };
        let shown = names
            .iter()
            .map(|name| match name.as_slice() {
                START_COLUMN => Shown::Start,
                CLASS_COLUMN => Shown::Class,
                _ => Shown::Field(
                    names
                        .iter()
                        .rposition(|held| held == name)
                        .unwrap_or_default(),
                ),
            })
            .collect();
        let class = names.iter().rposition(|held| held == CLASS_COLUMN);
        Ok(Header {
            names,
            shown,
            location,
            allele,
            feature,
            consequence,
            source,
            source_name,
            file,
            file_name,
            class,
        })
    }
}

/// The start of a Location: the text between its first colon and the first of the next
/// colon and the next `-`; the whole Location where it holds no colon.
pub fn start_of(location: &[u8]) -> &[u8] {
    let Some(colon) = memchr(b':', location) else {
        return location;
    };
    let coordinates = &location[colon + 1..];
    let end = coordinates
        .iter()
        .position(|byte| matches!(byte, b':' | b'-'))
        .unwrap_or(coordinates.len());
    &coordinates[..end]
}

/// Whether the value of the source column marks a row that the reference alone wrote.
pub fn is_reference_only(source: Option<&[u8]>) -> bool {
    matches!(source, Some(b"only_perl" | b"missing_in_rust"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::csv::Reader;

    fn header(line: &str) -> Result<Header, String> {
        let mut read = None;
        let mut reader = Reader::new();
        let mut sink = |record: Record<'_>| -> Result<(), String> {
            read = Some(Header::read(&record));
            Ok(())
        };
        reader.feed(line.as_bytes(), &mut sink)?;
        reader.finish(&mut sink)?;
        read.expect("one record")
    }

    #[test]
    fn the_tuple_comparator_s_columns_are_found_by_name() {
        let header =
            header("file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n")
                .unwrap();
        assert_eq!(
            (header.file, header.source, header.location, header.allele),
            (0, 1, 2, 3)
        );
        assert_eq!((header.feature, header.consequence), (4, 6));
        assert_eq!(
            (header.source_name, header.file_name),
            ("source", "file_name")
        );
        assert_eq!(header.class, None);
        assert_eq!(header.shown, (0..7).map(Shown::Field).collect::<Vec<_>>());
    }

    #[test]
    fn the_structural_variant_columns_are_taken_where_both_schemas_are_named() {
        let header = header(
            "source\tfile_name\tconsequence_set\tsource_file\tvariant_type\tdirection\tlocation\tallele\tfeature\tconsequence\n",
        )
        .unwrap();
        assert_eq!((header.file, header.source, header.consequence), (3, 5, 9));
        assert_eq!(
            (header.source_name, header.file_name),
            ("direction", "source_file")
        );
    }

    #[test]
    fn each_of_the_three_columns_falls_to_the_tuple_comparator_s_name_on_its_own() {
        let header =
            header("location\tallele\tfeature\tdirection\tfile_name\tconsequence\tsource\n")
                .unwrap();
        assert_eq!((header.source, header.file, header.consequence), (3, 4, 5));
        assert_eq!(
            (header.source_name, header.file_name),
            ("direction", "file_name")
        );
    }

    #[test]
    fn a_name_held_twice_reads_its_last_column() {
        let header = header(
            "location\tallele\tfeature\tsource\tfile_name\tconsequence_set\tallele\t_start\tmismatch_type\tlocation\tmismatch_type\n",
        )
        .unwrap();
        assert_eq!(
            (header.location, header.allele, header.class),
            (9, 6, Some(10))
        );
        assert_eq!(
            header.shown,
            [
                Shown::Field(9),
                Shown::Field(6),
                Shown::Field(2),
                Shown::Field(3),
                Shown::Field(4),
                Shown::Field(5),
                Shown::Field(6),
                Shown::Start,
                Shown::Class,
                Shown::Field(9),
                Shown::Class,
            ]
        );
    }

    #[test]
    fn the_missing_columns_are_listed_in_code_point_order_under_the_schema_they_fall_to() {
        let refused = |line: &str| header(line).unwrap_err();
        assert_eq!(
            refused("a\tb\n"),
            "Unsupported discordant TSV schema; missing required columns: \
             allele, consequence_set, feature, file_name, location, source"
        );
        assert_eq!(
            refused("direction\tsource_file\tconsequence\tsource\tfile_name\n"),
            "Unsupported discordant TSV schema; missing required columns: allele, feature, location"
        );
        assert_eq!(
            refused("location\tallele\tfeature\tdirection\tconsequence\tsource\n"),
            "Unsupported discordant TSV schema; missing required columns: file_name"
        );
        assert_eq!(
            refused("\u{feff}file_name\tsource\tlocation\tallele\tfeature\tconsequence_set\n"),
            "Unsupported discordant TSV schema; missing required columns: file_name"
        );
        assert_eq!(
            refused("\n"),
            "Unsupported discordant TSV schema; missing required columns: \
             allele, consequence_set, feature, file_name, location, source"
        );
    }

    #[test]
    fn the_start_of_a_location_lies_between_its_first_colon_and_the_next_colon_or_dash() {
        let cases = [
            ("21:25954688-25954689", "25954688"),
            ("chr1:5", "5"),
            ("1:100", "100"),
            ("a:b:c-d", "b"),
            ("a:b-c:d", "b"),
            ("no_colon", "no_colon"),
            ("no-colon", "no-colon"),
            ("", ""),
            (":", ""),
            ("x:", ""),
            (":-5", ""),
            ("1:-5-7", ""),
            ("1:5-", "5"),
            ("1::5", ""),
        ];
        for (location, start) in cases {
            assert_eq!(
                start_of(location.as_bytes()),
                start.as_bytes(),
                "{location:?}"
            );
        }
    }

    #[test]
    fn two_source_values_mark_a_reference_only_row() {
        assert!(is_reference_only(Some(b"only_perl")));
        assert!(is_reference_only(Some(b"missing_in_rust")));
        for other in [
            &b"only_rust"[..],
            b"extra_in_rust",
            b"",
            b"Only_perl",
            b"only_perl ",
        ] {
            assert!(!is_reference_only(Some(other)));
        }
        assert!(!is_reference_only(None));
    }
}
