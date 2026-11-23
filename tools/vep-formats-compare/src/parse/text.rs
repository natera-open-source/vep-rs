// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The default and tab parsers.

use crate::consts::JSON_OBJECT_FIELDS;
use crate::entry::{Entry, Fields};
use crate::error::{Error, Result};
use crate::input::{canonical_json_text_under, dedup_keys, line_slices, CANONICAL_CALLERS};
use crate::shared;
use std::borrow::Cow;

/// The call depth a parser runs at where a comparison parses an output. A nested value
/// is read as far as the interpreter's recursion limit (`pyjson::RECURSION_LIMIT`) less
/// the depth its reading runs at: `canonical_json_text` runs one call deeper than
/// `parse_default` (`CANONICAL_CALLERS`) and two deeper than `parse_tab`, so the tab
/// format refuses an object half a level sooner than the default format.
pub const COMPARISON_CALLERS: usize = CANONICAL_CALLERS - 1;

/// The default format: the columns of the `#` header line, the `Extra` column's
/// `key=value` pairs becoming fields of their own while `Extra` keeps its raw text, a
/// `JSON_OBJECT_FIELDS` value in canonical form in both.
pub fn parse_default(text: &str) -> Result<Vec<Entry>> {
    parse_default_under(text, COMPARISON_CALLERS)
}

/// `parse_default` called with `callers` frames of the interpreter's call stack below it
/// (`canonical_json_text_under`).
pub fn parse_default_under(text: &str, callers: usize) -> Result<Vec<Entry>> {
    parse_rows(line_slices(text), callers + 1)
}

/// `parse_default` over the lines of a text, as `text_lines` gives them.
pub fn parse_default_lines(lines: &[String]) -> Result<Vec<Entry>> {
    parse_rows(lines.iter().map(String::as_str), COMPARISON_CALLERS + 1)
}

/// The tab format: every column named in the `#` header line.
pub fn parse_tab(text: &str) -> Result<Vec<Entry>> {
    parse_tab_under(text, COMPARISON_CALLERS)
}

/// `parse_tab` called with `callers` frames of the interpreter's call stack below it:
/// `parse_default` stands one frame above it.
pub fn parse_tab_under(text: &str, callers: usize) -> Result<Vec<Entry>> {
    parse_default_under(text, callers + 1)
}

/// `parse_tab` over the lines of a text, as `text_lines` gives them.
pub fn parse_tab_lines(lines: &[String]) -> Result<Vec<Entry>> {
    parse_rows(lines.iter().map(String::as_str), COMPARISON_CALLERS + 2)
}

/// The entries of `lines`; `callers` counts the frames below `canonical_json_text`.
fn parse_rows<'a>(lines: impl Iterator<Item = &'a str>, callers: usize) -> Result<Vec<Entry>> {
    let mut cols: Vec<&str> = Vec::new();
    // The column whose pairs become fields: of several named `Extra`, the last.
    let mut extra_at = None;
    let mut vals: Vec<&str> = Vec::new();
    let mut entries = Vec::new();
    for line in lines {
        if let Some(header) = line.strip_prefix('#') {
            if !header.starts_with('#') {
                cols = header.split('\t').collect();
                extra_at = cols.iter().rposition(|col| *col == "Extra");
            }
            continue;
        }
        vals.clear();
        vals.extend(line.split('\t'));
        if vals.len() != cols.len() {
            return Err(Error::comparator(format!("column count on line: {line}")));
        }
        let pairs = extra_at
            .map(|at| vals[at])
            .filter(|extra| *extra != "-")
            .map_or(0, |extra| extra.split(';').count());
        let mut fields = Fields::with_capacity(cols.len() + pairs);
        for (col, val) in cols.iter().zip(&vals) {
            fields.insert(col.to_string(), val.to_string());
        }
        for name in JSON_OBJECT_FIELDS {
            if let Some(value) = fields.get_mut(*name) {
                *value = canonical_json_text_under(value, callers)?;
            }
        }
        if fields.get("Extra").is_some_and(|extra| extra != "-") {
            let extra = std::mem::take(&mut fields["Extra"]);
            // The pairs whose value is written in canonical form, by their place.
            let mut canonical: Vec<(usize, String)> = Vec::new();
            for (at, kv) in extra.split(';').enumerate() {
                match kv.split_once('=') {
                    Some((k, v)) if JSON_OBJECT_FIELDS.contains(&k) => {
                        let v = canonical_json_text_under(v, callers)?;
                        canonical.push((at, format!("{k}={v}")));
                        fields.insert(k.to_string(), v);
                    }
                    Some((k, v)) => {
                        fields.insert(k.to_string(), v.to_string());
                    }
                    None => {
                        fields.insert(kv.to_string(), String::new());
                    }
                }
            }
            fields["Extra"] = if canonical.is_empty() {
                extra
            } else {
                let mut canonical = canonical.into_iter().peekable();
                let pairs: Vec<Cow<'_, str>> = extra
                    .split(';')
                    .enumerate()
                    .map(
                        |(at, kv)| match canonical.next_if(|(place, _)| *place == at) {
                            Some((_, kv)) => Cow::Owned(kv),
                            None => Cow::Borrowed(kv),
                        },
                    )
                    .collect();
                pairs.join(";")
            };
        }
        let field = |name: &str| fields.get(name).map_or("", String::as_str);
        let (location, allele, feature) = (field("Location"), field("Allele"), field("Feature"));
        let (mut indel, mut insertion, mut reference) = (None, None, None);
        let uploaded = field("Uploaded_variation");
        if let Some((reference_allele, alternate_alleles)) = shared::uploaded_alleles(uploaded) {
            let alternates: Vec<&str> = alternate_alleles.iter().map(String::as_str).collect();
            indel = shared::allele_is_indel(&reference_allele, &alternates, allele, false);
            insertion = Some(shared::reference_allele_is_empty(
                &reference_allele,
                &alternates,
                false,
            ));
            if let Some(position) = shared::uploaded_position(uploaded)? {
                reference = shared::replaced_reference(
                    position,
                    &reference_allele,
                    &alternates,
                    Some(location),
                    allele,
                )?;
            }
        }
        let key = [location, allele, feature, field("Feature_type")].join("|");
        let (location, allele, feature, consequence) = (
            location.to_string(),
            allele.to_string(),
            feature.to_string(),
            field("Consequence").to_string(),
        );
        entries.push(Entry {
            key,
            record: None,
            location: Some(location),
            allele,
            feature,
            consequence,
            fields,
            indel,
            insertion,
            reference,
            line: None,
            allele_location: None,
        });
    }
    dedup_keys(&mut entries);
    Ok(entries)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::text_lines;

    const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

    fn row(uploaded: &str, location: &str, allele: &str, feature: &str, extra: &str) -> String {
        format!("{uploaded}\t{location}\t{allele}\tG1\t{feature}\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\t{extra}\n")
    }

    fn field<'a>(entry: &'a Entry, name: &str) -> Option<&'a str> {
        entry.fields.get(name).map(String::as_str)
    }

    fn names(entry: &Entry) -> Vec<&str> {
        entry.fields.keys().map(String::as_str).collect()
    }

    fn error_text(result: Result<Vec<Entry>>) -> String {
        match result.expect_err("the comparator's error") {
            Error::Comparator(text) => text,
        }
    }

    fn nested(depth: usize, leaf: &str) -> String {
        format!("{}{leaf}{}", "{\"a\":".repeat(depth), "}".repeat(depth))
    }

    #[test]
    fn the_extra_pairs_become_fields_after_the_columns_and_extra_keeps_its_text() {
        let text = format!(
            "{HEADER}{}",
            row(
                "rs1",
                "21:5",
                "C",
                "T1",
                "IMPACT=MODIFIER;STRAND=-1;FLAGS;A=b=c;E="
            )
        );
        let entries = parse_default(&text).unwrap();
        assert_eq!(entries.len(), 1);
        let entry = &entries[0];
        assert_eq!(entry.key, "21:5|C|T1|Transcript");
        assert_eq!(entry.record, None);
        assert_eq!(entry.location.as_deref(), Some("21:5"));
        assert_eq!((entry.allele.as_str(), entry.feature.as_str()), ("C", "T1"));
        assert_eq!(entry.consequence, "intron_variant");
        assert_eq!(
            (entry.indel, entry.insertion, entry.reference.as_deref()),
            (None, None, None)
        );
        assert_eq!((entry.line, entry.allele_location.as_deref()), (None, None));
        assert_eq!(
            field(entry, "Extra"),
            Some("IMPACT=MODIFIER;STRAND=-1;FLAGS;A=b=c;E=")
        );
        assert_eq!(field(entry, "IMPACT"), Some("MODIFIER"));
        assert_eq!(field(entry, "STRAND"), Some("-1"));
        assert_eq!(field(entry, "FLAGS"), Some(""));
        assert_eq!(field(entry, "A"), Some("b=c"));
        assert_eq!(field(entry, "E"), Some(""));
        assert_eq!(
            names(entry)[13..],
            ["Extra", "IMPACT", "STRAND", "FLAGS", "A", "E"]
        );
        // An Extra of `-` holds no pair; an empty one holds the pair of the empty name.
        let text = format!(
            "{HEADER}{}{}",
            row("rs1", "21:5", "C", "T1", "-"),
            row("rs1", "21:6", "C", "T1", "")
        );
        let entries = parse_default(&text).unwrap();
        assert_eq!(names(&entries[0]).len(), 14);
        assert_eq!(field(&entries[0], "Extra"), Some("-"));
        assert_eq!(names(&entries[1])[13..], ["Extra", ""]);
        assert_eq!(field(&entries[1], ""), Some(""));
    }

    #[test]
    fn a_json_object_field_is_canonical_as_a_column_and_as_an_extra_pair() {
        let perl = r#"{"type":"Allele","state":{"sequence":"C","type":"SequenceState"}}"#;
        let rust = r#"{"state":{"type":"SequenceState","sequence":"C"},"type":"Allele"}"#;
        let canonical = r#"{"state":{"sequence":"C","type":"SequenceState"},"type":"Allele"}"#;
        let extra = |value: &str| format!("IMPACT=MODIFIER;GA4GH_VRS={value};STRAND=1");
        let text = format!(
            "{HEADER}{}{}{}{}",
            row("rs1", "21:5", "C", "T1", &extra(perl)),
            row("rs1", "21:5", "C", "T2", &extra(rust)),
            row("rs1", "21:5", "C", "T3", &extra("-")),
            row(
                "rs1",
                "21:5",
                "C",
                "T4",
                "IMPACT=MODIFIER;GA4GH_VRS;STRAND=1"
            ),
        );
        let entries = parse_default(&text).unwrap();
        for entry in &entries[..2] {
            assert_eq!(field(entry, "GA4GH_VRS"), Some(canonical));
            assert_eq!(field(entry, "Extra"), Some(extra(canonical).as_str()));
            assert_eq!(field(entry, "STRAND"), Some("1"));
        }
        // A pair whose value is no object, and a bare key, are rebuilt as written.
        assert_eq!(
            field(&entries[2], "Extra"),
            Some("IMPACT=MODIFIER;GA4GH_VRS=-;STRAND=1")
        );
        assert_eq!(field(&entries[2], "GA4GH_VRS"), Some("-"));
        assert_eq!(
            field(&entries[3], "Extra"),
            Some("IMPACT=MODIFIER;GA4GH_VRS;STRAND=1")
        );
        assert_eq!(field(&entries[3], "GA4GH_VRS"), Some(""));
        let tab = format!("#Uploaded_variation\tLocation\tAllele\tGA4GH_VRS\nrs1\t21:5\tC\t{perl}\nrs1\t21:6\tC\t{rust}\nrs1\t21:7\tC\t-\n");
        let entries = parse_tab(&tab).unwrap();
        assert_eq!(field(&entries[0], "GA4GH_VRS"), Some(canonical));
        assert_eq!(field(&entries[1], "GA4GH_VRS"), Some(canonical));
        assert_eq!(field(&entries[2], "GA4GH_VRS"), Some("-"));
    }

    #[test]
    fn the_last_column_line_names_the_columns_and_a_row_of_another_count_is_refused() {
        let text = "## note\n#a\tb\n1\t2\n#b\ta\n3\t4\n";
        let entries = parse_default(text).unwrap();
        assert_eq!(names(&entries[0]), ["a", "b"]);
        assert_eq!(names(&entries[1]), ["b", "a"]);
        assert_eq!(field(&entries[1], "a"), Some("4"));
        assert_eq!(parse_default("").unwrap(), []);
        assert_eq!(parse_default(HEADER).unwrap(), []);
        // The column line `#` names one column, the empty name.
        let entries = parse_default("#\nx\n").unwrap();
        assert_eq!(names(&entries[0]), [""]);
        let refused = [
            ("#a\tb\n1\n", "column count on line: 1"),
            ("#a\tb\n1\t2\t3\n", "column count on line: 1\t2\t3"),
            ("#a\tb\n\n", "column count on line: "),
            ("1\t2\n#a\tb\n", "column count on line: 1\t2"),
            ("##a\tb\n1\t2\n", "column count on line: 1\t2"),
            ("#a\tb\n1\t2\n\n", "column count on line: "),
        ];
        for (text, error) in refused {
            assert_eq!(error_text(parse_default(text)), error, "{text:?}");
            assert_eq!(error_text(parse_tab(text)), error, "{text:?}");
        }
    }

    #[test]
    fn a_repeated_key_takes_a_counter_and_a_repeated_name_its_last_value() {
        let text = format!(
            "{HEADER}{}{}{}{}",
            row("rs1", "21:5", "C", "T1", "-"),
            row("rs1", "21:5", "C", "T1", "-"),
            row("rs1", "21:5", "C", "T2", "-"),
            row("rs1", "21:5", "C", "T1", "-"),
        );
        let keys: Vec<String> = parse_default(&text)
            .unwrap()
            .into_iter()
            .map(|e| e.key)
            .collect();
        assert_eq!(
            keys,
            [
                "21:5|C|T1|Transcript",
                "21:5|C|T1|Transcript#2",
                "21:5|C|T2|Transcript",
                "21:5|C|T1|Transcript#3"
            ]
        );
        // A column named twice keeps its first place and takes the later value.
        let entries = parse_default("#Location\tAllele\tLocation\n1:1\tA\t2:2\n").unwrap();
        assert_eq!(names(&entries[0]), ["Location", "Allele"]);
        assert_eq!(entries[0].key, "2:2|A||");
        // An Extra pair named as a column replaces the column's value, and so the key; a pair
        // named `Extra` is lost to the rebuilt text.
        let extra = "Extra=x;Location=9:9;Allele=G;Feature=T9;Feature_type=K;Consequence=c;A=1;A=2";
        let entries =
            parse_default(&format!("{HEADER}{}", row("rs1", "21:5", "C", "T1", extra))).unwrap();
        let entry = &entries[0];
        assert_eq!(entry.key, "9:9|G|T9|K");
        assert_eq!(entry.location.as_deref(), Some("9:9"));
        assert_eq!(entry.consequence, "c");
        assert_eq!(field(entry, "Extra"), Some(extra));
        assert_eq!(field(entry, "A"), Some("2"));
        assert_eq!(names(entry).len(), 15);
    }

    #[test]
    fn a_name_that_spells_its_alleles_gives_the_rows_shape_and_reference_bases() {
        let text = concat!(
            "#Uploaded_variation\tLocation\tAllele\n",
            "21_100_A/T\t21:100\tT\n",
            "21_100_A/AT\t21:100-101\tT\n",
            "21_100_AT/A\t21:101\t-\n",
            "rs100\t21:100\tT\n",
        );
        let entries = parse_default(text).unwrap();
        let read: Vec<_> = entries
            .iter()
            .map(|e| (e.indel, e.insertion, e.reference.as_deref()))
            .collect();
        assert_eq!(
            read,
            [
                (Some(false), Some(false), Some("A")),
                (Some(true), Some(true), Some("")),
                (Some(true), Some(false), Some("T")),
                (None, None, None),
            ]
        );
    }

    #[test]
    fn the_lines_of_a_text_parse_as_the_text() {
        let text = format!(
            "## note\r\n{}{}{}",
            HEADER.replace('\n', "\r\n"),
            row(
                "21_100_A/AT",
                "21:100-101",
                "T",
                "T1",
                "IMPACT=MODIFIER;GA4GH_VRS={\"b\":1,\"a\":2}"
            ),
            row("rs1", "21:5", "C", "T1", "IMPACT=MODIFIER").trim_end_matches('\n'),
        );
        let lines = text_lines(&text);
        let entries = parse_default(&text).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(field(&entries[0], "GA4GH_VRS"), Some(r#"{"a":2,"b":1}"#));
        assert_eq!(parse_default_lines(&lines).unwrap(), entries);
        assert_eq!(parse_tab(&text).unwrap(), entries);
        assert_eq!(parse_tab_lines(&lines).unwrap(), entries);
        assert_eq!(
            error_text(parse_default_lines(&["x".to_string()])),
            "column count on line: x"
        );
    }

    #[test]
    fn the_tab_parser_stands_one_frame_deeper_under_a_nested_object() {
        // An object the default parser still walks is one frame too deep under the tab
        // parser, which stands one call deeper.
        let text = format!("#Location\tGA4GH_VRS\n21:5\t{}\n", nested(495, "\"x\""));
        let lines = text_lines(&text);
        let refused =
            "a JSON object is nested too deep to read: {\"a\":{\"a\":{\"a\":{\"a\":{\"a\"...";
        assert_eq!(parse_default(&text).unwrap().len(), 1);
        assert_eq!(parse_default_lines(&lines).unwrap().len(), 1);
        assert_eq!(error_text(parse_tab(&text)), refused);
        assert_eq!(error_text(parse_tab_lines(&lines)), refused);
        assert_eq!(
            parse_tab_under(&text, COMPARISON_CALLERS - 1)
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            error_text(parse_default_under(&text, COMPARISON_CALLERS + 1)),
            refused
        );
        let extra = format!(
            "{HEADER}{}",
            row(
                "rs1",
                "21:5",
                "C",
                "T1",
                &format!("GA4GH_VRS={}", nested(1000, "1"))
            )
        );
        assert_eq!(error_text(parse_default(&extra)), refused);
    }
}
