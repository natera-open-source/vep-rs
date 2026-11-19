// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cell contract: `parse` from a `--tab` cell to a typed value and
//! `to_tab_text` back, per [`ColumnType`] (the module documentation of
//! [`super`] states the grammar of each type).

use super::{BoolSpelling, ColumnType};

/// The members of a `Positions` cell: a span and the `--total_length` total.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Positions {
    pub start: Option<i32>,
    pub end: Option<i32>,
    pub total: Option<i32>,
}

/// A stored cell value; NULL is the absence of one.
#[derive(Debug, Clone, PartialEq)]
pub enum CellValue {
    Text(String),
    Int32(i32),
    Int64(i64),
    Double(f64),
    Boolean(bool),
    /// A NULL element is a token that was not a number.
    DoubleList(Vec<Option<f64>>),
    /// A NULL element is a token that was neither `0` nor `1`.
    BooleanList(Vec<Option<bool>>),
    TextList(Vec<String>),
    Positions(Positions),
}

/// `parse(text)`: the value a non-NULL tab cell stores under `column_type`, or
/// `None` for a token outside the type's grammar, which is stored as NULL. The
/// caller decides NULL for a `-` or empty cell first.
pub fn parse_cell(column_type: ColumnType, text: &str) -> Option<CellValue> {
    Some(match column_type {
        ColumnType::Text => CellValue::Text(text.to_string()),
        ColumnType::Int32 => CellValue::Int32(text.parse().ok()?),
        ColumnType::Double => CellValue::Double(text.parse().ok()?),
        ColumnType::Boolean(spelling) => CellValue::Boolean(parse_bool(spelling, text)?),
        ColumnType::DoubleList => {
            CellValue::DoubleList(text.split(',').map(|t| t.parse().ok()).collect())
        }
        ColumnType::BooleanList => {
            CellValue::BooleanList(text.split(',').map(parse_flag).collect())
        }
        ColumnType::TextList => CellValue::TextList(text.split(',').map(str::to_string).collect()),
        ColumnType::Positions => CellValue::Positions(parse_positions(text)?),
    })
}

pub(super) fn parse_bool(spelling: BoolSpelling, text: &str) -> Option<bool> {
    match (spelling, text) {
        (BoolSpelling::Yes, "YES") | (BoolSpelling::One, "1") | (BoolSpelling::YesNo, "Y") => {
            Some(true)
        }
        (BoolSpelling::YesNo, "N") => Some(false),
        _ => None,
    }
}

/// A `LIST<BOOLEAN>` element: `0` is false and `1` is true.
pub(super) fn parse_flag(token: &str) -> Option<bool> {
    match token {
        "0" => Some(false),
        "1" => Some(true),
        _ => None,
    }
}

/// `[span][/total]`, the span `start`, `start-end`, `start-?`, `?-end` or `-`.
pub(super) fn parse_positions(text: &str) -> Option<Positions> {
    let (span, total) = match text.split_once('/') {
        Some((span, total)) => (span, Some(total.parse::<i32>().ok()?)),
        None => (text, None),
    };
    let member = |m: &str| -> Option<Option<i32>> {
        if m == "?" {
            Some(None)
        } else {
            m.parse::<i32>().ok().map(Some)
        }
    };
    let (start, end) = if span == "-" {
        (None, None)
    } else if let Some((start, end)) = span.split_once('-') {
        (member(start)?, member(end)?)
    } else {
        let position = member(span)?;
        (position, position)
    };
    Some(Positions { start, end, total })
}

/// `to_tab_text(value)`: the tab cell of a stored value; a NULL cell, which has
/// no value, prints `-`.
pub fn to_tab_text(column_type: ColumnType, value: &CellValue) -> String {
    match value {
        CellValue::Text(text) => text.clone(),
        CellValue::Int32(v) => v.to_string(),
        CellValue::Int64(v) => v.to_string(),
        CellValue::Double(v) => shortest_decimal(*v),
        CellValue::Boolean(v) => match column_type {
            ColumnType::Boolean(spelling) => bool_text(spelling, *v).to_string(),
            _ => flag_text(*v).to_string(),
        },
        CellValue::DoubleList(values) => join(values.iter().map(|v| match v {
            Some(v) => shortest_decimal(*v),
            None => null_tab_text().to_string(),
        })),
        CellValue::BooleanList(values) => join(values.iter().map(|v| match v {
            Some(v) => flag_text(*v).to_string(),
            None => null_tab_text().to_string(),
        })),
        CellValue::TextList(values) => values.join(","),
        CellValue::Positions(positions) => positions_text(positions),
    }
}

/// The tab cell of a NULL value.
pub const fn null_tab_text() -> &'static str {
    "-"
}

fn bool_text(spelling: BoolSpelling, value: bool) -> &'static str {
    match (spelling, value) {
        (BoolSpelling::Yes, true) => "YES",
        (BoolSpelling::One, true) => "1",
        (BoolSpelling::YesNo, true) => "Y",
        (BoolSpelling::YesNo, false) => "N",
        (BoolSpelling::Yes | BoolSpelling::One, false) => null_tab_text(),
    }
}

fn flag_text(value: bool) -> &'static str {
    if value {
        "1"
    } else {
        "0"
    }
}

fn positions_text(positions: &Positions) -> String {
    let member = |m: Option<i32>| match m {
        Some(v) => v.to_string(),
        None => "?".to_string(),
    };
    let mut text = match (positions.start, positions.end) {
        (None, None) => "-".to_string(),
        (Some(start), Some(end)) if start == end => start.to_string(),
        (start, end) => format!("{}-{}", member(start), member(end)),
    };
    if let Some(total) = positions.total {
        text.push('/');
        text.push_str(&total.to_string());
    }
    text
}

fn join(parts: impl Iterator<Item = String>) -> String {
    parts.collect::<Vec<_>>().join(",")
}

/// The shortest decimal that reads back to the same value: `0.0950` is `0.095`
/// and `1.0` is `1`.
fn shortest_decimal(value: f64) -> String {
    value.to_string()
}

/// The tab cell as the Parquet round trip prints it: a DOUBLE or
/// `LIST<DOUBLE>` cell with each number as its shortest round-trip decimal, an
/// INT32 cell as its integer (a `TSL` of `007` prints `7`) or NULL when the
/// text is no INT32, every other cell unchanged.
pub fn canonical_tab_text(column_type: ColumnType, text: &str) -> String {
    if text == null_tab_text() {
        return text.to_string();
    }
    let canonical_number = |t: &str| match t.parse::<f64>() {
        Ok(v) => shortest_decimal(v),
        Err(_) => t.to_string(),
    };
    match column_type {
        ColumnType::Double => canonical_number(text),
        ColumnType::DoubleList => join(text.split(',').map(canonical_number)),
        ColumnType::Int32 => match text.parse::<i32>() {
            Ok(v) => v.to_string(),
            Err(_) => null_tab_text().to_string(),
        },
        _ => text.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::super::column_type;
    use super::*;

    /// One tab cell per type and shape the tab format writes, keyed by the
    /// column that writes it.
    const FIXTURES: &[(&str, &str)] = &[
        ("Feature", "ENST00000456328"),
        ("Location", "21:25000100-25000101"),
        ("Existing_variation", "rs1,COSV12345"),
        ("STRAND", "-1"),
        ("DISTANCE", "4821"),
        ("ALLELE_NUM", "2"),
        ("TSL", "5"),
        ("TSL", "007"),
        ("TSL", "00"),
        ("HGVS_OFFSET", "-3"),
        ("OverlapBP", "272"),
        ("MOTIF_POS", "12"),
        ("MAX_AF", "0.0950"),
        ("MAX_AF", "1"),
        ("MAX_AF", "9.8e-05"),
        ("MOTIF_SCORE_CHANGE", "-0.012"),
        ("MOTIF_SCORE_CHANGE", "0.100"),
        ("CANONICAL", "YES"),
        ("PICK", "1"),
        ("GENCODE_PRIMARY", "1"),
        ("HIGH_INF_POS", "Y"),
        ("HIGH_INF_POS", "N"),
        ("AF", "0.0002"),
        ("gnomADe_AF", "0.0950,0.1"),
        ("gnomADg_NFE_AF", "0"),
        ("SOMATIC", "0,1"),
        ("PHENO", "1"),
        ("MANE", "MANE_Select,MANE_Plus_Clinical"),
        ("Consequence", "missense_variant,splice_region_variant"),
        ("CLIN_SIG", "benign;likely_benign,pathogenic"),
        ("PUBMED", "12345,67890"),
        ("FLAGS", "cds_start_NF,cds_end_NF"),
        ("TRANSCRIPTION_FACTORS", "CTCF,ZNF143"),
        ("EXON", "3/12"),
        ("EXON", "3-4/12"),
        ("INTRON", "7/11"),
        ("cDNA_position", "44"),
        ("cDNA_position", "44-45"),
        ("cDNA_position", "44/1200"),
        ("cDNA_position", "-/1200"),
        ("cDNA_position", "44-45/1200"),
        ("CDS_position", "146-?"),
        ("Protein_position", "?-52"),
        ("Protein_position", "49-?/300"),
    ];

    #[test]
    fn every_typed_cell_inverts_to_its_tab_text() {
        for (column, text) in FIXTURES {
            let column_type = column_type(column);
            let value =
                parse_cell(column_type, text).unwrap_or_else(|| panic!("{column} {text:?} parses"));
            assert_eq!(
                to_tab_text(column_type, &value),
                canonical_tab_text(column_type, text),
                "{column} {text:?}"
            );
        }
    }

    #[test]
    fn canon_rewrites_a_double_and_an_int32_cell() {
        assert_eq!(canonical_tab_text(ColumnType::Double, "0.0950"), "0.095");
        assert_eq!(canonical_tab_text(ColumnType::Double, "1.0"), "1");
        assert_eq!(
            canonical_tab_text(ColumnType::Double, "9.8e-05"),
            "0.000098"
        );
        assert_eq!(canonical_tab_text(ColumnType::Double, "-"), "-");
        assert_eq!(
            canonical_tab_text(ColumnType::DoubleList, "0.0950,0.1,1"),
            "0.095,0.1,1"
        );
        assert_eq!(canonical_tab_text(ColumnType::Text, "0.0950"), "0.0950");
        assert_eq!(canonical_tab_text(ColumnType::Int32, "007"), "7");
        assert_eq!(canonical_tab_text(ColumnType::Int32, "00"), "0");
        assert_eq!(canonical_tab_text(ColumnType::Int32, "-1"), "-1");
        assert_eq!(canonical_tab_text(ColumnType::Int32, "2147483648"), "-");
        assert_eq!(canonical_tab_text(ColumnType::Int32, "-"), "-");
    }

    #[test]
    fn a_token_outside_the_grammar_is_null() {
        assert_eq!(parse_cell(ColumnType::Int32, "plus"), None);
        assert_eq!(parse_cell(ColumnType::Int32, "1.5"), None);
        assert_eq!(parse_cell(ColumnType::Double, "n/a"), None);
        assert_eq!(
            parse_cell(ColumnType::Boolean(BoolSpelling::Yes), "1"),
            None
        );
        assert_eq!(
            parse_cell(ColumnType::Boolean(BoolSpelling::One), "YES"),
            None
        );
        assert_eq!(
            parse_cell(ColumnType::Boolean(BoolSpelling::YesNo), "YES"),
            None
        );
        assert_eq!(parse_cell(ColumnType::Positions, "a-b"), None);
        assert_eq!(parse_cell(ColumnType::Positions, "44/many"), None);
        assert_eq!(
            parse_cell(ColumnType::DoubleList, "0.1,x"),
            Some(CellValue::DoubleList(vec![Some(0.1), None]))
        );
        assert_eq!(
            parse_cell(ColumnType::BooleanList, "1,2"),
            Some(CellValue::BooleanList(vec![Some(true), None]))
        );
    }

    #[test]
    fn positions_members_follow_the_span_grammar() {
        let positions = |text: &str| match parse_cell(ColumnType::Positions, text) {
            Some(CellValue::Positions(p)) => p,
            other => panic!("{text}: {other:?}"),
        };
        let p = |start, end, total| Positions { start, end, total };
        assert_eq!(positions("44"), p(Some(44), Some(44), None));
        assert_eq!(positions("44-45"), p(Some(44), Some(45), None));
        assert_eq!(positions("44/1200"), p(Some(44), Some(44), Some(1200)));
        assert_eq!(positions("-/1200"), p(None, None, Some(1200)));
        assert_eq!(positions("146-?"), p(Some(146), None, None));
        assert_eq!(positions("?-52"), p(None, Some(52), None));
    }
}
