// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The parsers of the output formats and the reference side's text classes.

pub mod json;
pub mod parquet;
pub mod recoder;
pub mod text;
pub mod vcf;

use crate::classes::{CONSERVATION_SHORT_NAME_SPELLING, PC_DESCRIPTION_TEXT};
use crate::consts::{CONSERVATION_KEY_SPELLING, MISSPELT_CONSERVATION_KEY};
use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::input::read_text;
use crate::manifest::Documented;
use std::borrow::Cow;
use std::path::Path;
use std::sync::LazyLock;

pub const PC_DESCRIPTIONS: &[(&str, &str)] = &[
    (
        "Percentage of input variant covered by reference variant",
        "Percentage of the annotation record covered by the input variant",
    ),
    (
        "Percentage of minimum reciprocal overlap between input variant and reference variant",
        "Percentage of the input variant covered by the annotation record",
    ),
];

pub const PC_FIELD_SUFFIX: &str = "_PC";

/// `CONSERVATION_KEY_SPELLING` as the `regex` crate spells a replacement.
static CONSERVATION_KEY_REPLACEMENT: LazyLock<String> =
    LazyLock::new(|| CONSERVATION_KEY_SPELLING.replace("\\1", "${1}"));

/// The key a header line describes: `## <key> : <text>` in the default and tab formats
/// (joined at VEP_output.pm:159, :164 and :167 and at Tab.pm:153; a custom source's pairs
/// of key and text reach VEP_output.pm:167 and Tab.pm:149 from `get_custom_headers`,
/// OutputFactory.pm:2444-2512), `##INFO=<ID=<key>,...>` in the VCF; `None` for any
/// other line.
pub fn header_key(line: &str) -> Option<String> {
    header_key_slice(line).map(str::to_string)
}

fn header_key_slice(line: &str) -> Option<&str> {
    if let Some(described) = line.strip_prefix("## ") {
        return described.split_once(" : ").map(|(key, _)| key);
    }
    let named = line.strip_prefix("##INFO=<ID=")?;
    Some(named.split_once(',').map_or(named, |(key, _)| key))
}

/// `text` with every misspelt conservation key spelt as the candidate spells it, and the
/// count of the keys respelt.
fn conservation_keys_respelt(text: &str) -> (Cow<'_, str>, i64) {
    let mut count = 0;
    let respelt = MISSPELT_CONSERVATION_KEY.replace_all(text, |found: &regex::Captures<'_>| {
        count += 1;
        let mut spelling = String::new();
        found.expand(&CONSERVATION_KEY_REPLACEMENT, &mut spelling);
        spelling
    });
    (respelt, count)
}

/// A header line describing a `<key>_PC` field with one of the two `PC_DESCRIPTIONS`,
/// carrying the candidate's text; `None` for any other line.
fn pc_description_respelt(line: &str) -> Option<String> {
    if !header_key_slice(line)?.ends_with(PC_FIELD_SUFFIX) {
        return None;
    }
    PC_DESCRIPTIONS
        .iter()
        .find(|(perl_text, _)| line.contains(perl_text))
        .map(|(perl_text, candidate_text)| line.replacen(perl_text, candidate_text, 1))
}

/// The reference output with the text classes of `documented.masks` applied before it is
/// parsed: under `conservation_short_name_spelling` every `phlyoP<N>way` key spelt
/// `phyloP<N>way`, under `pc_description_text` every header line describing a `<key>_PC`
/// field with one of the two `PC_DESCRIPTIONS` carrying the candidate's text. A rewrite is
/// counted under its class (a key per occurrence, a header line per line); the entries read
/// no header line, so the second class changes what a header comparison would read and
/// nothing the entries compare.
pub fn reference_text_masked(text: &str, documented: &mut Documented) -> String {
    let mut text = Cow::Borrowed(text);
    if documented.masks.contains(CONSERVATION_SHORT_NAME_SPELLING) {
        let (respelt, count) = conservation_keys_respelt(&text);
        if count > 0 {
            documented
                .applied
                .add(CONSERVATION_SHORT_NAME_SPELLING, count);
            text = Cow::Owned(respelt.into_owned());
        }
    }
    if documented.masks.contains(PC_DESCRIPTION_TEXT) {
        let mut respelt = String::with_capacity(text.len());
        for (at, line) in text.split('\n').enumerate() {
            if at > 0 {
                respelt.push('\n');
            }
            match pc_description_respelt(line) {
                Some(line) => {
                    documented.applied.add(PC_DESCRIPTION_TEXT, 1);
                    respelt.push_str(&line);
                }
                None => respelt.push_str(line),
            }
        }
        text = Cow::Owned(respelt);
    }
    text.into_owned()
}

/// `reference_text_masked` over the lines of a text, as `text_lines` gives them, in place.
///
/// Both rewrites are local to a line: the misspelt key holds no line end, and a
/// description is replaced within its header line. The lines of a text without a `\r`,
/// joined again, are therefore `reference_text_masked` of that text, with the same counts.
pub fn reference_lines_masked(lines: &mut [String], documented: &mut Documented) {
    if documented.masks.contains(CONSERVATION_SHORT_NAME_SPELLING) {
        let mut count = 0;
        for line in lines.iter_mut() {
            let (respelt, found) = conservation_keys_respelt(line);
            if found > 0 {
                count += found;
                *line = respelt.into_owned();
            }
        }
        if count > 0 {
            documented
                .applied
                .add(CONSERVATION_SHORT_NAME_SPELLING, count);
        }
    }
    if documented.masks.contains(PC_DESCRIPTION_TEXT) {
        for line in lines.iter_mut() {
            if let Some(respelt) = pc_description_respelt(line) {
                documented.applied.add(PC_DESCRIPTION_TEXT, 1);
                *line = respelt;
            }
        }
    }
}

/// Parses the `--tab` output `path` a Parquet output is compared against, in the canonical
/// form read beside that output: a column `parquet` stores as strings stays as written,
/// where the `tab-canonical` format of [`parse_format`] reprints a column by its name
/// alone. `rewrite` is applied to the canonical text before it is parsed.
pub fn parse_tab_beside_parquet(
    path: &Path,
    parquet: &Path,
    rewrite: Option<&mut dyn FnMut(&str) -> String>,
) -> Result<Vec<Entry>> {
    let mut text = parquet::tab_beside_parquet(path, parquet)?;
    if let Some(rewrite) = rewrite {
        text = rewrite(&text);
    }
    text::parse_tab(&text)
}

/// Parses one output file in `fmt`; VCF and JSON need the input records. `parquet`
/// reads a Parquet output, and `tab-canonical` a `--tab` output in the canonical form a
/// Parquet output is compared against. `rewrite` is applied to the text before it is parsed
/// (the reference side's text classes, `reference_text_masked`), to the adapter's tab text
/// for the two Parquet forms; `None` leaves the text as it is read.
pub fn parse_format(
    fmt: &str,
    path: &Path,
    inputs: Option<&[String]>,
    rewrite: Option<&mut dyn FnMut(&str) -> String>,
) -> Result<Vec<Entry>> {
    let mut text = match fmt {
        "parquet" | "tab-canonical" => parquet::parquet_as_tab(path, true)?,
        _ => read_text(path)?,
    };
    if let Some(rewrite) = rewrite {
        text = rewrite(&text);
    }
    match fmt {
        "default" => return text::parse_default(&text),
        "tab" | "parquet" | "tab-canonical" => return text::parse_tab(&text),
        _ => {}
    }
    let Some(inputs) = inputs else {
        return Err(Error::comparator(format!(
            "--format {fmt} needs --input, the input VCF"
        )));
    };
    match fmt {
        "vcf" => vcf::parse_vcf(&text, inputs),
        _ => json::parse_json(&text, inputs),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECIPROCAL: &str =
        "Percentage of minimum reciprocal overlap between input variant and reference variant";
    const PLAIN: &str = "Percentage of input variant covered by reference variant";

    /// A `Documented` of no manifest whose mask set names `classes`.
    fn documented_with(classes: &[&str]) -> Documented {
        let masks = classes.iter().map(|class| class.to_string()).collect();
        Documented::from_manifest(None, &masks, None, None, None).unwrap()
    }

    #[test]
    fn the_lines_of_a_text_are_rewritten_and_counted_as_the_text_is() {
        // A description line stands first: the conservation keys are still counted first.
        let text = format!(
            "## a_PC : {PLAIN} from f\nphlyoP7way=1;phlyoP20way=2\n## phlyoP7way_PC : {RECIPROCAL}\nrow\n"
        );
        let both = [PC_DESCRIPTION_TEXT, CONSERVATION_SHORT_NAME_SPELLING];
        let mut whole = documented_with(&both);
        let masked = reference_text_masked(&text, &mut whole);
        assert_eq!(
            masked,
            "## a_PC : Percentage of the annotation record covered by the input variant from f\n\
             phyloP7way=1;phyloP20way=2\n\
             ## phyloP7way_PC : Percentage of the input variant covered by the annotation record\n\
             row\n"
        );
        let counted = [
            (CONSERVATION_SHORT_NAME_SPELLING, 3),
            (PC_DESCRIPTION_TEXT, 2),
        ];
        assert_eq!(whole.applied.iter().collect::<Vec<_>>(), counted);
        let mut lines = crate::input::text_lines(&text);
        let mut by_line = documented_with(&both);
        reference_lines_masked(&mut lines, &mut by_line);
        assert_eq!(lines, crate::input::text_lines(&masked));
        assert_eq!(by_line.applied.iter().collect::<Vec<_>>(), counted);
        // A class the masks do not name rewrites nothing and counts nothing.
        let mut unnamed = documented_with(&[]);
        assert_eq!(reference_text_masked(&text, &mut unnamed), text);
        let mut lines = crate::input::text_lines(&text);
        reference_lines_masked(&mut lines, &mut unnamed);
        assert_eq!(lines, crate::input::text_lines(&text));
        assert!(unnamed.applied.is_empty());
    }

    #[test]
    fn a_header_line_names_its_key_in_the_text_shape_and_in_the_vcf_shape() {
        let cases = [
            (
                "## IMPACT : Subjective impact classification",
                Some("IMPACT"),
            ),
            ("## a : b : c", Some("a")),
            ("##  : ", Some("")),
            ("## a :b", None),
            ("## a", None),
            ("## ", None),
            ("##a : b", None),
            ("## a : b\r", Some("a")),
            ("##INFO=<ID=CSQ,Number=.,Type=String", Some("CSQ")),
            ("##INFO=<ID=CSQ", Some("CSQ")),
            ("##INFO=<ID=", Some("")),
            ("##INFO=<ID=,x", Some("")),
            ("##INFO=<ID=a : b,c", Some("a : b")),
            ("## INFO=<ID=a,b", None),
            ("##INFO=<id=a,b", None),
            ("#INFO=<ID=a,b", None),
            ("#Uploaded_variation\tLocation", None),
            ("21_5_A/T\t21:5", None),
            ("", None),
        ];
        for (line, key) in cases {
            assert_eq!(header_key(line).as_deref(), key, "{line:?}");
        }
    }

    #[test]
    fn every_misspelt_conservation_key_is_respelt_and_counted() {
        let (text, count) = conservation_keys_respelt("phlyoP20way=1.5;phlyoP470way=-0.2");
        assert_eq!(
            (text.as_ref(), count),
            ("phyloP20way=1.5;phyloP470way=-0.2", 2)
        );
        // A key needs its number; the digits of any script are a number.
        let (text, count) = conservation_keys_respelt("phlyoPway phlyoP\u{663}way phyloP7way");
        assert_eq!(
            (text.as_ref(), count),
            ("phlyoPway phyloP\u{663}way phyloP7way", 1)
        );
        // A key holds no line end, so a text and its lines are respelt alike.
        let (text, count) = conservation_keys_respelt("phlyoP1\nway\nphlyoP\n1way");
        assert_eq!((text.as_ref(), count), ("phlyoP1\nway\nphlyoP\n1way", 0));
    }

    #[test]
    fn a_pc_header_line_takes_the_candidates_text_for_the_first_description_it_holds() {
        let respelt = |line: &str| pc_description_respelt(line);
        assert_eq!(
            respelt(&format!("## regions_PC : {RECIPROCAL} from [PATH]/regions.bed.gz")).unwrap(),
            "## regions_PC : Percentage of the input variant covered by the annotation record from [PATH]/regions.bed.gz"
        );
        assert_eq!(
            respelt(&format!("##INFO=<ID=plain_PC,Number=.,Type=String,Description=\"{PLAIN} from f\">")).unwrap(),
            "##INFO=<ID=plain_PC,Number=.,Type=String,Description=\"Percentage of the annotation record covered by the input variant from f\">"
        );
        // One description a line, the first of the two in the table's order, once.
        assert_eq!(
            respelt(&format!("## a_PC : {RECIPROCAL} then {PLAIN} and {PLAIN}")).unwrap(),
            format!("## a_PC : {RECIPROCAL} then Percentage of the annotation record covered by the input variant and {PLAIN}")
        );
        for line in [
            format!("## regions_min : {PLAIN}"),
            format!("## other : {PLAIN}"),
            format!("## a_PCx : {PLAIN}"),
            format!("## a_PC: {PLAIN}"),
            format!("#a_PC : {PLAIN}"),
            "## a_PC : another text".to_string(),
            format!("##INFO=<ID=x_PC {RECIPROCAL}"),
            format!("row\tregions_PC={PLAIN}"),
        ] {
            assert_eq!(respelt(&line), None, "{line:?}");
        }
    }
}
