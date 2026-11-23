// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What two readers of the regulatory GFF files the Ensembl site distributes
//! share: the annotator reads them at run time under `--regulatory_gff` and
//! the cache builder bakes them into a cache's regulatory shards, and both must
//! build the same feature from one record (Ensembl VEP's
//! `AnnotationSource/File/RegFeat.pm` at release/116.3). [`uri_unescape`] is
//! the decoding ensembl-io gives a GTF attribute value and a GFF3 type, for
//! every reader of either, [`decode_type`] the type column as text and
//! [`decode_attribute`] a GFF3 attribute value as text, for every reader of a
//! GFF3's ninth column. [`unindexed_coordinate`] names the start or end of a
//! line that no `tabix` indexes, which refuses a line of any type, and
//! [`index_columns`] is a line's columns as that index build ends them, at a
//! tab or at a NUL byte.

use std::borrow::Cow;

/// The GFF feature types a regulatory file yields (`RegFeat.pm:94-103`), each
/// a sequence ontology term printed as the row's BIOTYPE.
pub const REGULATORY_FEATURE_TYPES: [&str; 8] = [
    "promoter",
    "enhancer",
    "CTCF_binding_site",
    "open_chromatin_region",
    "promoter_flanking_region",
    "TF_binding_site",
    "regulatory_region",
    "EMAR",
];

/// The one type a motif file yields (`RegFeat.pm:118-120`).
pub const MOTIF_FEATURE_TYPES: [&str; 1] = ["TF_binding_site"];

/// The types whose `extended_start` and `extended_end` attributes
/// `--extended_promoters` widens (`RegFeat.pm:106-108`).
pub const EXTENDED_BOUNDS_TYPES: [&str; 1] = ["promoter"];

/// The text a record of an unsupported type is warned about with, once per
/// type and file (`RegFeat.pm:318-322`).
pub fn unsupported_type_warning(feature_type: &str, file: &str) -> String {
    format!("Ignoring '{feature_type}' feature type from {file}; not a supported regulatory feature type")
}

/// The text a promoter's `extended_start` or `extended_end` draws under the
/// extended-promoter option when its value is not a whole number: the bound
/// is ignored and the promoter keeps its own. Ensembl VEP compares the text
/// as a number (`RegFeat.pm:363-364`), so a value such as `x` reads as 0 there
/// and widens the promoter to the chromosome's start.
pub fn ignored_bound_warning(
    attribute: &str,
    value: &str,
    feature_type: &str,
    stable_id: &str,
    file: &str,
) -> String {
    format!("Ignoring {attribute} '{value}' of {feature_type} {stable_id} in {file}; not a whole number")
}

/// The text two distinct records of one file sharing an identifier draw, once
/// per identifier and file, Ensembl VEP's own words (`RegFeat.pm:487-491`).
pub fn duplicate_identifier_warning(stable_id: &str, file: &str) -> String {
    format!("Duplicate regulatory feature identifier '{stable_id}' in {file}; identifiers must be unique or annotations will be lost")
}

/// The identifier of a regulatory record (`RegFeat.pm:575-584`): its `ID`,
/// else its `Name`, an empty value passed over, else the chromosome as the
/// file spells it (`AnnotationSource.pm:179,198` maps the input's name onto
/// the file's before the region reaches `RegFeat.pm:372,583`) with the
/// feature's coordinates and type, so a record without either still yields a
/// named feature that two reads agree on.
pub fn feature_id<'a>(
    id: Option<&'a str>,
    name: Option<&'a str>,
    chr: &str,
    start: u64,
    end: u64,
    feature_type: &str,
) -> Cow<'a, str> {
    let named = [id, name]
        .into_iter()
        .flatten()
        .find(|value| !value.is_empty());
    match named {
        Some(value) => Cow::Borrowed(value),
        None => Cow::Owned(format!("{chr}:{start}-{end}:{feature_type}")),
    }
}

/// A GFF3 attribute value as text, its escapes decoded as URI::Escape's
/// `uri_unescape` decodes them, the decoding the format names for its ninth
/// column: the value's bytes with every `%` followed by two hexadecimal
/// digits, in either letter case, replaced by the byte they spell, in one pass
/// (`%252C` reads `%2C`; a `%` not so followed stays written), then read as
/// text, a byte sequence that is no UTF-8 once decoded standing as the
/// replacement character. A UTF-8 value without `%` is borrowed as it is.
/// Ensembl VEP's parser stores the value as written (ensembl-io
/// `GFF3::get_attributes`), so only a value holding an escape, which the files
/// Ensembl distributes do not, reads otherwise there.
pub fn decode_attribute(value: &[u8]) -> Cow<'_, str> {
    unescape_to_text(value)
}

/// The text a build-time reader of a GFF3 or GTF refuses a line with when an
/// attribute value's escapes spell bytes that are not UTF-8 once decoded
/// ([`uri_unescape`]): `format` is the file's format as the reader names it
/// (`GFF3`, `GTF`), `key` and `value` the pair as the reader holds it before
/// decoding, a GFF3 value as written and a GTF value with its quotes removed.
pub fn undecodable_attribute_error(format: &str, key: &str, value: &str) -> String {
    format!("{format} attribute {key} {value:?} decodes to text that is not UTF-8")
}

/// `text` with every `%` followed by two hexadecimal digits replaced by the
/// byte they spell and everything else as written: URI::Escape's
/// `uri_unescape`, the `decode_string` of ensembl-io's `GXF.pm` (lines
/// 303-306) that a GTF attribute value (`GTF.pm` line 126) and a GFF3 type
/// (`GFF3.pm` line 163) go through. The result is bytes, as Ensembl VEP's is,
/// since an escape can spell a byte that is no UTF-8 on its own; a text without
/// `%` is returned as it is.
pub fn uri_unescape(text: &str) -> Cow<'_, [u8]> {
    unescape_bytes(text.as_bytes())
}

/// [`uri_unescape`] over bytes: the column as a file holds it, before any
/// reading as text, so an escape can complete a multi-byte character whose
/// first byte stands raw in the column, as it does when Ensembl VEP decodes the
/// bytes.
fn unescape_bytes(column: &[u8]) -> Cow<'_, [u8]> {
    if !column.contains(&b'%') {
        return Cow::Borrowed(column);
    }
    let mut decoded = Vec::with_capacity(column.len());
    let mut rest = column;
    while let Some((&byte, tail)) = rest.split_first() {
        let escaped = match tail.get(..2) {
            Some(hex) if byte == b'%' && hex.iter().all(u8::is_ascii_hexdigit) => {
                std::str::from_utf8(hex)
                    .ok()
                    .and_then(|hex| u8::from_str_radix(hex, 16).ok())
            }
            _ => None,
        };
        match escaped {
            Some(value) => {
                decoded.push(value);
                rest = &tail[2..];
            }
            None => {
                decoded.push(byte);
                rest = tail;
            }
        }
    }
    Cow::Owned(decoded)
}

/// [`unescape_bytes`] read as text: a byte sequence that is no UTF-8 once
/// decoded stands as the replacement character, since the readers hold the
/// text where Ensembl VEP holds the bytes; a column that is UTF-8 and holds no
/// `%` is borrowed as it is.
fn unescape_to_text(column: &[u8]) -> Cow<'_, str> {
    match unescape_bytes(column) {
        Cow::Borrowed(bytes) => String::from_utf8_lossy(bytes),
        Cow::Owned(bytes) => Cow::Owned(String::from_utf8_lossy(&bytes).into_owned()),
    }
}

/// A GFF3 type column as `get_type` reads it (`GFF3.pm:161-164`): the column's
/// bytes with every `%` followed by two hexadecimal digits replaced by the
/// byte they spell, then read as text (`unescape_to_text`), which no
/// accepted type matches in either engine when the bytes are no UTF-8.
pub fn decode_type(column: &[u8]) -> Cow<'_, str> {
    unescape_to_text(column)
}

/// The transcription factors of a motif record's `transcription_factor`
/// attribute, `split(/,/, $value || '')` (`RegFeat.pm:436-438`) on the value
/// as written: none for an absent value or one false to Perl (empty or `0`),
/// else the comma-separated names with the empty ones at the end dropped,
/// each name then decoded ([`decode_attribute`]), so an escaped comma stays
/// inside its name and a `%30` lists the name `0`.
pub fn transcription_factors(value: Option<&[u8]>) -> Vec<Cow<'_, str>> {
    let mut names = value.unwrap_or(b"");
    if names == b"0" {
        return Vec::new();
    }
    while let Some(kept) = names.strip_suffix(b",") {
        names = kept;
    }
    if names.is_empty() {
        return Vec::new();
    }
    names
        .split(|byte| *byte == b',')
        .map(decode_attribute)
        .collect()
}

/// Whether `byte` ends a column to the index build of Ensembl VEP's release:
/// a tab or a NUL byte (htslib 1.9 `tbx.c:89`, the test of `tbx_parse1`'s
/// loop). An annotation reader ends a column at a tab alone (ensembl-io
/// `ColumnBasedParser.pm:250`).
fn ends_an_index_column(byte: &u8) -> bool {
    matches!(byte, b'\t' | 0)
}

/// The columns of `line` as the index build ends them: at a tab or at a NUL
/// byte.
pub fn index_columns(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(ends_an_index_column)
}

/// The first of a GFF line's `start` and `end` that no `tabix` indexes, as
/// the column's name and its text. `line` is the line as the file writes it,
/// whose fourth and fifth columns, each ended at a tab or at a NUL byte
/// ([`index_columns`]), the index build reads with `strtol` in base 0
/// (htslib 1.9 `tbx.c:95` and `:105`): it fails a column no number leads
/// (`:96`, `:106`) and an end that reads below zero (`:142`), and takes a
/// start that does as the first base (`:99`). `strtol` reads on past the tab
/// of an empty or blank column into the next and stops at a NUL byte, where
/// the text it is handed ends, so a column whose blanks run into a NUL byte
/// has no number. `None` for a line whose two columns the build parses,
/// whatever else they spell (`100x`, a start of `-5`, `0x1A`, an empty start
/// before a numeric end), and for a line too short to hold them.
pub fn unindexed_coordinate(line: &[u8]) -> Option<(&'static str, &[u8])> {
    let mut rest = line.splitn(4, ends_an_index_column).nth(3);
    for (name, fails_below_zero) in [("start", false), ("end", true)] {
        let text = rest?;
        let (value, after) = match text.iter().position(ends_an_index_column) {
            Some(end) => (&text[..end], Some(&text[end + 1..])),
            None => (text, None),
        };
        if !leads_a_number(text) || (fails_below_zero && reads_below_zero(text)) {
            return Some((name, value));
        }
        rest = after;
    }
    None
}

/// `text` past the blanks `strtol` skips before a number: the tab between
/// two columns is one, a NUL byte is none.
fn after_blanks(text: &[u8]) -> &[u8] {
    let blanks = text
        .iter()
        .take_while(|byte| b" \t\n\x0b\x0c\r".contains(byte))
        .count();
    &text[blanks..]
}

/// Whether a number leads `text`, a column and what follows it on its line,
/// as `strtol` reads it: any blanks, the tab that closes an empty or blank
/// column among them, one sign, then a digit.
fn leads_a_number(text: &[u8]) -> bool {
    let text = after_blanks(text);
    let unsigned = match text.split_first() {
        Some((b'+' | b'-', unsigned)) => unsigned,
        _ => text,
    };
    unsigned.first().is_some_and(u8::is_ascii_digit)
}

/// Whether the number that leads `text` reads below zero as `strtol` with
/// base 0 reads it: a minus sign, then digits, one of them not `0`, in the
/// base a `0x` or `0` prefix sets. `-9`, `-0x1A` and `-010` do; `-0`, `-08`
/// and `-0x` read zero.
fn reads_below_zero(text: &[u8]) -> bool {
    let Some(digits) = after_blanks(text).strip_prefix(b"-") else {
        return false;
    };
    let hexadecimal = digits
        .strip_prefix(b"0x")
        .or_else(|| digits.strip_prefix(b"0X"));
    let (digits, radix) = match hexadecimal {
        Some(hexadecimal) => (hexadecimal, 16),
        None if digits.starts_with(b"0") => (digits, 8),
        None => (digits, 10),
    };
    digits
        .iter()
        .take_while(|digit| char::from(**digit).is_digit(radix))
        .any(|digit| *digit != b'0')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn id_falls_through_an_empty_id_to_the_name_then_the_coordinates() {
        assert_eq!(
            feature_id(Some("ENSR1"), Some("n"), "21", 1, 2, "enhancer"),
            "ENSR1"
        );
        assert_eq!(
            feature_id(Some(""), Some("named"), "21", 1, 2, "enhancer"),
            "named"
        );
        assert_eq!(
            feature_id(None, None, "21", 1000, 1200, "enhancer"),
            "21:1000-1200:enhancer"
        );
        assert_eq!(
            feature_id(None, Some(""), "chr21", 5, 5, "EMAR"),
            "chr21:5-5:EMAR"
        );
    }

    /// `decode_attribute` reads a value as `uri_unescape` does: every `%XX` in
    /// either letter case is its byte, a non-reserved escape (`%3E`, `%30`)
    /// included, decoded in one pass (`%252C` is `%2C`); a `%` not followed by
    /// two hexadecimal digits stays written; an escape completes a character
    /// whose first byte stands raw; a byte that is no UTF-8 once decoded reads
    /// as the replacement character; a UTF-8 value without `%` is borrowed.
    #[test]
    fn attribute_escapes_decode_once_in_either_case() {
        assert!(matches!(decode_attribute(b"plain"), Cow::Borrowed("plain")));
        assert_eq!(decode_attribute(b"a%2Cb%3Bc%3Dd%26e%20f"), "a,b;c=d&e f");
        assert_eq!(decode_attribute(b"a%2cb%3bc%3dd%26e%3Ef"), "a,b;c=d&e>f");
        assert_eq!(decode_attribute(b"%30"), "0");
        assert_eq!(decode_attribute(b"%252C"), "%2C");
        assert_eq!(decode_attribute(b"a%2Gb%2"), "a%2Gb%2");
        assert_eq!(decode_attribute(b"caf%C3%A9"), "caf\u{E9}");
        assert_eq!(decode_attribute(b"caf\xC3%A9"), "caf\u{E9}");
        assert_eq!(decode_attribute(b"caf%FF"), "caf\u{FFFD}");
    }

    /// `uri_unescape` as URI::Escape runs it: `%XX` in either letter case is
    /// the byte, a `%` not followed by two hexadecimal digits stays, the
    /// scan does not revisit a decoded `%`, and a high byte comes out as the
    /// byte, not as a character.
    #[test]
    fn uri_unescape_decodes_every_two_digit_escape_to_its_byte() {
        assert!(matches!(uri_unescape("plain"), Cow::Borrowed(b"plain")));
        assert_eq!(uri_unescape("A%2CB%2cC").as_ref(), b"A,B,C");
        assert_eq!(uri_unescape("%30").as_ref(), b"0");
        assert_eq!(uri_unescape("a%2Gb%2").as_ref(), b"a%2Gb%2");
        assert_eq!(uri_unescape("%%41").as_ref(), b"%A");
        assert_eq!(uri_unescape("%252C").as_ref(), b"%2C");
        assert_eq!(uri_unescape("%C3").as_ref(), b"\xC3");
        assert_eq!(uri_unescape("%C3%A9").as_ref(), "\u{E9}".as_bytes());
    }

    /// `decode_type` decodes every escape of the column before reading it as
    /// text, so an escape completes a character whose first byte stands raw in
    /// the column, a byte that is no UTF-8 once decoded reads as the
    /// replacement character, and a UTF-8 column without `%` is borrowed.
    #[test]
    fn type_column_decodes_to_text() {
        assert!(matches!(
            decode_type(b"CTCF_binding_site"),
            Cow::Borrowed(_)
        ));
        assert_eq!(decode_type(b"CTCF%5Fbinding%5Fsite"), "CTCF_binding_site");
        assert_eq!(decode_type(b"a%2Gb%2"), "a%2Gb%2");
        assert_eq!(decode_type(b"caf%FF"), "caf\u{FFFD}");
        assert_eq!(decode_type(b"caf\xFF"), "caf\u{FFFD}");
        assert_eq!(decode_type(b"caf%CA%80"), "caf\u{280}");
        assert_eq!(decode_type(b"mix\xCA%80"), "mix\u{280}");
    }

    /// The list is split on the commas written, so an escaped comma stays
    /// inside its name once decoded, and each name is decoded after the split;
    /// a value false to Perl as written, `0`, lists nothing where `%30` lists
    /// the name `0`.
    #[test]
    fn transcription_factors_split_as_perl_does() {
        assert_eq!(
            transcription_factors(Some(b"HNF4A,NR2F1,RXRA")),
            ["HNF4A", "NR2F1", "RXRA"]
        );
        assert_eq!(transcription_factors(Some(b"A,,B,")), ["A", "", "B"]);
        assert_eq!(transcription_factors(Some(b"A%2cB,C%3E")), ["A,B", "C>"]);
        assert_eq!(transcription_factors(Some(b"%30")), ["0"]);
        assert!(transcription_factors(Some(b"0")).is_empty());
        assert!(transcription_factors(Some(b"")).is_empty());
        assert!(transcription_factors(None).is_empty());
        assert_eq!(
            unsupported_type_warning("exon", "f.gff3.gz"),
            "Ignoring 'exon' feature type from f.gff3.gz; not a supported regulatory feature type"
        );
    }

    /// The two texts both readers of a file write: the ignored bound names the
    /// attribute, its value, the record's type and identifier and the file;
    /// the collision is Ensembl VEP 116.3's text (`RegFeat.pm:489-490`).
    #[test]
    fn bound_and_identifier_warnings_name_the_record_and_the_file() {
        assert_eq!(
            ignored_bound_warning("extended_start", "x", "promoter", "ENSR_PROMOTER", "f.gff3.gz"),
            "Ignoring extended_start 'x' of promoter ENSR_PROMOTER in f.gff3.gz; not a whole number"
        );
        assert_eq!(
            duplicate_identifier_warning("ENSR_TWICE", "f.gff3.gz"),
            "Duplicate regulatory feature identifier 'ENSR_TWICE' in f.gff3.gz; identifiers must be unique or annotations will be lost"
        );
    }

    /// A line whose start is `start` and whose end is `end`, a score of `.`
    /// after them.
    fn line(start: &str, end: &str) -> Vec<u8> {
        format!("20\tx\texon\t{start}\t{end}\t.\t+\t.\tID=e").into_bytes()
    }

    /// The index build fails a start or an end no number leads and an end
    /// that reads below zero, and the first of the two it fails is named with
    /// its text. An empty column is read on into the next, so an empty start
    /// fails before an end that is no number, an empty end before the score
    /// `.` and before a score below zero; the columns are the line's as
    /// written, so one that opens with a tab holds its type where the build
    /// reads a start. Under the lines stand the two readings they rest on:
    /// no number leads a column that is empty, blank, a sign alone, two
    /// signs, a sign and a blank or other text, on its own and before the
    /// tabs of empty columns, nor a `.` behind a tab, nor a `-` whose digit
    /// stands in the next column; and a number reads below zero when a minus
    /// sign leads digits of which one is not `0`, in the base a `0x` or `0`
    /// prefix sets, behind blanks and a tab too (`tbx.c:142`).
    #[test]
    fn a_column_no_number_leads_and_an_end_below_zero_are_not_indexed() {
        for (start, end, unindexed) in [
            ("abc", "1100", ("start", "abc")),
            ("1000", "xyz", ("end", "xyz")),
            ("", "xyz", ("start", "")),
            ("1000", "", ("end", "")),
            ("1000", "-", ("end", "-")),
            ("1000", "-9", ("end", "-9")),
            ("abc", "-9", ("start", "abc")),
            ("-", "1100", ("start", "-")),
            ("+-5", "1100", ("start", "+-5")),
            ("1000", "-0x1A", ("end", "-0x1A")),
            ("1000", "-0X1a", ("end", "-0X1a")),
            ("1000", "-010", ("end", "-010")),
            ("1000", " -9", ("end", " -9")),
        ] {
            assert_eq!(
                unindexed_coordinate(&line(start, end)),
                Some((unindexed.0, unindexed.1.as_bytes())),
                "{start:?} {end:?}"
            );
        }
        assert_eq!(
            unindexed_coordinate(b"20\tx\texon\t1000\t\t-3\t+\t.\tID=e"),
            Some(("end", b"".as_slice()))
        );
        assert_eq!(
            unindexed_coordinate(b"\t20\tx\texon\t1000\t1100\t.\t+\t.\tID=e"),
            Some(("start", b"exon".as_slice()))
        );
        assert_eq!(
            unindexed_coordinate(b"20\tx\texon\t1000\txyz"),
            Some(("end", b"xyz".as_slice()))
        );
        for column in ["", " ", "abc", "x7", "-", "+", "--7", "+-7", ".7", "- 7"] {
            assert!(!leads_a_number(column.as_bytes()), "{column:?}");
            assert!(
                !leads_a_number(format!("\t{column}\t\t").as_bytes()),
                "{column:?}"
            );
        }
        assert!(!leads_a_number(b"\t.\t+\t."));
        assert!(!leads_a_number(b"-\t7"));
        for text in [
            "-9", "-1.5", "-010", "-0x1A", "-0X1a", " \t-9", "\t-3\t+", "-9x",
        ] {
            assert!(reads_below_zero(text.as_bytes()), "{text:?}");
        }
    }

    /// Every other spelling is one the index build parses, so nothing of it
    /// is named: a start of 0, an end before its start, trailing text, a start
    /// below zero, a blank, a hexadecimal or octal prefix, a `+`, zero
    /// padding, an empty or blank start before its end, an end that reads
    /// zero after its minus sign, an empty end before a score that is a
    /// number, a byte that is no UTF-8 after the digits, and a line too short
    /// to hold the columns. Under the lines stand the two readings they rest
    /// on: a number leads a column when, after any blanks and one sign, a
    /// digit follows, a reading that runs on past the tab of an empty or
    /// blank column into the next, so an empty start is led by its end's
    /// number and an empty end by a score that is one; and no number reads
    /// below zero that reads zero behind its minus sign, has none, or that
    /// no number leads.
    #[test]
    fn a_number_led_start_and_end_are_indexed_whatever_else_they_spell() {
        for (start, end) in [
            ("1000", "1100"),
            ("0", "1100"),
            ("1100", "1000"),
            ("1000x", "1100"),
            ("-5", "1100"),
            (" 1000", "1100"),
            ("\x0b1000", "1100"),
            ("0x3E8", "1100"),
            ("+1000", "1100"),
            ("001000", "001100"),
            ("1000", "1100x"),
            ("1000.5", "1100.5"),
            ("", "1100"),
            (" ", "1100"),
            ("1000", "-0"),
            ("1000", "-00"),
            ("1000", "-08"),
            ("1000", "-0x"),
            ("1000", "-0xg"),
            ("1000", "-0.5"),
            ("1000", "+9"),
        ] {
            assert_eq!(
                unindexed_coordinate(&line(start, end)),
                None,
                "{start:?} {end:?}"
            );
        }
        for written in [
            b"20\tx\texon\t1000\t\t7\t+\t.\tID=e".as_slice(),
            b"20\tx\texon\t1000\xFF\t1100\xFF\t.\t+\t.\tID=e",
            b"20\tx\texon\t1000",
            b"20\tx\texon",
            b"",
        ] {
            assert_eq!(unindexed_coordinate(written), None, "{written:?}");
        }
        for column in [
            "7", "+7", "-7", " 7", "\x0b7", "7x", "0x1A", "007", "7.5", "7 ",
        ] {
            assert!(leads_a_number(column.as_bytes()), "{column:?}");
            assert!(
                leads_a_number(format!("\t \t{column}\tx").as_bytes()),
                "{column:?}"
            );
        }
        assert!(leads_a_number(b"\t1100\t.\t+"));
        assert!(leads_a_number(b"\t7\t+\t."));
        for text in [
            "-0", "-00", "-08", "-0x", "-0xg", "-0.5", "9", "+9", "0", "", "-", "- 9", "--9", "abc",
        ] {
            assert!(!reads_below_zero(text.as_bytes()), "{text:?}");
        }
    }

    /// The index build ends a column at a tab and at a NUL byte, so a NUL
    /// byte inside a tab-ended column is one more column to it, and a line
    /// whose columns NUL bytes end has the columns of the line with tabs in
    /// their place.
    #[test]
    fn a_column_ends_at_a_tab_or_at_a_nul_byte_to_the_index_build() {
        let columns = |line: &'static [u8]| index_columns(line).collect::<Vec<_>>();
        assert_eq!(columns(b"a\0\tb\0"), [&b"a"[..], b"", b"b", b""]);
        assert_eq!(
            columns(b"21\0x\0exon\x001000\x001100"),
            columns(b"21\tx\texon\t1000\t1100")
        );
        assert_eq!(columns(b"21\tx\tex\0on\t1000").len(), 5);
        assert_eq!(columns(b""), [b""]);
    }

    /// A NUL byte before the end column moves the columns the index build
    /// reads: one in the sequence name or the source, or inside the type,
    /// puts the type or its tail where a start is read, and the line is
    /// named for it; one that leads a start or an end leaves an empty
    /// column whose reading stops at the byte, as blanks that run into one
    /// do, where a tab is read on past. A line whose NUL bytes leave a
    /// number in the fourth and fifth columns is one the build parses: its
    /// columns NUL bytes end, its type, start or end ends in one, or the
    /// byte stands among a start's or an end's digits or behind the end.
    #[test]
    fn a_nul_byte_ends_a_column_where_a_start_and_an_end_are_read() {
        for (line, unindexed) in [
            (
                &b"20\0\tx\texon\t1000\t1100\t.\t+\t.\tID=e"[..],
                ("start", &b"exon"[..]),
            ),
            (b"2\x000\tx\texon\t1000\t1100", ("start", b"exon")),
            (b"20\tx\0\texon\t1000\t1100", ("start", b"exon")),
            (b"20\tx\0y\texon\t1000\t1100", ("start", b"exon")),
            (b"20\tx\tex\0on\t1000\t1100", ("start", b"on")),
            (b"20\tx\texon\t\x001000\t1100", ("start", b"")),
            (b"20\tx\texon\t \t\x001100", ("start", b" ")),
            (b"20\tx\texon\t1000\t\x001100", ("end", b"")),
            (b"20\tx\texon\t1000\t \x005", ("end", b" ")),
            (b"20\tx\texon\t1000\0\0", ("end", b"")),
            (b"20\0x\0exon\0abc\x001100", ("start", b"abc")),
            (b"20\0x\0exon\x001000\0-9", ("end", b"-9")),
            (b"20\0x\0exon\0", ("start", b"")),
        ] {
            assert_eq!(unindexed_coordinate(line), Some(unindexed), "{line:?}");
        }
        for line in [
            &b"20\0x\0exon\x001000\x001100\0.\0+\0.\0ID=e"[..],
            b"20\0x\0exon\x001000",
            b"20\tx\texon\0\t1000\t1100\t.\t+\t.\tID=e",
            b"20\tx\texon\t1000\0\t1100\t.\t+\t.\tID=e",
            b"20\tx\texon\t1000\t1100\0\t.\t+\t.\tID=e",
            b"20\tx\texon\t10\x0000\t1100\t.\t+\t.\tID=e",
            b"20\tx\texon\t1000\t11\x0000\t.\t+\t.\tID=e",
            b"20\tx\texon\t1000\t1100\t.\0\t+\t.\tID=e\0",
            b"20\0x\0exon\0\t1100\0.",
            b"20\0\0\0\t7\texon\t1000\t1100\t.\t+\t.\tID=e",
            b"\0\0\x007",
        ] {
            assert_eq!(unindexed_coordinate(line), None, "{line:?}");
        }
    }
}
