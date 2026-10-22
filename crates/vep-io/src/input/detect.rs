// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Input format detection from the first data line, as Ensembl VEP's
//! `Parser.pm` `detect_format` (release/116.2, lines 368-435) does over
//! `check_format` (ensembl-variation `Utils/VEP.pm`, release/116, lines 395-450).
//!
//! `#` lines and lines with no field are skipped. The delimiter starts as a tab
//! and carries over from one line to the next (Parser.pm:387): a line without
//! the current delimiter moves it to a run of spaces when the line has two in a
//! row, else to whitespace when it has a space, so a whitespace-only line
//! before the first record can set the delimiter that record is split on. The
//! first data line is split on it and classified by the shape of its fields:
//! one field is a region, an SPDI, a CAid, an HGVS notation or an identifier,
//! in that order; several fields are a VCF record when the second is a number,
//! the fourth a reference allele and a fifth exists, else an Ensembl record
//! when the second and third are numbers and the fourth names an allele pair
//! or a variant class. A line of neither shape detects nothing.

use std::borrow::Cow;
use std::fmt;
use std::io::{self, BufRead};

use super::{read_until_bounded, Delimiter};

/// The input formats Ensembl VEP's `%FORMAT_MAP` names (Parser.pm:88-96), in
/// the lower-case spelling `--format` takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    Vcf,
    Ensembl,
    Id,
    Hgvs,
    Region,
    Spdi,
    Caid,
}

impl InputFormat {
    /// The format's name as `--format` takes it and `detect_format` returns it.
    pub fn as_str(self) -> &'static str {
        match self {
            InputFormat::Vcf => "vcf",
            InputFormat::Ensembl => "ensembl",
            InputFormat::Id => "id",
            InputFormat::Hgvs => "hgvs",
            InputFormat::Region => "region",
            InputFormat::Spdi => "spdi",
            InputFormat::Caid => "caid",
        }
    }

    /// The format a `--format` value names, case-insensitively; `None` for `guess`
    /// and for anything else.
    pub fn from_flag(value: &str) -> Option<Self> {
        match value.to_ascii_lowercase().as_str() {
            "vcf" => Some(InputFormat::Vcf),
            "ensembl" => Some(InputFormat::Ensembl),
            "id" => Some(InputFormat::Id),
            "hgvs" => Some(InputFormat::Hgvs),
            "region" => Some(InputFormat::Region),
            "spdi" => Some(InputFormat::Spdi),
            "caid" => Some(InputFormat::Caid),
            _ => None,
        }
    }
}

impl fmt::Display for InputFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// What `detect_format` read from the input: the format of the first data line,
/// if one was classified, every byte consumed up to and including that line, to
/// be replayed ahead of the rest of the input, and the delimiter the line was
/// split on, which Ensembl VEP keeps for the parser it then builds
/// (Parser.pm:414-415 store it; :159 hands it on).
pub struct Detection {
    pub format: Option<InputFormat>,
    pub consumed: Vec<u8>,
    pub delimiter: Delimiter,
}

/// Reads past `#` lines and lines with no field, then classifies the first data
/// line as [`detect_line`] does, except that a whitespace-only line before it
/// can change the delimiter it is split on; an input with no data line detects
/// nothing. The bytes read are returned with the verdict, so the caller can put
/// them back in front of the reader (`detect_format` seeks its handle back).
pub fn detect_format<R: BufRead>(reader: &mut R) -> io::Result<Detection> {
    let mut consumed = Vec::new();
    let mut line = Vec::new();
    let mut line_number = 0;
    // Kept across lines, as Perl's `$delimiter` is (Parser.pm:387).
    let mut delimiter = Delimiter::Tab;
    loop {
        line.clear();
        line_number += 1;
        if read_until_bounded(reader, b'\n', &mut line, line_number)? == 0 {
            return Ok(Detection {
                format: None,
                consumed,
                delimiter,
            });
        }
        consumed.extend_from_slice(&line);
        if line.starts_with(b"#") {
            continue;
        }
        let cleaned_bytes = without_vertical_whitespace(line.strip_suffix(b"\n").unwrap_or(&line));
        let cleaned = String::from_utf8_lossy(&cleaned_bytes);
        delimiter = delimiter.chosen_for(&cleaned);
        let fields = delimiter.split(&cleaned);
        if fields.is_empty() {
            continue;
        }
        return Ok(Detection {
            format: check_format(&fields),
            consumed,
            delimiter,
        });
    }
}

/// The format of one data line, with its line end already removed: `None` for a
/// line with no field or of no known shape.
pub fn detect_line(line: &str) -> Option<InputFormat> {
    let cleaned_bytes = without_vertical_whitespace(line.as_bytes());
    let cleaned = String::from_utf8_lossy(&cleaned_bytes);
    let fields = Delimiter::Tab.chosen_for(&cleaned).split(&cleaned);
    if fields.is_empty() {
        return None;
    }
    check_format(&fields)
}

/// `s/\r|(?>\v|\x0D\x0A).+//` (Parser.pm:398) on the line's bytes, as Perl runs
/// it on the bytes it read: at the first vertical whitespace byte (`\v` on a
/// byte string is 0x0A to 0x0D and 0x85), a `\r` is removed by itself and any
/// other one together with everything after it, when something follows.
fn without_vertical_whitespace(line: &[u8]) -> Cow<'_, [u8]> {
    let is_vertical = |b: &u8| matches!(b, 0x0A..=0x0D | 0x85);
    let Some(first_vertical) = line.iter().position(is_vertical) else {
        return Cow::Borrowed(line);
    };
    let (before, from_vertical) = line.split_at(first_vertical);
    match from_vertical.strip_prefix(b"\r") {
        Some([]) => Cow::Borrowed(before),
        Some(after) => Cow::Owned([before, after].concat()),
        None if from_vertical.len() > 1 => Cow::Borrowed(before),
        None => Cow::Borrowed(line),
    }
}

/// `check_format` (Utils/VEP.pm:395-450): the format a line's fields name, or
/// `None`. The seven branches are tested in Ensembl VEP's order.
pub fn check_format(fields: &[&str]) -> Option<InputFormat> {
    let [single] = fields else {
        return check_multi_column(fields);
    };
    if is_region(single) {
        Some(InputFormat::Region)
    } else if single.matches(':').count() >= 2 {
        // `^(.*?\:){2}([^\:]+|)$`: two or more colons.
        Some(InputFormat::Spdi)
    } else if is_caid(single) {
        Some(InputFormat::Caid)
    } else if is_hgvs(single) {
        Some(InputFormat::Hgvs)
    } else {
        Some(InputFormat::Id)
    }
}

/// The VCF and Ensembl branches, for a line not of exactly one field; fewer
/// than two fields detect nothing.
fn check_multi_column(fields: &[&str]) -> Option<InputFormat> {
    let field = |i: usize| fields.get(i).copied();
    let (Some(first), Some(second)) = (field(0), field(1)) else {
        return None;
    };
    let has_word_char = first
        .bytes()
        .any(|b| b.is_ascii_alphanumeric() || b == b'_');
    if !has_word_char || !all_digits(second) {
        return None;
    }
    // `$line[3] =~ /^[ACGTN\-\.]+$/i && exists $line[4]`.
    if let (Some(reference), Some(_alt)) = (field(3), field(4)) {
        let is_reference =
            !reference.is_empty() && reference.bytes().all(|b| is_allele_byte(b) || b == b'.');
        if is_reference {
            return Some(InputFormat::Vcf);
        }
    }
    // `$line[2] =~ /^\d+$/ && $line[3] =~ /([a-z]{2,})|([ACGTN-]+\/[ACGTN-]+)/i`.
    if let (Some(end), Some(allele)) = (field(2), field(3)) {
        if all_digits(end) && (has_two_letters(allele) || has_allele_pair(allele)) {
            return Some(InputFormat::Ensembl);
        }
    }
    None
}

fn all_digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// `[a-z]{2,}` under `/i`, unanchored.
fn has_two_letters(s: &str) -> bool {
    s.as_bytes()
        .windows(2)
        .any(|w| w[0].is_ascii_alphabetic() && w[1].is_ascii_alphabetic())
}

fn is_allele_byte(b: u8) -> bool {
    matches!(
        b.to_ascii_uppercase(),
        b'A' | b'C' | b'G' | b'T' | b'N' | b'-'
    )
}

/// `[ACGTN-]+\/[ACGTN-]+` under `/i`, unanchored: a `/` with an allele byte on
/// each side.
fn has_allele_pair(s: &str) -> bool {
    s.as_bytes()
        .windows(3)
        .any(|w| is_allele_byte(w[0]) && w[1] == b'/' && is_allele_byte(w[2]))
}

/// `_valid_region_regex` (Utils/VEP.pm:391), the expression the region parser
/// reads a line by.
fn is_region(s: &str) -> bool {
    super::region::match_region(s).is_some()
}

/// `^CA\d{1,}$` under `/i`.
fn is_caid(s: &str) -> bool {
    let Some((prefix, digits)) = s.split_at_checked(2) else {
        return false;
    };
    prefix.eq_ignore_ascii_case("CA") && all_digits(digits)
}

/// `^([^\:]+)\:.*?([cgmrp]?)\.?([\*\-0-9]+.*)$`: a non-empty reference before the
/// first colon and a `*`, `-` or digit somewhere after it; the optional groups
/// in between can match nothing.
fn is_hgvs(s: &str) -> bool {
    match s.split_once(':') {
        Some((reference, rest)) if !reference.is_empty() => rest
            .bytes()
            .any(|b| b == b'*' || b == b'-' || b.is_ascii_digit()),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Read};

    /// The format cases of Ensembl VEP's `t/Parser.t` (release/116.2, lines
    /// 240-302), one line each, with the format the test there expects: the
    /// Windows-newline case is the first line of its `t/testdata/input/windows.vcf`
    /// without the `\n`, and the `incomplete` line detects nothing. The last case
    /// is the SPDI example of `check_format` (Utils/VEP.pm:407), which Parser.t
    /// has no case for.
    #[test]
    fn check_format_seven_rules() {
        let cases = [
            (
                "21\t25587759\ttest\tG\tC,<DEL:*>\t.\t.\t.",
                Some(InputFormat::Vcf),
            ),
            (
                "21\t25587758\tsv_dup\tT\t.\t.\t.\tSVTYPE=DUP;END=25587769",
                Some(InputFormat::Vcf),
            ),
            (
                "21\t25587758\tsv_dup\tT\t<DUP>\t.\t.\t.",
                Some(InputFormat::Vcf),
            ),
            ("1\t884108\t.\tT\tTCCCTTTGGC\r", Some(InputFormat::Vcf)),
            ("rs699", Some(InputFormat::Id)),
            ("21:g.25585733C>T", Some(InputFormat::Hgvs)),
            ("ENST00000352957.8:c.991G>A", Some(InputFormat::Hgvs)),
            ("ENSP00000284967.6:p.Ala331Thr", Some(InputFormat::Hgvs)),
            (
                "21 25587759 25587759 C/A + test",
                Some(InputFormat::Ensembl),
            ),
            (
                "21   25587759 25587759     C/A  +  test",
                Some(InputFormat::Ensembl),
            ),
            (
                "21 25587759 25587759 C/A/G + test",
                Some(InputFormat::Ensembl),
            ),
            (
                "21 25587759 25587769 DUP + test",
                Some(InputFormat::Ensembl),
            ),
            ("21:25587759-25587759:1/A", Some(InputFormat::Region)),
            ("21:25587759-25587759:1/DUP", Some(InputFormat::Region)),
            ("chr21:25587759-25587759:1/A", Some(InputFormat::Region)),
            ("CHR_21:25587759-25587759:1/A", Some(InputFormat::Region)),
            ("21:25587759-25587759/A", Some(InputFormat::Region)),
            ("21:25587759-25587759:+1/A", Some(InputFormat::Region)),
            ("21:25587759-25587759:-1/A", Some(InputFormat::Region)),
            ("21:25587759-25587759:1:A", Some(InputFormat::Region)),
            ("21 25587759 25587769", None),
            ("NC_000016.10:68684738:G:A", Some(InputFormat::Spdi)),
        ];
        for (line, expected) in cases {
            assert_eq!(detect_line(line), expected, "{line:?}");
        }
    }

    /// `^CA\d{1,}$` under `/i` (Utils/VEP.pm:412-415): the ClinGen allele
    /// identifier shape, in either case; a CAid with a letter after the digits or
    /// nothing after `CA` is an identifier.
    #[test]
    fn caid_line_is_detected() {
        assert_eq!(detect_line("CA9985736"), Some(InputFormat::Caid));
        assert_eq!(detect_line("ca9985736"), Some(InputFormat::Caid));
        assert_eq!(detect_line("CA9985736X"), Some(InputFormat::Id));
        assert_eq!(detect_line("CA"), Some(InputFormat::Id));
        assert_eq!(check_format(&["CA123"]), Some(InputFormat::Caid));
    }

    /// The Ensembl branch runs after the VCF branch, so a line whose fourth field
    /// is a bare reference allele with a fifth field behind it is a VCF record
    /// whatever its third field; the two-letter class and the allele pair are
    /// matched anywhere in the fourth field, and a line of two or three fields,
    /// or one whose second field is not a number, detects nothing.
    #[test]
    fn multi_column_branches_follow_check_format() {
        assert_eq!(
            detect_line("21\t100\t100\tA\t+\tid"),
            Some(InputFormat::Vcf)
        );
        assert_eq!(detect_line("21\t100\t100\tA/G"), Some(InputFormat::Ensembl));
        assert_eq!(
            detect_line("21\t100\t100\tx-a/g-y\t1"),
            Some(InputFormat::Ensembl)
        );
        assert_eq!(detect_line("21\t100"), None);
        assert_eq!(detect_line("21\tpos\t100\tA/G\t1"), None);
        assert_eq!(detect_line("\t21\t100\t100\tA/G"), None);
        assert_eq!(detect_line("  21 100 100 A/G +"), None);
        assert_eq!(detect_line(""), None);
        assert_eq!(detect_line("   "), None);
        assert_eq!(detect_line("\t\t"), None);
    }

    /// The single-field branches in their order: a region with two colons is a
    /// region, not an SPDI; a reference with no number after its colon is an
    /// identifier; a bare colon is not an HGVS reference.
    #[test]
    fn single_field_branches_follow_check_format() {
        assert_eq!(
            detect_line("21:100-100:1:A"),
            Some(InputFormat::Region),
            "region before spdi"
        );
        assert_eq!(detect_line("21:100-100:10/A"), Some(InputFormat::Spdi));
        assert_eq!(detect_line("foo:bar"), Some(InputFormat::Id));
        assert_eq!(detect_line(":123"), Some(InputFormat::Id));
        assert_eq!(detect_line("LRG_1:g.8463-1C>T"), Some(InputFormat::Hgvs));
        assert_eq!(
            detect_line("ENST00000366667:r.803C>T"),
            Some(InputFormat::Hgvs)
        );
    }

    /// The delimiter is read from the line (Parser.pm:401-412): tab first; two
    /// spaces in a row make runs of spaces the delimiter, keeping a leading
    /// empty field; a single space splits on whitespace with leading spaces
    /// ignored and a trailing vertical tab a separator, as Perl's `\s` makes it;
    /// trailing empty fields go.
    #[test]
    fn fields_split_as_detect_format_splits() {
        fn split(line: &str) -> Vec<&str> {
            Delimiter::Tab.chosen_for(line).split(line)
        }
        assert_eq!(split("a\tb c\t"), ["a", "b c"]);
        assert_eq!(split(" a  b "), ["", "a", "b"]);
        assert_eq!(split(" a b "), ["a", "b"]);
        assert_eq!(split("21 100 id A \x0B"), ["21", "100", "id", "A"]);
        assert_eq!(split("a  b\x0B"), ["a", "b\x0B"]);
        assert_eq!(split("rs699"), ["rs699"]);
        assert!(split("").is_empty());
    }

    /// The delimiter changes only on a line that does not contain the current
    /// one, and it outlives the line (Parser.pm:387, 401-412): a line of two
    /// spaces before the first record makes runs of spaces the delimiter, so a
    /// tab-delimited record after it is one field and reads as an identifier,
    /// while a single-space line leaves whitespace splitting, which still
    /// splits the tabs.
    #[test]
    fn delimiter_outlives_the_line_as_the_loop_variable_does() {
        assert_eq!(Delimiter::Tab.chosen_for("  "), Delimiter::SpaceRun);
        assert_eq!(
            Delimiter::SpaceRun.chosen_for("21\t100\tid\tA\tG"),
            Delimiter::SpaceRun
        );
        assert_eq!(
            Delimiter::Whitespace.chosen_for("21\t100\tid\tA\tG"),
            Delimiter::Whitespace
        );
        assert_eq!(Delimiter::SpaceRun.chosen_for("a b"), Delimiter::SpaceRun);
        assert_eq!(Delimiter::Tab.chosen_for("a b"), Delimiter::Whitespace);

        let record = "21\t100\tid\tA\tG\t.\t.\t.\n";
        let two_spaces = format!("  \n{record}");
        let detection = detect_format(&mut Cursor::new(two_spaces.as_bytes())).unwrap();
        assert_eq!(detection.format, Some(InputFormat::Id));
        assert_eq!(detection.delimiter, Delimiter::SpaceRun);
        let one_space = format!(" \n{record}");
        let detection = detect_format(&mut Cursor::new(one_space.as_bytes())).unwrap();
        assert_eq!(detection.format, Some(InputFormat::Vcf));
        assert_eq!(detection.delimiter, Delimiter::Whitespace);
        assert_eq!(detect_line(record.trim_end()), Some(InputFormat::Vcf));
        let detection = detect_format(&mut Cursor::new(record.as_bytes())).unwrap();
        assert_eq!(detection.delimiter, Delimiter::Tab);
    }

    /// Parser.pm:398 strips one `\r` wherever it stands, and a vertical-tab or
    /// form-feed character together with the rest of the line when text
    /// follows it; a vertical whitespace character at the end stays.
    #[test]
    fn vertical_whitespace_is_stripped_as_the_substitution_does() {
        let stripped = |line: &str| without_vertical_whitespace(line.as_bytes()).into_owned();
        assert_eq!(stripped("a\tb\r"), b"a\tb");
        assert_eq!(stripped("a\rb\r"), b"ab\r");
        assert_eq!(stripped("a\x0Bb\tc"), b"a");
        assert_eq!(stripped("a\x0B"), b"a\x0B");
        assert_eq!(stripped("plain"), b"plain");
        // The byte 0x85 is vertical whitespace, inside a UTF-8 sequence too; the
        // multi-byte line and paragraph separators are not.
        assert_eq!(stripped("a\u{C5}b"), b"a\xC3");
        assert_eq!(stripped("a\u{2028}b"), "a\u{2028}b".as_bytes());
    }

    /// Header lines and empty lines are read past, the first data line decides,
    /// and every byte read comes back for replay; an input with no data line
    /// detects nothing and is consumed whole.
    #[test]
    fn detect_format_reads_through_the_first_data_line_and_returns_the_bytes() {
        let text = "##fileformat=VCFv4.2\n#CHROM\tPOS\n\n21\t100\tid\tA\tG\t.\t.\t.\n21\t200\tid2\tA\tG\t.\t.\t.\n";
        let mut reader = Cursor::new(text.as_bytes());
        let detection = detect_format(&mut reader).unwrap();
        assert_eq!(detection.format, Some(InputFormat::Vcf));
        assert_eq!(
            String::from_utf8(detection.consumed).unwrap(),
            "##fileformat=VCFv4.2\n#CHROM\tPOS\n\n21\t100\tid\tA\tG\t.\t.\t.\n"
        );
        let mut rest = String::new();
        reader.read_to_string(&mut rest).unwrap();
        assert_eq!(rest, "21\t200\tid2\tA\tG\t.\t.\t.\n");

        let mut header_only = Cursor::new(b"##fileformat=VCFv4.2\n#CHROM\tPOS\n".as_slice());
        let detection = detect_format(&mut header_only).unwrap();
        assert_eq!(detection.format, None);
        assert_eq!(detection.consumed, b"##fileformat=VCFv4.2\n#CHROM\tPOS\n");

        let mut no_newline = Cursor::new(b"rs699".as_slice());
        let detection = detect_format(&mut no_newline).unwrap();
        assert_eq!(detection.format, Some(InputFormat::Id));
        assert_eq!(detection.consumed, b"rs699");
    }

    #[test]
    fn flag_values_name_formats() {
        for format in [
            InputFormat::Vcf,
            InputFormat::Ensembl,
            InputFormat::Id,
            InputFormat::Hgvs,
            InputFormat::Region,
            InputFormat::Spdi,
            InputFormat::Caid,
        ] {
            assert_eq!(InputFormat::from_flag(format.as_str()), Some(format));
            assert_eq!(format.to_string(), format.as_str());
        }
        assert_eq!(InputFormat::from_flag("VCF"), Some(InputFormat::Vcf));
        assert_eq!(InputFormat::from_flag("guess"), None);
    }
}
