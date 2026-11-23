// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One value percent-decoded with its line breaks left encoded, the halves of a
//! `word(score)` part, and the pairs and counts of one cell.

use super::*;
use crate::columns::place_columns;
use crate::map::read_field_map;

fn decoded(text: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    decode_into(&mut out, text, &mut Vec::new());
    out
}

#[test]
fn an_escape_of_two_hexadecimal_digits_is_its_byte_and_any_other_percent_stays() {
    for (text, expected) in [
        (&b"plain"[..], &b"plain"[..]),
        (b"", b""),
        (b"Likely%20benign%3B%3D%25%22", b"Likely benign;=%\""),
        (b"%2C%7c%26", b",|&"),
        (b"50%2520off", b"50%20off"),
        (b"%", b"%"),
        (b"%4", b"%4"),
        (b"%zz", b"%zz"),
        (b"%4g%g4", b"%4g%g4"),
        (b"%%41", b"%A"),
        (b"%%%", b"%%%"),
        (b"100%", b"100%"),
        (b"a%4", b"a%4"),
        (b"%41%", b"A%"),
        (b"caf%C3%A9", b"caf\xc3\xa9"),
        (b"esc%FFaped", b"esc\xffaped"),
        (b"byte\xff", b"byte\xff"),
        (b"%00", b"\x00"),
    ] {
        assert_eq!(
            decoded(text),
            expected,
            "{:?}",
            String::from_utf8_lossy(text)
        );
    }
}

#[test]
fn a_tab_and_every_line_breaking_character_stay_encoded() {
    for (text, expected) in [
        (&b"a%09b"[..], &b"a%09b"[..]),
        (b"%0A%0a", b"%0A%0A"),
        (b"%0D", b"%0D"),
        (b"raw\x0bbreak", b"raw%0Bbreak"),
        (b"\x0c\x1c\x1d\x1e", b"%0C%1C%1D%1E"),
        (b"\x1f\x08\x00", b"\x1f\x08\x00"),
        (b"%1F%08", b"\x1f\x08"),
        (b"one\xc2%85two", b"one%C2%85two"),
        (b"\xc2\x85", b"%C2%85"),
        (b"%C2%85", b"%C2%85"),
        (b"%c2%85", b"%C2%85"),
        (b"\xc2\x86\xc2", b"\xc2\x86\xc2"),
        (b"%85\xc2", b"\x85\xc2"),
        (b"\xe2\x80\xaa\xe2\x80", b"\xe2\x80\xaa\xe2\x80"),
        (b"\xe2\xc2\x85", b"\xe2%C2%85"),
    ] {
        assert_eq!(
            decoded(text),
            expected,
            "{:?}",
            String::from_utf8_lossy(text)
        );
    }
}

/// Each byte as `%` and two uppercase hexadecimal digits, as `urllib.parse.quote` writes a
/// character that is not safe.
fn quoted(bytes: &[u8]) -> Vec<u8> {
    bytes
        .iter()
        .flat_map(|byte| format!("%{byte:02X}").into_bytes())
        .collect()
}

#[test]
fn the_line_and_paragraph_separators_stay_encoded_over_their_bytes() {
    for separator in ["\u{2028}", "\u{2029}"] {
        let raw = separator.as_bytes();
        let escape = quoted(raw);
        // Raw, escaped in either case, and half raw and half escaped.
        assert_eq!(decoded(raw), escape, "{separator:?}");
        assert_eq!(decoded(&escape), escape, "{separator:?}");
        assert_eq!(
            decoded(&escape.to_ascii_lowercase()),
            escape,
            "{separator:?}"
        );
        let mixed = [&raw[..1], &escape[3..6], &raw[2..], b"x"].concat();
        assert_eq!(decoded(&mixed), [escape.as_slice(), b"x"].concat());
        // A cut-short sequence before it stays as its bytes.
        let after = [b"\xf0\x9f", raw].concat();
        assert_eq!(decoded(&after), [b"\xf0\x9f", escape.as_slice()].concat());
    }
}

#[test]
fn a_word_score_part_splits_at_its_last_opening_parenthesis() {
    for (part, expected) in [
        (
            &b"deleterious(0.012)"[..],
            Some((&b"deleterious"[..], &b"0.012"[..])),
        ),
        (b"a(b)(0.5)", Some((b"a(b)", b"0.5"))),
        (b"(0.012)", Some((b"", b"0.012"))),
        (b"tolerated()", Some((b"tolerated", b""))),
        (b"()", Some((b"", b""))),
        (b"a(b)c)", Some((b"a", b"b)c"))),
        (b"a)(b)", Some((b"a)", b"b"))),
        (b"deleterious", None),
        (b"", None),
        (b")", None),
        (b"(", None),
        (b")(", None),
        (b"a(b", None),
        (b"a(b) ", None),
    ] {
        assert_eq!(
            word_score(part),
            expected,
            "{:?}",
            String::from_utf8_lossy(part)
        );
    }
}

const MAP: &str = "\
fastvep_column\tpart\tvep_field\tsource\tkind\tnote\tlist_separator\treference_set
FV_A\t1\tA_ONE\tS\ttext\tn\t-\t-
FV_A\t2\tA_LIST\tS\ttext\tn\t|\t-
FV_A\t4-6\tA_RANGE\tS\ttext\tn\t-\t-
FV_B\t1.score\tB_SCORE\tS\tnumber\tn\t-\t-
FV_B\t2.word\tB_WORD\tS\ttext\tn\t-\t-
FV_B\t2.score\tB_TWO\tS\tnumber\tn\t-\t-
";

/// The pairs of one cell as `(field, value)` text, then its other counts.
type Read = (Vec<(String, String)>, Vec<(u64, Half)>, u64, u64, bool);

fn cell(column: &str, cell: &[u8], kept: Option<&[&str]>) -> Read {
    let kept: Option<Vec<Vec<u8>>> =
        kept.map(|names| names.iter().map(|name| name.as_bytes().to_vec()).collect());
    let map = read_field_map("map.tsv", MAP.as_bytes(), kept.as_deref()).unwrap();
    let header: Vec<Vec<u8>> = ["v", "Location", "Allele", "Gene", "Feature", "FV_A", "FV_B"]
        .iter()
        .map(|name| name.as_bytes().to_vec())
        .collect();
    let placed = place_columns(Some(&header), &map, false).unwrap();
    let column = placed
        .columns
        .iter()
        .find(|placed| placed.name == column)
        .unwrap();
    let mut decoded = Decoded::default();
    decode_cell(cell, column, &mut decoded, &mut Scratch::default());
    let pairs = decoded
        .pairs
        .iter()
        .map(|(slot, tail)| {
            let tail = String::from_utf8_lossy(tail);
            let (field, value) = tail
                .strip_suffix('\n')
                .and_then(|tail| tail.split_once('\t'))
                .expect("a field, a tab, a value and a line feed");
            assert_eq!(placed.fields[*slot], field);
            (field.to_string(), value.to_string())
        })
        .collect();
    (
        pairs,
        decoded.dropped,
        decoded.deselected,
        decoded.unshaped,
        decoded.several,
    )
}

fn pairs(list: &[(&str, &str)]) -> Vec<(String, String)> {
    list.iter()
        .map(|(field, value)| (field.to_string(), value.to_string()))
        .collect()
}

#[test]
fn a_cell_gives_one_pair_per_populated_part_a_rule_reads() {
    assert_eq!(
        cell("FV_A", b"G|one|a&b%26c&|skipped|x||z|extra", None),
        (
            pairs(&[("A_ONE", "one"), ("A_LIST", "a|b&c|"), ("A_RANGE", "x||z")]),
            vec![(3, Half::Whole), (7, Half::Whole)],
            0,
            0,
            false
        )
    );
    // A range that runs past the entry joins the parts that are there; an entry of the
    // dropped part alone, and parts that are all empty, write nothing.
    assert_eq!(
        cell("FV_A", b"G||||x|y", None).0,
        pairs(&[("A_RANGE", "x|y")])
    );
    assert_eq!(cell("FV_A", b"G||||", None).0, pairs(&[]));
    assert_eq!(cell("FV_A", b"G|||||", None).0, pairs(&[]));
    assert_eq!(cell("FV_A", b"G", None), (pairs(&[]), vec![], 0, 0, false));
    assert_eq!(
        cell("FV_A", b"|%0A|%7C", None).0,
        pairs(&[("A_ONE", "%0A"), ("A_LIST", "|")])
    );
}

#[test]
fn several_entries_write_each_value_once_and_count_every_drop() {
    assert_eq!(
        cell("FV_A", b"G|one||3,G|one|two|3,|two", None),
        (
            pairs(&[("A_ONE", "one"), ("A_LIST", "two"), ("A_ONE", "two")]),
            vec![(3, Half::Whole), (3, Half::Whole)],
            0,
            0,
            true
        )
    );
    assert_eq!(cell("FV_A", b",", None), (pairs(&[]), vec![], 0, 0, true));
}

#[test]
fn a_half_is_read_from_a_part_of_that_shape_and_the_other_half_is_dropped() {
    assert_eq!(
        cell(
            "FV_B",
            b"G|deleterious(0.012)|probably%5Fdamaging(0%2E998)",
            None
        ),
        (
            pairs(&[
                ("B_SCORE", "0.012"),
                ("B_WORD", "probably_damaging"),
                ("B_TWO", "0.998")
            ]),
            vec![(1, Half::Word)],
            0,
            0,
            false
        )
    );
    // A part of another shape is counted once however many rules read it; an empty half
    // writes nothing and is no drop.
    assert_eq!(
        cell("FV_B", b"G|plain|plain,G|(0.5)|w()", None),
        (
            pairs(&[("B_SCORE", "0.5"), ("B_WORD", "w")]),
            vec![],
            0,
            2,
            true
        )
    );
}

#[test]
fn a_value_fields_leaves_out_is_counted_and_not_written() {
    assert_eq!(
        cell("FV_A", b"G|one|two||x,G|one", Some(&["A_LIST"])),
        (pairs(&[("A_LIST", "two")]), vec![], 3, 0, true)
    );
    // A part that is not `word(score)` is no value of a rule the run leaves out.
    assert_eq!(
        cell("FV_B", b"G|plain|w(1)", Some(&["B_TWO"])),
        (pairs(&[("B_TWO", "1")]), vec![], 1, 1, false)
    );
}
