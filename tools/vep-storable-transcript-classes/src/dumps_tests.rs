// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Each expectation is what CPython 3.14's `re` and `str` give for the same text.

use super::{key_of, load, peptide_of, unquote_integer, NoLineEnd};

/// `re.sub(r"=> '(-?\d+)'(,?)$", r"=> \1\2", line)`.
#[test]
fn an_integer_in_quotes_at_the_end_of_a_line_loses_them() {
    let sub = |line: &str| unquote_integer(line).unwrap_or_else(|| line.to_string());
    assert_eq!(sub("  'start' => '12',"), "  'start' => 12,");
    assert_eq!(sub("  'start' => '12'"), "  'start' => 12");
    assert_eq!(sub("'a' => '-5',"), "'a' => -5,");
    assert_eq!(sub("'a' => '0012'"), "'a' => 0012");
    assert_eq!(sub("=> '1'"), "=> 1");
    assert_eq!(sub("'a' => 'x => '7',"), "'a' => 'x => 7,");
    assert_eq!(sub("'k' => '12', # => '3'"), "'k' => '12', # => 3");
    // A decimal digit of any script is a digit.
    assert_eq!(sub("'a' => '\u{663}',"), "'a' => \u{663},");
    for line in [
        "'a' => '--5',",
        "'a' => '5' ,",
        "'a' => '',",
        "'a' => '5x',",
        "'a' => '5',,",
        "x=> '1', ",
        "'a' => '-',",
        "'a' =>  '5',",
        "'\u{b2}' => '\u{b2}',",
        "'a' => 5,",
        "",
    ] {
        assert_eq!(unquote_integer(line), None, "{line:?}");
    }
}

/// `re.search(r"'peptide' => '([A-Z*]*)'", line)` over the lines in turn.
#[test]
fn the_peptide_is_the_first_one_a_line_holds_whole() {
    let peptide = |lines: &[&str]| {
        let lines: Vec<String> = lines.iter().map(|line| (*line).to_string()).collect();
        peptide_of(&lines).map(str::to_string)
    };
    assert_eq!(
        peptide(&["x", "'peptide' => 'MKG*',"]).as_deref(),
        Some("MKG*")
    );
    assert_eq!(
        peptide(&["'peptide' => 'MkG', 'peptide' => 'AB'"]).as_deref(),
        Some("AB")
    );
    assert_eq!(peptide(&["'peptide' => ''"]).as_deref(), Some(""));
    assert_eq!(peptide(&["'peptide' => 'ABC"]), None);
    assert_eq!(
        peptide(&["'peptide' => undef", "  'peptide' => 'ZZ'"]).as_deref(),
        Some("ZZ")
    );
    assert_eq!(
        peptide(&["'peptide' => 'A' 'peptide' => 'B'"]).as_deref(),
        Some("A")
    );
    assert_eq!(peptide(&[]), None);
}

/// `re.search(r"'([A-Za-z_0-9]+)' =>", line)`.
#[test]
fn the_key_is_the_first_word_in_quotes_before_an_arrow() {
    assert_eq!(key_of("-  'start' => 5,"), Some("start"));
    assert_eq!(key_of("+'a b' => 1, 'c' => 2"), Some("c"));
    assert_eq!(key_of("-''x' =>"), Some("x"));
    assert_eq!(key_of("+'\u{e9}' => 1 'k_9' =>"), Some("k_9"));
    assert_eq!(key_of("-'a'=> 1"), None);
    assert_eq!(key_of("+'a' =>"), Some("a"));
    assert_eq!(key_of("-'' => 'b' => "), Some("b"));
    assert_eq!(key_of("+ 'it''s' =>"), Some("s"));
    assert_eq!(key_of("-'A'b' =>"), Some("b"));
    assert_eq!(key_of("+ no key here"), None);
}

#[test]
fn a_dump_is_its_blocks_by_id() {
    let dump = "before the first block\n### ENST1 \n{\n  'dbID' => '5',\n  'start' => '12',\n  'name' => 'x'\n}\n\n### ENST2\r\nline\r### ENST1\n{}\n";
    let blocks = load(dump).expect("blocks");
    assert_eq!(blocks.len(), 2);
    // A later block of an id replaces an earlier one; `\r\n` and `\r` end a line.
    assert_eq!(blocks["ENST1"], ["{}", ""]);
    assert_eq!(blocks["ENST2"], ["line", ""]);
    let first = load("### ENST1 \n{\n  'dbID' => '5',\n  '_gene_version' => 7,\n  'start' => '12',\n  'name' => 'x'\n}\n\n").expect("blocks");
    assert_eq!(
        first["ENST1"],
        ["{", "  'start' => 12,", "  'name' => 'x'", "}", "", ""]
    );
    // The id loses the white space `str.strip()` takes off, and no other character.
    let spaced = load("###  \u{1c}\u{85}\u{3000}id \u{b} \nbody").expect("blocks");
    assert_eq!(spaced["id"], ["body"]);
    let kept = load("### \u{feff}id\u{200b}\n").expect("blocks");
    assert!(kept.contains_key("\u{feff}id\u{200b}"));
    assert!(load("no marker at all").expect("blocks").is_empty());
    assert_eq!(load("### an id without a line end"), Err(NoLineEnd));
}
