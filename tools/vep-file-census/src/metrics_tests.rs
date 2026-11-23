// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

// int(" +0_07 "), int("-0"), int("1_"), int("1.5"), int("")
#[test]
fn seconds_are_read_as_int_reads_them() {
    for (text, expected) in [
        ("12", Some("12")),
        (" +0_07 ", Some("7")),
        ("-0", Some("0")),
        ("-12", Some("-12")),
        ("000", Some("0")),
        ("99999999999999999999999", Some("99999999999999999999999")),
        ("\u{2003}5\n", Some("5")),
        ("1_", None),
        ("_1", None),
        ("1__2", None),
        ("1.5", None),
        ("", None),
        ("+", None),
        ("0x10", None),
    ] {
        assert_eq!(py_int(text).as_deref(), expected, "{text:?}");
    }
}

// json.dumps('café "q"\\\x7f\x01\U0001d11e')
#[test]
fn a_key_is_written_as_json_dumps_writes_it_under_ensure_ascii() {
    assert_eq!(json_string("missense_variant"), "\"missense_variant\"");
    assert_eq!(
        json_string("caf\u{e9} \"q\"\\\u{7f}\u{1}\u{1d11e}\n\t"),
        "\"caf\\u00e9 \\\"q\\\"\\\\\\u007f\\u0001\\ud834\\udd1e\\n\\t\""
    );
}

#[test]
fn a_cell_is_stripped_of_what_str_strip_strips() {
    assert_eq!(py_strip(" \t\u{a0}\u{1f}x y\u{3000}\u{85}\r\n"), "x y");
    assert_eq!(py_strip("\u{200b}x"), "\u{200b}x");
    assert_eq!(py_strip(" "), "");
}

#[test]
fn the_extra_column_is_the_first_one_named_extra_after_stripping() {
    assert_eq!(extra_column("#A\tB\t Extra \tExtra"), Some(2));
    assert_eq!(extra_column("###Extra\tB"), Some(0));
    assert_eq!(extra_column("#A\tB"), None);
    assert_eq!(extra_column("#"), None);
}

fn counted(rows: &[&str], extra: Option<usize>) -> Counts {
    let mut counts = Counts::default();
    for row in rows {
        counts.row(row, extra);
    }
    counts
}

#[test]
fn a_row_counts_its_terms_its_existing_variation_and_the_keys_it_fills() {
    let row = "v\tloc\tA\tG\tF\tT\t a , ,b,a\t-\t-\t-\t-\t-\trs1\tCLIN_SIG=x;AF=;AF=0.1;HGVSc=;XAF=1;MAX_AF=2";
    let counts = counted(&[row], Some(13));
    assert_eq!(counts.rows, 1);
    assert_eq!(counts.terms, [("a".to_string(), 2), ("b".to_string(), 1)]);
    assert_eq!(counts.existing, 1);
    // CLIN_SIG, HGVSc, HGVSp, then the nine frequency fields.
    assert_eq!(counts.keys, [1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0]);
    // Under seven columns a row is counted and nothing of it is read.
    let counts = counted(&["a\tb\tc\td\te\tf"], None);
    assert_eq!(
        (counts.rows, counts.terms.len(), counts.existing),
        (1, 0, 0)
    );
    // No line named the columns: the fourteenth is the Extra column; a named column
    // past the row's end leaves the fourteenth as well.
    assert_eq!(counted(&[row], None).keys[0], 1);
    assert_eq!(counted(&[row], Some(20)).keys[0], 1);
    assert_eq!(counted(&[row], Some(12)).keys, [0; 12]);
    // A dash and blanks are no existing variation.
    let dash = "v\tloc\tA\tG\tF\tT\tx\t-\t-\t-\t-\t-\t-\t-";
    let blank = "v\tloc\tA\tG\tF\tT\tx\t-\t-\t-\t-\t-\t \t-";
    assert_eq!(counted(&[dash, blank], Some(13)).existing, 0);
}

#[test]
fn counts_merged_from_runs_keep_the_order_the_terms_were_first_met() {
    let rows = [
        "v\tl\tA\tG\tF\tT\tzeta\t-\t-\t-\t-\t-\t-\t-",
        "v\tl\tA\tG\tF\tT\talpha,beta\t-\t-\t-\t-\t-\t-\t-",
        "v\tl\tA\tG\tF\tT\tbeta,zeta\t-\t-\t-\t-\t-\t-\t-",
        "v\tl\tA\tG\tF\tT\tgamma\t-\t-\t-\t-\t-\t-\t-",
    ];
    let whole = counted(&rows, None);
    let mut merged = counted(&rows[..1], None);
    merged.merge(counted(&rows[1..3], None));
    merged.merge(counted(&rows[3..], None));
    assert_eq!(merged.terms, whole.terms);
    // The most frequent first, equal counts in the order first met.
    let report = whole.report("4");
    let zeta = report.find("\"zeta\": 2").unwrap();
    let beta = report.find("\"beta\": 2").unwrap();
    let alpha = report.find("\"alpha\": 1").unwrap();
    let gamma = report.find("\"gamma\": 1").unwrap();
    assert!(zeta < beta && beta < alpha && alpha < gamma);
    assert!(report.ends_with("  \"elapsed_seconds\": 4\n}\n"));
    assert!(Counts::default()
        .report("0")
        .contains("\"consequence_counts\": {},\n"));
}
