// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The read of the two identifier shapes: the bounds of each position a shape is read in,
//! the bytes that tell whether data can hold one, and the fabricated form. Every identifier
//! here is made up, and is written as its letters and its digits apart.

use super::*;

/// A made-up identifier: its letters, then its digits.
fn made_up(letters: &str, digits: &str) -> String {
    format!("{letters}{digits}")
}

/// The identifiers the reader finds in `text`, in order.
fn read(text: &str) -> Vec<&str> {
    let places = RESTRICTED_IDENTIFIER.find_iter(text);
    places.into_iter().map(|place| &text[place]).collect()
}

#[test]
fn the_first_shape_is_read_where_it_stands_as_a_word() {
    let first = made_up("COSV", "12");
    let cases = [
        ("{}", true),
        (" {} ", true),
        ("({}).", true),
        ("\u{e9}{}\u{e9}", true),
        // The two characters of an escaped tab or line feed stand before a word.
        ("\\t{}\\t", true),
        ("\\n{}", true),
        ("\\r{}", false),
        ("t{}", false),
        ("x{}", false),
        ("_{}", false),
        ("9{}", false),
        ("{}x", false),
        ("{}_", false),
    ];
    for (shape, found) in cases {
        let text = shape.replace("{}", &first);
        let expected: Vec<&str> = if found { vec![&first] } else { Vec::new() };
        assert_eq!(read(&text), expected, "{shape}");
    }
    // The type letter is a capital, and a digit follows it.
    for other in [
        made_up("COSv", "12"),
        made_up("COS", "12"),
        "COSV".to_string(),
    ] {
        assert_eq!(read(&other), Vec::<&str>::new(), "{other}");
    }
}

#[test]
fn the_second_shape_is_read_as_a_list_member_and_as_a_json_id() {
    let second = made_up("CM", "654321");
    let cases = [
        // After a separator, up to the next separator or the end of the field or line.
        ("&{}", true),
        (",{}", true),
        ("|{}", true),
        (";{}", true),
        (";{}\n", true),
        (";{}&", true),
        (";{},", true),
        (";{}|", true),
        (";{};", true),
        (";{}\t", true),
        (";{}\r\n", true),
        (";{}\"", true),
        (";{}\\t", true),
        (";{} ", false),
        (";{}=", false),
        (";{}.1", false),
        (";{}78", false),
        // Opening a line or a field, before a separator that is no `;`.
        ("{}&", true),
        ("x\n{},", true),
        ("={}|", true),
        ("\t{}&", true),
        ("\"{},", true),
        ("\\t{}&", true),
        ("\\n{}|", true),
        ("{};", false),
        ("{}", false),
        ("\t{}\t", false),
        (" {}&", false),
        (":{}&", false),
        ("x{}&", false),
        // A JSON `id`, with or without a space after the colon.
        ("\"id\":\"{}\"", true),
        ("\"id\": \"{}\"", true),
        ("\"id\" : \"{}\"", false),
        ("\"name\":\"{}\"", false),
        ("\"id\":\"{}", false),
    ];
    for (shape, found) in cases {
        let text = shape.replace("{}", &second);
        let expected: Vec<&str> = if found { vec![&second] } else { Vec::new() };
        assert_eq!(read(&text), expected, "{shape}");
    }
    // A letter pair of the table and a count of digits the pair takes, no more and no fewer.
    for (letters, digits, found) in [
        ("CD", "123456", true),
        ("CI", "123456", true),
        ("CN", "123456", false),
        ("CM", "12345", false),
        ("CM", "1234567", true),
        ("CM", "12345678", false),
        ("CP", "12345", true),
        ("CP", "1234567", false),
        ("BM", "1234567", true),
        ("BM", "123456", false),
        ("HM", "1234", true),
        ("HM", "12345", false),
        ("HM", "123456", true),
        ("HX", "1234", false),
        ("HX", "123456", true),
    ] {
        let text = format!("rs1,{}", made_up(letters, digits));
        assert_eq!(read(&text).len(), usize::from(found), "{letters} {digits}");
    }
}

#[test]
fn identifiers_are_given_in_the_order_of_the_text_and_never_overlap() {
    let (first, second) = (made_up("COSV", "12"), made_up("CM", "654321"));
    let zeroed = made_up("CM", "000000");
    let text = format!("{second},{first}|{zeroed}&{second}");
    assert_eq!(read(&text), [&second, &first, &zeroed, &second]);
    assert_eq!(
        restricted_identifiers(&text),
        BTreeSet::from([first.clone(), second.clone()])
    );
    assert!(is_fabricated(&zeroed) && is_fabricated("COSV") && !is_fabricated(&second));
}

#[test]
fn the_bytes_of_either_shape_tell_that_data_can_hold_an_identifier() {
    let second = made_up("CM", "654321");
    let holding = [
        made_up("COSV", "1"),
        format!("x{}x", made_up("COSV", "1")),
        made_up("CD", "123456"),
        made_up("HM", "1234"),
        made_up("BM", "1234567"),
        // The second shape is held to no bound on its left.
        format!("x{second}"),
        format!("{second}\t"),
        format!("{second}."),
        format!("{second}.x"),
    ];
    for data in holding {
        assert!(holds_identifier_bytes(data.as_bytes()), "{data}");
    }
    let without = [
        "COSV".to_string(),
        made_up("COSv", "1"),
        made_up("CM", "12345"),
        made_up("CN", "123456"),
        "C".to_string(),
        String::new(),
        // A word character or a version after the second shape.
        format!("x{second}78x"),
        format!("{second}x"),
        format!("{second}_"),
        format!("{second}.1"),
    ];
    for data in without {
        assert!(!holds_identifier_bytes(data.as_bytes()), "{data}");
    }
}

#[test]
fn the_second_shape_is_read_alone_in_the_cell_of_the_known_variants() {
    let second = made_up("CD", "123456");
    let zeroed = made_up("CD", "000000");
    let table = |cell: &str| format!("#Name\tExisting_variation\tOther\nrs1\t{cell}\t{cell}\n");
    assert_eq!(
        lone_cell_identifiers(&table(&second)),
        BTreeSet::from([second.clone()])
    );
    // Read in the cell, the identifier is replaced in the other column too.
    assert_eq!(
        replace_restricted_identifiers(&table(&second), &BTreeSet::new()),
        (Cow::Owned(table(&zeroed)), 2)
    );
    // The fabricated form, a cell with more than the identifier, a row short of the column
    // and a table without the column hold none.
    let none = [
        table(&zeroed),
        table(&format!("{second},rs2")),
        format!("#Name\tExisting_variation\n{second}\n"),
        format!("#Name\tOther\nrs1\t{second}\n"),
        format!("rs1\t{second}\n"),
    ];
    for text in none {
        assert_eq!(lone_cell_identifiers(&text), BTreeSet::new(), "{text:?}");
    }
}

#[test]
fn a_named_identifier_is_replaced_where_it_stands_as_a_word() {
    let second = made_up("CM", "654321");
    let zeroed = made_up("CM", "000000");
    let named = BTreeSet::from([second.clone()]);
    let cases = [
        ("{}", true),
        ("\t{}\t", true),
        ("\"{}\"", true),
        ("see {} here", true),
        ("\\t{}\\n", true),
        ("{}.", true),
        ("{}.x", true),
        ("{}.1", false),
        (":{}", false),
        ("x{}", false),
        ("{}x", false),
        ("{}78", false),
        ("_{}", false),
    ];
    for (shape, replaced) in cases {
        let text = shape.replace("{}", &second);
        let expected = if replaced {
            shape.replace("{}", &zeroed)
        } else {
            text.clone()
        };
        let (made, count) = replace_restricted_identifiers(&text, &named);
        assert_eq!(made, expected, "{shape}");
        assert_eq!(count, usize::from(replaced), "{shape}");
    }
    // A text that holds no identifier in another form comes back as it is, and not copied.
    let kept = format!("rs1,{zeroed}\t{second}.1");
    let (made, count) = replace_restricted_identifiers(&kept, &named);
    assert!(matches!(made, Cow::Borrowed(_)) && count == 0);
    // A first-shape identifier named to the call is not read where the reader leaves it.
    let first = made_up("COSV", "12");
    let glued = format!("x{first}");
    assert_eq!(
        fabricated_identifiers(&glued, &BTreeSet::from([first])),
        glued
    );
}
