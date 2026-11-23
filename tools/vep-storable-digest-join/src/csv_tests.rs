// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Each expectation is what CPython 3.14's `csv` module gives for the same text or row.

use super::{records, write_record, NoColumn, Table};

fn read(text: &str) -> Vec<Vec<String>> {
    records(text).expect("records")
}

fn rows(expected: &[&[&str]]) -> Vec<Vec<String>> {
    expected
        .iter()
        .map(|row| row.iter().map(|field| (*field).to_string()).collect())
        .collect()
}

/// `list(csv.reader(io.StringIO(text, newline="")))`.
#[test]
fn records_as_the_reader_gives_them() {
    assert_eq!(read("a,b\n"), rows(&[&["a", "b"]]));
    assert_eq!(read("a,b"), rows(&[&["a", "b"]]));
    assert_eq!(read(""), rows(&[]));
    assert_eq!(read("\n"), rows(&[&[]]));
    assert_eq!(
        read("a,,c\r\n\r\nx\n"),
        rows(&[&["a", "", "c"], &[], &["x"]])
    );
    assert_eq!(read("x\ry\r\nz"), rows(&[&["x"], &["y"], &["z"]]));
    assert_eq!(read("a\n\n\nb\n"), rows(&[&["a"], &[], &[], &["b"]]));
    assert_eq!(read(" a , b \n"), rows(&[&[" a ", " b "]]));
}

#[test]
fn quoted_fields_as_the_reader_gives_them() {
    assert_eq!(
        read("\"q,1\",\"w\"\"x\",e\"f\"g\n"),
        rows(&[&["q,1", "w\"x", "e\"f\"g"]])
    );
    assert_eq!(
        read("\"multi\nline\",2\nnext\n"),
        rows(&[&["multi\nline", "2"], &["next"]])
    );
    assert_eq!(read("\"open,\nstill"), rows(&[&["open,\nstill"]]));
    assert_eq!(
        read("a,\"b\"c,\"d\"e\"f\"\n"),
        rows(&[&["a", "bc", "de\"f\""]])
    );
    assert_eq!(read("\"a\"\n\"\"\n,\n"), rows(&[&["a"], &[""], &["", ""]]));
    assert_eq!(
        read("a,\"b\r\nc\"\r\nd\n"),
        rows(&[&["a", "b\r\nc"], &["d"]])
    );
    assert_eq!(read("\"a\" ,b\n"), rows(&[&["a ", "b"]]));
}

#[test]
fn a_field_past_the_limit_is_refused() {
    let long = "x".repeat(131_072);
    assert_eq!(read(&long), vec![vec![long.clone()]]);
    assert_eq!(
        records(&(long + "y")),
        Err("field larger than field limit (131072)".to_string())
    );
}

/// `csv.writer(out, lineterminator="\n").writerow(row)`.
#[test]
fn a_record_as_the_writer_writes_it() {
    let written = |fields: &[Option<&str>]| {
        let mut out = String::new();
        write_record(&mut out, fields);
        out
    };
    assert_eq!(written(&[Some("a"), Some("b")]), "a,b\n");
    assert_eq!(
        written(&[Some("a,b"), Some("c\"d"), Some("e\nf"), Some("g\rh")]),
        "\"a,b\",\"c\"\"d\",\"e\nf\",\"g\rh\"\n"
    );
    assert_eq!(written(&[Some("")]), "\"\"\n");
    assert_eq!(written(&[None]), "\"\"\n");
    assert_eq!(written(&[Some(""), Some("")]), ",\n");
    assert_eq!(written(&[None, Some("x")]), ",x\n");
    assert_eq!(written(&[]), "\n");
    assert_eq!(written(&[Some(" sp "), Some("it's")]), " sp ,it's\n");
    assert_eq!(written(&[Some("\t")]), "\t\n");
}

/// `csv.DictReader`: the first record names the columns, an empty record is passed over,
/// a short record lacks its last values and a long one loses what is past the names.
#[test]
fn a_table_as_the_dict_reader_gives_it() {
    let table = Table::read("h1,h2\n1\n1,2,3\n\n4,5\n").expect("a table");
    assert_eq!(table.header, ["h1", "h2"]);
    assert_eq!(
        table.rows,
        vec![
            vec![Some("1".to_string()), None],
            vec![Some("1".to_string()), Some("2".to_string())],
            vec![Some("4".to_string()), Some("5".to_string())],
        ]
    );
    assert_eq!(table.get(&table.rows[0], "h2"), Ok(None));
    assert_eq!(table.get(&table.rows[2], "h2"), Ok(Some("5")));
    assert_eq!(table.get(&table.rows[2], "h3"), Err(NoColumn));

    let empty = Table::read("").expect("a table");
    assert!(empty.header.is_empty() && empty.rows.is_empty());

    // An empty first line names no column.
    let unnamed = Table::read("\n1,2\n").expect("a table");
    assert!(unnamed.header.is_empty());
    assert_eq!(unnamed.rows, vec![Vec::<Option<String>>::new()]);

    // Of two columns with one name the last one's value is the name's.
    let twice = Table::read("a,a\n1,2\n3\n").expect("a table");
    assert_eq!(twice.get(&twice.rows[0], "a"), Ok(Some("2")));
    assert_eq!(twice.get(&twice.rows[1], "a"), Ok(None));
}
