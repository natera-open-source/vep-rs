// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::dump::dump_records;
use crate::testing::Scratch;

const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures");

/// `tests/fixtures/expected_records_dump.jsonl` holds what `tomllib` reads of the
/// documents of `tests/fixtures/toml/`: every value and, for every table, its keys in
/// order.
#[test]
fn the_fixture_documents_are_read_as_tomllib_reads_them() {
    let directory = format!("{FIXTURES}/toml");
    let repros = load_reproductions(None, &directory).expect("the fixture documents");
    let names: Vec<&str> = repros.iter().map(|r| r.name.as_str()).collect();
    assert_eq!(names, ["crlf_lines", "empty", "every_shape", "table_order"]);
    let expected = std::fs::read_to_string(format!("{FIXTURES}/expected_records_dump.jsonl"))
        .expect("the expected dump");
    assert_eq!(dump_records(&repros), expected);
}

#[test]
fn a_table_stands_where_the_document_first_names_it() {
    let document = b"z = 1\n[b.c]\nk = 1\n[a]\nk = 2\n[b]\nlater = true\n[b.a]\nk = 3\n";
    let value = parse_toml(document).unwrap();
    let keys = |value: &Value| -> Vec<String> {
        value.items().iter().map(|(key, _)| key.clone()).collect()
    };
    assert_eq!(keys(&value), ["z", "b", "a"]);
    assert_eq!(keys(value.section("b")), ["c", "later", "a"]);
}

#[test]
fn a_document_that_is_not_toml_is_refused_in_one_line() {
    let refusal = parse_toml(b"entry = [\n").unwrap_err();
    assert!(!refusal.contains('\n'), "{refusal:?}");
    assert!(
        refusal.ends_with(')') && refusal.contains(" (at "),
        "{refusal:?}"
    );
    assert_eq!(
        parse_toml(b"entry = \"caf\xe9\"\n").unwrap_err(),
        "'utf-8' codec can't decode byte 0xe9 in position 12: invalid continuation byte"
    );
    assert_eq!(
        parse_toml(b"a = 1\r\r\nb = 2\n").unwrap_err(),
        "a carriage return that no line feed follows (at line 1, column 6)"
    );
    assert_eq!(
        parse_toml(b"a = 1\r\nb = 2\r\n").unwrap().py_repr(),
        "{'a': 1, 'b': 2}"
    );
    assert_eq!(parse_toml(b"").unwrap().py_repr(), "{}");
}

#[test]
fn the_records_of_a_directory_are_loaded_in_name_order_or_as_named() {
    let scratch = Scratch::new("records");
    for name in ["b.toml", "a.toml", ".hidden.toml", "x.TOML", "notes.txt"] {
        scratch.write(name, b"entry = 1\n");
    }
    let names = |repros: Vec<Reproduction>| -> Vec<String> {
        repros.into_iter().map(|repro| repro.name).collect()
    };
    let directory = scratch.0.as_str();
    let all = load_reproductions(None, directory).unwrap();
    assert_eq!(all[1].path, format!("{directory}/a.toml"));
    assert_eq!(all[1].data.py_repr(), "{'entry': 1}");
    assert_eq!(names(all), [".hidden", "a", "b"]);
    let named: Vec<String> = ["b", "a", "b"].iter().map(|n| n.to_string()).collect();
    assert_eq!(
        names(load_reproductions(Some(&named), directory).unwrap()),
        ["b", "a"]
    );
    assert_eq!(
        names(load_reproductions(Some(&[]), directory).unwrap()),
        [".hidden", "a", "b"]
    );
    let missing: Vec<String> = ["c", "a", "d"].iter().map(|n| n.to_string()).collect();
    match load_reproductions(Some(&missing), directory) {
        Err(LoadError::Reproduction(message)) => {
            assert_eq!(
                message,
                format!("no reproduction for c, d under {directory}")
            );
        }
        other => panic!("{other:?}"),
    }
    let empty = scratch.mkdir("empty");
    for directory in [empty, scratch.path("absent")] {
        match load_reproductions(None, &directory) {
            Err(LoadError::Reproduction(message)) => {
                assert_eq!(message, format!("no reproduction under {directory}"));
            }
            other => panic!("{other:?}"),
        }
    }
    scratch.write("bad/bad.toml", b"entry = [\n");
    match load_reproductions(None, &scratch.path("bad")) {
        Err(LoadError::Reproduction(message)) => {
            let prefix = format!("{}: ", scratch.path("bad/bad.toml"));
            assert!(message.starts_with(&prefix), "{message}");
        }
        other => panic!("{other:?}"),
    }
    scratch.mkdir("directory/is_a_directory.toml");
    match load_reproductions(None, &scratch.path("directory")) {
        Err(LoadError::Os(message)) => assert_eq!(
            message,
            format!(
                "[Errno 21] Is a directory: '{}'",
                scratch.path("directory/is_a_directory.toml")
            )
        ),
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_record_holds_its_own_run_and_one_for_each_table_of_its_also_list() {
    let document = b"entry = 'x'\nprepare = ['p']\n[input]\ntext = 'own'\n[run]\nflags = ['--a']\n\
        [cache]\nassembly = 'GRCh38'\n[[also]]\nnotes = 'n'\n[also.run]\nflags = ['--b']\n\
        [also.expected]\nvep_rs = ['v']\n[[also]]\n[also.input]\ntext = 'second'\n";
    let record = Reproduction {
        name: "x".to_string(),
        path: "/nowhere/x.toml".to_string(),
        data: parse_toml(document).unwrap(),
    };
    let runs = record.runs();
    assert_eq!(runs.len(), 3);
    for run in &runs {
        assert_eq!(
            (run.name.as_str(), run.path.as_str()),
            ("x", "/nowhere/x.toml")
        );
        assert!(!run.data.has("also"));
        assert_eq!(run.data.get("prepare").unwrap().py_repr(), "['p']");
        assert_eq!(run.cache().py_repr(), "{'assembly': 'GRCh38'}");
    }
    // The record's own run reads the record's tables.
    assert_eq!(runs[0].input().py_repr(), "{'text': 'own'}");
    assert_eq!(runs[0].run().py_repr(), "{'flags': ['--a']}");
    assert_eq!(runs[0].expected().py_repr(), "{}");
    // A further run reads its own tables and, for one it leaves out, the record's.
    assert_eq!(runs[1].input().py_repr(), "{'text': 'own'}");
    assert_eq!(runs[1].run().py_repr(), "{'flags': ['--b']}");
    assert_eq!(runs[1].expected().py_repr(), "{'vep_rs': ['v']}");
    assert!(!runs[1].data.has("notes"));
    assert_eq!(runs[2].input().py_repr(), "{'text': 'second'}");
    assert_eq!(runs[2].run().py_repr(), "{'flags': ['--a']}");

    let plain = Reproduction {
        name: "y".to_string(),
        path: "/nowhere/y.toml".to_string(),
        data: parse_toml(b"entry = 'y'\n[run]\nflags = []\n").unwrap(),
    };
    let runs = plain.runs();
    assert_eq!(runs.len(), 1);
    assert_eq!(runs[0].data, plain.data);
}
