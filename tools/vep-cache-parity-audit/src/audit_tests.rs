// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The converter's rules against what perl answers for each path.
//!
//! `tests/fixtures/rules/cases.tsv` holds, for paths of both sides, what perl 5.34.0
//! gives: an `S` line is a Storable path, the JSON path it maps to, the path above its
//! last key and the text of the first rule that accounts for its absence; a `J` line
//! is a JSON path and the text of the rule that produces it. A line feed is written
//! `\n`, the empty text `<empty>` and no result `<none>`.

use super::{contig_of, derived_rule, dropped_rule, mapped_path, parent};

fn text(shown: &str) -> Vec<u8> {
    match shown {
        "<empty>" => Vec::new(),
        shown => shown.replace("\\n", "\n").into_bytes(),
    }
}

fn result(shown: &str) -> Option<Vec<u8>> {
    (shown != "<none>").then(|| text(shown))
}

#[test]
fn the_rules_answer_as_the_perl_patterns_do() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/rules/cases.tsv"
    );
    let cases = std::fs::read_to_string(path).expect("the cases");
    let (mut storable, mut json) = (0, 0);
    for line in cases.lines() {
        let columns: Vec<&str> = line.split('\t').collect();
        match columns.as_slice() {
            ["S", path, mapped, above, rule] => {
                let path = text(path);
                assert_eq!(mapped_path(&path), text(mapped), "mapped: {line}");
                assert_eq!(
                    parent(&path).map(<[u8]>::to_vec),
                    result(above),
                    "parent: {line}"
                );
                assert_eq!(
                    dropped_rule(&path).map(|rule| rule.as_bytes().to_vec()),
                    result(rule),
                    "rule: {line}"
                );
                storable += 1;
            }
            ["J", path, rule] => {
                assert_eq!(
                    derived_rule(&text(path)).map(|rule| rule.as_bytes().to_vec()),
                    result(rule),
                    "rule: {line}"
                );
                json += 1;
            }
            _ => panic!("a line of another shape: {line:?}"),
        }
    }
    assert_eq!((storable, json), (60, 12));
}

/// `m{/(\w+)/[^/]+\.gz$}` on each name, `unknown` where it finds nothing.
#[test]
fn the_contig_is_the_directory_of_word_characters_above_the_shard() {
    let contig = |file: &str| String::from_utf8(contig_of(file.as_bytes())).expect("ASCII");
    assert_eq!(contig("/cache/116_GRCh38/21/1-1000000.gz"), "21");
    assert_eq!(contig("cache/MT/1-16569.gz"), "MT");
    assert_eq!(contig("cache/GL000191.1/1-1000000.gz"), "unknown");
    assert_eq!(contig("cache/1-b/1-1000000.gz"), "unknown");
    assert_eq!(contig("1/1-1000000.gz"), "unknown");
    assert_eq!(contig("./1/1-1000000.gz"), "1");
    assert_eq!(contig("cache/1/.gz"), "unknown");
    assert_eq!(contig("cache/1/x.gz\n"), "1");
    assert_eq!(contig("cache/1/x.json"), "unknown");
    assert_eq!(contig("x.gz"), "unknown");
}
