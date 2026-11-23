// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

fn rule_counts(rule: Rule, data: &[u8]) -> u64 {
    rule.count(data)
}

// printf 'a\n#b\n\n#c' | grep -c '^#'; ... | grep -vc '^#'
#[test]
fn a_line_starts_with_the_text_or_does_not_and_the_last_line_needs_no_newline() {
    let data = b"a\n#b\n\n#c";
    assert_eq!(rule_counts(Rule::Starting(b"#".to_vec()), data), 2);
    assert_eq!(rule_counts(Rule::NotStarting(b"#".to_vec()), data), 2);
    assert_eq!(
        rule_counts(Rule::NotStarting(b"##".to_vec()), b"##x\n#y\n##\n"),
        1
    );
    assert_eq!(rule_counts(Rule::NotStarting(b"#".to_vec()), b""), 0);
    assert_eq!(rule_counts(Rule::Starting(b"#".to_vec()), b"#"), 1);
    // Every line starts with the empty text.
    assert_eq!(rule_counts(Rule::Starting(Vec::new()), data), 4);
    assert_eq!(rule_counts(Rule::NotStarting(Vec::new()), data), 0);
}

// printf '#a\0b\nc\n' | grep -vc '^#'   -> 2
// printf 'x\0\0\0y\0' | grep -vc '^#'   -> 4
// printf 'y\0y\nz\n' | grep -c y        -> 2
#[test]
fn a_zero_byte_ends_a_line_as_grep_reads_a_file_that_holds_one() {
    assert_eq!(
        rule_counts(Rule::NotStarting(b"#".to_vec()), b"#a\0b\nc\n"),
        2
    );
    assert_eq!(
        rule_counts(Rule::Starting(b"#".to_vec()), b"a\0#b\n#c\n"),
        2
    );
    assert_eq!(
        rule_counts(Rule::NotStarting(b"#".to_vec()), b"x\0\0\0y\0"),
        4
    );
    assert_eq!(
        rule_counts(Rule::Containing(b"y".to_vec()), b"y\0y\nz\n"),
        2
    );
}

// printf 'a_b a_b\nc\na_b' | grep -c 'a_b'
#[test]
fn a_line_that_holds_the_text_twice_is_counted_once() {
    assert_eq!(
        rule_counts(Rule::Containing(b"a_b".to_vec()), b"a_b a_b\nc\na_b"),
        2
    );
    assert_eq!(
        rule_counts(Rule::Containing(b"a_b".to_vec()), b"a_\nb\n"),
        0
    );
    assert_eq!(rule_counts(Rule::Containing(Vec::new()), b"a\n\nb"), 3);
}

// awk -F'\t' '!/^#/ && index(";" $8, ";CSQ=") {c++} END {print c + 0}'
#[test]
fn the_eighth_column_holds_the_key_at_its_start_or_after_a_semicolon() {
    let rule = Rule::InfoField(b"CSQ".to_vec());
    let held = [
        "1\t2\t3\t4\t5\t6\t7\tCSQ=x",
        "1\t2\t3\t4\t5\t6\t7\tAF=1;CSQ=x\tGT",
        "1\t2\t3\t4\t5\t6\t7\tAF=1;CSQ=",
    ];
    let not_held = [
        "#1\t2\t3\t4\t5\t6\t7\tCSQ=x",
        "1\t2\t3\t4\t5\t6\t7\tXCSQ=x",
        "1\t2\t3\t4\t5\t6\t7\tCSQ",
        "1\t2\t3\t4\t5\t6\tCSQ=x",
        "1\t2\t3\t4\t5\t6\t7\t.\tCSQ=x",
        "",
    ];
    for row in held {
        assert_eq!(rule_counts(rule.clone(), row.as_bytes()), 1, "{row}");
    }
    for row in not_held {
        assert_eq!(rule_counts(rule.clone(), row.as_bytes()), 0, "{row}");
    }
    let all = [held.join("\n"), not_held.join("\n")].join("\n");
    assert_eq!(rule_counts(rule, all.as_bytes()), 3);
}

// awk -F'\t' '!/^#/ && index(";" $NF ";", ";IMPACT=") {c++} END {print c + 0}'
#[test]
fn the_last_column_holds_the_key_at_its_start_or_after_a_semicolon() {
    let rule = Rule::ExtraKey(b"IMPACT".to_vec());
    assert_eq!(rule_counts(rule.clone(), b"a\tb\tIMPACT=LOW\n"), 1);
    assert_eq!(rule_counts(rule.clone(), b"a\tb\tAF=1;IMPACT=\r\n"), 1);
    assert_eq!(rule_counts(rule.clone(), b"IMPACT=LOW"), 1);
    assert_eq!(rule_counts(rule.clone(), b"a\tIMPACT=LOW\tb\n"), 0);
    assert_eq!(rule_counts(rule.clone(), b"a\tXIMPACT=LOW\n"), 0);
    assert_eq!(rule_counts(rule.clone(), b"#a\tIMPACT=LOW\n"), 0);
    // A zero byte ends no line of an awk pass.
    assert_eq!(rule_counts(rule, b"a\0\tb\nc\tIMPACT=1\0\n"), 1);
}

#[test]
fn the_counts_of_a_file_are_the_sums_over_any_cut_at_its_newlines() {
    let rules = [
        Rule::Starting(b"#".to_vec()),
        Rule::NotStarting(b"#".to_vec()),
        Rule::Containing(b"=x".to_vec()),
        Rule::InfoField(b"CSQ".to_vec()),
        Rule::ExtraKey(b"CSQ".to_vec()),
    ];
    let data = b"##h\n#c\n1\t2\t3\t4\t5\t6\t7\tCSQ=x\n\n#late\n1\t2\t3\t4\t5\t6\t7\tAF=1;CSQ=x\tCSQ=y\nlast\tCSQ=x";
    let whole = counts(&rules, data);
    assert_eq!(whole, [3, 4, 3, 2, 3]);
    let newlines: Vec<usize> = memchr::memchr_iter(b'\n', data).collect();
    for cut in newlines {
        let (head, tail) = data.split_at(cut + 1);
        let summed: Vec<u64> = counts(&rules, head)
            .iter()
            .zip(counts(&rules, tail))
            .map(|(a, b)| a + b)
            .collect();
        assert_eq!(summed, whole, "cut after byte {cut}");
    }
    assert!(rules[0].is_grep() && rules[1].is_grep() && rules[2].is_grep());
    assert!(!rules[3].is_grep() && !rules[4].is_grep());
}
