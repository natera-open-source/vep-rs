// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

fn sum(data: &[u8]) -> (u64, String) {
    let made = LineSum::of(data, None);
    (made.lines, made.hex())
}

// The sums are those of
//   sum(int.from_bytes(hashlib.sha256(line).digest(), "big") for line in lines) % 2**256
// printed as format(total, "064x").
const A: &str = "ca978112ca1bbdcafac231b39a23dc4da786eff8147c4e72b9807785afee48bb";
const A_B: &str = "08bb6928ca5517152e4b8118ff058d82334469f89d507abd84f46634858a4958";
const A_A: &str = "952f022594377b95f58463673447b89b4f0ddff028f89ce57300ef0b5fdc9176";
const A_B_C_D: &str = "4fe4d39fb795a880a08985622da33938d0b1eff09fb7f27ba1767adf5cda3202";
const ONE_TO_NINE: &str = "6f0465314894ca0bcb995477018997988633ae3f2a76e1cd33363d48b7fe1ea3";

#[test]
fn the_lines_of_a_run_sum_the_same_in_any_order() {
    assert_eq!(sum(b"a\n"), (1, A.to_string()));
    assert_eq!(sum(b"a\nb\n"), (2, A_B.to_string()));
    assert_eq!(sum(b"b\na\n"), (2, A_B.to_string()));
    assert_eq!(sum(b"d\nb\na\nc\n"), (4, A_B_C_D.to_string()));
    assert_eq!(sum(b"c\nd\na\nb\n"), (4, A_B_C_D.to_string()));
}

#[test]
fn a_line_that_stands_twice_is_summed_twice() {
    assert_eq!(sum(b"a\na\n"), (2, A_A.to_string()));
    assert_ne!(sum(b"a\na\n"), sum(b"a\n"));
}

#[test]
fn nothing_sums_to_zero() {
    assert_eq!(sum(b""), (0, "0".repeat(64)));
    assert_eq!(sum(b"\n\r\n\r"), (0, "0".repeat(64)));
}

#[test]
fn a_line_ends_at_a_newline_a_return_or_both_and_the_last_needs_neither() {
    assert_eq!(sum(b"a\r\nb\rc\n\nd"), (4, A_B_C_D.to_string()));
    assert_eq!(sum(b"a\r\r\nb"), (2, A_B.to_string()));
    // The bytes of a line are digested as they stand: a blank is part of it.
    assert_ne!(sum(b"a \nb\n"), sum(b"a\nb\n"));
}

#[test]
fn a_line_that_starts_with_the_text_is_left_out() {
    let made = LineSum::of(b"##meta\n#names\na\n#late\nb", Some(b"#"));
    assert_eq!((made.lines, made.hex()), (2, A_B.to_string()));
    let made = LineSum::of(b"##meta\n#names\na\n", Some(b"##"));
    assert_eq!(made.lines, 2);
    // Every line starts with the empty text.
    assert_eq!(LineSum::of(b"a\nb\n", Some(b"")), LineSum::default());
}

#[test]
fn the_sum_wraps_at_two_to_the_256() {
    // The nine digests add up to more than four times 2**256.
    assert_eq!(
        sum(b"1\n2\n3\n4\n5\n6\n7\n8\n9\n"),
        (9, ONE_TO_NINE.to_string())
    );
    let mut made = LineSum {
        lines: 0,
        sum: [0, u64::MAX, u64::MAX, u64::MAX],
    };
    made.add([0, 0, 0, 1]);
    assert_eq!(made.sum, [1, 0, 0, 0]);
    let mut made = LineSum {
        lines: 0,
        sum: [u64::MAX; 4],
    };
    made.add([0, 0, 0, 2]);
    assert_eq!(made.sum, [0, 0, 0, 1]);
}

#[test]
fn runs_merge_into_the_sum_of_their_lines() {
    let mut total = LineSum::of(b"3\n1\n2\n", None);
    total.merge(LineSum::of(b"9\n8\n", None));
    total.merge(LineSum::default());
    total.merge(LineSum::of(b"4\n5\n6\n7", None));
    assert_eq!((total.lines, total.hex()), (9, ONE_TO_NINE.to_string()));
}
