// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The sort against Perl's own.
//!
//! `tests/fixtures/perlsort/cases.txt` holds one line a case, written by perl 5.34.0:
//! `<n> <seed>:` and the order `sort { rand(1) <=> rand(1) } 0 .. $n - 1` leaves after
//! `srand($seed)`, its places joined by commas. The seed is the text given to `srand`, with a space written `\s`,
//! a line feed `\n` and the empty text `<empty>`.

use super::{shuffled, sorted, srand_seed, Drand48};

#[test]
fn the_randomising_sort_leaves_perls_order() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/perlsort/cases.txt"
    );
    let text = std::fs::read_to_string(path).expect("the cases");
    let mut cases = 0;
    for line in text.lines() {
        let (head, order) = line.split_once(':').expect("a case");
        let (n, seed) = head.split_once(' ').expect("a length and a seed");
        let n: usize = n.parse().expect("a length");
        let seed = match seed {
            "<empty>" => String::new(),
            seed => seed.replace("\\s", " ").replace("\\n", "\n"),
        };
        let expected: Vec<usize> = order
            .split(',')
            .filter(|place| !place.is_empty())
            .map(|place| place.parse().expect("a place"))
            .collect();
        assert_eq!(expected.len(), n, "{line}");
        assert_eq!(
            shuffled(n, srand_seed(seed.as_bytes())),
            expected,
            "n {n}, seed {seed:?}"
        );
        cases += 1;
    }
    assert_eq!(cases, 126);
}

/// `perl -e 'srand(1); printf "%.17g\n", rand(1) for 1 .. 3'`.
#[test]
fn rand_gives_perls_numbers() {
    let mut random = Drand48::seeded(1);
    assert_eq!(random.rand(), 0.041630344771878214);
    assert_eq!(random.rand(), 0.45449244472862915);
    assert_eq!(random.rand(), 0.8348172181669149);
}

#[test]
fn a_comparison_that_orders_sorts() {
    for n in [0usize, 1, 2, 3, 15, 16, 17, 33, 100, 257] {
        // Keys in a scrambled order with repeats: equal keys keep their order.
        let keys: Vec<usize> = (0..n).map(|place| (place * 7919 + 13) % 31).collect();
        let order = sorted(n, |a, b| match keys[a].cmp(&keys[b]) {
            std::cmp::Ordering::Less => -1,
            std::cmp::Ordering::Equal => 0,
            std::cmp::Ordering::Greater => 1,
        });
        let mut expected: Vec<usize> = (0..n).collect();
        expected.sort_by_key(|place| keys[*place]);
        assert_eq!(order, expected, "{n} elements");
    }
}

/// The number `srand` takes from a text: what `grok_number` reads at its head.
#[test]
fn a_seed_is_the_integer_its_text_opens_with() {
    assert_eq!(srand_seed(b"1"), 1);
    assert_eq!(srand_seed(b"0"), 0);
    assert_eq!(srand_seed(b" 9 "), 9);
    assert_eq!(srand_seed(b"-3"), 3);
    assert_eq!(srand_seed(b"+6"), 6);
    assert_eq!(srand_seed(b"1.9"), 1);
    assert_eq!(srand_seed(b".5"), 0);
    assert_eq!(srand_seed(b"5\n"), 5);
    assert_eq!(srand_seed(b"4294967296"), 0);
    assert_eq!(srand_seed(b"4294967297"), 1);
    assert_eq!(srand_seed(b"18446744073709551615"), u32::MAX);
    for text in [
        &b"abc"[..],
        b"",
        b"12abc",
        b"0x10",
        b"1_000",
        b"18446744073709551616",
        b"-",
        b".",
        b"1e3",
        b"1e",
    ] {
        assert_eq!(
            srand_seed(text),
            u32::MAX,
            "{:?}",
            String::from_utf8_lossy(text)
        );
    }
}
