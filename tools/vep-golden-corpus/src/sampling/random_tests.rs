// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The generator held to `tests/fixtures/random/draws.tsv`, CPython's draws for the seeds
//! and ranges `select` uses and for a spread of others. The table writes a number of ten
//! digits or more with `_` between its groups of three, as `int()` reads it.

use super::Random;
use crate::sampling::reservoir_step;

const DRAWS: &str = include_str!("../../tests/fixtures/random/draws.tsv");

struct Case {
    seed: i64,
    kind: &'static str,
    argument: &'static str,
    values: Vec<u64>,
}

/// A number of the table, without the `_` between its groups of three.
fn number<T: std::str::FromStr>(text: &str) -> T {
    let digits = text.replace('_', "");
    digits.parse().ok().expect("a number of the table")
}

fn numbers(text: &str, separator: char) -> Vec<u64> {
    text.split(separator)
        .filter(|item| !item.is_empty())
        .map(number)
        .collect()
}

fn cases() -> Vec<Case> {
    let mut lines = DRAWS.lines();
    assert_eq!(lines.next(), Some("seed\tkind\targument\tvalues"));
    lines
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let [seed, kind, argument, values] = columns[..] else {
                panic!("a row of four columns: {line}");
            };
            Case {
                seed: number(seed),
                kind,
                argument,
                values: numbers(values, ','),
            }
        })
        .collect()
}

/// What the generator draws for one row of the table.
fn drawn(case: &Case) -> Vec<u64> {
    let mut rng = Random::new(case.seed);
    match case.kind {
        "words" => {
            for _ in 0..number::<usize>(case.argument) {
                rng.getrandbits(32);
            }
            (0..4).map(|_| rng.getrandbits(32)).collect()
        }
        "randrange" => {
            let n: u64 = number(case.argument);
            (0..16).map(|_| rng.randrange(n)).collect()
        }
        "shuffle" => {
            let mut items: Vec<u64> = (0..number(case.argument)).collect();
            rng.shuffle(&mut items);
            items
        }
        "reservoir" => {
            let bounds = numbers(case.argument, '/');
            let (capacity, offers) = (bounds[0] as usize, bounds[1] as usize);
            let mut bucket: Vec<u64> = Vec::new();
            for offered in 1..=offers {
                reservoir_step(&mut rng, &mut bucket, offered as u64 - 1, offered, capacity);
            }
            bucket
        }
        "mixed" => {
            let ranges = numbers(case.argument, '/');
            let mut values: Vec<u64> = ranges.iter().map(|n| rng.randrange(*n)).collect();
            let mut items: Vec<u64> = (0..9).collect();
            rng.shuffle(&mut items);
            values.extend(items);
            values.extend(ranges.iter().map(|n| rng.randrange(*n)));
            values
        }
        other => panic!("a kind of the table: {other}"),
    }
}

#[test]
fn every_draw_of_the_table_is_the_generators() {
    let cases = cases();
    assert_eq!(cases.len(), 1078);
    for case in &cases {
        assert_eq!(
            drawn(case),
            case.values,
            "seed {} {} {}",
            case.seed,
            case.kind,
            case.argument
        );
    }
}

#[test]
fn the_table_holds_every_kind_of_draw_for_every_seed() {
    let cases = cases();
    let seeds: std::collections::BTreeSet<i64> = cases.iter().map(|case| case.seed).collect();
    assert_eq!(seeds.len(), 22);
    assert!(seeds.contains(&i64::MIN) && seeds.contains(&i64::MAX) && seeds.contains(&0));
    for kind in ["words", "randrange", "shuffle", "reservoir", "mixed"] {
        let with_kind = cases.iter().filter(|case| case.kind == kind);
        let of_kind: std::collections::BTreeSet<i64> = with_kind.map(|case| case.seed).collect();
        assert_eq!(of_kind, seeds, "{kind}");
    }
}

#[test]
fn a_seed_and_its_negative_seed_one_generator() {
    let (mut plus, mut minus) = (Random::new(7), Random::new(-7));
    for _ in 0..700 {
        assert_eq!(plus.getrandbits(32), minus.getrandbits(32));
    }
}

#[test]
fn a_draw_of_more_than_32_bits_takes_two_words_low_half_first() {
    let mut words = Random::new(1);
    let (low, high) = (words.getrandbits(32), words.getrandbits(32));
    assert_eq!(Random::new(1).getrandbits(64), low | (high << 32));
    assert_eq!(Random::new(1).getrandbits(40), low | ((high >> 24) << 32));
    assert_eq!(Random::new(1).getrandbits(33), low | ((high >> 31) << 32));
    assert_eq!(Random::new(1).getrandbits(1), low >> 31);
}

#[test]
fn a_range_of_one_draws_a_bit_until_it_is_zero() {
    let mut rng = Random::new(3);
    let mut bits = Random::new(3);
    for _ in 0..50 {
        assert_eq!(rng.randrange(1), 0);
        while bits.getrandbits(1) != 0 {}
    }
    assert_eq!(rng.getrandbits(32), bits.getrandbits(32));
}
