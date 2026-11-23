// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The digest of a sequence of lines.
//!
//! A sequence of lines `l[0] .. l[n-1]` has the digest `sum(h(l[i]) * X^(n-1-i))` modulo
//! the prime `2^127 - 1`, `h` a 128-bit hash of a line's text. The digest of two sequences
//! joined is the first one's times `X^len(second)` plus the second one's, so a group's
//! digest does not depend on where its lines were cut into runs or blocks.

use std::sync::OnceLock;
use twox_hash::XxHash3_128;

const PRIME: u128 = (1 << 127) - 1;
const X: u128 = 0x2545_f491_4f6c_dd1d_9e37_79b9_7f4a_7c15;

fn reduce(value: u128) -> u128 {
    let folded = (value & PRIME) + (value >> 127);
    if folded >= PRIME {
        folded - PRIME
    } else {
        folded
    }
}

/// `a + b` modulo the prime, both below it.
fn add(a: u128, b: u128) -> u128 {
    let sum = a + b;
    if sum >= PRIME {
        sum - PRIME
    } else {
        sum
    }
}

/// `a * b` modulo the prime, both below it.
fn mul(a: u128, b: u128) -> u128 {
    let (a1, a0) = ((a >> 64) as u64, a as u64);
    let (b1, b0) = ((b >> 64) as u64, b as u64);
    let low = u128::from(a0) * u128::from(b0);
    // Each cross product is below 2^127, a high half being below 2^63.
    let middle = u128::from(a1) * u128::from(b0) + u128::from(a0) * u128::from(b1);
    let high = u128::from(a1) * u128::from(b1);
    // 2^127 is 1 modulo the prime, so 2^128 is 2.
    let mut sum = reduce(low);
    sum = add(sum, reduce(u128::from(middle as u64) << 64));
    sum = add(sum, reduce((middle >> 64) << 1));
    add(sum, reduce(high << 1))
}

/// `X^(2^k)` for every `k` below 64.
fn squares() -> &'static [u128; 64] {
    static SQUARES: OnceLock<[u128; 64]> = OnceLock::new();
    SQUARES.get_or_init(|| {
        let mut squares = [reduce(X); 64];
        for k in 1..64 {
            squares[k] = mul(squares[k - 1], squares[k - 1]);
        }
        squares
    })
}

fn power(mut exponent: u64) -> u128 {
    let squares = squares();
    let mut product = 1;
    while exponent != 0 {
        product = mul(product, squares[exponent.trailing_zeros() as usize]);
        exponent &= exponent - 1;
    }
    product
}

/// The digest of no line.
pub(crate) const EMPTY: u128 = 0;

/// The hash of a line's text.
pub(crate) fn of_line(text: &[u8]) -> u128 {
    XxHash3_128::oneshot(text)
}

/// The digest of the lines of `digest` followed by a line of the hash `line`.
pub(crate) fn push(digest: u128, line: u128) -> u128 {
    add(mul(digest, squares()[0]), reduce(line))
}

/// The digest of the lines of `left` followed by the `lines` lines of `right`.
pub(crate) fn join(left: u128, right: u128, lines: u64) -> u128 {
    if left == EMPTY {
        return right;
    }
    add(mul(reduce(left), power(lines)), right)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `a * b` modulo the prime by doubling and adding.
    fn mul_slowly(a: u128, mut b: u128) -> u128 {
        let (mut product, mut addend) = (0, a);
        while b != 0 {
            if b & 1 == 1 {
                product = add(product, addend);
            }
            addend = add(addend, addend);
            b >>= 1;
        }
        product
    }

    /// A sequence of values below the prime, from a 64-bit generator.
    fn values(count: usize) -> Vec<u128> {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        (0..count)
            .map(|_| reduce(u128::from(next()) << 64 | u128::from(next())))
            .collect()
    }

    #[test]
    fn a_product_is_the_product_modulo_the_prime() {
        let mut edges = vec![
            0,
            1,
            2,
            PRIME - 1,
            PRIME - 2,
            1 << 126,
            (1 << 64) - 1,
            1 << 64,
        ];
        edges.extend(values(40));
        for &a in &edges {
            for &b in &edges {
                assert_eq!(mul(a, b), mul_slowly(a, b), "{a:#x} * {b:#x}");
            }
        }
        assert_eq!(reduce(PRIME), 0);
        assert_eq!(reduce(u128::MAX), 1);
    }

    #[test]
    fn a_power_is_the_repeated_product() {
        let mut product = 1;
        for exponent in 0..200u64 {
            assert_eq!(power(exponent), product);
            product = mul(product, reduce(X));
        }
        assert_eq!(power(1 << 40), squares()[40]);
        assert_eq!(
            power((1 << 40) + 3),
            mul(squares()[40], mul(squares()[1], squares()[0]))
        );
    }

    #[test]
    fn a_digest_is_the_same_wherever_its_lines_are_cut() {
        let lines: Vec<String> = (0..37)
            .map(|i| format!("row {i}\tvalue {}", i * 7))
            .collect();
        let digest = |lines: &[String]| {
            lines
                .iter()
                .fold(EMPTY, |digest, line| push(digest, of_line(line.as_bytes())))
        };
        let whole = digest(&lines);
        for cut in 0..=lines.len() {
            let (left, right) = lines.split_at(cut);
            assert_eq!(
                join(digest(left), digest(right), right.len() as u64),
                whole,
                "cut at {cut}"
            );
        }
        let (a, rest) = lines.split_at(5);
        let (b, c) = rest.split_at(20);
        let ab = join(digest(a), digest(b), b.len() as u64);
        assert_eq!(join(ab, digest(c), c.len() as u64), whole);
    }

    #[test]
    fn lines_in_another_order_or_another_text_have_another_digest() {
        let digest = |lines: &[&str]| {
            lines
                .iter()
                .fold(EMPTY, |digest, line| push(digest, of_line(line.as_bytes())))
        };
        let base = digest(&["a", "b", "c"]);
        assert_ne!(base, digest(&["a", "c", "b"]));
        assert_ne!(base, digest(&["a", "b", "d"]));
        assert_ne!(base, digest(&["a", "b"]));
        assert_ne!(base, digest(&["a", "b", "c", ""]));
        assert_ne!(digest(&["ab", "c"]), digest(&["a", "bc"]));
        assert_ne!(digest(&[""]), EMPTY);
    }
}
