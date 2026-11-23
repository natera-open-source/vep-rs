// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The lines two sequences differ in, as CPython's `difflib.unified_diff(a, b, n=0)`
//! marks them with `-` and `+`: the matching blocks of `difflib.SequenceMatcher` with
//! its default `autojunk`, and every line outside them.

use std::collections::HashMap;

/// The matcher over two sequences of items, each item a number: equal items of the
/// two sequences carry equal numbers.
struct Matcher<'a> {
    a: &'a [u32],
    b: &'a [u32],
    /// For each item, the places of `b` that hold it, in rising order. An item that
    /// `autojunk` finds popular has none.
    b2j: Vec<Vec<u32>>,
}

impl<'a> Matcher<'a> {
    fn new(a: &'a [u32], b: &'a [u32], items: usize) -> Self {
        let mut b2j: Vec<Vec<u32>> = vec![Vec::new(); items];
        for (index, item) in b.iter().enumerate() {
            b2j[*item as usize].push(index as u32);
        }
        // An item that fills more than one place in a hundred of a sequence of 200 or
        // more is popular and takes no part in the search.
        let n = b.len();
        if n >= 200 {
            let limit = n / 100 + 1;
            for places in &mut b2j {
                if places.len() > limit {
                    places.clear();
                }
            }
        }
        Matcher { a, b, b2j }
    }

    /// The longest block `a[i..i + k] == b[j..j + k]` within the two ranges, the
    /// earliest in `a` and then in `b` among the longest.
    fn find_longest_match(
        &self,
        alo: usize,
        ahi: usize,
        blo: usize,
        bhi: usize,
    ) -> (usize, usize, usize) {
        let (mut besti, mut bestj, mut bestsize) = (alo, blo, 0usize);
        // The length of the longest match that ends with a[i - 1] and b[j].
        let mut j2len: HashMap<usize, usize> = HashMap::new();
        let mut newj2len: HashMap<usize, usize> = HashMap::new();
        for i in alo..ahi {
            newj2len.clear();
            for &j in &self.b2j[self.a[i] as usize] {
                let j = j as usize;
                if j < blo {
                    continue;
                }
                if j >= bhi {
                    break;
                }
                let before = if j == 0 {
                    0
                } else {
                    j2len.get(&(j - 1)).copied().unwrap_or(0)
                };
                let k = before + 1;
                newj2len.insert(j, k);
                if k > bestsize {
                    besti = i + 1 - k;
                    bestj = j + 1 - k;
                    bestsize = k;
                }
            }
            std::mem::swap(&mut j2len, &mut newj2len);
        }
        // The best match is widened over equal items on each side: popular items took
        // no part in the search.
        while besti > alo && bestj > blo && self.a[besti - 1] == self.b[bestj - 1] {
            besti -= 1;
            bestj -= 1;
            bestsize += 1;
        }
        while besti + bestsize < ahi
            && bestj + bestsize < bhi
            && self.a[besti + bestsize] == self.b[bestj + bestsize]
        {
            bestsize += 1;
        }
        (besti, bestj, bestsize)
    }

    /// The matching blocks in rising order, adjacent ones joined, without the closing
    /// block of no length.
    fn matching_blocks(&self) -> Vec<(usize, usize, usize)> {
        let mut queue = vec![(0, self.a.len(), 0, self.b.len())];
        let mut blocks = Vec::new();
        while let Some((alo, ahi, blo, bhi)) = queue.pop() {
            let (i, j, k) = self.find_longest_match(alo, ahi, blo, bhi);
            if k > 0 {
                blocks.push((i, j, k));
                if alo < i && blo < j {
                    queue.push((alo, i, blo, j));
                }
                if i + k < ahi && j + k < bhi {
                    queue.push((i + k, ahi, j + k, bhi));
                }
            }
        }
        blocks.sort_unstable();
        let mut joined: Vec<(usize, usize, usize)> = Vec::new();
        for (i, j, k) in blocks {
            match joined.last_mut() {
                Some((i1, j1, k1)) if *i1 + *k1 == i && *j1 + *k1 == j => *k1 += k,
                _ => joined.push((i, j, k)),
            }
        }
        joined
    }
}

/// The lines of `a` and of `b` outside every matching block, by their places: what
/// `unified_diff(a, b, n=0)` writes with `-` and with `+`.
pub fn changed_lines(a: &[&str], b: &[&str]) -> (Vec<usize>, Vec<usize>) {
    let mut numbers: HashMap<&str, u32> = HashMap::new();
    let mut number = |line| {
        let next = numbers.len() as u32;
        *numbers.entry(line).or_insert(next)
    };
    let a_items: Vec<u32> = a.iter().map(|line| number(line)).collect();
    let b_items: Vec<u32> = b.iter().map(|line| number(line)).collect();
    let matcher = Matcher::new(&a_items, &b_items, numbers.len());
    let (mut removed, mut added) = (Vec::new(), Vec::new());
    let (mut i, mut j) = (0, 0);
    for (ai, bj, size) in matcher
        .matching_blocks()
        .into_iter()
        .chain([(a.len(), b.len(), 0)])
    {
        removed.extend(i..ai);
        added.extend(j..bj);
        i = ai + size;
        j = bj + size;
    }
    (removed, added)
}

#[cfg(test)]
#[path = "difflib_tests.rs"]
mod tests;
