// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The record groups of the two sides and which of them are kept back.
//!
//! A group is every data line of a side that shares a Location. A group is kept back when
//! no list of the comparator's names it and its lines are the same text in the same order
//! on both sides; the comparator pairs entries by a key that starts with the Location and
//! numbers a repeated key in line order, so such a group pairs each line with itself.
//!
//! The comparator also reads the order of the keys both sides share, and reports the first
//! place where the two orders part. Two more conditions therefore hold for what is kept
//! back, and a group that breaks one is handed over: the kept-back lines stand in one
//! order on both sides, and every handed-over group with lines on both sides has the same
//! number of kept-back lines before each of its lines, on either side. The first place
//! where the shared keys of the remainder part is then that place in the whole outputs,
//! moved by the kept-back lines before it.

use crate::input::{Line, Side, HEADER};
use crate::layout::{Layout, Whole};
use crate::scan::{scan, Force, Scan};
use std::collections::HashMap;

const NO_GROUP: u32 = u32::MAX;
/// Rounds of the second order condition before every line is handed over: each round
/// hands over the kept-back groups that stand inside a handed-over group's span.
const ORDER_ROUNDS: usize = 64;

pub struct Split {
    /// Why every line is handed over, when it is.
    pub whole: Option<Whole>,
    /// Per line of each side, whether it is kept back.
    pub kept_expected: Vec<bool>,
    pub kept_actual: Vec<bool>,
    pub groups: u64,
    pub groups_handed_over: u64,
}

impl Split {
    fn whole(reason: Whole, expected: &Side, actual: &Side) -> Split {
        Split {
            whole: Some(reason),
            kept_expected: vec![false; expected.lines.len()],
            kept_actual: vec![false; actual.lines.len()],
            groups: 0,
            groups_handed_over: 0,
        }
    }
}

/// A side's lines by group: `order[starts[g]..starts[g + 1]]` are group `g`'s line
/// indices in file order.
struct Grouped {
    group: Vec<u32>,
    starts: Vec<u32>,
    order: Vec<u32>,
}

fn assign(lines: &[Line], location: &[u64], ids: &mut HashMap<u64, u32>) -> Vec<u32> {
    let mut group = Vec::with_capacity(lines.len());
    let mut last: Option<(u64, u32)> = None;
    for (line, &hash) in lines.iter().zip(location) {
        if line.flags & HEADER != 0 {
            group.push(NO_GROUP);
            continue;
        }
        let id = match last {
            Some((h, id)) if h == hash => id,
            _ => {
                let next = ids.len() as u32;
                *ids.entry(hash).or_insert(next)
            }
        };
        last = Some((hash, id));
        group.push(id);
    }
    group
}

fn grouped(group: Vec<u32>, groups: usize) -> Grouped {
    let mut starts = vec![0u32; groups + 1];
    for &g in group.iter().filter(|&&g| g != NO_GROUP) {
        starts[g as usize + 1] += 1;
    }
    for g in 0..groups {
        starts[g + 1] += starts[g];
    }
    let mut next = starts.clone();
    let mut order = vec![0u32; starts[groups] as usize];
    for (i, &g) in group.iter().enumerate() {
        if g != NO_GROUP {
            order[next[g as usize] as usize] = i as u32;
            next[g as usize] += 1;
        }
    }
    Grouped {
        group,
        starts,
        order,
    }
}

impl Grouped {
    fn lines(&self, g: usize) -> &[u32] {
        &self.order[self.starts[g] as usize..self.starts[g + 1] as usize]
    }
}

/// Per group, whether it is the same text in the same order on both sides and unnamed by
/// the comparator's lists.
fn equal_groups(
    expected: &Side,
    actual: &Side,
    ge: &Grouped,
    ga: &Grouped,
    forced: &[bool],
    threads: usize,
) -> Vec<bool> {
    let groups = forced.len();
    let mut kept = vec![false; groups];
    let per_worker = groups.div_ceil(threads.max(1)).max(1);
    std::thread::scope(|scope| {
        for (n, slots) in kept.chunks_mut(per_worker).enumerate() {
            scope.spawn(move || {
                for (k, slot) in slots.iter_mut().enumerate() {
                    let g = n * per_worker + k;
                    let (le, la) = (ge.lines(g), ga.lines(g));
                    *slot = !forced[g]
                        && le.len() == la.len()
                        && le.iter().zip(la).all(|(&i, &j)| {
                            expected.text(&expected.lines[i as usize])
                                == actual.text(&actual.lines[j as usize])
                        });
                }
            });
        }
    });
    kept
}

/// The first order condition: walks the kept-back lines of both sides together and hands
/// over the two groups met wherever the sides name different groups.
fn one_order(ge: &[u32], ga: &[u32], kept: &mut [bool]) {
    let is_kept = |kept: &[bool], g: u32| g != NO_GROUP && kept[g as usize];
    let (mut i, mut j) = (0, 0);
    loop {
        while i < ge.len() && !is_kept(kept, ge[i]) {
            i += 1;
        }
        while j < ga.len() && !is_kept(kept, ga[j]) {
            j += 1;
        }
        if i == ge.len() || j == ga.len() {
            break;
        }
        if ge[i] == ga[j] {
            i += 1;
            j += 1;
        } else {
            kept[ge[i] as usize] = false;
            kept[ga[j] as usize] = false;
        }
    }
    for &g in ge[i..].iter().chain(&ga[j..]) {
        if g != NO_GROUP {
            kept[g as usize] = false;
        }
    }
}

/// One round of the second order condition: hands over every kept-back group with a line
/// between two lines of one handed-over group that both sides hold. Returns whether any
/// group was handed over.
fn close_spans(ge: &[u32], ga: &[u32], two_sided: &[bool], kept: &mut [bool]) -> bool {
    let groups = kept.len();
    let mut low = vec![u32::MAX; groups];
    let mut high = vec![0u32; groups];
    let mut kept_lines = 0u32;
    for side in [ge, ga] {
        let mut before = 0u32;
        for &g in side.iter().filter(|&&g| g != NO_GROUP) {
            let g = g as usize;
            if kept[g] {
                before += 1;
            } else if two_sided[g] {
                low[g] = low[g].min(before);
                high[g] = high[g].max(before);
            }
        }
        kept_lines = before;
    }
    // `spans[k]` is the change, at kept-back line `k`, in the number of spans open over it.
    let mut spans = vec![0i32; kept_lines as usize + 1];
    let mut any = false;
    for g in 0..groups {
        if !kept[g] && two_sided[g] && low[g] < high[g] {
            spans[low[g] as usize] += 1;
            spans[high[g] as usize] -= 1;
            any = true;
        }
    }
    if !any {
        return false;
    }
    let mut inside = vec![false; groups];
    let (mut open, mut k) = (0i32, 0usize);
    for &g in ge.iter().filter(|&&g| g != NO_GROUP) {
        if kept[g as usize] {
            open += spans[k];
            k += 1;
            if open > 0 {
                inside[g as usize] = true;
            }
        }
    }
    for g in 0..groups {
        if inside[g] {
            kept[g] = false;
        }
    }
    true
}

pub fn split(
    expected: &Side,
    actual: &Side,
    force: &Force,
    unordered: bool,
    threads: usize,
) -> Result<Split, String> {
    let (le, la) = match (Layout::of(expected), Layout::of(actual)) {
        (Ok(le), Ok(la)) => (le, la),
        (Err(reason), Ok(_)) | (Ok(_), Err(reason)) => {
            return Ok(Split::whole(reason, expected, actual))
        }
        (Err(a), Err(b)) => return Ok(Split::whole(a.min(b), expected, actual)),
    };
    if le.names != la.names {
        return Ok(Split::whole(Whole::ColumnLinesDiffer, expected, actual));
    }
    // A rewrite of the reference's column line would rename a column on one side alone.
    if force
        .texts
        .iter()
        .any(|text| memchr::memmem::find(le.line, text).is_some())
    {
        return Ok(Split::whole(Whole::ListedTextInKey, expected, actual));
    }
    let half = threads.div_ceil(2);
    let (se, sa): (Scan, Scan) = std::thread::scope(|scope| {
        let se = scope.spawn(|| scan(expected, &le, force, true, half));
        let sa = scope.spawn(|| scan(actual, &la, force, false, half));
        (
            se.join().expect("reference scan thread"),
            sa.join().expect("candidate scan thread"),
        )
    });
    let met = se.whole | sa.whole;
    if let Some(reason) = Whole::ALL.into_iter().find(|w| met & (1 << *w as u32) != 0) {
        return Ok(Split::whole(reason, expected, actual));
    }

    let mut ids: HashMap<u64, u32> = HashMap::new();
    let group_e = assign(&expected.lines, &se.location, &mut ids);
    let group_a = assign(&actual.lines, &sa.location, &mut ids);
    let groups = ids.len();
    drop(ids);
    let mut forced = vec![false; groups];
    for (group, scan) in [(&group_e, &se), (&group_a, &sa)] {
        for (&g, &f) in group.iter().zip(&scan.forced) {
            if f {
                forced[g as usize] = true;
            }
        }
    }
    drop((se, sa));
    let ge = grouped(group_e, groups);
    let ga = grouped(group_a, groups);
    let mut kept = equal_groups(expected, actual, &ge, &ga, &forced, threads);
    let mut whole = None;
    if !unordered {
        let two_sided: Vec<bool> = (0..groups)
            .map(|g| !ge.lines(g).is_empty() && !ga.lines(g).is_empty())
            .collect();
        one_order(&ge.group, &ga.group, &mut kept);
        let mut rounds = 0;
        while close_spans(&ge.group, &ga.group, &two_sided, &mut kept) {
            rounds += 1;
            if rounds == ORDER_ROUNDS {
                whole = Some(Whole::OrderUnsettled);
                kept.fill(false);
                break;
            }
        }
    }
    let line_kept = |group: &[u32]| -> Vec<bool> {
        group
            .iter()
            .map(|&g| g != NO_GROUP && kept[g as usize])
            .collect()
    };
    Ok(Split {
        whole,
        kept_expected: line_kept(&ge.group),
        kept_actual: line_kept(&ga.group),
        groups: groups as u64,
        groups_handed_over: kept.iter().filter(|k| !**k).count() as u64,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Groups by letter, one line each: which stay kept back under the order conditions
    /// when `handed` are handed over from the start.
    fn order(expected: &str, actual: &str, handed: &str) -> String {
        let id = |c: char| c as u32 - 'a' as u32;
        let ge: Vec<u32> = expected.chars().map(id).collect();
        let ga: Vec<u32> = actual.chars().map(id).collect();
        let groups = ge.iter().chain(&ga).max().map_or(0, |m| *m as usize + 1);
        let mut kept: Vec<bool> = (0..groups as u32)
            .map(|g| {
                let c = char::from(b'a' + g as u8);
                let count = |side: &[u32]| side.iter().filter(|x| **x == g).count();
                !handed.contains(c) && count(&ge) > 0 && count(&ge) == count(&ga)
            })
            .collect();
        let two_sided: Vec<bool> = (0..groups as u32)
            .map(|g| ge.contains(&g) && ga.contains(&g))
            .collect();
        one_order(&ge, &ga, &mut kept);
        while close_spans(&ge, &ga, &two_sided, &mut kept) {}
        (0..groups as u32)
            .filter(|g| kept[*g as usize])
            .map(|g| char::from(b'a' + g as u8))
            .collect()
    }

    #[test]
    fn groups_in_one_order_are_all_kept_back() {
        assert_eq!(order("abcd", "abcd", ""), "abcd");
        assert_eq!(order("abcd", "abcd", "b"), "acd");
    }

    #[test]
    fn groups_that_change_places_are_handed_over() {
        assert_eq!(order("abcd", "acbd", ""), "ad");
        // d ends the reference side and starts the candidate side: b and c stand inside its span.
        assert_eq!(order("abcd", "dabc", ""), "");
    }

    #[test]
    fn a_kept_back_group_inside_a_handed_over_group_s_span_is_handed_over() {
        // Group b stands apart on the reference side, with c between its lines.
        assert_eq!(order("abcbd", "abbcd", "b"), "ad");
        // The same on both sides: c is still between two lines of a handed-over group.
        assert_eq!(order("abcbd", "abcbd", "b"), "ad");
        assert_eq!(order("abcbd", "abcbd", ""), "abcd");
    }

    #[test]
    fn a_handed_over_group_at_two_places_hands_over_what_stands_between() {
        // b is before c on one side and after d on the other.
        assert_eq!(order("abcde", "acdbe", "b"), "ae");
        // A chain: handing c over leaves it around d on the candidate side.
        assert_eq!(order("abcbdce", "abbcdce", "b"), "ae");
    }

    #[test]
    fn a_one_sided_group_moves_nothing() {
        assert_eq!(order("axbc", "abcy", ""), "abc");
    }
}
