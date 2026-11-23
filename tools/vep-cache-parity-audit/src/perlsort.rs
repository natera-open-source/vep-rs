// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Perl's `sort` with a comparison block, `srand` and `rand`, as perl 5.34 runs them:
//! the order `sort { rand(1) <=> rand(1) } @list` leaves a list in after `srand($seed)`.
//!
//! The sort is the merge sort of `pp_sort.c`. It first cuts the list into runs by
//! comparing neighbours (`dynprep`: pairs, or longer runs where eight pairs in a row
//! compare alike), then merges runs two by two. With a comparison that answers at
//! random the outcome depends on which elements are compared in which order, so both
//! steps follow the C code comparison for comparison.

/// The number of pairs with the same sense that justifies looking for a longer run.
const PTHRESH: usize = 8;
/// The number of elements that must compare low in a row before a merge steps in
/// doubling strides.
const RTHRESH: usize = 6;

/// Perl's random number generator: the 48-bit linear congruence of `drand48`.
#[derive(Debug, Clone)]
pub struct Drand48 {
    state: u64,
}

impl Drand48 {
    /// `srand($seed)`: the seed's low 32 bits above the fixed low 16 bits `0x330e`.
    pub fn seeded(seed: u32) -> Self {
        Drand48 {
            state: 0x330e + (u64::from(seed) << 16),
        }
    }

    /// `rand(1)`: the next number in [0, 1).
    pub fn rand(&mut self) -> f64 {
        self.state = (self.state.wrapping_mul(0x5_deec_e66d).wrapping_add(0xb)) & 0xffff_ffff_ffff;
        self.state as f64 / (1u64 << 48) as f64
    }
}

/// The number `srand` seeds with for the text of its argument: the integer the text
/// opens with (an optional sign and a fraction are read and dropped), cut to 32 bits; a
/// text that is no number, a number with an exponent and an integer past 64 bits seed
/// with all bits set.
pub fn srand_seed(text: &[u8]) -> u32 {
    const NO_NUMBER: u32 = u32::MAX;
    let is_space = |byte: &u8| matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c);
    let mut rest = text;
    while rest.first().is_some_and(is_space) {
        rest = &rest[1..];
    }
    if let Some(after) = rest.strip_prefix(b"-").or_else(|| rest.strip_prefix(b"+")) {
        rest = after;
    }
    let digits = rest.iter().take_while(|byte| byte.is_ascii_digit()).count();
    let mut value: u64 = 0;
    for digit in &rest[..digits] {
        value = match value
            .checked_mul(10)
            .and_then(|value| value.checked_add(u64::from(digit - b'0')))
        {
            Some(value) => value,
            None => return NO_NUMBER,
        };
    }
    rest = &rest[digits..];
    let mut fraction = 0;
    if let Some(after) = rest.strip_prefix(b".") {
        fraction = after
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count();
        rest = &after[fraction..];
    }
    if digits == 0 && fraction == 0 {
        return NO_NUMBER;
    }
    // A number with an exponent is no integer to `srand`, whatever its value.
    if matches!(rest.first(), Some(b'e' | b'E')) {
        return NO_NUMBER;
    }
    while rest.first().is_some_and(is_space) {
        rest = &rest[1..];
    }
    if !rest.is_empty() {
        return NO_NUMBER;
    }
    value as u32
}

/// The two arrays the sort works in. A slot holds an element, or the place where the
/// run that starts at the slot ends: the sort keeps its runs' ends in the array it is
/// not reading elements from.
struct Lists {
    slots: [Vec<usize>; 2],
}

const BASE: usize = 0;
const AUX: usize = 1;

/// Cuts `base` into runs and writes the end of each run at its start in `aux`.
/// Returns the number of runs.
fn dynprep(lists: &mut Lists, nmemb: usize, cmp: &mut impl FnMut(usize, usize) -> i32) -> usize {
    let mut runs = 0;
    let last = nmemb;
    let mut b = 0;
    let mut p2 = 0;
    let greater = |lists: &Lists, cmp: &mut dyn FnMut(usize, usize) -> i32, x: usize, y: usize| {
        cmp(lists.slots[BASE][x], lists.slots[BASE][y]) > 0
    };
    let mut sense = greater(lists, cmp, b, b + 1);
    while b < last {
        // The pairs from b on that compare with the prevailing sense end at t.
        let mut p = b + 2;
        let mut t = p;
        loop {
            p += 1;
            if p >= last {
                break;
            }
            if greater(lists, cmp, t, p) != sense {
                break;
            }
            p += 1;
            t = p;
        }
        let mut q = b;
        loop {
            let mut r = b + 2 * PTHRESH;
            p = r;
            if r >= t {
                // Too short to look for a long run.
                p = t;
                r = t;
            } else {
                // Walk back from r: do the pairs join into one run down to q?
                loop {
                    if greater(lists, cmp, p - 1, p) != sense {
                        break;
                    }
                    p -= 2;
                    if p <= q {
                        break;
                    }
                }
                if p <= q {
                    // b through r is a run: push it on towards t.
                    p = r;
                    q = r;
                    loop {
                        p += 2;
                        if p < t && greater(lists, cmp, p - 1, p) == sense {
                            q = p;
                        } else {
                            break;
                        }
                    }
                    p = q + 2;
                    r = p;
                }
            }
            if q > b {
                // A run of more than two at b.
                let mut savep = p;
                q += 2;
                p = q;
                // Pick up a single element left at the end.
                if p == t && t + 1 == last && greater(lists, cmp, p - 1, p) == sense {
                    savep = last;
                    r = last;
                    p = last;
                    q = last;
                }
                let next = p2 + (p - b);
                lists.slots[AUX][p2] = next;
                p2 = next;
                runs += 1;
                if sense {
                    // A falling run is turned round.
                    loop {
                        p -= 1;
                        if b >= p {
                            break;
                        }
                        lists.slots[BASE].swap(b, p);
                        b += 1;
                    }
                }
                p = savep;
            }
            while q < p {
                // Single pairs.
                let next = p2 + 2;
                lists.slots[AUX][p2] = next;
                p2 = next;
                runs += 1;
                if sense {
                    lists.slots[BASE].swap(q, q + 1);
                }
                q += 2;
            }
            b = p;
            if b == t && t + 1 == last {
                lists.slots[AUX][p2] = p2 + 1;
                runs += 1;
                b += 1;
            }
            q = r;
            if b >= t {
                break;
            }
        }
        sense = !sense;
    }
    runs
}

/// Sorts the numbers `0..nmemb` with `cmp` as Perl's `sort` does and returns them in
/// their final order. `cmp(a, b)` answers for the elements `a` and `b` as a sort block
/// does: below zero, zero or above.
pub fn sorted(nmemb: usize, mut cmp: impl FnMut(usize, usize) -> i32) -> Vec<usize> {
    let mut lists = Lists {
        slots: [(0..nmemb).collect(), vec![0; nmemb]],
    };
    if nmemb <= 1 {
        return lists.slots[BASE].clone();
    }
    // One entry a level: the offset of the first of two runs, and how many runs are
    // still to be combined into one there.
    let mut stack: Vec<(usize, usize)> = vec![(0, 0); 64];
    let mut sp = 0;
    let mut level: usize = 0;
    let runs = dynprep(&mut lists, nmemb, &mut cmp);
    stack[0] = (0, runs);
    let mut offset = 0;
    let which = [BASE, AUX, BASE];
    'sort: loop {
        let mut runs = stack[sp].1;
        if runs == 0 {
            // Both runs of this level are built: merge them, and go up a level for as
            // long as the level above is waiting for its second run.
            let mut list1 = which[level & 1];
            let mut list2 = which[(level & 1) + 1];
            loop {
                offset = stack[sp].0;
                let mut f1 = offset;
                let p1 = offset;
                let mut tp2 = offset;
                let l1 = lists.slots[list2][offset];
                let mut f2 = l1;
                let l2 = lists.slots[list2][l1];
                offset = l2;
                while f1 < l1 && f2 < l2 {
                    // The head that goes second is q; b and t close in on its place
                    // among the elements of the other run.
                    let (q, mut b, mut t, sense) =
                        if cmp(lists.slots[list1][f1], lists.slots[list1][f2]) <= 0 {
                            (f2, f1, l1, -1)
                        } else {
                            (f1, f2, l2, 0)
                        };
                    let mut low = |lists: &Lists, p: usize| {
                        cmp(lists.slots[list1][q], lists.slots[list1][p]) <= sense
                    };
                    let mut stride = 1;
                    let mut run = 0;
                    loop {
                        let p = b + stride;
                        if p >= t {
                            let p = t - 1;
                            if p > b && low(&lists, p) {
                                t = p;
                            } else {
                                b = p;
                            }
                            break;
                        } else if low(&lists, p) {
                            t = p;
                            break;
                        }
                        b = p;
                        run += 1;
                        if run >= RTHRESH {
                            stride += stride;
                        }
                    }
                    b += 1;
                    while b < t {
                        let p = b + (t - b - 1) / 2;
                        if low(&lists, p) {
                            t = p;
                        } else {
                            b = p + 1;
                        }
                    }
                    // Every element below q goes first, then q.
                    let (from, head) = if q == f1 {
                        (&mut f2, &mut f1)
                    } else {
                        (&mut f1, &mut f2)
                    };
                    loop {
                        lists.slots[list2][tp2] = lists.slots[list1][*from];
                        tp2 += 1;
                        *from += 1;
                        if *from >= t {
                            break;
                        }
                    }
                    lists.slots[list2][tp2] = lists.slots[list1][*head];
                    tp2 += 1;
                    *head += 1;
                }
                // What is left of one run.
                let (mut from, end) = if f1 == l1 { (f2, l2) } else { (f1, l1) };
                while from < end {
                    lists.slots[list2][tp2] = lists.slots[list1][from];
                    tp2 += 1;
                    from += 1;
                }
                lists.slots[list1][p1] = tp2;
                level -= 1;
                if level == 0 {
                    break 'sort;
                }
                sp -= 1;
                std::mem::swap(&mut list1, &mut list2);
                runs = stack[sp].1;
                if runs != 0 {
                    break;
                }
            }
        }
        stack[sp].1 = 0;
        // More than two runs become two at the level below, each of about half of them;
        // the second half waits on the stack.
        while runs > 2 {
            level += 1;
            sp += 1;
            stack[sp] = (offset, runs / 2);
            runs -= runs / 2;
        }
        let in_aux = level & 1 == 1;
        if runs == 1 {
            if in_aux {
                // The run belongs in aux and stands in base: copy it.
                let start = offset;
                let end = lists.slots[AUX][offset];
                offset = end;
                for at in start..end.max(start + 1) {
                    lists.slots[AUX][at] = lists.slots[BASE][at];
                }
                lists.slots[BASE][start] = end;
            } else if level == 0 {
                break 'sort;
            }
        } else {
            level += 1;
            sp += 1;
            stack[sp] = (offset, 0);
            if !in_aux {
                // The merged runs belong in base: the merge reads them from aux.
                let start = offset;
                let first_end = lists.slots[AUX][offset];
                offset = first_end;
                let second_end = lists.slots[AUX][first_end];
                for at in start..second_end.max(start + 1) {
                    lists.slots[AUX][at] = lists.slots[BASE][at];
                }
                lists.slots[BASE][start] = first_end;
                lists.slots[BASE][first_end] = second_end;
            }
        }
    }
    lists.slots[BASE].clone()
}

/// The order `sort { rand(1) <=> rand(1) }` leaves `nmemb` elements in after
/// `srand($seed)`: the places of the original list, in their final order.
pub fn shuffled(nmemb: usize, seed: u32) -> Vec<usize> {
    let mut random = Drand48::seeded(seed);
    sorted(nmemb, |_, _| {
        let left = random.rand();
        let right = random.rand();
        match left.partial_cmp(&right) {
            Some(std::cmp::Ordering::Less) => -1,
            Some(std::cmp::Ordering::Greater) => 1,
            _ => 0,
        }
    })
}

#[cfg(test)]
#[path = "perlsort_tests.rs"]
mod tests;
