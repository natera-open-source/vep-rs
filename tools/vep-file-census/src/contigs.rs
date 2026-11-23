// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The records of a VCF by contig, and an equal share of every contig.
//!
//! `Tally` is this `awk` pass: the first column of every line that does not start with
//! `#`, in the order first met, with how many lines carry it.
//!
//! ```text
//! awk -F'\t' '!/^#/ {if (!($1 in c)) o[++k] = $1; c[$1]++} END {for (i = 1; i <= k; i++) print o[i] "\t" c[o[i]]}'
//! ```
//!
//! `Slice` is the pass that follows it, over the same file and that tally: every header
//! line, and of each contig the records at which `floor(i * q / c)` moves, `i` the
//! record's number on its contig from 1, `c` the contig's records and `q` its share of
//! the `n` records asked for over `k` contigs: `floor(n / k)`, one more for the first
//! `n mod k` contigs of the tally, and no more than `c`. The products and quotients are
//! those of `awk`, in double precision.
//!
//! ```text
//! awk -F'\t' -v n=N -v k=K 'FNR == NR {idx[$1] = FNR; c[$1] = $2; next} /^#/ {print; next}
//!     {if (!($1 in q)) {q[$1] = int(n / k) + (idx[$1] <= n % k); if (q[$1] > c[$1]) q[$1] = c[$1]}
//!      i = ++seen[$1]; if (int(i * q[$1] / c[$1]) > int((i - 1) * q[$1] / c[$1])) print}' TALLY -
//! ```

use memchr::memchr;
use std::collections::HashMap;

/// The first tab-separated column of `line`.
fn first_column(line: &[u8]) -> &[u8] {
    &line[..memchr(b'\t', line).unwrap_or(line.len())]
}

/// The lines of `data`, a run of whole lines, without their newlines.
fn rows(data: &[u8]) -> impl Iterator<Item = &[u8]> {
    data.split_inclusive(|byte| *byte == b'\n')
        .map(|line| line.strip_suffix(b"\n").unwrap_or(line))
}

/// Contigs in the order first met, each with a count.
#[derive(Default)]
pub struct Tally {
    index: HashMap<Vec<u8>, usize>,
    pub contigs: Vec<(Vec<u8>, u64)>,
}

impl Tally {
    fn slot(&mut self, contig: &[u8]) -> usize {
        match self.index.get(contig) {
            Some(at) => *at,
            None => {
                let at = self.contigs.len();
                self.index.insert(contig.to_vec(), at);
                self.contigs.push((contig.to_vec(), 0));
                at
            }
        }
    }

    pub fn add(&mut self, contig: &[u8], count: u64) {
        let at = self.slot(contig);
        self.contigs[at].1 += count;
    }

    /// The tally of `data`, a run of whole lines.
    pub fn of(data: &[u8]) -> Tally {
        let mut tally = Tally::default();
        // The slot of the line before, which most lines share.
        let mut last: Option<usize> = None;
        for line in rows(data) {
            if line.first() == Some(&b'#') {
                continue;
            }
            let contig = first_column(line);
            let at = match last {
                Some(at) if tally.contigs[at].0 == contig => at,
                _ => tally.slot(contig),
            };
            tally.contigs[at].1 += 1;
            last = Some(at);
        }
        tally
    }

    /// Adds `other`, the tally of what follows in the file.
    pub fn merge(&mut self, other: Tally) {
        for (contig, count) in other.contigs {
            self.add(&contig, count);
        }
    }
}

/// One contig of a slice: the count its tally line states, as written and as a number,
/// its place in the tally from 1, its share once a record names it, whether that share is
/// the stated count itself, and the records of it seen so far.
struct Share {
    stated: Vec<u8>,
    records: f64,
    index: usize,
    quota: Option<f64>,
    capped: bool,
    seen: f64,
}

/// What `awk` reads of a count: the leading decimal number of the cell, zero without one.
fn number(cell: &[u8]) -> f64 {
    let text = String::from_utf8_lossy(cell);
    let text = text.trim_start_matches([' ', '\t']);
    let digits = text
        .char_indices()
        .take_while(|(at, c)| c.is_ascii_digit() || *c == '.' || *at == 0 && matches!(c, '+' | '-'))
        .count();
    text[..digits].parse().unwrap_or(0.0)
}

/// How `awk` prints a number that holds a whole value.
fn whole(value: f64) -> String {
    format!("{value:.0}")
}

pub struct Slice {
    records: f64,
    contigs: f64,
    shares: HashMap<Vec<u8>, Share>,
    /// The contigs in the order the file's records first name them.
    order: Vec<Vec<u8>>,
}

impl Slice {
    /// A slice of `records` records over the contigs of `tally`, the text of a tally:
    /// `contig<TAB>count` lines.
    pub fn new(records: u64, tally: &[u8]) -> Slice {
        let mut shares = HashMap::new();
        let mut lines = 0;
        for (at, line) in rows(tally).enumerate() {
            lines = at + 1;
            let contig = first_column(line);
            let stated = line.get(contig.len() + 1..).map_or(&[][..], first_column);
            shares.insert(
                contig.to_vec(),
                Share {
                    stated: stated.to_vec(),
                    records: number(stated),
                    index: at + 1,
                    quota: None,
                    capped: false,
                    seen: 0.0,
                },
            );
        }
        Slice {
            records: records as f64,
            contigs: lines as f64,
            shares,
            order: Vec::new(),
        }
    }

    /// Appends to `out` the lines of `data`, a run of whole lines, that the slice keeps,
    /// each ended by a newline.
    pub fn take(&mut self, data: &[u8], out: &mut Vec<u8>) {
        for line in rows(data) {
            if line.first() == Some(&b'#') {
                out.extend_from_slice(line);
                out.push(b'\n');
                continue;
            }
            let contig = first_column(line);
            if !self.shares.contains_key(contig) {
                // A contig the tally does not hold has no record to give.
                self.shares.insert(
                    contig.to_vec(),
                    Share {
                        stated: Vec::new(),
                        records: 0.0,
                        index: 0,
                        quota: None,
                        capped: false,
                        seen: 0.0,
                    },
                );
            }
            let share = self.shares.get_mut(contig).expect("inserted above");
            let quota = match share.quota {
                Some(quota) => quota,
                None => {
                    let extra = (share.index as f64) <= self.records % self.contigs;
                    let even = (self.records / self.contigs).trunc() + f64::from(u8::from(extra));
                    share.capped = even > share.records;
                    let quota = even.min(share.records);
                    share.quota = Some(quota);
                    self.order.push(contig.to_vec());
                    quota
                }
            };
            share.seen += 1.0;
            let at = share.seen;
            if (at * quota / share.records).trunc() > ((at - 1.0) * quota / share.records).trunc() {
                out.extend_from_slice(line);
                out.push(b'\n');
            }
        }
    }

    /// One line a contig the file's records named, in that order: the label, the contig,
    /// its stated count and its share.
    pub fn table(&self, label: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        for contig in &self.order {
            let share = &self.shares[contig];
            out.extend_from_slice(label);
            out.push(b'\t');
            out.extend_from_slice(contig);
            out.push(b'\t');
            out.extend_from_slice(&share.stated);
            out.push(b'\t');
            // A share cut down to the contig's records is that cell as the tally wrote it.
            if share.capped {
                out.extend_from_slice(&share.stated);
            } else {
                out.extend_from_slice(whole(share.quota.unwrap_or(0.0)).as_bytes());
            }
            out.push(b'\n');
        }
        out
    }
}

#[cfg(test)]
#[path = "contigs_tests.rs"]
mod tests;
