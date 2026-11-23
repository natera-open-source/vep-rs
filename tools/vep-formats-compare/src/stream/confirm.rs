// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The text of the kept-back lines of the two sides, held against each other as the two
//! sides are read.
//!
//! Where the kept-back lines stand in one order on both sides, the n-th kept-back line of
//! the reference and the n-th of the candidate are lines of one group, and a pair that
//! differs names a group two digests called equal.

use memchr::{memchr, memchr_iter};
use std::sync::mpsc::Receiver;
use std::sync::Mutex;

/// The most buffers kept for the next blocks.
const SPARE: usize = 64;

/// Buffers of kept-back text handed back for the next block.
pub(crate) struct Spare(Mutex<Vec<Vec<u8>>>);

impl Spare {
    pub(crate) fn new() -> Spare {
        Spare(Mutex::new(Vec::new()))
    }

    /// An empty buffer.
    pub(crate) fn take(&self) -> Vec<u8> {
        let spare = self.0.lock().expect("text buffers").pop();
        spare.unwrap_or_default()
    }

    fn give(&self, mut buffer: Vec<u8>) {
        buffer.clear();
        let mut spare = self.0.lock().expect("text buffers");
        if spare.len() < SPARE {
            spare.push(buffer);
        }
    }
}

/// The kept-back lines of one block, each ended by `\n`.
pub(crate) struct KeptText {
    pub(crate) text: Vec<u8>,
}

struct Stream<'a> {
    chunks: Receiver<KeptText>,
    spare: &'a Spare,
    chunk: Vec<u8>,
    at: usize,
    /// The lines of the chunks before `chunk`.
    lines_before: u64,
}

impl<'a> Stream<'a> {
    fn new(chunks: Receiver<KeptText>, spare: &'a Spare) -> Stream<'a> {
        Stream {
            chunks,
            spare,
            chunk: Vec::new(),
            at: 0,
            lines_before: 0,
        }
    }

    /// The text still to compare, `None` at the end of the side.
    fn rest(&mut self) -> Option<&[u8]> {
        while self.at == self.chunk.len() {
            self.lines_before += memchr_iter(b'\n', &self.chunk).count() as u64;
            let next = self.chunks.recv().ok()?.text;
            self.spare.give(std::mem::replace(&mut self.chunk, next));
            self.at = 0;
        }
        Some(&self.chunk[self.at..])
    }

    /// The ordinal of the line that holds the byte `ahead` bytes on.
    fn line(&self, ahead: usize) -> u64 {
        self.lines_before + memchr_iter(b'\n', &self.chunk[..self.at + ahead]).count() as u64
    }

    /// Moves `ahead` bytes on and then past the end of the line there.
    fn past_line(&mut self, ahead: usize) {
        self.at += ahead;
        self.at +=
            memchr(b'\n', &self.chunk[self.at..]).map_or(self.chunk.len() - self.at, |i| i + 1);
    }

    fn drain(self) {
        for _ in self.chunks {}
    }
}

/// What holding the two sides' kept-back lines against each other found.
pub(crate) struct Confirmed {
    /// The ordinals, among the reference's kept-back lines, of the lines whose text is
    /// not the candidate's.
    pub(crate) differing: Vec<u64>,
    /// Whether one side had kept-back lines left when the other ended.
    pub(crate) uneven: bool,
}

pub(crate) fn confirm(
    expected: Receiver<KeptText>,
    actual: Receiver<KeptText>,
    spare: &Spare,
) -> Confirmed {
    let (mut e, mut a) = (Stream::new(expected, spare), Stream::new(actual, spare));
    let mut confirmed = Confirmed {
        differing: Vec::new(),
        uneven: false,
    };
    loop {
        let (rest_e, rest_a) = match (e.rest(), a.rest()) {
            (Some(rest_e), Some(rest_a)) => (rest_e, rest_a),
            (None, None) => break,
            _ => {
                confirmed.uneven = true;
                break;
            }
        };
        let shared = rest_e.len().min(rest_a.len());
        if rest_e[..shared] == rest_a[..shared] {
            e.at += shared;
            a.at += shared;
            continue;
        }
        let ahead = rest_e
            .iter()
            .zip(rest_a)
            .position(|(x, y)| x != y)
            .expect("the two texts differ");
        confirmed.differing.push(e.line(ahead));
        e.past_line(ahead);
        a.past_line(ahead);
    }
    e.drain();
    a.drain();
    confirmed
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc::channel;

    fn confirmed(expected: &[&str], actual: &[&str]) -> (Vec<u64>, bool) {
        let feed = |chunks: &[&str]| {
            let (sender, receiver) = channel();
            for chunk in chunks {
                sender
                    .send(KeptText {
                        text: chunk.as_bytes().to_vec(),
                    })
                    .unwrap();
            }
            receiver
        };
        let confirmed = confirm(feed(expected), feed(actual), &Spare::new());
        (confirmed.differing, confirmed.uneven)
    }

    #[test]
    fn the_same_lines_cut_into_other_chunks_are_the_same() {
        assert_eq!(
            confirmed(&["a\nbc\n", "d\n"], &["a\n", "bc\nd\n"]),
            (vec![], false)
        );
        assert_eq!(confirmed(&[], &[]), (vec![], false));
        assert_eq!(confirmed(&["", "a\n"], &["a\n", ""]), (vec![], false));
    }

    #[test]
    fn each_line_that_differs_is_named_by_its_ordinal() {
        assert_eq!(
            confirmed(&["a\nbc\n", "d\ne\n"], &["a\n", "bx\nd\n", "f\n"]),
            (vec![1, 3], false)
        );
        // A line that is the start of the other side's line, either way.
        assert_eq!(confirmed(&["ab\nc\n"], &["abc\nc\n"]), (vec![0], false));
        assert_eq!(
            confirmed(&["abc\n", "c\n"], &["ab\n", "c\n"]),
            (vec![0], false)
        );
        assert_eq!(confirmed(&["\nx\n\n"], &["\ny\nz\n"]), (vec![1, 2], false));
    }

    #[test]
    fn a_side_that_ends_first_is_reported() {
        assert_eq!(confirmed(&["a\n", "b\n"], &["a\n"]), (vec![], true));
        assert_eq!(confirmed(&[], &["a\n"]), (vec![], true));
    }
}
