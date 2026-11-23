// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What one input becomes, line by line.
//!
//! A line ends at `\n`, `\r\n` or a lone `\r`, as universal newlines read it, and is
//! decoded as UTF-8 with `errors="replace"`. A line that starts with `#` is a header line
//! and is written as it is, its line end as `\n`. Every other line is a variant line: its
//! first column goes through the contig naming mode and the canonical-contig filter, and
//! the line is dropped or written with that column and a `\n`.
//!
//! `process` takes a run of whole lines. `Collector` takes the input in pieces cut
//! anywhere, holds back the line a piece ends inside, and hands on what is to be written
//! in input order; a piece whose whole lines another thread has already processed is
//! joined to its neighbours at its first and last line end.

use crate::contigs::ContigSet;
use memchr::{memchr, memchr2, memrchr};
use std::io;
use std::ops::Range;
use std::sync::Arc;

/// The contig naming a file's variant lines are given.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Keep,
    Strip,
    Add,
}

impl Mode {
    pub fn name(self) -> &'static str {
        match self {
            Mode::Keep => "keep",
            Mode::Strip => "strip",
            Mode::Add => "add",
        }
    }
}

/// What decides a variant line of one file.
pub struct Rules {
    pub mode: Mode,
    pub contigs: Option<Arc<ContigSet>>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub headers: u64,
    pub variants: u64,
    pub dropped: u64,
    pub chr_stripped: bool,
}

impl Counts {
    pub fn add(&mut self, other: &Counts) {
        self.headers += other.headers;
        self.variants += other.variants;
        self.dropped += other.dropped;
        self.chr_stripped |= other.chr_stripped;
    }
}

/// Whether a contig name starts with `chr` in any letter case.
pub fn has_chr_prefix(name: &[u8]) -> bool {
    name.len() >= 3 && name[..3].eq_ignore_ascii_case(b"chr")
}

/// What a variant line's first column comes to.
struct Outcome {
    /// The column to write where it is not the one read.
    renamed: Option<Vec<u8>>,
    keep: bool,
    stripped: bool,
}

impl Rules {
    /// The naming mode first; then, under a contig set, a name that still starts with
    /// `chr` and is not in the set loses the prefix, and a name not in the set drops the
    /// line.
    fn outcome(&self, chrom: &[u8]) -> Outcome {
        let prefixed = has_chr_prefix(chrom);
        let mut name: Vec<u8> = match self.mode {
            Mode::Strip if prefixed => chrom[3..].to_vec(),
            Mode::Add if !prefixed => [b"chr", chrom].concat(),
            _ => chrom.to_vec(),
        };
        let mut stripped = false;
        let mut keep = true;
        if let Some(contigs) = &self.contigs {
            if has_chr_prefix(&name) && !contigs.contains(&name) {
                name.drain(..3);
                stripped = true;
            }
            keep = contigs.contains(&name);
        }
        Outcome {
            renamed: (name != chrom).then_some(name),
            keep,
            stripped,
        }
    }
}

/// What a run of lines is written as.
#[derive(Debug, PartialEq, Eq)]
pub enum Text {
    /// The first bytes of the run, as they are.
    Same(usize),
    Rewritten(Vec<u8>),
}

/// Processes a run of whole lines; only its last line may lack a line end, and a `\r` that
/// ends the run is a line end. Under `cap`, stops before the first variant line met once
/// that many are written, and says so.
pub fn process(unit: &[u8], rules: &Rules, cap: Option<u64>, counts: &mut Counts) -> (Text, bool) {
    if std::str::from_utf8(unit).is_ok() {
        return process_decoded(unit, rules, cap, counts);
    }
    // Line ends are ASCII, so the replacement characters fall inside the lines.
    let decoded = String::from_utf8_lossy(unit).into_owned().into_bytes();
    match process_decoded(&decoded, rules, cap, counts) {
        (Text::Same(kept), stopped) => {
            let mut decoded = decoded;
            decoded.truncate(kept);
            (Text::Rewritten(decoded), stopped)
        }
        rewritten => rewritten,
    }
}

fn process_decoded(
    src: &[u8],
    rules: &Rules,
    cap: Option<u64>,
    counts: &mut Counts,
) -> (Text, bool) {
    let has_cr = memchr(b'\r', src).is_some();
    let mut out: Vec<u8> = Vec::new();
    // `out` is in use once a line is dropped or written otherwise than read; until then,
    // and between such lines, the bytes from `run_start` on are written as they are.
    let mut rewritten = false;
    let mut run_start = 0;
    let mut last: Option<(&[u8], Outcome)> = None;
    let mut stopped = false;
    let mut pos = 0;
    while pos < src.len() {
        // The line's text, where the next line starts, and whether its line end is `\n`.
        let (text_end, next, plain) = if has_cr {
            match memchr2(b'\n', b'\r', &src[pos..]) {
                Some(i) if src[pos + i] == b'\n' => (pos + i, pos + i + 1, true),
                Some(i) if src.get(pos + i + 1) == Some(&b'\n') => (pos + i, pos + i + 2, false),
                Some(i) => (pos + i, pos + i + 1, false),
                None => (src.len(), src.len(), false),
            }
        } else {
            match memchr(b'\n', &src[pos..]) {
                Some(i) => (pos + i, pos + i + 1, true),
                None => (src.len(), src.len(), false),
            }
        };
        let text = &src[pos..text_end];
        let ended = next > text_end;
        let mut rewrite = |parts: &[&[u8]]| {
            if !rewritten {
                rewritten = true;
                out.reserve(src.len());
            }
            out.extend_from_slice(&src[run_start..pos]);
            for part in parts {
                out.extend_from_slice(part);
            }
            run_start = next;
        };
        if text.first() == Some(&b'#') {
            counts.headers += 1;
            if ended && !plain {
                rewrite(&[text, b"\n"]);
            }
            pos = next;
            continue;
        }
        if cap.is_some_and(|cap| counts.variants >= cap) {
            stopped = true;
            break;
        }
        let tab = memchr(b'\t', text).unwrap_or(text.len());
        let chrom = &text[..tab];
        if last.as_ref().map(|(name, _)| *name) != Some(chrom) {
            last = Some((chrom, rules.outcome(chrom)));
        }
        let outcome = &last.as_ref().expect("set above").1;
        counts.chr_stripped |= outcome.stripped;
        if !outcome.keep {
            counts.dropped += 1;
            rewrite(&[]);
        } else {
            counts.variants += 1;
            match &outcome.renamed {
                Some(name) => rewrite(&[name, &text[tab..], b"\n"]),
                None if !plain => rewrite(&[text, b"\n"]),
                None => {}
            }
        }
        pos = next;
    }
    if rewritten {
        out.extend_from_slice(&src[run_start..pos]);
        (Text::Rewritten(out), stopped)
    } else {
        (Text::Same(pos), stopped)
    }
}

/// The end of the last line of `bytes` whose line end is certain: after the last `\n`,
/// else after the last `\r` that a byte follows.
fn last_line_end(bytes: &[u8]) -> Option<usize> {
    if let Some(at) = memrchr(b'\n', bytes) {
        return Some(at + 1);
    }
    let before_last = bytes.len().checked_sub(1)?;
    memrchr(b'\r', &bytes[..before_last]).map(|at| at + 1)
}

/// Bytes to write: a buffer, or a part of one.
pub enum Chunk {
    Whole(Vec<u8>),
    Part(Vec<u8>, Range<usize>),
}

impl Chunk {
    pub fn bytes(&self) -> &[u8] {
        match self {
            Chunk::Whole(bytes) => bytes,
            Chunk::Part(bytes, range) => &bytes[range.clone()],
        }
    }
}

/// A piece of the input whose whole lines are processed: the bytes read, the end of the
/// first line end in them, the end of the last, what the lines between the two are
/// written as where that is not those bytes, and their counts.
pub struct Split {
    pub piece: Vec<u8>,
    pub head_end: usize,
    pub tail_start: usize,
    pub middle: Option<Vec<u8>>,
    pub counts: Counts,
}

/// Processes the lines of `piece` that start and end inside it; None where it holds no
/// `\n` to cut at.
pub fn split(piece: Vec<u8>, rules: &Rules) -> Result<Split, Vec<u8>> {
    let Some(first) = memchr(b'\n', &piece) else {
        return Err(piece);
    };
    let (head_end, tail_start) = (first + 1, memrchr(b'\n', &piece).expect("found above") + 1);
    let mut counts = Counts::default();
    let middle = match process(&piece[head_end..tail_start], rules, None, &mut counts).0 {
        Text::Same(_) => None,
        Text::Rewritten(text) => Some(text),
    };
    Ok(Split {
        piece,
        head_end,
        tail_start,
        middle,
        counts,
    })
}

/// Joins the pieces of one input in order and hands what is to be written to `emit`.
pub struct Collector<'a, E> {
    rules: &'a Rules,
    cap: Option<u64>,
    emit: E,
    /// The end of what was taken, which may be part of a line that goes on.
    held: Vec<u8>,
    pub counts: Counts,
}

impl<'a, E: FnMut(Chunk) -> io::Result<()>> Collector<'a, E> {
    pub fn new(rules: &'a Rules, cap: Option<u64>, emit: E) -> Collector<'a, E> {
        Collector {
            rules,
            cap,
            emit,
            held: Vec::new(),
            counts: Counts::default(),
        }
    }

    /// Processes a run of whole lines and emits it; true where the cap stopped it.
    fn unit(&mut self, unit: Vec<u8>) -> io::Result<bool> {
        let (text, stopped) = process(&unit, self.rules, self.cap, &mut self.counts);
        match text {
            Text::Same(0) => {}
            Text::Same(kept) => (self.emit)(Chunk::Part(unit, 0..kept))?,
            Text::Rewritten(text) if text.is_empty() => {}
            Text::Rewritten(text) => (self.emit)(Chunk::Whole(text))?,
        }
        Ok(stopped)
    }

    /// Takes the next bytes of the input, cut anywhere; true where the cap stopped the
    /// input, after which nothing more is taken.
    pub fn raw(&mut self, mut bytes: Vec<u8>) -> io::Result<bool> {
        if !self.held.is_empty() {
            // The held line ends at the first `\n` of the new bytes, if they hold one.
            let joined = memchr(b'\n', &bytes).map_or(bytes.len(), |at| at + 1);
            self.held.extend_from_slice(&bytes[..joined]);
            bytes.drain(..joined);
            if let Some(cut) = last_line_end(&self.held) {
                let rest = self.held.split_off(cut);
                let unit = std::mem::replace(&mut self.held, rest);
                if self.unit(unit)? {
                    return Ok(true);
                }
            }
            if bytes.is_empty() {
                return Ok(false);
            }
        }
        match last_line_end(&bytes) {
            Some(cut) => {
                self.held = bytes[cut..].to_vec();
                bytes.truncate(cut);
                self.unit(bytes)
            }
            None => {
                self.held = bytes;
                Ok(false)
            }
        }
    }

    /// Takes the next piece of the input, its whole lines already processed. Not for an
    /// input under a cap, whose lines are counted in order.
    pub fn split(&mut self, split: Split) -> io::Result<()> {
        let Split {
            mut piece,
            head_end,
            tail_start,
            middle,
            counts,
        } = split;
        self.held.extend_from_slice(&piece[..head_end]);
        let unit = std::mem::take(&mut self.held);
        self.unit(unit)?;
        self.held = piece[tail_start..].to_vec();
        match middle {
            Some(text) if text.is_empty() => {}
            Some(text) => (self.emit)(Chunk::Whole(text))?,
            None if head_end == tail_start => {}
            None => {
                piece.truncate(tail_start);
                (self.emit)(Chunk::Part(piece, head_end..tail_start))?;
            }
        }
        self.counts.add(&counts);
        Ok(())
    }

    /// Ends the input: the held bytes are its last line.
    pub fn finish(&mut self) -> io::Result<()> {
        let unit = std::mem::take(&mut self.held);
        self.unit(unit).map(|_| ())
    }
}

#[cfg(test)]
#[path = "lines_tests.rs"]
mod tests;
