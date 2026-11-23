// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Digest lines as `sha256sum`, `md5sum` and `xxh128sum` print them.
//!
//! A file digested whole is read on one thread, and the files of a run are read at once on
//! as many threads as the run has; their lines are written in the order the names were
//! given. A body (a file inflated, a file without the lines that start with a text, the
//! bytes after a file's leading lines) is digested as the digest program reads it from a
//! pipe, so its line names standard input: the inflating and the dropping of lines run
//! across the threads and one thread digests in file order. A body digested in chunks is
//! the digest of its chunks' digests, one hexadecimal digest a line in chunk order, the
//! chunks read and digested across the threads.

use crate::pipeline::{lines, ordered, pieces, Ends, Reading};
use memchr::memmem;
use sha2::{Digest, Sha256};
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::FileExt;
use std::sync::mpsc::sync_channel;
use std::sync::Arc;
use std::thread;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Algorithm {
    Sha256,
    Md5,
    Xxh128,
}

impl Algorithm {
    pub fn parse(word: &str) -> Option<Algorithm> {
        match word {
            "sha256" => Some(Algorithm::Sha256),
            "md5" => Some(Algorithm::Md5),
            "xxh128" => Some(Algorithm::Xxh128),
            _ => None,
        }
    }

    /// The name the digest program prints for its standard input.
    fn stdin_name(self) -> &'static [u8] {
        match self {
            Algorithm::Xxh128 => b"stdin",
            _ => b"-",
        }
    }

    /// Whether the digest program escapes a name that holds a backslash or a newline.
    fn escapes(self) -> bool {
        !matches!(self, Algorithm::Xxh128)
    }
}

enum State {
    Sha256(Box<Sha256>),
    Md5(Box<md5::Context>),
    Xxh128(Box<twox_hash::XxHash3_128>),
}

impl State {
    fn new(algorithm: Algorithm) -> State {
        match algorithm {
            Algorithm::Sha256 => State::Sha256(Box::new(Sha256::new())),
            Algorithm::Md5 => State::Md5(Box::new(md5::Context::new())),
            Algorithm::Xxh128 => State::Xxh128(Box::new(twox_hash::XxHash3_128::new())),
        }
    }

    fn update(&mut self, bytes: &[u8]) {
        match self {
            State::Sha256(state) => state.update(bytes),
            State::Md5(state) => state.consume(bytes),
            State::Xxh128(state) => state.write(bytes),
        }
    }

    fn hex(self) -> String {
        match self {
            State::Sha256(state) => state
                .finalize()
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect(),
            State::Md5(state) => format!("{:x}", state.finalize()),
            State::Xxh128(state) => format!("{:032x}", state.finish_128()),
        }
    }
}

/// What of a file is digested.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Body {
    /// Every byte.
    Whole,
    /// The lines that do not start with the text, each ended by a newline, as
    /// `grep -a -v '^TEXT'` writes them.
    WithoutLinesStarting(Vec<u8>),
    /// The bytes after the file's leading lines that start with the text.
    AfterLeadingLinesStarting(Vec<u8>),
}

#[derive(Clone, Debug)]
pub struct Plan {
    pub algorithms: Vec<Algorithm>,
    pub gunzip: bool,
    pub body: Body,
    pub chunk_bytes: Option<u64>,
}

impl Plan {
    /// Whether the digest program would read the bytes from a pipe.
    fn piped(&self) -> bool {
        self.gunzip || self.body != Body::Whole || self.chunk_bytes.is_some()
    }
}

/// How one name came out.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Status {
    Digested,
    /// Digested, with bytes after the last gzip member that start no member.
    Warned,
    /// Digested, a body without a line: `grep -v` leaves with status 1 where it writes none.
    NoLine,
    Failed,
}

/// The lines of one name and what went wrong with it.
pub struct Outcome {
    pub lines: Vec<u8>,
    pub status: Status,
    pub complaint: Option<String>,
}

fn line(hex: &str, name: &[u8], escapes: bool) -> Vec<u8> {
    let mut out = Vec::with_capacity(hex.len() + name.len() + 4);
    let escaped = escapes && name.iter().any(|byte| matches!(byte, b'\\' | b'\n'));
    if escaped {
        out.push(b'\\');
    }
    out.extend_from_slice(hex.as_bytes());
    out.extend_from_slice(b"  ");
    for byte in name {
        match byte {
            b'\\' if escaped => out.extend_from_slice(b"\\\\"),
            b'\n' if escaped => out.extend_from_slice(b"\\n"),
            other => out.push(*other),
        }
    }
    out.push(b'\n');
    out
}

/// How much of a file one read asks for.
const READ_BYTES: usize = 1 << 20;

/// Digests what `source` reads under every algorithm, the second and further ones on
/// threads of their own.
fn digest_stream(mut source: impl Read, algorithms: &[Algorithm]) -> io::Result<Vec<String>> {
    let mut first = State::new(algorithms[0]);
    if algorithms.len() == 1 {
        let mut buf = vec![0u8; READ_BYTES];
        loop {
            match source.read(&mut buf) {
                Ok(0) => return Ok(vec![first.hex()]),
                Ok(count) => first.update(&buf[..count]),
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                Err(error) => return Err(error),
            }
        }
    }
    thread::scope(|scope| {
        let mut feeds = Vec::new();
        let mut others = Vec::new();
        for algorithm in &algorithms[1..] {
            let (feed, fed) = sync_channel::<Arc<Vec<u8>>>(4);
            feeds.push(feed);
            others.push(scope.spawn(move || {
                let mut state = State::new(*algorithm);
                for bytes in fed {
                    state.update(&bytes);
                }
                state.hex()
            }));
        }
        let read = (|| -> io::Result<()> {
            loop {
                let mut bytes = Vec::with_capacity(READ_BYTES);
                let count = (&mut source)
                    .take(READ_BYTES as u64)
                    .read_to_end(&mut bytes)?;
                if count == 0 {
                    return Ok(());
                }
                let bytes = Arc::new(bytes);
                for feed in &feeds {
                    feed.send(Arc::clone(&bytes))
                        .map_err(|_| crate::pipeline::lost_thread())?;
                }
                first.update(&bytes);
            }
        })();
        drop(feeds);
        let mut digests = vec![first.hex()];
        for other in others {
            digests.push(other.join().map_err(|_| crate::pipeline::lost_thread())?);
        }
        read.map(|()| digests)
    })
}

/// A file of this many bytes or more is read ahead of its digest.
const READ_AHEAD_BYTES: u64 = 32 << 20;
/// The threads that read ahead of one file's digest.
const READ_AHEAD_THREADS: usize = 2;

/// Digests a file whole. One digest of a long file is fed by threads that read ahead of
/// it, so that the digesting thread does nothing but digest.
fn digest_file(file: &File, algorithms: &[Algorithm]) -> io::Result<Vec<String>> {
    let long = file
        .metadata()
        .is_ok_and(|meta| meta.is_file() && meta.len() >= READ_AHEAD_BYTES);
    if !long || algorithms.len() > 1 {
        return digest_stream(file, algorithms);
    }
    let mut state = State::new(algorithms[0]);
    pieces(
        file,
        Reading::Plain { from: 0 },
        READ_AHEAD_THREADS,
        |bytes| bytes,
        |bytes| {
            state.update(&bytes);
            Ok(())
        },
    )?;
    Ok(vec![state.hex()])
}

/// The bytes of the leading lines of `file` that start with `text`.
fn leading_lines(file: &File, text: &[u8]) -> io::Result<u64> {
    let mut reader = BufReader::with_capacity(1 << 16, file);
    let mut total = 0u64;
    loop {
        if reader.fill_buf()?.is_empty() {
            return Ok(total);
        }
        let mut start = vec![0u8; text.len()];
        let mut filled = 0;
        while filled < start.len() {
            let count = reader.read(&mut start[filled..])?;
            if count == 0 {
                break;
            }
            filled += count;
        }
        if start[..filled] != *text {
            return Ok(total);
        }
        total += filled as u64;
        total += reader.skip_until(b'\n')? as u64;
    }
}

/// `data` without its lines that start with `text`; a last line without a newline gains
/// one.
fn without_lines_starting(data: &[u8], text: &[u8]) -> Vec<u8> {
    let mut kept = Vec::with_capacity(data.len());
    let mut needle = Vec::with_capacity(text.len() + 1);
    needle.push(b'\n');
    needle.extend_from_slice(text);
    let finder = memmem::Finder::new(&needle);
    // Where the next line to look at starts, and where the lines still to keep start.
    let mut at = 0;
    let mut keep_from = 0;
    while at < data.len() {
        let dropped = if data[at..].starts_with(text) {
            Some(at)
        } else {
            finder.find(&data[at..]).map(|found| at + found + 1)
        };
        let Some(start) = dropped else {
            break;
        };
        kept.extend_from_slice(&data[keep_from..start]);
        let end = memchr::memchr(b'\n', &data[start..]).map_or(data.len(), |nl| start + nl + 1);
        at = end;
        keep_from = end;
    }
    kept.extend_from_slice(&data[keep_from..]);
    if kept.last().is_some_and(|byte| *byte != b'\n') {
        kept.push(b'\n');
    }
    kept
}

struct Bodies<'a> {
    plan: &'a Plan,
    threads: usize,
}

impl Bodies<'_> {
    /// The digests of the bytes after `offset` of a plain file, whole.
    fn rest(&self, file: &File, offset: u64) -> io::Result<Vec<String>> {
        let mut states: Vec<State> = self
            .plan
            .algorithms
            .iter()
            .map(|a| State::new(*a))
            .collect();
        pieces(
            file,
            Reading::Plain { from: offset },
            self.threads,
            |bytes| bytes,
            |bytes| {
                for state in &mut states {
                    state.update(&bytes);
                }
                Ok(())
            },
        )?;
        Ok(states.into_iter().map(State::hex).collect())
    }

    /// The digest of the digests of the chunks of the bytes after `offset` of a plain
    /// file, under one algorithm.
    fn chunked(
        &self,
        file: &File,
        offset: u64,
        chunk: u64,
        algorithm: Algorithm,
    ) -> io::Result<String> {
        let size = file.metadata()?.len();
        let body = size.saturating_sub(offset);
        let count = body.div_ceil(chunk);
        let mut fold = State::new(algorithm);
        let one = |index: u64| -> io::Result<String> {
            let mut state = State::new(algorithm);
            let mut at = offset + index * chunk;
            let end = (at + chunk).min(size);
            let mut buf = vec![0u8; ((end - at) as usize).min(4 << 20)];
            while at < end {
                let want = ((end - at) as usize).min(buf.len());
                file.read_exact_at(&mut buf[..want], at)?;
                state.update(&buf[..want]);
                at += want as u64;
            }
            Ok(state.hex())
        };
        ordered(
            self.threads,
            |emit| {
                for index in 0..count {
                    if !emit(index) {
                        break;
                    }
                }
            },
            one,
            |hex: io::Result<String>| {
                fold.update(hex?.as_bytes());
                fold.update(b"\n");
                Ok(())
            },
        )?;
        Ok(fold.hex())
    }

    fn digests(&self, file: &File) -> (Vec<String>, Status, Option<io::Error>) {
        let plan = self.plan;
        let reading = if plan.gunzip {
            Reading::Gunzip
        } else {
            Reading::Plain { from: 0 }
        };
        let mut states: Vec<State> = plan.algorithms.iter().map(|a| State::new(*a)).collect();
        let mut fed = false;
        let mut feed = |bytes: Vec<u8>| -> io::Result<()> {
            fed |= !bytes.is_empty();
            for state in &mut states {
                state.update(&bytes);
            }
            Ok(())
        };
        let finish = match &plan.body {
            Body::Whole => pieces(file, reading, self.threads, |bytes| bytes, &mut feed),
            Body::WithoutLinesStarting(text) => lines(
                file,
                reading,
                self.threads,
                Ends::Newline,
                |data| without_lines_starting(data, text),
                &mut feed,
            ),
            Body::AfterLeadingLinesStarting(_) => unreachable!("read from an offset"),
        };
        let digests = states.into_iter().map(State::hex).collect();
        let no_line = !fed && matches!(plan.body, Body::WithoutLinesStarting(_));
        match finish {
            Ok(_) if no_line => (digests, Status::NoLine, None),
            Ok(finish) if finish.garbage => (digests, Status::Warned, None),
            Ok(_) => (digests, Status::Digested, None),
            Err(error) => (digests, Status::Failed, Some(error)),
        }
    }
}

fn open(name: &[u8]) -> io::Result<File> {
    File::open(OsStr::from_bytes(name))
}

fn shown(name: &[u8]) -> String {
    String::from_utf8_lossy(name).into_owned()
}

/// The lines of a file digested whole, as the digest program names it.
fn whole(name: &[u8], algorithms: &[Algorithm]) -> Outcome {
    let digests = if name == b"-" {
        digest_stream(io::stdin().lock(), algorithms)
    } else {
        open(name).and_then(|file| digest_file(&file, algorithms))
    };
    match digests {
        Ok(digests) => Outcome {
            lines: algorithms
                .iter()
                .zip(&digests)
                .flat_map(|(algorithm, hex)| line(hex, name, algorithm.escapes()))
                .collect(),
            status: Status::Digested,
            complaint: None,
        },
        Err(error) => Outcome {
            lines: Vec::new(),
            status: Status::Failed,
            complaint: Some(format!("{}: {error}", shown(name))),
        },
    }
}

/// The lines of a body of one file, named as the digest program names standard input.
fn body(name: &[u8], bodies: &Bodies<'_>) -> Outcome {
    let plan = bodies.plan;
    let named = |digests: &[String]| -> Vec<u8> {
        plan.algorithms
            .iter()
            .zip(digests)
            .flat_map(|(algorithm, hex)| line(hex, algorithm.stdin_name(), false))
            .collect()
    };
    let failed = |lines: Vec<u8>, error: io::Error| Outcome {
        lines,
        status: Status::Failed,
        complaint: Some(format!("{}: {error}", shown(name))),
    };
    let empty: Vec<String> = plan
        .algorithms
        .iter()
        .map(|algorithm| State::new(*algorithm).hex())
        .collect();
    let file = match open(name) {
        Ok(file) => file,
        // The digest program reads nothing from the pipe and digests that.
        Err(error) => return failed(named(&empty), error),
    };
    let offset = match &plan.body {
        Body::AfterLeadingLinesStarting(text) => match leading_lines(&file, text) {
            Ok(offset) => Some(offset),
            Err(error) => return failed(named(&empty), error),
        },
        _ if plan.chunk_bytes.is_some() => Some(0),
        _ => None,
    };
    let Some(offset) = offset else {
        let (digests, status, error) = bodies.digests(&file);
        return Outcome {
            lines: named(&digests),
            status,
            complaint: error.map(|error| format!("{}: {error}", shown(name))),
        };
    };
    let digests = match plan.chunk_bytes {
        Some(chunk) => plan
            .algorithms
            .iter()
            .map(|algorithm| bodies.chunked(&file, offset, chunk, *algorithm))
            .collect::<io::Result<Vec<String>>>(),
        None => bodies.rest(&file, offset),
    };
    match digests {
        Ok(digests) => Outcome {
            lines: named(&digests),
            status: Status::Digested,
            complaint: None,
        },
        Err(error) => failed(Vec::new(), error),
    }
}

/// Digests `names` in order and writes their lines to `out`; each name's status, in order.
/// `on_no_name` asks for the one line a digest program started on no name writes of its
/// empty standard input.
pub fn run(
    plan: &Plan,
    names: &[Vec<u8>],
    on_no_name: bool,
    threads: usize,
    out: &mut impl Write,
    complain: &mut impl FnMut(&str),
) -> io::Result<Vec<Status>> {
    if names.is_empty() && on_no_name {
        for algorithm in &plan.algorithms {
            let hex = State::new(*algorithm).hex();
            out.write_all(&line(&hex, algorithm.stdin_name(), false))?;
        }
        return Ok(Vec::new());
    }
    let mut statuses = Vec::with_capacity(names.len());
    let mut write = |outcome: Outcome| -> io::Result<()> {
        out.write_all(&outcome.lines)?;
        if let Some(complaint) = &outcome.complaint {
            complain(complaint);
        }
        statuses.push(outcome.status);
        Ok(())
    };
    if plan.piped() {
        let bodies = Bodies { plan, threads };
        for name in names {
            write(body(name, &bodies))?;
        }
    } else {
        ordered(
            threads.min(names.len().max(1)),
            |emit| {
                for name in names {
                    if !emit(name) {
                        break;
                    }
                }
            },
            |name: &Vec<u8>| whole(name, &plan.algorithms),
            &mut write,
        )?;
    }
    Ok(statuses)
}

#[cfg(test)]
#[path = "digest_tests.rs"]
mod tests;
