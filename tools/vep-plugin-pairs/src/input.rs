// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One input as the passes read it: its digest as stored, a plain file to read in blocks
//! that each start at a line, and whether its lines are default-format rows or tuples.
//!
//! A gzip input is unpacked once into the scratch, every member of it, so that each pass
//! reads any block of it without unpacking what comes before.

use crate::rules::{Refusal, Res};
use crate::score::interrupt;
use flate2::bufread::GzDecoder;
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::OnceLock;

pub const DEFAULT_COLUMNS: usize = 14;
pub const TUPLE_COLUMNS: usize = 5;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Default,
    Tuples,
}

impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Kind::Default => "default",
            Kind::Tuples => "tuples",
        }
    }
}

pub struct Input {
    /// The path as the command line gave it, for the messages.
    pub shown: String,
    pub tuple_file_allowed: bool,
    pub kind: Kind,
    pub file: File,
    /// Where each block starts, and the file's length last.
    pub bounds: Vec<u64>,
    /// The sha256 of the file as stored, once it is taken.
    digest: OnceLock<io::Result<String>>,
}

/// Digests every byte read through it.
struct Digested<R> {
    inner: R,
    digest: Sha256,
}

impl<R: Read> Read for Digested<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let count = self.inner.read(buffer)?;
        self.digest.update(&buffer[..count]);
        Ok(count)
    }
}

fn hex_digest(digest: Sha256) -> String {
    crate::order::hex(&digest.finalize())
}

/// Unpacks every member of the gzip file into `plain` and returns the digest of the file as
/// stored. Zero bytes after a member are padding, as Python's `gzip` reads them.
fn unpack(shown: &str, packed: File, plain: &Path) -> Res<String> {
    let broken =
        |error: io::Error| Refusal(format!("{shown}: cannot read the gzip stream: {error}"));
    let mut stored = BufReader::with_capacity(
        1 << 20,
        Digested {
            inner: packed,
            digest: Sha256::new(),
        },
    );
    let mut unpacked = BufWriter::with_capacity(1 << 20, File::create(plain)?);
    loop {
        let mut member = GzDecoder::new(&mut stored);
        let mut block = vec![0u8; 1 << 20];
        loop {
            let count = member.read(&mut block).map_err(broken)?;
            if count == 0 {
                break;
            }
            unpacked.write_all(&block[..count])?;
        }
        loop {
            let padding = stored.fill_buf().map_err(broken)?;
            let zeros = padding.iter().take_while(|&&byte| byte == 0).count();
            if zeros == 0 {
                break;
            }
            stored.consume(zeros);
        }
        if stored.fill_buf().map_err(broken)?.is_empty() {
            break;
        }
    }
    unpacked.flush()?;
    Ok(hex_digest(stored.into_inner().digest))
}

/// What the column count of the first data line makes of a file that may be a tuple file;
/// a file with no data line is a default-format output where it holds a `#` line and a
/// tuple file otherwise.
fn kind_of(file: &File) -> io::Result<Kind> {
    let mut lines = BufReader::with_capacity(1 << 16, file);
    let mut line = Vec::new();
    let mut headed = false;
    loop {
        line.clear();
        if lines.read_until(b'\n', &mut line)? == 0 {
            return Ok(if headed { Kind::Default } else { Kind::Tuples });
        }
        while line.last().is_some_and(|byte| b"\r\n".contains(byte)) {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }
        if line.starts_with(b"#") {
            headed = true;
            continue;
        }
        let columns = memchr::memchr_iter(b'\t', &line).count() + 1;
        return Ok(if columns == TUPLE_COLUMNS && !headed {
            Kind::Tuples
        } else {
            Kind::Default
        });
    }
}

/// The start of every block: the first line that starts at or after each multiple of
/// `block_bytes`, each once.
fn block_bounds(file: &File, length: u64, block_bytes: u64) -> io::Result<Vec<u64>> {
    let mut bounds = vec![0];
    let mut window = vec![0u8; 1 << 16];
    let mut target = block_bytes;
    while target < length {
        // A line starts at `target` when the byte before it ends one.
        let mut at = target - 1;
        let start = loop {
            let count = file.read_at(&mut window, at)?;
            if count == 0 {
                break length;
            }
            if let Some(found) = memchr::memchr(b'\n', &window[..count]) {
                break at + found as u64 + 1;
            }
            at += count as u64;
        };
        if start >= length {
            break;
        }
        if start > *bounds.last().unwrap() {
            bounds.push(start);
        }
        target = (start / block_bytes + 1) * block_bytes;
    }
    bounds.push(length);
    Ok(bounds)
}

impl Input {
    /// Opens `path`, unpacking a gzip file to `plain`, which is digested as it is unpacked;
    /// a plain file's digest is taken by `digest`.
    pub fn open(
        path: &Path,
        tuple_file_allowed: bool,
        plain: &Path,
        block_bytes: u64,
    ) -> Res<Input> {
        let shown = path.display().to_string();
        let named = |error: io::Error| Refusal(format!("{shown}: {error}"));
        let stored = File::open(path).map_err(named)?;
        let mut magic = [0u8; 2];
        let packed = stored.read_at(&mut magic, 0).map_err(named)? == 2 && magic == [0x1f, 0x8b];
        let (file, digest) = if packed {
            let digest = unpack(&shown, stored, plain)?;
            (File::open(plain)?, OnceLock::from(Ok(digest)))
        } else {
            (stored, OnceLock::new())
        };
        let kind = if tuple_file_allowed {
            kind_of(&file).map_err(named)?
        } else {
            Kind::Default
        };
        let length = file.metadata().map_err(named)?.len();
        let bounds = block_bounds(&file, length, block_bytes).map_err(named)?;
        Ok(Input {
            shown,
            tuple_file_allowed,
            kind,
            file,
            bounds,
            digest,
        })
    }

    pub fn blocks(&self) -> usize {
        self.bounds.len() - 1
    }

    pub fn length(&self) -> u64 {
        *self.bounds.last().unwrap()
    }

    /// Reads block `index` into `buffer`.
    pub fn read_block(&self, index: usize, buffer: &mut Vec<u8>) -> io::Result<()> {
        let (start, end) = (self.bounds[index], self.bounds[index + 1]);
        buffer.resize((end - start) as usize, 0);
        self.file.read_exact_at(buffer, start)
    }

    /// Whether the sha256 of the file as stored is still to be taken.
    pub fn digest_is_due(&self) -> bool {
        self.digest.get().is_none()
    }

    /// Takes the sha256 of a plain file as stored, unless it is known. The reading ends with
    /// no digest taken once `stop` is set or a stop signal arrived.
    pub fn digest(&self, stop: &AtomicBool) {
        if !self.digest_is_due() {
            return;
        }
        let mut digest = Sha256::new();
        let mut block = vec![0u8; 1 << 20];
        let mut at = 0;
        let taken = loop {
            if stop.load(Ordering::Relaxed) || interrupt::requested() {
                return;
            }
            match self.file.read_at(&mut block, at) {
                Ok(0) => break Ok(hex_digest(digest)),
                Ok(count) => {
                    digest.update(&block[..count]);
                    at += count as u64;
                }
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                Err(error) => break Err(error),
            }
        };
        let _ = self.digest.set(taken);
    }

    /// The sha256 of the file as stored, taken here where it is still due.
    pub fn sha256(&self) -> Res<&str> {
        self.digest(&AtomicBool::new(false));
        match self.digest.get() {
            Some(Ok(digest)) => Ok(digest),
            Some(Err(error)) => Err(Refusal(format!("{}: {error}", self.shown))),
            // A stop signal ended the reading, and the run ends on the signal.
            None => Err(Refusal(format!("{}: the digest was not taken", self.shown))),
        }
    }
}

/// The paths of one input's files in the scratch.
pub fn scratch_file(scratch: &Path, input: usize, part: &str) -> PathBuf {
    scratch.join(format!("input{input}.{part}"))
}

/// The lines of a block, each without its line end: every trailing `\r` and `\n`, as the
/// scorer strips them.
pub fn lines(block: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut rest = block;
    std::iter::from_fn(move || {
        if rest.is_empty() {
            return None;
        }
        let end = memchr::memchr(b'\n', rest).unwrap_or(rest.len());
        let mut line = &rest[..end];
        rest = &rest[(end + 1).min(rest.len())..];
        while let [kept @ .., b'\r'] = line {
            line = kept;
        }
        Some(line)
    })
}
