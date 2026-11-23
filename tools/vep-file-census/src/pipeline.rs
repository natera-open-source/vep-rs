// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One file from its bytes to what a pass makes of them, across threads.
//!
//! `ordered` runs a source on one thread, hands every item it emits to a worker thread
//! and gives the calling thread what the workers made, in the order the source emitted.
//! `pieces` is that over a file: a plain file is cut into ranges that the worker threads
//! read themselves, a BGZF file into runs of whole blocks that the worker threads inflate,
//! and any other gzip file is inflated member by member on the source's thread. `lines`
//! is `pieces` for a pass over lines: the whole lines inside a piece are processed on the
//! worker threads, and the line that crosses from one piece to the next is joined and
//! processed on the calling thread. No more of a file is held than the pieces in flight
//! and the line being joined.

use crate::gzip::{bgzf_block_in_bound, bgzf_block_len, inflate_blocks, Members, BGZF_HEADER};
use flate2::Decompress;
use memchr::{memchr, memchr2, memrchr, memrchr2};
use std::fs::File;
use std::io::{self, BufReader, Cursor, Read};
use std::os::unix::fs::FileExt;
use std::sync::mpsc::{channel, sync_channel, Receiver, Sender};
use std::sync::Mutex;
use std::thread;

/// The environment variable that bounds the worker threads; unset or empty, the host's
/// available cores are the bound. An environment variable and not an option: the command
/// line names the pass, and the thread count changes no byte of it.
pub const THREADS_ENV: &str = "VEP_FILE_CENSUS_THREADS";

/// The worker threads of this run.
pub fn thread_bound() -> Result<usize, String> {
    match std::env::var(THREADS_ENV) {
        Ok(value) if !value.is_empty() => match value.parse::<usize>() {
            Ok(count) if count >= 1 => Ok(count),
            _ => Err(format!(
                "{THREADS_ENV} takes a whole number of 1 or more, not '{value}'"
            )),
        },
        _ => Ok(thread::available_parallelism().map_or(1, |count| count.get())),
    }
}

pub fn lost_thread() -> io::Error {
    io::Error::other("a thread of the run ended early")
}

/// Runs `source` on a thread of its own, `work` on `threads` threads over every item the
/// source emits, and `each` on the calling thread over what `work` made, in the order the
/// source emitted. The source's emitter answers false once the calling thread takes no
/// more. Returns what the source returned, or the failure of `each`.
pub fn ordered<I, T, R>(
    threads: usize,
    source: impl FnOnce(&mut dyn FnMut(I) -> bool) -> R + Send,
    work: impl Fn(I) -> T + Sync,
    mut each: impl FnMut(T) -> io::Result<()>,
) -> io::Result<R>
where
    I: Send,
    T: Send,
    R: Send,
{
    let threads = threads.max(1);
    let (jobs, queue) = sync_channel::<(I, Sender<T>)>(threads);
    let queue = Mutex::new(queue);
    let (slots, taken) = sync_channel::<Receiver<T>>(threads * 2);
    thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| loop {
                let job = queue.lock().expect("no job runs under the lock").recv();
                match job {
                    Ok((item, done)) => {
                        let _ = done.send(work(item));
                    }
                    Err(_) => return,
                }
            });
        }
        let reader = scope.spawn(move || {
            let mut emit = |item: I| -> bool {
                let (done, slot) = channel();
                slots.send(slot).is_ok() && jobs.send((item, done)).is_ok()
            };
            source(&mut emit)
        });
        let taken = taken;
        let mut failure = None;
        for slot in &taken {
            let made = match slot.recv() {
                Ok(made) => made,
                Err(_) => {
                    failure = Some(lost_thread());
                    break;
                }
            };
            if let Err(error) = each(made) {
                failure = Some(error);
                break;
            }
        }
        // The source's next emit fails once the queue of slots is gone.
        drop(taken);
        let result = reader.join().map_err(|_| lost_thread());
        match failure {
            Some(error) => Err(error),
            None => result,
        }
    })
}

/// How a file is read.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Reading {
    /// The file's own bytes, from an offset on.
    Plain { from: u64 },
    /// What `gzip -dc` writes of it.
    Gunzip,
    /// What `zcat -f` writes of it: inflated where it starts as a gzip file, else as it is.
    Auto,
}

/// How a reading ended where it did not fail.
#[derive(Clone, Copy, Default, Debug)]
pub struct Finish {
    /// Bytes after the last gzip member that start no member: `gzip` warns and exits 2.
    pub garbage: bool,
}

/// The BGZF blocks of one piece; they inflate to about two mebibytes.
const BLOCKS_PER_PIECE: usize = 32;
/// The bytes of one piece of a file that is not BGZF.
const PIECE_BYTES: usize = 2 << 20;
/// How much of a gzip file the source asks for at once.
const READ_BYTES: usize = 1 << 20;

enum Raw {
    Bytes(Vec<u8>),
    Blocks(Vec<u8>),
    Range(u64, usize),
}

/// Fills `buf` from `offset` of `file`; the bytes read where the file ends first.
fn read_range(file: &File, mut offset: u64, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match file.read_at(&mut buf[filled..], offset) {
            Ok(0) => break,
            Ok(count) => {
                filled += count;
                offset += count as u64;
            }
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(filled)
}

/// Emits what `stream` reads, in pieces, to its end; what it read before a failure is
/// emitted before the failure is returned.
fn stream(mut stream: impl Read, emit: &mut dyn FnMut(Raw) -> bool) -> io::Result<()> {
    loop {
        let mut bytes = Vec::with_capacity(PIECE_BYTES);
        let read = (&mut stream)
            .take(PIECE_BYTES as u64)
            .read_to_end(&mut bytes);
        if !bytes.is_empty() && !emit(Raw::Bytes(bytes)) {
            return Ok(());
        }
        if read? == 0 {
            return Ok(());
        }
    }
}

/// A file read from its first byte on, wherever its own position stands; a file that
/// cannot be read at an offset, a pipe, is read from where it stands.
struct FromStart<'a> {
    file: &'a File,
    offset: Option<u64>,
}

impl Read for FromStart<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if let Some(offset) = self.offset {
            match self.file.read_at(buf, offset) {
                Ok(count) => {
                    self.offset = Some(offset + count as u64);
                    return Ok(count);
                }
                Err(error) if error.raw_os_error() == Some(ESPIPE) && offset == 0 => {
                    self.offset = None;
                }
                Err(error) => return Err(error),
            }
        }
        let mut file = self.file;
        file.read(buf)
    }
}

/// `errno` of a read at an offset from a file that has none.
const ESPIPE: i32 = 29;

/// Emits a gzip file: its leading BGZF blocks in runs, and from the first byte that does
/// not start a whole BGZF block on, what the member reader reads. A member under a block's
/// header that announces more than a block's bound starts no block.
fn gzip(file: &File, pass_through: bool, emit: &mut dyn FnMut(Raw) -> bool) -> io::Result<Finish> {
    let mut source = FromStart {
        file,
        offset: Some(0),
    };
    // The bytes of the file that are read and still to be emitted.
    let mut unread: Vec<u8> = Vec::new();
    let mut first = true;
    let mut at_end = false;
    loop {
        let (mut end, mut blocks) = (0, 0);
        while blocks < BLOCKS_PER_PIECE {
            let length = bgzf_block_len(&unread[end..]);
            let whole = length.filter(|length| end + length <= unread.len());
            if let Some(length) =
                whole.filter(|length| bgzf_block_in_bound(&unread[end..end + length]))
            {
                end += length;
                blocks += 1;
            } else if whole.is_some()
                || at_end
                || length.is_none() && unread.len() - end >= BGZF_HEADER
            {
                break;
            } else {
                let before = unread.len();
                (&mut source)
                    .take(READ_BYTES as u64)
                    .read_to_end(&mut unread)?;
                at_end = unread.len() == before;
            }
        }
        if end == 0 {
            let rest = BufReader::with_capacity(READ_BYTES, Cursor::new(unread).chain(source));
            let mut members = if first {
                Members::new(rest, pass_through)
            } else {
                Members::after_member(rest, pass_through)
            };
            stream(&mut members, emit)?;
            return Ok(Finish {
                garbage: members.ended_in_garbage(),
            });
        }
        let rest = unread.split_off(end);
        let run = std::mem::replace(&mut unread, rest);
        if !emit(Raw::Blocks(run)) {
            return Ok(Finish::default());
        }
        first = false;
    }
}

fn read(file: &File, reading: Reading, emit: &mut dyn FnMut(Raw) -> bool) -> io::Result<Finish> {
    let reading = match reading {
        Reading::Auto => {
            let mut magic = [0u8; 2];
            // A file that cannot be read at an offset is read as a stream.
            match read_range(file, 0, &mut magic) {
                Ok(2) if magic != [0x1f, 0x8b] => Reading::Plain { from: 0 },
                _ => Reading::Auto,
            }
        }
        other => other,
    };
    match reading {
        Reading::Plain { from } => {
            let size = file
                .metadata()
                .map_or(0, |meta| if meta.is_file() { meta.len() } else { 0 });
            if size == 0 {
                // No length to cut into ranges: a pipe, a device, or nothing at all.
                stream(file, emit)?;
                return Ok(Finish::default());
            }
            let mut offset = from;
            while offset < size {
                let length = (size - offset).min(PIECE_BYTES as u64) as usize;
                if !emit(Raw::Range(offset, length)) {
                    break;
                }
                offset += length as u64;
            }
            Ok(Finish::default())
        }
        Reading::Gunzip => gzip(file, false, emit),
        Reading::Auto => gzip(file, true, emit),
    }
}

/// Hands every piece of `file` to `work` on the worker threads and what it made to `each`
/// on the calling thread, in file order. What was read before a failure reaches `each`
/// before the failure is returned.
pub fn pieces<T: Send>(
    file: &File,
    reading: Reading,
    threads: usize,
    work: impl Fn(Vec<u8>) -> T + Sync,
    mut each: impl FnMut(T) -> io::Result<()>,
) -> io::Result<Finish> {
    let made = |raw: Raw| -> (T, Option<io::Error>) {
        match raw {
            Raw::Bytes(bytes) => (work(bytes), None),
            Raw::Blocks(blocks) => {
                let mut bytes = Vec::new();
                let failure =
                    inflate_blocks(&blocks, &mut Decompress::new(false), &mut bytes).err();
                (work(bytes), failure)
            }
            Raw::Range(offset, length) => {
                let mut bytes = vec![0u8; length];
                let (filled, failure) = match read_range(file, offset, &mut bytes) {
                    Ok(filled) => (filled, None),
                    Err(error) => (0, Some(error)),
                };
                bytes.truncate(filled);
                (work(bytes), failure)
            }
        }
    };
    ordered(
        threads,
        |emit| read(file, reading, emit),
        made,
        |(value, failure)| {
            each(value)?;
            failure.map_or(Ok(()), Err)
        },
    )?
}

/// Which bytes end a line.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Ends {
    /// `\n`.
    Newline,
    /// `\n` and `\r`, as universal newlines read a file whose empty lines are not counted.
    NewlineOrReturn,
}

impl Ends {
    fn first(self, bytes: &[u8]) -> Option<usize> {
        match self {
            Ends::Newline => memchr(b'\n', bytes),
            Ends::NewlineOrReturn => memchr2(b'\n', b'\r', bytes),
        }
    }

    fn last(self, bytes: &[u8]) -> Option<usize> {
        match self {
            Ends::Newline => memrchr(b'\n', bytes),
            Ends::NewlineOrReturn => memrchr2(b'\n', b'\r', bytes),
        }
    }
}

/// A piece cut where its lines end: the end of the line it starts inside, what was made
/// of its whole lines, and the start of the line it ends inside. A piece that ends no
/// line is all `head`.
struct Cut<P> {
    head: Vec<u8>,
    part: Option<P>,
    tail: Vec<u8>,
}

fn cut<P>(mut bytes: Vec<u8>, ends: Ends, part: &impl Fn(&[u8]) -> P) -> Cut<P> {
    let Some(first) = ends.first(&bytes) else {
        return Cut {
            head: bytes,
            part: None,
            tail: Vec::new(),
        };
    };
    let last = ends.last(&bytes).expect("a line ends in the piece");
    let made = part(&bytes[first + 1..=last]);
    let tail = bytes[last + 1..].to_vec();
    bytes.truncate(first + 1);
    Cut {
        head: bytes,
        part: Some(made),
        tail,
    }
}

/// Hands the lines of `file` to `part` in runs of whole lines, each run ending where a
/// line ends but the file's last, and what `part` made to `each` in file order. The lines
/// read before a failure, the last one cut short among them, reach `each` before the
/// failure is returned.
pub fn lines<P: Send>(
    file: &File,
    reading: Reading,
    threads: usize,
    ends: Ends,
    part: impl Fn(&[u8]) -> P + Sync,
    mut each: impl FnMut(P) -> io::Result<()>,
) -> io::Result<Finish> {
    let mut carry: Vec<u8> = Vec::new();
    let finish = pieces(
        file,
        reading,
        threads,
        |bytes| cut(bytes, ends, &part),
        |piece: Cut<P>| {
            carry.extend_from_slice(&piece.head);
            if let Some(made) = piece.part {
                each(part(&carry))?;
                each(made)?;
                carry = piece.tail;
            }
            Ok(())
        },
    );
    if !carry.is_empty() {
        each(part(&carry))?;
    }
    finish
}

#[cfg(test)]
#[path = "pipeline_tests.rs"]
mod tests;
