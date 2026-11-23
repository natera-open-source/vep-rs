// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A BGZF file inflated on several threads. BGZF is gzip whose members are blocks of at most
//! 64 KiB that each state their own length, so the blocks of a file can be cut apart before
//! they are inflated. The bytes are the ones the member-by-member reader gives: a block that
//! is not of that shape, is cut short or does not inflate to what its trailer states is left
//! to that reader, from the block's own place on.

use super::{open_file, Compressed, Gzip, Member, ReadError, CHUNK};
use flate2::{Crc, Decompress, FlushDecompress, Status};
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};

/// The most threads that inflate the blocks of one file.
const MOST_THREADS: usize = 16;

/// The blocks handed to each thread ahead of the block being read.
const AHEAD: usize = 4;

/// The bytes of a block header: the gzip header with one extra field, `BC`, whose two bytes
/// are the length of the whole block less one.
const HEADER: usize = 18;

/// The bytes of a block trailer: the checksum and the length of what the block inflates to.
const TRAILER: usize = 8;

/// The most bytes a block inflates to.
const MOST_BYTES: usize = 64 * 1024;

/// What a thread makes of a block: its bytes, or nothing for a block that does not inflate
/// to what its trailer states.
type Inflated = Option<Vec<u8>>;

/// One thread that inflates blocks: the blocks go in and their bytes come out in order.
struct Inflater {
    blocks: SyncSender<Vec<u8>>,
    bytes: Receiver<Inflated>,
}

/// Inflates the deflate data of a block held to its trailer.
fn inflate_block(inflate: &mut Decompress, block: &[u8]) -> Inflated {
    let (data, trailer) = block.split_at(block.len() - TRAILER);
    let stored = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
    let length = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]) as usize;
    if length > MOST_BYTES {
        return None;
    }
    inflate.reset(false);
    let mut bytes = vec![0u8; length.max(1)];
    let status = inflate
        .decompress(data, &mut bytes, FlushDecompress::Finish)
        .ok()?;
    let whole = status == Status::StreamEnd
        && inflate.total_in() as usize == data.len()
        && inflate.total_out() as usize == length;
    if !whole {
        return None;
    }
    bytes.truncate(length);
    let mut crc = Crc::new();
    crc.update(&bytes);
    (crc.sum() == stored).then_some(bytes)
}

/// Starts a thread that inflates the blocks it is sent; `None` where the host starts none.
fn inflater() -> Option<Inflater> {
    let (blocks, taken) = sync_channel::<Vec<u8>>(AHEAD);
    let (given, bytes) = sync_channel::<Inflated>(AHEAD);
    let work = move || {
        let mut inflate = Decompress::new(false);
        for block in taken {
            if given.send(inflate_block(&mut inflate, &block)).is_err() {
                break;
            }
        }
    };
    let started = std::thread::Builder::new().spawn(work);
    started.ok().map(|_| Inflater { blocks, bytes })
}

/// What a read of the blocks gives.
pub(super) enum Step {
    /// Bytes were read into the buffer.
    Bytes(usize),
    /// The reading goes on member by member from this place of the file, behind a member
    /// when one was read before it.
    Members { offset: u64, after_member: bool },
}

/// A BGZF file whose blocks are inflated on several threads and read in the file's order.
pub(super) struct Bgzf {
    file: Option<BufReader<File>>,
    threads: usize,
    inflaters: Vec<Inflater>,
    /// The place of the block that is handed to a thread next.
    next: u64,
    /// The places of the blocks handed out and not read back, in order.
    ahead: VecDeque<u64>,
    handed: usize,
    taken: usize,
    /// Whether the bytes at `next` are a block to hand out.
    more: bool,
    /// The bytes of the block being read, and how many of them are read.
    block: Vec<u8>,
    read: usize,
}

impl Bgzf {
    /// The reader of `path` on at most `threads` threads; `None` for a file whose place
    /// cannot be set, which is read member by member alone.
    pub(super) fn open(path: &str, threads: usize) -> Result<Option<Bgzf>, ReadError> {
        let mut file = open_file(path)?;
        let regular = file.metadata().is_ok_and(|found| found.is_file());
        if !regular || file.stream_position().is_err() {
            return Ok(None);
        }
        Ok(Some(Bgzf {
            file: Some(BufReader::with_capacity(CHUNK, file)),
            threads: threads.min(MOST_THREADS),
            inflaters: Vec::new(),
            next: 0,
            ahead: VecDeque::new(),
            handed: 0,
            taken: 0,
            more: true,
            block: Vec::new(),
            read: 0,
        }))
    }

    /// The whole block at `next`, without its header; `None` for bytes that are no block of
    /// the shape, or that the file ends in or cannot be read.
    fn next_block(&mut self) -> Option<Vec<u8>> {
        let file = self.file.as_mut()?;
        let mut header = [0u8; HEADER];
        file.read_exact(&mut header).ok()?;
        let shaped =
            header[..4] == [0x1f, 0x8b, 8, 4] && header[10..16] == [6, 0, b'B', b'C', 2, 0];
        let length = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;
        if !shaped || length < HEADER + TRAILER {
            return None;
        }
        let mut block = vec![0u8; length - HEADER];
        file.read_exact(&mut block).ok()?;
        self.next += length as u64;
        Some(block)
    }

    /// Hands blocks to the threads until each has its share ahead or the file gives no
    /// further block.
    fn hand_out(&mut self) {
        while self.more && self.ahead.len() < self.threads * AHEAD {
            let place = self.next;
            let Some(block) = self.next_block() else {
                self.more = false;
                break;
            };
            if self.inflaters.len() < self.threads {
                match inflater() {
                    Some(started) => self.inflaters.push(started),
                    None => self.threads = self.inflaters.len(),
                }
            }
            let sent = match self.inflaters.len() {
                0 => false,
                count => self.inflaters[self.handed % count]
                    .blocks
                    .send(block)
                    .is_ok(),
            };
            if !sent {
                // The block is read member by member, from its own place.
                self.next = place;
                self.more = false;
                break;
            }
            self.ahead.push_back(place);
            self.handed += 1;
        }
    }

    /// Reads into `into`, which has room.
    pub(super) fn read(&mut self, into: &mut [u8]) -> Step {
        loop {
            if self.read < self.block.len() {
                let count = into.len().min(self.block.len() - self.read);
                into[..count].copy_from_slice(&self.block[self.read..self.read + count]);
                self.read += count;
                return Step::Bytes(count);
            }
            self.hand_out();
            let after_member = self.taken > 0;
            let Some(offset) = self.ahead.pop_front() else {
                return Step::Members {
                    offset: self.next,
                    after_member,
                };
            };
            let inflated = self.inflaters[self.taken % self.inflaters.len()]
                .bytes
                .recv();
            self.taken += 1;
            match inflated {
                Ok(Some(bytes)) => {
                    self.block = bytes;
                    self.read = 0;
                }
                _ => {
                    return Step::Members {
                        offset,
                        after_member,
                    }
                }
            }
        }
    }

    /// The member-by-member reader of the file from `offset` on; behind a member the zero
    /// bytes before the next one are skipped first.
    pub(super) fn members(&mut self, offset: u64, after_member: bool) -> Result<Gzip, ReadError> {
        self.inflaters.clear();
        let Some(reader) = self.file.take() else {
            return Err(ReadError::Tool(
                "the file is read to its failure".to_string(),
            ));
        };
        let mut file = reader.into_inner();
        file.seek(SeekFrom::Start(offset))
            .map_err(|error| super::read_error(&error))?;
        let mut gzip = Gzip {
            compressed: Compressed {
                file,
                bytes: vec![0; CHUNK],
                start: 0,
                end: 0,
            },
            member: Member::Header,
            inflate: Decompress::new(false),
            crc: Crc::new(),
        };
        if after_member {
            gzip.skip_padding()?;
        }
        Ok(gzip)
    }
}

#[cfg(test)]
#[path = "bgzf_tests.rs"]
mod tests;
