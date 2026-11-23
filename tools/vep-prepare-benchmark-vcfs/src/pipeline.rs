// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One input from its file to its prepared file and digest, across threads.
//!
//! A reader thread cuts the input into pieces and queues them in input order. The pieces
//! of a BGZF file are runs of whole blocks, inflated on the worker threads; any other gzip
//! file is inflated member by member on the reader thread, and a plain file is read as it
//! is. The whole lines of each piece are processed on the worker threads. The calling
//! thread takes the pieces in order, joins the lines that cross from one piece to the next,
//! writes the prepared file and hands every chunk it has written to a thread that digests
//! it. Under a cap on the variant lines the lines are counted in order, so the calling
//! thread processes them itself and the input is read no further than the cap needs.

use crate::gzip::{bgzf_block_len, inflate_blocks, GzipMembers};
use crate::lines::{split, Chunk, Collector, Counts, Rules, Split};
use crate::paths::Style;
use flate2::Decompress;
use memchr::{memchr, memchr2};
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Cursor, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{channel, sync_channel, Receiver, Sender, SyncSender};
use std::sync::{Arc, Mutex};
use std::thread;

/// What a worker thread keeps from one job to the next.
pub struct Worker {
    inflater: Option<Decompress>,
}

type Job = Box<dyn FnOnce(&mut Worker) + Send + 'static>;

/// The worker threads every input of a run shares. A thread is started when a job comes
/// and no thread waits for one, up to the bound, so a run of small inputs starts none.
pub struct Pool {
    jobs: Option<Sender<Job>>,
    queue: Arc<Mutex<Receiver<Job>>>,
    waiting: Arc<AtomicUsize>,
    threads: Mutex<Vec<thread::JoinHandle<()>>>,
    most: usize,
}

impl Pool {
    pub fn new(most: usize) -> Pool {
        let (jobs, queue) = channel::<Job>();
        Pool {
            jobs: Some(jobs),
            queue: Arc::new(Mutex::new(queue)),
            waiting: Arc::new(AtomicUsize::new(0)),
            threads: Mutex::new(Vec::new()),
            most: most.max(1),
        }
    }

    fn start_worker(&self) -> thread::JoinHandle<()> {
        let queue = Arc::clone(&self.queue);
        let waiting = Arc::clone(&self.waiting);
        thread::spawn(move || {
            let mut worker = Worker { inflater: None };
            loop {
                waiting.fetch_add(1, Ordering::Relaxed);
                let job = queue.lock().expect("no job runs under the lock").recv();
                waiting.fetch_sub(1, Ordering::Relaxed);
                match job {
                    Ok(job) => job(&mut worker),
                    Err(_) => return,
                }
            }
        })
    }

    fn submit(&self, job: Job) {
        if self.waiting.load(Ordering::Relaxed) == 0 {
            let mut threads = self.threads.lock().expect("no panic under the lock");
            if threads.len() < self.most {
                threads.push(self.start_worker());
            }
        }
        self.jobs
            .as_ref()
            .expect("the pool takes jobs until it is dropped")
            .send(job)
            .expect("the queue lives as long as the pool");
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        self.jobs = None;
        let threads = std::mem::take(&mut *self.threads.lock().expect("no panic under the lock"));
        for thread in threads {
            let _ = thread.join();
        }
    }
}

/// A piece of one input, as the calling thread takes it.
enum Piece {
    /// Bytes whose lines are not processed.
    Raw(Vec<u8>),
    Split(Split),
    /// BGZF blocks that are not what the member reader would read: the input is read
    /// again from this offset, member by member.
    Reread(u64),
}

/// A place in the queue of one input's pieces: the piece, or where a worker thread will
/// leave it.
enum Slot {
    Ready(io::Result<Piece>),
    Pending(Receiver<Piece>),
}

/// Why a reader stops early.
enum Halt {
    /// The calling thread takes no more pieces.
    Closed,
    Failed(io::Error),
}

impl From<io::Error> for Halt {
    fn from(error: io::Error) -> Halt {
        Halt::Failed(error)
    }
}

/// The BGZF blocks of one piece; they inflate to about two mebibytes.
const BLOCKS_PER_PIECE: usize = 32;
/// The bytes of one piece of an input that is not BGZF.
const PIECE_BYTES: usize = 2 << 20;
/// The bytes of one piece under a cap on the variant lines.
const CAPPED_PIECE_BYTES: usize = 1 << 18;
/// How much of the file a reader asks for at once.
const READ_BYTES: usize = 1 << 20;
/// An input of up to this many bytes is prepared on the calling thread alone.
const SMALL_INPUT_BYTES: u64 = 1 << 20;
/// The bytes of one piece of such an input.
const SMALL_PIECE_BYTES: usize = 1 << 16;

struct Reader<'a> {
    queue: SyncSender<Slot>,
    pool: &'a Pool,
    rules: Arc<Rules>,
    /// Whether the lines are processed in order on the calling thread.
    in_order: bool,
}

impl Reader<'_> {
    fn send(&self, slot: Slot) -> Result<(), Halt> {
        self.queue.send(slot).map_err(|_| Halt::Closed)
    }

    /// Queues bytes of the input: as they are where the lines are processed in order,
    /// else through a worker thread that processes their whole lines.
    fn bytes(&self, bytes: Vec<u8>) -> Result<(), Halt> {
        if self.in_order {
            return self.send(Slot::Ready(Ok(Piece::Raw(bytes))));
        }
        let (done, slot) = channel();
        let rules = Arc::clone(&self.rules);
        self.pool.submit(Box::new(move |_| {
            let _ = done.send(processed(bytes, &rules));
        }));
        self.send(Slot::Pending(slot))
    }

    /// Queues a run of whole BGZF blocks that starts at `offset` of the file.
    fn blocks(&self, blocks: Vec<u8>, offset: u64) -> Result<(), Halt> {
        let (done, slot) = channel();
        let rules = Arc::clone(&self.rules);
        self.pool.submit(Box::new(move |worker| {
            let mut bytes = Vec::new();
            let inflater = worker
                .inflater
                .get_or_insert_with(|| Decompress::new(false));
            let piece = if inflate_blocks(&blocks, inflater, &mut bytes) {
                processed(bytes, &rules)
            } else {
                Piece::Reread(offset)
            };
            let _ = done.send(piece);
        }));
        self.send(Slot::Pending(slot))
    }

    /// Queues what `stream` reads, in pieces of `piece_bytes`, to its end; what it read
    /// before a failure is queued before the failure.
    fn stream(&self, mut stream: impl Read, piece_bytes: usize) -> Result<(), Halt> {
        loop {
            let mut bytes = Vec::with_capacity(piece_bytes);
            let read = (&mut stream)
                .take(piece_bytes as u64)
                .read_to_end(&mut bytes);
            if !bytes.is_empty() {
                self.bytes(bytes)?;
            }
            if read? == 0 {
                return Ok(());
            }
        }
    }

    /// Queues a gzip file: its leading BGZF blocks in runs, and from the first byte that
    /// does not start a whole BGZF block on, what the member reader reads.
    fn gzip(&self, mut source: File) -> Result<(), Halt> {
        if self.in_order {
            let members = GzipMembers::new(BufReader::with_capacity(READ_BYTES, source));
            return self.stream(members, CAPPED_PIECE_BYTES);
        }
        // The bytes of the file from `offset` on that are read and still to be queued.
        let mut unread: Vec<u8> = Vec::new();
        let mut offset = 0u64;
        let mut at_end = false;
        loop {
            let (mut end, mut blocks) = (0, 0);
            while blocks < BLOCKS_PER_PIECE {
                let length = bgzf_block_len(&unread[end..]);
                if length.is_some_and(|length| end + length <= unread.len()) {
                    end += length.expect("checked above");
                    blocks += 1;
                } else if at_end
                    || length.is_none() && unread.len() - end >= crate::gzip::BGZF_HEADER
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
                let members = if offset == 0 {
                    GzipMembers::new(rest)
                } else {
                    GzipMembers::after_member(rest)
                };
                return self.stream(members, PIECE_BYTES);
            }
            let rest = unread.split_off(end);
            let run = std::mem::replace(&mut unread, rest);
            self.blocks(run, offset)?;
            offset += end as u64;
        }
    }

    fn run(self, source: File, gzip: bool) {
        let result = if gzip {
            self.gzip(source)
        } else if self.in_order {
            self.stream(source, CAPPED_PIECE_BYTES)
        } else {
            self.stream(source, PIECE_BYTES)
        };
        if let Err(Halt::Failed(error)) = result {
            let _ = self.queue.send(Slot::Ready(Err(error)));
        }
    }
}

fn processed(bytes: Vec<u8>, rules: &Rules) -> Piece {
    match split(bytes, rules) {
        Ok(split) => Piece::Split(split),
        Err(bytes) => Piece::Raw(bytes),
    }
}

/// What one prepared file comes to.
pub struct Prepared {
    pub counts: Counts,
    pub sha256: String,
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn lost_thread() -> io::Error {
    io::Error::other("a thread of the run ended early")
}

/// Hands what `stream` reads to `collector` in pieces of `piece_bytes`, to its end; true
/// where the cap stopped the input. What was read before a failure is handed on before the
/// failure.
fn feed<E: FnMut(Chunk) -> io::Result<()>>(
    collector: &mut Collector<'_, E>,
    mut stream: impl Read,
    piece_bytes: usize,
) -> io::Result<bool> {
    loop {
        let mut bytes = Vec::with_capacity(piece_bytes);
        let read = (&mut stream)
            .take(piece_bytes as u64)
            .read_to_end(&mut bytes);
        if !bytes.is_empty() && collector.raw(bytes)? {
            return Ok(true);
        }
        if read? == 0 {
            return Ok(false);
        }
    }
}

/// How the inputs of a run are prepared: the rules of one file, the cap on its variant
/// lines, the worker threads and the pieces that may be read and still unwritten.
pub struct Preparation<'a> {
    pub rules: &'a Arc<Rules>,
    pub cap: Option<u64>,
    pub pool: &'a Pool,
    pub in_flight: usize,
}

impl Preparation<'_> {
    /// Prepares `src` into `dst`: opens the input, then creates the output, then reads.
    /// A small input is prepared on the calling thread alone, any other across threads.
    pub fn file(&self, src: &Path, dst: &Path, gzip: bool) -> io::Result<Prepared> {
        let source = File::open(src)?;
        let small = source
            .metadata()
            .is_ok_and(|meta| meta.len() <= SMALL_INPUT_BYTES);
        let target = File::create(dst)?;
        if small {
            self.on_one_thread(source, target, gzip)
        } else {
            self.across_threads(source, target, src, gzip)
        }
    }

    /// Reads, processes, writes and digests the input piece by piece on the calling thread.
    fn on_one_thread(&self, source: File, target: File, gzip: bool) -> io::Result<Prepared> {
        let mut writer = BufWriter::with_capacity(SMALL_PIECE_BYTES, target);
        let mut sha256 = Sha256::new();
        let emit = |chunk: Chunk| -> io::Result<()> {
            writer.write_all(chunk.bytes())?;
            sha256.update(chunk.bytes());
            Ok(())
        };
        let mut collector = Collector::new(self.rules, self.cap, emit);
        let stopped = if gzip {
            let mut members = GzipMembers::new(BufReader::with_capacity(SMALL_PIECE_BYTES, source));
            feed(&mut collector, &mut members, SMALL_PIECE_BYTES)?
        } else {
            let mut source = source;
            feed(&mut collector, &mut source, SMALL_PIECE_BYTES)?
        };
        if !stopped {
            collector.finish()?;
        }
        let counts = collector.counts;
        drop(collector);
        writer.flush()?;
        Ok(Prepared {
            counts,
            sha256: hex(&sha256.finalize()),
        })
    }

    /// Reads the input on one thread, processes its pieces on the worker threads, joins
    /// and writes them on the calling thread and digests on a fourth. `source` is `src`
    /// opened; `src` is opened again where BGZF blocks have to be read member by member.
    fn across_threads(
        &self,
        source: File,
        target: File,
        src: &Path,
        gzip: bool,
    ) -> io::Result<Prepared> {
        let (rules, cap) = (self.rules, self.cap);
        let mut writer = BufWriter::with_capacity(1 << 20, target);
        let (queue, slots) = sync_channel::<Slot>(self.in_flight.max(1));
        let (written, to_digest) = sync_channel::<Chunk>(8);
        thread::scope(|scope| {
            let reader = Reader {
                queue,
                pool: self.pool,
                rules: Arc::clone(rules),
                in_order: cap.is_some(),
            };
            scope.spawn(move || reader.run(source, gzip));
            let digest = scope.spawn(move || {
                let mut sha256 = Sha256::new();
                for chunk in to_digest {
                    sha256.update(chunk.bytes());
                }
                hex(&sha256.finalize())
            });
            // The queue and the digest's channel close when this closure returns, which
            // ends both threads on any path out of it.
            let counts = (move || -> io::Result<Counts> {
                let emit = |chunk: Chunk| -> io::Result<()> {
                    writer.write_all(chunk.bytes())?;
                    written.send(chunk).map_err(|_| lost_thread())
                };
                let mut collector = Collector::new(rules, cap, emit);
                let mut reread = None;
                let mut stopped = false;
                for slot in &slots {
                    let piece = match slot {
                        Slot::Ready(piece) => piece?,
                        Slot::Pending(piece) => piece.recv().map_err(|_| lost_thread())?,
                    };
                    match piece {
                        Piece::Raw(bytes) => stopped = collector.raw(bytes)?,
                        Piece::Split(split) => collector.split(split)?,
                        Piece::Reread(offset) => reread = Some(offset),
                    }
                    if stopped || reread.is_some() {
                        break;
                    }
                }
                drop(slots);
                if let Some(offset) = reread {
                    let mut source = File::open(src)?;
                    source.seek(SeekFrom::Start(offset))?;
                    let source = BufReader::with_capacity(READ_BYTES, source);
                    let mut members = if offset == 0 {
                        GzipMembers::new(source)
                    } else {
                        GzipMembers::after_member(source)
                    };
                    feed(&mut collector, &mut members, PIECE_BYTES)?;
                }
                if !stopped {
                    collector.finish()?;
                }
                let counts = collector.counts;
                drop(collector);
                writer.flush()?;
                Ok(counts)
            })()?;
            let sha256 = digest.join().map_err(|_| lost_thread())?;
            Ok(Prepared { counts, sha256 })
        })
    }
}

/// The contig naming of the first variant line of `src` whose first column is not empty;
/// unknown where it holds none.
pub fn first_variant_style(src: &Path, gzip: bool) -> io::Result<Style> {
    let source = File::open(src)?;
    let mut reader: Box<dyn Read> = if gzip {
        Box::new(GzipMembers::new(BufReader::with_capacity(1 << 16, source)))
    } else {
        Box::new(source)
    };
    let style_of = |line: &[u8]| {
        if line.first() == Some(&b'#') {
            return None;
        }
        let chrom = &line[..memchr(b'\t', line).unwrap_or(line.len())];
        (!chrom.is_empty()).then(|| Style::of(chrom))
    };
    let mut buf = vec![0u8; 1 << 16];
    // The text of the line being read, and whether a `\n` that follows belongs to the
    // `\r` that ended the line before it.
    let mut line: Vec<u8> = Vec::new();
    let mut after_cr = false;
    loop {
        let read = reader.read(&mut buf)?;
        if read == 0 {
            break;
        }
        let mut rest = &buf[..read];
        if after_cr && rest[0] == b'\n' {
            rest = &rest[1..];
        }
        after_cr = false;
        while !rest.is_empty() {
            let Some(end) = memchr2(b'\n', b'\r', rest) else {
                line.extend_from_slice(rest);
                break;
            };
            line.extend_from_slice(&rest[..end]);
            if let Some(style) = style_of(&line) {
                return Ok(style);
            }
            line.clear();
            let ended_by_cr = rest[end] == b'\r';
            rest = &rest[end + 1..];
            if ended_by_cr {
                match rest.first() {
                    Some(b'\n') => rest = &rest[1..],
                    Some(_) => {}
                    None => after_cr = true,
                }
            }
        }
    }
    Ok(style_of(&line).unwrap_or(Style::Unknown))
}

#[cfg(test)]
#[path = "pipeline_tests.rs"]
mod tests;
