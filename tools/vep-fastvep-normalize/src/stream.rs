// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The run from input to output: one reader cuts the input into runs of whole lines, the
//! rewriting threads take the runs as they come, and one writer puts the rewritten runs out
//! in input order.
//!
//! A rewriting thread is started when a run is ready and none waits for one, up to the
//! thread bound, so a run whose reading or writing is the slower part keeps few threads.
//! What is held at once follows the threads started and not the input: a run is about
//! [`RUN_BYTES`] long, two runs a rewriting thread and two more are in work at most, and
//! only a line longer than a run makes a run longer.

use crate::rewrite::{rewrite, Counts};
use memchr::{memrchr, memrchr2};
use std::collections::BTreeMap;
use std::io::{self, Read, Write};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{mpsc, Condvar, Mutex, MutexGuard, PoisonError};

/// The bytes read from the input at a time.
pub const RUN_BYTES: usize = 1 << 20;

/// Why a run stopped.
#[derive(Debug)]
pub enum Failure {
    /// The input could not be read on. The output holds the rewrite of the runs read before.
    Read(io::Error),
    /// The output took no more. It holds the rewrite of the runs written before.
    Write(io::Error),
    /// The input is not UTF-8, first at this offset. The output holds the rewrite of every
    /// line before the line of that byte.
    NotUtf8(u64),
}

/// The bytes of one run. What the buffer has held stays initialised, so a read into it
/// writes over old bytes and nothing is cleared between runs.
#[derive(Default)]
struct Run {
    bytes: Vec<u8>,
    len: usize,
}

impl Run {
    fn filled(&self) -> &[u8] {
        &self.bytes[..self.len]
    }

    /// The `want` bytes after the filled ones.
    fn room(&mut self, want: usize) -> &mut [u8] {
        let end = self.len + want;
        if self.bytes.len() < end {
            self.bytes.resize(end, 0);
        }
        &mut self.bytes[self.len..end]
    }

    fn refill(&mut self, bytes: &[u8]) {
        self.len = 0;
        self.room(bytes.len()).copy_from_slice(bytes);
        self.len = bytes.len();
    }
}

/// Reads until `buffer` is full or the input ends, and returns the bytes read.
fn read_whole(input: &mut impl Read, buffer: &mut [u8]) -> io::Result<usize> {
    let mut read = 0;
    while read < buffer.len() {
        match input.read(&mut buffer[read..]) {
            Ok(0) => break,
            Ok(n) => read += n,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(read)
}

/// The offset after the last line ending of `bytes` that is known to be whole: after the
/// last `\n`, or, in bytes without one, after the last `\r` that has a byte after it, which
/// shows the `\r` to end a line by itself.
fn end_of_last_line(bytes: &[u8]) -> Option<usize> {
    if let Some(at) = memrchr(b'\n', bytes) {
        return Some(at + 1);
    }
    let (_, before_last) = bytes.split_last()?;
    memrchr(b'\r', before_last).map(|at| at + 1)
}

/// Cuts the input into runs of whole lines.
struct Reader<R> {
    input: R,
    run_bytes: usize,
    /// What was read past the end of the last run: the start of the next.
    carry: Vec<u8>,
    /// The offset in the input of the next run.
    offset: u64,
    ended: bool,
}

impl<R: Read> Reader<R> {
    fn new(input: R, run_bytes: usize) -> Reader<R> {
        Reader {
            input,
            run_bytes: run_bytes.max(1),
            carry: Vec::new(),
            offset: 0,
            ended: false,
        }
    }

    /// Fills `run` with the next run and returns the offset of its first byte in the input,
    /// or `None` when the input holds no more. A run ends after a line ending, the last run
    /// where the input ends.
    fn next_run(&mut self, run: &mut Run) -> io::Result<Option<u64>> {
        run.refill(&self.carry);
        self.carry.clear();
        let mut want = self.run_bytes;
        while !self.ended {
            let read = read_whole(&mut self.input, run.room(want))?;
            run.len += read;
            if read < want {
                self.ended = true;
            } else if let Some(cut) = end_of_last_line(run.filled()) {
                self.carry.extend_from_slice(&run.filled()[cut..]);
                run.len = cut;
                break;
            } else {
                // One line fills the run: hold twice as much and look again.
                want = run.len;
            }
        }
        if run.len == 0 {
            return Ok(None);
        }
        let offset = self.offset;
        self.offset += run.len as u64;
        Ok(Some(offset))
    }
}

/// Rewrites one run into `out`. A run that is not UTF-8 is rewritten up to the line of its
/// first byte outside UTF-8 and refused with the offset of that byte in the input; a line
/// ending is one byte of ASCII, so the runs are UTF-8 exactly when the input is.
fn rewrite_run(run: &[u8], offset: u64, out: &mut Vec<u8>) -> Result<Counts, u64> {
    match std::str::from_utf8(run) {
        Ok(_) => Ok(rewrite(run, out)),
        Err(error) => {
            let valid = error.valid_up_to();
            let line = memrchr2(b'\n', b'\r', &run[..valid]).map_or(0, |at| at + 1);
            rewrite(&run[..line], out);
            Err(offset + valid as u64)
        }
    }
}

/// The buffers of one run on its way through: the bytes read and the bytes to write.
#[derive(Default)]
struct Buffers {
    run: Run,
    out: Vec<u8>,
}

/// The buffer sets in work at once. The reader takes a set for each run and waits while
/// `limit` sets are out; the writer gives a set back once its run is written. The set that
/// came back last goes out first, so the sets in memory are the ones in use. The set of the
/// first run, which is read before the pool is made, counts as out from the start.
struct Pool {
    state: Mutex<PoolState>,
    given_back: Condvar,
}

#[derive(Default)]
struct PoolState {
    idle: Vec<Buffers>,
    out: usize,
    limit: usize,
    closed: bool,
}

impl Pool {
    fn new(limit: usize) -> Pool {
        Pool {
            state: Mutex::new(PoolState {
                out: 1,
                limit,
                ..PoolState::default()
            }),
            given_back: Condvar::new(),
        }
    }

    /// The sets a number of rewriting threads keep in work: one in each thread, one waiting
    /// for each, one being read and one being written.
    fn limit_for(rewriting: usize) -> usize {
        2 * rewriting + 2
    }

    fn raise_limit(&self, limit: usize) {
        let mut state = self.state();
        state.limit = state.limit.max(limit);
    }

    fn state(&self) -> MutexGuard<'_, PoolState> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// A set to fill, or `None` once the pool is closed.
    fn take(&self) -> Option<Buffers> {
        let mut state = self.state();
        loop {
            if state.closed {
                return None;
            }
            if let Some(buffers) = state.idle.pop() {
                state.out += 1;
                return Some(buffers);
            }
            if state.out < state.limit {
                state.out += 1;
                return Some(Buffers::default());
            }
            state = self
                .given_back
                .wait(state)
                .unwrap_or_else(PoisonError::into_inner);
        }
    }

    fn give(&self, buffers: Buffers) {
        let mut state = self.state();
        state.idle.push(buffers);
        state.out -= 1;
        drop(state);
        self.given_back.notify_one();
    }

    /// Ends every wait on the pool: no further run is read.
    fn close(&self) {
        self.state().closed = true;
        self.given_back.notify_all();
    }
}

/// Stops the run when the thread that holds it panics, so that no thread waits on that one
/// and the panic reaches the caller.
struct StopOnPanic<'a> {
    pool: &'a Pool,
    stopped: &'a AtomicBool,
}

impl Drop for StopOnPanic<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.stopped.store(true, Ordering::Relaxed);
            self.pool.close();
        }
    }
}

struct Job {
    index: u64,
    offset: u64,
    buffers: Buffers,
}

struct Done {
    index: u64,
    buffers: Buffers,
    result: Result<Counts, u64>,
}

/// Rewrites `input` into `output` and counts its rows.
///
/// `threads` bounds the threads that rewrite; one more reads and one more writes. With a
/// bound of one, or an input of one run, the whole run stays on the calling thread. Every
/// bound and every `run_bytes` write the same bytes.
pub fn normalise<R: Read, W: Write + Send>(
    input: R,
    output: W,
    threads: usize,
    run_bytes: usize,
) -> Result<Counts, Failure> {
    let mut reader = Reader::new(input, run_bytes);
    let mut buffers = Buffers::default();
    let Some(offset) = reader.next_run(&mut buffers.run).map_err(Failure::Read)? else {
        return Ok(Counts::default());
    };
    if threads <= 1 || reader.ended {
        in_turn(reader, output, buffers, offset)
    } else {
        across_threads(reader, output, threads, buffers, offset)
    }
}

/// Reads, rewrites and writes one run after the other on the calling thread.
fn in_turn<R: Read, W: Write>(
    mut reader: Reader<R>,
    mut output: W,
    mut buffers: Buffers,
    mut offset: u64,
) -> Result<Counts, Failure> {
    let mut counts = Counts::default();
    loop {
        let rewritten = rewrite_run(buffers.run.filled(), offset, &mut buffers.out);
        output.write_all(&buffers.out).map_err(Failure::Write)?;
        counts += rewritten.map_err(Failure::NotUtf8)?;
        match reader.next_run(&mut buffers.run).map_err(Failure::Read)? {
            Some(next) => offset = next,
            None => return Ok(counts),
        }
    }
}

/// Reads on the calling thread, rewrites on up to `threads` threads and writes on one more,
/// the runs in input order.
fn across_threads<R: Read, W: Write + Send>(
    mut reader: Reader<R>,
    mut output: W,
    threads: usize,
    first: Buffers,
    offset: u64,
) -> Result<Counts, Failure> {
    let pool = &Pool::new(Pool::limit_for(1));
    let stopped = &AtomicBool::new(false);
    // The rewriting threads that wait for a run.
    let waiting_for_a_run = &AtomicUsize::new(0);
    let (jobs, taken) = mpsc::channel::<Job>();
    let taken = &Mutex::new(taken);
    let (finished, done) = mpsc::channel::<Done>();
    std::thread::scope(|scope| {
        let writer = scope.spawn(move || {
            let _guard = StopOnPanic { pool, stopped };
            let mut waiting = BTreeMap::new();
            let mut next = 0u64;
            let mut counts = Counts::default();
            let mut failure = None;
            for Done {
                index,
                buffers,
                result,
            } in done
            {
                waiting.insert(index, (buffers, result));
                while let Some((buffers, result)) = waiting.remove(&next) {
                    next += 1;
                    if failure.is_none() {
                        failure = match (output.write_all(&buffers.out), result) {
                            (Err(error), _) => Some(Failure::Write(error)),
                            (Ok(()), Err(offset)) => Some(Failure::NotUtf8(offset)),
                            (Ok(()), Ok(run)) => {
                                counts += run;
                                None
                            }
                        };
                        if failure.is_some() {
                            stopped.store(true, Ordering::Relaxed);
                            pool.close();
                        }
                    }
                    pool.give(buffers);
                }
            }
            failure.map_or(Ok(counts), Err)
        });

        let rewrite_runs = |finished: mpsc::Sender<Done>| {
            let _guard = StopOnPanic { pool, stopped };
            loop {
                waiting_for_a_run.fetch_add(1, Ordering::Relaxed);
                let job = taken.lock().unwrap_or_else(PoisonError::into_inner).recv();
                waiting_for_a_run.fetch_sub(1, Ordering::Relaxed);
                let Ok(Job {
                    index,
                    offset,
                    mut buffers,
                }) = job
                else {
                    break;
                };
                if stopped.load(Ordering::Relaxed) {
                    continue;
                }
                let result = rewrite_run(buffers.run.filled(), offset, &mut buffers.out);
                let done = Done {
                    index,
                    buffers,
                    result,
                };
                if finished.send(done).is_err() {
                    break;
                }
            }
        };

        let mut read_failure = None;
        let mut rewriting = 0;
        let mut index = 0u64;
        let mut job = Some(Job {
            index,
            offset,
            buffers: first,
        });
        while let Some(ready) = job.take() {
            // A run is ready and no thread waits for one: one more rewrites, up to the bound.
            if rewriting < threads && waiting_for_a_run.load(Ordering::Relaxed) == 0 {
                let finished = finished.clone();
                scope.spawn(move || rewrite_runs(finished));
                rewriting += 1;
                pool.raise_limit(Pool::limit_for(rewriting));
            }
            if jobs.send(ready).is_err() {
                break;
            }
            index += 1;
            let Some(mut buffers) = pool.take() else {
                break;
            };
            match reader.next_run(&mut buffers.run) {
                Ok(Some(offset)) => {
                    job = Some(Job {
                        index,
                        offset,
                        buffers,
                    });
                }
                Ok(None) => {}
                Err(error) => read_failure = Some(error),
            }
        }
        drop(jobs);
        drop(finished);

        let written = writer
            .join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic));
        // A run that stopped the writer lies before the point where the reader stopped.
        match (written, read_failure) {
            (Err(failure), _) => Err(failure),
            (Ok(_), Some(error)) => Err(Failure::Read(error)),
            (Ok(counts), None) => Ok(counts),
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Rows and `#` lines under every line ending, a row longer than its neighbours, an empty
    /// line and a last row without an ending.
    const MIXED: &[u8] = b"## fastVEP output\r\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n\
v1\t21:1\tG\t-\t-\tTranscript\tintergenic_variant\r\n\
v2\t21:2\tT\tENSG1\tENST1\tTranscript\tmissense_variant\tIMPACT=MODERATE;SYMBOL=GENE1;BIOTYPE=protein_coding\r\
v3\t21:3\tC\t-\t-\tTranscript\r\r\n\n\
v4\t21:4\tA\t-\t-\tTranscript\tintergenic_variant\tx\ty\n\
#late\rv5\t21:5\tG\t-\t-\tTranscript";

    const THREADS: [usize; 4] = [1, 2, 3, 8];

    fn whole(input: &[u8]) -> (Vec<u8>, Counts) {
        let mut out = Vec::new();
        let counts = rewrite(input, &mut out);
        (out, counts)
    }

    fn run(input: &[u8], threads: usize, run_bytes: usize) -> (Result<Counts, Failure>, Vec<u8>) {
        let mut out = Vec::new();
        let result = normalise(input, &mut out, threads, run_bytes);
        (result, out)
    }

    fn held_for_every_run_length(input: &[u8]) {
        let (expected, counts) = whole(input);
        for run_bytes in 1..=input.len() + 2 {
            for threads in THREADS {
                let (result, out) = run(input, threads, run_bytes);
                assert_eq!(
                    result.unwrap(),
                    counts,
                    "{threads} threads, runs of {run_bytes}"
                );
                assert_eq!(out, expected, "{threads} threads, runs of {run_bytes}");
            }
        }
    }

    #[test]
    fn every_run_length_and_thread_count_writes_the_bytes_of_one_run() {
        assert_eq!(
            whole(MIXED).1,
            Counts {
                rows: 7,
                rewritten: 4
            }
        );
        held_for_every_run_length(MIXED);
    }

    #[test]
    fn inputs_of_one_line_ending_or_none_are_cut_the_same() {
        for input in [
            &b"a\tb\tc\t-\t-\tTranscript\rd\re\tf\tg\t-\t-\tTranscript\th\r"[..],
            b"a\tb\tc\t-\t-\tTranscript\r\nd\r\n\r\ne\r\n",
            b"a\tb\tc\t-\t-\tTranscript\tone row and no line ending",
            b"#a\tb\tc\t-\t-\tTranscript\tone hash line and no line ending",
            b"\r",
            b"\n",
            b"\r\n",
            b"x",
        ] {
            held_for_every_run_length(input);
        }
    }

    #[test]
    fn an_empty_input_writes_nothing() {
        for threads in THREADS {
            let (result, out) = run(b"", threads, 4);
            assert_eq!(result.unwrap(), Counts::default());
            assert!(out.is_empty());
        }
    }

    #[test]
    fn a_long_input_goes_through_the_buffer_pool_in_order() {
        let mut input = Vec::new();
        for row in 0..20_000u32 {
            let cells = if row % 3 == 0 {
                "-\t-\tTranscript"
            } else {
                "ENSG1\tENST1\tTranscript"
            };
            let ending = ["\n", "\r\n", "\n", "\r"][(row % 4) as usize];
            input.extend_from_slice(
                format!("v{row}\t21:{row}\tG\t{cells}\tterm{}{ending}", row % 7).as_bytes(),
            );
        }
        let (expected, counts) = whole(&input);
        assert_eq!(
            counts,
            Counts {
                rows: 20_000,
                rewritten: 6_667
            }
        );
        for threads in THREADS {
            for run_bytes in [64, 1000, 4096, RUN_BYTES] {
                let (result, out) = run(&input, threads, run_bytes);
                assert_eq!(
                    result.unwrap(),
                    counts,
                    "{threads} threads, runs of {run_bytes}"
                );
                assert!(out == expected, "{threads} threads, runs of {run_bytes}");
            }
        }
    }

    #[test]
    fn runs_end_after_a_line_ending_and_never_inside_one() {
        for run_bytes in 1..=MIXED.len() + 2 {
            let mut reader = Reader::new(MIXED, run_bytes);
            let mut run = Run::default();
            let mut read = Vec::new();
            while let Some(offset) = reader.next_run(&mut run).unwrap() {
                assert_eq!(offset as usize, read.len());
                read.extend_from_slice(run.filled());
                if read.len() < MIXED.len() {
                    let last = *read.last().unwrap();
                    assert!(last == b'\n' || last == b'\r', "runs of {run_bytes}");
                    assert!(
                        !(last == b'\r' && MIXED[read.len()] == b'\n'),
                        "runs of {run_bytes}"
                    );
                }
            }
            assert_eq!(read, MIXED);
        }
    }

    #[test]
    fn a_byte_outside_utf_8_stops_the_run_at_its_offset_after_the_lines_before_its_line() {
        let mut cases = Vec::new();
        for at in [0, 1, 40, MIXED.len() / 2, MIXED.len() - 1, MIXED.len()] {
            let mut input = MIXED.to_vec();
            input.insert(at, 0xff);
            cases.push((input, at));
        }
        // Between the two bytes of one line ending, a sequence the input ends inside, and a
        // surrogate, which UTF-8 does not hold.
        let crlf = MIXED.iter().position(|byte| *byte == b'\r').unwrap() + 1;
        let mut input = MIXED.to_vec();
        input.insert(crlf, 0xff);
        cases.push((input, crlf));
        cases.push(([MIXED, b"\n\xe2\x86"].concat(), MIXED.len() + 1));
        cases.push(([b"a\tb\n\xed\xa0\x80\n", MIXED].concat(), 4));
        for (input, at) in cases {
            let line = input[..at]
                .iter()
                .rposition(|byte| *byte == b'\n' || *byte == b'\r')
                .map_or(0, |ending| ending + 1);
            let (before, _) = whole(&input[..line]);
            for run_bytes in [1, 2, 7, 64, 1000] {
                for threads in THREADS {
                    let (result, out) = run(&input, threads, run_bytes);
                    match result {
                        Err(Failure::NotUtf8(offset)) => assert_eq!(offset as usize, at),
                        other => {
                            panic!("{other:?} at {at}, {threads} threads, runs of {run_bytes}")
                        }
                    }
                    assert_eq!(
                        out, before,
                        "at {at}, {threads} threads, runs of {run_bytes}"
                    );
                }
            }
        }
    }

    /// Takes `room` bytes and refuses the write that would pass it.
    struct Full {
        room: usize,
        written: Vec<u8>,
    }

    impl Write for &mut Full {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if self.written.len() + bytes.len() > self.room {
                return Err(io::Error::other("no room"));
            }
            self.written.extend_from_slice(bytes);
            Ok(bytes.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn an_output_that_refuses_a_write_stops_the_run() {
        let (expected, _) = whole(MIXED);
        for run_bytes in [1, 7, 64, 1000] {
            for threads in THREADS {
                let mut full = Full {
                    room: expected.len() / 2,
                    written: Vec::new(),
                };
                match normalise(MIXED, &mut full, threads, run_bytes) {
                    Err(Failure::Write(error)) => assert_eq!(error.to_string(), "no room"),
                    other => panic!("{other:?}, {threads} threads, runs of {run_bytes}"),
                }
                assert!(expected.starts_with(&full.written));
            }
        }
    }

    /// Gives its bytes a few at a time, then fails.
    struct Broken<'a> {
        bytes: &'a [u8],
    }

    impl Read for Broken<'_> {
        fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
            if self.bytes.is_empty() {
                return Err(io::Error::other("input gone"));
            }
            let n = buffer.len().min(self.bytes.len()).min(5);
            buffer[..n].copy_from_slice(&self.bytes[..n]);
            self.bytes = &self.bytes[n..];
            Ok(n)
        }
    }

    #[test]
    fn an_input_that_fails_stops_the_run_after_the_runs_read_before() {
        let (expected, _) = whole(MIXED);
        for run_bytes in [1, 7, 64, 1000] {
            for threads in THREADS {
                let mut out = Vec::new();
                match normalise(Broken { bytes: MIXED }, &mut out, threads, run_bytes) {
                    Err(Failure::Read(error)) => assert_eq!(error.to_string(), "input gone"),
                    other => panic!("{other:?}, {threads} threads, runs of {run_bytes}"),
                }
                assert!(expected.starts_with(&out));
            }
        }
    }

    #[test]
    fn a_read_that_is_interrupted_is_taken_up_again() {
        struct Interrupted<'a> {
            bytes: &'a [u8],
            interrupt: bool,
        }
        impl Read for Interrupted<'_> {
            fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
                self.interrupt = !self.interrupt;
                if self.interrupt {
                    return Err(io::ErrorKind::Interrupted.into());
                }
                let n = buffer.len().min(self.bytes.len()).min(3);
                buffer[..n].copy_from_slice(&self.bytes[..n]);
                self.bytes = &self.bytes[n..];
                Ok(n)
            }
        }
        let (expected, counts) = whole(MIXED);
        let mut out = Vec::new();
        let input = Interrupted {
            bytes: MIXED,
            interrupt: false,
        };
        assert_eq!(normalise(input, &mut out, 2, 16).unwrap(), counts);
        assert_eq!(out, expected);
    }

    #[test]
    fn the_pool_lends_its_limit_and_no_set_after_it_closes() {
        let pool = Pool::new(2);
        let first = Buffers::default();
        let _second = pool.take().unwrap();
        std::thread::scope(|scope| {
            let third = scope.spawn(|| pool.take().is_some());
            pool.give(first);
            assert!(third.join().unwrap());
            pool.raise_limit(1);
            pool.raise_limit(3);
            assert!(pool.take().is_some());
            let fifth = scope.spawn(|| pool.take().is_some());
            pool.close();
            assert!(!fifth.join().unwrap());
        });
        assert!(pool.take().is_none());
        assert_eq!(Pool::limit_for(1), 4);
    }
}
