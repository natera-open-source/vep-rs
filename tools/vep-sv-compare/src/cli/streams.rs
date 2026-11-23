// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Stdout and stderr as Python's `print` writes them: what a stream holds back, when it writes, and how a run ends whose stream cannot be written.
//!
//! Stdout that is no terminal is held in two layers, as the interpreter's text writer and
//! its buffered writer hold it. The text of the writes is gathered until it reaches
//! [`TEXT_CHUNK_BYTES`], and the chunks are gathered in a block that is written when the
//! next chunk does not fit: the first byte reaches the file after about 128 KiB of output,
//! or when the run ends. Stdout on a terminal and stderr everywhere are written line by
//! line. Where one file takes both streams, a line of stderr therefore stands before the
//! stdout text that was held when the line was written.
//!
//! A write that fails ends nothing by itself. The bytes the file did not take stay held and
//! every later flush tries them again; a run whose last flush of a stream fails ends with
//! [`FLUSH_FAILED`]. A write past the file-size limit is such a write and no end of the
//! process: the binary ignores the limit's signal from its start.

use std::fs::File;
use std::io::{self, ErrorKind, IsTerminal, Write};
use std::os::fd::AsFd;
use std::os::unix::fs::MetadataExt;

use memchr::memchr2;

use crate::text::os_failure;

/// The bytes of text a stream gathers before it hands them to its block.
pub const TEXT_CHUNK_BYTES: usize = 8192;

/// The smallest block of a stream, and the block of a file that states no size of its own.
pub const MIN_BLOCK_BYTES: usize = 128 * 1024;

/// The largest block of a stream.
pub const MAX_BLOCK_BYTES: usize = 8192 * 1024;

/// The exit code of a run that ends with stdout or stderr holding bytes it cannot write.
pub const FLUSH_FAILED: u8 = 120;

/// The block of a stream on a file whose preferred size of a write is `blksize` bytes.
pub fn block_bytes(blksize: u64) -> usize {
    match usize::try_from(blksize) {
        Ok(0) => MIN_BLOCK_BYTES,
        Ok(blksize) => blksize.clamp(MIN_BLOCK_BYTES, MAX_BLOCK_BYTES),
        Err(_) => MAX_BLOCK_BYTES,
    }
}

/// The name of the exception Python raises for a failed write.
pub fn py_error_class(error: &io::Error) -> &'static str {
    match error.kind() {
        ErrorKind::BrokenPipe => "BrokenPipeError",
        ErrorKind::ConnectionReset => "ConnectionResetError",
        ErrorKind::ConnectionAborted => "ConnectionAbortedError",
        ErrorKind::ConnectionRefused => "ConnectionRefusedError",
        ErrorKind::WouldBlock => "BlockingIOError",
        ErrorKind::PermissionDenied => "PermissionError",
        ErrorKind::TimedOut => "TimeoutError",
        _ => "OSError",
    }
}

/// A stream that holds its text as a Python text stream over a buffered writer does.
pub struct Held<W: Write> {
    raw: W,
    /// The text of the writes that waits to be handed to the block.
    text: Vec<u8>,
    /// The chunks of text the file has not taken.
    block: Vec<u8>,
    chunk_bytes: usize,
    block_bytes: usize,
    /// Whether a write that holds a line end is written at once.
    by_line: bool,
}

impl<W: Write> Held<W> {
    /// A stream over `raw` with a block of `block_bytes`, written line by line when
    /// `by_line` is set.
    pub fn new(raw: W, block_bytes: usize, by_line: bool) -> Held<W> {
        Held::with_sizes(raw, TEXT_CHUNK_BYTES, block_bytes, by_line)
    }

    fn with_sizes(raw: W, chunk_bytes: usize, block_bytes: usize, by_line: bool) -> Held<W> {
        Held {
            raw,
            text: Vec::new(),
            block: Vec::new(),
            chunk_bytes,
            block_bytes,
            by_line,
        }
    }

    /// Takes one write of text, as one call of a Python text stream's `write`. The text is
    /// handed to the block once a chunk of it is held, a text of a chunk or more after the
    /// text held before it; a stream written line by line writes what it holds after a text
    /// with a line feed or a carriage return.
    pub fn write_text(&mut self, text: &[u8]) -> io::Result<()> {
        if text.len() >= self.chunk_bytes {
            self.hand_on()?;
        }
        self.text.extend_from_slice(text);
        let ends_line = self.by_line && memchr2(b'\n', b'\r', text).is_some();
        if ends_line || self.text.len() >= self.chunk_bytes {
            self.hand_on()?;
        }
        if ends_line {
            self.write_block()?;
        }
        Ok(())
    }

    /// Writes a line as `print` does: the text, then the line feed, each a write of its own.
    pub fn print(&mut self, line: &str) -> io::Result<()> {
        self.write_text(line.as_bytes())?;
        self.write_text(b"\n")
    }

    /// Hands the text held to the block as one chunk. A chunk the block has no room for
    /// waits until the block is written, and is lost when that fails; a chunk of the
    /// block's size or more is then written by itself while a block or more of it is
    /// left, and what is left of it is held.
    fn hand_on(&mut self) -> io::Result<()> {
        if self.text.is_empty() {
            return Ok(());
        }
        let chunk = std::mem::take(&mut self.text);
        let room = self.block_bytes.saturating_sub(self.block.len());
        if chunk.len() <= room && chunk.len() < self.block_bytes {
            self.block.extend_from_slice(&chunk);
            return Ok(());
        }
        self.write_block()?;
        self.block = chunk;
        self.write_while_holding(self.block_bytes.max(1))
    }

    /// Writes the block to the file. The bytes the file does not take stay held.
    fn write_block(&mut self) -> io::Result<()> {
        self.write_while_holding(1)
    }

    /// Writes the block to the file, one write after another, while `least` bytes or more
    /// of it are left. The bytes the file does not take stay held.
    fn write_while_holding(&mut self, least: usize) -> io::Result<()> {
        let mut written = 0;
        let result = loop {
            if self.block.len() - written < least {
                break Ok(());
            }
            match self.raw.write(&self.block[written..]) {
                Ok(0) => break Err(ErrorKind::WriteZero.into()),
                Ok(taken) => written += taken,
                Err(error) if error.kind() == ErrorKind::Interrupted => {}
                Err(error) => break Err(error),
            }
        };
        self.block.drain(..written);
        result
    }

    /// Writes everything held to the file.
    pub fn flush(&mut self) -> io::Result<()> {
        self.hand_on()?;
        self.write_block()
    }
}

/// Each call of `write` is one write of text and takes all of it.
impl<W: Write> Write for Held<W> {
    fn write(&mut self, text: &[u8]) -> io::Result<usize> {
        self.write_text(text).map(|()| text.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Held::flush(self)
    }
}

/// The two standard streams of a run.
pub struct Streams<O: Write, E: Write> {
    pub out: Held<O>,
    pub err: Held<E>,
}

impl<O: Write, E: Write> Streams<O, E> {
    /// Ends the run as the interpreter ends one that exits with `code`, and gives the exit
    /// code of the process. Stdout is flushed, and a failure of that is written to stderr
    /// in the interpreter's words; then stderr is flushed. Either failure makes the exit
    /// code [`FLUSH_FAILED`], whatever `code` is.
    pub fn close(&mut self, code: u8) -> u8 {
        let mut code = code;
        if let Err(error) = self.out.flush() {
            code = FLUSH_FAILED;
            let words = format!(
                "Exception ignored while flushing sys.stdout:\n{}: {}",
                py_error_class(&error),
                os_failure(&error, None)
            );
            // A stderr that cannot take the words fails its own flush below.
            let _ = self.err.print(&words);
        }
        if self.err.flush().is_err() {
            code = FLUSH_FAILED;
        }
        code
    }
}

/// The file behind a standard stream and its preferred size of a write: a second
/// descriptor of the same open file, so that a write reaches the file when it is made.
/// A stream whose descriptor cannot be duplicated is written through its own handle.
fn file_of<S: AsFd + Write + 'static>(stream: S) -> (Box<dyn Write>, u64) {
    match stream.as_fd().try_clone_to_owned() {
        Ok(descriptor) => {
            let file = File::from(descriptor);
            let blksize = file.metadata().map_or(0, |metadata| metadata.blksize());
            (Box::new(file), blksize)
        }
        Err(_) => (Box::new(stream), 0),
    }
}

impl Streams<Box<dyn Write>, Box<dyn Write>> {
    /// The streams of the process: stdout held in blocks unless it is a terminal, stderr
    /// written line by line.
    pub fn of_process() -> Streams<Box<dyn Write>, Box<dyn Write>> {
        let on_terminal = io::stdout().is_terminal();
        let (out, out_blksize) = file_of(io::stdout());
        let (err, err_blksize) = file_of(io::stderr());
        Streams {
            out: Held::new(out, block_bytes(out_blksize), on_terminal),
            err: Held::new(err, block_bytes(err_blksize), true),
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::cell::{Cell, RefCell};
    use std::rc::Rc;

    /// The system error number of a write to a pipe or socket whose reader has gone.
    pub(crate) const EPIPE: i32 = 32;

    /// The system error number of a write to a device with no space left.
    const ENOSPC: i32 = 28;

    /// The writes a file took, in order, each under the name of the stream that made it.
    type Record = Vec<(&'static str, Vec<u8>)>;

    /// The texts written to a stream, and the sizes of the writes it makes for them before
    /// and at its flush.
    type Case<'a> = (&'a [&'a [u8]], &'a [usize], &'a [usize]);

    /// A file that records every write it takes, in order, under the name of its stream;
    /// while `refuses` holds a system error number it takes none and fails with that error.
    #[derive(Clone, Default)]
    pub(crate) struct Sink {
        pub(crate) writes: Rc<RefCell<Record>>,
        pub(crate) refuses: Rc<Cell<Option<i32>>>,
        name: &'static str,
    }

    impl Sink {
        /// This file as another stream writes it: one record of writes for both.
        pub(crate) fn named(&self, name: &'static str) -> Sink {
            Sink {
                name,
                ..self.clone()
            }
        }

        /// The sizes of the writes taken.
        fn sizes(&self) -> Vec<usize> {
            self.writes
                .borrow()
                .iter()
                .map(|(_, bytes)| bytes.len())
                .collect()
        }

        /// Everything written, in order.
        pub(crate) fn text(&self) -> String {
            let writes = self.writes.borrow();
            let bytes: Vec<u8> = writes
                .iter()
                .flat_map(|(_, bytes)| bytes.iter().copied())
                .collect();
            String::from_utf8(bytes).expect("UTF-8 text")
        }
    }

    impl Write for Sink {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if let Some(code) = self.refuses.get() {
                return Err(io::Error::from_raw_os_error(code));
            }
            self.writes.borrow_mut().push((self.name, bytes.to_vec()));
            Ok(bytes.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// The sizes of the writes a stream of these sizes makes for `texts`, before and at its
    /// flush.
    fn written(
        chunk_bytes: usize,
        block_bytes: usize,
        by_line: bool,
        texts: &[&[u8]],
    ) -> (Vec<usize>, Vec<usize>) {
        let sink = Sink::default();
        let mut held = Held::with_sizes(sink.clone(), chunk_bytes, block_bytes, by_line);
        for text in texts {
            held.write_text(text).unwrap();
        }
        let before = sink.sizes();
        held.flush().unwrap();
        let after = sink.sizes()[before.len()..].to_vec();
        (before, after)
    }

    #[test]
    fn a_stream_that_is_no_terminal_writes_in_blocks_of_whole_chunks() {
        // 4,000 lines as `print` writes them, the text of line i being (i * 37) % 200 bytes.
        let texts: Vec<Vec<u8>> = (0..4000).map(|i| vec![b'x'; (i * 37) % 200]).collect();
        let sink = Sink::default();
        let mut held = Held::new(sink.clone(), MIN_BLOCK_BYTES, false);
        for text in &texts {
            held.write_text(text).unwrap();
            held.write_text(b"\n").unwrap();
        }
        assert_eq!(sink.sizes(), [123_968, 124_102, 124_193]);
        held.flush().unwrap();
        assert_eq!(sink.sizes(), [123_968, 124_102, 124_193, 29_737]);
        let expected: Vec<u8> = texts
            .iter()
            .flat_map(|text| text.iter().copied().chain(*b"\n"))
            .collect();
        assert_eq!(sink.text().into_bytes(), expected);
    }

    #[test]
    fn text_is_handed_on_by_the_chunk_and_written_by_the_block() {
        let x = |count: usize| vec![b'x'; count];
        let cases: [Case; 10] = [
            (&[b"abc", b"\n", b"defg", b"\n"], &[], &[9]),
            (&[b"abcdefgh", b"ij"], &[], &[10]),
            (&[b"abc", &x(9), b"y"], &[], &[13]),
            (&[&x(8), &x(8), &x(8), &x(8), b"y"], &[], &[32, 1]),
            (&[&x(8), &x(8), &x(8), &x(8), &x(8)], &[32], &[8]),
            (&[b"ab", &x(32), b"c"], &[2, 32], &[1]),
            (&[b"ab", &x(40), b"c"], &[2, 40], &[1]),
            (&[&x(16), &x(16), b"z"], &[], &[32, 1]),
            (&[b"", b"", b"abc", b""], &[], &[3]),
            (&[&x(7), &x(7), &x(7), &x(7), &x(7), &x(7)], &[28], &[14]),
        ];
        for (texts, before, after) in cases {
            let sizes: Vec<usize> = texts.iter().map(|text| text.len()).collect();
            assert_eq!(
                written(8, 32, false, texts),
                (before.to_vec(), after.to_vec()),
                "{sizes:?}"
            );
        }
    }

    #[test]
    fn a_stream_written_by_line_writes_after_every_line_end() {
        let cases: [Case; 3] = [
            (&[b"ab", b"\n", b"cd\nef", b"\n", b"gh"], &[3, 5, 1], &[2]),
            (&[b"ab\rcd", b"ef", b"\n"], &[5, 3], &[]),
            (&[&[b'x'; 40], b"\n"], &[40, 1], &[]),
        ];
        for (texts, before, after) in cases {
            assert_eq!(
                written(8, 32, true, texts),
                (before.to_vec(), after.to_vec()),
                "{texts:?}"
            );
        }
    }

    #[test]
    fn print_is_the_text_and_then_the_line_feed() {
        let sink = Sink::default();
        let mut held = Held::with_sizes(sink.clone(), 4, 32, false);
        // The text reaches a chunk and is handed on without its line feed.
        held.print("abcd").unwrap();
        held.print("").unwrap();
        held.flush().unwrap();
        assert_eq!(sink.text(), "abcd\n\n");
        assert_eq!(sink.sizes(), [6]);
        let mut by_write = Held::with_sizes(sink.clone(), 4, 32, false);
        by_write.write_all(b"ab").unwrap();
        Write::flush(&mut by_write).unwrap();
        assert_eq!(sink.sizes(), [6, 2]);
    }

    #[test]
    fn bytes_a_file_does_not_take_stay_held_and_are_tried_again() {
        let sink = Sink::default();
        let mut held = Held::with_sizes(sink.clone(), 8, 32, true);
        sink.refuses.set(Some(EPIPE));
        let refused = held.print("first").unwrap_err();
        assert_eq!(refused.kind(), ErrorKind::BrokenPipe);
        assert_eq!(
            held.print("second").unwrap_err().raw_os_error(),
            Some(EPIPE)
        );
        assert_eq!(held.flush().unwrap_err().raw_os_error(), Some(EPIPE));
        assert!(sink.writes.borrow().is_empty());
        sink.refuses.set(None);
        held.flush().unwrap();
        assert_eq!(sink.text(), "first\nsecond\n");
        assert_eq!(sink.sizes(), [13]);
    }

    #[test]
    fn a_chunk_is_lost_where_the_block_before_it_cannot_be_written() {
        let sink = Sink::default();
        let mut held = Held::with_sizes(sink.clone(), 8, 16, false);
        held.write_text(b"01234567").unwrap();
        held.write_text(b"89abcdef").unwrap();
        assert!(sink.writes.borrow().is_empty());
        sink.refuses.set(Some(ENOSPC));
        assert_eq!(
            held.write_text(b"ghijklmn").unwrap_err().raw_os_error(),
            Some(ENOSPC)
        );
        sink.refuses.set(None);
        held.write_text(b"op").unwrap();
        held.flush().unwrap();
        assert_eq!(sink.text(), "0123456789abcdefop");
    }

    /// A file that takes at most three bytes a write, and fails once between two of them
    /// with an error that asks for the write again.
    struct Trickle {
        taken: Vec<u8>,
        interrupted: bool,
    }

    impl Write for Trickle {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if self.taken.len() == 3 && !self.interrupted {
                self.interrupted = true;
                return Err(ErrorKind::Interrupted.into());
            }
            let count = bytes.len().min(3);
            self.taken.extend_from_slice(&bytes[..count]);
            Ok(count)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn a_block_is_written_whole_over_short_and_interrupted_writes() {
        let trickle = Trickle {
            taken: Vec::new(),
            interrupted: false,
        };
        let mut held = Held::with_sizes(trickle, 8, 32, false);
        held.print("0123456789").unwrap();
        held.flush().unwrap();
        assert_eq!(held.raw.taken, b"0123456789\n");
        assert!(held.raw.interrupted);
    }

    #[test]
    fn a_text_of_the_block_s_size_is_written_while_a_block_of_it_is_left() {
        let trickle = Trickle {
            taken: Vec::new(),
            interrupted: false,
        };
        let mut held = Held::with_sizes(trickle, 8, 32, false);
        held.write_text(&[b'x'; 40]).unwrap();
        // Three writes of three bytes leave 31, which are held as any block is.
        assert_eq!(held.raw.taken.len(), 9);
        held.flush().unwrap();
        assert_eq!(held.raw.taken, [b'x'; 40]);
    }

    #[test]
    fn a_block_is_sized_by_the_file_between_128_kib_and_8_mib() {
        assert_eq!(block_bytes(0), 131_072);
        assert_eq!(block_bytes(512), 131_072);
        assert_eq!(block_bytes(4096), 131_072);
        assert_eq!(block_bytes(131_072), 131_072);
        assert_eq!(block_bytes(1 << 20), 1 << 20);
        assert_eq!(block_bytes(8192 * 1024), 8192 * 1024);
        assert_eq!(block_bytes(1 << 26), 8192 * 1024);
        assert_eq!(block_bytes(u64::MAX), 8192 * 1024);
    }

    #[test]
    fn a_failed_write_is_named_as_python_names_its_exception() {
        let class = |code: i32| py_error_class(&io::Error::from_raw_os_error(code));
        assert_eq!(class(EPIPE), "BrokenPipeError");
        assert_eq!(class(ENOSPC), "OSError");
        let of_kind = |kind: ErrorKind| py_error_class(&kind.into());
        assert_eq!(of_kind(ErrorKind::ConnectionReset), "ConnectionResetError");
        assert_eq!(of_kind(ErrorKind::WouldBlock), "BlockingIOError");
        assert_eq!(of_kind(ErrorKind::PermissionDenied), "PermissionError");
        assert_eq!(of_kind(ErrorKind::TimedOut), "TimeoutError");
        assert_eq!(of_kind(ErrorKind::WriteZero), "OSError");
        assert_eq!(of_kind(ErrorKind::Other), "OSError");
    }

    /// Streams over two files with one record of writes, both of the given kinds.
    fn two_streams(out_by_line: bool) -> (Streams<Sink, Sink>, Sink, Sink) {
        let out = Sink::default().named("out");
        let err = Sink {
            refuses: Rc::default(),
            ..out.named("err")
        };
        let streams = Streams {
            out: Held::new(out.clone(), MIN_BLOCK_BYTES, out_by_line),
            err: Held::new(err.clone(), MIN_BLOCK_BYTES, true),
        };
        (streams, out, err)
    }

    /// The streams that made each write of a record, in order.
    fn order(sink: &Sink) -> Vec<&'static str> {
        sink.writes.borrow().iter().map(|(name, _)| *name).collect()
    }

    #[test]
    fn a_run_whose_streams_take_everything_keeps_its_exit_code() {
        for code in [0, 1, 2] {
            let (mut streams, out, _) = two_streams(false);
            streams.out.print("a line").unwrap();
            streams.err.print("a refusal").unwrap();
            assert_eq!(streams.close(code), code);
            // Stdout was held until the run ended; stderr was written at once.
            assert_eq!(order(&out), ["err", "out"]);
            assert_eq!(out.text(), "a refusal\na line\n");
        }
        let (mut streams, out, _) = two_streams(true);
        streams.out.print("a line").unwrap();
        streams.err.print("a refusal").unwrap();
        assert_eq!(streams.close(1), 1);
        assert_eq!(out.text(), "a line\na refusal\n");
    }

    #[test]
    fn a_stdout_that_cannot_be_flushed_ends_the_run_with_120_and_two_lines_on_stderr() {
        for (code, class, words) in [
            (EPIPE, "BrokenPipeError", "Broken pipe"),
            (ENOSPC, "OSError", "No space left on device"),
        ] {
            for exit in [0, 1] {
                let (mut streams, out, err) = two_streams(false);
                streams.out.print("a line").unwrap();
                out.refuses.set(Some(code));
                assert_eq!(streams.close(exit), FLUSH_FAILED);
                assert_eq!(
                    err.text(),
                    format!(
                        "Exception ignored while flushing sys.stdout:\n{class}: [Errno {code}] {words}\n"
                    )
                );
            }
        }
        // A stdout that holds nothing has nothing to fail on.
        let (mut streams, out, err) = two_streams(false);
        out.refuses.set(Some(EPIPE));
        assert_eq!(streams.close(1), 1);
        assert_eq!(err.text(), "");
    }

    #[test]
    fn a_stderr_that_cannot_be_flushed_ends_the_run_with_120() {
        let (mut streams, out, err) = two_streams(false);
        streams.out.print("a line").unwrap();
        err.refuses.set(Some(EPIPE));
        assert!(streams.err.print("a refusal").is_err());
        assert_eq!(streams.close(1), FLUSH_FAILED);
        assert_eq!(out.text(), "a line\n");
        // Both gone: the words about stdout have nowhere to go.
        let (mut streams, out, err) = two_streams(false);
        streams.out.print("a line").unwrap();
        out.refuses.set(Some(EPIPE));
        err.refuses.set(Some(EPIPE));
        assert_eq!(streams.close(0), FLUSH_FAILED);
        assert!(out.writes.borrow().is_empty());
        // A stderr that was never written to fails nothing.
        let (mut streams, _, err) = two_streams(false);
        err.refuses.set(Some(EPIPE));
        assert_eq!(streams.close(0), 0);
    }
}
