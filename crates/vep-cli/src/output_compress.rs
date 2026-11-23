// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Compressed output for `--compress_output`.
//!
//! Ensembl VEP pipes the output file, standard output included, through the
//! external `gzip -c` or `bgzip -c` (`BaseRunner.pm` `get_output_file_handle`).
//! vep-rs encodes in process instead: the decompressed bytes are the same
//! output, the compressed bytes are this encoder's, and the output file keeps the
//! name it was given. The BGZF encoder deflates its blocks on worker threads
//! (`--compression_threads`); the gzip encoder, one member, runs on the thread
//! that writes.

use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::sync::{Arc, Mutex};

use flate2::write::GzEncoder;
use noodles::bgzf;

/// The encoders `--compress_output` names (`Config.pm` `%VALID`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputCompression {
    /// One gzip member, as `gzip -c` writes.
    Gzip,
    /// A BGZF stream, as `bgzip -c` writes: gzip members of at most 64 KiB, each
    /// carrying the `BC` extra subfield, closed by the empty end-of-file member,
    /// which `bgzip`, `tabix` and `samtools` read.
    Bgzip,
}

impl OutputCompression {
    /// The encoder a `--compress_output` value names; `None` for a value outside
    /// the flag's enumeration.
    pub fn parse(value: &str) -> Option<Self> {
        match value {
            "gzip" => Some(OutputCompression::Gzip),
            "bgzip" => Some(OutputCompression::Bgzip),
            _ => None,
        }
    }

    /// Opens the encoder over `inner`: the writer the output goes through, and
    /// the handle whose [`CompressedOutput::finish`] writes the stream's trailer
    /// once every line is written. A handle dropped unfinished still closes the
    /// stream, without reporting a failure.
    ///
    /// `workers` is the BGZF encoder's deflate thread count: above one, the
    /// blocks are deflated on that many worker threads and written in order by
    /// one more; at one (or zero), on the thread that writes. Both cut the text
    /// into blocks at the same fill and deflate each block with the same call,
    /// so every count writes the same bytes. The gzip encoder is one member and
    /// runs on the writing thread whatever the count.
    pub fn open(
        self,
        inner: Box<dyn Write + Send>,
        workers: usize,
    ) -> (Box<dyn Write + Send>, CompressedOutput) {
        let encoder = match (self, NonZeroUsize::new(workers)) {
            (OutputCompression::Gzip, _) => {
                Encoder::Gzip(GzEncoder::new(inner, flate2::Compression::default()))
            }
            (OutputCompression::Bgzip, Some(workers)) if workers.get() > 1 => {
                Encoder::BgzipWorkers(bgzf::io::MultithreadedWriter::with_worker_count(
                    workers, inner,
                ))
            }
            (OutputCompression::Bgzip, _) => Encoder::Bgzip(bgzf::io::Writer::new(inner)),
        };
        let shared = Arc::new(Mutex::new(Some(encoder)));
        (
            Box::new(CompressedWriter {
                shared: Arc::clone(&shared),
                failure: None,
            }),
            CompressedOutput(shared),
        )
    }
}

enum Encoder {
    Gzip(GzEncoder<Box<dyn Write + Send>>),
    Bgzip(bgzf::io::Writer<Box<dyn Write + Send>>),
    BgzipWorkers(bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>),
}

impl Write for Encoder {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Encoder::Gzip(e) => e.write(buf),
            Encoder::Bgzip(e) => e.write(buf),
            Encoder::BgzipWorkers(e) => e.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Encoder::Gzip(e) => e.flush(),
            Encoder::Bgzip(e) => e.flush(),
            Encoder::BgzipWorkers(e) => e.flush(),
        }
    }
}

impl Encoder {
    /// Writes the stream's trailer and hands back the writer beneath.
    fn finish(self) -> io::Result<Box<dyn Write + Send>> {
        match self {
            Encoder::Gzip(e) => e.finish(),
            Encoder::Bgzip(e) => e.finish(),
            Encoder::BgzipWorkers(mut e) => e.finish(),
        }
    }
}

type Shared = Arc<Mutex<Option<Encoder>>>;

fn lock(shared: &Shared) -> std::sync::MutexGuard<'_, Option<Encoder>> {
    shared
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// The writer the output pipeline holds. A write after the stream is finished
/// is refused; a write after a failed write or flush repeats that failure, so
/// a buffer above the writer that retries its bytes reports the output's own
/// error.
struct CompressedWriter {
    shared: Shared,
    /// The kind and text of the error that closed the encoder.
    failure: Option<(io::ErrorKind, String)>,
}

/// Drops the encoder after `error` and hands back the record of it.
fn close(encoder: &mut Option<Encoder>, error: &io::Error) -> (io::ErrorKind, String) {
    // The worker-thread encoder has joined its frame writer by the time it
    // reports an error and panics on a later flush or finish, so it is dropped
    // here and never called again.
    *encoder = None;
    (error.kind(), error.to_string())
}

impl CompressedWriter {
    /// The error a write meets once the encoder is gone.
    fn closed(&self) -> io::Error {
        match &self.failure {
            Some((kind, message)) => io::Error::new(*kind, message.clone()),
            None => io::Error::other("the compressed output is finished"),
        }
    }
}

impl Write for CompressedWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut guard = lock(&self.shared);
        let Some(encoder) = guard.as_mut() else {
            return Err(self.closed());
        };
        let written = encoder.write(buf);
        if let Err(error) = &written {
            self.failure = Some(close(&mut guard, error));
        }
        written
    }

    fn flush(&mut self) -> io::Result<()> {
        let mut guard = lock(&self.shared);
        let Some(encoder) = guard.as_mut() else {
            return Ok(());
        };
        let flushed = encoder.flush();
        if let Err(error) = &flushed {
            self.failure = Some(close(&mut guard, error));
        }
        flushed
    }
}

/// The handle that closes a compressed output.
pub struct CompressedOutput(Shared);

impl CompressedOutput {
    /// Writes the stream's trailer (the gzip member's end, or BGZF's last block
    /// and end-of-file member) and flushes the file beneath; a second call does
    /// nothing, as does a call after a failed write.
    pub fn finish(&self) -> io::Result<()> {
        let Some(encoder) = lock(&self.0).take() else {
            return Ok(());
        };
        encoder.finish()?.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    /// A writer whose bytes the test reads back after the pipeline is done with it.
    #[derive(Clone, Default)]
    struct Sink(Arc<Mutex<Vec<u8>>>);

    impl Write for Sink {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// The counts the tests open the BGZF encoder with: the single-threaded
    /// encoder (one, and zero, which reads as one) and worker pools of several
    /// sizes.
    const WORKER_COUNTS: [usize; 5] = [0, 1, 2, 3, 8];

    /// `text` through the encoder `kind` with `workers`, written in `chunk`-byte
    /// pieces and finished twice (the second call is a no-op).
    fn compress_in_chunks(
        kind: OutputCompression,
        workers: usize,
        text: &[u8],
        chunk: usize,
    ) -> Vec<u8> {
        let sink = Sink::default();
        let (mut writer, output) = kind.open(Box::new(sink.clone()), workers);
        for piece in text.chunks(chunk) {
            writer.write_all(piece).unwrap();
        }
        writer.flush().unwrap();
        drop(writer);
        output.finish().unwrap();
        output.finish().unwrap();
        let bytes = sink.0.lock().unwrap().clone();
        bytes
    }

    fn compress(kind: OutputCompression, workers: usize, text: &[u8]) -> Vec<u8> {
        compress_in_chunks(kind, workers, text, text.len().max(1))
    }

    /// Output-shaped text wide enough for several BGZF blocks: `lines` tab-separated
    /// rows of about 70 bytes, so 5,000 rows fill five blocks and part of a sixth.
    fn rows(lines: usize) -> Vec<u8> {
        (0..lines)
            .flat_map(|i| {
                format!(
                    "rs{i}\t21:{}\tA\tENSG00000{i:05}\tENST00000{i:05}\tTranscript\tmissense_variant\n",
                    10_000 + i
                )
                .into_bytes()
            })
            .collect()
    }

    /// The 28-byte empty BGZF member (`bgzip`'s end-of-file marker).
    const BGZF_EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// Every worker count writes a BGZF stream (gzip members carrying the `BC`
    /// extra subfield, RFC 1952 `FEXTRA` with subfield id `B` `C`) that a BGZF
    /// reader and a multi-member gzip reader both decompress to the text written,
    /// whether the text arrives in one write or in pieces smaller than a block.
    #[test]
    fn bgzip_stream_decompresses_to_the_plain_output_at_every_worker_count() {
        let text = rows(5_000);
        for workers in WORKER_COUNTS {
            for chunk in [text.len(), 4_096, 100] {
                let bgzip = compress_in_chunks(OutputCompression::Bgzip, workers, &text, chunk);
                assert_eq!(&bgzip[..4], &[0x1f, 0x8b, 0x08, 0x04], "{workers} workers");
                assert_eq!(
                    &bgzip[12..14],
                    b"BC",
                    "{workers} workers: the BGZF extra subfield"
                );
                let mut back = Vec::new();
                bgzf::io::Reader::new(&bgzip[..])
                    .read_to_end(&mut back)
                    .unwrap();
                assert_eq!(back, text, "{workers} workers, {chunk}-byte writes");
                let mut back = Vec::new();
                flate2::read::MultiGzDecoder::new(&bgzip[..])
                    .read_to_end(&mut back)
                    .unwrap();
                assert_eq!(
                    back, text,
                    "{workers} workers, {chunk}-byte writes: a BGZF stream is also a multi-member gzip stream"
                );
            }
        }
    }

    /// Every worker count ends the stream with the empty end-of-file member, once,
    /// after the last data block; an empty output is the end-of-file member alone.
    #[test]
    fn bgzip_stream_ends_with_the_eof_member() {
        let text = rows(5_000);
        for workers in WORKER_COUNTS {
            let bgzip = compress(OutputCompression::Bgzip, workers, &text);
            assert_eq!(
                &bgzip[bgzip.len() - 28..],
                &BGZF_EOF,
                "{workers} workers: the BGZF end-of-file member"
            );
            assert_ne!(
                &bgzip[bgzip.len() - 56..bgzip.len() - 28],
                &BGZF_EOF,
                "{workers} workers: one end-of-file member"
            );
            let empty = compress(OutputCompression::Bgzip, workers, b"");
            assert_eq!(empty, BGZF_EOF, "{workers} workers: an empty output");
        }
    }

    /// The uncompressed size (`ISIZE`) of each BGZF block of `stream`, walked
    /// header to header: a block opens with the gzip magic and the `FEXTRA`
    /// flag, carries `BC` at offset 12 and its own length minus one as `BSIZE`
    /// at offset 16, little-endian, and ends with its `ISIZE`; the walk must end
    /// exactly at the stream's end.
    fn bgzf_block_sizes(stream: &[u8]) -> Vec<usize> {
        let mut sizes = Vec::new();
        let mut at = 0;
        while at < stream.len() {
            let block = &stream[at..];
            assert_eq!(
                &block[..4],
                &[0x1f, 0x8b, 0x08, 0x04],
                "block header at {at}"
            );
            assert_eq!(&block[12..14], b"BC", "block header at {at}");
            let len = usize::from(u16::from_le_bytes([block[16], block[17]])) + 1;
            assert!(len <= block.len(), "the block at {at} runs past the stream");
            let trailer: [u8; 4] = block[len - 4..len].try_into().unwrap();
            sizes.push(u32::from_le_bytes(trailer) as usize);
            at += len;
        }
        sizes
    }

    /// The single-threaded encoder and every worker pool cut the text into blocks
    /// at the same fill and deflate each block with the same call, so the
    /// compressed bytes are identical whatever the count and however the text is
    /// split across writes.
    #[test]
    fn bgzip_stream_is_byte_identical_across_worker_counts() {
        let text = rows(5_000);
        let reference = compress(OutputCompression::Bgzip, 1, &text);
        let sizes = bgzf_block_sizes(&reference);
        let (eof, blocks) = sizes.split_last().unwrap();
        assert_eq!(*eof, 0, "the end-of-file member closes the stream");
        assert!(
            blocks.len() >= 6,
            "several blocks, {} written",
            blocks.len()
        );
        assert!(
            blocks[..blocks.len() - 1].iter().all(|&n| n == blocks[0]),
            "every block but the last is filled to the same size: {blocks:?}"
        );
        assert_eq!(
            blocks.iter().sum::<usize>(),
            text.len(),
            "the blocks carry the whole text"
        );
        for workers in WORKER_COUNTS {
            for chunk in [text.len(), 4_096, 100] {
                let bgzip = compress_in_chunks(OutputCompression::Bgzip, workers, &text, chunk);
                assert!(
                    bgzip == reference,
                    "{workers} workers, {chunk}-byte writes: the compressed bytes differ from the single-threaded encoder's"
                );
            }
        }
    }

    /// `gzip -c` writes one member with no extra field; the gzip encoder ignores
    /// the worker count and writes the same member at every value, decompressed
    /// whole by a single-member reader.
    #[test]
    fn gzip_member_is_unchanged() {
        let text = rows(5_000);
        let reference = compress_in_chunks(OutputCompression::Gzip, 1, &text, 4_096);
        assert_eq!(&reference[..2], &[0x1f, 0x8b]);
        assert_eq!(reference[3] & 0x04, 0, "gzip -c writes no FEXTRA field");
        let mut back = Vec::new();
        flate2::read::MultiGzDecoder::new(&reference[..])
            .read_to_end(&mut back)
            .unwrap();
        assert_eq!(back, text);
        // One member: the first member alone decodes to the whole text.
        let mut first_member = Vec::new();
        flate2::read::GzDecoder::new(&reference[..])
            .read_to_end(&mut first_member)
            .unwrap();
        assert_eq!(first_member, text);
        for workers in WORKER_COUNTS {
            let gzip = compress_in_chunks(OutputCompression::Gzip, workers, &text, 4_096);
            assert!(
                gzip == reference,
                "{workers} workers: the gzip member differs"
            );
        }
        assert_ne!(reference, compress(OutputCompression::Bgzip, 1, &text));
    }

    /// The encoders the write-pattern tests open: gzip, and BGZF on the writing
    /// thread and on a worker pool.
    const ENCODERS: [(OutputCompression, usize); 3] = [
        (OutputCompression::Gzip, 1),
        (OutputCompression::Bgzip, 1),
        (OutputCompression::Bgzip, 4),
    ];

    /// `text` through the encoder as the output sink writes a run: in pieces of
    /// the sizes `sizes` cycles through, beneath a buffer of `buffer` bytes when
    /// one is given, flushed once every byte is written, then finished.
    fn compress_as_the_sink_writes(
        kind: OutputCompression,
        workers: usize,
        text: &[u8],
        sizes: &[usize],
        buffer: Option<usize>,
    ) -> Vec<u8> {
        let sink = Sink::default();
        let (writer, output) = kind.open(Box::new(sink.clone()), workers);
        let mut writer: Box<dyn Write> = match buffer {
            Some(bytes) => Box::new(io::BufWriter::with_capacity(bytes, writer)),
            None => writer,
        };
        let mut rest = text;
        for &size in sizes.iter().cycle() {
            if rest.is_empty() {
                break;
            }
            let (piece, tail) = rest.split_at(size.min(rest.len()));
            writer.write_all(piece).unwrap();
            rest = tail;
        }
        writer.flush().unwrap();
        drop(writer);
        output.finish().unwrap();
        let bytes = sink.0.lock().unwrap().clone();
        bytes
    }

    /// The compressed bytes follow the text alone. The output sink hands the
    /// encoder each rendered chunk whole, or flushes a buffer into it wherever a
    /// chunk overflows the buffer, so the sizes of the writes follow the task
    /// split and the buffer size; its one `flush` call comes when every byte is
    /// written. Each encoder writes the same stream under every such pattern:
    /// the gzip encoder's deflate stream takes no flush before its end, and the
    /// BGZF encoders cut the text at the block fill.
    #[test]
    fn a_compressed_stream_is_the_same_bytes_under_any_write_sizes_and_flushes() {
        let text = rows(5_000);
        let sizes: [&[usize]; 6] = [
            &[text.len()],
            &[65_536],
            &[4_096],
            &[100],
            &[7],
            &[1, 13, 70, 911, 4_097, 70_000],
        ];
        let buffers = [None, Some(1), Some(4_096), Some(4 << 20)];
        for (kind, workers) in ENCODERS {
            let reference = compress(kind, workers, &text);
            for pieces in sizes {
                for buffer in buffers {
                    let stream = compress_as_the_sink_writes(kind, workers, &text, pieces, buffer);
                    assert!(
                        stream == reference,
                        "{kind:?} with {workers} workers, writes of {pieces:?} bytes beneath a buffer of {buffer:?}: the compressed bytes differ from one write's"
                    );
                }
            }
        }
    }

    /// A `flush` call inside the stream closes a deflate block (gzip) or a BGZF
    /// block where it lands, so the compressed bytes follow such calls and the
    /// decompressed bytes do not. With a flush after every 4,096-byte write, and
    /// with one flush halfway, each encoder's stream differs from the unflushed
    /// one and still decompresses to the text: the gzip stream as one member, the
    /// BGZF stream as blocks closed by one end-of-file member.
    #[test]
    fn a_flush_inside_a_compressed_stream_leaves_the_decompressed_bytes_as_they_are() {
        let text = rows(5_000);
        for (kind, workers) in ENCODERS {
            let reference = compress(kind, workers, &text);
            for chunk in [4_096, text.len() / 2 + 1] {
                let sink = Sink::default();
                let (mut writer, output) = kind.open(Box::new(sink.clone()), workers);
                for piece in text.chunks(chunk) {
                    writer.write_all(piece).unwrap();
                    writer.flush().unwrap();
                }
                drop(writer);
                output.finish().unwrap();
                let stream = sink.0.lock().unwrap().clone();
                let pattern =
                    format!("{kind:?} with {workers} workers, a flush every {chunk} bytes");
                assert!(
                    stream != reference,
                    "{pattern}: a flush inside the stream closes a block there"
                );
                let mut back = Vec::new();
                flate2::read::MultiGzDecoder::new(&stream[..])
                    .read_to_end(&mut back)
                    .unwrap();
                assert!(back == text, "{pattern}: the decompressed bytes");
                match kind {
                    OutputCompression::Gzip => {
                        let mut member = Vec::new();
                        flate2::read::GzDecoder::new(&stream[..])
                            .read_to_end(&mut member)
                            .unwrap();
                        assert!(member == text, "{pattern}: one member");
                    }
                    OutputCompression::Bgzip => {
                        let sizes = bgzf_block_sizes(&stream);
                        assert_eq!(sizes.iter().sum::<usize>(), text.len(), "{pattern}");
                        assert_eq!(
                            sizes.iter().filter(|&&size| size == 0).count(),
                            1,
                            "{pattern}: one end-of-file member"
                        );
                        assert_eq!(&stream[stream.len() - 28..], &BGZF_EOF, "{pattern}");
                    }
                }
            }
        }
    }

    /// A finished output refuses further bytes, so a late write is an error and
    /// not a silent loss, at every worker count.
    #[test]
    fn writes_after_finish_fail() {
        for (kind, workers) in [
            (OutputCompression::Gzip, 1),
            (OutputCompression::Bgzip, 1),
            (OutputCompression::Bgzip, 4),
        ] {
            let sink = Sink::default();
            let (mut writer, output) = kind.open(Box::new(sink), workers);
            writer.write_all(b"x\n").unwrap();
            output.finish().unwrap();
            assert!(writer.write_all(b"y\n").is_err(), "{kind:?} {workers}");
            assert!(writer.flush().is_ok(), "{kind:?} {workers}");
        }
    }

    /// A writer that fails after `limit` bytes.
    struct Failing {
        written: usize,
        limit: usize,
    }

    impl Write for Failing {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.written + buf.len() > self.limit {
                return Err(io::Error::new(io::ErrorKind::StorageFull, "no space left"));
            }
            self.written += buf.len();
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// Once the output beneath has failed, the encoder reports the failure to the
    /// writer that meets it, and every later write repeats that failure, kind and
    /// text, while a flush and the finish report nothing further: the pipeline
    /// completes its output after a failed write, retrying the bytes its buffer
    /// holds, and the worker-thread encoder has shut down by then.
    #[test]
    fn a_failed_write_closes_the_encoder_and_later_writes_repeat_it() {
        let text = rows(20_000);
        for (kind, workers) in [
            (OutputCompression::Gzip, 1),
            (OutputCompression::Bgzip, 1),
            (OutputCompression::Bgzip, 4),
        ] {
            let failing = Failing {
                written: 0,
                limit: 1_000,
            };
            let (mut writer, output) = kind.open(Box::new(failing), workers);
            let failed = text
                .chunks(4_096)
                .find_map(|piece| writer.write_all(piece).err())
                .or_else(|| writer.flush().err())
                .unwrap_or_else(|| {
                    panic!("{kind:?} {workers}: the failing output was not reported")
                });
            assert_eq!(
                failed.kind(),
                io::ErrorKind::StorageFull,
                "{kind:?} {workers}"
            );
            for _ in 0..2 {
                let again = writer
                    .write_all(b"x")
                    .expect_err("a write after the failure");
                assert_eq!(again.kind(), failed.kind(), "{kind:?} {workers}");
                assert_eq!(again.to_string(), failed.to_string(), "{kind:?} {workers}");
            }
            assert!(writer.flush().is_ok(), "{kind:?} {workers}");
            assert!(output.finish().is_ok(), "{kind:?} {workers}");
        }
    }

    #[test]
    fn parse_names_the_two_encoders() {
        assert_eq!(
            OutputCompression::parse("gzip"),
            Some(OutputCompression::Gzip)
        );
        assert_eq!(
            OutputCompression::parse("bgzip"),
            Some(OutputCompression::Bgzip)
        );
        assert_eq!(OutputCompression::parse("zip"), None);
    }
}
