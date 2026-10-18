// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Compressed output for `--compress_output`.
//!
//! Ensembl VEP pipes the output file, standard output included, through the
//! external `gzip -c` or `bgzip -c` (`BaseRunner.pm` `get_output_file_handle`).
//! vep-rs encodes in process instead: the decompressed bytes are the same
//! output, the compressed bytes are this encoder's, and the output file keeps the
//! name it was given.

use std::io::{self, Write};
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
    pub fn open(self, inner: Box<dyn Write + Send>) -> (Box<dyn Write + Send>, CompressedOutput) {
        let encoder = match self {
            OutputCompression::Gzip => {
                Encoder::Gzip(GzEncoder::new(inner, flate2::Compression::default()))
            }
            OutputCompression::Bgzip => Encoder::Bgzip(bgzf::io::Writer::new(inner)),
        };
        let shared = Arc::new(Mutex::new(Some(encoder)));
        (
            Box::new(CompressedWriter(Arc::clone(&shared))),
            CompressedOutput(shared),
        )
    }
}

enum Encoder {
    Gzip(GzEncoder<Box<dyn Write + Send>>),
    Bgzip(bgzf::io::Writer<Box<dyn Write + Send>>),
}

impl Write for Encoder {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Encoder::Gzip(e) => e.write(buf),
            Encoder::Bgzip(e) => e.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Encoder::Gzip(e) => e.flush(),
            Encoder::Bgzip(e) => e.flush(),
        }
    }
}

type Shared = Arc<Mutex<Option<Encoder>>>;

fn lock(shared: &Shared) -> std::sync::MutexGuard<'_, Option<Encoder>> {
    shared
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// The writer the output pipeline holds; it fails once the stream is finished.
struct CompressedWriter(Shared);

impl Write for CompressedWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match lock(&self.0).as_mut() {
            Some(encoder) => encoder.write(buf),
            None => Err(io::Error::other("the compressed output is finished")),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match lock(&self.0).as_mut() {
            Some(encoder) => encoder.flush(),
            None => Ok(()),
        }
    }
}

/// The handle that closes a compressed output.
pub struct CompressedOutput(Shared);

impl CompressedOutput {
    /// Writes the stream's trailer (the gzip member's end, or BGZF's last block
    /// and end-of-file member) and flushes the file beneath; a second call does
    /// nothing.
    pub fn finish(&self) -> io::Result<()> {
        let Some(encoder) = lock(&self.0).take() else {
            return Ok(());
        };
        let mut inner = match encoder {
            Encoder::Gzip(e) => e.finish()?,
            Encoder::Bgzip(e) => e.finish()?,
        };
        inner.flush()
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

    fn compress(kind: OutputCompression, text: &[u8]) -> Vec<u8> {
        let sink = Sink::default();
        let (mut writer, output) = kind.open(Box::new(sink.clone()));
        writer.write_all(text).unwrap();
        writer.flush().unwrap();
        drop(writer);
        output.finish().unwrap();
        output.finish().unwrap();
        let bytes = sink.0.lock().unwrap().clone();
        bytes
    }

    /// `bgzip -c` writes gzip members carrying the `BC` extra subfield (RFC 1952
    /// `FEXTRA`, subfield id `B` `C`) and ends the file with the 28-byte empty
    /// member; `gzip -c` writes one member with no extra field. Both decompress to
    /// the bytes written, and the BGZF stream is read by a BGZF reader.
    #[test]
    fn bgzip_output_is_bgzf_and_gzip_output_is_plain_gzip() {
        let text: Vec<u8> = (0..3_000)
            .flat_map(|i| format!("21\t{}\t.\tA\tG\t.\t.\t.\n", 10_000 + i).into_bytes())
            .collect();

        let gzip = compress(OutputCompression::Gzip, &text);
        assert_eq!(&gzip[..2], &[0x1f, 0x8b]);
        assert_eq!(gzip[3] & 0x04, 0, "gzip -c writes no FEXTRA field");
        let mut back = Vec::new();
        flate2::read::MultiGzDecoder::new(&gzip[..])
            .read_to_end(&mut back)
            .unwrap();
        assert_eq!(back, text);
        // One member: the first member alone decodes to the whole text.
        let mut first_member = Vec::new();
        flate2::read::GzDecoder::new(&gzip[..])
            .read_to_end(&mut first_member)
            .unwrap();
        assert_eq!(first_member, text);

        let bgzip = compress(OutputCompression::Bgzip, &text);
        assert_eq!(&bgzip[..4], &[0x1f, 0x8b, 0x08, 0x04]);
        assert_eq!(&bgzip[12..14], b"BC", "the BGZF extra subfield");
        let eof: [u8; 28] = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        assert_eq!(
            &bgzip[bgzip.len() - 28..],
            &eof,
            "the BGZF end-of-file member"
        );
        let mut back = Vec::new();
        bgzf::io::Reader::new(&bgzip[..])
            .read_to_end(&mut back)
            .unwrap();
        assert_eq!(back, text);
        let mut back = Vec::new();
        flate2::read::MultiGzDecoder::new(&bgzip[..])
            .read_to_end(&mut back)
            .unwrap();
        assert_eq!(
            back, text,
            "a BGZF stream is also a multi-member gzip stream"
        );
        assert_ne!(gzip, bgzip);
    }

    /// A finished output refuses further bytes, so a late write is an error and
    /// not a silent loss.
    #[test]
    fn writes_after_finish_fail() {
        let sink = Sink::default();
        let (mut writer, output) = OutputCompression::Gzip.open(Box::new(sink));
        writer.write_all(b"x\n").unwrap();
        output.finish().unwrap();
        assert!(writer.write_all(b"y\n").is_err());
        assert!(writer.flush().is_ok());
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
