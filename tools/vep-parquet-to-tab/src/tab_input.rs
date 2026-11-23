// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--tab` output, plain or gzip, read as a stream of blocks of whole lines.
//!
//! The file is strict UTF-8. A line ends at `\n` alone, and every `\r` before it is
//! dropped, so a cell holding another separator character stays in its row. A `##` line is
//! dropped, the one `#` line is the header, an empty line is dropped, and every other line
//! is a data line.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use std::sync::mpsc::sync_channel;

use memchr::{memchr_iter, memrchr};
use rayon::prelude::*;

use crate::gunzip::{GzipReader, ReadFailure, GZIP_MAGIC};
use crate::pyfmt::{os_error_line, utf8_error_text};
use crate::Failure;

/// Bytes read at a time; a block ends at the last line end inside them.
const BLOCK_BYTES: usize = 16 << 20;

enum Source {
    Plain(File),
    Gzip(Box<GzipReader<BufReader<File>>>),
}

impl Source {
    fn read(&mut self, out: &mut [u8]) -> Result<usize, ReadFailure> {
        match self {
            Self::Plain(file) => Ok(file.read(out)?),
            Self::Gzip(reader) => reader.read(out),
        }
    }
}

/// A run of whole lines of the input, `offset` bytes into its decoded data.
struct Block {
    index: usize,
    offset: u64,
    bytes: Vec<u8>,
}

/// What one block holds beside its data lines.
struct Scan {
    block: Block,
    /// The first sequence that is not UTF-8, as the codec reports it.
    decode_error: Option<String>,
    /// The first two header lines of the block, without the `#`.
    headers: Vec<Vec<u8>>,
}

/// A `--tab` output read whole: its header and its blocks of lines.
pub struct TabInput {
    /// The header's cells.
    pub columns: Vec<Vec<u8>>,
    blocks: Vec<Vec<u8>>,
}

/// The lines of a block, each without its line end and without the `\r`s before it.
fn lines(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut start = 0;
    let mut ends = memchr_iter(b'\n', bytes);
    let mut done = false;
    std::iter::from_fn(move || {
        if done {
            return None;
        }
        let end = ends.next().unwrap_or_else(|| {
            done = true;
            bytes.len()
        });
        let mut line = &bytes[start..end];
        start = end + 1;
        while let [rest @ .., b'\r'] = line {
            line = rest;
        }
        Some(line)
    })
}

/// The data lines of a block.
pub fn data_lines(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    lines(bytes).filter(|line| !line.is_empty() && line[0] != b'#')
}

fn scan(block: Block) -> Scan {
    let decode_error = std::str::from_utf8(&block.bytes)
        .err()
        .map(|error| utf8_error_text(&block.bytes, block.offset, &error));
    let headers = lines(&block.bytes)
        .filter(|line| line.starts_with(b"#") && !line.starts_with(b"##"))
        .take(2)
        .map(|line| line[1..].to_vec())
        .collect();
    Scan {
        block,
        decode_error,
        headers,
    }
}

/// Reads the source to its end and sends it on as blocks that end at a line end; the last
/// block ends where the data does.
fn read_blocks(mut source: Source, send: impl Fn(Block) -> bool) -> Result<(), ReadFailure> {
    let mut carried = Vec::new();
    let mut offset = 0u64;
    let mut index = 0;
    let mut at_end = false;
    while !at_end {
        let mut bytes = std::mem::take(&mut carried);
        let mut filled = bytes.len();
        bytes.resize(filled + BLOCK_BYTES, 0);
        while filled < bytes.len() {
            match source.read(&mut bytes[filled..])? {
                0 => {
                    at_end = true;
                    break;
                }
                n => filled += n,
            }
        }
        bytes.truncate(filled);
        if !at_end {
            match memrchr(b'\n', &bytes) {
                Some(last) => carried = bytes.split_off(last + 1),
                // One line longer than a read: keep reading it.
                None => {
                    carried = bytes;
                    continue;
                }
            }
        }
        if bytes.is_empty() {
            continue;
        }
        let length = bytes.len() as u64;
        if !send(Block {
            index,
            offset,
            bytes,
        }) {
            return Ok(());
        }
        offset += length;
        index += 1;
    }
    Ok(())
}

impl TabInput {
    /// Reads the `--tab` output at `path`, which `display` names in messages. The failures
    /// come in the order a reader of the whole file meets them: a broken gzip stream, then
    /// the first byte sequence that is not UTF-8, then the header.
    pub fn read(path: &Path, display: &[u8]) -> Result<Self, Failure> {
        let open = || File::open(path).map_err(|e| Failure::Reported(os_error_line(&e, display)));
        let mut head = [0u8; 2];
        let mut file = open()?;
        let mut filled = 0;
        while filled < head.len() {
            match file
                .read(&mut head[filled..])
                .map_err(|e| Failure::io(&e))?
            {
                0 => break,
                n => filled += n,
            }
        }
        let source = if head[..filled] == GZIP_MAGIC {
            Source::Gzip(Box::new(GzipReader::new(BufReader::with_capacity(
                1 << 20,
                open()?,
            ))))
        } else {
            Source::Plain(open()?)
        };
        let (sender, receiver) = sync_channel::<Block>(rayon::current_num_threads());
        let (read, mut scans) = std::thread::scope(|scope| {
            let reader =
                scope.spawn(move || read_blocks(source, |block| sender.send(block).is_ok()));
            let scans: Vec<Scan> = receiver.into_iter().par_bridge().map(scan).collect();
            (reader.join().expect("the reading thread"), scans)
        });
        read.map_err(|failure| match failure {
            ReadFailure::Io(error) => Failure::io(&error),
            ReadFailure::Gzip(message) => Failure::message(message),
            ReadFailure::CorruptDeflate => Failure::corrupt_deflate(display),
        })?;
        scans.sort_by_key(|scan| scan.block.index);
        if let Some(error) = scans.iter().find_map(|scan| scan.decode_error.as_ref()) {
            return Err(Failure::message(error.clone()));
        }
        let mut headers = scans.iter().flat_map(|scan| scan.headers.iter());
        let header = headers.next();
        if header.is_some() && headers.next().is_some() {
            return Err(Failure::named(display, "two header lines"));
        }
        let Some(header) = header else {
            return Err(Failure::named(
                display,
                "no # header line; a tab output needs its header",
            ));
        };
        let columns = header.split(|b| *b == b'\t').map(<[u8]>::to_vec).collect();
        Ok(Self {
            columns,
            blocks: scans.into_iter().map(|scan| scan.block.bytes).collect(),
        })
    }

    /// The blocks of the input, each a run of whole lines.
    pub fn into_blocks(self) -> Vec<Vec<u8>> {
        self.blocks
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all(bytes: &[u8]) -> Vec<&[u8]> {
        lines(bytes).collect()
    }

    #[test]
    fn a_line_ends_at_a_newline_alone_and_loses_the_carriage_returns_before_it() {
        assert_eq!(all(b"a\nb\r\nc\r\r\n"), [&b"a"[..], b"b", b"c", b""]);
        assert_eq!(all(b"a\rb\n"), [&b"a\rb"[..], b""]);
        assert_eq!(all(b"last"), [&b"last"[..]]);
        assert_eq!(all(b""), [&b""[..]]);
        assert_eq!(all(b"\r\r"), [&b""[..]]);
    }

    #[test]
    fn data_lines_drop_the_header_the_meta_lines_and_the_empty_lines() {
        let text = b"## meta\n#Feature\tMAX_AF\r\na\t1\n\n\r\n##x\nb\t2";
        let data: Vec<&[u8]> = data_lines(text).collect();
        assert_eq!(data, [&b"a\t1"[..], b"b\t2"]);
    }

    #[test]
    fn a_scan_finds_the_header_lines_and_the_first_decode_error() {
        let found = scan(Block {
            index: 0,
            offset: 100,
            bytes: b"##m\n#a\tb\r\nrow\n#c\n#d\nx\xff\n".to_vec(),
        });
        assert_eq!(found.headers, [b"a\tb".to_vec(), b"c".to_vec()]);
        assert_eq!(
            found.decode_error.as_deref(),
            Some("'utf-8' codec can't decode byte 0xff in position 121: invalid start byte")
        );
    }
}
