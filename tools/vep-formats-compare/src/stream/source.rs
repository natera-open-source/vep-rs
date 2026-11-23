// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One output as a stream that can be read from its start more than once, cut into
//! blocks that each end at a line end.
//!
//! The blocks of a text depend on its bytes alone, so two readings of one output cut it
//! at the same places. An output is a file, or a text another reader writes anew at each
//! reading (`Text`).

use super::gzip::{Failure, Gunzip, Open, Opened};
use super::lines::cut;
use crate::argparse::path_text;
use crate::error::{Error, Result};
use crate::report::strerror;
use std::fs::File;
use std::io::{self, Cursor, Read};
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// The bytes a block is filled to before it is cut at its last line end.
pub(crate) const BLOCK: usize = 1 << 22;

/// Block buffers handed back for the next block.
pub(crate) struct Buffers {
    spare: Mutex<Vec<Vec<u8>>>,
    size: usize,
}

impl Buffers {
    pub(crate) fn new(size: usize) -> Buffers {
        Buffers {
            spare: Mutex::new(Vec::new()),
            size,
        }
    }

    /// A buffer of `size` bytes of any content.
    fn take(&self) -> Vec<u8> {
        let spare = self.spare.lock().expect("buffer pool").pop();
        spare.unwrap_or_else(|| vec![0; self.size])
    }

    pub(crate) fn give(&self, buffer: Vec<u8>) {
        // A buffer grown for a long line is dropped.
        if buffer.len() == self.size {
            self.spare.lock().expect("buffer pool").push(buffer);
        }
    }
}

/// One piece of an output: whole lines, the last one with its end unless the piece ends
/// the output.
pub(crate) struct Block {
    buffer: Vec<u8>,
    length: usize,
    /// The place of the block's first byte in the output's text.
    pub(crate) offset: u64,
    /// The block's ordinal.
    pub(crate) index: u64,
}

impl Block {
    pub(crate) fn bytes(&self) -> &[u8] {
        &self.buffer[..self.length]
    }

    pub(crate) fn into_buffer(self) -> Vec<u8> {
        self.buffer
    }
}

/// A text written anew at each reading, every reading giving the same bytes.
pub trait Text: Send + Sync {
    /// The text from its start.
    fn open(&self) -> io::Result<Box<dyn Read + Send>>;

    /// The error of a reading that failed with `detail`.
    fn failed(&self, detail: &str) -> Error;
}

/// Where an output's bytes are read from.
enum Origin {
    /// A regular file, opened anew for each reading.
    File,
    /// The bytes of a file that can be read once, as the file holds them.
    Memory(Arc<Vec<u8>>),
    /// A text its own reader writes at each reading.
    Written(Arc<dyn Text>),
}

pub(crate) struct Source {
    path: PathBuf,
    gzip: bool,
    origin: Origin,
}

impl Source {
    /// The output at `path`, gunzipped when the name ends in `.gz`. A regular file is
    /// read from the disk at each reading; any other file (a pipe, a device) is read
    /// once, here, and held whole in memory.
    pub(crate) fn open(path: &Path) -> Result<Source> {
        let mut source = Source {
            path: path.to_path_buf(),
            gzip: path.as_os_str().as_bytes().ends_with(b".gz"),
            origin: Origin::File,
        };
        let regular = std::fs::metadata(path)
            .map_err(|error| source.failed(&strerror(&error)))?
            .is_file();
        if !regular {
            let mut bytes = Vec::new();
            File::open(path)
                .and_then(|mut file| file.read_to_end(&mut bytes))
                .map_err(|error| source.failed(&strerror(&error)))?;
            source.origin = Origin::Memory(Arc::new(bytes));
        }
        Ok(source)
    }

    /// The output `text` writes, read as it is written: no reading gunzips it.
    pub(crate) fn written(text: Arc<dyn Text>) -> Source {
        Source {
            path: PathBuf::new(),
            gzip: false,
            origin: Origin::Written(text),
        }
    }

    /// The error of an output that cannot be read, worded as `input::read_text` words it;
    /// a written text words its own.
    pub(crate) fn failed(&self, detail: &str) -> Error {
        if let Origin::Written(text) = &self.origin {
            return text.failed(detail);
        }
        Error::comparator(format!("cannot read {}: {detail}", path_text(&self.path)))
    }

    /// A reading of the output from its start.
    pub(crate) fn blocks<'a>(&self, buffers: &'a Buffers) -> Result<Blocks<'a>> {
        let open: Open = match &self.origin {
            Origin::File => {
                let path = self.path.clone();
                Box::new(move || Ok(Box::new(File::open(&path)?)))
            }
            Origin::Memory(bytes) => {
                let bytes = Arc::clone(bytes);
                Box::new(move || Ok(Box::new(Cursor::new(SharedBytes(Arc::clone(&bytes))))))
            }
            Origin::Written(text) => {
                let text = Arc::clone(text);
                Box::new(move || text.open())
            }
        };
        let input = if self.gzip {
            Gunzip::new(open).map(|gunzip| Input::Gzip(Box::new(gunzip)))
        } else {
            open().map(Input::Plain)
        };
        Ok(Blocks {
            input: input.map_err(|error| self.failed(&strerror(&error)))?,
            buffers,
            carry: Vec::new(),
            offset: 0,
            index: 0,
            ended: false,
        })
    }
}

struct SharedBytes(Arc<Vec<u8>>);

impl AsRef<[u8]> for SharedBytes {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

enum Input {
    Plain(Opened),
    Gzip(Box<Gunzip>),
}

impl Input {
    /// Reads onto `out[*filled..]` until `out` is full or the text ends; `false` once
    /// the text has ended.
    fn read(&mut self, out: &mut [u8], filled: &mut usize) -> std::result::Result<bool, Failure> {
        match self {
            Input::Gzip(gunzip) => gunzip.read(out, filled),
            Input::Plain(file) => {
                while *filled < out.len() {
                    match file.read(&mut out[*filled..]) {
                        Ok(0) => return Ok(false),
                        Ok(read) => *filled += read,
                        Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                        Err(error) => return Err(Failure::System(error)),
                    }
                }
                Ok(true)
            }
        }
    }
}

/// One reading of an output, block by block.
pub(crate) struct Blocks<'a> {
    input: Input,
    buffers: &'a Buffers,
    /// The start of the line the last block was cut before.
    carry: Vec<u8>,
    offset: u64,
    index: u64,
    ended: bool,
}

impl Blocks<'_> {
    /// The next block, `None` after the last; the detail of a failure is worded as CPython
    /// words the exception.
    pub(crate) fn next(&mut self) -> std::result::Result<Option<Block>, String> {
        if self.ended && self.carry.is_empty() {
            return Ok(None);
        }
        let mut target = self.buffers.size;
        while target <= self.carry.len() {
            target *= 2;
        }
        let mut buffer = self.buffers.take();
        if buffer.len() < target {
            buffer.resize(target, 0);
        }
        let mut filled = self.carry.len();
        buffer[..filled].copy_from_slice(&self.carry);
        self.carry.clear();
        let length = loop {
            if !self.ended {
                let more = self
                    .input
                    .read(&mut buffer[..target], &mut filled)
                    .map_err(|failure| match failure {
                        Failure::System(error) => strerror(&error),
                        Failure::Text(text) => text,
                    })?;
                self.ended = !more;
            }
            if self.ended {
                break filled;
            }
            if let Some(at) = cut(&buffer[..filled]) {
                self.carry.extend_from_slice(&buffer[at..filled]);
                break at;
            }
            // A line longer than the block.
            target *= 2;
            buffer.resize(target, 0);
        };
        if length == 0 {
            self.buffers.give(buffer);
            return Ok(None);
        }
        let block = Block {
            buffer,
            length,
            offset: self.offset,
            index: self.index,
        };
        self.offset += length as u64;
        self.index += 1;
        Ok(Some(block))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stream::lines::lines;

    fn scratch(name: &str, bytes: &[u8]) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-compare-source-{}-{name}",
            std::process::id()
        ));
        std::fs::write(&path, bytes).unwrap();
        path
    }

    /// The blocks of `path` read with buffers of `size` bytes.
    fn blocks_of(path: &Path, size: usize) -> Vec<Vec<u8>> {
        let buffers = Buffers::new(size);
        let source = Source::open(path).unwrap();
        let mut blocks = source.blocks(&buffers).unwrap();
        let mut read = Vec::new();
        let mut offset = 0;
        while let Some(block) = blocks.next().unwrap() {
            assert_eq!(block.offset, offset);
            assert_eq!(block.index, read.len() as u64);
            offset += block.bytes().len() as u64;
            read.push(block.bytes().to_vec());
            buffers.give(block.into_buffer());
        }
        read
    }

    #[test]
    fn blocks_end_at_line_ends_and_hold_every_byte() {
        let mut text = Vec::new();
        for i in 0..500 {
            let end = ["\n", "\r\n", "\r"][i % 3];
            text.extend_from_slice(format!("row{i}\t{}{end}", "v".repeat(i % 40)).as_bytes());
        }
        text.extend_from_slice(b"a last line without its end");
        let path = scratch("lines", &text);
        let whole: Vec<Vec<u8>> = lines(&text).map(|line| text[line].to_vec()).collect();
        for size in [16, 64, 100, 4096, BLOCK] {
            let blocks = blocks_of(&path, size);
            assert_eq!(blocks.concat(), text, "blocks of {size} bytes");
            let cut: Vec<Vec<u8>> = blocks
                .iter()
                .flat_map(|block| lines(block).map(|line| block[line].to_vec()))
                .collect();
            assert_eq!(cut, whole, "blocks of {size} bytes");
            assert_eq!(blocks_of(&path, size), blocks, "a second reading");
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn a_line_longer_than_a_block_is_one_block() {
        let mut text = b"short\n".to_vec();
        text.extend(std::iter::repeat_n(b'x', 1000));
        text.extend_from_slice(b"\nafter\n");
        let path = scratch("long", &text);
        let blocks = blocks_of(&path, 64);
        assert_eq!(blocks.concat(), text);
        assert!(blocks.iter().any(|block| block.len() > 1000));
        let empty = scratch("empty", b"");
        assert!(blocks_of(&empty, 64).is_empty());
        for path in [path, empty] {
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn a_file_that_cannot_be_read_is_named_with_the_system_s_words() {
        let absent = std::env::temp_dir().join("vep-formats-compare-source-absent.txt");
        let message = |path: &Path| match Source::open(path) {
            Err(Error::Comparator(text)) => text,
            Ok(source) => {
                let buffers = Buffers::new(64);
                let mut blocks = source.blocks(&buffers).unwrap();
                match blocks.next() {
                    Err(detail) => match source.failed(&detail) {
                        Error::Comparator(text) => text,
                    },
                    Ok(_) => panic!("{} was read", path.display()),
                }
            }
        };
        assert_eq!(
            message(&absent),
            format!(
                "cannot read {}: No such file or directory",
                absent.display()
            )
        );
        let directory = std::env::temp_dir();
        assert_eq!(
            message(&directory),
            format!("cannot read {}: Is a directory", directory.display())
        );
    }
}
