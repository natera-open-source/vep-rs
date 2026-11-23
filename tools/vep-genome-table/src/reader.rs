// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Opens a file as it is or, when it starts with the gzip magic, through a decoder that
//! reads every member, so a block-compressed file reads whole.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// A buffered reader over the file's bytes, decompressed when the file is gzip.
pub fn open(path: &Path) -> io::Result<Box<dyn BufRead>> {
    buffered(File::open(path)?)
}

/// A buffered reader over `source`'s bytes, decompressed when they start with the gzip
/// magic. The magic is read until both bytes are held, since one read of a pipe may yield
/// the first byte alone.
fn buffered<R: Read + 'static>(mut source: R) -> io::Result<Box<dyn BufRead>> {
    let mut head = [0u8; GZIP_MAGIC.len()];
    let mut held = 0;
    while held < head.len() {
        match source.read(&mut head[held..]) {
            Ok(0) => break,
            Ok(count) => held += count,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    let whole = BufReader::new(Cursor::new(head[..held].to_vec()).chain(source));
    if head[..held] == GZIP_MAGIC {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(whole))))
    } else {
        Ok(Box::new(whole))
    }
}

/// The file's lines without their newlines, the file read as UTF-8. A newline that ends
/// the file ends its last line and opens no empty one.
pub fn lines(path: &Path) -> io::Result<Vec<String>> {
    let mut text = String::new();
    open(path)?.read_to_string(&mut text)?;
    let mut lines: Vec<String> = text.split('\n').map(str::to_owned).collect();
    if lines.last().is_some_and(String::is_empty) {
        lines.pop();
    }
    Ok(lines)
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::write::GzEncoder;
    use flate2::Compression;

    use super::*;

    /// A source that yields one byte a read.
    struct ByteAtATime(Cursor<Vec<u8>>);

    impl Read for ByteAtATime {
        fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
            let one = buffer.len().min(1);
            self.0.read(&mut buffer[..one])
        }
    }

    fn read_whole(bytes: &[u8]) -> Vec<u8> {
        let mut whole = Vec::new();
        buffered(ByteAtATime(Cursor::new(bytes.to_vec())))
            .unwrap()
            .read_to_end(&mut whole)
            .unwrap();
        whole
    }

    #[test]
    fn a_gzip_source_that_yields_its_first_byte_alone_is_decompressed() {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(b"one\ntwo\n").unwrap();
        let compressed = encoder.finish().unwrap();
        assert_eq!(read_whole(&compressed), b"one\ntwo\n");
    }

    #[test]
    fn a_source_without_the_whole_magic_is_read_as_it_is() {
        for bytes in [
            &b""[..],
            &b"\x1f"[..],
            &b"\x1fone\n"[..],
            &b"one\ntwo\n"[..],
        ] {
            assert_eq!(read_whole(bytes), bytes);
        }
    }
}
