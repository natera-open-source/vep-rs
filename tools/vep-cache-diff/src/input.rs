// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How an input file is opened and cut into records: gunzipped member by member as Python's
//! `gzip` module reads a file, and split at the line ends Python's universal-newlines text
//! mode recognises.

use flate2::bufread::GzDecoder;
use memchr::memchr2;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

/// The bytes a record reader takes from its stream at a time.
const READ_BUF: usize = 4 << 20;
/// The compressed bytes a gzip reader takes from its file at a time.
const GZIP_BUF: usize = 1 << 18;

/// Opens `path` as `Path.exists()` and then `open()` see it: `None` where nothing answers
/// under the name (a link whose target is gone among them), an error where the file is
/// there and cannot be opened.
fn open_existing(path: &Path) -> io::Result<Option<File>> {
    match File::open(path) {
        Ok(file) => Ok(Some(file)),
        Err(_) if std::fs::metadata(path).is_err() => Ok(None),
        Err(e) => Err(e),
    }
}

/// The whole text of a transcript shard, gunzipped when its name ends in `.gz`; `None`
/// where the file does not exist.
pub fn read_shard(path: &Path) -> io::Result<Option<Vec<u8>>> {
    let Some(mut file) = open_existing(path)? else {
        return Ok(None);
    };
    let size = file.metadata().map_or(0, |m| m.len() as usize);
    let mut text = Vec::new();
    if path.as_os_str().as_bytes().ends_with(b".gz") {
        text.reserve(size.saturating_mul(4));
        GzMembers::new(BufReader::with_capacity(GZIP_BUF, file)).read_to_end(&mut text)?;
    } else {
        text.reserve_exact(size + 1);
        file.read_to_end(&mut text)?;
    }
    Ok(Some(text))
}

/// A gzip file as a stream of its members' bytes; `None` where the file does not exist.
pub fn open_gzip(path: &Path) -> io::Result<Option<GzMembers<BufReader<File>>>> {
    Ok(open_existing(path)?.map(|file| GzMembers::new(BufReader::with_capacity(GZIP_BUF, file))))
}

/// The members of a gzip file one after the other, as Python's `gzip` module reads them: a
/// zero-byte file is an empty stream, zero bytes after a member are padding, and anything
/// else after a member has to be a further member.
pub struct GzMembers<R: BufRead> {
    state: State<R>,
}

enum State<R: BufRead> {
    Start(R),
    Member(Box<GzDecoder<R>>),
    Padding(R),
    End,
}

impl<R: BufRead> GzMembers<R> {
    pub fn new(reader: R) -> GzMembers<R> {
        GzMembers {
            state: State::Start(reader),
        }
    }
}

impl<R: BufRead> Read for GzMembers<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            match std::mem::replace(&mut self.state, State::End) {
                State::End => return Ok(0),
                State::Start(mut reader) => {
                    if !reader.fill_buf()?.is_empty() {
                        self.state = State::Member(Box::new(GzDecoder::new(reader)));
                    }
                }
                State::Padding(mut reader) => loop {
                    let available = reader.fill_buf()?;
                    if available.is_empty() {
                        break;
                    }
                    let zeros = available.iter().take_while(|b| **b == 0).count();
                    let only_zeros = zeros == available.len();
                    reader.consume(zeros);
                    if !only_zeros {
                        self.state = State::Member(Box::new(GzDecoder::new(reader)));
                        break;
                    }
                },
                State::Member(mut member) => match member.read(buf) {
                    Ok(0) => self.state = State::Padding(member.into_inner()),
                    Ok(n) => {
                        self.state = State::Member(member);
                        return Ok(n);
                    }
                    Err(e) => {
                        if e.kind() == io::ErrorKind::Interrupted {
                            self.state = State::Member(member);
                        }
                        return Err(e);
                    }
                },
            }
        }
    }
}

/// Calls `f` with each record of `reader` without its terminator. `\n`, `\r\n` and a lone
/// `\r` each end a record, as Python's universal-newlines text mode reads them; a final
/// record needs no terminator.
pub fn for_each_record<R: Read>(
    reader: R,
    f: impl FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<()> {
    for_each_record_buffered(reader, READ_BUF, f)
}

fn for_each_record_buffered<R: Read>(
    mut reader: R,
    capacity: usize,
    mut f: impl FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let mut buf = vec![0u8; capacity.max(2)];
    let mut start = 0;
    let mut end = 0;
    let mut eof = false;
    loop {
        while let Some(i) = memchr2(b'\n', b'\r', &buf[start..end]) {
            let pos = start + i;
            if buf[pos] == b'\r' && pos + 1 == end && !eof {
                break;
            }
            f(&buf[start..pos])?;
            start = pos + 1;
            if buf[pos] == b'\r' && start < end && buf[start] == b'\n' {
                start += 1;
            }
        }
        if eof {
            if start < end {
                f(&buf[start..end])?;
            }
            return Ok(());
        }
        if start > 0 {
            buf.copy_within(start..end, 0);
            end -= start;
            start = 0;
        }
        if end == buf.len() {
            buf.resize(buf.len() * 2, 0);
        }
        match reader.read(&mut buf[end..]) {
            Ok(0) => eof = true,
            Ok(n) => end += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    fn records(bytes: &[u8], capacity: usize) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        for_each_record_buffered(bytes, capacity, |r| {
            out.push(r.to_vec());
            Ok(())
        })
        .unwrap();
        out
    }

    #[test]
    fn records_end_at_lf_crlf_and_lone_cr() {
        let data = b"a\nb\r\nc\rd\r\n\r\ne";
        let expected: Vec<&[u8]> = vec![b"a", b"b", b"c", b"d", b"", b"e"];
        for capacity in [2, 3, 4, 5, 7, 4096] {
            assert_eq!(records(data, capacity), expected, "capacity {capacity}");
        }
        assert_eq!(records(b"", 4096), Vec::<Vec<u8>>::new());
        assert_eq!(records(b"\r", 4096), vec![b"".to_vec()]);
        assert_eq!(records(b"x\n", 4096), vec![b"x".to_vec()]);
    }

    fn member(text: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(text).unwrap();
        encoder.finish().unwrap()
    }

    fn gunzip(bytes: &[u8]) -> io::Result<Vec<u8>> {
        let mut out = Vec::new();
        // A three-byte buffer puts member ends, padding and headers across refills.
        GzMembers::new(BufReader::with_capacity(3, bytes)).read_to_end(&mut out)?;
        Ok(out)
    }

    #[test]
    fn every_member_is_read_and_zero_padding_after_a_member_is_skipped() {
        assert_eq!(gunzip(b"").unwrap(), b"");
        assert_eq!(gunzip(&member(b"one\n")).unwrap(), b"one\n");
        let two = [member(b"one\n"), member(b""), member(b"two\n")].concat();
        assert_eq!(gunzip(&two).unwrap(), b"one\ntwo\n");
        let padded = [member(b"one\n"), vec![0; 700], member(b"two\n"), vec![0; 5]].concat();
        assert_eq!(gunzip(&padded).unwrap(), b"one\ntwo\n");
    }

    #[test]
    fn bytes_that_are_no_member_are_refused() {
        assert!(gunzip(b"plain text").is_err());
        assert!(gunzip(&[0, 0, 0, 0]).is_err());
        let trailing = [member(b"one\n"), b"x".to_vec()].concat();
        assert!(gunzip(&trailing).is_err());
        let whole = member(b"a record that is cut short\n");
        assert!(gunzip(&whole[..whole.len() - 6]).is_err());
        let mut corrupt = whole.clone();
        let last = corrupt.len() - 5;
        corrupt[last] ^= 0xff;
        assert!(gunzip(&corrupt).is_err());
    }
}
