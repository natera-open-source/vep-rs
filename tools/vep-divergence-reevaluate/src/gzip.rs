// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A gzip file as CPython's `gzip` module reads it: the members one after another, the zero
//! bytes that pad the file after a member skipped, and the file ended wherever a member
//! would start and no byte is left, a zero-byte file included. Any other byte where a
//! member would start is a corrupt stream, zero bytes before the first member among them.

use flate2::bufread::GzDecoder;
use std::io::{self, BufRead, Read};

enum State<R> {
    /// Where a member starts or the file ends.
    Between(R),
    Member(Box<GzDecoder<R>>),
    Ended,
}

pub struct GzMembers<R> {
    state: State<R>,
}

impl<R: BufRead> GzMembers<R> {
    pub fn new(raw: R) -> GzMembers<R> {
        GzMembers {
            state: State::Between(raw),
        }
    }
}

fn skip_zero_bytes(raw: &mut impl BufRead) -> io::Result<()> {
    loop {
        let available = raw.fill_buf()?;
        let zeros = available.iter().take_while(|byte| **byte == 0).count();
        let more = zeros > 0 && zeros == available.len();
        raw.consume(zeros);
        if !more {
            return Ok(());
        }
    }
}

impl<R: BufRead> Read for GzMembers<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            // The state is put back before any return, an error's included, so a read
            // repeated after an interrupted one finds the stream where it stood.
            match std::mem::replace(&mut self.state, State::Ended) {
                State::Ended => return Ok(0),
                State::Between(mut raw) => match raw.fill_buf().map(<[u8]>::is_empty) {
                    Ok(true) => return Ok(0),
                    Ok(false) => self.state = State::Member(Box::new(GzDecoder::new(raw))),
                    Err(e) => {
                        self.state = State::Between(raw);
                        return Err(e);
                    }
                },
                State::Member(mut member) => match member.read(buf) {
                    Ok(0) => {
                        let mut raw = member.into_inner();
                        let skipped = skip_zero_bytes(&mut raw);
                        self.state = State::Between(raw);
                        skipped?;
                    }
                    read => {
                        self.state = State::Member(member);
                        return read;
                    }
                },
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    fn member(data: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    }

    /// Reads `file` whole through a buffer of `capacity` bytes, so a member's end and the
    /// padding after it fall at every place in the buffer.
    fn read_all(file: &[u8], capacity: usize) -> io::Result<Vec<u8>> {
        let mut out = Vec::new();
        GzMembers::new(io::BufReader::with_capacity(capacity, file)).read_to_end(&mut out)?;
        Ok(out)
    }

    #[test]
    fn members_are_read_one_after_another_and_zero_padding_after_a_member_is_skipped() {
        let files: [(Vec<u8>, &[u8]); 7] = [
            (Vec::new(), b""),
            (member(b""), b""),
            (member(b"one\n"), b"one\n"),
            ([member(b"one\n"), member(b"two\n")].concat(), b"one\ntwo\n"),
            ([member(b"one\n"), vec![0; 512]].concat(), b"one\n"),
            (
                [member(b"on"), vec![0; 7], member(b"e\ntwo\n"), vec![0]].concat(),
                b"one\ntwo\n",
            ),
            ([member(b""), vec![0; 3], member(b"")].concat(), b""),
        ];
        for (file, expected) in &files {
            for capacity in [1, 2, 3, 5, 8, 64, 4096] {
                assert_eq!(
                    read_all(file, capacity).unwrap(),
                    *expected,
                    "capacity {capacity}, {} bytes",
                    file.len()
                );
            }
        }
    }

    #[test]
    fn any_other_byte_where_a_member_would_start_is_a_corrupt_stream() {
        let good = member(b"one\n");
        let files = [
            vec![0u8; 4],
            [vec![0], good.clone()].concat(),
            [good.clone(), b"garbage".to_vec()].concat(),
            [good.clone(), vec![0; 9], b"x".to_vec()].concat(),
            [good.clone(), vec![0x1f]].concat(),
            good[..good.len() - 9].to_vec(),
            good[..good.len() - 1].to_vec(),
            b"#plain text\n".to_vec(),
        ];
        for file in &files {
            for capacity in [1, 4, 4096] {
                assert!(
                    read_all(file, capacity).is_err(),
                    "capacity {capacity}, {} bytes",
                    file.len()
                );
            }
        }
    }
}
