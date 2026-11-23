// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The lines of an input VCF, plain or gzipped, as `tools/vep-outputs-compare` reads them:
//! `\n`, `\r\n` and a lone `\r` each end a line, and a last line needs no end.

use crate::argparse::path_text;
use crate::error::{Error, Result};
use crate::input::{byte_stream, text_lines_of, Stopped};
use crate::report::strerror;
use std::io;
use std::path::Path;

/// Why the reading of a file ends before its last line.
enum Stop {
    /// The file cannot be read: why, in the interpreter's words.
    Input(String),
    /// The visitor's own error.
    Visitor(Error),
}

impl From<io::Error> for Stop {
    fn from(error: io::Error) -> Stop {
        Stop::Input(strerror(&error))
    }
}

/// Whether the file is read as gzip: its name ends in `.gz` or `.bgz` after a character
/// that is no `.`, the suffix the interpreter's `Path.suffix` gives. The name is read by its
/// bytes, so one that is not UTF-8 has its suffix as any other.
fn gzipped(path: &Path) -> bool {
    let Some(name) = path.file_name().map(|name| name.as_encoded_bytes()) else {
        return false;
    };
    name.iter()
        .rposition(|byte| *byte == b'.')
        .is_some_and(|dot| {
            name[..dot].iter().any(|byte| *byte != b'.') && matches!(&name[dot..], b".gz" | b".bgz")
        })
}

fn read_lines(
    path: &Path,
    visit: &mut dyn FnMut(&[u8]) -> Result<()>,
) -> std::result::Result<(), Stop> {
    let mut reader = byte_stream(path, gzipped(path))?;
    text_lines_of(reader.as_mut(), &mut |line, _, _| visit(line)).map_err(|stopped| match stopped {
        Stopped::Read(error) => Stop::from(error),
        Stopped::Visit(error) => Stop::Visitor(error),
    })
}

/// Calls `visit` with each line of `path` without its end, the file gunzipped where its
/// name ends in `.gz` or `.bgz`. A file that cannot be opened or read is the comparator's
/// error `cannot read <path>: <why>`.
pub(super) fn for_each_line(path: &Path, mut visit: impl FnMut(&[u8]) -> Result<()>) -> Result<()> {
    read_lines(path, &mut visit).map_err(|stop| match stop {
        Stop::Input(why) => Error::comparator(format!("cannot read {}: {why}", path_text(path))),
        Stop::Visitor(error) => error,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::OsStr;
    use std::os::unix::ffi::OsStrExt;

    #[test]
    fn a_name_is_gzip_by_its_last_suffix_alone() {
        for (name, expected) in [
            ("input.vcf.gz", true),
            ("input.vcf.bgz", true),
            ("input.vcf", false),
            ("input.GZ", false),
            (".gz", false),
            ("..gz", false),
            ("a..gz", true),
            ("dir.gz/input.vcf", false),
            ("input.gz/", true),
        ] {
            assert_eq!(gzipped(Path::new(name)), expected, "{name}");
        }
        for (name, expected) in [
            (&b"in\xffput.vcf.gz"[..], true),
            (b"\xff.bgz", true),
            (b"in\xffput.vcf", false),
            (b".\xff.gz", true),
        ] {
            let path = Path::new(OsStr::from_bytes(name));
            assert_eq!(gzipped(path), expected, "{name:?}");
        }
    }

    /// A gzip member the interpreter made of a small VCF, as hexadecimal digits.
    const MEMBER: &str = "1f8b08000000000002ff558e5b0b82301886afdffd0ac1cb86aca969175d8c95265896adee45ec00453024fcf96dcd82d8c5c7fb3cfb0ebe7fbeddbbf3533f9a7e7192d92b0a38f1e5baae36d85507144bd4ab0ca254d81f4589ac28d5aa46b1cd2a9260ca18f41442e6102a47601f69af3a01b7864318aee8d724082d0e6d8382a0b6a8d1a56e58649ce991ff3476748473841d879e5998cb1ff42611631eccf2a6c5a5fdf1013a755f4d88edb8b989c484d9e77ce624a5e8865e37783c7547dea7f5b2d918010000";

    fn member() -> Vec<u8> {
        (0..MEMBER.len())
            .step_by(2)
            .map(|at| u8::from_str_radix(&MEMBER[at..at + 2], 16).unwrap())
            .collect()
    }

    /// The lines of a file of `bytes` under a gzip name, or the comparator's error with
    /// the path left out.
    fn read(name: &str, bytes: &[u8]) -> std::result::Result<Vec<Vec<u8>>, String> {
        let root = crate::cache::testing::scratch(name);
        let path = root.join("input.vcf.gz");
        std::fs::write(&path, bytes).unwrap();
        let mut lines = Vec::new();
        let read = for_each_line(&path, |line| {
            lines.push(line.to_vec());
            Ok(())
        });
        let shown = path_text(&path);
        std::fs::remove_dir_all(root).unwrap();
        match read {
            Ok(()) => Ok(lines),
            Err(Error::Comparator(text)) => Err(text
                .strip_prefix(&format!("cannot read {shown}: "))
                .unwrap_or(&text)
                .to_string()),
        }
    }

    #[test]
    fn a_broken_gzip_input_is_refused_in_the_interpreters_words() {
        let whole = member();
        let lines = read("lines-whole", &whole).unwrap();
        assert!(lines.len() > 1 && lines[0].starts_with(b"#"));
        let with = |at: usize, byte: u8| {
            let mut bytes = whole.clone();
            bytes[at] = byte;
            bytes
        };
        // A byte of the deflate stream that still inflates, to another text.
        assert_eq!(
            read("lines-stream", &with(100, 0x35)).unwrap_err(),
            "CRC check failed 0xd9b2f5a7 != 0x277fa2ad"
        );
        assert_eq!(
            read("lines-checksum", &with(195, 0x4d)).unwrap_err(),
            "CRC check failed 0xd94df5a7 != 0xd9b2f5a7"
        );
        assert_eq!(
            read("lines-length", &with(199, 0xff)).unwrap_err(),
            "Incorrect length of data produced"
        );
        // The bits of the flag byte the format reserves are read past.
        assert_eq!(read("lines-reserved", &with(3, 0xe0)).unwrap(), lines);
        assert_eq!(
            read("lines-method", &with(2, 7)).unwrap_err(),
            "Unknown compression method"
        );
        assert_eq!(
            read("lines-cut", &whole[..whole.len() - 3]).unwrap_err(),
            "Compressed file ended before the end-of-stream marker was reached"
        );
        let after = |more: &[u8]| [whole.as_slice(), more].concat();
        assert_eq!(
            read("lines-after", &after(b"no")).unwrap_err(),
            "Not a gzipped file (b'no')"
        );
        assert_eq!(
            read("lines-padded", &after(&[0, 0, 0x1f])).unwrap_err(),
            "Not a gzipped file (b'\\x1f')"
        );
        assert_eq!(read("lines-padding", &after(&[0, 0, 0])).unwrap(), lines);
        assert_eq!(
            read("lines-twice", &after(&whole)).unwrap().len(),
            2 * lines.len()
        );
        assert!(read("lines-empty", b"").unwrap().is_empty());
    }

    #[test]
    fn the_lines_before_a_wrong_trailer_are_read_before_the_file_is_refused() {
        let mut bytes = member();
        bytes[195] = 0x4d;
        let root = crate::cache::testing::scratch("lines-order");
        let path = root.join("input.vcf.gz");
        std::fs::write(&path, bytes).unwrap();
        let mut seen = 0;
        let refused = for_each_line(&path, |_| {
            seen += 1;
            Ok(())
        });
        assert!(refused.is_err());
        assert_eq!(seen, read("lines-order-whole", &member()).unwrap().len());
        // The visitor's own error on such a line is the one the reading ends with.
        let stopped = for_each_line(&path, |_| Err(Error::comparator("the visitor's")));
        assert_eq!(stopped.unwrap_err().to_string(), "the visitor's");
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn the_bytes_of_a_file_that_is_no_gzip_file_are_written_as_the_interpreter_writes_them() {
        for (bytes, shown) in [
            (&b"ab"[..], "b'ab'"),
            (&[0x1f, 0x00], "b'\\x1f\\x00'"),
            (b"'\\", "b\"'\\\\\""),
            (b"'\"", "b'\\'\"'"),
            (b"\t\xff", "b'\\t\\xff'"),
        ] {
            assert_eq!(
                read("lines-magic", bytes).unwrap_err(),
                format!("Not a gzipped file ({shown})")
            );
        }
    }
}
