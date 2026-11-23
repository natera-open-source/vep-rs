// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The digest a table row carries in its `sha256` column: the hash of a default-format
//! output with every line that starts with `##` left out.

use std::io::{self, BufRead};
use std::path::Path;

use sha2::{Digest, Sha256};

use crate::reader;

const META_PREFIX: &[u8] = b"##";

/// The digest of the lines `reader` yields, as 64 lower-case hex digits. A line that
/// starts with `##` is left out with its newline; every other byte is hashed as it stands.
pub fn digest_reader<R: BufRead>(mut reader: R) -> io::Result<String> {
    let mut hasher = Sha256::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if !line.starts_with(META_PREFIX) {
            hasher.update(&line);
        }
    }
    let mut hex = String::with_capacity(64);
    for byte in hasher.finalize() {
        hex.push_str(&format!("{byte:02x}"));
    }
    Ok(hex)
}

/// The digest of an output file, plain or gzip.
pub fn digest_file(path: &Path) -> io::Result<String> {
    digest_reader(reader::open(path)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    // The digest of no bytes.
    const EMPTY: &str = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";

    #[test]
    fn nothing_but_meta_lines_hashes_as_no_bytes() {
        assert_eq!(digest_reader(&b""[..]).unwrap(), EMPTY);
        assert_eq!(digest_reader(&b"## one\n## two\n"[..]).unwrap(), EMPTY);
        assert_eq!(digest_reader(&b"## no newline"[..]).unwrap(), EMPTY);
    }

    #[test]
    fn a_single_hash_line_and_an_unterminated_line_are_hashed() {
        let with_meta = digest_reader(&b"## meta\n#Uploaded_variation\tLocation\nrow"[..]).unwrap();
        let without = digest_reader(&b"#Uploaded_variation\tLocation\nrow"[..]).unwrap();
        assert_eq!(with_meta, without);
        assert_ne!(with_meta, EMPTY);
        let terminated = digest_reader(&b"#Uploaded_variation\tLocation\nrow\n"[..]).unwrap();
        assert_ne!(terminated, without);
    }
}
