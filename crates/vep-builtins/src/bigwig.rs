// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The one way a bigWig is opened for reading.
//!
//! `BigWigRead::open_file` prints its own `Error when opening` line on standard
//! error beside the error it returns, a line `--quiet` cannot silence and the
//! refusal built from the error repeats. [`open_bigwig`] opens the file itself
//! and hands the library the open handle, so standard error carries the
//! caller's message alone.

use std::fs::File;
use std::path::Path;

use bigtools::utils::file::reopen::ReopenableFile;
use bigtools::{BigWigRead, BigWigReadOpenError};

/// Opens the bigWig at `path`, writing nothing to standard error. The error is
/// the one `BigWigRead::open_file` returns for the same file, the open error
/// of a file that cannot be opened and the library's for one that is no
/// bigWig, so a message built from it reads the same.
pub fn open_bigwig(
    path: impl AsRef<Path>,
) -> Result<BigWigRead<ReopenableFile>, BigWigReadOpenError> {
    let path = path.as_ref();
    let file = File::open(path)?;
    BigWigRead::open(ReopenableFile {
        path: path.to_owned(),
        file,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A file that is no bigWig fails with the library's own error for the
    /// shape: the end-of-file error on a file shorter than the header.
    #[test]
    fn a_text_file_is_refused_with_the_open_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.bw");
        std::fs::write(&path, "not a bigWig\n").unwrap();
        let error = open_bigwig(&path).err().expect("a text file is no bigWig");
        assert_eq!(error.to_string(), "failed to fill whole buffer");
    }

    #[test]
    fn a_missing_file_is_refused_with_the_io_error() {
        let dir = tempfile::tempdir().unwrap();
        let error = open_bigwig(dir.path().join("none.bw"))
            .err()
            .expect("a missing file does not open");
        assert!(
            matches!(error, BigWigReadOpenError::IoError(ref io) if io.kind() == std::io::ErrorKind::NotFound),
            "{error}"
        );
    }
}
