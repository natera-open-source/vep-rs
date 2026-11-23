// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The one way a bigWig is opened for reading.
//!
//! `BigWigRead::open_file` prints its own `Error when opening` line on standard
//! error beside the error it returns, a line `--quiet` cannot silence and the
//! refusal built from the error repeats. `open_bigwig` opens the file itself
//! and hands the library the open handle, so standard error carries the
//! caller's message alone.
//!
//! `BigWigFile` is a file opened that way for a caller outside the crate:
//! its methods return standard-library types, so the caller names no type of
//! the library.

use std::fs::File;
use std::path::Path;

use bigtools::utils::file::reopen::ReopenableFile;
use bigtools::{BigWigRead, BigWigReadOpenError};

use crate::PluginError;

/// Opens the bigWig at `path`, writing nothing to standard error. The error is
/// the one `BigWigRead::open_file` returns for the same file, the open error
/// of a file that cannot be opened and the library's for one that is no
/// bigWig, so a message built from it reads the same.
pub(crate) fn open_bigwig(
    path: impl AsRef<Path>,
) -> Result<BigWigRead<ReopenableFile>, BigWigReadOpenError> {
    let path = path.as_ref();
    let file = File::open(path)?;
    BigWigRead::open(ReopenableFile {
        path: path.to_owned(),
        file,
    })
}

/// A bigWig opened through `open_bigwig`. Each method returns what the
/// library's reader holds or yields for the same question, copied out of the
/// reader's records into standard-library types and changed in nothing.
#[doc(hidden)]
pub struct BigWigFile {
    reader: BigWigRead<ReopenableFile>,
}

impl BigWigFile {
    /// Opens the bigWig at `path`. A file that does not open is
    /// [`PluginError::Init`] holding the text of the error `open_bigwig`
    /// returns for it.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, PluginError> {
        open_bigwig(path)
            .map(|reader| Self { reader })
            .map_err(|error| PluginError::Init(error.to_string()))
    }

    /// The file's chromosomes, each with its length in bases, in the reader's
    /// order.
    pub fn chromosome_lengths(&self) -> Vec<(String, u32)> {
        self.reader
            .chroms()
            .iter()
            .map(|chromosome| (chromosome.name.clone(), chromosome.length))
            .collect()
    }

    /// The values the reader yields for the 0-based, half-open span
    /// `start..end` of `chromosome`, each as `(start, end, value)`. A
    /// chromosome the file lacks and a read that fails are
    /// [`PluginError::Run`] holding the text of the reader's error.
    pub fn values(
        &mut self,
        chromosome: &str,
        start: u32,
        end: u32,
    ) -> Result<Vec<(u32, u32, f32)>, PluginError> {
        self.reader
            .get_interval(chromosome, start, end)
            .map_err(read_error)?
            .map(|value| {
                value
                    .map(|value| (value.start, value.end, value.value))
                    .map_err(read_error)
            })
            .collect()
    }

    /// The reduction level of each zoom level the file's header names, in the
    /// header's order.
    pub fn zoom_reduction_levels(&self) -> Vec<u32> {
        self.reader
            .info()
            .zoom_headers
            .iter()
            .map(|level| level.reduction_level)
            .collect()
    }

    /// The summary records the reader yields from the zoom level
    /// `reduction_level` for the 0-based, half-open span `start..end` of
    /// `chromosome`, each as `(start, end, bases covered, sum)`: the record's
    /// span, how many of its bases carry a value and the sum of those values.
    /// A level or a chromosome the file lacks and a read that fails are
    /// [`PluginError::Run`] holding the text of the reader's error.
    pub fn zoom_sums(
        &mut self,
        chromosome: &str,
        start: u32,
        end: u32,
        reduction_level: u32,
    ) -> Result<Vec<(u32, u32, u64, f64)>, PluginError> {
        self.reader
            .get_zoom_interval(chromosome, start, end, reduction_level)
            .map_err(read_error)?
            .map(|record| {
                record
                    .map(|record| {
                        (
                            record.start,
                            record.end,
                            record.summary.bases_covered,
                            record.summary.sum,
                        )
                    })
                    .map_err(read_error)
            })
            .collect()
    }
}

/// A reader's error as [`PluginError::Run`] holding the error's own text.
fn read_error(error: impl std::fmt::Display) -> PluginError {
    PluginError::Run(error.to_string())
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

    /// `tests/data/cons.bw` holds two chromosomes, values on both and one
    /// zoom level of 640 bases. Every reading of `BigWigFile` is the one the
    /// library's reader gives for the same file, an error included.
    #[test]
    fn a_bigwig_file_reads_what_the_reader_gives_for_the_same_file() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cons.bw");
        let mut reader = open_bigwig(&path).unwrap();
        let mut file = BigWigFile::open(&path).unwrap();

        let chromosomes: Vec<(String, u32)> = reader
            .chroms()
            .iter()
            .map(|chromosome| (chromosome.name.clone(), chromosome.length))
            .collect();
        assert_eq!(file.chromosome_lengths(), chromosomes);
        let mut by_name = chromosomes.clone();
        by_name.sort();
        assert_eq!(
            by_name,
            [
                ("CM000684.2".to_string(), 50_818_468),
                ("chr21".to_string(), 46_709_983)
            ]
        );

        let levels: Vec<u32> = reader
            .info()
            .zoom_headers
            .iter()
            .map(|level| level.reduction_level)
            .collect();
        assert_eq!(levels, [640]);
        assert_eq!(file.zoom_reduction_levels(), levels);

        for (name, length) in &chromosomes {
            for (start, end) in [(0, *length), (100, 103), (0, 1)] {
                let values: Vec<(u32, u32, f32)> = reader
                    .get_interval(name, start, end)
                    .unwrap()
                    .map(|value| {
                        let value = value.unwrap();
                        (value.start, value.end, value.value)
                    })
                    .collect();
                assert_eq!(file.values(name, start, end).unwrap(), values);
                let sums: Vec<(u32, u32, u64, f64)> = reader
                    .get_zoom_interval(name, start, end, 640)
                    .unwrap()
                    .map(|record| {
                        let record = record.unwrap();
                        (
                            record.start,
                            record.end,
                            record.summary.bases_covered,
                            record.summary.sum,
                        )
                    })
                    .collect();
                assert_eq!(file.zoom_sums(name, start, end, 640).unwrap(), sums);
            }
        }
        assert_eq!(
            file.values("chr21", 0, 46_709_983).unwrap(),
            [
                (99, 100, 10.0),
                (100, 101, 20.0),
                (101, 102, 30.0),
                (102, 103, 0.0),
                (105, 106, -1.5),
                (199, 210, 7.25),
                (300, 302, 0.1)
            ]
        );
        assert!(!file
            .zoom_sums("chr21", 0, 46_709_983, 640)
            .unwrap()
            .is_empty());

        let no_chromosome = reader.get_interval("none", 0, 1).err().unwrap();
        assert_eq!(
            file.values("none", 0, 1).unwrap_err().to_string(),
            PluginError::Run(no_chromosome.to_string()).to_string()
        );
        let no_level = reader.get_zoom_interval("chr21", 0, 1, 7).err().unwrap();
        assert_eq!(
            file.zoom_sums("chr21", 0, 1, 7).unwrap_err().to_string(),
            PluginError::Run(no_level.to_string()).to_string()
        );
        let missing = path.with_file_name("none.bw");
        let no_file = open_bigwig(&missing).err().unwrap();
        assert_eq!(
            BigWigFile::open(&missing).err().unwrap().to_string(),
            PluginError::Init(no_file.to_string()).to_string()
        );
    }
}
