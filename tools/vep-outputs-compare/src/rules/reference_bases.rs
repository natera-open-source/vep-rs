// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The bases of a reference sequence file, asked one position at a time by the one row
//! rule that reads past a transcript's last base.
//!
//! A FASTA that is not gzipped and has its `.fai` index beside it is read through the
//! index, a base a read. Any other FASTA, gzipped (`.gz`, `.bgz`) or without an index, is
//! held whole, which suits the reference file of a corpus and not that of a genome.

use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};

/// One sequence of an indexed FASTA, as a line of the index states it: its length, the
/// offset of its first base, the bases a line holds and the bytes a line holds.
#[derive(Debug, Clone, Copy)]
struct Indexed {
    length: i64,
    offset: i64,
    line_bases: i64,
    line_bytes: i64,
}

#[derive(Debug)]
enum Held {
    Indexed(File, HashMap<String, Indexed>),
    Whole(HashMap<String, Vec<u8>>),
}

/// The reference sequences of one FASTA file.
#[derive(Debug)]
pub struct ReferenceBases {
    held: Held,
}

fn named(path: &Path, reason: impl std::fmt::Display) -> io::Error {
    io::Error::other(format!("--reference-fasta: {}: {reason}", path.display()))
}

/// The sequences of a FASTA's text by the first word of their header lines, without their
/// line ends.
fn whole(text: &[u8]) -> HashMap<String, Vec<u8>> {
    let mut sequences: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current: Option<String> = None;
    for line in text.split(|byte| *byte == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if let Some(header) = line.strip_prefix(b">") {
            let name = header
                .split(|byte| byte.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            let name = String::from_utf8_lossy(name).into_owned();
            sequences.entry(name.clone()).or_default();
            current = Some(name);
        } else if let Some(bases) = current.as_ref().and_then(|name| sequences.get_mut(name)) {
            bases.extend_from_slice(line);
        }
    }
    sequences
}

/// The sequences a FASTA index names; `None` for a line that is not five tab-separated
/// cells of which the last four are whole numbers and the fourth is above 0.
fn index(text: &str) -> Option<HashMap<String, Indexed>> {
    let mut sequences = HashMap::new();
    for line in text.lines().filter(|line| !line.is_empty()) {
        let cells: Vec<&str> = line.split('\t').collect();
        let number = |at: usize| cells.get(at)?.parse::<i64>().ok().filter(|n| *n >= 0);
        let indexed = Indexed {
            length: number(1)?,
            offset: number(2)?,
            line_bases: number(3).filter(|bases| *bases > 0)?,
            line_bytes: number(4)?,
        };
        sequences.insert(cells[0].to_string(), indexed);
    }
    Some(sequences)
}

impl ReferenceBases {
    /// The sequences of the FASTA at `path`; a failure names the file it is about.
    pub fn open(path: &Path) -> io::Result<ReferenceBases> {
        let gzipped = matches!(
            path.extension().and_then(|extension| extension.to_str()),
            Some("gz") | Some("bgz")
        );
        let mut index_path = path.as_os_str().to_os_string();
        index_path.push(".fai");
        let index_path = PathBuf::from(index_path);
        let file = File::open(path).map_err(|error| named(path, error))?;
        if !gzipped && index_path.is_file() {
            let text =
                std::fs::read_to_string(&index_path).map_err(|error| named(&index_path, error))?;
            let sequences = index(&text).ok_or_else(|| named(&index_path, "is no FASTA index"))?;
            // A base the index states past the file's end cannot be read, and a base that
            // is not read is answered as one outside its sequence.
            let size = file.metadata().map_err(|error| named(path, error))?.len();
            let within_the_file = |indexed: &Indexed| {
                let before = i128::from(indexed.length) - 1;
                let (line_bases, line_bytes) = (
                    i128::from(indexed.line_bases),
                    i128::from(indexed.line_bytes),
                );
                let last = i128::from(indexed.offset)
                    + before / line_bases * line_bytes
                    + before % line_bases;
                before < 0 || last < i128::from(size)
            };
            if !sequences.values().all(within_the_file) {
                return Err(named(
                    &index_path,
                    "states a sequence longer than the FASTA holds",
                ));
            }
            return Ok(ReferenceBases {
                held: Held::Indexed(file, sequences),
            });
        }
        let mut text = Vec::new();
        let read = if gzipped {
            MultiGzDecoder::new(file).read_to_end(&mut text)
        } else {
            let mut file = file;
            file.read_to_end(&mut text)
        };
        read.map_err(|error| named(path, error))?;
        Ok(ReferenceBases {
            held: Held::Whole(whole(&text)),
        })
    }

    /// Whether the file names `sequence`: a header line's first word, or an index line's
    /// first cell.
    pub fn holds(&self, sequence: &str) -> bool {
        match &self.held {
            Held::Whole(sequences) => sequences.contains_key(sequence),
            Held::Indexed(_, sequences) => sequences.contains_key(sequence),
        }
    }

    /// The base at the 1-based `position` of `sequence`, in upper case; `None` for a
    /// sequence the file does not name, a position outside it and a byte that is no letter.
    pub fn base(&self, sequence: &str, position: i64) -> Option<u8> {
        let base = match &self.held {
            Held::Whole(sequences) => {
                let at = usize::try_from(position.checked_sub(1)?).ok()?;
                *sequences.get(sequence)?.get(at)?
            }
            Held::Indexed(file, sequences) => {
                let indexed = sequences.get(sequence)?;
                if !(1..=indexed.length).contains(&position) {
                    return None;
                }
                let before = i128::from(position - 1);
                let (line_bases, line_bytes) = (
                    i128::from(indexed.line_bases),
                    i128::from(indexed.line_bytes),
                );
                let at = i128::from(indexed.offset)
                    + before / line_bases * line_bytes
                    + before % line_bases;
                let mut base = [0u8];
                file.read_exact_at(&mut base, u64::try_from(at).ok()?)
                    .ok()?;
                base[0]
            }
        };
        base.is_ascii_alphabetic()
            .then_some(base.to_ascii_uppercase())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TempDir;
    use std::io::Write;

    /// Two made-up sequences, the first on three lines of four bases.
    const FASTA: &str = ">7 made up\nACGT\nTGCA\nacgn\n>MT\nGG-A\n";

    fn written(name: &str, bytes: &[u8]) -> (TempDir, PathBuf) {
        let scratch = TempDir::create(&std::env::temp_dir(), "reference_bases_test_").unwrap();
        let path = scratch.path.join(name);
        std::fs::write(&path, bytes).unwrap();
        (scratch, path)
    }

    fn read_alike(reference: &ReferenceBases) {
        let bases = |sequence: &str, positions: std::ops::RangeInclusive<i64>| -> Vec<Option<u8>> {
            positions
                .map(|position| reference.base(sequence, position))
                .collect()
        };
        assert_eq!(
            bases("7", 0..=13),
            [
                None,
                Some(b'A'),
                Some(b'C'),
                Some(b'G'),
                Some(b'T'),
                Some(b'T'),
                Some(b'G'),
                Some(b'C'),
                Some(b'A'),
                Some(b'A'),
                Some(b'C'),
                Some(b'G'),
                Some(b'N'),
                None
            ]
        );
        assert_eq!(
            bases("MT", 1..=5),
            [Some(b'G'), Some(b'G'), None, Some(b'A'), None]
        );
        assert_eq!(reference.base("8", 1), None);
        assert!(reference.holds("7") && reference.holds("MT") && !reference.holds("8"));
        assert_eq!(reference.base("7", -3), None);
        assert_eq!(reference.base("7", i64::MAX), None);
    }

    #[test]
    fn a_base_is_read_alike_from_a_plain_a_gzipped_and_an_indexed_file() {
        let (_plain_scratch, plain) = written("reference.fa", FASTA.as_bytes());
        read_alike(&ReferenceBases::open(&plain).unwrap());

        let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip.write_all(FASTA.as_bytes()).unwrap();
        let (_gz_scratch, gz) = written("reference.fa.gz", &gzip.finish().unwrap());
        read_alike(&ReferenceBases::open(&gz).unwrap());

        // The index of the plain file: a sequence's length, the offset of its first base,
        // the bases and the bytes of a line.
        let (indexed_scratch, indexed) = written("reference.fa", FASTA.as_bytes());
        std::fs::write(
            indexed_scratch.path.join("reference.fa.fai"),
            "7\t12\t11\t4\t5\nMT\t4\t30\t4\t5\n",
        )
        .unwrap();
        let reference = ReferenceBases::open(&indexed).unwrap();
        assert!(matches!(reference.held, Held::Indexed(..)));
        read_alike(&reference);
    }

    #[test]
    fn a_file_that_cannot_be_read_and_an_index_that_is_none_are_named() {
        let scratch = TempDir::create(&std::env::temp_dir(), "reference_bases_test_").unwrap();
        let absent = scratch.path.join("absent.fa");
        let refused = ReferenceBases::open(&absent).unwrap_err().to_string();
        assert!(
            refused.starts_with(&format!("--reference-fasta: {}: ", absent.display())),
            "{refused}"
        );
        let (indexed_scratch, indexed) = written("reference.fa", FASTA.as_bytes());
        let index = indexed_scratch.path.join("reference.fa.fai");
        for no_index in ["7\t12\t11\t0\t5\n", "7\t12\n", "7\ttwelve\t11\t4\t5\n"] {
            std::fs::write(&index, no_index).unwrap();
            assert_eq!(
                ReferenceBases::open(&indexed).unwrap_err().to_string(),
                format!("--reference-fasta: {}: is no FASTA index", index.display()),
                "{no_index:?}"
            );
        }
    }

    /// The index states nine bases of a sequence whose one line, the file's last, holds
    /// four.
    #[test]
    fn an_index_that_states_a_sequence_longer_than_the_file_holds_is_named() {
        let (indexed_scratch, indexed) = written("reference.fa", FASTA.as_bytes());
        let index = indexed_scratch.path.join("reference.fa.fai");
        std::fs::write(&index, "7\t12\t11\t4\t5\nMT\t9\t30\t4\t5\n").unwrap();
        assert_eq!(
            ReferenceBases::open(&indexed).unwrap_err().to_string(),
            format!(
                "--reference-fasta: {}: states a sequence longer than the FASTA holds",
                index.display()
            )
        );
        // A sequence the index states with no base lies in any file.
        std::fs::write(&index, "7\t12\t11\t4\t5\nMT\t4\t30\t4\t5\n\t0\t0\t4\t5\n").unwrap();
        assert!(ReferenceBases::open(&indexed).is_ok());
    }
}
